!===============================================================================
! PROGRAM TDGP3D_MPI
!
! Multi-chain, parallel-tempered 3D TDGP sampler -- the 3D counterpart of
! tdgp_main.f90, with the same two-level MPI topology (TDGP_MPI_TOPOLOGY's
! SETUP_TEAMS) and control flow: rank 0 is the PT coordinator (swaps, sample
! logging); ranks 1..N_chains*nproc_per_team form N_chains teams, one per chain.
! In each team, team rank 0 is the chain driver (RJMCMC + GP) and every rank
! of the team, the driver included, solves a share of the frequencies in each
! forward call (tdgp3d_lib_interface.f90).
!
! Run with:  mpirun -np <1+N_chains*nproc_per_team> ./tdgp3d_mpi
! The run folder needs tdgp.settings (see tdgp.settings.example), MT3D_V6.dat
! (the forward mesh) and the ModEM data file named by dataFile.
!===============================================================================
PROGRAM TDGP3D_MPI
USE MPI
USE TDGP_DATA_TYPE
USE TDGP3D_COM
USE TDGP_MPI_TOPOLOGY,    ONLY: SETUP_TEAMS
USE TDGP3D_RJMCMC,        ONLY: BUILD_STARTING_MODEL3D, INIT_CHAIN_DRIVER3D, INIT_CHAIN_WORKER3D, &
                              & CHAIN_STEP_DRIVER3D, CHAIN_STEP_WORKER3D
USE TDGP3D_PT,            ONLY: SEND_NODEOBJ3D, RECV_NODEOBJ3D, TEMPSWP_MH3D, MAKE_BETA_LADDER3D, &
                              & PARALLEL_SEED3D, SAVESAMPLE3D, OPEN_RESPONSES3D
USE TDGP3D_LIB_INTERFACE, ONLY: tdgp3d_mt3d_init_team, tdgp3d_mt3d_get_grid, tdgp3d_mt3d_get_stations, &
                              & tdgp3d_mt3d_get_mesh, tdgp3d_mt3d_get_data
IMPLICIT NONE

INTEGER, PARAMETER :: PRINT_EVERY = 10   ! coordinator progress line every PRINT_EVERY iterations

TYPE(nodeobj) :: obj, objnew           ! this rank's own chain (drivers only)
TYPE(nodeobj) :: objm(2)               ! coordinator's swap-pair receive buffers
INTEGER :: rank, nproc, ierr
INTEGER :: ndat, nbuf(2)
INTEGER :: iunit_summary, iunit_ensemble, move_type
INTEGER(KIND=8) :: imcmc
REAL(RP), ALLOCATABLE :: beta_pt(:)
INTEGER :: ic, ic_start, isource1, isource2
INTEGER(KIND=8) :: ncswap, ncswapprop
LOGICAL :: swapped
INTEGER :: nsta, nfreq
REAL(RP), ALLOCATABLE :: xsta(:), ysta(:), freq(:)
REAL(RP) :: xlo, xhi, ylo, yhi, zbot, t0
REAL(RP) :: cold_logL = 0._RP, cold_chi = 0._RP, cold_acc = 0._RP, cold_m0 = 0._RP   ! latest state of the beta = 1 chain,
INTEGER  :: cold_k = 0, cold_fail = 0                               ! for the progress line

CALL MPI_INIT(ierr)
CALL MPI_COMM_RANK(MPI_COMM_WORLD, rank, ierr)
CALL MPI_COMM_SIZE(MPI_COMM_WORLD, nproc, ierr)

! Every rank reads the settings file itself:
CALL READ_TDGP_SETTINGS('tdgp.settings')

IF (kmax > NODEMX) THEN
   IF (rank == 0) WRITE(*,'(A,I0,A,I0,A)') 'FATAL: kmax=', kmax, ' exceeds the compile-time ceiling NODEMX=', &
                                          & NODEMX, ' (tdgp3d_com.f90). Lower kmax or raise NODEMX and rebuild.'
   CALL MPI_BARRIER(MPI_COMM_WORLD, ierr)
   CALL MPI_ABORT(MPI_COMM_WORLD, 1, ierr)
END IF
IF (N_chains < 1 .OR. nproc_per_team < 1) THEN
   IF (rank == 0) WRITE(*,'(A)') 'FATAL: need N_chains >= 1 and nproc_per_team >= 1 in tdgp.settings.'
   CALL MPI_BARRIER(MPI_COMM_WORLD, ierr)
   CALL MPI_ABORT(MPI_COMM_WORLD, 1, ierr)
END IF

! --- Two-level topology: team_comm (this rank's chain team, if any) and pt_comm
!     (coordinator + every chain's driver). Collective over MPI_COMM_WORLD. ---
CALL SETUP_TEAMS(nproc_per_team, N_chains, rank, team_comm, pt_comm, team_rank, pt_rank, is_driver)
! is_driver is TRUE for the PT coordinator (rank 0) and every chain driver
! (team_rank == 0); pt_rank == 0 is the coordinator.

IF (rank > 0) THEN
   ! Every non-coordinator rank is part of exactly one team: set up its forward code.
   CALL tdgp3d_mt3d_init_team(team_comm, TRIM(dataFile), nstar, ndat, rank == 1)
END IF

! --- Every team loads the same mesh and data, so chain 1's driver (pt_rank 1)
!     tells the PT coordinator the sizes: ---
IF (is_driver) THEN
   IF (pt_rank == 1) nbuf = (/ nstar, ndat /)
   CALL MPI_BCAST(nbuf, 2, MPI_INTEGER, 1, pt_comm, ierr)
   nstar = nbuf(1)
   ndat  = nbuf(2)
END IF

! --- RNG seeds and the temperature ladder (coordinator + drivers only): ---
IF (is_driver) THEN
   CALL PARALLEL_SEED3D(pt_comm, pt_rank, N_chains+1, TRIM(outputBase)//'_seed')
   ALLOCATE(beta_pt(N_chains))
   CALL MAKE_BETA_LADDER3D(N_chains, dTlog, beta_pt)
END IF

!===============================================================================
IF (pt_rank == 0 .AND. is_driver) THEN
!===================================================================== PT coordinator (global rank 0)
   CALL ALLOC_NODEOBJ(objm(1), ndat, nstar)
   CALL ALLOC_NODEOBJ(objm(2), ndat, nstar)

   OPEN(NEWUNIT=iunit_summary,  FILE=TRIM(outputBase)//'_summary.txt',  STATUS='REPLACE', ACTION='WRITE')
   OPEN(NEWUNIT=iunit_ensemble, FILE=TRIM(outputBase)//'_ensemble.txt', STATUS='REPLACE', ACTION='WRITE')
   CALL OPEN_RESPONSES3D(TRIM(outputBase)//'_responses.txt')

   ncswap     = 0_8
   ncswapprop = 0_8

   WRITE(*,'(A,I0,A,I0,A)') 'TDGP3D_MPI: coordinating ', N_chains, ' chains for ', NCHAIN, ' iterations...'

   t0 = MPI_WTIME()
   DO imcmc = 1,NCHAIN

      IF (MOD(imcmc,2_8) == 1) THEN
         ic_start = 1
      ELSE
         ! Even iterations: chain 1 gets a passthrough round, pairs (2,3)(4,5)...
         ic_start = 2
         CALL RECV_NODEOBJ3D(objm(1), 1, pt_comm)
         CALL SEND_NODEOBJ3D(objm(1), 1, pt_comm)
         IF (beta_pt(1) == 1._RP) CALL SAVE_COLD(objm(1), 1)
      END IF

      ic = ic_start
      DO WHILE (ic+1 <= N_chains)
         isource1 = ic
         isource2 = ic+1

         CALL RECV_NODEOBJ3D(objm(1), isource1, pt_comm)
         CALL RECV_NODEOBJ3D(objm(2), isource2, pt_comm)
         objm(1)%beta = beta_pt(isource1)
         objm(2)%beta = beta_pt(isource2)

         ncswapprop = ncswapprop + 1_8
         CALL TEMPSWP_MH3D(objm(1), objm(2), swapped)
         IF (swapped) ncswap = ncswap + 1_8

         CALL SEND_NODEOBJ3D(objm(1), isource1, pt_comm)
         CALL SEND_NODEOBJ3D(objm(2), isource2, pt_comm)

         IF (beta_pt(isource1) == 1._RP) CALL SAVE_COLD(objm(1), isource1)
         IF (beta_pt(isource2) == 1._RP) CALL SAVE_COLD(objm(2), isource2)

         ic = ic + 2
      END DO

      IF (ic <= N_chains) THEN   ! trailing unpaired chain
         CALL RECV_NODEOBJ3D(objm(1), ic, pt_comm)
         CALL SEND_NODEOBJ3D(objm(1), ic, pt_comm)
         IF (beta_pt(ic) == 1._RP) CALL SAVE_COLD(objm(1), ic)
      END IF

      IF (MOD(imcmc,INT(PRINT_EVERY,8)) == 0) THEN
         WRITE(*,'(A,I10,A,F6.3,A,ES13.5,A,F9.4,A,I5,A,F6.3,A,I6,A,F9.2)') ' imcmc=', imcmc, &
            & '  swap acc=', REAL(ncswap,RP)/MAX(1._RP,REAL(ncswapprop,RP)), &
            & '  beta=1: logL=', cold_logL, '  RMS=', SQRT(cold_chi/REAL(MAX(ndat,1),RP)), '  k=', cold_k, &
            & '  acc=', cold_acc, '  fwd not converged=', cold_fail, &
            & '  s/iteration=', (MPI_WTIME()-t0)/REAL(PRINT_EVERY,RP)
         IF (sample_m0) WRITE(*,'(A,F9.4,A,F10.1,A)') '            beta=1: m0=', cold_m0, ' (', 10._RP**cold_m0, ' ohm-m)'
         t0 = MPI_WTIME()
      END IF

   END DO

   CLOSE(iunit_summary)
   CLOSE(iunit_ensemble)
   WRITE(*,'(A)') 'TDGP3D_MPI: done.'

ELSE IF (is_driver) THEN
!===================================================================== chain drivers (pt_rank 1..N_chains)
   IF (ALLOCATED(xstar)) DEALLOCATE(xstar, ystar, zstar)
   ALLOCATE(xstar(nstar), ystar(nstar), zstar(nstar))
   CALL tdgp3d_mt3d_get_grid(nstar, xstar, ystar, zstar)
   CALL tdgp3d_mt3d_get_mesh(xlo, xhi, ylo, yhi, zbot, logrho_bg)

   ! Domain for the GP node locations, from the data and the forward mesh. Every
   ! driver derives the same domain; chain 1's driver reports it and writes
   ! <outputBase>_domain.txt:
   CALL tdgp3d_mt3d_get_stations(nsta, xsta, ysta, nfreq, freq)
   CALL SET_DOMAIN3D(nsta, xsta, ysta, nfreq, freq, xlo, xhi, ylo, yhi, zbot, pt_rank == 1)
   IF (pt_rank == 1) CALL REPORT_DOMAIN3D(TRIM(outputBase)//'_domain.txt')

   CALL ALLOC_NODEOBJ(obj,    ndat, nstar)
   CALL ALLOC_NODEOBJ(objnew, ndat, nstar)
   CALL tdgp3d_mt3d_get_data(ndat, obj%Dobs, obj%Dunc)
   obj%beta = beta_pt(pt_rank)

   CALL BUILD_STARTING_MODEL3D(obj)
   CALL INIT_CHAIN_DRIVER3D(obj)

   DO imcmc = 1,NCHAIN

      CALL CHAIN_STEP_DRIVER3D(obj, objnew, obj%beta, move_type)

      ! PT swap round: exactly one send-then-receive with the coordinator every
      ! iteration (the coordinator's pairing gives every chain one RECV+SEND).
      CALL SEND_NODEOBJ3D(obj, 0, pt_comm)
      CALL RECV_NODEOBJ3D(obj, 0, pt_comm)
      ! obj%beta is never sent/received -- it stays beta_pt(pt_rank).

   END DO

ELSE
!===================================================================== forward-only ranks of each team
   CALL INIT_CHAIN_WORKER3D()

   DO imcmc = 1,NCHAIN
      CALL CHAIN_STEP_WORKER3D()
   END DO

END IF
!===============================================================================

CALL MPI_BARRIER(MPI_COMM_WORLD, ierr)
CALL MPI_FINALIZE(ierr)

CONTAINS

   SUBROUTINE SAVE_COLD(o, chain_id)
   ! Logs the beta = 1 chain's sample and keeps its state for the progress line.
   TYPE(nodeobj), INTENT(IN) :: o
   INTEGER,       INTENT(IN) :: chain_id
   CALL SAVESAMPLE3D(o, chain_id, imcmc, iunit_summary, iunit_ensemble)
   cold_logL = o%logL
   cold_chi  = o%chiSquaredMisfit
   cold_k    = o%k
   cold_acc  = o%acceptance
   cold_fail = o%nfwdfail
   cold_m0   = o%m0
   END SUBROUTINE SAVE_COLD

END PROGRAM TDGP3D_MPI
