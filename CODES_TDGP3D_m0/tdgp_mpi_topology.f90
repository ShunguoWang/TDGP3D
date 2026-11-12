!===============================================================================
! MODULE TDGP_MPI_TOPOLOGY
!
! Builds the two-level MPI rank topology the TDGP sampler needs: N_chains
! independent MARE2DEM manager+worker teams (one per parallel-tempering
! chain), plus a separate, small communicator connecting the PT coordinator
! to each chain's driver for temperature-swap bookkeeping.
!
! Global rank 0            = PT coordinator (never joins a team; mirrors the
!                             1D RJMCMC code's rank-0 role -- see
!                             rjmcmc_thick_resistivity_daterr.f90's
!                             PROGRAM RJMCMC_PLANE).
! Global ranks 1..N_chains*nproc_per_team
!                           = split into N_chains contiguous teams of
!                             nproc_per_team ranks each. Within a team, the
!                             first rank (team_rank==0) is that chain's
!                             "driver": it runs the RJMCMC/GP logic (see
!                             TDGP_RJMCMC) and acts as MARE2DEM's manager for
!                             its team (see tdgp_lib_interface.f90). The
!                             other ranks in the team are unmodified MARE2DEM
!                             workers running no TDGP-specific code at all.
!
! Required total rank count: 1 + N_chains*nproc_per_team.
!
! Two communicators come out of SETUP_TEAMS:
!   team_comm -- this rank's MARE2DEM team (undefined/unused on the PT
!                coordinator). Adopted into mare2dem_mpi's module variable
!                mcomm by tdgp_mare2dem_init_team (tdgp_lib_interface.f90).
!   pt_comm   -- a THIRD, separate communicator carrying only the PT
!                coordinator (pt_rank 0) and the N_chains drivers (pt_rank
!                1..N_chains, in team order) -- used for TEMPSWP_MH2D /
!                SAVESAMPLE2D / the staggered-pairing loop (TDGP_PT,
!                tdgp_main.f90), so worker-rank MPI traffic can never collide
!                with the coordinator's MPI_RECV(...,MPI_ANY_TAG,...) calls.
!                A dedicated communicator is used here rather than reusing
!                MPI_COMM_WORLD (as the 1D code could, having no workers to
!                worry about) for exactly this reason.
!===============================================================================
MODULE TDGP_MPI_TOPOLOGY
USE MPI
IMPLICIT NONE
PRIVATE

PUBLIC :: SETUP_TEAMS

CONTAINS

!-------------------------------------------------------------------------------
SUBROUTINE SETUP_TEAMS(nproc_per_team, N_chains, rank, team_comm, pt_comm, team_rank, pt_rank, is_driver)
! Must be called collectively by EVERY rank in MPI_COMM_WORLD, including the
! PT coordinator (rank 0) -- MPI_COMM_SPLIT is collective over its parent
! communicator even for ranks passing MPI_UNDEFINED.
!-------------------------------------------------------------------------------
IMPLICIT NONE
INTEGER, INTENT(IN)  :: nproc_per_team, N_chains, rank
INTEGER, INTENT(OUT) :: team_comm, pt_comm, team_rank, pt_rank
LOGICAL, INTENT(OUT) :: is_driver

INTEGER :: team_id, pt_color, ierr, nproc, nproc_expected

CALL MPI_COMM_SIZE(MPI_COMM_WORLD, nproc, ierr)
nproc_expected = 1 + N_chains*nproc_per_team
IF (nproc /= nproc_expected) THEN
   IF (rank == 0) THEN
      WRITE(*,'(A)')            'SETUP_TEAMS: rank count mismatch.'
      WRITE(*,'(A,I0)')         '  mpirun -np gave: ', nproc
      WRITE(*,'(A,I0,A,I0,A,I0)') '  required: 1 + N_chains*nproc_per_team = 1 + ', &
                                 & N_chains, '*', nproc_per_team, ' = ', nproc_expected
      WRITE(*,'(A)')            '  (1 PT coordinator + N_chains teams of nproc_per_team ranks each)'
   END IF
   CALL MPI_BARRIER(MPI_COMM_WORLD, ierr)
   CALL MPI_ABORT(MPI_COMM_WORLD, 1, ierr)
END IF

! --- Team communicator: PT coordinator (rank 0) joins no team ---
team_id = MPI_UNDEFINED
IF (rank > 0) team_id = (rank-1)/nproc_per_team
CALL MPI_COMM_SPLIT(MPI_COMM_WORLD, team_id, rank, team_comm, ierr)

team_rank = -1
IF (rank > 0) CALL MPI_COMM_RANK(team_comm, team_rank, ierr)

is_driver = (rank == 0) .OR. (team_rank == 0)

! --- PT communicator: coordinator + every team's driver, ordered by global
!     rank so pt_rank 0 = coordinator, pt_rank 1..N_chains = chain drivers in
!     team order (global ranks 1, 1+nproc_per_team, 1+2*nproc_per_team, ...) ---
pt_color = MPI_UNDEFINED
IF (is_driver) pt_color = 1
CALL MPI_COMM_SPLIT(MPI_COMM_WORLD, pt_color, rank, pt_comm, ierr)

pt_rank = -1
IF (is_driver) CALL MPI_COMM_RANK(pt_comm, pt_rank, ierr)

END SUBROUTINE SETUP_TEAMS

END MODULE TDGP_MPI_TOPOLOGY
