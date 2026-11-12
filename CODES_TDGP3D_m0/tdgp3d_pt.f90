!===============================================================================
! MODULE TDGP3D_PT: parallel-tempering support for the 3D TDGP sampler -- the 3D
! counterpart of TDGP_PT (tdgp_pt.f90): chain state transfer over pt_comm,
! temperature-swap acceptance, the geometric temperature ladder, RNG seeding,
! and posterior sample logging. Same logic as the 2D code, with the nodes'
! x coordinate, the nfwdfail counter and the background m0 added to the chain state.
!===============================================================================
MODULE TDGP3D_PT
USE MPI
USE TDGP_DATA_TYPE
USE TDGP3D_COM
IMPLICIT NONE
PRIVATE

PUBLIC :: SEND_NODEOBJ3D, RECV_NODEOBJ3D, TEMPSWP_MH3D, MAKE_BETA_LADDER3D, PARALLEL_SEED3D, SAVESAMPLE3D
PUBLIC :: OPEN_RESPONSES3D

INTEGER, PARAMETER :: TAG_SWAP = 100   ! base tag; SEND/RECV use TAG_SWAP+0..+6

INTEGER, PARAMETER :: RESP_EVERY = 10  ! the predicted data are saved for every RESP_EVERY-th iteration
INTEGER, SAVE      :: iunit_resp = 0   ! response file, opened by OPEN_RESPONSES3D
LOGICAL, SAVE      :: resp_open = .FALSE.

CONTAINS

!-------------------------------------------------------------------------------
SUBROUTINE SEND_NODEOBJ3D(obj, dest, comm)
! Sends the swappable part of a nodeobj to rank dest on comm: the 7 integer
! fields, the 4 scalar reals (logL, chi2, acceptance, background m0), then
! x/y/z/logrho (dim NODEMX each) and Dpred (dim obj%NDAT). Not sent: Dobs/Dunc
! (the same on every rank), logrho_star (recomputed on the next proposal) and
! beta (stays with the chain).
!-------------------------------------------------------------------------------
IMPLICIT NONE
TYPE(nodeobj), INTENT(IN) :: obj
INTEGER,       INTENT(IN) :: dest, comm

INTEGER  :: ibuf(7), ierr
REAL(RP) :: rbuf(4)

ibuf = (/ obj%k, obj%NDAT, obj%ireject_bd, obj%iaccept_bd, obj%iaccept, obj%ireject, obj%nfwdfail /)
CALL MPI_SEND(ibuf, 7, MPI_INTEGER, dest, TAG_SWAP, comm, ierr)

rbuf = (/ obj%logL, obj%chiSquaredMisfit, obj%acceptance, obj%m0 /)
CALL MPI_SEND(rbuf, 4, MPI_DOUBLE_PRECISION, dest, TAG_SWAP+1, comm, ierr)

CALL MPI_SEND(obj%x,      NODEMX,   MPI_DOUBLE_PRECISION, dest, TAG_SWAP+2, comm, ierr)
CALL MPI_SEND(obj%y,      NODEMX,   MPI_DOUBLE_PRECISION, dest, TAG_SWAP+3, comm, ierr)
CALL MPI_SEND(obj%z,      NODEMX,   MPI_DOUBLE_PRECISION, dest, TAG_SWAP+4, comm, ierr)
CALL MPI_SEND(obj%logrho, NODEMX,   MPI_DOUBLE_PRECISION, dest, TAG_SWAP+5, comm, ierr)
CALL MPI_SEND(obj%Dpred,  obj%NDAT, MPI_DOUBLE_PRECISION, dest, TAG_SWAP+6, comm, ierr)

END SUBROUTINE SEND_NODEOBJ3D

!-------------------------------------------------------------------------------
SUBROUTINE RECV_NODEOBJ3D(obj, source, comm)
! Receives what SEND_NODEOBJ3D sent. obj must already be allocated
! (ALLOC_NODEOBJ) with the correct NDAT.
!-------------------------------------------------------------------------------
IMPLICIT NONE
TYPE(nodeobj), INTENT(INOUT) :: obj
INTEGER,       INTENT(IN)    :: source, comm

INTEGER  :: ibuf(7), ierr, mpistat(MPI_STATUS_SIZE)
REAL(RP) :: rbuf(4)

CALL MPI_RECV(ibuf, 7, MPI_INTEGER, source, TAG_SWAP, comm, mpistat, ierr)
obj%k          = ibuf(1)
obj%NDAT       = ibuf(2)
obj%ireject_bd = ibuf(3)
obj%iaccept_bd = ibuf(4)
obj%iaccept    = ibuf(5)
obj%ireject    = ibuf(6)
obj%nfwdfail   = ibuf(7)

CALL MPI_RECV(rbuf, 4, MPI_DOUBLE_PRECISION, source, TAG_SWAP+1, comm, mpistat, ierr)
obj%logL             = rbuf(1)
obj%chiSquaredMisfit = rbuf(2)
obj%acceptance       = rbuf(3)
obj%m0               = rbuf(4)

CALL MPI_RECV(obj%x,      NODEMX,   MPI_DOUBLE_PRECISION, source, TAG_SWAP+2, comm, mpistat, ierr)
CALL MPI_RECV(obj%y,      NODEMX,   MPI_DOUBLE_PRECISION, source, TAG_SWAP+3, comm, mpistat, ierr)
CALL MPI_RECV(obj%z,      NODEMX,   MPI_DOUBLE_PRECISION, source, TAG_SWAP+4, comm, mpistat, ierr)
CALL MPI_RECV(obj%logrho, NODEMX,   MPI_DOUBLE_PRECISION, source, TAG_SWAP+5, comm, mpistat, ierr)
CALL MPI_RECV(obj%Dpred,  obj%NDAT, MPI_DOUBLE_PRECISION, source, TAG_SWAP+6, comm, mpistat, ierr)

END SUBROUTINE RECV_NODEOBJ3D

!-------------------------------------------------------------------------------
SUBROUTINE TEMPSWP_MH3D(obj1, obj2, swapped)
! Metropolis-Hastings temperature swap between two chains' current models:
! alpha = min(1, exp((beta2-beta1)*(logL1-logL2))). After an accepted swap only
! the model state is exchanged -- each chain keeps its own beta.
!-------------------------------------------------------------------------------
IMPLICIT NONE
TYPE(nodeobj), INTENT(INOUT) :: obj1, obj2
LOGICAL,       INTENT(OUT)   :: swapped

REAL(RP) :: logPLratio, ran_uni

CALL RANDOM_NUMBER(ran_uni)
logPLratio = (obj2%beta - obj1%beta) * (obj1%logL - obj2%logL)

swapped = (ran_uni < EXP(logPLratio))
IF (swapped) CALL SWAP_MODEL_STATE3D(obj1, obj2)

END SUBROUTINE TEMPSWP_MH3D

!-------------------------------------------------------------------------------
SUBROUTINE SWAP_MODEL_STATE3D(obj1, obj2)
! Exchanges exactly the fields SEND/RECV_NODEOBJ3D transfer.
!-------------------------------------------------------------------------------
IMPLICIT NONE
TYPE(nodeobj), INTENT(INOUT) :: obj1, obj2

INTEGER  :: itmp
REAL(RP) :: rtmp
REAL(RP), ALLOCATABLE :: atmp(:)

itmp = obj1%k;          obj1%k          = obj2%k;          obj2%k          = itmp
itmp = obj1%ireject_bd; obj1%ireject_bd = obj2%ireject_bd; obj2%ireject_bd = itmp
itmp = obj1%iaccept_bd; obj1%iaccept_bd = obj2%iaccept_bd; obj2%iaccept_bd = itmp
itmp = obj1%iaccept;    obj1%iaccept    = obj2%iaccept;    obj2%iaccept    = itmp
itmp = obj1%ireject;    obj1%ireject    = obj2%ireject;    obj2%ireject    = itmp
itmp = obj1%nfwdfail;   obj1%nfwdfail   = obj2%nfwdfail;   obj2%nfwdfail   = itmp

rtmp = obj1%logL;             obj1%logL             = obj2%logL;             obj2%logL             = rtmp
rtmp = obj1%chiSquaredMisfit; obj1%chiSquaredMisfit = obj2%chiSquaredMisfit; obj2%chiSquaredMisfit = rtmp
rtmp = obj1%acceptance;       obj1%acceptance       = obj2%acceptance;       obj2%acceptance       = rtmp
rtmp = obj1%m0;               obj1%m0               = obj2%m0;               obj2%m0               = rtmp

atmp = obj1%x;      obj1%x      = obj2%x;      obj2%x      = atmp
atmp = obj1%y;      obj1%y      = obj2%y;      obj2%y      = atmp
atmp = obj1%z;      obj1%z      = obj2%z;      obj2%z      = atmp
atmp = obj1%logrho; obj1%logrho = obj2%logrho; obj2%logrho = atmp
atmp = obj1%Dpred;  obj1%Dpred  = obj2%Dpred;  obj2%Dpred  = atmp

END SUBROUTINE SWAP_MODEL_STATE3D

!-------------------------------------------------------------------------------
SUBROUTINE MAKE_BETA_LADDER3D(N_chains_in, dTlog_in, beta_pt)
! Geometric ladder: beta_pt(1) = 1 (the true posterior), beta_pt(i) = 1/dTlog**(i-1).
!-------------------------------------------------------------------------------
IMPLICIT NONE
INTEGER,  INTENT(IN)  :: N_chains_in
REAL(RP), INTENT(IN)  :: dTlog_in
REAL(RP), INTENT(OUT) :: beta_pt(N_chains_in)

INTEGER :: i

DO i = 1,N_chains_in
   beta_pt(i) = 1._RP / dTlog_in**(i-1)
END DO

END SUBROUTINE MAKE_BETA_LADDER3D

!-------------------------------------------------------------------------------
SUBROUTINE PARALLEL_SEED3D(pt_comm_in, pt_rank_in, n_pt_ranks, iseedfile_base)
! Same as the 2D code's PARALLEL_SEED2D: the PT coordinator draws one seed
! vector per pt_comm rank, broadcasts them, and every rank seeds its own RNG
! from its column (and writes it to <iseedfile_base>_<pt_rank>.log).
!-------------------------------------------------------------------------------
IMPLICIT NONE
INTEGER,          INTENT(IN) :: pt_comm_in, pt_rank_in, n_pt_ranks
CHARACTER(LEN=*), INTENT(IN) :: iseedfile_base

INTEGER :: iseedsize, i, ierr, iunit
INTEGER, ALLOCATABLE  :: iseeds(:,:), iseed(:)
REAL(RP), ALLOCATABLE :: rseeds(:,:)
REAL(RP) :: junk
CHARACTER(LEN=16) :: rankstr

CALL RANDOM_SEED()
CALL RANDOM_SEED(SIZE=iseedsize)
ALLOCATE(iseed(iseedsize), iseeds(iseedsize,n_pt_ranks))

IF (pt_rank_in == 0) THEN
   ALLOCATE(rseeds(iseedsize,n_pt_ranks))
   CALL RANDOM_NUMBER(rseeds)
   iseeds = -NINT(rseeds*1.0e6_RP)
   DEALLOCATE(rseeds)
END IF

DO i = 1,iseedsize
   CALL MPI_BCAST(iseeds(i,:), n_pt_ranks, MPI_INTEGER, 0, pt_comm_in, ierr)
END DO

iseed = iseeds(:,pt_rank_in+1)

WRITE(rankstr,'(I0)') pt_rank_in
OPEN(NEWUNIT=iunit, FILE=TRIM(iseedfile_base)//'_'//TRIM(rankstr)//'.log', STATUS='REPLACE', ACTION='WRITE')
WRITE(iunit,*) 'pt_rank ', pt_rank_in, ' seed ', iseed
CLOSE(iunit)

CALL RANDOM_SEED(PUT=iseed)

DO i = 1,100   ! burn-in draws
   CALL RANDOM_NUMBER(junk)
END DO

DEALLOCATE(iseeds, iseed)

END SUBROUTINE PARALLEL_SEED3D

!-------------------------------------------------------------------------------
SUBROUTINE SAVESAMPLE3D(obj, chain_id, imcmc, iunit_summary, iunit_ensemble)
! Appends one row to the summary file (imcmc, chain, logL, chi-squared misfit,
! k, acceptance, nfwdfail; with sample_m0 also m0) and one row to the ensemble
! file (imcmc, chain, k, then k node quadruples x, y, z, log10 rho; with
! sample_m0 then m0, so a row has 3+4k+1 values). Call only for the chain at
! beta = 1. Also saves the predicted data if OPEN_RESPONSES3D was called (see
! SAVERESPONSE3D).
!-------------------------------------------------------------------------------
IMPLICIT NONE
TYPE(nodeobj),   INTENT(IN) :: obj
INTEGER,         INTENT(IN) :: chain_id
INTEGER(KIND=8), INTENT(IN) :: imcmc
INTEGER,         INTENT(IN) :: iunit_summary, iunit_ensemble

INTEGER :: i

IF (sample_m0) THEN
   WRITE(iunit_summary,'(I12,1X,I6,1X,ES16.8,1X,ES16.8,1X,I6,1X,F8.4,1X,I10,1X,ES16.8)') &
       & imcmc, chain_id, obj%logL, obj%chiSquaredMisfit, obj%k, obj%acceptance, obj%nfwdfail, obj%m0
ELSE
   WRITE(iunit_summary,'(I12,1X,I6,1X,ES16.8,1X,ES16.8,1X,I6,1X,F8.4,1X,I10)') &
       & imcmc, chain_id, obj%logL, obj%chiSquaredMisfit, obj%k, obj%acceptance, obj%nfwdfail
END IF
FLUSH(iunit_summary)

WRITE(iunit_ensemble,'(I12,1X,I6,1X,I6)',ADVANCE='NO') imcmc, chain_id, obj%k
DO i = 1,obj%k
   WRITE(iunit_ensemble,'(4(1X,ES14.6))',ADVANCE='NO') obj%x(i), obj%y(i), obj%z(i), obj%logrho(i)
END DO
IF (sample_m0) WRITE(iunit_ensemble,'(1X,ES16.8)',ADVANCE='NO') obj%m0
WRITE(iunit_ensemble,*)
FLUSH(iunit_ensemble)

CALL SAVERESPONSE3D(obj, imcmc)

END SUBROUTINE SAVESAMPLE3D

!-------------------------------------------------------------------------------
SUBROUTINE OPEN_RESPONSES3D(filename)
! Opens the file in which SAVERESPONSE3D saves the predicted data of the cold
! chain. Call once, before the first SAVESAMPLE3D.
!-------------------------------------------------------------------------------
IMPLICIT NONE
CHARACTER(LEN=*), INTENT(IN) :: filename

OPEN(NEWUNIT=iunit_resp, FILE=TRIM(filename), STATUS='REPLACE', ACTION='WRITE')
resp_open = .TRUE.

END SUBROUTINE OPEN_RESPONSES3D

!-------------------------------------------------------------------------------
SUBROUTINE SAVERESPONSE3D(obj, imcmc)
! For every RESP_EVERY-th iteration appends one row to the response file: the
! iteration number, then the predicted data (real, imaginary part of each row of
! the ModEM data file, in its order, units and sign convention).
!-------------------------------------------------------------------------------
IMPLICIT NONE
TYPE(nodeobj),   INTENT(IN) :: obj
INTEGER(KIND=8), INTENT(IN) :: imcmc

IF (.NOT. resp_open) RETURN
IF (MOD(imcmc, INT(RESP_EVERY,8)) /= 0) RETURN

WRITE(iunit_resp,'(I12,*(1X,ES13.5))') imcmc, obj%Dpred(1:obj%NDAT)
FLUSH(iunit_resp)

END SUBROUTINE SAVERESPONSE3D

END MODULE TDGP3D_PT
