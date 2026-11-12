!===============================================================================
! MODULE TDGP3D_RJMCMC: the trans-dimensional Metropolis-Hastings-Green engine
! for one 3D TDGP chain -- the 3D counterpart of TDGP_RJMCMC (tdgp_rjmcmc.f90),
! with nodes at (x, y, z) instead of (y, z) and the V6 3D forward code
! (tdgp3d_lib_interface.f90) instead of MARE2DEM.
!
! As in the 2D code, every chain step proposes and evaluates exactly ONE thing
! (a birth, a death, or one node's position or value; with sample_m0 also the
! background m0, the GP prior mean and forward half-space), and the forward call is
! COLLECTIVE over the chain's team: CHAIN_STEP_DRIVER3D does the local
! propose + prior-bounds + GP step, broadcasts a "forward call coming" flag over
! team_comm, and only then calls the forward code; CHAIN_STEP_WORKER3D is the
! matching call on every other rank of the team, once per iteration.
!
! One addition to the 2D logic: the forward code's iteration may not converge
! within iteration_max_number (MT3D_V6.dat) for a proposed model. Such a
! proposal is rejected (its misfit is not reliable) and counted in nfwdfail.
!===============================================================================
MODULE TDGP3D_RJMCMC
USE TDGP_DATA_TYPE
USE TDGP3D_COM
IMPLICIT NONE
PRIVATE

PUBLIC :: BUILD_STARTING_MODEL3D
PUBLIC :: INIT_CHAIN_DRIVER3D, INIT_CHAIN_WORKER3D
PUBLIC :: CHAIN_STEP_DRIVER3D, CHAIN_STEP_WORKER3D

CONTAINS

!-------------------------------------------------------------------------------
SUBROUTINE BUILD_STARTING_MODEL3D(obj)
! A valid starting model: kmin nodes at random positions respecting the dmin
! minimum-separation prior, each at the midpoint of the log-resistivity prior
! range (as in the 2D code). Purely local -- no MPI, no forward call.
!-------------------------------------------------------------------------------
IMPLICIT NONE
TYPE(nodeobj), INTENT(INOUT) :: obj

INTEGER  :: i, j, itry
REAL(RP) :: u, xtry, ytry, ztry
LOGICAL  :: tooClose
INTEGER, PARAMETER :: MAXTRY = 10000

obj%k  = kmin
obj%m0 = logrho_bg          ! the background of MT3D_V6.dat; the starting value if sample_m0

DO i = 1,obj%k
   tooClose = .true.
   DO itry = 1,MAXTRY
      CALL RANDOM_NUMBER(u); xtry = xmin + (xmax-xmin)*u
      CALL RANDOM_NUMBER(u); ytry = ymin + (ymax-ymin)*u
      CALL RANDOM_NUMBER(u); ztry = zmin + (zmax-zmin)*u
      tooClose = .false.
      DO j = 1,i-1
         IF ( DIST3(xtry,ytry,ztry, obj%x(j),obj%y(j),obj%z(j)) < dmin ) THEN
            tooClose = .true.
            EXIT
         END IF
      END DO
      IF (.NOT. tooClose) EXIT
   END DO
   ! If still tooClose after MAXTRY attempts, the domain/dmin/kmin combination is
   ! infeasible -- the node is placed anyway and INIT_CHAIN_DRIVER3D's GP check
   ! aborts the job cleanly if the kernel matrix is singular.
   obj%x(i)      = xtry
   obj%y(i)      = ytry
   obj%z(i)      = ztry
   obj%logrho(i) = 0.5_RP*(logrho_min+logrho_max)
END DO

END SUBROUTINE BUILD_STARTING_MODEL3D

!-------------------------------------------------------------------------------
SUBROUTINE INIT_CHAIN_DRIVER3D(obj)
! Call ONCE, before the MCMC loop, on the team driver only: GP mean of the
! starting model and one forward call (collective -- every other rank of the
! team calls INIT_CHAIN_WORKER3D once, at the same time) to set obj%logL.
!-------------------------------------------------------------------------------
USE MPI
USE TDGP3D_GP,            ONLY: GP_PREDICT_MEAN3D
USE TDGP3D_LIB_INTERFACE, ONLY: tdgp3d_mt3d_forward
IMPLICIT NONE
TYPE(nodeobj), INTENT(INOUT) :: obj

INTEGER :: ierr_gp, ierr
LOGICAL :: converged

CALL GP_PREDICT_MEAN3D(obj%x(1:obj%k), obj%y(1:obj%k), obj%z(1:obj%k), obj%logrho(1:obj%k), obj%k, &
                     & xstar, ystar, zstar, nstar, lambda_x, lambda_y, lambda_z, sigma_m, obj%m0, &
                     & obj%logrho_star, ierr_gp, lambda_increase)

IF (ierr_gp /= 0) THEN
   ! MPI_ABORT rather than STOP: a STOP here would leave this team's other
   ! ranks (and every other chain) waiting forever for a forward call.
   WRITE(*,*) 'INIT_CHAIN_DRIVER3D: FATAL -- starting model gives a non-positive-definite GP kernel matrix.'
   CALL MPI_ABORT(MPI_COMM_WORLD, 1, ierr)
END IF

CALL tdgp3d_mt3d_forward(obj%logrho_star, obj%chiSquaredMisfit, obj%Dpred, converged, logrho_bg_model=obj%m0)
obj%logL = -0.5_RP * obj%chiSquaredMisfit

IF (.NOT. converged) WRITE(*,'(A,I0,A)') ' INIT_CHAIN_DRIVER3D: WARNING -- chain ', pt_rank, ': the forward '// &
   & 'iteration did not converge for the starting model (iteration_max_number in MT3D_V6.dat); continuing'

END SUBROUTINE INIT_CHAIN_DRIVER3D

!-------------------------------------------------------------------------------
SUBROUTINE INIT_CHAIN_WORKER3D()
! Matching worker-side call for INIT_CHAIN_DRIVER3D. Workers pass empty arrays:
! only the driver's model is used and only the driver gets the predicted data.
!-------------------------------------------------------------------------------
USE TDGP3D_LIB_INTERFACE, ONLY: tdgp3d_mt3d_forward
IMPLICIT NONE
REAL(RP) :: rhoDummy(0), respDummy(0), misfitDummy
LOGICAL  :: convDummy

CALL tdgp3d_mt3d_forward(rhoDummy, misfitDummy, respDummy, convDummy)

END SUBROUTINE INIT_CHAIN_WORKER3D

!-------------------------------------------------------------------------------
SUBROUTINE CHAIN_STEP_DRIVER3D(obj, objnew, beta_mh, move_type)
! Call once per outer RJMCMC iteration, on the team driver (team rank 0) only.
! Every OTHER rank of the team must call CHAIN_STEP_WORKER3D exactly once for
! each call made here.
!-------------------------------------------------------------------------------
USE MPI
USE TDGP3D_LIB_INTERFACE, ONLY: tdgp3d_mt3d_forward
IMPLICIT NONE
TYPE(nodeobj), INTENT(INOUT) :: obj, objnew
REAL(RP),      INTENT(IN)    :: beta_mh
INTEGER,       INTENT(OUT)   :: move_type

LOGICAL :: ready, converged
INTEGER :: ierr

CALL PROPOSE_MH3D(obj, objnew, move_type, ready)

! Every rank of the team must agree on whether a forward call is coming this
! iteration -- CHAIN_STEP_WORKER3D issues the matching broadcast:
CALL MPI_BCAST(ready, 1, MPI_LOGICAL, 0, team_comm, ierr)

IF (ready) THEN
   CALL tdgp3d_mt3d_forward(objnew%logrho_star, objnew%chiSquaredMisfit, objnew%Dpred, converged, &
                          & logrho_bg_model=objnew%m0)
   IF (converged) THEN
      CALL DECIDE_MH3D(obj, objnew, beta_mh, move_type)
   ELSE
      obj%nfwdfail = obj%nfwdfail + 1
      CALL REJECT_MH3D(obj, objnew, move_type)
   END IF
ELSE
   CALL REJECT_MH3D(obj, objnew, move_type)
END IF

END SUBROUTINE CHAIN_STEP_DRIVER3D

!-------------------------------------------------------------------------------
SUBROUTINE CHAIN_STEP_WORKER3D()
! Call once per outer RJMCMC iteration on every non-driver rank of the team,
! in lockstep with that team's CHAIN_STEP_DRIVER3D call.
!-------------------------------------------------------------------------------
USE MPI
USE TDGP3D_LIB_INTERFACE, ONLY: tdgp3d_mt3d_forward
IMPLICIT NONE
LOGICAL  :: ready, convDummy
INTEGER  :: ierr
REAL(RP) :: rhoDummy(0), respDummy(0), misfitDummy

CALL MPI_BCAST(ready, 1, MPI_LOGICAL, 0, team_comm, ierr)

IF (ready) CALL tdgp3d_mt3d_forward(rhoDummy, misfitDummy, respDummy, convDummy)

END SUBROUTINE CHAIN_STEP_WORKER3D

!-------------------------------------------------------------------------------
SUBROUTINE PROPOSE_MH3D(obj, objnew, move_type, ready)
! Chooses a move (birth 1/3, death 1/3, stay 1/3, clamped at kmin/kmax as in
! the 2D code), applies it, checks the prior bounds and (only if valid) computes
! the GP mean at the earth nodes of the forward mesh. Entirely local -- no MPI.
! A stay perturbs one node's position or value (1/2 each); with sample_m0 it
! perturbs the background m0, one node's position or one node's value (1/3 each).
!
! ready = .TRUE. only if the proposal is within the prior AND the GP solve
! succeeded; objnew%logrho_star is only valid when ready is .TRUE.
!-------------------------------------------------------------------------------
USE TDGP3D_GP, ONLY: GP_PREDICT_MEAN3D
IMPLICIT NONE
TYPE(nodeobj), INTENT(INOUT) :: obj, objnew
INTEGER,       INTENT(OUT)   :: move_type
LOGICAL,       INTENT(OUT)   :: ready

REAL(RP) :: u_bd, u_move, u_m0
LOGICAL  :: outside
INTEGER  :: inode, ierr_gp

objnew = obj

! --- choose move type ---
CALL RANDOM_NUMBER(u_bd)
move_type = 1  ! stay
IF (kmin /= kmax) THEN
   IF (obj%k >= kmax) THEN
      IF (u_bd <= 1._RP/3._RP) move_type = 3   ! death only -- no room to grow
   ELSE IF (obj%k <= kmin) THEN
      IF (u_bd <= 1._RP/3._RP) move_type = 2   ! birth only -- no room to shrink
   ELSE
      IF (u_bd <= 1._RP/3._RP) THEN
         move_type = 2
      ELSE IF (u_bd > 2._RP/3._RP) THEN
         move_type = 3
      END IF
   END IF
END IF

outside = .false.

SELECT CASE (move_type)
CASE (2)   ! birth
   CALL BIRTH3D(obj, objnew, outside)
CASE (3)   ! death
   CALL DEATH3D(obj, objnew, outside)
CASE DEFAULT   ! stay: pick one active node, perturb either its position or its value (or, sample_m0, m0)
   u_m0 = 1._RP
   IF (sample_m0) CALL RANDOM_NUMBER(u_m0)    ! no extra draw without sample_m0: the same chain as before
   IF (u_m0 <= 1._RP/3._RP) THEN
      CALL PROPOSAL_M03D(obj, objnew)
   ELSE IF (obj%k < 1) THEN
      outside = .true.
   ELSE
      CALL RANDOM_NUMBER(u_move)
      inode = 1 + INT(u_move*REAL(obj%k,RP))
      inode = MIN(MAX(inode,1),obj%k)
      CALL RANDOM_NUMBER(u_move)
      IF (u_move <= 0.5_RP) THEN
         CALL PROPOSAL_POS3D(obj, objnew, inode)
      ELSE
         CALL PROPOSAL_VAL3D(obj, objnew, inode)
      END IF
   END IF
END SELECT

IF (.NOT. outside) CALL CHECKBOUNDS3D(objnew, outside)

IF (outside) THEN
   ready = .false.
   RETURN
END IF

CALL GP_PREDICT_MEAN3D(objnew%x(1:objnew%k), objnew%y(1:objnew%k), objnew%z(1:objnew%k), &
                     & objnew%logrho(1:objnew%k), objnew%k, xstar, ystar, zstar, nstar, &
                     & lambda_x, lambda_y, lambda_z, sigma_m, objnew%m0, objnew%logrho_star, ierr_gp, &
                     & lambda_increase)

! ierr_gp /= 0: the node-node kernel matrix is not positive definite (two nodes
! almost on top of one another, which the dmin check should normally prevent) --
! an invalid proposal, like a prior-bounds violation:
ready = (ierr_gp == 0)

END SUBROUTINE PROPOSE_MH3D

!-------------------------------------------------------------------------------
SUBROUTINE DECIDE_MH3D(obj, objnew, beta_mh, move_type)
! Tempered Metropolis-Hastings-Green accept/reject, given objnew with a freshly
! computed (converged) chiSquaredMisfit.
!-------------------------------------------------------------------------------
IMPLICIT NONE
TYPE(nodeobj), INTENT(INOUT) :: obj, objnew
REAL(RP),      INTENT(IN)    :: beta_mh
INTEGER,       INTENT(IN)    :: move_type

REAL(RP) :: logPLratio, u_accept

objnew%logL = -0.5_RP * objnew%chiSquaredMisfit
logPLratio  = (objnew%logL - obj%logL) * beta_mh

CALL RANDOM_NUMBER(u_accept)
IF (u_accept >= EXP(logPLratio)) THEN   ! reject
   objnew = obj
   CALL TALLY_REJECT(obj, move_type)
ELSE                                     ! accept
   obj = objnew
   CALL TALLY_ACCEPT(obj, move_type)
END IF

CALL UPDATE_ACCEPTANCE(obj)

END SUBROUTINE DECIDE_MH3D

!-------------------------------------------------------------------------------
SUBROUTINE REJECT_MH3D(obj, objnew, move_type)
! Rejection without a usable forward result (prior violation, GP failure, or a
! forward iteration that did not converge).
!-------------------------------------------------------------------------------
IMPLICIT NONE
TYPE(nodeobj), INTENT(INOUT) :: obj, objnew
INTEGER,       INTENT(IN)    :: move_type

objnew = obj
CALL TALLY_REJECT(obj, move_type)
CALL UPDATE_ACCEPTANCE(obj)

END SUBROUTINE REJECT_MH3D

!-------------------------------------------------------------------------------
SUBROUTINE TALLY_REJECT(obj, move_type)
IMPLICIT NONE
TYPE(nodeobj), INTENT(INOUT) :: obj
INTEGER,       INTENT(IN)    :: move_type
IF (move_type == 1) THEN
   obj%ireject = obj%ireject + 1
ELSE
   obj%ireject_bd = obj%ireject_bd + 1
END IF
END SUBROUTINE TALLY_REJECT

SUBROUTINE TALLY_ACCEPT(obj, move_type)
IMPLICIT NONE
TYPE(nodeobj), INTENT(INOUT) :: obj
INTEGER,       INTENT(IN)    :: move_type
IF (move_type == 1) THEN
   obj%iaccept = obj%iaccept + 1
ELSE
   obj%iaccept_bd = obj%iaccept_bd + 1
END IF
END SUBROUTINE TALLY_ACCEPT

SUBROUTINE UPDATE_ACCEPTANCE(obj)
IMPLICIT NONE
TYPE(nodeobj), INTENT(INOUT) :: obj
obj%acceptance = REAL(obj%iaccept+obj%iaccept_bd,RP) / &
               & MAX(1._RP, REAL(obj%iaccept+obj%iaccept_bd+obj%ireject+obj%ireject_bd,RP))
END SUBROUTINE UPDATE_ACCEPTANCE

!-------------------------------------------------------------------------------
SUBROUTINE BIRTH3D(obj, objnew, outside)
! Proposes one new node at a uniformly random position in the domain (uniform in
! x, y and z), with log10 resistivity drawn from the uniform prior: proposal and
! prior ratios cancel and the Jacobian is 1, as in the 2D code. Appended at index
! k+1 (node order carries no meaning).
!-------------------------------------------------------------------------------
IMPLICIT NONE
TYPE(nodeobj), INTENT(IN)    :: obj
TYPE(nodeobj), INTENT(INOUT) :: objnew
LOGICAL,       INTENT(OUT)   :: outside

REAL(RP) :: u, xnew, ynew, znew, logrhonew
INTEGER  :: i

outside  = .false.
objnew   = obj
objnew%k = obj%k + 1

CALL RANDOM_NUMBER(u); xnew = xmin + (xmax-xmin)*u
CALL RANDOM_NUMBER(u); ynew = ymin + (ymax-ymin)*u
CALL RANDOM_NUMBER(u); znew = zmin + (zmax-zmin)*u
CALL RANDOM_NUMBER(u); logrhonew = logrho_min + (logrho_max-logrho_min)*u

! Reject if the new node would be closer than dmin to any existing node:
DO i = 1,obj%k
   IF ( DIST3(xnew,ynew,znew, obj%x(i),obj%y(i),obj%z(i)) < dmin ) THEN
      outside = .true.
      RETURN
   END IF
END DO

objnew%x(objnew%k)      = xnew
objnew%y(objnew%k)      = ynew
objnew%z(objnew%k)      = znew
objnew%logrho(objnew%k) = logrhonew

END SUBROUTINE BIRTH3D

!-------------------------------------------------------------------------------
SUBROUTINE DEATH3D(obj, objnew, outside)
! Proposes removing one uniformly chosen active node; its slot is backfilled with
! the last active node (swap-remove), as in the 2D code.
!-------------------------------------------------------------------------------
IMPLICIT NONE
TYPE(nodeobj), INTENT(IN)    :: obj
TYPE(nodeobj), INTENT(INOUT) :: objnew
LOGICAL,       INTENT(OUT)   :: outside

REAL(RP) :: u
INTEGER  :: idel

outside = .false.
objnew  = obj

IF (obj%k <= 0) THEN
   outside = .true.
   RETURN
END IF

CALL RANDOM_NUMBER(u)
idel = 1 + INT(u*REAL(obj%k,RP))
idel = MIN(MAX(idel,1),obj%k)

objnew%k = obj%k - 1

IF (idel /= obj%k) THEN
   objnew%x(idel)      = obj%x(obj%k)
   objnew%y(idel)      = obj%y(obj%k)
   objnew%z(idel)      = obj%z(obj%k)
   objnew%logrho(idel) = obj%logrho(obj%k)
END IF

END SUBROUTINE DEATH3D

!-------------------------------------------------------------------------------
SUBROUTINE PROPOSAL_POS3D(obj, objnew, inode)
! Cauchy-perturbs node inode's (x, y, z) position.
!-------------------------------------------------------------------------------
IMPLICIT NONE
TYPE(nodeobj), INTENT(IN)    :: obj
TYPE(nodeobj), INTENT(INOUT) :: objnew
INTEGER,       INTENT(IN)    :: inode

REAL(RP) :: u

objnew = obj
CALL RANDOM_NUMBER(u); objnew%x(inode) = obj%x(inode) + pertsd_x*TAN(PI2*(u-0.5_RP))
CALL RANDOM_NUMBER(u); objnew%y(inode) = obj%y(inode) + pertsd_y*TAN(PI2*(u-0.5_RP))
CALL RANDOM_NUMBER(u); objnew%z(inode) = obj%z(inode) + pertsd_z*TAN(PI2*(u-0.5_RP))

END SUBROUTINE PROPOSAL_POS3D

!-------------------------------------------------------------------------------
SUBROUTINE PROPOSAL_VAL3D(obj, objnew, inode)
! Cauchy-perturbs node inode's log10-resistivity.
!-------------------------------------------------------------------------------
IMPLICIT NONE
TYPE(nodeobj), INTENT(IN)    :: obj
TYPE(nodeobj), INTENT(INOUT) :: objnew
INTEGER,       INTENT(IN)    :: inode

REAL(RP) :: u

objnew = obj
CALL RANDOM_NUMBER(u); objnew%logrho(inode) = obj%logrho(inode) + pertsd_logrho*TAN(PI2*(u-0.5_RP))

END SUBROUTINE PROPOSAL_VAL3D

!-------------------------------------------------------------------------------
SUBROUTINE PROPOSAL_M03D(obj, objnew)
! Cauchy-perturbs the background m0 (sample_m0). The Cauchy scale is drawn
! log-uniformly between pertsd_logrho/100 and pertsd_logrho for every proposal:
! the data fix the background much more tightly than a node value (a 0.005-decade
! change of a half-space raises chi2 by ~18 in the sphere test), so the large
! steps let m0 travel and the small ones let it mix. A mixture of symmetric
! steps is symmetric: the proposal ratio is 1, and with the uniform prior the
! acceptance is the likelihood ratio, as for a node value.
!-------------------------------------------------------------------------------
IMPLICIT NONE
TYPE(nodeobj), INTENT(IN)    :: obj
TYPE(nodeobj), INTENT(INOUT) :: objnew

REAL(RP) :: u, scale

objnew = obj
CALL RANDOM_NUMBER(u); scale = pertsd_logrho*10._RP**(-2._RP*u)
CALL RANDOM_NUMBER(u); objnew%m0 = obj%m0 + scale*TAN(PI2*(u-0.5_RP))

END SUBROUTINE PROPOSAL_M03D

!-------------------------------------------------------------------------------
SUBROUTINE CHECKBOUNDS3D(obj, outside)
! Every active node's (x, y, z, logrho) against the prior bounds, and every
! pair of nodes against the minimum separation dmin (Euclidean in x, y, z).
! With sample_m0, also the background m0 against logrho_min..logrho_max.
!-------------------------------------------------------------------------------
IMPLICIT NONE
TYPE(nodeobj), INTENT(IN)  :: obj
LOGICAL,       INTENT(OUT) :: outside

INTEGER :: i,j

outside = .false.

IF (sample_m0) THEN
   IF (obj%m0 < logrho_min .OR. obj%m0 > logrho_max) THEN
      outside = .true.; RETURN
   END IF
END IF

DO i = 1,obj%k
   IF (obj%x(i) < xmin .OR. obj%x(i) > xmax) THEN
      outside = .true.; RETURN
   END IF
   IF (obj%y(i) < ymin .OR. obj%y(i) > ymax) THEN
      outside = .true.; RETURN
   END IF
   IF (obj%z(i) < zmin .OR. obj%z(i) > zmax) THEN
      outside = .true.; RETURN
   END IF
   IF (obj%logrho(i) < logrho_min .OR. obj%logrho(i) > logrho_max) THEN
      outside = .true.; RETURN
   END IF
END DO

DO i = 1,obj%k-1
   DO j = i+1,obj%k
      IF ( DIST3(obj%x(i),obj%y(i),obj%z(i), obj%x(j),obj%y(j),obj%z(j)) < dmin ) THEN
         outside = .true.; RETURN
      END IF
   END DO
END DO

END SUBROUTINE CHECKBOUNDS3D

!-------------------------------------------------------------------------------
PURE REAL(RP) FUNCTION DIST3(x1,y1,z1,x2,y2,z2) RESULT(d)
IMPLICIT NONE
REAL(RP), INTENT(IN) :: x1,y1,z1,x2,y2,z2
d = SQRT((x1-x2)**2 + (y1-y2)**2 + (z1-z2)**2)
END FUNCTION DIST3

END MODULE TDGP3D_RJMCMC
