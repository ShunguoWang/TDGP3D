!===============================================================================
! MODULE TDGP3D_COM: shared state for the 3D TDGP trans-dimensional Bayesian
! sampler -- the 3D counterpart of TDGP_COM (tdgp_com.f90): the nodeobj type (one
! chain's complete state), prior bounds, GP hyperparameters, parallel-tempering
! schedule, MPI topology handles, and the tdgp.settings namelist reader.
!
! Coordinates are those of the 3D forward mesh (MT3D_V6.dat): x north, y east,
! z depth [m], positive down, z = 0 at the surface.
!===============================================================================
MODULE TDGP3D_COM
USE TDGP_DATA_TYPE
IMPLICIT NONE

! --- Trans-dimensional node-count ceiling (compile-time max; active count is
!     tracked separately by nodeobj%k) ---
INTEGER, PARAMETER :: NODEMX = 10000

! --- nodeobj: complete state of one MCMC chain ---
TYPE :: nodeobj
   INTEGER  :: k = 0, NDAT = 0                                     ! # active GP nodes (trans-D), # data
   INTEGER  :: ireject_bd = 0, iaccept_bd = 0, iaccept = 0, ireject = 0
   INTEGER  :: nfwdfail = 0                     ! proposals rejected because the forward iteration did not converge
   REAL(RP), ALLOCATABLE :: x(:), y(:)          ! dim NODEMX -- node position [m]
   REAL(RP), ALLOCATABLE :: z(:)                ! dim NODEMX -- node depth [m]
   REAL(RP), ALLOCATABLE :: logrho(:)           ! dim NODEMX -- node log10-resistivity
   REAL(RP), ALLOCATABLE :: logrho_star(:)      ! dim nstar  -- GP mean at the earth nodes of the forward mesh
   REAL(RP), ALLOCATABLE :: Dobs(:), Dunc(:)    ! dim NDAT   -- observed data + uncertainty
   REAL(RP), ALLOCATABLE :: Dpred(:)            ! dim NDAT   -- predicted data
   REAL(RP) :: beta = 1._RP, logL = 0._RP, chiSquaredMisfit = 0._RP
   REAL(RP) :: acceptance = 0._RP
   REAL(RP) :: m0 = 0._RP                       ! log10 rho of the background: the GP prior mean and the forward
                                                ! code's half-space (fixed at logrho_bg unless sample_m0)
END TYPE nodeobj

! --- Prior bounds (overridable via tdgp.settings) ---
INTEGER  :: kmin = 8, kmax = 150
REAL(RP) :: logrho_min = -1._RP, logrho_max = 5._RP  ! log10(ohm-m) prior bounds (0.1 -- 100000 ohm-m)

! --- Domain for the GP node locations (nodes can only be located inside it). Not in tdgp.settings:
!     SET_DOMAIN3D derives it from the data and the forward mesh. ---
REAL(RP) :: xmin = -3000._RP, xmax = 3000._RP        ! north [m]
REAL(RP) :: ymin = -3000._RP, ymax = 3000._RP        ! east [m]
REAL(RP) :: zmin = 0._RP,     zmax = 3000._RP        ! depth [m]
REAL(RP) :: rho_ref = 100._RP                        ! [ohm-m] resistivity of the skin depth that sets the domain depth
REAL(RP) :: xy_margin_frac = 0.5_RP                  ! x and y margin on each side = max(xy_margin_frac*(receiver span), xy_margin_min)
REAL(RP) :: xy_margin_min  = 5000._RP                ! [m] floor on that margin
REAL(RP) :: z_depth_min    = 0._RP                   ! [m] domain depth zmax-zmin if > 0 (replaces the skin depth)

! --- GP hyperparameters ---
REAL(RP) :: lambda_x = 1000._RP, lambda_y = 1000._RP, lambda_z = 500._RP   ! kernel length scales [m]
! lambda_increase > 0: all three length scales grow with the depth d below the surface,
! lambda + lambda_increase*d [m per m of depth], so deeper nodes reach farther in x, y and z
! (kernel in TDGP3D_GP). 0 (default) = constant.
REAL(RP) :: lambda_increase = 0._RP
REAL(RP), PARAMETER :: sigma_m   = 0.05_RP           ! nugget, log10 rho units
REAL(RP), PARAMETER :: DMIN_FRAC = 0.25_RP           ! minimum node separation = DMIN_FRAC*min(lambda_x, lambda_y, lambda_z)
REAL(RP) :: dmin = 0._RP                             ! [m], in x, y and z; set by READ_TDGP_SETTINGS
! GP prior mean: log10 of the background resistivity of the forward mesh (1/background_conductivity
! in MT3D_V6.dat) -- the model is the background wherever no node is near. Not a setting: every
! chain driver sets it at startup from the forward mesh.
REAL(RP) :: logrho_bg = 2._RP
! sample_m0 = .true.: the background m0 (GP prior mean and forward half-space, nodeobj%m0) is a parameter of
! every model, sampled with the nodes from logrho_bg on (uniform prior logrho_min..logrho_max).
! .false. (default): m0 = logrho_bg throughout, as before.
LOGICAL  :: sample_m0 = .false.

! --- Proposal widths (Cauchy scale) ---
REAL(RP) :: pertsd_x = 0._RP, pertsd_y = 0._RP, pertsd_z = 0._RP, pertsd_logrho = 0._RP
REAL(RP) :: pertsdsc = 30._RP
! pertsd_* = (range of x, y, z, log10 rho)/pertsdsc, set by UPDATE_PROPOSAL_WIDTHS --
! larger pertsdsc => narrower proposals => higher acceptance, slower mixing.

! --- Parallel tempering (geometric ladder; chain 1 is the only one at beta = 1) ---
REAL(RP) :: dTlog = 1.15_RP
INTEGER  :: N_chains = 8

! --- MPI topology (filled in by TDGP_MPI_TOPOLOGY%SETUP_TEAMS; team_comm = MPI_COMM_WORLD in tdgp3d_single) ---
INTEGER :: nproc_per_team = 6
INTEGER :: team_comm = -1, pt_comm = -1
INTEGER :: team_rank = -1, pt_rank = -1
LOGICAL :: is_driver = .false.

! --- The earth nodes of the forward mesh (z >= 0), where the GP mean is evaluated -- the 3D
!     counterpart of the MARE2DEM free-region centroids. Filled once on every chain driver by
!     tdgp3d_mt3d_get_grid (tdgp3d_lib_interface.f90). ---
INTEGER :: nstar = 0
REAL(RP), ALLOCATABLE :: xstar(:), ystar(:), zstar(:)

! --- Run length / I/O ---
INTEGER(KIND=8)    :: NCHAIN = 100000_8
CHARACTER(LEN=256) :: dataFile   = 'data.dat'      ! ModEM data file; the forward mesh is always MT3D_V6.dat
CHARACTER(LEN=256) :: outputBase = 'tdgp3d_out'

CONTAINS

!-------------------------------------------------------------------------------
SUBROUTINE READ_TDGP_SETTINGS(filename)
! Reads run parameters from a Fortran namelist file (default values above
! apply to anything the file omits, or if the file itself is missing). An
! unknown or malformed setting stops the run.
!-------------------------------------------------------------------------------
IMPLICIT NONE
CHARACTER(LEN=*), INTENT(IN) :: filename
INTEGER :: iunit, ios
CHARACTER(LEN=256) :: msg

NAMELIST /tdgp_settings/ kmin, kmax, logrho_min, logrho_max, rho_ref, &
                        & xy_margin_frac, xy_margin_min, z_depth_min, &
                        & lambda_x, lambda_y, lambda_z, lambda_increase, pertsdsc, &
                        & dTlog, N_chains, nproc_per_team, &
                        & NCHAIN, dataFile, outputBase, sample_m0

OPEN(NEWUNIT=iunit, FILE=TRIM(filename), STATUS='OLD', ACTION='READ', IOSTAT=ios)
IF (ios /= 0) THEN
   WRITE(*,*) 'READ_TDGP_SETTINGS: could not open ', TRIM(filename), ' -- using built-in defaults'
ELSE
   READ(iunit, NML=tdgp_settings, IOSTAT=ios, IOMSG=msg)
   CLOSE(iunit)
   IF (ios /= 0) THEN
      WRITE(*,*) 'READ_TDGP_SETTINGS: cannot read &tdgp_settings from ', TRIM(filename), ': ', TRIM(msg)
      WRITE(*,*) '  Unknown or misspelled setting? See tdgp.settings.example for the complete list.'
      STOP 1
   END IF
END IF

IF (lambda_increase < 0._RP) THEN
   WRITE(*,*) 'READ_TDGP_SETTINGS: lambda_increase must be >= 0'
   STOP 1
END IF

dmin = DMIN_FRAC*MIN(lambda_x, lambda_y, lambda_z)
CALL UPDATE_PROPOSAL_WIDTHS()      ! recomputed by SET_DOMAIN3D once the domain is known

END SUBROUTINE READ_TDGP_SETTINGS

!-------------------------------------------------------------------------------
SUBROUTINE UPDATE_PROPOSAL_WIDTHS()
! pertsd_* = (prior range)/pertsdsc.
!-------------------------------------------------------------------------------
IMPLICIT NONE

pertsd_x      = (xmax-xmin)/pertsdsc
pertsd_y      = (ymax-ymin)/pertsdsc
pertsd_z      = (zmax-zmin)/pertsdsc
pertsd_logrho = (logrho_max-logrho_min)/pertsdsc

END SUBROUTINE UPDATE_PROPOSAL_WIDTHS

!-------------------------------------------------------------------------------
SUBROUTINE SET_DOMAIN3D(nsta, xsta, ysta, nfreq, freq, xlo, xhi, ylo, yhi, zbot, report)
! Sets xmin/xmax/ymin/ymax/zmin/zmax, the domain in which the GP nodes may be
! located, i.e. where the model is allowed to depart from the background:
!
!   x, y : from the first MT receiver to the last, plus a margin on each side of
!          max( xy_margin_frac*(receiver span), xy_margin_min ), each direction on
!          its own, and never outside the forward mesh (xlo..xhi, ylo..yhi)
!   z    : from the surface (z = 0) down z_depth_min [m] if it is given (> 0),
!          otherwise one skin depth of the lowest frequency in the data,
!          503.3*sqrt(rho_ref/f_min) [m]; never below the bottom of the forward mesh (zbot)
!
! Receiver depths are not used. Identical on every chain driver. report = .true.
! prints what was derived (call it with .true. on one driver only).
!-------------------------------------------------------------------------------
USE MPI
IMPLICIT NONE
INTEGER,  INTENT(IN) :: nsta, nfreq
REAL(RP), INTENT(IN) :: xsta(nsta), ysta(nsta), freq(nfreq)
REAL(RP), INTENT(IN) :: xlo, xhi, ylo, yhi, zbot
LOGICAL,  INTENT(IN) :: report

REAL(RP), PARAMETER :: SKIN_M = 503.3_RP     ! skin depth [m] = 503.3*sqrt(rho/f)
REAL(RP) :: x_lo, x_hi, y_lo, y_hi, margx, margy, depth, fmin
INTEGER  :: ierr

IF (nsta < 1 .OR. nfreq < 1) THEN
   IF (report) WRITE(*,'(A)') ' TDGP3D: FATAL -- the domain needs the MT receivers and frequencies of the data file'
   CALL MPI_ABORT(MPI_COMM_WORLD, 1, ierr)
   RETURN
END IF

x_lo  = MINVAL(xsta);  x_hi = MAXVAL(xsta)
y_lo  = MINVAL(ysta);  y_hi = MAXVAL(ysta)
margx = MAX(xy_margin_frac*(x_hi-x_lo), xy_margin_min)
margy = MAX(xy_margin_frac*(y_hi-y_lo), xy_margin_min)
fmin  = MINVAL(freq)
depth = SKIN_M*SQRT(rho_ref/fmin)
IF (z_depth_min > 0._RP) depth = z_depth_min          ! given in tdgp.settings: replaces the skin depth

xmin = MAX(x_lo-margx, xlo)
xmax = MIN(x_hi+margx, xhi)
ymin = MAX(y_lo-margy, ylo)
ymax = MIN(y_hi+margy, yhi)
zmin = 0._RP
zmax = MIN(zmin + depth, zbot)

CALL UPDATE_PROPOSAL_WIDTHS()

IF (report) THEN
   WRITE(*,'(A)') ' TDGP3D: domain for the GP node locations:'
   WRITE(*,'(A,I0,A,2F12.1,A,F10.1,A)') '   x: ', nsta, ' MT receivers at', x_lo, x_hi, ' m, margin', margx, ' m on each side'
   WRITE(*,'(A,2F12.1,A,F10.1,A)')      '   y:    MT receivers at', y_lo, y_hi, ' m, margin', margy, ' m on each side'
   IF (xmin > x_lo-margx .OR. xmax < x_hi+margx .OR. ymin > y_lo-margy .OR. ymax < y_hi+margy) &
      & WRITE(*,'(A)') '   (x/y limited to the forward mesh)'
   IF (z_depth_min > 0._RP) THEN
      WRITE(*,'(A,F10.1,A)') '   z: 0 down z_depth_min =', depth, ' m (tdgp.settings; the skin depth is not used)'
   ELSE
      WRITE(*,'(A,ES10.3,A,F7.1,A,F10.1,A)') '   z: 0 down one skin depth at the lowest frequency (', fmin, ' Hz, rho_ref =', &
         & rho_ref, ' ohm-m) =', depth, ' m'
   END IF
   IF (zmax < zmin + depth) WRITE(*,'(A)') '   (z limited to the bottom of the forward mesh)'
END IF

END SUBROUTINE SET_DOMAIN3D

!-------------------------------------------------------------------------------
SUBROUTINE REPORT_DOMAIN3D(filename)
! Prints the domain actually in use and how many of the forward mesh's earth
! nodes lie inside it (outside it, and far from every node, the model is the
! background resistivity), and writes the domain to <filename> as a namelist so
! post-processing can read what this run used. Call on ONE chain driver only,
! after xstar/ystar/zstar are filled.
!-------------------------------------------------------------------------------
IMPLICIT NONE
CHARACTER(LEN=*), INTENT(IN) :: filename
INTEGER  :: iunit, nin, ios
REAL(RP) :: nlx, nly, nlz

nin = 0
IF (nstar > 0) THEN
   nin = COUNT( xstar(1:nstar) >= xmin .AND. xstar(1:nstar) <= xmax .AND. &
              & ystar(1:nstar) >= ymin .AND. ystar(1:nstar) <= ymax .AND. &
              & zstar(1:nstar) >= zmin .AND. zstar(1:nstar) <= zmax )
END IF

WRITE(*,'(A)')                ' TDGP3D domain in use (GP nodes can only be located inside it):'
WRITE(*,'(A,2F14.1,A)')       '   x = ', xmin, xmax, ' m'
WRITE(*,'(A,2F14.1,A)')       '   y = ', ymin, ymax, ' m'
WRITE(*,'(A,2F14.1,A)')       '   z = ', zmin, zmax, ' m'
WRITE(*,'(A,I0,A,I0,A)')      '   earth nodes of the forward mesh inside the domain: ', nin, ' of ', nstar, &
                             & ' (the GP mean is evaluated at all of them)'
WRITE(*,'(A,F8.3,A,F10.1,A)') ' TDGP3D GP prior mean (background of MT3D_V6.dat): log10(rho) =', logrho_bg, &
                             & ' (', 10._RP**logrho_bg, ' ohm-m)'
IF (logrho_bg < logrho_min .OR. logrho_bg > logrho_max) &
   & WRITE(*,'(A)') ' WARNING: the background resistivity is outside the prior range logrho_min..logrho_max'
IF (sample_m0) WRITE(*,'(A)') '   sample_m0: the background m0 is sampled with the nodes, starting from this value '// &
   & '(uniform prior logrho_min..logrho_max)'

IF (lambda_increase > 0._RP) THEN
   WRITE(*,'(A,F7.3,A)') ' TDGP3D length scales grow with depth: lambda + lambda_increase*depth, lambda_increase =', &
      & lambda_increase, ' m/m'
   WRITE(*,'(A,3F9.1,A,3F9.1,A)') '   lambda_x, lambda_y, lambda_z =', lambda_x, lambda_y, lambda_z, ' m at the surface ->', &
      & (/ lambda_x, lambda_y, lambda_z /) + lambda_increase*MAX(zmax, 0._RP), ' m at the bottom of the domain'
END IF
nlx = (xmax-xmin)/lambda_x
nly = (ymax-ymin)/lambda_y
nlz = (zmax-zmin)/lambda_z
WRITE(*,'(A,F6.1,A,F6.1,A,F6.1,A,I0,A,I0,A)') ' TDGP3D domain size in GP length scales:', nlx, ' (x) x', nly, ' (y) x', &
   & nlz, ' (z) = ', NINT(nlx*nly*nlz), '   (kmax = ', kmax, ')'

OPEN(NEWUNIT=iunit, FILE=TRIM(filename), STATUS='REPLACE', ACTION='WRITE', IOSTAT=ios)
IF (ios == 0) THEN
   WRITE(iunit,'(A)')        '! domain for the GP node locations used by this run (written by REPORT_DOMAIN3D in tdgp3d_com.f90)'
   WRITE(iunit,'(A)')        '&tdgp_domain'
   WRITE(iunit,'(A,ES16.8)') 'xmin = ', xmin
   WRITE(iunit,'(A,ES16.8)') 'xmax = ', xmax
   WRITE(iunit,'(A,ES16.8)') 'ymin = ', ymin
   WRITE(iunit,'(A,ES16.8)') 'ymax = ', ymax
   WRITE(iunit,'(A,ES16.8)') 'zmin = ', zmin
   WRITE(iunit,'(A,ES16.8)') 'zmax = ', zmax
   WRITE(iunit,'(A,ES16.8)') 'lambda_increase = ', lambda_increase
   WRITE(iunit,'(A,ES16.8)') 'logrho_bg = ', logrho_bg
   WRITE(iunit,'(A)')        '/'
   CLOSE(iunit)
ELSE
   WRITE(*,*) 'REPORT_DOMAIN3D: could not write ', TRIM(filename)
END IF

END SUBROUTINE REPORT_DOMAIN3D

!-------------------------------------------------------------------------------
SUBROUTINE ALLOC_NODEOBJ(obj, ndat, nFree)
! Allocates a nodeobj's arrays: node arrays to the compile-time NODEMX ceiling,
! data arrays to ndat, and the GP mean at the earth nodes to nFree (= nstar).
! Only called on chain drivers and the PT coordinator.
!-------------------------------------------------------------------------------
IMPLICIT NONE
TYPE(nodeobj), INTENT(INOUT) :: obj
INTEGER,       INTENT(IN)    :: ndat, nFree

IF (ALLOCATED(obj%x))           DEALLOCATE(obj%x)
IF (ALLOCATED(obj%y))           DEALLOCATE(obj%y)
IF (ALLOCATED(obj%z))           DEALLOCATE(obj%z)
IF (ALLOCATED(obj%logrho))      DEALLOCATE(obj%logrho)
IF (ALLOCATED(obj%logrho_star)) DEALLOCATE(obj%logrho_star)
IF (ALLOCATED(obj%Dobs))        DEALLOCATE(obj%Dobs)
IF (ALLOCATED(obj%Dunc))        DEALLOCATE(obj%Dunc)
IF (ALLOCATED(obj%Dpred))       DEALLOCATE(obj%Dpred)

ALLOCATE(obj%x(NODEMX), obj%y(NODEMX), obj%z(NODEMX), obj%logrho(NODEMX))
ALLOCATE(obj%logrho_star(MAX(nFree,1)))
ALLOCATE(obj%Dobs(MAX(ndat,1)), obj%Dunc(MAX(ndat,1)), obj%Dpred(MAX(ndat,1)))

obj%x           = 0._RP
obj%y           = 0._RP
obj%z           = 0._RP
obj%logrho      = 0._RP
obj%logrho_star = 0._RP
obj%Dobs        = 0._RP
obj%Dunc        = 0._RP
obj%Dpred       = 0._RP
obj%k           = 0
obj%NDAT        = ndat
obj%ireject_bd  = 0
obj%iaccept_bd  = 0
obj%iaccept     = 0
obj%ireject     = 0
obj%nfwdfail    = 0
obj%beta        = 1._RP
obj%logL        = 0._RP
obj%chiSquaredMisfit = 0._RP
obj%acceptance  = 0._RP
obj%m0          = logrho_bg

END SUBROUTINE ALLOC_NODEOBJ

END MODULE TDGP3D_COM
