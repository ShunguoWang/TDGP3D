!===============================================================================
! MODULE TDGP3D_LIB_INTERFACE
!
! The V6 3D MT forward code (CODES_MPIV6: mixed spatial-spectral GS-FFT solver
! with variable dz) as the forward operator of the 3D TDGP sampler -- the 3D
! counterpart of tdgp_lib_interface.f90 (MARE2DEM teams in the 2D code).
!
! Each PT chain has a team of ranks (team_comm, from TDGP_MPI_TOPOLOGY's
! SETUP_TEAMS). Every rank of a team holds the whole forward mesh and solves
! its share of the 2*NF polarization solves -- the x polarization of frequency K
! on team rank MOD(K-1, team size), the y polarization on MOD(K-1+NF, team size),
! so a team of up to 2*NF ranks is useful (with a team size that divides NF, both
! polarizations of a frequency are on one rank) -- with OpenMP threads inside
! each solve, as in the V6 driver, except that the team driver (team rank 0)
! solves as well instead of only collecting results. The scattered surface fields
! are then summed over the team (MPI_ALLREDUCE), and the rank of frequency K's x
! polarization computes the data of frequency K.
!
! Calling pattern (collective over the team, as in the 2D interface):
!   (1) every team rank calls tdgp3d_mt3d_init_team ONCE. It reads the forward
!       mesh from MT3D_V6.dat (V6's READ_PARA and MODEL) and the observed data
!       from the ModEM file (TDGP3D_MODEM). The data's periods replace the
!       freqmin/freqmax/nfreq lines of MT3D_V6.dat; its mesh_master lines are
!       not used (every frequency is solved on the MT3D_V6.dat mesh).
!   (2) every team rank calls tdgp3d_mt3d_forward once per forward solve. Only
!       the driver's input -- log10 resistivity at the earth nodes of the mesh
!       -- is used; the predicted data (driver only), the chi-squared misfit and
!       the convergence flag come back.
!
! Model: the forward code solves for the field scattered by the anomalous
! conductivity DR = sigma - sigma_background at every earth node (z >= 0), where
! sigma_background is the half-space conductivity of MT3D_V6.dat
! (background_conductivity) and air is 1e-8 S/m (V6's MODEL). The GP mean
! (TDGP3D_GP) has the background as its prior mean, so DR = 0 away from the
! GP nodes. With sample_m0 the background is each model's m0: the driver passes
! it with the model and SET_BACKGROUND applies it before the solve.
!
! SOLVE_FREQ repeats the V6 driver's SOLVE_FREQ_SURFACE_FIELDS
! (WFEM-3D_F_hmhz_fft_OMP_2018_V_02.f90), an internal procedure there that cannot
! be called from outside, with the same constants, without printing each
! iteration, and reporting whether the iteration converged within
! iteration_max_number. IMPEDANCE_AT evaluates V6's COMPACT_IMPEDANCE formulas at
! an MT site, with the scattered surface fields interpolated bilinearly between
! grid nodes (exact for a site on a node).
!===============================================================================
MODULE TDGP3D_LIB_INTERFACE
USE MPI
USE TDGP_DATA_TYPE
USE TDGP3D_MODEM
USE MODPAR1,         ONLY: NX, NY, NZ, KX, KY, DX, DY, JF, NF, FRE, SIGMA0, ITERMAX, NNA, NNB, kx_shift, ky_shift
USE MODPAR2,         ONLY: X, Y, Z
USE MODPAR3,         ONLY: VI, VF
USE MODPAR4,         ONLY: DZ, R0
USE MODPAR5,         ONLY: DR
USE PARAME1,         ONLY: WU, WE, WU0
USE GMTR1,           ONLY: G8
USE MT3D_A,          ONLY: Ax, Ay, Az
USE MT3D_E0_X,       ONLY: Exx, Exy, Exz
USE csem_jx_1D,      ONLY: EFxx, EFxy, EFxz, HFxx, HFxy, HFxz, EFxx1, EFxy1, EFxz1, HFxx1, HFxy1, HFxz1, &
                         & E0xx, E0xy, E0xz
USE MESHMASTER_CTRL, ONLY: MESHMASTER_ENABLE
USE INTERP_UTIL,     ONLY: BILINEAR2D_C
IMPLICIT NONE
PRIVATE

PUBLIC :: tdgp3d_mt3d_init_team, tdgp3d_mt3d_forward
PUBLIC :: tdgp3d_mt3d_get_grid, tdgp3d_mt3d_get_stations, tdgp3d_mt3d_get_mesh, tdgp3d_mt3d_get_data

CHARACTER(LEN=*), PARAMETER :: MESHFILE = 'MT3D_V6.dat'   ! the file name V6's READ_PARA reads

INTEGER, SAVE :: fcomm = -1, frank = 0, fsize = 1
INTEGER, SAVE :: nze = 0, nstar = 0          ! earth nodes per (x,y) column (N = JF..NZ), and in total
REAL(RP),    ALLOCATABLE, SAVE :: rbuf(:)    ! (nstar) log10 rho at the earth nodes, broadcast from the driver
REAL(RP),    ALLOCATABLE, SAVE :: dloc(:), dall(:)   ! (ndat) this rank's / the team's predicted data
COMPLEX(RP), ALLOCATABLE, TARGET, SAVE :: sfield(:,:,:,:)   ! (NX,NY,10,NF) scattered surface fields of every frequency
COMPLEX(RP), POINTER, SAVE :: e1x(:,:), e1y(:,:), h1x(:,:), h1y(:,:), h1z(:,:)   ! scattered surface fields,
COMPLEX(RP), POINTER, SAVE :: e2x(:,:), e2y(:,:), h2x(:,:), h2y(:,:), h2z(:,:)   ! x and y polarization (NX,NY) of
                                                                                   ! the frequency set by POINT_FIELDS

! Constants of the V6 driver, with its types and expressions, so SOLVE_FREQ gives the same numbers:
COMPLEX(KIND=8), SAVE :: C0, PAI
REAL(KIND=8),    SAVE :: FMU, EPS

EXTERNAL :: READ_PARA, SZ_ALLOCATE, MODEL, GS_FFT_WAVE_NUMBER, GMATR1, GMATR2, csem_x_1d
EXTERNAL :: fxfyfz0_x, fxfyfz_x, fxfyfz0_y, fxfyfz_y, MT3DA, GaussFFT_2D

CONTAINS

!-------------------------------------------------------------------------------
SUBROUTINE tdgp3d_mt3d_init_team(team_comm_in, dataFileName, nstar_out, ndat_out, report)
! One-time setup of this rank's forward code, collective over the team. Every
! rank reads the data file and MT3D_V6.dat itself (both are small). report =
! .true. prints a summary (pass .true. on one rank of the whole job only).
!
! Outputs (every rank): nstar_out = number of earth nodes of the forward mesh
! (the length of the GP mean vector), ndat_out = number of data.
!
! Note: V6's MODEL prints the vertical grid and writes model_data.dat (its own
! example model, not used here) when it runs on MPI_COMM_WORLD rank 0, i.e. only
! in tdgp3d_single; in tdgp3d_mpi, rank 0 is the PT coordinator.
!-------------------------------------------------------------------------------
IMPLICIT NONE
INTEGER,          INTENT(IN)  :: team_comm_in
CHARACTER(LEN=*), INTENT(IN)  :: dataFileName
INTEGER,          INTENT(OUT) :: nstar_out, ndat_out
LOGICAL,          INTENT(IN)  :: report

INTEGER  :: ierr, ios, K, is, noff
LOGICAL  :: ex
REAL(RP) :: dsite, fmax, skin
CHARACTER(LEN=512) :: msg

fcomm = team_comm_in
CALL MPI_COMM_RANK(fcomm, frank, ierr)
CALL MPI_COMM_SIZE(fcomm, fsize, ierr)

! --- observed data (ModEM) ---
CALL READ_MODEM_DATA(TRIM(dataFileName), ios, msg)
IF (ios /= 0) CALL FATAL(TRIM(msg))

! --- forward mesh (V6) ---
INQUIRE(FILE=MESHFILE, EXIST=ex)
IF (.NOT. ex) CALL FATAL('the forward mesh file '//MESHFILE//' is not in the run folder')
CALL READ_PARA()

NF = nper                         ! the data's periods replace the frequencies of MT3D_V6.dat
CALL SZ_ALLOCATE()
DO K = 1,NF
   FRE(K) = 1._RP/period(K)
END DO
CALL MODEL()                      ! X, Y, Z, DZ, R0 and the dz-derivative weights (its example DR is replaced below)
DR = 0._RP
CALL GS_FFT_WAVE_NUMBER()
CALL GMATR1()

C0  = (0.,1.)
PAI = ATAN2(0.,-1.)
FMU = 4.*PAI*1.e-7
EPS = 8.854*1.e-12

! --- MT sites must lie inside the mesh, on its surface ---
noff = 0
DO is = 1,nsite
   IF (site_x(is) < X(1) .OR. site_x(is) > X(NX) .OR. site_y(is) < Y(1) .OR. site_y(is) > Y(NY)) THEN
      WRITE(msg,'(A,A,A,2F12.1,A)') 'site ', TRIM(site_code(is)), ' at x, y =', site_x(is), site_y(is), &
                                  & ' m is outside the forward mesh of '//MESHFILE
      CALL FATAL(TRIM(msg))
   END IF
   dsite = MINVAL(ABS(X(1:NX)-site_x(is))) + MINVAL(ABS(Y(1:NY)-site_y(is)))
   IF (dsite > 1.e-3_RP) noff = noff + 1
END DO

nze   = NZ - JF + 1
nstar = nze*NX*NY
ALLOCATE(rbuf(nstar), dloc(ndat), dall(ndat))
ALLOCATE(sfield(NX,NY,10,NF))

nstar_out = nstar
ndat_out  = ndat

IF (report) THEN
   fmax = MAXVAL(FRE(1:NF))
   skin = 503.3_RP*SQRT(1._RP/(SIGMA0*fmax))
   WRITE(*,'(A)')                    ' TDGP3D forward: V6 3D MT code, mesh from '//MESHFILE
   WRITE(*,'(A,I0,A,I0,A,I0,A,2F9.1,A)') '   NX x NY x NZ = ', NX, ' x ', NY, ' x ', NZ, ',  dx, dy =', DX, DY, ' m'
   WRITE(*,'(A,2F12.1,A,2F12.1,A,F12.1,A)') '   x =', X(1), X(NX), ' m,  y =', Y(1), Y(NY), ' m,  bottom z =', Z(NZ), ' m'
   WRITE(*,'(A,ES11.4,A,F10.2,A,I0)') '   background conductivity', SIGMA0, ' S/m (', 1._RP/SIGMA0, &
                                     & ' ohm-m),  iteration_max_number = ', ITERMAX
   WRITE(*,'(A,I0,A,I0,A,I0,A,I0,A)') ' TDGP3D data: ', nsite, ' sites, ', nper, ' periods, ', nrow, ' data rows = ', &
                                     & ndat, ' data (real and imaginary parts)'
   WRITE(*,'(A,*(1X,ES10.3))')       '   frequencies [Hz]:', FRE(1:NF)
   WRITE(*,'(A,I0,A,I0,A)') '   the frequencies of '//MESHFILE//' are not used; each chain solves them on ', fsize, &
                     & ' rank(s) of the team: x polarization of frequency K on team rank MOD(K-1, n), '// &
                     & 'y polarization on MOD(K-1+', NF, ', n), n = ranks per team'
   IF (fsize > 2*NF) WRITE(*,'(A)') ' WARNING: more ranks per team than twice the number of periods -- '// &
                                   & 'some ranks have nothing to solve'
   IF (MESHMASTER_ENABLE) WRITE(*,'(A)') '   the mesh_master lines of '//MESHFILE//' are not used: '// &
                                        & 'every frequency is solved on the '//MESHFILE//' mesh'
   IF (noff > 0) WRITE(*,'(A,I0,A)') '   ', noff, ' site(s) not on a grid node: their surface fields are interpolated bilinearly'
   WRITE(*,'(A,F10.1,A,F6.1,A)') '   skin depth in the background at the highest frequency:', skin, ' m = ', &
                                & skin/MAX(DX,DY), ' horizontal cells'
END IF

END SUBROUTINE tdgp3d_mt3d_init_team

!-------------------------------------------------------------------------------
SUBROUTINE tdgp3d_mt3d_forward(logrho_star, chiSquaredMisfit, Dpred, converged, logrho_bg_model)
! One forward solve, collective over the team.
!
! logrho_star      : (nstar) log10 resistivity at the earth nodes, in the order of
!                    tdgp3d_mt3d_get_grid -- used on the team driver only (other
!                    ranks may pass an empty array).
! logrho_bg_model  : optional, team driver only: log10 resistivity of this model's
!                    half-space background (its m0, sampled if sample_m0). Absent:
!                    the background stays as it is (MT3D_V6.dat's, unless changed).
! chiSquaredMisfit : sum(((dobs - Dpred)/dunc)**2), on every rank.
! Dpred            : (ndat) predicted data, in the order of the data rows (real,
!                    imaginary) and in their units and sign convention. Filled if
!                    its size is ndat (the driver passes obj%Dpred, other ranks may
!                    pass an empty array).
! converged        : .false. if the iteration of any frequency did not converge
!                    within iteration_max_number (MT3D_V6.dat), on every rank.
!-------------------------------------------------------------------------------
IMPLICIT NONE
REAL(RP), INTENT(IN)  :: logrho_star(:)
REAL(RP), INTENT(OUT) :: chiSquaredMisfit
REAL(RP), INTENT(OUT) :: Dpred(:)
LOGICAL,  INTENT(OUT) :: converged
REAL(RP), INTENT(IN), OPTIONAL :: logrho_bg_model

INTEGER  :: ierr, K, I, J, N, is
LOGICAL  :: conv_loc, conv_k, dox, doy
REAL(RP) :: bg(1)

IF (frank == 0) rbuf = logrho_star(1:nstar)
CALL MPI_BCAST(rbuf, nstar, MPI_DOUBLE_PRECISION, 0, fcomm, ierr)

! this model's half-space background, from the driver (every rank sets it before DR is formed):
bg(1) = LOG10(1._RP/SIGMA0)
IF (frank == 0 .AND. PRESENT(logrho_bg_model)) bg(1) = logrho_bg_model
CALL MPI_BCAST(bg, 1, MPI_DOUBLE_PRECISION, 0, fcomm, ierr)
CALL SET_BACKGROUND(bg(1))

! anomalous conductivity at every earth node (air: 0):
DR = 0._RP
DO J = 1,NY
   DO I = 1,NX
      is = ((J-1)*NX + (I-1))*nze
      DO N = JF,NZ
         DR(N,I,J) = 10._RP**(-rbuf(is+N-JF+1)) - R0(N)
      END DO
   END DO
END DO

dloc     = 0._RP
sfield   = (0._RP, 0._RP)
conv_loc = .TRUE.
DO K = 1,NF
   dox = MOD(K-1,    fsize) == frank           ! x polarization of frequency K on this rank
   doy = MOD(K-1+NF, fsize) == frank           ! y polarization
   IF (.NOT. (dox .OR. doy)) CYCLE
   CALL csem_x_1d(K)
   CALL POINT_FIELDS(K)
   CALL SOLVE_FREQ(K, dox, doy, conv_k)
   conv_loc = conv_loc .AND. conv_k
END DO

! the two polarizations of a frequency may have been solved on different ranks:
CALL MPI_ALLREDUCE(MPI_IN_PLACE, sfield, SIZE(sfield), MPI_DOUBLE_COMPLEX, MPI_SUM, fcomm, ierr)
DO K = 1,NF
   IF (MOD(K-1,fsize) /= frank) CYCLE
   CALL csem_x_1d(K)                           ! frequency K's primary fields (a later solve may have replaced them)
   CALL POINT_FIELDS(K)
   CALL SITE_RESPONSES(K)
END DO

CALL MPI_ALLREDUCE(dloc, dall, ndat, MPI_DOUBLE_PRECISION, MPI_SUM, fcomm, ierr)
CALL MPI_ALLREDUCE(conv_loc, converged, 1, MPI_LOGICAL, MPI_LAND, fcomm, ierr)

chiSquaredMisfit = SUM(((dobs(1:ndat)-dall)/dunc(1:ndat))**2)
IF (SIZE(Dpred) == ndat) Dpred = dall

END SUBROUTINE tdgp3d_mt3d_forward

!-------------------------------------------------------------------------------
SUBROUTINE SET_BACKGROUND(logrho_bg_new)
! Sets the forward code's half-space background to log10 rho = logrho_bg_new:
! SIGMA0 (background_conductivity of MT3D_V6.dat) and the earth part of the
! background profile, R0(JF:NZ), as MODEL does at startup. Everything that
! depends on the background reads R0 during each solve: DR here, the background
! fields (csem_x_1d), the depth operator (MT3DA) and the contraction (fxfyfz_x/y).
! Does nothing if the background is unchanged, so a fixed background gives
! exactly the numbers of the code without this routine.
!-------------------------------------------------------------------------------
IMPLICIT NONE
REAL(RP), INTENT(IN) :: logrho_bg_new
REAL(RP) :: s

s = 10._RP**(-logrho_bg_new)
IF (ABS(s - SIGMA0) <= 1.e-12_RP*SIGMA0) RETURN
SIGMA0    = s
R0(JF:NZ) = s

END SUBROUTINE SET_BACKGROUND

!-------------------------------------------------------------------------------
SUBROUTINE SOLVE_FREQ(K, dox, doy, conv)
! The polarizations of frequency K given by dox (x) and doy (y) on the current
! model (DR): the scattered surface fields e1x..h1z (x) and e2x..h2z (y).
! Repeats SOLVE_FREQ_SURFACE_FIELDS of the V6 driver; conv = .false. if a
! polarization's iteration did not converge. Needs csem_x_1d(K) (the primary
! fields) and POINT_FIELDS(K) first.
!-------------------------------------------------------------------------------
IMPLICIT NONE
INTEGER, INTENT(IN)  :: K
LOGICAL, INTENT(IN)  :: dox, doy
LOGICAL, INTENT(OUT) :: conv
LOGICAL :: convx, convy
INTEGER :: I, J

convx = .TRUE.
convy = .TRUE.

!-------------------------- X polarization ---------------------
IF (dox) THEN
   CALL GMATR2()
   E0xx(:,:,:) = EFxx(:,:,:)
   E0xy(:,:,:) = EFxy(:,:,:)
   E0xz(:,:,:) = EFxz(:,:,:)
   WU  = C0*2.*PAI*FRE(K)*FMU
   We  = C0*2.*PAI*FRE(K)*EPS
   Wu0 = 2.*PAI*FRE(K)*FMU
   Exx(:,:,:) = 0.
   Exy(:,:,:) = 0.
   Exz(:,:,:) = 0.
   CALL fxfyfz0_x()
   CALL ITERATE(1, convx)
   DO I = 1,NX
      DO J = 1,NY
         e1x(I,J) = Exx(JF,I,J)
         e1y(I,J) = Exy(JF,I,J)
      END DO
   END DO
   CALL SURFACE_H(h1x, h1y, h1z)
END IF

!-------------------------- Y polarization ---------------------
IF (doy) THEN
   CALL GMATR2()
   E0xx(:,:,:) = EFxx1(:,:,:)
   E0xy(:,:,:) = EFxy1(:,:,:)
   E0xz(:,:,:) = EFxz1(:,:,:)
   WU  = C0*2.*PAI*FRE(K)*FMU
   We  = C0*2.*PAI*FRE(K)*EPS
   Wu0 = 2.*PAI*FRE(K)*FMU
   Exx(:,:,:) = 0.
   Exy(:,:,:) = 0.
   Exz(:,:,:) = 0.
   CALL fxfyfz0_y()
   CALL ITERATE(2, convy)
   DO I = 1,NX
      DO J = 1,NY
         e2x(I,J) = Exx(JF,I,J)
         e2y(I,J) = Exy(JF,I,J)
      END DO
   END DO
   CALL SURFACE_H(h2x, h2y, h2z)
END IF

conv = convx .AND. convy

END SUBROUTINE SOLVE_FREQ

!-------------------------------------------------------------------------------
SUBROUTINE POINT_FIELDS(K)
! Points e1x..h2z at frequency K's scattered surface fields in sfield.
!-------------------------------------------------------------------------------
IMPLICIT NONE
INTEGER, INTENT(IN) :: K

e1x => sfield(:,:,1,K)
e1y => sfield(:,:,2,K)
h1x => sfield(:,:,3,K)
h1y => sfield(:,:,4,K)
h1z => sfield(:,:,5,K)
e2x => sfield(:,:,6,K)
e2y => sfield(:,:,7,K)
h2x => sfield(:,:,8,K)
h2y => sfield(:,:,9,K)
h2z => sfield(:,:,10,K)

END SUBROUTINE POINT_FIELDS

!-------------------------------------------------------------------------------
SUBROUTINE ITERATE(ipol, conv)
! The V6 iteration for one polarization (ipol = 1: x, 2: y) with its stopping
! rules: |sum of Ex| < 1e-20, or a relative change of it below 1e-4. conv =
! .false. if neither happened within ITERMAX iterations.
!-------------------------------------------------------------------------------
IMPLICIT NONE
INTEGER, INTENT(IN)  :: ipol
LOGICAL, INTENT(OUT) :: conv
INTEGER      :: ITER
REAL(KIND=8) :: s, s1

conv = .FALSE.
s1   = 0.
DO ITER = 1,ITERMAX
   CALL MT3DA()
   IF (ipol == 1) THEN
      CALL fxfyfz_x(s)
   ELSE
      CALL fxfyfz_y(s)
   END IF
   IF (s .LT. 1.E-20) THEN
      conv = .TRUE.
      EXIT
   END IF
   IF (ABS(s1-s)/s .LT. 0.0001) THEN
      conv = .TRUE.
      EXIT
   END IF
   s1 = s
END DO

END SUBROUTINE ITERATE

!-------------------------------------------------------------------------------
SUBROUTINE SURFACE_H(hx, hy, hz)
! Scattered magnetic field at the surface node JF from the potentials of the
! last MT3DA call -- the three "Surface H" blocks of the V6 driver.
!-------------------------------------------------------------------------------
IMPLICIT NONE
COMPLEX(RP), INTENT(OUT) :: hx(NX,NY), hy(NX,NY), hz(NX,NY)
INTEGER :: I, J, KL, KP

DO KL=1,NNA
DO KP=1,NNB
DO I=1,KX
DO J=1,KY
      VI(KL,KP,I,J)=C0*ky_shift(KP,J)*Az(KL,KP,I,J,1) &
        -(G8(1,1)*Ay(KL,KP,I,J,1)+G8(1,2)*Ay(KL,KP,I,J,2)+G8(1,3)*Ay(KL,KP,I,J,3))/DZ((JF-1)/2+1)
END DO
END DO
END DO
END DO
CALL GaussFFT_2D(NX,NY,X(1),Y(1),DX,DY,NNA,X,Y,VF,VI,-1)
hx(:,:)=VF(:,:)/WU

DO KL=1,NNA
DO KP=1,NNB
DO I=1,KX
DO J=1,KY
      VI(KL,KP,I,J)= &
        (G8(1,1)*Ax(KL,KP,I,J,1)+G8(1,2)*Ax(KL,KP,I,J,2)+G8(1,3)*Ax(KL,KP,I,J,3))/DZ((JF-1)/2+1) &
        -C0*kx_shift(KL,I)*Az(KL,KP,I,J,1)
END DO
END DO
END DO
END DO
CALL GaussFFT_2D(NX,NY,X(1),Y(1),DX,DY,NNA,X,Y,VF,VI,-1)
hy(:,:)=VF(:,:)/WU

DO KL=1,NNA
DO KP=1,NNB
DO I=1,KX
DO J=1,KY
      VI(KL,KP,I,J)=C0*kx_shift(KL,I)*Ay(KL,KP,I,J,1) &
                    -C0*ky_shift(KP,J)*Ax(KL,KP,I,J,1)
END DO
END DO
END DO
END DO
CALL GaussFFT_2D(NX,NY,X(1),Y(1),DX,DY,NNA,X,Y,VF,VI,-1)
hz(:,:)=VF(:,:)/WU

END SUBROUTINE SURFACE_H

!-------------------------------------------------------------------------------
SUBROUTINE SITE_RESPONSES(K)
! Predicted data of every data row at period K, into dloc.
!-------------------------------------------------------------------------------
IMPLICIT NONE
INTEGER, INTENT(IN) :: K
LOGICAL     :: need(nsite)
COMPLEX(RP) :: zt(6,nsite), v
INTEGER     :: r, is

need = .FALSE.
DO r = 1,nrow
   IF (row_per(r) == K) need(row_site(r)) = .TRUE.
END DO
DO is = 1,nsite
   IF (need(is)) CALL IMPEDANCE_AT(site_x(is), site_y(is), zt(:,is))
END DO
DO r = 1,nrow
   IF (row_per(r) /= K) CYCLE
   v = MODEM_VALUE(zt(:,row_site(r)), r)
   dloc(2*r-1) = REAL(v, RP)
   dloc(2*r)   = AIMAG(v)
END DO

END SUBROUTINE SITE_RESPONSES

!-------------------------------------------------------------------------------
SUBROUTINE IMPEDANCE_AT(xq, yq, zt)
! Impedance and tipper at the surface point (xq, yq): V6's COMPACT_IMPEDANCE
! formulas, with total field = primary (1D, the same at every node) + scattered
! (interpolated bilinearly between grid nodes). zt = Zxx, Zxy, Zyx, Zyy [Ohm],
! Tzx, Tzy, time dependence exp(-i omega t).
!-------------------------------------------------------------------------------
IMPLICIT NONE
REAL(RP),    INTENT(IN)  :: xq, yq
COMPLEX(RP), INTENT(OUT) :: zt(6)
COMPLEX(RP) :: ex1, ey1, hx1, hy1, hz1, ex2, ey2, hx2, hy2, hz2, den, s

CALL BILINEAR2D_C(X,NX,Y,NY,e1x,xq,yq,s);  ex1 = EFxx(JF,1,1)  + s
CALL BILINEAR2D_C(X,NX,Y,NY,e1y,xq,yq,s);  ey1 = EFxy(JF,1,1)  + s
CALL BILINEAR2D_C(X,NX,Y,NY,h1x,xq,yq,s);  hx1 = HFxx(JF,1,1)  + s
CALL BILINEAR2D_C(X,NX,Y,NY,h1y,xq,yq,s);  hy1 = HFxy(JF,1,1)  + s
CALL BILINEAR2D_C(X,NX,Y,NY,h1z,xq,yq,s);  hz1 = HFxz(JF,1,1)  + s
CALL BILINEAR2D_C(X,NX,Y,NY,e2x,xq,yq,s);  ex2 = EFxx1(JF,1,1) + s
CALL BILINEAR2D_C(X,NX,Y,NY,e2y,xq,yq,s);  ey2 = EFxy1(JF,1,1) + s
CALL BILINEAR2D_C(X,NX,Y,NY,h2x,xq,yq,s);  hx2 = HFxx1(JF,1,1) + s
CALL BILINEAR2D_C(X,NX,Y,NY,h2y,xq,yq,s);  hy2 = HFxy1(JF,1,1) + s
CALL BILINEAR2D_C(X,NX,Y,NY,h2z,xq,yq,s);  hz2 = HFxz1(JF,1,1) + s

den = hx1*hy2 - hx2*hy1
zt(iZXX) = (ex1*hy2 - ex2*hy1)/den
zt(iZXY) = (ex2*hx1 - ex1*hx2)/den
zt(iZYX) = (ey1*hy2 - ey2*hy1)/den
zt(iZYY) = (ey2*hx1 - ey1*hx2)/den
zt(iTZX) = (hz1*hy2 - hz2*hy1)/den
zt(iTZY) = (hz2*hx1 - hz1*hx2)/den

END SUBROUTINE IMPEDANCE_AT

!-------------------------------------------------------------------------------
SUBROUTINE tdgp3d_mt3d_get_grid(nstar_in, xs, ys, zs)
! Positions of the earth nodes of the forward mesh (z >= 0) -- where the GP mean
! is evaluated -- in the order tdgp3d_mt3d_forward expects logrho_star in.
!-------------------------------------------------------------------------------
IMPLICIT NONE
INTEGER,  INTENT(IN)  :: nstar_in
REAL(RP), INTENT(OUT) :: xs(nstar_in), ys(nstar_in), zs(nstar_in)
INTEGER :: I, J, N, is, ierr

IF (nstar_in /= nstar) THEN
   WRITE(*,*) 'tdgp3d_mt3d_get_grid: nstar_in does not match the forward mesh -- ', nstar_in, nstar
   CALL MPI_ABORT(MPI_COMM_WORLD, 1, ierr)
END IF
DO J = 1,NY
   DO I = 1,NX
      is = ((J-1)*NX + (I-1))*nze
      DO N = JF,NZ
         xs(is+N-JF+1) = X(I)
         ys(is+N-JF+1) = Y(J)
         zs(is+N-JF+1) = Z(N)
      END DO
   END DO
END DO

END SUBROUTINE tdgp3d_mt3d_get_grid

!-------------------------------------------------------------------------------
SUBROUTINE tdgp3d_mt3d_get_stations(nsta, xsta, ysta, nfreq, freq)
! MT site positions [m] and frequencies [Hz] of the data, for SET_DOMAIN3D.
!-------------------------------------------------------------------------------
IMPLICIT NONE
INTEGER,               INTENT(OUT) :: nsta, nfreq
REAL(RP), ALLOCATABLE, INTENT(OUT) :: xsta(:), ysta(:), freq(:)

nsta  = nsite
nfreq = nper
ALLOCATE(xsta(nsta), ysta(nsta), freq(nfreq))
xsta = site_x(1:nsite)
ysta = site_y(1:nsite)
freq = 1._RP/period(1:nper)

END SUBROUTINE tdgp3d_mt3d_get_stations

!-------------------------------------------------------------------------------
SUBROUTINE tdgp3d_mt3d_get_mesh(xlo, xhi, ylo, yhi, zbot, logrho_background)
! Extent of the forward mesh [m] and log10 of its background resistivity.
!-------------------------------------------------------------------------------
IMPLICIT NONE
REAL(RP), INTENT(OUT) :: xlo, xhi, ylo, yhi, zbot, logrho_background

xlo  = X(1)
xhi  = X(NX)
ylo  = Y(1)
yhi  = Y(NY)
zbot = Z(NZ)
logrho_background = LOG10(1._RP/SIGMA0)

END SUBROUTINE tdgp3d_mt3d_get_mesh

!-------------------------------------------------------------------------------
SUBROUTINE tdgp3d_mt3d_get_data(ndat_in, d, sd)
! The observed data and their standard errors, in the order of the data rows.
!-------------------------------------------------------------------------------
IMPLICIT NONE
INTEGER,  INTENT(IN)  :: ndat_in
REAL(RP), INTENT(OUT) :: d(ndat_in), sd(ndat_in)

d  = dobs(1:ndat_in)
sd = dunc(1:ndat_in)

END SUBROUTINE tdgp3d_mt3d_get_data

!-------------------------------------------------------------------------------
SUBROUTINE FATAL(text)
! Stops the whole job: an error here happens on every rank of every team, but the
! PT coordinator never calls this module and would otherwise wait forever.
!-------------------------------------------------------------------------------
IMPLICIT NONE
CHARACTER(LEN=*), INTENT(IN) :: text
INTEGER :: ierr

IF (frank == 0) WRITE(*,'(A)') ' TDGP3D: FATAL -- '//TRIM(text)
CALL MPI_ABORT(MPI_COMM_WORLD, 1, ierr)

END SUBROUTINE FATAL

END MODULE TDGP3D_LIB_INTERFACE
