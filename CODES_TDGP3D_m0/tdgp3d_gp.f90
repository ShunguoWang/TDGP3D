!===============================================================================
! MODULE TDGP3D_GP: Gaussian-process interpolation from a sparse, trans-dimensional
! set of (x, y, z, log10-resistivity) nodes onto the nodes of the 3D forward
! mesh -- the 3D counterpart of TDGP_GP (tdgp_gp.f90).
!
! Same squared-exponential kernel and GP conditional mean as the 2D code, with one
! length scale per axis and a constant prior mean m0:
!
!   K(x,x') = exp( -0.5 * [ ((x1-x2)/lambda_x)^2 + ((y1-y2)/lambda_y)^2 + ((z1-z2)/lambda_z)^2 ] )
!   Km      = K(x,x) + sigma_m^2 * I                  (node-node kernel, with nugget)
!   mustar  = m0 + K(x,x*)^T * Km^-1 * (g(x) - m0)    (GP conditional mean on the grid)
!
! m0 is the log10 background resistivity of the forward mesh, so the model goes
! back to the background away from the nodes (the 2D code has m0 = 0).
!
! Optionally (lambda_increase > 0 in tdgp.settings) all three length scales grow
! with the depth below the surface, l_c(z) = lambda_c + lambda_increase*MAX(z,0),
! c = x, y, z, so deeper nodes reach farther in every direction. The kernel is then
! the non-stationary squared exponential of Gibbs (1997) and Paciorek & Schervish
! (2004), which stays positive definite and is the kernel above when
! lambda_increase = 0 (GP_KERNEL3D_INC):
!   K(a,b) = PROD_c sqrt(2 l_c(a) l_c(b)/(l_c(a)^2 + l_c(b)^2)) * exp(-SUM_c (c_a-c_b)^2/(l_c(a)^2 + l_c(b)^2))
!
! Km^-1 is never formed: Km is Cholesky-factorized (DPOTRF) and solved (DPOTRS),
! as in the 2D code. K(x,x*) is not stored either -- the forward mesh has many more
! nodes than a MARE2DEM inversion grid has free regions -- each grid value is
! summed directly over the k GP nodes (OpenMP over the grid points). With constant
! length scales, terms with a squared scaled distance above QCUT (kernel value
! below exp(-40) ~ 4e-18) are skipped.
!
! No MPI or forward-code dependency: tested standalone by tdgp3d_test_gp.f90.
!===============================================================================
MODULE TDGP3D_GP
USE TDGP_DATA_TYPE
IMPLICIT NONE
PRIVATE

PUBLIC :: GP_KERNEL3D, GP_PREDICT_MEAN3D

REAL(RP), PARAMETER :: QCUT = 80._RP

CONTAINS

!-------------------------------------------------------------------------------
PURE REAL(RP) FUNCTION GP_KERNEL3D(x1,y1,z1,x2,y2,z2,lambda_x,lambda_y,lambda_z) RESULT(kval)
!-------------------------------------------------------------------------------
IMPLICIT NONE
REAL(RP), INTENT(IN) :: x1,y1,z1,x2,y2,z2,lambda_x,lambda_y,lambda_z
REAL(RP) :: dx,dy,dz

dx = (x1-x2)/lambda_x
dy = (y1-y2)/lambda_y
dz = (z1-z2)/lambda_z
kval = EXP( -0.5_RP*(dx*dx+dy*dy+dz*dz) )

END FUNCTION GP_KERNEL3D

!-------------------------------------------------------------------------------
PURE REAL(RP) FUNCTION GP_KERNEL3D_INC(x1,y1,z1,x2,y2,z2,lambda_x,lambda_y,lambda_z,inc) RESULT(kval)
! The kernel with length scales lambda + inc*MAX(z,0) at each point (see the
! module header). Equals GP_KERNEL3D for inc = 0.
!-------------------------------------------------------------------------------
IMPLICIT NONE
REAL(RP), INTENT(IN) :: x1,y1,z1,x2,y2,z2,lambda_x,lambda_y,lambda_z,inc
REAL(RP) :: l1(3), l2(3), s(3)

l1 = (/ lambda_x, lambda_y, lambda_z /) + inc*MAX(z1, 0._RP)
l2 = (/ lambda_x, lambda_y, lambda_z /) + inc*MAX(z2, 0._RP)
s  = l1*l1 + l2*l2
kval = SQRT(PRODUCT(2._RP*l1*l2/s)) * EXP( -(x1-x2)**2/s(1) - (y1-y2)**2/s(2) - (z1-z2)**2/s(3) )

END FUNCTION GP_KERNEL3D_INC

!-------------------------------------------------------------------------------
SUBROUTINE GP_PREDICT_MEAN3D(x,y,z,g,k,xs,ys,zs,nstar,lambda_x,lambda_y,lambda_z,sigma_m,mean0,mustar,ierr, &
                           & lambda_increase)
! GP conditional mean mustar(1:nstar) at the points (xs,ys,zs) given k
! (node, value) pairs and the prior mean mean0. ierr = 0 on success; ierr =
! LAPACK INFO code (nonzero) if Km is not positive definite (e.g. two nodes
! at almost the same position). Optional lambda_increase > 0: the length scales
! grow with depth (GP_KERNEL3D_INC).
!-------------------------------------------------------------------------------
IMPLICIT NONE
INTEGER(IB), INTENT(IN)  :: k, nstar
REAL(RP),    INTENT(IN)  :: x(k), y(k), z(k), g(k), xs(nstar), ys(nstar), zs(nstar)
REAL(RP),    INTENT(IN)  :: lambda_x, lambda_y, lambda_z, sigma_m, mean0
REAL(RP),    INTENT(OUT) :: mustar(nstar)
INTEGER(IB), INTENT(OUT) :: ierr
REAL(RP),    INTENT(IN), OPTIONAL :: lambda_increase

REAL(RP), ALLOCATABLE :: Km(:,:), alpha(:,:)
REAL(RP) :: ax, ay, az, dx, dy, dz, q, s, inc
INTEGER  :: i, j, info
EXTERNAL :: DPOTRF, DPOTRS

inc = 0._RP
IF (PRESENT(lambda_increase)) inc = lambda_increase

ierr = 0
IF (k < 1) THEN          ! no nodes: the prior mean everywhere
   mustar = mean0
   RETURN
END IF

ALLOCATE(Km(k,k), alpha(k,1))

DO j = 1,k
   DO i = 1,k
      IF (inc > 0._RP) THEN
         Km(i,j) = GP_KERNEL3D_INC(x(i),y(i),z(i),x(j),y(j),z(j),lambda_x,lambda_y,lambda_z,inc)
      ELSE
         Km(i,j) = GP_KERNEL3D(x(i),y(i),z(i),x(j),y(j),z(j),lambda_x,lambda_y,lambda_z)
      END IF
   END DO
END DO
DO i = 1,k
   Km(i,i) = Km(i,i) + sigma_m*sigma_m
END DO

alpha(1:k,1) = g(1:k) - mean0

! Cholesky factorize Km = L*L^T, then solve Km*alpha = g - mean0:
CALL DPOTRF('L', k, Km, k, info)
IF (info == 0) CALL DPOTRS('L', k, 1, Km, k, alpha, k, info)
IF (info /= 0) THEN
   ierr   = info
   mustar = mean0
   DEALLOCATE(Km, alpha)
   RETURN
END IF

! mustar = mean0 + K(x,x*)^T * alpha, one grid point at a time:
IF (inc > 0._RP) THEN
   !$OMP PARALLEL DO PRIVATE(i,s) SCHEDULE(STATIC)
   DO j = 1,nstar
      s = 0._RP
      DO i = 1,k
         s = s + alpha(i,1)*GP_KERNEL3D_INC(x(i),y(i),z(i),xs(j),ys(j),zs(j),lambda_x,lambda_y,lambda_z,inc)
      END DO
      mustar(j) = mean0 + s
   END DO
   !$OMP END PARALLEL DO
   DEALLOCATE(Km, alpha)
   RETURN
END IF
ax = 1._RP/lambda_x
ay = 1._RP/lambda_y
az = 1._RP/lambda_z
!$OMP PARALLEL DO PRIVATE(i,dx,dy,dz,q,s) SCHEDULE(STATIC)
DO j = 1,nstar
   s = 0._RP
   DO i = 1,k
      dx = (x(i)-xs(j))*ax
      dy = (y(i)-ys(j))*ay
      dz = (z(i)-zs(j))*az
      q  = dx*dx + dy*dy + dz*dz
      IF (q < QCUT) s = s + alpha(i,1)*EXP(-0.5_RP*q)
   END DO
   mustar(j) = mean0 + s
END DO
!$OMP END PARALLEL DO

DEALLOCATE(Km, alpha)

END SUBROUTINE GP_PREDICT_MEAN3D

END MODULE TDGP3D_GP
