!===============================================================================
! MODULE TDGP_DATA_TYPE: shared precision/kind constants for the TDGP
! trans-dimensional Gaussian-process MT inversion codes (2D and 3D).
!===============================================================================
MODULE TDGP_DATA_TYPE
USE, INTRINSIC :: ISO_FORTRAN_ENV, ONLY: INT32, REAL64
IMPLICIT NONE

INTEGER, PARAMETER :: IB = INT32
INTEGER, PARAMETER :: RP = REAL64

REAL(RP), PARAMETER :: PI2 = 3.14159265358979323846264338327950288_RP

END MODULE TDGP_DATA_TYPE
