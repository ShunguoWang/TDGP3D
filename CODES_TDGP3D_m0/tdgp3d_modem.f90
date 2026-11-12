!===============================================================================
! MODULE TDGP3D_MODEM: the observed MT data, read from a ModEM data file.
!
! Supported ModEM data blocks:
!   Full_Impedance            components ZXX ZXY ZYX ZYY
!   Off_Diagonal_Impedance    components ZXY ZYX
!   Full_Vertical_Components  components TX TY (the tipper Tzx, Tzy)
! Each block has comment lines starting with '#', then six header lines starting with '>':
!   > block type
!   > sign convention: exp(-i\omega t) or exp(+i\omega t)
!   > units: [mV/km]/[nT], [V/m]/[T], [V/m]/[A/m] (= Ohm), or [] for the tipper
!   > orientation angle      (not used: the data must use the X/Y axes of the forward mesh)
!   > origin latitude, longitude          (not used)
!   > number of periods, number of sites  (not used: counted from the data rows)
! and then one data row per site, period and component:
!   Period(s) Code GG_Lat GG_Lon X(m) Y(m) Z(m) Component Real Imag Error
! X is north and Y east, in the coordinates of the forward mesh (MT3D_V6.dat);
! Z is not used (the forward mesh has a flat surface at z = 0). Every row gives two
! data, the real and the imaginary part, each with the standard error Error.
!
! The data are kept as they are in the file. The forward code's impedances (Ohm,
! exp(-i omega t)) are converted to the units and sign convention of each row by
! MODEM_VALUE, so the predicted data can be compared directly with the file.
!
! No MPI dependency: READ_MODEM_DATA returns an error code and message, and the
! caller aborts the run.
!===============================================================================
MODULE TDGP3D_MODEM
USE TDGP_DATA_TYPE
IMPLICIT NONE
PRIVATE

PUBLIC :: READ_MODEM_DATA, MODEM_VALUE

! component codes (row_comp):
INTEGER, PARAMETER, PUBLIC :: iZXX = 1, iZXY = 2, iZYX = 3, iZYY = 4, iTZX = 5, iTZY = 6

! --- MT sites ---
INTEGER,                        PUBLIC :: nsite = 0
CHARACTER(LEN=32), ALLOCATABLE, PUBLIC :: site_code(:)
REAL(RP),          ALLOCATABLE, PUBLIC :: site_x(:), site_y(:)   ! [m]

! --- periods, in the order they first appear in the file ---
INTEGER,               PUBLIC :: nper = 0
REAL(RP), ALLOCATABLE, PUBLIC :: period(:)                      ! [s]

! --- data rows: nrow rows give ndat = 2*nrow data (real, imaginary) ---
INTEGER,               PUBLIC :: nrow = 0, ndat = 0
INTEGER,  ALLOCATABLE, PUBLIC :: row_site(:), row_per(:), row_comp(:)
REAL(RP), ALLOCATABLE, PUBLIC :: row_conv(:)    ! impedance in the row's units = impedance in Ohm * row_conv
LOGICAL,  ALLOCATABLE, PUBLIC :: row_conj(:)    ! .TRUE. for exp(+i omega t) rows: complex conjugate of the forward code's value
REAL(RP), ALLOCATABLE, PUBLIC :: dobs(:), dunc(:)   ! (ndat): row r gives dobs(2r-1) = Real, dobs(2r) = Imag, both with dunc = Error

REAL(RP), PARAMETER :: MU0 = 4._RP*PI2*1.0e-7_RP

CONTAINS

!-------------------------------------------------------------------------------
SUBROUTINE READ_MODEM_DATA(filename, ierr, msg)
! Reads the sites, periods and data rows of a ModEM data file. ierr /= 0 (and a
! message in msg) if the file cannot be read or has something this code does not
! support.
!-------------------------------------------------------------------------------
IMPLICIT NONE
CHARACTER(LEN=*), INTENT(IN)  :: filename
INTEGER,          INTENT(OUT) :: ierr
CHARACTER(LEN=*), INTENT(OUT) :: msg

INTEGER  :: iu, ios, nmax, iline, nhead, btype, is, ip, icomp, j
LOGICAL  :: in_data, conj_blk
REAL(RP) :: conv_blk, per, rlat, rlon, xr, yr, zr, re, im, err
CHARACTER(LEN=1024) :: line
CHARACTER(LEN=512)  :: str
CHARACTER(LEN=64)   :: code, comp

ierr = 0
msg  = ''

OPEN(NEWUNIT=iu, FILE=TRIM(filename), STATUS='OLD', ACTION='READ', IOSTAT=ios)
IF (ios /= 0) THEN
   ierr = 1
   msg  = 'cannot open the ModEM data file '//TRIM(filename)
   RETURN
END IF

! --- pass 1: count the data rows ---
nmax = 0
DO
   READ(iu,'(A)',IOSTAT=ios) line
   IF (ios /= 0) EXIT
   line = ADJUSTL(line)
   IF (LEN_TRIM(line) == 0) CYCLE
   IF (line(1:1) == '#' .OR. line(1:1) == '>') CYCLE
   nmax = nmax + 1
END DO
IF (nmax == 0) THEN
   CLOSE(iu)
   ierr = 1
   msg  = 'no data rows in the ModEM data file '//TRIM(filename)
   RETURN
END IF

IF (ALLOCATED(site_code)) DEALLOCATE(site_code, site_x, site_y, period, row_site, row_per, row_comp, &
                                   & row_conv, row_conj, dobs, dunc)
ALLOCATE(site_code(nmax), site_x(nmax), site_y(nmax), period(nmax))
ALLOCATE(row_site(nmax), row_per(nmax), row_comp(nmax), row_conv(nmax), row_conj(nmax))
ALLOCATE(dobs(2*nmax), dunc(2*nmax))
nsite = 0
nper  = 0
nrow  = 0

! --- pass 2: parse ---
REWIND(iu)
nhead    = 0
in_data  = .FALSE.
btype    = 0
conj_blk = .FALSE.
conv_blk = 1._RP
iline    = 0
DO
   READ(iu,'(A)',IOSTAT=ios) line
   IF (ios /= 0) EXIT
   iline = iline + 1
   line  = ADJUSTL(line)
   IF (LEN_TRIM(line) == 0) CYCLE

   IF (line(1:1) == '#') THEN               ! comment line: a new block follows
      nhead   = 0
      in_data = .FALSE.
      CYCLE
   END IF

   IF (line(1:1) == '>') THEN               ! block header line
      IF (in_data) THEN                     ! a new block without '#' lines
         nhead   = 0
         in_data = .FALSE.
      END IF
      nhead = nhead + 1
      str   = LOWER(ADJUSTL(line(2:)))
      SELECT CASE (nhead)
      CASE (1)
         btype = 0
         IF (INDEX(str,'full_impedance') == 1)           btype = 1
         IF (INDEX(str,'off_diagonal_impedance') == 1)   btype = 2
         IF (INDEX(str,'full_vertical_components') == 1) btype = 3
         IF (btype == 0) THEN
            CALL FAIL('unsupported data type "'//TRIM(line(2:))//'" (supported: Full_Impedance, '// &
                    & 'Off_Diagonal_Impedance, Full_Vertical_Components)')
            RETURN
         END IF
      CASE (2)
         IF (INDEX(str,'+') > 0) THEN
            conj_blk = .TRUE.
         ELSE IF (INDEX(str,'-') > 0) THEN
            conj_blk = .FALSE.
         ELSE
            CALL FAIL('cannot read the sign convention "'//TRIM(line(2:))//'"')
            RETURN
         END IF
      CASE (3)
         IF (btype == 3) THEN
            conv_blk = 1._RP                       ! the tipper has no units
         ELSE IF (INDEX(str,'[mv/km]/[nt]') > 0) THEN
            conv_blk = 1._RP/(MU0*1000._RP)
         ELSE IF (INDEX(str,'[v/m]/[t]') > 0) THEN
            conv_blk = 1._RP/MU0
         ELSE IF (INDEX(str,'[v/m]/[a/m]') > 0 .OR. INDEX(str,'ohm') > 0) THEN
            conv_blk = 1._RP
         ELSE
            CALL FAIL('unknown impedance units "'//TRIM(line(2:))//'"')
            RETURN
         END IF
      CASE DEFAULT      ! orientation, origin, counts: not used
      END SELECT
      CYCLE
   END IF

   ! --- data row ---
   in_data = .TRUE.
   IF (btype == 0 .OR. nhead < 3) THEN
      CALL FAIL('data row before a complete block header')
      RETURN
   END IF
   READ(line,*,IOSTAT=ios) per, code, rlat, rlon, xr, yr, zr, comp, re, im, err
   IF (ios /= 0) THEN
      CALL FAIL('cannot read the data row')
      RETURN
   END IF

   icomp = 0
   SELECT CASE (TRIM(UPPER(comp)))
   CASE ('ZXX');        icomp = iZXX
   CASE ('ZXY');        icomp = iZXY
   CASE ('ZYX');        icomp = iZYX
   CASE ('ZYY');        icomp = iZYY
   CASE ('TX', 'TZX');  icomp = iTZX
   CASE ('TY', 'TZY');  icomp = iTZY
   END SELECT
   IF ( icomp == 0 .OR. (btype == 1 .AND. icomp > iZYY) .OR. &
      & (btype == 2 .AND. icomp /= iZXY .AND. icomp /= iZYX) .OR. (btype == 3 .AND. icomp < iTZX) ) THEN
      CALL FAIL('component '//TRIM(comp)//' does not belong to this data block')
      RETURN
   END IF
   IF (per <= 0._RP) THEN
      CALL FAIL('period <= 0')
      RETURN
   END IF
   IF (err <= 0._RP) THEN
      CALL FAIL('error <= 0 (give missing data a large error instead)')
      RETURN
   END IF

   ! site, by its code:
   is = 0
   DO j = 1,nsite
      IF (TRIM(site_code(j)) == TRIM(code)) THEN
         is = j
         EXIT
      END IF
   END DO
   IF (is == 0) THEN
      nsite = nsite + 1
      is    = nsite
      site_code(is) = code
      site_x(is)    = xr
      site_y(is)    = yr
   ELSE IF (ABS(xr-site_x(is)) > 0.01_RP .OR. ABS(yr-site_y(is)) > 0.01_RP) THEN
      CALL FAIL('site '//TRIM(code)//' appears with two different X, Y positions')
      RETURN
   END IF

   ! period (relative tolerance 1e-6):
   ip = 0
   DO j = 1,nper
      IF (ABS(per-period(j)) <= 1.0e-6_RP*period(j)) THEN
         ip = j
         EXIT
      END IF
   END DO
   IF (ip == 0) THEN
      nper = nper + 1
      ip   = nper
      period(ip) = per
   END IF

   nrow = nrow + 1
   row_site(nrow) = is
   row_per(nrow)  = ip
   row_comp(nrow) = icomp
   row_conv(nrow) = conv_blk
   row_conj(nrow) = conj_blk
   dobs(2*nrow-1) = re
   dobs(2*nrow)   = im
   dunc(2*nrow-1) = err
   dunc(2*nrow)   = err
END DO
CLOSE(iu)

ndat = 2*nrow

CONTAINS

   SUBROUTINE FAIL(text)
   CHARACTER(LEN=*), INTENT(IN) :: text
   WRITE(msg,'(A,A,I0,A,A)') TRIM(filename), ', line ', iline, ': ', TRIM(text)
   ierr = 1
   CLOSE(iu)
   END SUBROUTINE FAIL

END SUBROUTINE READ_MODEM_DATA

!-------------------------------------------------------------------------------
COMPLEX(RP) FUNCTION MODEM_VALUE(zt, r) RESULT(v)
! Value of data row r from the forward code's responses at the row's site and
! period: zt(1:6) = Zxx, Zxy, Zyx, Zyy [Ohm], Tzx, Tzy, all exp(-i omega t).
! Converted to the row's units and sign convention.
!-------------------------------------------------------------------------------
IMPLICIT NONE
COMPLEX(RP), INTENT(IN) :: zt(6)
INTEGER,     INTENT(IN) :: r

v = zt(row_comp(r))
IF (row_comp(r) <= iZYY) v = v*row_conv(r)
IF (row_conj(r)) v = CONJG(v)

END FUNCTION MODEM_VALUE

!-------------------------------------------------------------------------------
FUNCTION LOWER(s) RESULT(t)
IMPLICIT NONE
CHARACTER(LEN=*), INTENT(IN) :: s
CHARACTER(LEN=LEN(s)) :: t
INTEGER :: i
t = s
DO i = 1,LEN(s)
   IF (s(i:i) >= 'A' .AND. s(i:i) <= 'Z') t(i:i) = ACHAR(IACHAR(s(i:i))+32)
END DO
END FUNCTION LOWER

FUNCTION UPPER(s) RESULT(t)
IMPLICIT NONE
CHARACTER(LEN=*), INTENT(IN) :: s
CHARACTER(LEN=LEN(s)) :: t
INTEGER :: i
t = s
DO i = 1,LEN(s)
   IF (s(i:i) >= 'a' .AND. s(i:i) <= 'z') t(i:i) = ACHAR(IACHAR(s(i:i))-32)
END DO
END FUNCTION UPPER

END MODULE TDGP3D_MODEM
