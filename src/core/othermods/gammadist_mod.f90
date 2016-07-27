MODULE gammadist_mod
!module to create random deviates according to a Maxwell-Boltzmann Distribution
! (which is a Gamma Distribution, with shape 3/2 and scale kT)
! called directly by the main program (JT 2016)

USE GLOBAL

IMPLICIT NONE

  PRIVATE
  PUBLIC :: random_gamma
  PRIVATE:: random_gamma1, random_gamma2
  
 CONTAINS
!-----------------------------------------------------------------------------!   
FUNCTION random_gamma(s, b, first) 

! Adapted from Fortran 77 code from the book:
!     Dagpunar, J. 'Principles of random variate generation'
!     Clarendon Press, Oxford, 1988.   ISBN 0-19-852202-9

!     N.B. This version is in `double precision' and includes scaling

!     FUNCTION GENERATES A RANDOM GAMMA VARIATE.
!     CALLS EITHER random_gamma1 (S > 1.0)
!     OR random_exponential (S = 1.0)
!     OR random_gamma2 (S < 1.0).

!     S = SHAPE PARAMETER OF DISTRIBUTION (0 < REAL).
!     B = Scale parameter

IMPLICIT NONE
!INTEGER, PARAMETER  :: dp = SELECTED_REAL_KIND(12, 60)

REAL (num), INTENT(IN)  :: s, b
LOGICAL, INTENT(IN)    :: first
REAL (num)              :: random_gamma


IF (s <= zero) THEN
  WRITE(*, *) 'SHAPE PARAMETER VALUE MUST BE POSITIVE'
  STOP
END IF

IF (s >= one) THEN
  random_gamma = random_gamma1(s, first)
ELSE IF (s < one) THEN
  random_gamma = random_gamma2(s, first)
END IF

! Now scale the random variable
random_gamma = b * random_gamma

END FUNCTION random_gamma
!----------------------------------------------------------------------------!
FUNCTION random_gamma1(s, first)

! Adapted from Fortran 77 code from the book:
!     Dagpunar, J. 'Principles of random variate generation'
!     Clarendon Press, Oxford, 1988.   ISBN 0-19-852202-9

! FUNCTION GENERATES A RANDOM VARIATE IN [0,INFINITY) FROM
! A GAMMA DISTRIBUTION WITH DENSITY PROPORTIONAL TO GAMMA**(S-1)*EXP(-GAMMA),
! BASED UPON BEST'S T DISTRIBUTION METHOD

!     S = SHAPE PARAMETER OF DISTRIBUTION
!          (1.0 < REAL)

REAL (num), INTENT(IN)  :: s
LOGICAL, INTENT(IN)    :: first
REAL (num)              :: random_gamma1

!     Local variables
REAL (num)             :: d, r, g, f, x
REAL (num), SAVE       :: b, h
REAL (num), PARAMETER  :: sixty4 = 64.0_num, three = 3.0_num, pt75 = 0.75_num,  &
                         two = 2.0_num, half = 0.5_num

IF (s <= one) THEN
  WRITE(*, *) 'IMPERMISSIBLE SHAPE PARAMETER VALUE'
  STOP
END IF

IF (first) THEN                        ! Initialization, if necessary
  b = s - one
  h = SQRT(three*s - pt75)
END IF

DO
  CALL RANDOM_NUMBER(r)
  g = r - r*r
  IF (g <= zero) CYCLE
  f = (r - half)*h/SQRT(g)
  x = b + f
  IF (x <= zero) CYCLE
  CALL RANDOM_NUMBER(r)
  d = sixty4*g*(r*g)**2
  IF (d <= zero) EXIT
  IF (d*x < x - two*f*f) EXIT
  IF (LOG(d) < two*(b*LOG(x/b) - f)) EXIT
END DO
random_gamma1 = x

END FUNCTION random_gamma1
!----------------------------------------------------------------------------!
FUNCTION random_gamma2(s, first) 
! Adapted from Fortran 77 code from the book:
!     Dagpunar, J. 'Principles of random variate generation'
!     Clarendon Press, Oxford, 1988.   ISBN 0-19-852202-9

! FUNCTION GENERATES A RANDOM VARIATE IN [0,INFINITY) FROM
! A GAMMA DISTRIBUTION WITH DENSITY PROPORTIONAL TO
! GAMMA2**(S-1) * EXP(-GAMMA2),
! USING A SWITCHING METHOD.

!    S = SHAPE PARAMETER OF DISTRIBUTION
!          (REAL < 1.0)

REAL (num), INTENT(IN)  :: s
LOGICAL, INTENT(IN)     :: first
REAL (num)              :: random_gamma2

!     Local variables
REAL (num)             :: r, x, w
REAL (num), SAVE       :: a, p, c, uf, vr, d
REAL (num), SAVE       :: Vvsmall

vVsmall = EPSILON(1.0_num)

IF (s <= zero .OR. s >= one) THEN
  WRITE(*, *) 'SHAPE PARAMETER VALUE OUTSIDE PERMITTED RANGE'
  STOP
END IF

IF (first) THEN                        ! Initialization, if necessary
  a = one - s
  p = a/(a + s*EXP(-a))
  IF (s < vVsmall) THEN
    WRITE(*, *) 'SHAPE PARAMETER VALUE TOO SMALL'
    STOP
  END IF
  c = one/s
  uf = p*(vVsmall/a)**s
  vr = one - vVsmall
  d = a*LOG(a)
END IF

DO
  CALL RANDOM_NUMBER(r)
  IF (r >= vr) THEN
    CYCLE
  ELSE IF (r > p) THEN
    x = a - LOG((one - r)/(one - p))
    w = a*LOG(x)-d
  ELSE IF (r > uf) THEN
    x = a*(r/p)**c
    w = x
  ELSE
    random_gamma2 = zero
    RETURN
  END IF

  CALL RANDOM_NUMBER(r)
  IF (one-r <= w .AND. r > zero) THEN
    IF (r*(w + one) >= one) CYCLE
    IF (-LOG(r) <= w) CYCLE
  END IF
  EXIT
END DO

random_gamma2 = x

END FUNCTION random_gamma2


END MODULE gammadist_mod


