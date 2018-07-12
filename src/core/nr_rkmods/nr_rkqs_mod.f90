MODULE M_rkqsN 

   USE global
   USE M_rkckN, ONLY: RKCK
   IMPLICIT NONE

   PRIVATE
   PUBLIC :: RKQS

   CONTAINS

 SUBROUTINE RKQS (R,DRDT,VPAR,DVPARDT,T,HTRY,MU,EPS,RSCAL,HDID,HNEXT,T1,T2, uflag)

!##################################################################### 
 !this subroutine is the stepper and basically calls rkck to take one
 !step while monitoring truncation error to ensure accuracy. It adjusts
 !stepsize. Inputs are the what you'd expect. Outputs are new values of R,
 ! DRDT, T, hdid and hnext.
!###################################################################

 IMPLICIT NONE
 REAL(num), INTENT(INOUT), DIMENSION(3)	:: R
 REAL(num), INTENT(INOUT), DIMENSION(3)	:: DRDT
 REAL(num), INTENT(INOUT)		:: VPAR
 REAL(num), INTENT(INOUT)		:: DVPARDT
 REAL(num), INTENT(INOUT) 		:: T
 REAL(num), INTENT(IN)			:: EPS, HTRY,MU, T1, T2
 REAL(num), INTENT(IN), DIMENSION(4)	:: RSCAL
 REAL(num), INTENT(OUT)			:: HDID, HNEXT
 REAL(num), DIMENSION(4)		:: RERR 
 REAL(num), DIMENSION(3)		:: RTEMP
 REAL(num)				:: VPARTEMP
 REAL(num)				:: ERRMAX, H, HTEMP, TNEW, HNEW
 REAL(num), PARAMETER			:: SAFETY=0.9_num, PGROW=-0.2_num
 REAL(num), PARAMETER			:: PSHRINK=-0.25_num, ERRCON=1.89e-4
 INTEGER, INTENT(OUT)			:: uflag

uflag=0

 H=HTRY                   !Initial value for stepsize
 DO WHILE (.NOT.(Rlost)) 
  CALL RKCK(R,DRDT,VPAR,DVPARDT,T,H,MU,RTEMP,VPARTEMP,RERR,T1,T2)
  ERRMAX=maxval(abs(RERR(:)/RSCAL(:)))/EPS 
  !WRITE(*,1001)"time:",Tscl*(T-t1)
  !WRITE(*,1004)"RERR:",RERR
  !WRITE(*,1004)"RSCAL:",RSCAL
  !WRITE(*,1004)"abs(RE/RS):",abs(RERR(:)/RSCAL(:))
  !WRITE(*,1001)"ERRMAX:",ERRMAX
  !WRITE(*,*)'****'
  !1004 FORMAT (a,4E13.4)
  !1001 FORMAT (a,E13.4)
  
  
  IF (ERRMAX <= 1.0_num) EXIT
  HTEMP=SAFETY*H*(ERRMAX**PSHRINK)

  HNEW=sign(max(abs(HTEMP),0.1_num*abs(H)),H)
  !PRINT 1010,  abs(HTEMP), 0.1_num*abs(H), H, HNEW
  !1010 FORMAT ("abs(HTEMP)=",E11.4," 0.1|H|=",E11.4," H=",E11.4, " HNEW=",E11.4) 

  H=sign(max(abs(HTEMP),0.1_num*abs(H)),H)
  TNEW=T+H
  IF (TNEW == T) THEN
   PRINT *, 'STEPSIZE UNDERFLOW IN RKQS'
   PRINT *, 'Particle No.',pn
   PRINT *, 'T & TNEW',T,TNEW
   PRINT *, 'STOPPING.....'
   UFLAG=1
   EXIT
  END IF
 END DO
 IF (ERRMAX > ERRCON) THEN
  HNEXT=SAFETY*H*(ERRMAX**PGROW)
 ELSE
  HNEXT=5.0_num*H
 END IF
 HDID=H
 T=T+H
 R(:)=RTEMP(:)
 VPAR = VPARTEMP

 END SUBROUTINE RKQS

END MODULE M_rkqsN 
