MODULE M_rkckN 

 USE M_derivsN, ONLY: DERIVS
 USE global 
 IMPLICIT NONE

 PRIVATE
 PUBLIC :: RKCK

 CONTAINS

  SUBROUTINE RKCK (R,DRDT,VPAR,DVPARDT, T,H,MU,ROUT,VPAROUT,RERR,T1,T2)
 !####################################################################  
 !This part basically takes one Runga-Kutta step, given 6 position 
 !variables, t, and stepsize h. It calculates a fifth and fourth order
 !solution and takes the difference as the rerr (the error)
!#################################################################

  IMPLICIT NONE
  INTEGER :: I
  REAL(num) :: H, T,MU
  REAL(num), INTENT(IN)  :: VPAR, DVPARDT, T1, T2
  REAL(num), INTENT(OUT) :: VPAROUT
  REAL(num), DIMENSION(3), INTENT(IN) :: R, DRDT 
  REAL(num), DIMENSION(3), INTENT(OUT) :: ROUT
  REAL(num), DIMENSION(3) :: UGB, UC
  REAL(num), DIMENSION(4) :: RERR
  REAL(num), DIMENSION(3) :: AK2, AK3, AK4, AK5, AK6, RTEMP
  REAL(num) :: VK2, VK3, VK4, VK5, VK6, VPARTEMP
 ! REAL(num) :: A2, A3, A4, A5, A6, B21, B31, B32, B41, B42, B43
 ! REAL(num) :: B51, B52, B53, B54, B61, B62, B63, B64, B65
 ! REAL(num) :: C1, C3, C4, C6, DC1, DC3, DC4, DC5, DC6
  REAL(num), PARAMETER 			:: A2=0.2_num, A3=0.3_num, A4=0.6_num, A5=1.0_num, A6=0.875_num
  REAL(num), PARAMETER 			:: B21=0.2_num,B31=0.075_num,B32=0.225_num,B41=0.3_num,B42=-0.9_num,B43=1.2_num
  REAL(num), PARAMETER 			:: B51=-11.0_num/54.0_num,B52=2.5_num, B53=-70.0_num/27.0_num,B54=35.0_num/27.0_num
  REAL(num), PARAMETER 			:: B61=1631.0_num/55296.0_num, B62=175.0_num/512.0_num, B63=575.0_num/13824.0_num 
  REAL(num), PARAMETER			:: B64=44275.0_num/110592.0_num, B65=253.0_num/4096.0_num
  REAL(num), PARAMETER 			:: C1=37.0_num/378.0_num, C3=250.0_num/621.0_num, C4=125.0_num/594.0_num
  REAL(num), PARAMETER			:: C6=512.0_num/1771.0_num
  REAL(num), PARAMETER 			:: DC1=C1-2825.0_num/27648.0_num,DC3=C3-18575.0_num/48384.0_num 
  REAL(num), PARAMETER			:: DC4=C4-13525.0_num/55296.0_num,DC5=-277.0_num/14336.0_num,DC6=C6-0.25_num

  
  DO I = 1,3                               !first step
    RTEMP(I) = R(I) + B21*H*DRDT(I)          
  ENDDO
    VPARTEMP = VPAR + B21*H*DVPARDT  


  CALL DERIVS (T+A2*H, RTEMP, AK2, VPARTEMP,VK2,MU,T1,T2, UGB, UC)     !second step
  DO I = 1,3
    RTEMP(I) = R(I) + H *(B31*DRDT(I)+B32*AK2(I))
  ENDDO
    VPARTEMP = VPAR+ H *(B31*DVPARDT +B32*VK2)

  CALL DERIVS (T+A3*H, RTEMP, AK3,VPARTEMP,VK3,MU,T1,T2, UGB, UC)      !third step
  DO I = 1,3
    RTEMP(I) = R(I) + H*(B41*DRDT(I)+B42*AK2(I)+B43*AK3(I))
  ENDDO
    VPARTEMP = VPAR + H*(B41*DVPARDT+B42*VK2+B43*VK3)

  CALL DERIVS (T+A4*H, RTEMP, AK4,VPARTEMP,VK4,MU,T1,T2, UGB, UC)   !fourth step
  DO I = 1,3                  
   RTEMP(I) = R(I) + H*(B51*DRDT(I)+B52*AK2(I)+B53*AK3(I) &
                & + B54*AK4(I))
  ENDDO
   VPARTEMP = VPAR + H*(B51*DVPARDT+B52*VK2+B53*VK3 &
                & + B54*VK4)

 CALL DERIVS (T+A5*H, RTEMP, AK5, VPARTEMP, VK5,MU,T1,T2, UGB, UC)    !fifth step
  DO I = 1,3                  
   RTEMP(I) = R(I) + H*(B61*DRDT(I)+B62*AK2(I)+B63*AK3(I) &
                & + B64*AK4(I)+B65*AK5(I))
  ENDDO
   VPARTEMP = VPAR + H*(B61*DVPARDT+B62*VK2+B63*VK3 &
                & + B64*VK4+B65*VK5)

 CALL DERIVS (T+A6*H, RTEMP, AK6, VPARTEMP, VK6,MU,T1,T2, UGB, UC)     !sixth step
  DO I = 1,3           !accumulate increments with weights       
   ROUT(I) = R(I) + H*(C1*DRDT(I)+C3*AK3(I)+C4*AK4(I)+C6*AK6(I))
  ENDDO
   VPAROUT = VPAR + H*(C1*DVPARDT+C3*VK3+C4*VK4+C6*VK6)
   !estimate error as difference between 4th and 5th order methods
  DO I= 1,3
    RERR(I)=H*(DC1*DRDT(I)+DC3*AK3(I)+DC4*AK4(I)+DC5*AK5(I)+DC6*AK6(I))
  ENDDO
    RERR(4)=H*(DC1*DVPARDT+DC3*VK3+DC4*VK4+DC5*VK5+DC6*VK6) 
  RETURN

 END SUBROUTINE RKCK


END MODULE M_rkckN


