MODULE CMT_fields
  
  USE global
  USE M_products

  IMPLICIT NONE

  PRIVATE
  PUBLIC :: CMTFIELDS
  PRIVATE :: SUB1_X0Y0, dA0, SUB2_X0Y0, ddA0, SUB3_X0Y0
  
  CONTAINS 
!-----------------------------------------------------------------------------!   
SUBROUTINE CMTFIELDS(R,T,E,B,DBDX,DBDY,DBDZ,DBDT,DEDX,DEDY,DEDZ,DEDT,Vf)

 REAL(num)	:: X0,Y0,dX0dX,dX0dY,dY0dX,dY0dY,dX0dt,dY0dt
 REAL(num)	:: d2X0dXdY,d2X0dYdX,d2X0dX2,d2X0dY2
 REAL(num)	:: d2Y0dXdY,d2Y0dYdX,d2Y0dX2,d2Y0dY2
 REAL(num)	:: d2X0dtdX, d2X0dXdt, d2X0dtdY, d2X0dYdt
 REAL(num)	:: d2Y0dtdX, d2Y0dXdt, d2Y0dtdY, d2Y0dYdt
 REAL(num)	:: d2X0dt2, d2Y0dt2
 REAL(num)	:: dA0dX0, dA0dY0
 REAL(num)	:: d2A0dx02,d2A0dy02,d2A0dx0dy0,d2A0dy0dx0
 REAL(num)	:: DETERMINANT,der_det,xder_det,yder_det,zder_det
 REAL(num), DIMENSION(3), INTENT(OUT) :: B,E
 REAL(num), DIMENSION(3), INTENT(OUT) :: DBDX,DBDY,DBDZ,DBDT,DEDX,DEDY,DEDZ,DEDT
 REAL(num), DIMENSION(3), INTENT(IN) :: R
 REAL(num), DIMENSION(3)  :: Vf,dVfdt,dVfdx,dVfdy,dVfdz
 REAL(num), INTENT(IN) :: T
 REAL(num)	:: oneuponL=1./Lscl, oneuponT=1./Tscl
 REAL(num)	:: LuponTV=Lscl/Tscl/Vscl, oneuponTL=1./Tscl/Lscl, LoB=Lscl/Bscl, ToB=Tscl/Bscl


 !!!!!!!!!!! Magnetic Dipole !!!!!!!!!!
 ! md(1)=0.
 ! md(2)=0.
 ! md(3)=-1.
 ! RMOD2= DOT_PRODUCT(R,R)
 ! mdr  = DOT_PRODUCT(md,R)
 ! RMOD5= RMOD2**(5./2.)
 !E=0.
 !B=(3.*mdr*R-md*RMOD2)/RMOD5
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

 !!!!! Notice that X0 and Y0 and T represent dimensionless variables !!!!

 CALL SUB1_X0Y0(R,T,X0,Y0,dX0dX,dX0dY,dY0dX,dY0dY,dX0dt,dY0dt)
 CALL dA0(X0,Y0,dA0dX0, dA0dY0) 
 
! JT: B(1) and B(2) are equations 28 and 29 of Guiliani et al. (2005)
 B(1) = oneuponL *  ( dA0dX0 * dX0dY + dA0dY0 * dY0dY )
 B(2) = oneuponL * ( -(dA0dX0 * dX0dX + dA0dY0 * dY0dX) )
 B(3) =  dY0dY * Bzfinal					!??

 
 B=B/Bscl  !!!  B is made dimensionless

 !!!!! Velocity field is also in dimensionless units !!!!!
 !!!!! Notice the factors L/T/Vscl

 DETERMINANT = dX0dX * dY0dY - dX0dY * dY0dX
 Vf(1)= LuponTV*( -dX0dt*dY0dY + dX0dY*dY0dt )/DETERMINANT
 Vf(2)= LuponTV*( -dY0dt*dX0dX + dY0dX*dX0dt )/DETERMINANT
 Vf(3)= 0.0_num

 E = -CROSS(Vf,B)   ! This electric field is dimensionless
 ! print*,'E(3)=',E(3),'****'
 !!!! Notice that the dimensionless electric field E(3) can also be calculated 
 !!!! as
 ! E(3) = -(1.0/Tscl)*(dA0dX0*dX0dt + dA0dY0*dY0dt)/(Vscl*B0)
 !!!! see notebook for more information !!!!!!!!!!!!!!!!
 !print*,'E(3)=',E(3)

 !!!!!!! Derivatives !!!!!!!!!!!!!!!!!!

 CALL SUB2_X0Y0(R,T,d2X0dXdY,d2X0dYdX,d2X0dX2,d2X0dY2,d2Y0dXdY,d2Y0dYdX,d2Y0dX2,d2Y0dY2)
 CALL SUB3_X0Y0(R,T,d2X0dtdX,d2X0dXdt,d2X0dtdY,d2X0dYdt,d2Y0dtdX,d2Y0dXdt,d2Y0dtdY,d2Y0dYdt,d2X0dt2,d2Y0dt2 )
 CALL ddA0(X0,Y0,d2A0dx02,d2A0dy02,d2A0dx0dy0,d2A0dy0dx0)

 ! Now calculate derivative of Velocity field w.r.t. time

 !!! derivative of DETERMINANT is done with respect to normalized time !!!
 der_det=d2X0dtdX*dY0dY+dX0dX*d2Y0dtdY - d2X0dtdY*dY0dX-dX0dY*d2Y0dtdX
 xder_det=d2X0dX2*dY0dY+dX0dX*d2Y0dXdY-d2X0dXdY*dY0dX-dX0dY*d2Y0dX2
 yder_det=d2X0dXdY*dY0dY+dX0dX*d2Y0dY2-d2X0dY2*dY0dX-dX0dY*d2Y0dXdY
 zder_det=0.

 !!! Time derivatives of dmensionless Vf with respect to dimensionless time!!
 dVfdt(1)=LuponTV*((-d2X0dt2*dY0dY-dX0dt*d2Y0dtdY+d2X0dtdY*dY0dt+dX0dY*d2Y0dt2)/DETERMINANT+ &
          (-dX0dt*dY0dY+dX0dY*dY0dt)*(-der_det/DETERMINANT/DETERMINANT))
 dVfdt(2)=LuponTV*((-d2Y0dt2*dX0dX-dY0dt*d2X0dtdX+d2Y0dtdX*dX0dt+dY0dX*d2X0dt2)/DETERMINANT+ &
          (-dY0dt*dX0dX+dY0dX*dX0dt)*(-der_det/DETERMINANT/DETERMINANT))
 dVfdt(3)=0.
 
!---- SPACE DERIVATIVES: dimensionless quantities
 dVfdx(1)=LuponTV*((-d2X0dXdt*dY0dY-dX0dt*d2Y0dXdY+d2X0dXdY*dY0dt+dX0dY*d2Y0dXdt)/DETERMINANT+ &
          (-dX0dt*dY0dY + dX0dY*dY0dt )*(-xder_det/DETERMINANT/DETERMINANT))
 dVfdx(2)=LuponTV*((-d2Y0dX2*dX0dX-dY0dX*d2X0dX2+d2Y0dX2*dX0dt+dY0dX*d2X0dXdt)/DETERMINANT +&
         (-dY0dt*dX0dX+dY0dX*dX0dt)*(-xder_det/DETERMINANT/DETERMINANT))
 dVfdx(3)=0.
!----
 dVfdy(1)=LuponTV*((-d2X0dYdt*dY0dY-dX0dt*d2Y0dY2+d2X0dY2*dY0dt+dX0dY*d2Y0dYdt )/DETERMINANT+ &
          (-dX0dt*dY0dY+dX0dY*dY0dt)*(-yder_det/DETERMINANT/DETERMINANT))
 dVfdy(2)=LuponTV*((-d2Y0dYdt*dX0dX -dY0dt*d2X0dYdX+d2Y0dYdX*dX0dt+dY0dX*d2X0dYdt )/DETERMINANT + &
          (-dY0dt*dX0dX+dY0dX*dX0dt)*(-yder_det/DETERMINANT/DETERMINANT))
 dVfdy(3)= 0.
!----
 dVfdz(1)= 0.
 dVfdz(2)= 0.
 dVfdz(3)= 0.

!!!! Check: LHS and RHS should be equal !!!!!!
 !print*,'1-LHS=',(der_det+Vf(1)*xder_det+Vf(2)*yder_det+Vf(3)*zder_det)/DETERMINANT
 !print*,'2-RHS=',-(dVfdx(1)+dVfdy(2)+dVfdz(3))
 !print*,'****'

!!!! Magnetic field derivatives in dimensional form !!!
 !dBxdx
 dBdx(1)=oneuponL*oneuponL*( (d2A0dx02*dx0dx+d2A0dy0dx0*dy0dx)*dx0dy + dA0dx0*d2x0dxdy + &
                    (d2A0dx0dy0*dx0dx + d2A0dy02*dy0dx)*dy0dy + dA0dy0*d2y0dxdy )
 !dBxdy
 dBdy(1)=oneuponL*oneuponL*( (d2A0dx02*dx0dy+d2A0dy0dx0*dy0dy)*dx0dy + dA0dx0*d2x0dy2 +  &
                    (d2A0dx0dy0*dx0dy + d2A0dy02*dy0dy)*dy0dy + dA0dy0*d2y0dy2 )
 !dBxdz
 dBdz(1)=0.
 !dBydx
 dBdx(2)=-oneuponL*oneuponL*( (d2A0dx02*dx0dx+d2A0dy0dx0*dy0dx)*dx0dx + dA0dx0*d2x0dx2 + &
                     (d2A0dx0dy0*dx0dx + d2A0dy02*dy0dx)*dy0dx + dA0dy0*d2y0dx2 )
 !dBydy
 dBdy(2)=-oneuponL*oneuponL*( (d2A0dx02*dx0dy+d2A0dy0dx0*dy0dy)*dx0dx + dA0dx0*d2x0dydx + &
                      (d2A0dx0dy0*dx0dy + d2A0dy02*dy0dy)*dy0dx + dA0dy0*d2y0dydx )
 !dBydz
 dBdz(2)=0.
 !dBzdx
 dBdx(3)=oneuponL*d2y0dydx*Bzfinal
 !dBzdy
 dBdy(3)=oneuponL*d2y0dy2*Bzfinal
 !dBzdz
 dBdz(3)= 0.
 !dBxdt
 dBdt(1)=oneuponTL*(d2A0dx02 * dx0dt * dx0dy + d2A0dy0dx0 * dy0dt *dx0dy + &
                      d2A0dx0dy0 * dx0dt * dy0dy + d2A0dy02 * dy0dt* dy0dy +  &
                      dA0dx0 * d2x0dtdy + dA0dy0 * d2y0dtdy)
 !dBydt
 dBdt(2)=-oneuponTL*( d2A0dx02 * dx0dt * dx0dx + d2A0dy0dx0 * dy0dt *dx0dx + &
                      d2A0dx0dy0 * dx0dt * dy0dx + d2A0dy02 * dy0dt* dy0dx +   &
                      dA0dx0 * d2x0dtdx + dA0dy0 * d2y0dtdx)
 !dBzdt
 dBdt(3)=oneuponT*d2y0dtdy *Bzfinal

 !! Make derivatives dimensionless !!!!
 DBDX = LoB*DBDX 
 DBDY = LoB*DBDY
 DBDZ = LoB*DBDZ
 DBDT = ToB*DBDT 

 IF (T.eq.1.234) THEN
  PRINT *,DBDT
 ENDIF

 !!! The dimensionless dEdt is !!
 dEdt=-(cross(dVfdt,B)+cross(Vf,DBDT))
 !print*,'dEdt(3)=',dEdt(3),'***'

 !!! The dimensionless dEdx is !!
 dEdx=-(cross(dVfdx,B)+cross(Vf,DBDX)) 
! print*,'dEdx(3)=',dEdx(3)

 !!! The dimensionless dEdy is !!
 dEdy=-(cross(dVfdy,B)+cross(Vf,DBDY))
 !print*,'dEdy(3)=',dEdy(3)

 !!! The dimensionless dEdz is !!
 dEdz=-(cross(dVfdz,B)+cross(Vf,DBDZ))
 !print*,'dEdz(3)=',dEdz(3)


 !!!! Dimensionless derivatives of Ez can also be calculated as follows !!!
 !!!! (L/Vscl/B0) and (Tscl/Vscl/B0) reduce quantities to dimensionless.
 !!dEzdx
 !dEdx(3)=-(oneuponL/Tscl)*(L/Vscl/B0)* &
 !                 ((d2A0dx02*dx0dx+d2A0dy0dx0*dy0dx)*dx0dt+dA0dx0*d2x0dxdt+&
 !                 (d2A0dx0dy0*dx0dx+d2A0dy02*dy0dx)*dy0dt+dA0dy0*d2y0dxdt)
 !!dEzdy
 !dEdy(3)=-(oneuponL/Tscl)*(L/Vscl/B0)* &
 !                     ((d2A0dx02*dx0dy+d2A0dy0dx0*dy0dy)*dx0dt+dA0dx0*d2x0dydt+&
 !                     (d2A0dx0dy0*dx0dy+d2A0dy02*dy0dy)*dy0dt+dA0dy0*d2y0dydt)
 !!dEzdz
 !dEdz(3)=0.
 !
 !!dEdt
 !dEdt(3)=-oneuponT*oneuponT*(Tscl/Vscl/B0)* &
 !                    (d2A0dx02 * dx0dt*dx0dt + d2A0dy0dx0 *dy0dt*dx0dt + &
 !                     dA0dx0 * d2x0dt2 + &
 !                     d2A0dx0dy0 * dx0dt*dy0dt + d2A0dy02*dy0dt*dy0dt + &
 !                     dA0dy0 * d2y0dt2 )
 !print*,'DEDT(3)=',DEDT(3)
 !print*,'DEDX(3)=',DEDX(3)
 !print*,'DEDY(3)=',DEDY(3)
 !print*,'DEDZ(3)=',DEDZ(3)

 !print*,X0   ! OK
 !print*,Y0   ! OK
 !print*,dY0dY ! 0k
 !print*,dY0dX ! OK
 !print*,dY0dT ! OK
 !print*,dX0dX ! OK
 !print*,dX0dY ! OK
 !print*,dX0dt ! OK
 !print*,d2X0dTdX !OK
 !print*,d2X0dTdY !OK
 !print*,d2Y0dTdX !OK
 !print*,d2Y0dTdY !Ok
 !print*,d2X0dXdY !OK
 !print*,d2X0dYdX !OK
 !print*,d2X0dX2  !OK
 !print*,d2X0dY2  !OK
 !print*,d2X0dt2  !OK
 !print*,d2Y0dXdY  !Ok
 !print*,d2Y0dYdX  !OK
 !print*,d2Y0dX2   !OK
 !print*,d2Y0dY2   !OK
 !print*,d2Y0dt2   !OK 

 !print*,dBdy
 !print*,dBdt ! OK
 !print*,E    ! OK
 !print*,DEDT
 !print*,DEDX
 !print*,DEDY
 !print*,DEDZ

 !B(1)=0
 !B(2)=0
 !B(3)=B0*T*exp(-R(2))
 !DBDX=0; DBDY=0; DBDZ=0; 
 !DBDY(3)=(-oneuponL)*B(3)
 !DBDT(1)=0
 !DBDT(2)=0
 !DBDT(3)=B0*exp(-R(2))/Tscl
 !E=0
 !DEDX=0; DEDY=0; DEDZ=0; DEDT=0
 !Vf=0
 !dVfdx=0; dVfdy=0; dVfdz=0; dVfdt=0
 !
 !B=B/B0
 !DBDY=(L/B0)*DBDY
 !DBDT=(Tscl/B0)*DBDT
END SUBROUTINE CMTFIELDS
!-------------------------------------------------------------------------------------------!
SUBROUTINE SUB1_X0Y0(R,T,X0,Y0,dX0dX,dX0dY,dY0dX,dY0dY,dX0dt,dY0dt)
!JT+ reads in a given position vector at time t, 
! and calculates Lagrangian transformation, x0(=x_inf) and y0(=y_inf)
! based on eq. (36) of Guiliani et al (2005)

 REAL(num), DIMENSION(3), INTENT(IN)	:: R
 REAL(num), INTENT(IN)  		:: T
 REAL(num), INTENT(OUT)			:: X0, Y0, dX0dX, dX0dY, dY0dX, dY0dY, dX0dt, dY0dt
 REAL(num)				:: a=0.9,b=0.9

 !!!!!!  please notice Lv/L in the expressions  !!!!!

 X0=R(1)

 IF ( (1.+(R(2)/((cc*T)**esp))) .le.0.) THEN
  PRINT *, '(1.+(R(2)/((cc*T)**esp)))=',(1.+(R(2)/((cc*T)**esp)))
  PRINT *, 'R(2)=',R(2)
  PRINT *,'T=',T
 END IF

 Y0=(cc*T)**esp*log (1.+(R(2)/((cc*T)**esp)))*  &	
      (1.+tanh((R(2)-Lv/Lscl)*a))*0.5 +           &
       (1.-tanh((R(2)-Lv/Lscl)*b))*0.5*R(2)

!print*,X0,Y0

!JT: are these boundary conditions?
 dX0dX=1.
 dX0dY=0.
 dY0dX=0.
 dY0dY=0.5*(1.+tanh((R(2)-Lv/Lscl)*a))/(1.+R(2)/(cc*T)**esp)+        &
      0.5*(cc*T)**esp*log(1.+R(2)/(cc*T)**esp)*               &
      (1.-tanh((R(2)-Lv/Lscl)*a)**2)*a-0.5*                         &
      (1.-tanh((R(2)-Lv/Lscl)*b)**2)*b*R(2)+0.5-0.5*tanh((R(2)-Lv/Lscl)*b) 

 dX0dt=0.
 dY0dt=0.5*(cc*T)**esp*esp*log(1.+R(2)/(cc*T)**esp)*          &
      (1.+tanh((R(2)-Lv/Lscl)*a))/T-0.5*R(2)*esp*                  &
      (1.+tanh((R(2)-Lv/Lscl)*a))/(T*(1.+R(2)/(cc*T)**esp))          
 !**********************************************************
 ! X0=R(1)
 ! Y0=R(2) 

 ! dX0dX=1.
 ! dX0dY=0.
 ! dY0dX=0.
 ! dY0dY=1.

 !dX0dt=0.
 !dY0dt=0.
     
END SUBROUTINE SUB1_X0Y0
!-------------------------------------------------------------------------------------------!
SUBROUTINE SUB2_X0Y0(R,T,d2X0dXdY,d2X0dYdX,d2X0dX2,d2X0dY2,d2Y0dXdY,d2Y0dYdX,d2Y0dX2,d2Y0dY2)
 
 REAL(num), DIMENSION(3), INTENT(IN) :: R
 REAL(num), INTENT(IN)  	:: T
 REAL(num), INTENT(OUT) 	:: d2X0dXdY,d2X0dYdX,d2X0dX2,d2X0dY2
 REAL(num), INTENT(OUT) 	:: d2Y0dXdY,d2Y0dYdX,d2Y0dX2,d2Y0dY2
 REAL(num)              	:: y
 REAL(num)              	:: a=0.9,b=0.9

 y=R(2)
 !!!!  please notice Lv/L in the expressions !!!!
 
 d2y0dy2= -0.1e1 / (0.1e1 + y / (cc * t) ** esp) ** 2 * (0.1e1 + tan&
     &h((y - Lv/Lscl) * a)) / (cc * t) ** esp / 0.2e1 + 0.1e1 / (0.1e1 + y /&
     & (cc * t) ** esp) * (0.1e1 - tanh((y - Lv/Lscl) * a) ** 2) * a - (cc& 
     &* t) ** esp * log(0.1e1 + y / (cc * t) ** esp) * tanh((y - Lv/Lscl) * a&
     &) * (0.1e1 - tanh((y - Lv/Lscl) * a) ** 2) * a*a + tanh((y - Lv/Lscl)& 
     &* b) * (0.1e1 - tanh((y - Lv/Lscl) * b) ** 2) * b*b * y - (0.1e1& 
     &- tanh((y - Lv/Lscl) * b) ** 2) * b

 d2x0dxdy = 0.
 d2x0dydx = 0.
 d2x0dx2 = 0.
 d2x0dy2 = 0.
 d2y0dxdy = 0.
 d2y0dydx = 0.
 d2y0dx2 = 0.

!*********************************************************
!d2y0dy2=0.

END SUBROUTINE SUB2_X0Y0
!-------------------------------------------------------------------------------------------!
SUBROUTINE SUB3_X0Y0(R,T,d2X0dtdX,d2X0dXdt,d2X0dtdY,d2X0dYdt, &
                         d2Y0dtdX,d2Y0dXdt,d2Y0dtdY,d2Y0dYdt, &
                                           d2X0dt2,d2Y0dt2 )
 REAL(num), DIMENSION(3), INTENT(IN)	:: R
 REAL(num), INTENT(IN)  		:: T
 REAL(num), INTENT(OUT) 		:: d2X0dtdX,d2X0dXdt,d2X0dtdY,d2X0dYdt
 REAL(num), INTENT(OUT) 		:: d2Y0dtdX,d2Y0dXdt,d2Y0dtdY,d2Y0dYdt
 REAL(num), INTENT(OUT) 		:: d2X0dt2, d2Y0dt2
 REAL(num)             			:: y
 REAL(num)				:: a=0.9
 
y=R(2)

!!!!! Lv/L in the expressions !!!!!

d2x0dtdx =0.
d2x0dxdt=d2x0dtdx

d2x0dtdy =0.
d2x0dydt =d2x0dtdy

d2y0dtdx =0.
d2y0dxdt=d2y0dtdx

d2y0dtdy = 0.1e1 / (0.1e1 + y / (cc * t) ** esp) ** 2 * (0.1e1 + tanh &
     ((y - Lv/Lscl) * a)) * y / (cc * t) ** esp * esp / t / 0.2e1 + (cc * & 
     t) ** esp * esp / t * log(0.1e1 + y / (cc * t) ** esp) * (0.1e1 - &
     tanh((y - Lv/Lscl) * a) ** 2) * a / 0.2e1 - y * esp / t / (0.1e1 + y &
     / (cc * t) ** esp) * (0.1e1 - tanh((y - Lv/Lscl) * a) ** 2) * a / 0.2e1

d2y0dydt=d2y0dtdy

d2x0dt2 = 0.

d2y0dt2 = (cc * t) ** esp * esp ** 2 / t ** 2 * log(0.1e1 + y / (cc &
      * t) ** esp) * (0.1e1 + tanh((y - Lv/Lscl) * a)) / 0.2e1 - (cc * t) ** &
      esp * esp / t ** 2 * log(0.1e1 + y / (cc * t) ** esp) * (0.1e1 +  &
      tanh((y - Lv/Lscl) * a)) / 0.2e1 - esp ** 2 / t ** 2 * y / (0.1e1 + y  &
      / (cc * t) ** esp) * (0.1e1 + tanh((y - Lv/Lscl) * a)) / 0.2e1 + y *   &
     esp / t ** 2 / (0.1e1 + y / (cc * t) ** esp) * (0.1e1 + tanh((y -  &
     Lv/Lscl) * a)) / 0.2e1 - y ** 2 * esp ** 2 / t ** 2 / (0.1e1 + y / (cc  &
      * t) ** esp) ** 2 * (0.1e1 + tanh((y - Lv/Lscl) * a)) / (cc * t) **    &
     esp / 0.2e1

!******************************************************
!d2y0dtdy=0.
!d2y0dydt=0.
!d2y0dt2=0.

END SUBROUTINE SUB3_X0Y0
!-------------------------------------------------------------------------------------------!
SUBROUTINE dA0 (X0,Y0,dA0dX0, dA0dY0)
!+ returns differentials dA0/dx0 and dA0/dy0 for specific values of x0 and y0
! based on equation 31 of Guiliani et al (2005) 
!JT edited - all numbers forced to be real, and powers unrolled
 
 REAL(num), INTENT(IN)  :: X0,Y0
 REAL(num), INTENT(OUT) :: dA0dX0, dA0dY0


!!!!!!!!!!! These are the derivatives with respect to
!!!!!!!!!!! the dimensional variables, not to be used here
! dA0dX0 = 32.*(Y0+d)*c1*X0*L/  &
!         ((4.*X0**2-4.*X0*L+L**2+4*Y0**2+8.*Y0*d+4*d**2)*   &
!         (4.*X0**2+4*X0*L+L**2+4*Y0**2+8*Y0*d+4*d**2)) 
! dA0dY0 = -4.*c1*L*(4.*X0**2-L**2-4*Y0**2-8*Y0*d-4*d**2)/  & 
!         ((4*X0**2-4*X0*L+L**2+4*Y0**2+8*Y0*d+4*d**2)*     &
!         (4*X0**2+4*X0*L+L**2+4*Y0**2+8*Y0*d+4*d**2)) 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!dA0dX0 = 32.*c1*(Y0*L+d)*L**3*X0/  &
! ((4*L**2*X0**2-4*L**2*X0+L**2+4*Y0**2*L**2+8*Y0*L*d+4*d**2)* &
! (4*L**2*X0**2+4*L**2*X0+L**2+4*Y0**2*L**2+8*Y0*L*d+4*d**2))
!
!dA0dY0 = 4.*c1*L**2*  &
! (-4*L**2*X0**2+L**2+4*Y0**2*L**2+8*Y0*L*d+4*d**2)/  &
! ((4*L**2*X0**2-4*L**2*X0+L**2+4*Y0**2*L**2+8*Y0*L*d+4*d**2)*  &
! (4*L**2*X0**2+4*L**2*X0+L**2+4*Y0**2*L**2+8*Y0*L*d+4*d**2))

dA0dX0 = 32.*c1*(Y0*Lscl+d)*Lscl*Lscl*Lscl*X0/  &
 ((4.*Lscl*Lscl*X0*X0-4.*Lscl*Lscl*X0+Lscl*Lscl+4.*Y0*Y0*Lscl*Lscl+8.*Y0*Lscl*d+4.*d*d)* &
 (4.*Lscl*Lscl*X0*X0+4.*Lscl*Lscl*X0+Lscl*Lscl+4.*Y0*Y0*Lscl*Lscl+8.*Y0*Lscl*d+4.*d*d))

dA0dY0 = 4.*c1*Lscl*Lscl*  &
 (-4.*Lscl*Lscl*X0*X0+Lscl*Lscl+4.*Y0*Y0*Lscl*Lscl+8.*Y0*Lscl*d+4.*d*d)/  &
 ((4.*Lscl*Lscl*X0*X0-4.*Lscl*Lscl*X0+Lscl*Lscl+4*Y0*Y0*Lscl*Lscl+8.*Y0*Lscl*d+4.*d*d)*  &
 (4.*Lscl*Lscl*X0*X0+4.*Lscl*Lscl*X0+Lscl*Lscl+4.*Y0*Y0*Lscl*Lscl+8.*Y0*Lscl*d+4.*d*d))


END SUBROUTINE dA0
!-------------------------------------------------------------------------------------------!
SUBROUTINE ddA0 (X0,Y0,d2A0dx02, d2A0dy02, d2A0dx0dy0, d2A0dy0dx0 )
! second differentials of flux function A0,
! also based on Eq. 31 of Guiliani et al (2005)

!+ edit to make all numbers real and unroll powers 
 REAL(num), INTENT(OUT) :: d2A0dx02, d2A0dy02, d2A0dx0dy0, d2A0dy0dx0
 REAL(num), INTENT(IN)  :: X0,Y0

 d2A0dx02 = 32. * c1 * (y0 * Lscl + d) * Lscl*Lscl*Lscl / (4. * Lscl*Lscl * x0*x0 - & 
       4. * Lscl*Lscl * x0 + Lscl*Lscl + 4. * y0*y0 * Lscl*Lscl + 8. * y0 * Lscl * d + &
        4. * d*d) / (4. * Lscl*Lscl * x0*x0 + 4. * Lscl*Lscl * x0 + Lscl*Lscl +  &
       4. * y0*y0 * Lscl*Lscl + 8. * y0 * Lscl * d + 4. * d*d) - 32. * c1 * (y0 &
       * Lscl + d) * Lscl*Lscl*Lscl * x0 / (4. * Lscl*Lscl * x0*x0 - 4. * Lscl*Lscl * x0 &
        + Lscl*Lscl + 4. * y0*y0 * Lscl*Lscl + 8. * y0 * Lscl * d + 4. * d*d) ** &
        2 / (4. * Lscl*Lscl * x0*x0 + 4. * Lscl*Lscl * x0 + Lscl*Lscl + 4. * y0 **  & 
       2 * Lscl*Lscl + 8. * y0 * Lscl * d + 4. * d*d) * (8. * Lscl*Lscl * x0 - 4. * &
        Lscl*Lscl) - 32. * c1 * (y0 * Lscl + d) * Lscl*Lscl*Lscl * x0 / (4. * Lscl*Lscl * x0*x0 &
         - 4. * Lscl*Lscl * x0 + Lscl*Lscl + 4. * y0*y0 * Lscl*Lscl + 8. * y0 * &
        Lscl * d + 4. * d*d) / (4. * Lscl*Lscl * x0*x0 + 4. * Lscl*Lscl * x0 + Lscl*Lscl + &
        4. * y0*y0 * Lscl*Lscl + 8. * y0 * Lscl * d + 4. * d*d) ** 2 * &
       (8. * Lscl*Lscl * x0 + 4. * Lscl*Lscl)

d2A0dy02 = 4. * c1 * Lscl*Lscl * (8. * y0 * Lscl*Lscl + 8. * Lscl * d) / (4. * Lscl*Lscl &
      * x0*x0 - 4. * Lscl*Lscl * x0 + Lscl*Lscl + 4. * y0*y0 * Lscl*Lscl + &
       8. * y0 * Lscl * d + 4. * d*d) / (4. * Lscl*Lscl * x0*x0 + 4. * Lscl*Lscl &
       * x0 + Lscl*Lscl + 4. * y0*y0 * Lscl*Lscl + 8. * y0 * Lscl * d + 4. * d*d) & 
       - 4. * c1 * Lscl*Lscl * (-4. * Lscl*Lscl * x0*x0 + Lscl*Lscl + 4. * y0*y0  &
       * Lscl*Lscl + 8. * y0 * Lscl * d + 4. * d*d) / (4. * Lscl*Lscl * x0*x0 &
       - 4. * Lscl*Lscl * x0 + Lscl*Lscl + 4. * y0*y0 * Lscl*Lscl + 8. * y0 * Lscl *  &
      d + 4. * d*d) ** 2 / (4. * Lscl*Lscl * x0*x0 + 4. * Lscl*Lscl * x0 + Lscl*Lscl &
       + 4. * y0*y0 * Lscl*Lscl + 8. * y0 * Lscl * d + 4. * d*d) * (8. * &
       y0 * Lscl*Lscl + 8. * Lscl * d) - 4. * c1 * Lscl*Lscl * (-4. * Lscl*Lscl * x0*x0 &
       + Lscl*Lscl + 4. * y0*y0 * Lscl*Lscl + 8. * y0 * Lscl * d + 4. * d*d)  &
      / (4. * Lscl*Lscl * x0*x0 - 4. * Lscl*Lscl * x0 + Lscl*Lscl + 4. * y0*y0 * &
       Lscl*Lscl + 8. * y0 * Lscl * d + 4. * d*d) / (4. * Lscl*Lscl * x0*x0 + 4. &
       * Lscl*Lscl * x0 + Lscl*Lscl + 4. * y0*y0 * Lscl*Lscl + 8. * y0 * Lscl * d +  &
      4. * d*d) ** 2 * (8. * y0 * Lscl*Lscl + 8. * Lscl * d)

d2A0dx0dy0 = 32. * c1 * Lscl ** 4. * x0 / (4. * Lscl*Lscl * x0*x0 - 4. * Lscl*Lscl &
     * x0 + Lscl*Lscl + 4. * y0*y0 * Lscl*Lscl + 8. * y0 * Lscl * d + 4. * d*d) &
      / (4. * Lscl*Lscl * x0*x0 + 4. * Lscl*Lscl * x0 + Lscl*Lscl + 4. * y0*y0  &
      * Lscl*Lscl + 8. * y0 * Lscl * d + 4. * d*d) - 32. * c1 * (y0 * Lscl + d &
     ) * Lscl*Lscl*Lscl * x0 / (4. * Lscl*Lscl * x0*x0 - 4. * Lscl*Lscl * x0 + Lscl*Lscl &
      + 4. * y0*y0 * Lscl*Lscl + 8. * y0 * Lscl * d + 4. * d*d) ** 2 / (4. * &
      Lscl*Lscl * x0*x0 + 4. * Lscl*Lscl * x0 + Lscl*Lscl + 4. * y0*y0 * Lscl*Lscl  &
     + 8. * y0 * Lscl * d + 4. * d*d) * (8. * y0 * Lscl*Lscl + 8. * Lscl * d) - &
      32. * c1 * (y0 * Lscl + d) * Lscl*Lscl*Lscl * x0 / (4. * Lscl*Lscl * x0*x0 - 4.  &
     * Lscl*Lscl * x0 + Lscl*Lscl + 4. * y0*y0 * Lscl*Lscl + 8. * y0 * Lscl * d + 4. &
      * d*d) / (4. * Lscl*Lscl * x0*x0 + 4. * Lscl*Lscl * x0 + Lscl*Lscl + 4.  &
     * y0*y0 * Lscl*Lscl + 8. * y0 * Lscl * d + 4. * d*d) ** 2 * (8. * y0 * &
      Lscl*Lscl + 8. * Lscl * d)

d2A0dy0dx0 = -32. * c1 * Lscl ** 4. * x0 / (4. * Lscl*Lscl * x0*x0 - 4. * Lscl*Lscl &
      * x0 + Lscl*Lscl + 4. * y0*y0 * Lscl*Lscl + 8. * y0 * Lscl * d + 4. * d*d) &
       / (4. * Lscl*Lscl * x0*x0 + 4. * Lscl*Lscl * x0 + Lscl*Lscl + 4. * y0*y0 &
      * Lscl*Lscl + 8. * y0 * Lscl * d + 4. * d*d) - 4. * c1 * Lscl*Lscl * (- &
     4. * Lscl*Lscl * x0*x0 + Lscl*Lscl + 4. * y0*y0 * Lscl*Lscl + 8. * y0 * Lscl  &
     * d + 4. * d*d) / (4. * Lscl*Lscl * x0*x0 - 4. * Lscl*Lscl * x0 + Lscl*Lscl &
      + 4. * y0*y0 * Lscl*Lscl + 8. * y0 * Lscl * d + 4. * d*d) ** 2 / (4. &
      * Lscl*Lscl * x0*x0 + 4. * Lscl*Lscl * x0 + Lscl*Lscl + 4. * y0*y0 * Lscl*Lscl &
      + 8. * y0 * Lscl * d + 4. * d*d) * (8. * Lscl*Lscl * x0 - 4. * Lscl*Lscl &
     ) - 4. * c1 * Lscl*Lscl * (-4. * Lscl*Lscl * x0*x0 + Lscl*Lscl + 4. * y0*y0   &
     * Lscl*Lscl + 8. * y0 * Lscl * d + 4. * d*d) / (4. * Lscl*Lscl * x0*x0  &
     - 4. * Lscl*Lscl * x0 + Lscl*Lscl + 4. * y0*y0 * Lscl*Lscl + 8. * y0 * Lscl * d &
      + 4. * d*d) / (4. * Lscl*Lscl * x0*x0 + 4. * Lscl*Lscl * x0 + Lscl*Lscl  &
     + 4. * y0*y0 * Lscl*Lscl + 8. * y0 * Lscl * d + 4. * d*d) ** 2 * (8. *  &
     Lscl*Lscl * x0 + 4. * Lscl*Lscl)


END SUBROUTINE ddA0
!-------------------------------------------------------------------------------------------!
!SUBROUTINE ddA0 (X0,Y0,d2A0dx02, d2A0dy02, d2A0dx0dy0, d2A0dy0dx0 )
!
! REAL, INTENT(OUT) :: d2A0dx02, d2A0dy02, d2A0dx0dy0, d2A0dy0dx0
! REAL, INTENT(IN)  :: X0,Y0
!
! d2A0dx02 = 32 * c1 * (y0 * L + d) * L ** 3 / (4 * L ** 2 * x0 ** 2 - & 
!     4 * L ** 2 * x0 + L ** 2 + 4 * y0 ** 2 * L ** 2 + 8 * y0 * L * d + &
!      4 * d ** 2) / (4 * L ** 2 * x0 ** 2 + 4 * L ** 2 * x0 + L ** 2 +  &
!     4 * y0 ** 2 * L ** 2 + 8 * y0 * L * d + 4 * d ** 2) - 32 * c1 * (y&
!    &0 * L + d) * L ** 3 * x0 / (4 * L ** 2 * x0 ** 2 - 4 * L ** 2 * x0 &
!      + L ** 2 + 4 * y0 ** 2 * L ** 2 + 8 * y0 * L * d + 4 * d ** 2) ** &
!      2 / (4 * L ** 2 * x0 ** 2 + 4 * L ** 2 * x0 + L ** 2 + 4 * y0 **  & 
!     2 * L ** 2 + 8 * y0 * L * d + 4 * d ** 2) * (8 * L ** 2 * x0 - 4 * &
!      L ** 2) - 32 * c1 * (y0 * L + d) * L ** 3 * x0 / (4 * L ** 2 * x0 &
!      ** 2 - 4 * L ** 2 * x0 + L ** 2 + 4 * y0 ** 2 * L ** 2 + 8 * y0 * &
!      L * d + 4 * d ** 2) / (4 * L ** 2 * x0 ** 2 + 4 * L ** 2 * x0 + L &
!      ** 2 + 4 * y0 ** 2 * L ** 2 + 8 * y0 * L * d + 4 * d ** 2) ** 2 * &
!      (8 * L ** 2 * x0 + 4 * L ** 2)
!
!d2A0dy02 = 4 * c1 * L ** 2 * (8 * y0 * L ** 2 + 8 * L * d) / (4 * L &
!     ** 2 * x0 ** 2 - 4 * L ** 2 * x0 + L ** 2 + 4 * y0 ** 2 * L ** 2 + &
!      8 * y0 * L * d + 4 * d ** 2) / (4 * L ** 2 * x0 ** 2 + 4 * L ** 2 &
!      * x0 + L ** 2 + 4 * y0 ** 2 * L ** 2 + 8 * y0 * L * d + 4 * d ** & 
!     2) - 4 * c1 * L ** 2 * (-4 * L ** 2 * x0 ** 2 + L ** 2 + 4 * y0 ** &
!      2 * L ** 2 + 8 * y0 * L * d + 4 * d ** 2) / (4 * L ** 2 * x0 ** 2 &
!      - 4 * L ** 2 * x0 + L ** 2 + 4 * y0 ** 2 * L ** 2 + 8 * y0 * L *  &
!     d + 4 * d ** 2) ** 2 / (4 * L ** 2 * x0 ** 2 + 4 * L ** 2 * x0 + L &
!      ** 2 + 4 * y0 ** 2 * L ** 2 + 8 * y0 * L * d + 4 * d ** 2) * (8 * &
!      y0 * L ** 2 + 8 * L * d) - 4 * c1 * L ** 2 * (-4 * L ** 2 * x0 ** &
!      2 + L ** 2 + 4 * y0 ** 2 * L ** 2 + 8 * y0 * L * d + 4 * d ** 2)  &
!     / (4 * L ** 2 * x0 ** 2 - 4 * L ** 2 * x0 + L ** 2 + 4 * y0 ** 2 * &
!      L ** 2 + 8 * y0 * L * d + 4 * d ** 2) / (4 * L ** 2 * x0 ** 2 + 4 &
!      * L ** 2 * x0 + L ** 2 + 4 * y0 ** 2 * L ** 2 + 8 * y0 * L * d +  &
!     4 * d ** 2) ** 2 * (8 * y0 * L ** 2 + 8 * L * d)
!
!d2A0dx0dy0 = 32 * c1 * L ** 4 * x0 / (4 * L ** 2 * x0 ** 2 - 4 * L ** &
!     2 * x0 + L ** 2 + 4 * y0 ** 2 * L ** 2 + 8 * y0 * L * d + 4 * d ** &
!      2) / (4 * L ** 2 * x0 ** 2 + 4 * L ** 2 * x0 + L ** 2 + 4 * y0 ** &
!      2 * L ** 2 + 8 * y0 * L * d + 4 * d ** 2) - 32 * c1 * (y0 * L + d &
!     ) * L ** 3 * x0 / (4 * L ** 2 * x0 ** 2 - 4 * L ** 2 * x0 + L ** 2 &
!      + 4 * y0 ** 2 * L ** 2 + 8 * y0 * L * d + 4 * d ** 2) ** 2 / (4 * &
!      L ** 2 * x0 ** 2 + 4 * L ** 2 * x0 + L ** 2 + 4 * y0 ** 2 * L **  &
!     2 + 8 * y0 * L * d + 4 * d ** 2) * (8 * y0 * L ** 2 + 8 * L * d) - &
!      32 * c1 * (y0 * L + d) * L ** 3 * x0 / (4 * L ** 2 * x0 ** 2 - 4  &
!     * L ** 2 * x0 + L ** 2 + 4 * y0 ** 2 * L ** 2 + 8 * y0 * L * d + 4 &
!      * d ** 2) / (4 * L ** 2 * x0 ** 2 + 4 * L ** 2 * x0 + L ** 2 + 4  &
!     * y0 ** 2 * L ** 2 + 8 * y0 * L * d + 4 * d ** 2) ** 2 * (8 * y0 * &
!      L ** 2 + 8 * L * d)
!
!d2A0dy0dx0 = -32 * c1 * L ** 4 * x0 / (4 * L ** 2 * x0 ** 2 - 4 * L ** &
!      2 * x0 + L ** 2 + 4 * y0 ** 2 * L ** 2 + 8 * y0 * L * d + 4 * d ** &
!      2) / (4 * L ** 2 * x0 ** 2 + 4 * L ** 2 * x0 + L ** 2 + 4 * y0 **&
!      2 * L ** 2 + 8 * y0 * L * d + 4 * d ** 2) - 4 * c1 * L ** 2 * (- &
!     4 * L ** 2 * x0 ** 2 + L ** 2 + 4 * y0 ** 2 * L ** 2 + 8 * y0 * L  &
!     * d + 4 * d ** 2) / (4 * L ** 2 * x0 ** 2 - 4 * L ** 2 * x0 + L ** &
!      2 + 4 * y0 ** 2 * L ** 2 + 8 * y0 * L * d + 4 * d ** 2) ** 2 / (4 &
!      * L ** 2 * x0 ** 2 + 4 * L ** 2 * x0 + L ** 2 + 4 * y0 ** 2 * L **&
!      2 + 8 * y0 * L * d + 4 * d ** 2) * (8 * L ** 2 * x0 - 4 * L ** 2 &
!     ) - 4 * c1 * L ** 2 * (-4 * L ** 2 * x0 ** 2 + L ** 2 + 4 * y0 **  &
!     2 * L ** 2 + 8 * y0 * L * d + 4 * d ** 2) / (4 * L ** 2 * x0 ** 2  &
!     - 4 * L ** 2 * x0 + L ** 2 + 4 * y0 ** 2 * L ** 2 + 8 * y0 * L * d &
!      + 4 * d ** 2) / (4 * L ** 2 * x0 ** 2 + 4 * L ** 2 * x0 + L ** 2  &
!     + 4 * y0 ** 2 * L ** 2 + 8 * y0 * L * d + 4 * d ** 2) ** 2 * (8 *  &
!     L ** 2 * x0 + 4 * L ** 2)
!
!
!END SUBROUTINE ddA0

END MODULE CMT_fields

