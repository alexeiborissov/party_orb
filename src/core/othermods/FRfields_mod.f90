MODULE FR_fields
  
  USE global
  USE M_products

  IMPLICIT NONE
  
  PRIVATE
  PUBLIC :: FRFIELDS

  CONTAINS 
!------------------------------------------------------------------------------!    
SUBROUTINE FRFIELDS(R,T,E,B,DBDX,DBDY,DBDZ,DBDT,DEDX,DEDY,DEDZ,DEDT,Vf)
! FR fields - attempt to use setup of Hesse et al, ApJ (2005) for kinematic flux rope eruption experiment

 REAL(num), DIMENSION(3),INTENT(OUT)	:: B,E
 REAL(num), DIMENSION(3),INTENT(OUT)	:: DBDX,DBDY,DBDZ,DBDT,DEDX,DEDY,DEDZ,DEDT
 REAL(num), DIMENSION(3),INTENT(IN)	:: R
 REAL(num), DIMENSION(3)  		:: Vf, dBedx,dBedy,dBedz
 REAL(num), INTENT(IN) 			:: T
 REAL(num)				:: opzsq, opxsq, omzsq, epsilon, depsilondt
  

! first define some simple quantities:
 opzsq=(1.0_num+R(3)*R(3)/Lz/Lz)
 omzsq=(1.0_num-R(3)*R(3)/Lz/Lz)
 opxsq=(1.0_num+R(1)*R(1)/Lx/Lx)
 
 epsilon=T/tau_n
 depsilondt=1.0_num/tau	!<- could replace with 1/tau_n and divide E, DEDXYZ by tscl 
 !depsilondtbar=tscl/tau
  
  !full B=Be+B1
 B(1)=B0*0.2_num
 B(2)=-B0*(epsilon*omzsq/opzsq/opzsq/opxsq+1.0_num)
 B(3)=B0*R(2)/Ly
 
 ! non-dimensionalised
 B=B/Bscl

 ! no velocity?
 Vf(1)=0.0_num 
 Vf(2)=0.0_num
 Vf(3)=0.0_num

 ! Electric fields
 !E(1)=0.0_num
 E(1)=B0*lscl*R(3)/opzsq/opxsq*depsilondt	! assume depsilon/dt=1
 !E(1)=-R(3)/opzsq/omzsp	! assume depsilon/dt=1
 E(2)=0.0_num
 E(3)=0.0_num
 
 !STOP
! non-dimensionalised 
 E=E/Escl

 !---------DERIVATIVES-------!
 ! calculate the spatial derivatives of individual components - in this case, only the x derivatives matter, but make sure only the z component is assigned..
 dBdx(1)=0.0_num
 dBdx(2)=2.0_num*B0*epsilon*omzsq/opzsq/opzsq/opxsq*R(1)/Lx/Lx
 dBdx(3)=0.0_num
 dBdy(1)=0.0_num
 dBdy(2)=0.0_num
 dBdy(3)=B0/Ly
 dBdz(1)=0.0_num
 dBdz(2)=2.0_num*B0*epsilon*R(3)/opzsq/opzsq/opzsq/opxsq/Lz/Lz*(3.0_num-R(2)*R(2)/Lz/Lz)
 dBdz(3)=0.0_num
 
 ! not sure quite how to deal with time dependence in here yet. Or normalisation.
 !dBdx = 0.0_num
 !dBdy = 0.0_num
 !dBdz = 0.0_num
 DBDX = dBdx/Bscl		! lengths normalised, fields NOT 
 DBDY = dBdy/Bscl
 DBDZ = dBdz/Bscl

 !dBdt= Bf/tau
 dBdt(1)= 0.0_num
 !dBdt(2)= 0.0_num
 dBdt(2)= -depsilondt*B0*omzsq/opzsq/opzsq/opxsq
 dBdt(3)= 0.0_num		
! non-dimensionalised
 DBDT= dBdt*Tscl/Bscl

! E is purely in z, therefore derivatives only vary in z too
 !dEdx(1)=0.0_num
 dEdx(1)=-2.0_num*B0*R(3)*R(1)/Lx/Lx/opzsq/opxsq/opxsq*depsilondt
 dEdx(2)=0.0_num
 dEdx(3)=0.0_num
 DEDX=DEDX/Escl*lscl			! again, all lengths in problem ALREADY normalised, so 1/Escl
 
 dEdy(1)=0.0_num
 dEdy(2)=0.0_num
 dEdy(3)=0.0_num
 DEDY=DEDY/Escl*lscl

 !dEdz(1)=0.0_num
 dEdz(1)=B0*omzsq/opzsq/opzsq/opxsq*depsilondt
 dEdz(2)=0.0_num
 dEdz(3)=0.0_num
 DEDZ=dEdz/Escl*lscl

! no time dependence
 dEdt=0.0_num
 DEDT= dEdt*Tscl/Escl

END SUBROUTINE FRFIELDS

END MODULE FR_fields

