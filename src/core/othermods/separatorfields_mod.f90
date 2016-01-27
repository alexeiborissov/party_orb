MODULE sep_fields
  
  USE global
  USE M_products

  IMPLICIT NONE
  
  PRIVATE
  PUBLIC :: SEPFIELDS

  CONTAINS 
!------------------------------------------------------------------------------!    
SUBROUTINE SEPFIELDS(R,T,E,B,DBDX,DBDY,DBDZ,DBDT,DEDX,DEDY,DEDZ,DEDT,Vf)
! re-written by JT to use other fields
! These are the fields specified by Wilmot-Smith and Hornig (2011)


 REAL(num), DIMENSION(3),INTENT(OUT)	:: B,E
 REAL(num), DIMENSION(3),INTENT(OUT)	:: DBDX,DBDY,DBDZ,DBDT,DEDX,DEDY,DEDZ,DEDT
 REAL(num), DIMENSION(3),INTENT(IN)	:: R
 REAL(num), DIMENSION(3)  		:: Be, Bf, Vf, dBedx,dBedy,dBedz,dBfdx,dBfdy,dBfdz
 REAL(num), INTENT(IN) 			:: T!, T1, T2
 REAL(num)				:: expfring, oan=L/a, oln=L/ll
  
 ! Be is the equilibrium B-field
 !Be(1)=b0*oL*oL*R(1)*(R(3)-3.0_num*z0)
 !Be(2)=b0*oL*oL*R(2)*(R(3)+3.0_num*z0) 
 !Be(3)=0.5_num*b0*oL*oL*(z0*z0-R(3)*R(3)+R(1)*R(1)+R(2)*R(2))
 Be(1)=b0*R(1)*(R(3)-3.0_num*z0*oL)
 Be(2)=b0*R(2)*(R(3)+3.0_num*z0*oL) 
 Be(3)=0.5_num*b0*(z0*z0*oL*oL-R(3)*R(3)+R(1)*R(1)+R(2)*R(2))
 ! now the general form of the flux ring
 
 !fring=b1*a*exp(-(R(1)-xc)*(R(1)-xc)*oan*oan-(R(2)-yc)*(R(2)-yc)*oan*oan-(R(3)-xc)*(R(3)-xc)*oln*oln)
 expfring=exp(-(R(1)-xc)*(R(1)-xc)*oan*oan-(R(2)-yc)*(R(2)-yc)*oan*oan-(R(3)-xc)*(R(3)-xc)*oln*oln)
 
 !Components of the perturbed field caused by flux ring
 Bf(1)=-2.0_num*b1*(R(2)-yc)*oan*expfring
 Bf(2)=2.0_num*b1*(R(1)-xc)*oan*expfring
 Bf(3)=0.0_num
 
 !full B=Be+B1
 !B(1)=Be(1)+(T-T1)/(T2-T1)*Bf(1)			!T/tau needs fixing here.. made Tau=T2
 !B(2)=Be(2)+(T-T1)/(T2-T1)*Bf(2)
 !B(3)=Be(3)+(T-T1)/(T2-T1)*Bf(3)
 B(1)=Be(1)+(T-myT1)/tau_n*Bf(1)			!not convinced this is the correct time dependent form...
 B(2)=Be(2)+(T-myT1)/tau_n*Bf(2)
 B(3)=Be(3)+(T-myT1)/tau_n*Bf(3)
 !print*, T, T1, tau
 
 ! non-dimensionalised
 B=B/Bscl

 ! no velocity?
 Vf(1)=0.0_num 
 Vf(2)=0.0_num
 Vf(3)=0.0_num

 ! Electric fields
 E(1)=0.0_num
 E(2)=0.0_num
 E(3)=-b1*a*expfring/tau
 !E(3)=-fring/Tscl
 !E(3)=0.0_num
 
 !print*, B1*a
 !print*, E
 !print*, Escl
 !print*, E/Escl
 !STOP
! non-dimensionalised 
 E=E/Escl

 !---------DERIVATIVES-------!
 ! calculate the spatial derivatives of individual components..
 dBedx(1)=b0*(R(3)-3.0_num*z0*oL)
 dBedx(2)=0.0_num
 dBedx(3)=b0*R(1)
 
 dBedy(1)=0.0_num
 dBedy(2)=b0*(R(3)+3.0_num*z0*oL)
 dBedy(3)=b0*R(2)
 
 dBedz(1)=b0*R(1)
 dBedz(2)=b0*R(2)
 dBedz(3)=-b0*R(3)

 !dBfdx(1)=Bf(1)*(-2.0_num)*(R(1)-xc)*oa*oa				!original - wrong (a unnormalised)
 !dBfdx(2)=Bf(2)*(-2.0_num)*(R(1)-xc)*oa*oa+2.0_num*fring*oa*oa
 !dBfdx(3)=0.0_num 
 !dBfdx(1)=-2.0_num*R(1)*oan*oan*Bf(1)					! written in terms of original B
 !dBfdx(2)=Bf(2)/R(1)-2.0_num*R(1)*oan*oan*Bf(2)
 !dBfdx(3)=0.0_num
 dBfdx(1)=4.0_num*b1*(R(2)-yc)*(R(1)-xc)*oan*oan*oan*expfring
 dBfdx(2)=(2.0_num*b1*oan-4.0_num*b1*(R(1)-xc)*(R(1)-xc)*oan*oan*oan)*expfring
 dBfdx(3)=0.0_num

 !dBfdy(1)=Bf(1)*(-2.0_num)*(R(2)-yc)*oa*oa-2.0_num*fring*oa*oa
 !dBfdy(2)=Bf(2)*(-2.0_num)*(R(2)-yc)*oa*oa
 !dBfdy(3)=0.0_num
 dBfdy(1)=(-2.0_num*b1*oan+4.0_num*b1*(R(2)-yc)*(R(2)-yc)*oan*oan*oan)*expfring
 dBfdy(2)=-4.0_num*b1*(R(2)-yc)*(R(1)-xc)*oan*oan*oan*expfring
 dBfdy(3)=0.0_num
 
 !dBfdz(1)=Bf(1)*2.0_num*(R(3)-zc)*oa*oa
 !dBfdz(2)=Bf(2)*2.0_num*(R(3)-zc)*oa*oa		! are these really oa's?
 !dBfdz(3)=0.0_num
 dBfdz(1)=4.0_num*b1*(R(3)-zc)*(R(2)-yc)*oan*oln*oln
 dBfdz(2)=-4.0_num*b1*(R(3)-zc)*(R(1)-xc)*oan*oln*oln
 dBfdz(3)=0.0_num
 
 ! not sure quite how to deal with time dependence in here yet. Or normalisation.
 dBdx = dBedx+(T-myT1)/tau_n*dBfdx
 dBdy = dBedy+(T-myT1)/tau_n*dBfdy
 dBdz = dBedz+(T-myT1)/tau_n*dBfdz
 !dBdx = dBedx+(T-T1)/(T2-T1)*dBfdx
 !dBdy = dBedy+(T-T1)/(T2-T1)*dBfdy
 !dBdz = dBedz+(T-T1)/(T2-T1)*dBfdz
! non-dimensionalised
! DBDX = dBdx*Lscl/Bscl
! DBDY = dBdy*Lscl/Bscl
! DBDZ = dBdz*Lscl/Bscl
 DBDX = dBdx/Bscl		! lengths normalised, fields NOT 
 DBDY = dBdy/Bscl
 DBDZ = dBdz/Bscl

 
 !dBdt= Bf/tau
 dBdt= 0.0_num		
! non-dimensionalised
 DBDT= dBdt*Tscl/Bscl

! E is purely in z, therefore derivatives only vary in z too
 dEdx(1)=0.0_num
 dEdx(2)=0.0_num
 !dEdx(3)=0.0_num
 dEdx(3)=2.0_num*(R(1)-xc)*oan*oan*(b1*a/tau)*expfring
 !DEDX=DEDX*Lscl/Escl
 DEDX=DEDX/Escl			! again, all lengths in problem ALREADY normalised, so 1/Escl
 
 dEdy(1)=0.0_num
 dEdy(2)=0.0_num
 !dEdy(3)=0.0_num
 dEdy(3)=2.0_num*(R(2)-yc)*oan*oan*(b1*a/tau)*expfring
 !DEDY=DEDY*Lscl/Escl
 DEDY=DEDY/Escl

 dEdz(1)=0.0_num
 dEdz(2)=0.0_num
 !dEdz(3)=0.0_num
 dEdz(3)=2.0_num*(R(3)-zc)*oln*oln*(b1*a/tau)*expfring
 !DEDZ=dEdz*Lscl/Escl
 DEDZ=dEdz/Escl

! no time dependence
 dEdt(1)=0.0_num
 dEdt(2)=0.0_num
 dEdt(3)=0.0_num
 DEDT= dEdt*Tscl/Escl

END SUBROUTINE SEPFIELDS

END MODULE sep_fields

