MODULE test_fields
  
  USE global
  USE M_products

  IMPLICIT NONE
  
  PRIVATE
  PUBLIC :: TESTFIELDS

  CONTAINS 
!------------------------------------------------------------------------------!    
SUBROUTINE TESTFIELDS(R,T,E,B,DBDX,DBDY,DBDZ,DBDT,DEDX,DEDY,DEDZ,DEDT,Vf)
! testfields - routines to test how the relativistic and nonrelativistic versions conform.

 REAL(num), DIMENSION(3),INTENT(OUT)	:: B,E
 REAL(num), DIMENSION(3),INTENT(OUT)	:: DBDX,DBDY,DBDZ,DBDT,DEDX,DEDY,DEDZ,DEDT
 REAL(num), DIMENSION(3),INTENT(IN)	:: R
 REAL(num), DIMENSION(3)  		:: Vf, dBedx,dBedy,dBedz
 REAL(num), INTENT(IN) 			:: T
 REAL(num)				:: ntemp, ontemp 
  

 ntemp=64.0_num
 ontemp=1.0_num/ntemp 
  
  !full B=Be+B1
 B(1)=0.0_num
 B(2)=0.0_num
 B(3)=B0
 !B(3)=B0*exp(R(1))

 !print*, T, T1, tau
 
 ! non-dimensionalised
 B=B/Bscl

 ! no velocity?
 Vf(1)=0.0_num 
 Vf(2)=0.0_num
 Vf(3)=0.0_num

 ! Electric fields
 E(1)=0.0_num
 !E(1)=E0
 E(2)=0.0_num
 E(3)=0.0_num
 
 !print*, B1*a
 !print*, E
 !print*, Escl
 !print*, E/Escl
 !STOP
! non-dimensionalised 
 E=E/Escl

 !---------DERIVATIVES-------!
 ! calculate the spatial derivatives of individual components - in this case, only the x derivatives matter, but make sure only the z component is assigned..
 dBedx(1)=0.0_num
 dBedx(2)=0.0_num
 dBedx(3)=0.0_num
 !dBedx(3)=b0*exp(R(1))
 
 ! not sure quite how to deal with time dependence in here yet. Or normalisation.
 dBdx = dBedx
 dBdy = 0.0_num
 dBdz = 0.0_num
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
 dEdx(3)=0.0_num
 DEDX=DEDX/Escl			! again, all lengths in problem ALREADY normalised, so 1/Escl
 
 dEdy(1)=0.0_num
 dEdy(2)=0.0_num
 dEdy(3)=0.0_num
 DEDY=DEDY/Escl

 dEdz(1)=0.0_num
 dEdz(2)=0.0_num
 dEdz(3)=0.0_num
 DEDZ=dEdz/Escl

! no time dependence
 dEdt(1)=0.0_num
 dEdt(2)=0.0_num
 dEdt(3)=0.0_num
 DEDT= dEdt*Tscl/Escl

END SUBROUTINE TESTFIELDS

END MODULE test_fields

