PROGRAM reltest		! Relativistic Particle Code: Main Program (JT Dec 2015)
!+ sets up initial field conditions, particle grid and hands over to rkdriver

USE GLOBAL
USE lare_functions
USE mpi_routines
USE M_DRIVERR, ONLY: RKDRIVE
USE bourdin_fields, ONLY: bour_ini, bour_fini
USE M_products, ONLY: DOT, CROSS
USE M_fields, ONLY: FIELDS

IMPLICIT NONE

 INTEGER :: NOK, NBAD
 INTEGER :: NKEEP,time_no
 INTEGER :: pos_no_x,pos_no_y,pos_no_z,pos_no_alpha,pos_no_ekin
 INTEGER :: pnmax
 INTEGER, DIMENSION(3) :: pos_no_r
 
 REAL(num), DIMENSION(3) :: gds, lbox
 REAL(num), DIMENSION(3) :: RSTART,RSTARTKEEP!, R1,R2
 REAL(num), DIMENSION(NKEEPMAX) :: TT 
 REAL(num), DIMENSION(NKEEPMAX,3) :: S, TOTAL
 REAL(num), DIMENSION(3) :: Et,Bt,DBDXt,DBDYt,DBDZt,DBDTt,DEDXt,DEDYt,DEDZt,DEDTt,Vft
 CHARACTER(LEN=35)	 :: finfile
 
 !welcome screen
 WRITE(*,*) "====RELATIVISTIC particle code====="
 WRITE(*,*) "using ", FMOD, " fields module."
 bourdinflag=.FALSE.
 l3dflag=.FALSE.
 l2dflag=.FALSE.
 analyticalflag=.FALSE.

 CALL read_param
  
 !The initial setup and comparison of the grid and environment depends on which module has been chosen.
  
  IF ((str_cmp(FMOD, "L3D")).OR.(str_cmp(FMOD, "l3d"))) THEN
   ndims=3
   allocate(dims(ndims))
   l3dflag=.TRUE.
   CALL MPI_INIT(errcode)
   CALL mpi_initialise      ! mpi_routines.f90 
   IF (evenlarefield) THEN
    PRINT*, '..evaluating particle array against even lare grid..'
    IF (((R1(1)/lscl).le.(-1.0_num*x_end)).OR.((R2(1)/lscl).ge.(x_end)))  THEN
      WRITE(*,*) '..particles not within lare x extent '
      gflag=.true.
    ENDIF
    IF (((R1(2)/lscl).le.(-1.0_num*y_end)).OR.((R2(2)/lscl).ge.(y_end))) THEN
      WRITE(*,*) '..particles not within lare y extent '
      gflag=.true.
    ENDIF
    IF (((R1(3)/lscl).le.(-1.0_num*z_end)).OR.((R2(3)/lscl).ge.(z_end))) THEN
      WRITE(*,*) '..particles not within lare z extent '
      gflag=.true.
    ENDIF
    IF (gflag) THEN
      WRITE(*,*) 'terminating: particle grid out of bounds set in Lare.'
      STOP
    ELSE 
      WRITE(*,*) 'fine!'
    ENDIF 
   ELSE
    PRINT*, '..evaluating particle array against ODD lare grid..'
    IF (((R1(1)/lscl).le.xe(1)).OR.((R2(1)/lscl).ge.xe(2)))  THEN
      WRITE(*,*) '..particles not within lare x extent '
      gflag=.true.
    ENDIF
    IF (((R1(2)/lscl).le.ye(1)).OR.((R2(2)/lscl).ge.ye(2)))  THEN
      WRITE(*,*) '..particles not within lare y extent '
      gflag=.true.
    ENDIF
    IF (((R1(3)/lscl).le.ze(1)).OR.((R2(3)/lscl).ge.ze(2)))  THEN
      WRITE(*,*) '..particles not within lare z extent '
      gflag=.true.
    ENDIF
    IF (gflag) THEN
      WRITE(*,*) 'terminating: particle grid out of bounds set in Lare.'
      STOP
    ELSE 
      WRITE(*,*) 'fine!'
    ENDIF 
   ENDIF
  ELSE IF ((str_cmp(FMOD, "L2D")).OR.(str_cmp(FMOD, "l2d"))) THEN
   l2dflag=.TRUE.
   ndims=2
   c_ndims=ndims
   
   allocate(dims(ndims))
   CALL MPI_INIT(errcode)
   CALL mpi_initialise_2d
   IF ((R1(3).NE.R2(3)).OR.(RSTEPS(3).GT.1)) THEN
    PRINT*, 'Using Lare2d data - z position must be single valued!!'
    STOP
   ENDIF
   IF (evenlarefield) THEN
    PRINT*, '..evaluating particle array against even lare grid..'
    IF (((R1(1)/lscl).le.(-1.0_num*x_end)).OR.((R2(1)/lscl).ge.(x_end)))  THEN
      WRITE(*,*) '..particles not within lare x extent '
      gflag=.true.
    ENDIF
    IF (((R1(2)/lscl).le.(-1.0_num*y_end)).OR.((R2(2)/lscl).ge.(y_end))) THEN
      WRITE(*,*) '..particles not within lare y extent '
      gflag=.true.
    ENDIF
    IF (gflag) THEN
      WRITE(*,*) 'terminating: particle grid out of bounds set in Lare.'
      STOP
    ELSE 
      WRITE(*,*) 'fine!'
    ENDIF 
   ELSE
    PRINT*, '..evaluating particle array against ODD lare grid..'
    IF (((R1(1)/lscl).le.xe(1)).OR.((R2(1)/lscl).ge.xe(2)))  THEN
      WRITE(*,*) '..particles not within lare x extent '
      gflag=.true.
    ENDIF
    IF (((R1(2)/lscl).le.ye(1)).OR.((R2(2)/lscl).ge.ye(2)))  THEN
      WRITE(*,*) '..particles not within lare y extent '
      gflag=.true.
    ENDIF
    IF (gflag) THEN
      WRITE(*,*) 'terminating: particle grid out of bounds set in Lare.'
      STOP
    ELSE 
      WRITE(*,*) 'fine!'
    ENDIF 
   ENDIF
  ELSE IF ((str_cmp(FMOD, "SEP")).OR.(str_cmp(FMOD, "sep"))) THEN
    analyticalflag=.TRUE.
    PRINT*, '..evaluating particle array against analytical field bounds..'
    IF (((R1(1)/lscl).le.(-1.0_num*x_end)).OR.((R2(1)/lscl).ge.(x_end)))  THEN
      WRITE(*,*) '..particles not within lare x extent '
      gflag=.true.
    ENDIF
    IF (((R1(2)/lscl).le.(-1.0_num*y_end)).OR.((R2(2)/lscl).ge.(y_end))) THEN
      WRITE(*,*) '..particles not within lare y extent '
      gflag=.true.
    ENDIF
    IF (((R1(3)/lscl).le.(-1.0_num*z_end)).OR.((R2(3)/lscl).ge.(z_end))) THEN
      WRITE(*,*) '..particles not within lare z extent '
      gflag=.true.
    ENDIF
    IF (gflag) THEN
      WRITE(*,*) 'terminating: particle grid out of bounds.'
      STOP
    ELSE 
      WRITE(*,*) 'fine!'
    ENDIF  
  ELSE IF ((str_cmp(FMOD, "CMT")).OR.(str_cmp(FMOD, "cmt"))) THEN
      !CMT setup?
  ELSE IF ((str_cmp(FMOD, "TEST")).OR.(str_cmp(FMOD, "test"))) THEN
      !CMT setup?
  ELSE IF ((str_cmp(FMOD, "BOR")).OR.(str_cmp(FMOD, "bor"))) THEN
   bourdinflag=.TRUE.
   CALL bour_ini      ! read in data
   xee=(/myx(6),myx(nx-5)/)
   yee=(/myy(6),myy(ny-5)/)
   zee=(/myz(6),myz(nz-5)/)
   print*, '----'
   PRINT*, '..evaluating particle array against BOURDIN grid..'
   IF ((R1(1)/lscl.le.xee(1)).OR.(R2(1)/lscl.ge.xee(2)))  THEN
    WRITE(*,*) '..particles not within x extent '
    gflag=.true.
   ENDIF
    IF ((R1(2)/lscl.le.yee(1)).OR.(R2(2)/lscl.ge.yee(2)))  THEN
    WRITE(*,*) '..particles not within y extent '
    gflag=.true.
   ENDIF
   IF ((R1(3)/lscl.le.zee(1)).OR.(R2(3)/lscl.ge.zee(2)))  THEN
    WRITE(*,*) '..particles not within z extent '
    gflag=.true.
   ENDIF
   IF (gflag) THEN
    WRITE(*,*) 'terminating: particle grid out of bounds of BOURDIN grid.'
    STOP
   ELSE 
    WRITE(*,*) 'fine!'
   ENDIF 
  ELSE
   PRINT*, "incorrect module selection, choose from:"
   PRINT*, "['lare','sep','CMT','test','bour']"
   STOP
  END IF
   
! STOP
   
  DO i=1,3
   IF (RSTEPS(i).EQ.1) THEN 	;! if rsteps=1, 1/(rsteps-1)=1/0!!
    gds(i)=1.0_num 
   ELSE 
    gds(i)=1.0_num/REAL(RSTEPS(i)-1)
   ENDIF
  ENDDO
  lbox=(/R2(1)-R1(1),R2(2)-R1(2),R2(3)-R1(3)/)

  nparticles=RSTEPS(1)*RSTEPS(2)*RSTEPS(3)*(AlphaSteps-1)*EkinSteps
  dalpha = (AlphaMax-Alphamin)/(Alphasteps - 1.0d0) !added by S.Oskoui
  
  !Adjust T2 to use loop value.
  !T2=time_no*1.0_num

  T1Keep=T1
  T2Keep=T2
  
  writervs=1
 
  IF (JTo4.eq.1) WRITE(finfile,"(A,'finishr.tmp')"),dlocR
  IF (JTo4.eq.1)  open(49,file=finfile,recl=1024,status='unknown')
  PRINT*, ''
  ! restart our calculation for certain particles within a given grid?
  ! can use to divvy up the same grid to different CPUs.
  IF (p_restart) THEN	! p_restart_no is the place a calculation starts again at
   PRINT*, '--RESTARTING p_grid--'
   pn=p_restart_no-1
   pos_no_x = pn / (RSTEPS(2)*RSTEPS(3))  
   pos_no_y = MOD(pn / RSTEPS(3), RSTEPS(2))
   pos_no_z = MOD(pn,RSTEPS(3))
  ELSE
   PRINT*, '--starting particle grid from beginning--'   
   pn=0
   pos_no_x=0
   pos_no_y=0
   pos_no_z=0
  ENDIF
  
  IF (p_stop) THEN	!p_stop_no is the particle number where we stop calculating
   pnmax=p_stop_no-1
  ELSE			! otherwise iterate to total number of expected particles
   pnmax=nparticles
  ENDIF

  
  DO WHILE (pos_no_x .LE. RSTEPS(1)-1)
   DO WHILE (pos_no_y .LE. RSTEPS(2)-1)
    DO WHILE ((pos_no_z .LE. RSTEPS(3)-1).AND.(pn .LE. pnmax))
     DO pos_no_alpha = 2, AlphaSteps,1
      DO pos_no_ekin = 1, EkinSteps,1

       pos_no_r = (/pos_no_x,pos_no_y, pos_no_z/)
       print*, tempr
       IF (RANDOMISE_R) THEN
        CALL init_random_seed()
        CALL RANDOM_NUMBER(tempr)
        RSTART   = R1+tempr*lbox	!randomise position in bounds set in input
       ELSE
        RSTART   = R1+lbox*(pos_no_r*1.0_num)*gds
       ENDIF
       
       pn= pn + 1
       
       !call progress(pn,nparticles) ! generate the progress bar.
       
       IF (JTo4.eq.1) write(49,"(I4)",advance='no'), pn	   
	   
       IF (nparticles.gt.1000) THEN 
        print 1000, pn,nparticles, RSTART     
        1000 format ("particle no. ",i4,"/",i4, ", R=(",ES9.2,",",ES9.2,",",ES9.2,")")
       ELSE 
        print 1001, pn,nparticles, RSTART     
        1001 format ("particle no. ",i3,"/",i3, ", R=(",ES9.2,",",ES9.2,",",ES9.2,")")
       ENDIF    

       T1=T1Keep
       T2=T2Keep

       IF (RANDOMISE_A) THEN
        alpha = Alphamin+dalpha*(pos_no_alpha -1)*tempa	!added by S.Oskoui
       ELSE
        alpha = Alphamin+dalpha*(pos_no_alpha -1)	!added by S.Oskoui
       ENDIF
       alpha = alpha*Pi/180.0d0				! RADEG: added by S.Oskoui

       
       !pos_no_ekin starts from 0, if started from 1 then (stepekin-1)
       IF (RANDOMISE_E) THEN
        Ekin=EKinLow+(EKinHigh-EKinLow)*pos_no_ekin/(EkinSteps*1.0d0)*tempe 
       ELSE
        Ekin=EKinLow+(EKinHigh-EKinLow)*pos_no_ekin/(EkinSteps*1.0d0)
       ENDIF
       !alpha = pi/(no of steps+1) if fullangle is 1 (ie, steps from >=0 to >Pi (but not including Pi))
       !alpha = pi/2/(no of steps) if fullangle is 0 (steps from 0 to Pi/2 inclusive)
  
       RSTARTKEEP=RSTART
       USTARTKEEP=USTART
       GAMMASTARTKEEP=GAMMASTART
       
       !PRINT*,'Normalising:'
       RSTART=RSTART/Lscl
       RSTARTKEEP=RSTARTKEEP/Lscl
       T1=T1/Tscl
       T2=T2/Tscl
       
       ! WARNING passing in dimensional Ekin into mu calc
       CALL JTMUcalc(MU,USTART,GAMMASTART,Ekin,Alpha,RSTART,T1,T2)
              
       Ekin = Ekin*AQ/Ekscl        
       
       !Call the rk sophisticated driver, which then works out the arrays for the
       !time steps and positions.
       CALL RKDRIVE(RSTART,USTART,GAMMASTART,MU,T1,T2,EPS,H1,NOK,NBAD,TT,S,TOTAL)
        
       NKEEP = (NOK +NBAD)/NSTORE

      ! CALL WRITE_ENDTIME(RSTART,T2,MU,VPARSTART)
      
      END DO
     END DO
      pos_no_z=pos_no_z+1
    END DO
    pos_no_y=pos_no_y+1
    pos_no_z=0
   END DO
   pos_no_x=pos_no_x+1
   pos_no_y=0
  END DO
 CLOSE(49)
 !CALL MAKEFILE(time_no)
  
 !END DO
 
 IF ((str_cmp(FMOD, "LARE")).OR.(str_cmp(FMOD, "lare"))) THEN	!forget arrays at end
  CALL mpi_close                     ! mpi_routines.f90
  CALL MPI_FINALIZE(errcode)
 ENDIF

 IF ((str_cmp(FMOD, "BOUR")).OR.(str_cmp(FMOD, "bour"))) THEN
  CALL bour_fini      		! deallocate stuff, leave everything nice and tidy
 ENDIF 

!------------------------------------------------------------------------------!
 Contains
!------------------------------------------------------------------------------!
SUBROUTINE JTMUcalc(mu,USTART,GAMMASTART, Ekin,alpha,RSTART,T1,T2)

  REAL(num), DIMENSION(3),INTENT(IN) :: RSTART
  REAL(num), INTENT(IN) :: T1,T2, Ekin, Alpha				! 
  REAL(num), INTENT(OUT) :: mu, gammastart, Ustart
  REAL(num), DIMENSION(3) :: B,El,a2,a3,a4,a5,a6,a7,a8,a9,a10,ue, RT
  !REAL(num) :: magB,vtot,vperp, vparstart!,Erest
  REAL(num) :: modB,vtot, gamma
 
 !calculate B, E, V at this point/time:
 CALL FIELDS(RSTART,T1,El,B,a2,a3,a4,a5,a6,a7,a8,a9,a10,T1,T2)

 !calculate magnitude of B
 modB=B(1)*B(1)+B(2)*B(2)+B(3)*B(3)
 modB=sqrt(modB)
 RT=RSTART

 !print*, 'RT=', RT
 !PRINT*, 'B=', B
 !STOP
 !print*, modB

 !Erest = (M*c*c)*oneuponAQ
 
 ! E X B drift
 ue=cross(El,B)/dot(B,B)  !*0.5
 
 gamma=Ekin/m/c/c*AQ+1.0d0
 
 ! vtot^2=vpar^2+vperp^2+UE^2!!, UE should NOT come out of Gamma!
 !vtot=c/gamma*sqrt(gamma*gamma-1)			! I *think* vtot is dimensional
 !vtot=sqrt(c*c/gamma/gamma*(gamma*gamma-1)-vscl*vscl*dot(ue,ue))
 vtot=sqrt(c*c-c/gamma*c/gamma-vscl*vscl*dot(ue,ue))			! THRELFALL ET AL 2015
 ! vtot=sqrt(c*c-c/gamma*c/gamma-vce*vce*dot(ue,ue))			! THRELFALL ET AL 2015
 
 !Ustart=(c*sqrt(gamma*gamma-1))*cos(alpha)
 Ustart=gamma*vtot*cos(alpha)

 
 mu=0.5_num*m*vtot*vtot*sin(alpha)*sin(alpha)*gamma*gamma/modB		
 
 USTART=Ustart/vscl					! hence, non-dimensionalising..
 GAMMASTART=gamma
 mu=mu/m/vscl/vscl					! no bscl normalising factor - using normalised B's already!
 !PRINT*, 'vtot=', vtot
 !PRINT*, 'gamma=', gamma
 !PRINT*, 'mu=', mu
 !PRINT*, 'modB=', modB
  !calculate mu
  !mu=0.5_num*vperp*vperp/magB*gammastart*gammastart
 
 ! WRITE (19,*) RStart,T1,Ekin,Alpha, mu*magB, 0.5_num*vparstart*vparstart
  !WRITE (19,*) vtot,vperp,vparstart,El,B,magB,mu

END SUBROUTINE

END PROGRAM reltest