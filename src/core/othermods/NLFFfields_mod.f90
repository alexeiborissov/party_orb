MODULE NLFF_fields
  
  USE global
  USE M_products
  USE lare_functions

  IMPLICIT NONE

  PRIVATE
  PUBLIC :: NLFFF_ini, NLFFF_fini, NLFFFIELDS

  CONTAINS 
!------------------------------------------------------------------------------!
SUBROUTINE NLFFF_ini
 
 INTEGER				:: ii
 CHARACTER(len=21)			:: nom='ar11437_run20_256'
 CHARACTER(len=21)			:: nomB='_ftn2idlB',nomJ='_ftn2idlJ'
 CHARACTER(len=21)			:: nomG='_ftn2idlgrid' 
 REAL(num), DIMENSION(:), ALLOCATABLE  	:: griddat
 REAL(sp), DIMENSION(:, :, :), ALLOCATABLE  :: mydat
   
 CHARACTER(LEN = 20+data_dir_max_length):: locB, locJ, gridloc
 LOGICAL				:: ONLY=.FALSE.

 gridloc=trim(adjustl(sloc))//trim(nom)//trim(nomG)//filetype3
 locB=trim(adjustl(sloc))//trim(nom)//trim(nomB)//filetype3
 locJ=trim(adjustl(sloc))//trim(nom)//trim(nomJ)//filetype3
 
 nx=nx_global
 ny=ny_global
 nz=nz_global
 
 ALLOCATE(griddat(1:nx+1))	
 ALLOCATE(myx(1:nx+1))
 ALLOCATE(myy(1:ny+1))
 ALLOCATE(myz(1:nz+1))
 
 !read grid
 PRINT*, 'reading grid'
 OPEN(24,file=gridloc, access='stream', status='old')
 READ(24) griddat
 myx(1:nx+1)=griddat(1:nx+1)
 READ(24) griddat
 myy(1:ny+1)=griddat(1:nx+1)
 READ(24) griddat
 myz(1:nz+1)=griddat(1:nx+1)
 CLOSE(24)

 DEALLOCATE(griddat)
 !WRITE(*,'(a,f15.5,a,f15.5,a,e15.5,a)') "xr=[", minval(myx),"->", maxval(myx),"]*", lscl, "m"
 !WRITE(*,'(a,f15.5,a,f15.5,a,e15.5,a)') "yr=[", minval(myy),"->", maxval(myy),"]*", lscl, "m"
 !WRITE(*,'(a,f15.5,a,f15.5,a,e15.5,a)') "zr=[", minval(myz),"->", maxval(myz),"]*", lscl, "m"
 print*, minval(myx),"->", maxval(myx)
 print*, minval(myy),"->", maxval(myy)
 print*, minval(myz),"->", maxval(myz) 
 !NORMALISE NOW:
 myx=myx*1e6/lscl
 myy=myy*1e6/lscl
 myz=myz*1e6/lscl	!the values are already in megameters - *10^6 to be in meters, but divide by arbitrary new lengthscale
 

! read fields
 write (*,"(a)",advance="no") 'reading NLFF field variables: '
 ALLOCATE(bx(1:nx+1,1:ny+1,1:nz+1,1))
 ALLOCATE(mydat(1:nx+1,1:ny+1,1:nz+1))
 OPEN(25,file=locB, access='stream', status='old')
 READ(25) mydat
 !bx(1:nx+1,1:ny+1,1:nz+1,1)=mydat
 DO ii=1,nx+1
  bx(ii,1:ny+1,1:nz+1,1)=stagger_bx(DBLE(mydat(ii,1:ny+1,1:nz+1)))*1E-4
 ENDDO  
 ALLOCATE(by(1:nx+1,1:ny+1,1:nz+1,1))
 READ(25) mydat
 !by(1:nx+1,1:ny+1,1:nz+1,1)=mydat
 DO ii=1,ny+1
  by(1:nx,ii,1:nz+1,1)=stagger_by(DBLE(mydat(1:nx+1,ii,1:nz+1)))*1E-4
 ENDDO  
 ALLOCATE(bz(1:nx+1,1:ny+1,1:nz+1,1)) 
 READ(25) mydat
 !bz(1:nx+1,1:ny+1,1:nz+1,1)=mydat
 DO ii=1,nz
  bz(1:nx+1,1:ny+1,ii,1)=stagger_bz(DBLE(mydat(1:nx+1,1:ny+1,ii)))*1E-4
 ENDDO 
 CLOSE(25)
 
 !reading in j
 ALLOCATE(jx(1:nx+1,1:ny+1,1:nz+1,1))
 ALLOCATE(jy(1:nx+1,1:ny+1,1:nz+1,1))
 ALLOCATE(jz(1:nx+1,1:ny+1,1:nz+1,1)) 
 OPEN(26,file=locJ, access='stream', status='old')
 READ(26) mydat
 !jx(1:nx+1,1:ny+1,1:nz+1,1)=mydat
 DO ii=1,nx+1
  jx(ii,1:ny+1,1:nz+1,1)=stagger_j(DBLE(mydat(ii,1:ny+1,1:nz+1)))*1E-4*1e-6
 ENDDO 
 READ(26) mydat
 !jy(1:nx+1,1:ny+1,1:nz+1,1)=mydat
 DO ii=1,ny+1
  jy(1:nx,ii,1:nz+1,1)=stagger_j(DBLE(mydat(1:nx+1,ii,1:nz+1)))*1E-4*1e-6
 ENDDO 
 READ(26) mydat
 !jz(1:nx+1,1:ny+1,1:nz+1,1)=mydat
 DO ii=1,nz
  jz(1:nx+1,1:ny+1,ii,1)=stagger_j(DBLE(mydat(1:nx+1,1:ny+1,ii)))*1E-4*1e-6
 ENDDO 
 CLOSE(26)
 
 IF (FIELDDUMP) THEN
   ! quick way to check the fields read in are the same as those seen in Lare
    PRINT*, 'DUMPING MAGNETIC AND VELOCITY FIELD DATA:'
    OPEN(34, file=trim(adjustl(dloc))//'bx2idl.dat', form="unformatted")
    WRITE(34) bx(1:nx+1,1:ny+1,1:nz+1,1)
    CLOSE(34)
    OPEN(35, file=trim(adjustl(dloc))//'by2idl.dat', form="unformatted")
    WRITE(35) by(1:nx+1,1:ny+1,1:nz+1,1)
    CLOSE(35)
    OPEN(36, file=trim(adjustl(dloc))//'bz2idl.dat', form="unformatted")
    WRITE(36) bz(1:nx+1,1:ny+1,1:nz+1,1)
    CLOSE(36)
    OPEN(37, file=trim(adjustl(dloc))//'jx2idl.dat', form="unformatted")
    WRITE(37) jx(1:nx+1,1:ny+1,1:nz+1,1)
    CLOSE(37)
    OPEN(38, file=trim(adjustl(dloc))//'jy2idl.dat', form="unformatted")
    WRITE(38) jy(1:nx+1,1:ny+1,1:nz+1,1)
    CLOSE(38)
    OPEN(39, file=trim(adjustl(dloc))//'jz2idl.dat', form="unformatted")
    WRITE(39) jz(1:nx+1,1:ny+1,1:nz+1,1)
    CLOSE(39)
    PRINT*, 'DONE. TERMINATING.'  
    STOP
 ENDIF
 
 DEALLOCATE(mydat) 
 write (*,"(a)") 'done!' 
  !currently do not have E field derivs
 ! can easily include this, but for now calculate E-field derivs on the fly like for Lare3d..
 
END SUBROUTINE NLFFF_ini
!------------------------------------------------------------------------------!
SUBROUTINE NLFFF_fini

    DEALLOCATE(bx)
    DEALLOCATE(by)
    DEALLOCATE(bz)
    DEALLOCATE(Jx)
    DEALLOCATE(Jy)
    DEALLOCATE(Jz)
    DEALLOCATE(myx)
    DEALLOCATE(myy)
    DEALLOCATE(myz)

END SUBROUTINE NLFFF_fini
!------------------------------------------------------------------------------!
FUNCTION N3d(R,T)
! function to calculate local quantities at the particle position
! (slimmed down version of f3d with only local interpolation!).
! OCT2014: updated to include higher order derivs to increase accuracy

   REAL(num), DIMENSION(3), INTENT(IN)		:: R		!actual position
   REAL(num), INTENT(IN)			:: T
   REAL(num), DIMENSION(36)			:: N3d
   REAL(num)					:: temp,  modj
   REAL(num), DIMENSION(3)			:: dg, odg, coffset
   REAL(num), DIMENSION(:), ALLOCATABLE		:: dgt, odgt
   REAL(num), DIMENSION(:), ALLOCATABLE		:: bxt, byt, bzt,vxt, vyt, vzt, Ext, Eyt, Ezt, jxt, jyt, jzt
   REAL(num), DIMENSION(:), ALLOCATABLE		:: dbxdxt,dbxdyt,dbxdzt,dbydxt,dbydyt,dbydzt,dbzdxt,dbzdyt,dbzdzt
   REAL(num), DIMENSION(:), ALLOCATABLE		:: dExdxt,dExdyt,dExdzt,dEydxt,dEydyt,dEydzt,dEzdxt,dEzdyt,dEzdzt
   REAL(num), DIMENSION(:, :, :), ALLOCATABLE	:: mbx, mby, mbz, mEx, mEy, mEz, mjx, mjy, mjz, mvx, mvy, mvz
   REAL(num), DIMENSION(:, :, :), ALLOCATABLE	:: dmbxdx,dmbxdy,dmbxdz,dmbydx,dmbydy,dmbydz,dmbzdx,dmbzdy,dmbzdz
   REAL(num), DIMENSION(:, :, :), ALLOCATABLE	:: dmExdx,dmExdy,dmExdz,dmEydx,dmEydy,dmEydz,dmEzdx,dmEzdy,dmEzdz
   REAL(num), DIMENSION(:, :, :), ALLOCATABLE	:: meta
   INTEGER, DIMENSION(4) 			:: l		! set at value it could never reach!
   INTEGER					:: jjx, jjy, jjz,jjt, rpt
   LOGICAL					:: fxflag=.FALSE., fyflag=.FALSE., fzflag=.FALSE.
 
   temp=0.0_num
   dg=(/myx(2)-myx(1),myy(2)-myy(1),myz(2)-myz(1)/)	! grid spacing
   odg=(/1.0_num/dg(1),1.0_num/dg(2),1.0_num/dg(3)/)	! one over grid spacing
   l=(/-nx,-ny,-nz,-nframes/)					! initial value of l, set to silly value (as >=0 triggers flag)   
    
          
  ! --STEP ONE-- !
  ! first, locate the (x,y,z) position on the grid
  ! particle found between (jjx,jjy,jjz) and (jjx+1,jjy+1,jjz+1)   
     
   DO jjx=5,nx-5
    IF ((R(1).ge.myx(jjx)).and.(R(1).lt.myx(jjx+1))) THEN 
      l(1)=jjx
      EXIT
    ENDIF
   ENDDO
   !IF (l(1).eq.(-nx)) fxflag=.TRUE.
   DO jjy=5,ny-5
    IF ((R(2).ge.myy(jjy)).and.(R(2).lt.myy(jjy+1))) THEN
      l(2)=jjy
      EXIT
    ENDIF
   ENDDO
   !IF (l(2).eq.(-ny)) fyflag=.TRUE.
   DO jjz=5,nz-5
    IF ((R(3).ge.myz(jjz)).and.(R(3).lt.myz(jjz+1))) THEN
      l(3)=jjz
      EXIT
    ENDIF
   ENDDO
   IF (l(3).eq.(-nz)) print*, R(3), myz(1), myz(5)

! No guarantee we have more than one frame. IF we have one, this routine doesn't bother interpolating in time
  IF (nframes.gt.1) THEN
   ALLOCATE(dgt(nframes-1),odgt(nframes-1))
   dgt=ltimes(2:nframes)-ltimes(1:nframes-1)
   odgt=1.0_num/dgt
   DO jjt=1,nframes
    IF ((T.GE.ltimes(jjt)).AND.((T.LT.ltimes(jjt+1)))) THEN
      l(4)=jjt 
      !print*, (T-ltimes(l(4)))*odgt(l(4))
      EXIT
    !ELSE
    !  PRINT *, 'CANNOT FIND TIME IN LARE TIME RANGE'
    ENDIF
   ENDDO
   rpt=1
   ALLOCATE(bxt(2), byt(2), bzt(2), Ext(2), Eyt(2), Ezt(2), jxt(2), jyt(2), jzt(2))
   ALLOCATE(dbxdxt(2),dbxdyt(2),dbxdzt(2),dbydxt(2),dbydyt(2),dbydzt(2),dbzdxt(2),dbzdyt(2),dbzdzt(2))
   ALLOCATE(dExdxt(2),dExdyt(2),dExdzt(2),dEydxt(2),dEydyt(2),dEydzt(2),dEzdxt(2),dEzdyt(2),dEzdzt(2))
  ELSE
   l(4)=1
   rpt=0
   ALLOCATE(bxt(1), byt(1), bzt(1), Ext(1), Eyt(1), Ezt(1), jxt(1), jyt(1), jzt(1))
   ALLOCATE(dbxdxt(1),dbxdyt(1),dbxdzt(1),dbydxt(1),dbydyt(1),dbydzt(1),dbzdxt(1),dbzdyt(1),dbzdzt(1))
   ALLOCATE(dExdxt(1),dExdyt(1),dExdzt(1),dEydxt(1),dEydyt(1),dEydzt(1),dEzdxt(1),dEzdyt(1),dEzdzt(1))
  !what happens if there is one frame to work with?
  ENDIF 


   coffset=(/(R(1)-myx(l(1)))*odg(1),(R(2)-myy(l(2)))*odg(2),(R(3)-myz(l(3)))*odg(3)/)

   
   ! --STEP TWO-- !
   ! begin interpolation of basic variables we already have, e.g. bx, by, bz, vx, vy, vz
  DO it=0,rpt	! NEED TO REPEAT FOR INTERPOLATION BETWEEN FRAMES, JT DEC 2015
  
   !temp=linterp3d(R(1)-myx(l(1)),R(2)-myy(l(2)),R(3)-myz(l(3)),&	!1   
   temp=linterp3d(coffset(1),coffset(2),coffset(3), &
   bx(l(1),l(2),l(3),l(4)+it),bx(l(1)+1,l(2),l(3),l(4)+it),bx(l(1),l(2)+1,l(3),l(4)+it),&
   bx(l(1)+1,l(2)+1,l(3),l(4)+it),bx(l(1),l(2),l(3)+1,l(4)+it),bx(l(1)+1,l(2),l(3)+1,l(4)+it), &
   bx(l(1),l(2)+1,l(3)+1,l(4)+it),bx(l(1)+1,l(2)+1,l(3)+1,l(4)+it))
   bxt(it+1)=temp
   
   !temp=linterp3d(R(1)-myx(l(1)),R(2)-myy(l(2)),R(3)-myz(l(3)),&	!2
   temp=linterp3d(coffset(1),coffset(2),coffset(3), &
   by(l(1),l(2),l(3),l(4)+it),by(l(1)+1,l(2),l(3),l(4)+it),by(l(1),l(2)+1,l(3),l(4)+it),&
   by(l(1)+1,l(2)+1,l(3),l(4)+it),by(l(1),l(2),l(3)+1,l(4)+it),by(l(1)+1,l(2),l(3)+1,l(4)+it), &
   by(l(1),l(2)+1,l(3)+1,l(4)+it),by(l(1)+1,l(2)+1,l(3)+1,l(4)+it))
   byt(it+1)=temp

   !temp=linterp3d(R(1)-myx(l(1)),R(2)-myy(l(2)),R(3)-myz(l(3)),&	!3
   temp=linterp3d(coffset(1),coffset(2),coffset(3), &
   bz(l(1),l(2),l(3),l(4)+it),bz(l(1)+1,l(2),l(3),l(4)+it),bz(l(1),l(2)+1,l(3),l(4)+it),&
   bz(l(1)+1,l(2)+1,l(3),l(4)+it),bz(l(1),l(2),l(3)+1,l(4)+it),bz(l(1)+1,l(2),l(3)+1,l(4)+it), &
   bz(l(1),l(2)+1,l(3)+1,l(4)+it),bz(l(1)+1,l(2)+1,l(3)+1,l(4)+it))
   bzt(it+1)=temp
   
   ! --STEP THREE--!
   ! create temporary arrays around target cell and calculate derivs using extra cells
   ! - how many cells depend on which finite difference routine used
   ! - need db/dr->j and etadj/dr for E-derivs
   ! JT OCT2014: commented out O(h) derivs, and insert O(h^2)
   
   ! begin with temporary B arrays:
   ALLOCATE(mbx(-4:5,-4:5,-4:5))
   ALLOCATE(mby(-4:5,-4:5,-4:5))
   ALLOCATE(mbz(-4:5,-4:5,-4:5))

   mbx=bx(l(1)-4:l(1)+5,l(2)-4:l(2)+5,l(3)-4:l(3)+5,l(4)+it)
   mby=by(l(1)-4:l(1)+5,l(2)-4:l(2)+5,l(3)-4:l(3)+5,l(4)+it)
   mbz=bz(l(1)-4:l(1)+5,l(2)-4:l(2)+5,l(3)-4:l(3)+5,l(4)+it)


   ! calculate x-derivs first..
   ALLOCATE(dmbxdx(-2:3,-4:5,-4:5))
   ALLOCATE(dmbydx(-2:3,-4:5,-4:5))
   ALLOCATE(dmbzdx(-2:3,-4:5,-4:5))

   DO iz=-4,5   
    DO iy=-4,5
     DO ix=-2,3
      dmbxdx(ix,iy,iz)=(-mbx(ix+2,iy,iz)+8.0_num*mbx(ix+1,iy,iz)-8.0_num*mbx(ix-1,iy,iz)+mbx(ix-2,iy,iz))*odg(1)*oneotwelve
      dmbydx(ix,iy,iz)=(-mby(ix+2,iy,iz)+8.0_num*mby(ix+1,iy,iz)-8.0_num*mby(ix-1,iy,iz)+mby(ix-2,iy,iz))*odg(1)*oneotwelve
      dmbzdx(ix,iy,iz)=(-mbz(ix+2,iy,iz)+8.0_num*mbz(ix+1,iy,iz)-8.0_num*mbz(ix-1,iy,iz)+mbz(ix-2,iy,iz))*odg(1)*oneotwelve
     END DO
    END DO
   END DO
   
   ! ..followed by y-derivs..
   ALLOCATE(dmbxdy(-4:5,-2:3,-4:5))
   ALLOCATE(dmbydy(-4:5,-2:3,-4:5))
   ALLOCATE(dmbzdy(-4:5,-2:3,-4:5))

   DO iz=-4,5   
    DO iy=-2,3
     DO ix=-4,5
      dmbxdy(ix,iy,iz)=(-mbx(ix,iy+2,iz)+8.0_num*mbx(ix,iy+1,iz)-8.0_num*mbx(ix,iy-1,iz)+mbx(ix,iy-2,iz))*odg(2)*oneotwelve
      dmbydy(ix,iy,iz)=(-mby(ix,iy+2,iz)+8.0_num*mby(ix,iy+1,iz)-8.0_num*mby(ix,iy-1,iz)+mby(ix,iy-2,iz))*odg(2)*oneotwelve
      dmbzdy(ix,iy,iz)=(-mbz(ix,iy+2,iz)+8.0_num*mbz(ix,iy+1,iz)-8.0_num*mbz(ix,iy-1,iz)+mbz(ix,iy-2,iz))*odg(2)*oneotwelve
     END DO
    END DO
   END DO

   ! .. and finally z derivs.
   ALLOCATE(dmbxdz(-4:5,-4:5,-2:3))
   ALLOCATE(dmbydz(-4:5,-4:5,-2:3))
   ALLOCATE(dmbzdz(-4:5,-4:5,-2:3))

   DO iz=-2,3
    DO iy=-4,5   
     DO ix=-4,5
      dmbxdz(ix,iy,iz)=oneotwelve*(-mbx(ix,iy,iz+2)+8.0_num*mbx(ix,iy,iz+1)-8.0_num*mbx(ix,iy,iz-1)+mbx(ix,iy,iz-2))*odg(3)
      dmbydz(ix,iy,iz)=oneotwelve*(-mby(ix,iy,iz+2)+8.0_num*mby(ix,iy,iz+1)-8.0_num*mby(ix,iy,iz-1)+mby(ix,iy,iz-2))*odg(3)
      dmbzdz(ix,iy,iz)=oneotwelve*(-mbz(ix,iy,iz+2)+8.0_num*mbz(ix,iy,iz+1)-8.0_num*mbz(ix,iy,iz-1)+mbz(ix,iy,iz-2))*odg(3)
     ENDDO
    ENDDO
   ENDDO

   ! now calculate temporary currents and also local eta values
   ALLOCATE(mjx(-2:3,-2:3,-2:3))
   ALLOCATE(mjy(-2:3,-2:3,-2:3))
   ALLOCATE(mjz(-2:3,-2:3,-2:3))
   ALLOCATE(mEx(-2:3,-2:3,-2:3))
   ALLOCATE(mEy(-2:3,-2:3,-2:3))
   ALLOCATE(mEz(-2:3,-2:3,-2:3))
   ALLOCATE(meta(-2:3,-2:3,-2:3))
   
   !but Steph and Duncan have already calculated the local currents.
   
   mjx=jx(l(1)-2:l(1)+3,l(2)-2:l(2)+3,l(3)-2:l(3)+3,l(4)+it)
   mjy=jy(l(1)-2:l(1)+3,l(2)-2:l(2)+3,l(3)-2:l(3)+3,l(4)+it)
   mjz=jz(l(1)-2:l(1)+3,l(2)-2:l(2)+3,l(3)-2:l(3)+3,l(4)+it)
   
   DO iz=-2,3  
    DO iy=-2,3
     DO ix=-2,3
      !mjx(ix,iy,iz)=dmbzdy(ix,iy,iz)-dmbydz(ix,iy,iz)
      !mjy(ix,iy,iz)=dmbxdz(ix,iy,iz)-dmbzdx(ix,iy,iz)
      !mjz(ix,iy,iz)=dmbydx(ix,iy,iz)-dmbxdy(ix,iy,iz)
      modj=(mjx(ix,iy,iz)*mjx(ix,iy,iz)+mjy(ix,iy,iz)*mjy(ix,iy,iz)+mjz(ix,iy,iz)*mjz(ix,iy,iz))**0.5_num 
      !!OCT 2014: replace step function with smooth eta ramp using tanh profile:
      !IF (modj.gt.jcrit) THEN
      ! meta(ix,iy,iz)=eta
      !ELSE
      ! meta(ix,iy,iz)=0.0_num
      !ENDIF
      meta(ix,iy,iz)=0.5_num*(tanh((modj-jcrit)/rwidth)+1.0_num)*eta+etabkg          
     END DO
    END DO
   END DO

   ! now interpolate derivs to correct point:
   !linterp3d(dx,dy,dz,f000,f100,f010,f110,f001,f101,f011,f111)    
   !dbxdxt=linterp3d(R(1)-myx(l(1)),R(2)-myy(l(2)),R(3)-myz(l(3)),&	!22
   dbxdxt(it+1)=linterp3d(coffset(1),coffset(2),coffset(3), &
    dmbxdx(0,0,0), dmbxdx(1,0,0), dmbxdx(0,1,0), dmbxdx(1,1,0),&
    dmbxdx(0,0,1), dmbxdx(1,0,1), dmbxdx(0,1,1), dmbxdx(1,1,1))
   !dbydxt=linterp3d(R(1)-myx(l(1)),R(2)-myy(l(2)),R(3)-myz(l(3)),&	!23
   dbydxt(it+1)=linterp3d(coffset(1),coffset(2),coffset(3), &
    dmbydx(0,0,0), dmbydx(1,0,0), dmbydx(0,1,0), dmbydx(1,1,0),&
    dmbydx(0,0,1), dmbydx(1,0,1), dmbydx(0,1,1), dmbydx(1,1,1))
   !dbzdxt=linterp3d(R(1)-myx(l(1)),R(2)-myy(l(2)),R(3)-myz(l(3)),&	!24
   dbzdxt(it+1)=linterp3d(coffset(1),coffset(2),coffset(3), &
    dmbzdx(0,0,0), dmbzdx(1,0,0), dmbzdx(0,1,0), dmbzdx(1,1,0),&
    dmbzdx(0,0,1), dmbzdx(1,0,1), dmbzdx(0,1,1), dmbzdx(1,1,1))
   !dbxdyt=linterp3d(R(1)-myx(l(1)),R(2)-myy(l(2)),R(3)-myz(l(3)),&	!25
   dbxdyt(it+1)=linterp3d(coffset(1),coffset(2),coffset(3), &
    dmbxdy(0,0,0), dmbxdy(1,0,0), dmbxdy(0,1,0), dmbxdy(1,1,0),&
    dmbxdy(0,0,1), dmbxdy(1,0,1), dmbxdy(0,1,1), dmbxdy(1,1,1))
   !dbydyt=linterp3d(R(1)-myx(l(1)),R(2)-myy(l(2)),R(3)-myz(l(3)),&	!26
   dbydyt(it+1)=linterp3d(coffset(1),coffset(2),coffset(3), &
    dmbydy(0,0,0), dmbydy(1,0,0), dmbydy(0,1,0), dmbydy(1,1,0),&
    dmbydy(0,0,1), dmbydy(1,0,1), dmbydy(0,1,1), dmbydy(1,1,1))
   !dbzdyt=linterp3d(R(1)-myx(l(1)),R(2)-myy(l(2)),R(3)-myz(l(3)),&	!27
   dbzdyt(it+1)=linterp3d(coffset(1),coffset(2),coffset(3), &
    dmbzdy(0,0,0), dmbzdy(1,0,0), dmbzdy(0,1,0), dmbzdy(1,1,0),&
    dmbzdy(0,0,1), dmbzdy(1,0,1), dmbzdy(0,1,1), dmbzdy(1,1,1))            
   !dbxdzt=linterp3d(R(1)-myx(l(1)),R(2)-myy(l(2)),R(3)-myz(l(3)),&	!28
   dbxdzt(it+1)=linterp3d(coffset(1),coffset(2),coffset(3), &
    dmbxdz(0,0,0), dmbxdz(1,0,0), dmbxdz(0,1,0), dmbxdz(1,1,0),&
    dmbxdz(0,0,1), dmbxdz(1,0,1), dmbxdz(0,1,1), dmbxdz(1,1,1))
   !dbydzt=linterp3d(R(1)-myx(l(1)),R(2)-myy(l(2)),R(3)-myz(l(3)),&	!29
   dbydzt(it+1)=linterp3d(coffset(1),coffset(2),coffset(3), &
    dmbydz(0,0,0), dmbydz(1,0,0), dmbydz(0,1,0), dmbydz(1,1,0),&
    dmbydz(0,0,1), dmbydz(1,0,1), dmbydz(0,1,1), dmbydz(1,1,1))
   !dbzdzt=linterp3d(R(1)-myx(l(1)),R(2)-myy(l(2)),R(3)-myz(l(3)),&	!30
   dbzdzt(it+1)=linterp3d(coffset(1),coffset(2),coffset(3), &
    dmbzdz(0,0,0), dmbzdz(1,0,0), dmbzdz(0,1,0), dmbzdz(1,1,0),&
    dmbzdz(0,0,1), dmbzdz(1,0,1), dmbzdz(0,1,1), dmbzdz(1,1,1))
   
   DEALLOCATE(dmbxdx)
   DEALLOCATE(dmbydx)
   DEALLOCATE(dmbzdx)
   DEALLOCATE(dmbxdy)
   DEALLOCATE(dmbydy)
   DEALLOCATE(dmbzdy)
   DEALLOCATE(dmbxdz)
   DEALLOCATE(dmbydz)
   DEALLOCATE(dmbzdz)

   temp=linterp3d(coffset(1),coffset(2),coffset(3), &
   jx(l(1),l(2),l(3),l(4)+it),jx(l(1)+1,l(2),l(3),l(4)+it),jx(l(1),l(2)+1,l(3),l(4)+it),&
   jx(l(1)+1,l(2)+1,l(3),l(4)+it),jx(l(1),l(2),l(3)+1,l(4)+it),jx(l(1)+1,l(2),l(3)+1,l(4)+it), &
   jx(l(1),l(2)+1,l(3)+1,l(4)+it),jx(l(1)+1,l(2)+1,l(3)+1,l(4)+it))
   jxt(it+1)=temp
   
   !temp=linterp3d(R(1)-myx(l(1)),R(2)-myy(l(2)),R(3)-myz(l(3)),&	!2
   temp=linterp3d(coffset(1),coffset(2),coffset(3), &
   jy(l(1),l(2),l(3),l(4)+it),jy(l(1)+1,l(2),l(3),l(4)+it),jy(l(1),l(2)+1,l(3),l(4)+it),&
   jy(l(1)+1,l(2)+1,l(3),l(4)+it),jy(l(1),l(2),l(3)+1,l(4)+it),jy(l(1)+1,l(2),l(3)+1,l(4)+it), &
   jy(l(1),l(2)+1,l(3)+1,l(4)+it),jy(l(1)+1,l(2)+1,l(3)+1,l(4)+it))
   jyt(it+1)=temp

   !temp=linterp3d(R(1)-myx(l(1)),R(2)-myy(l(2)),R(3)-myz(l(3)),&	!3
   temp=linterp3d(coffset(1),coffset(2),coffset(3), &
   jz(l(1),l(2),l(3),l(4)+it),jz(l(1)+1,l(2),l(3),l(4)+it),jz(l(1),l(2)+1,l(3),l(4)+it),&
   jz(l(1)+1,l(2)+1,l(3),l(4)+it),jz(l(1),l(2),l(3)+1,l(4)+it),jz(l(1)+1,l(2),l(3)+1,l(4)+it), &
   jz(l(1),l(2)+1,l(3)+1,l(4)+it),jz(l(1)+1,l(2)+1,l(3)+1,l(4)+it))
   jzt(it+1)=temp

   !! interpolate j too at this stage:
   !!mjx(0:3,0:3,0:3)	<-x=1-2,y=1-2,z=1-2
   !!temp=linterp3d(R(1)-myx(l(1)),R(2)-myy(l(2)),R(3)-myz(l(3)),&	!7
   !temp=linterp3d(coffset(1),coffset(2),coffset(3), &
   ! mjx(0,0,0), mjx(1,0,0), mjx(0,1,0), mjx(1,1,0), &
   ! mjx(0,0,1), mjx(1,0,1), mjx(0,1,1), mjx(1,1,1))
   !jxt(it+1)=temp
   !!temp=linterp3d(R(1)-myx(l(1)),R(2)-myy(l(2)),R(3)-myz(l(3)),&	!8
   !temp=linterp3d(coffset(1),coffset(2),coffset(3), &
   ! mjy(0,0,0), mjy(1,0,0), mjy(0,1,0), mjy(1,1,0), &
   ! mjy(0,0,1), mjy(1,0,1), mjy(0,1,1), mjy(1,1,1))
   !jyt(it+1)=temp
   !!temp=linterp3d(R(1)-myx(l(1)),R(2)-myy(l(2)),R(3)-myz(l(3)),&	!9
   !temp=linterp3d(coffset(1),coffset(2),coffset(3), &
   ! mjz(0,0,0), mjz(1,0,0), mjz(0,1,0), mjz(1,1,0), &
   ! mjz(0,0,1), mjz(1,0,1), mjz(0,1,1), mjz(1,1,1))
   !jzt(it+1)=temp
   
   ! --STEP FOUR--
   ! Ohm's Law: E=vxB-etaJ OR E=-etaJ => v=0 IN THIS CASE
   DO iz=-2,3 
    DO iy=-2,3
     DO ix=-2,3
      mEx(ix,iy,iz)=meta(ix,iy,iz)*mjx(ix,iy,iz)
      mEy(ix,iy,iz)=meta(ix,iy,iz)*mjy(ix,iy,iz)
      mEz(ix,iy,iz)=meta(ix,iy,iz)*mjz(ix,iy,iz)
     ENDDO
    ENDDO
   ENDDO
     
     
   DEALLOCATE(mbx)
   DEALLOCATE(mby)
   DEALLOCATE(mbz)
   DEALLOCATE(mjx)
   DEALLOCATE(mjy)
   DEALLOCATE(mjz)    
   DEALLOCATE(meta)
   
   ! --STEP FIVE-- !
   ! finally, calculate derivatives of electric field:

   ALLOCATE(dmExdx(0:1,-2:3,-2:3))
   ALLOCATE(dmEydx(0:1,-2:3,-2:3))
   ALLOCATE(dmEzdx(0:1,-2:3,-2:3))
   DO iz=-2,3   
    DO iy=-2,3
     DO ix=0,1
      dmExdx(ix,iy,iz)=oneotwelve*(-mEx(ix+2,iy,iz)+8.0_num*mEx(ix+1,iy,iz)-8.0_num*mEx(ix-1,iy,iz)+mEx(ix-2,iy,iz))*odg(1)
      dmEydx(ix,iy,iz)=oneotwelve*(-mEy(ix+2,iy,iz)+8.0_num*mEy(ix+1,iy,iz)-8.0_num*mEy(ix-1,iy,iz)+mEy(ix-2,iy,iz))*odg(1)
      dmEzdx(ix,iy,iz)=oneotwelve*(-mEz(ix+2,iy,iz)+8.0_num*mEz(ix+1,iy,iz)-8.0_num*mEz(ix-1,iy,iz)+mEz(ix-2,iy,iz))*odg(1)
     ENDDO
    ENDDO
   ENDDO
   
   ALLOCATE(dmExdy(-2:3,0:1,-2:3))
   ALLOCATE(dmEydy(-2:3,0:1,-2:3))
   ALLOCATE(dmEzdy(-2:3,0:1,-2:3))  
   DO iz=-2,3   
    DO iy=0,1
     DO ix=-2,3  
      dmExdy(ix,iy,iz)=oneotwelve*(-mEx(ix,iy+2,iz)+8.0_num*mEx(ix,iy+1,iz)-8.0_num*mEx(ix,iy-1,iz)+mEx(ix,iy-2,iz))*odg(2)
      dmEydy(ix,iy,iz)=oneotwelve*(-mEy(ix,iy+2,iz)+8.0_num*mEy(ix,iy+1,iz)-8.0_num*mEy(ix,iy-1,iz)+mEy(ix,iy-2,iz))*odg(2)
      dmEzdy(ix,iy,iz)=oneotwelve*(-mEz(ix,iy+2,iz)+8.0_num*mEz(ix,iy+1,iz)-8.0_num*mEz(ix,iy-1,iz)+mEz(ix,iy-2,iz))*odg(2)
     ENDDO
    ENDDO
   ENDDO   
   
   ALLOCATE(dmExdz(-2:3,-2:3,0:1))
   ALLOCATE(dmEydz(-2:3,-2:3,0:1))
   ALLOCATE(dmEzdz(-2:3,-2:3,0:1))   
   DO iz=0,1   
    DO iy=-2,3
     DO ix=-2,3
      dmExdz(ix,iy,iz)=oneotwelve*(-mEx(ix,iy,iz+2)+8.0_num*mEx(ix,iy,iz+1)-8.0_num*mEx(ix,iy,iz-1)+mEx(ix,iy,iz-2))*odg(3)
      dmEydz(ix,iy,iz)=oneotwelve*(-mEy(ix,iy,iz+2)+8.0_num*mEy(ix,iy,iz+1)-8.0_num*mEy(ix,iy,iz-1)+mEy(ix,iy,iz-2))*odg(3)
      dmEzdz(ix,iy,iz)=oneotwelve*(-mEz(ix,iy,iz+2)+8.0_num*mEz(ix,iy,iz+1)-8.0_num*mEz(ix,iy,iz-1)+mEz(ix,iy,iz-2))*odg(3)
     ENDDO
    ENDDO
   ENDDO
   
   !dmExdx(1:2,0:3,0:3)	<-x=1-2,y=1-2,z=1-2
   !linterp3d(dx,dy,dz,f000,f100,f010,f110,f001,f101,f011,f111)
         
   !temp=linterp3d(R(1)-myx(l(1)),R(2)-myy(l(2)),R(3)-myz(l(3)),&	!10
   temp=linterp3d(coffset(1),coffset(2),coffset(3), &
    dmExdx(0,0,0), dmExdx(1,0,0), dmExdx(0,1,0), dmExdx(1,1,0),&
    dmExdx(0,0,1), dmExdx(1,0,1), dmExdx(0,1,1), dmExdx(1,1,1))
   dExdxt(it+1)=temp
   !temp=linterp3d(R(1)-myx(l(1)),R(2)-myy(l(2)),R(3)-myz(l(3)),&	!11
   temp=linterp3d(coffset(1),coffset(2),coffset(3), &
    dmEydx(0,0,0), dmEydx(1,0,0), dmEydx(0,1,0), dmEydx(1,1,0),&
    dmEydx(0,0,1), dmEydx(1,0,1), dmEydx(0,1,1), dmEydx(1,1,1))
   dEydxt(it+1)=temp
   !temp=linterp3d(R(1)-myx(l(1)),R(2)-myy(l(2)),R(3)-myz(l(3)),&	!12
   temp=linterp3d(coffset(1),coffset(2),coffset(3), &
    dmEzdx(0,0,0), dmEzdx(1,0,0), dmEzdx(0,1,0), dmEzdx(1,1,0),&
    dmEzdx(0,0,1), dmEzdx(1,0,1), dmEzdx(0,1,1), dmEzdx(1,1,1))
   dEzdxt(it+1)=temp
   !temp=linterp3d(R(1)-myx(l(1)),R(2)-myy(l(2)),R(3)-myz(l(3)),&	!13
   temp=linterp3d(coffset(1),coffset(2),coffset(3), &
    dmExdy(0,0,0), dmExdy(1,0,0), dmExdy(0,1,0), dmExdy(1,1,0),&
    dmExdy(0,0,1), dmExdy(1,0,1), dmExdy(0,1,1), dmExdy(1,1,1))
   dExdyt(it+1)=temp
   !temp=linterp3d(R(1)-myx(l(1)),R(2)-myy(l(2)),R(3)-myz(l(3)),&	!14
   temp=linterp3d(coffset(1),coffset(2),coffset(3), &
    dmEydy(0,0,0), dmEydy(1,0,0), dmEydy(0,1,0), dmEydy(1,1,0),&
    dmEydy(0,0,1), dmEydy(1,0,1), dmEydy(0,1,1), dmEydy(1,1,1))
   dEydyt(it+1)=temp
   !temp=linterp3d(R(1)-myx(l(1)),R(2)-myy(l(2)),R(3)-myz(l(3)),&	!15
   temp=linterp3d(coffset(1),coffset(2),coffset(3), &
    dmEzdy(0,0,0), dmEzdy(1,0,0), dmEzdy(0,1,0), dmEzdy(1,1,0),&
    dmEzdy(0,0,1), dmEzdy(1,0,1), dmEzdy(0,1,1), dmEzdy(1,1,1))
   dEzdyt(it+1)=temp
   !temp=linterp3d(R(1)-myx(l(1)),R(2)-myy(l(2)),R(3)-myz(l(3)),&	!16
   temp=linterp3d(coffset(1),coffset(2),coffset(3), &
    dmExdz(0,0,0), dmExdz(1,0,0), dmExdz(0,1,0), dmExdz(1,1,0),&
    dmExdz(0,0,1), dmExdz(1,0,1), dmExdz(0,1,1), dmExdz(1,1,1))
   dExdzt(it+1)=temp
   !temp=linterp3d(R(1)-myx(l(1)),R(2)-myy(l(2)),R(3)-myz(l(3)),&	!17
   temp=linterp3d(coffset(1),coffset(2),coffset(3), &
    dmEydz(0,0,0), dmEydz(1,0,0), dmEydz(0,1,0), dmEydz(1,1,0),&
    dmEydz(0,0,1), dmEydz(1,0,1), dmEydz(0,1,1), dmEydz(1,1,1))
   dEydzt(it+1)=temp
   !temp=linterp3d(R(1)-myx(l(1)),R(2)-myy(l(2)),R(3)-myz(l(3)),&	!18
   temp=linterp3d(coffset(1),coffset(2),coffset(3), &
    dmEzdz(0,0,0), dmEzdz(1,0,0), dmEzdz(0,1,0), dmEzdz(1,1,0),&
    dmEzdz(0,0,1), dmEzdz(1,0,1), dmEzdz(0,1,1), dmEzdz(1,1,1))
   dEzdzt(it+1)=temp
   
   DEALLOCATE(dmExdx)
   DEALLOCATE(dmEydx)
   DEALLOCATE(dmEzdx)
   DEALLOCATE(dmExdy)
   DEALLOCATE(dmEydy)
   DEALLOCATE(dmEzdy)
   DEALLOCATE(dmExdz)
   DEALLOCATE(dmEydz)
   DEALLOCATE(dmEzdz)
   
   !mEx(0:3,0:3,0:3)	<-x=1-2,y=1-2,z=1-2
   !linterp3d(dx,dy,dz,f000,f100,f010,f110,f001,f101,f011,f111)
         
   !temp=linterp3d(R(1)-myx(l(1)),R(2)-myy(l(2)),R(3)-myz(l(3)),&	!19
   temp=linterp3d(coffset(1),coffset(2),coffset(3), &
    mEx(0,0,0), mEx(1,0,0), mEx(0,1,0), mEx(1,1,0),&
    mEx(0,0,1), mEx(1,0,1), mEx(0,1,1), mEx(1,1,1))
   Ext(it+1)=temp 
   !temp=linterp3d(R(1)-myx(l(1)),R(2)-myy(l(2)),R(3)-myz(l(3)),&	!20
   temp=linterp3d(coffset(1),coffset(2),coffset(3), &
    mEy(0,0,0), mEy(1,0,0), mEy(0,1,0), mEy(1,1,0),&
    mEy(0,0,1), mEy(1,0,1), mEy(0,1,1), mEy(1,1,1))
   Eyt(it+1)=temp
   !temp=linterp3d(R(1)-myx(l(1)),R(2)-myy(l(2)),R(3)-myz(l(3)),&	!21
   temp=linterp3d(coffset(1),coffset(2),coffset(3), &
    mEz(0,0,0), mEz(1,0,0), mEz(0,1,0), mEz(1,1,0),&
    mEz(0,0,1), mEz(1,0,1), mEz(0,1,1), mEz(1,1,1))
   Ezt(it+1)=temp 
    
   DEALLOCATE(mEx)
   DEALLOCATE(mEy)
   DEALLOCATE(mEz)
   
  END DO
   
   !we can try and sort other variables later (e.g. density/temp) when we're happy

   ! its possible to run lare3d in real units - lets make sure things are normalised if this happens

   IF (nframes.gt.1) THEN

    N3d(1)=linterp1d((T-ltimes(l(4)))*odgt(l(4)),bxt(1),bxt(2))
    N3d(2)=linterp1d((T-ltimes(l(4)))*odgt(l(4)),byt(1),byt(2))
    N3d(3)=linterp1d((T-ltimes(l(4)))*odgt(l(4)),bzt(1),bzt(2))
    N3d(4)=0.0_num
    N3d(5)=0.0_num
    N3d(6)=0.0_num
    N3d(7)=linterp1d((T-ltimes(l(4)))*odgt(l(4)),Ext(1),Ext(2))
    N3d(8)=linterp1d((T-ltimes(l(4)))*odgt(l(4)),Eyt(1),Eyt(2))
    N3d(9)=linterp1d((T-ltimes(l(4)))*odgt(l(4)),Ezt(1),Ezt(2))   
    N3d(10)=linterp1d((T-ltimes(l(4)))*odgt(l(4)),jxt(1),jxt(2))
    N3d(11)=linterp1d((T-ltimes(l(4)))*odgt(l(4)),jyt(1),jyt(2))
    N3d(12)=linterp1d((T-ltimes(l(4)))*odgt(l(4)),jzt(1),jzt(2))
    N3d(13)=linterp1d((T-ltimes(l(4)))*odgt(l(4)),dbxdxt(1),dbxdxt(2))
    N3d(14)=linterp1d((T-ltimes(l(4)))*odgt(l(4)),dbydxt(1),dbydxt(2))
    N3d(15)=linterp1d((T-ltimes(l(4)))*odgt(l(4)),dbzdxt(1),dbzdxt(2))
    N3d(16)=linterp1d((T-ltimes(l(4)))*odgt(l(4)),dbxdyt(1),dbxdyt(2))
    N3d(17)=linterp1d((T-ltimes(l(4)))*odgt(l(4)),dbydyt(1),dbydyt(2))
    N3d(18)=linterp1d((T-ltimes(l(4)))*odgt(l(4)),dbzdyt(1),dbzdyt(2))
    N3d(19)=linterp1d((T-ltimes(l(4)))*odgt(l(4)),dbxdzt(1),dbxdzt(2))
    N3d(20)=linterp1d((T-ltimes(l(4)))*odgt(l(4)),dbydzt(1),dbydzt(2))
    N3d(21)=linterp1d((T-ltimes(l(4)))*odgt(l(4)),dbzdzt(1),dbzdzt(2))
    N3d(22)=linterp1d((T-ltimes(l(4)))*odgt(l(4)),dExdxt(1),dExdxt(2))
    N3d(23)=linterp1d((T-ltimes(l(4)))*odgt(l(4)),dEydxt(1),dEydxt(2))
    N3d(24)=linterp1d((T-ltimes(l(4)))*odgt(l(4)),dEzdxt(1),dEzdxt(2))
    N3d(25)=linterp1d((T-ltimes(l(4)))*odgt(l(4)),dExdyt(1),dExdyt(2))
    N3d(26)=linterp1d((T-ltimes(l(4)))*odgt(l(4)),dEydyt(1),dEydyt(2))
    N3d(27)=linterp1d((T-ltimes(l(4)))*odgt(l(4)),dEzdyt(1),dEzdyt(2))
    N3d(28)=linterp1d((T-ltimes(l(4)))*odgt(l(4)),dExdzt(1),dExdzt(2))
    N3d(29)=linterp1d((T-ltimes(l(4)))*odgt(l(4)),dEydzt(1),dEydzt(2))
    N3d(30)=linterp1d((T-ltimes(l(4)))*odgt(l(4)),dEzdzt(1),dEzdzt(2))
    N3d(31)=(bxt(2)-bxt(1))*odgt(l(4))
    N3d(32)=(byt(2)-byt(1))*odgt(l(4))
    N3d(33)=(bzt(2)-bzt(1))*odgt(l(4))
    N3d(34)=(Ext(2)-Ext(1))*odgt(l(4))	! as we have multiple snapshots, we can calculate time derivs finally!
    N3d(35)=(Eyt(2)-Eyt(1))*odgt(l(4))
    N3d(36)=(Ezt(2)-Ezt(1))*odgt(l(4))
    DEALLOCATE(dgt,odgt)
    
   ELSE
   
    N3d(1)=bxt(1)
    N3d(2)=byt(1)
    N3d(3)=bzt(1)
    N3d(4)=0.0_num
    N3d(5)=0.0_num
    N3d(6)=0.0_num
    N3d(7)=Ext(1)
    N3d(8)=Eyt(1)
    N3d(9)=Ezt(1)   
    N3d(10)=jxt(1)
    N3d(11)=jyt(1)
    N3d(12)=jzt(1)
    N3d(13)=dbxdxt(1)
    N3d(14)=dbydxt(1)
    N3d(15)=dbzdxt(1)
    N3d(16)=dbxdyt(1)
    N3d(17)=dbydyt(1)
    N3d(18)=dbzdyt(1)
    N3d(19)=dbxdzt(1)
    N3d(20)=dbydzt(1)
    N3d(21)=dbzdzt(1)
    N3d(22)=dExdxt(1)
    N3d(23)=dEydxt(1)
    N3d(24)=dEzdxt(1)
    N3d(25)=dExdyt(1)
    N3d(26)=dEydyt(1)
    N3d(27)=dEzdyt(1)
    N3d(28)=dExdzt(1)
    N3d(29)=dEydzt(1)
    N3d(30)=dEzdzt(1)
    N3d(31:36)=0.0_num
   ENDIF
   
   DEALLOCATE(bxt, byt, bzt, Ext, Eyt, Ezt, jxt, jyt, jzt)
   DEALLOCATE(dbxdxt,dbxdyt,dbxdzt,dbydxt,dbydyt,dbydzt,dbzdxt,dbzdyt,dbzdzt)
   DEALLOCATE(dExdxt,dExdyt,dExdzt,dEydxt,dEydyt,dEydzt,dEzdxt,dEzdyt,dEzdzt)

   RETURN

END FUNCTION N3d
!------------------------------------------------------------------------------!    
SUBROUTINE NLFFFIELDS(R,T,E,B,DBDX,DBDY,DBDZ,DBDT,DEDX,DEDY,DEDZ,DEDT,Vf)

 REAL(num), DIMENSION(3),INTENT(OUT)	:: B,E
 REAL(num), DIMENSION(3),INTENT(OUT)	:: DBDX,DBDY,DBDZ,DBDT,DEDX,DEDY,DEDZ,DEDT
 REAL(num), DIMENSION(3),INTENT(IN)	:: R
 REAL(num), DIMENSION(3)  		:: Vf
 REAL(num), INTENT(IN) 			:: T
 REAL(num), DIMENSION(36)		:: iquants

! For the NLFFF code, each of the returned variables will be in REAL units:
! normalise these quantities so the orbits can also work in non-dimensional quantities

   iquants=N3d(R,T)
   B	=iquants(1:3)/Bscl
   Vf	=0.0_num
   E	=iquants(7:9)/Escl
   DBDX	=iquants(10:12)/Bscl*lscl
   DBDY	=iquants(13:15)/Bscl*lscl
   DBDZ	=iquants(16:18)/Bscl*lscl
   DBDT	=0.0_num
   DEDX	=iquants(19:21)/Escl*lscl
   DEDY	=iquants(22:24)/Escl*lscl
   DEDZ	=iquants(25:27)/Escl*lscl
   DEDT	=0.0_num

END SUBROUTINE NLFFFIELDS

END MODULE NLFF_fields

