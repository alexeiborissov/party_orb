MODULE l2dc_fields
! Module contains several routines which set up the global lare grid and variable arrays.
  
  USE GLOBAL
  USE iocommon
  USE iocontrol
  USE input
  USE input_cartesian
  USE lare_functions
  !USE sep_fields

  IMPLICIT NONE
  
  PRIVATE
  PUBLIC :: L2DCGRID, L2DCINIFIELDS, L2DFIELDS

  CONTAINS 
!-------------------------------------------------------------------  
SUBROUTINE L2DCGRID
! subroutine to read in ONLY the Lare grid of a *.cfd file
! features several horrible hacks of the lare3d code. Sorry Tony Arber!

 !CHARACTER(len=7) 				:: fmt1='(I4.4)', istring 	! format descriptor
 !CHARACTER(LEN = 20+data_dir_max_length) 	:: cfdloc
 CHARACTER(LEN = 20) :: name, class, mesh_name, mesh_class
 !INTEGER, DIMENSION(3) :: dims
 REAL(num) :: time_d
 !REAL(num), DIMENSION(3) :: extent
 !REAL(num), DIMENSION(3) :: stagger
 INTEGER:: nblocks, type, nd, sof, snap
 
 ALLOCATE(extent(2), stagger(2))
  
!  WRITE (istring,fmt1) mysnap 			! converting integer to string using an 'internal file'
!  cfdloc=trim(adjustl(sloc))//trim(istring)//filetype1		! store new filename
  
  print*, 'reading grid from:', cfdloc
  
! BEGIN HACK OF LAREXD READ:  

   !CALL cfd_open(cfdloc, rank, comm, MPI_MODE_RDONLY)
   cfd_comm = comm
   cfd_rank = rank
   cfd_mode = MPI_MODE_RDONLY
   CALL cfd_open_read(cfdloc)
   nblocks = cfd_get_nblocks()

    DO ix = 1, nblocks
      CALL cfd_get_next_block_info_all(name, class, type)
      !IF (rank == 0) PRINT *, ix,nblocks, name, class, type
      !print*, name
      
      IF (type == TYPE_SNAPSHOT) THEN	!if type=3
        CALL cfd_get_snapshot(time_d, snap)
        time = time_d
      END IF

      IF (type == TYPE_MESH) THEN	!if type=1

        CALL cfd_get_2d_cartesian_grid_all(myx, myy)

        !CALL cfd_skip_block()
      ELSE IF (type == TYPE_MESH_VARIABLE) THEN		!if type=2
        CALL cfd_get_common_meshtype_metadata_all(type, nd, sof)

        IF (nd /= DIMENSION_2D) THEN
          IF (rank == 0) PRINT *, "Non 2D Dataset found in input file, ", &
              "ignoring and continuing."
          CALL cfd_skip_block()
          CYCLE
        END IF

        IF (type /= VAR_CARTESIAN) THEN
          IF (rank == 0) PRINT *, "Non - Cartesian variable block found ", &
              "in file, ignoring and continuing"
          CALL cfd_skip_block()
          CYCLE
        END IF

        ! Should be at end of block, but force the point anyway
        CALL cfd_skip_block()
      ELSE
        ! Unknown block, just skip it
        CALL cfd_skip_block()
      END IF
    END DO

    CALL cfd_close()
    CALL MPI_BARRIER(comm, errcode)
    DEALLOCATE(extent, stagger)


END SUBROUTINE L2DCGRID
!------------------------------------------------------------------------------! 
SUBROUTINE L2DCINIFIELDS
! designed to read in ALL field variables from Lare3d cfd files. Apologies to T. Arber.

    INTEGER				:: ii
    CHARACTER(LEN = 20) 		:: name, class, mesh_name, mesh_class
    !INTEGER				:: nblocks, type, nd, sof, snap
    INTEGER				:: nblocks, type, nd, sof, snap
    REAL(num) :: time_d
    REAL(num), DIMENSION(:, :), ALLOCATABLE :: data
  
   !ALLOCATE(extent(2), stagger(2))
 
 ! WRITE(istring,fmt1) mysnap 			! converting integer to string using an 'internal file'
 ! cfdloc=trim(adjustl(sloc))//trim(istring)//filetype1		! store new filename

  !print*, 'adjustl(sloc):', adjustl(sloc)
  !print*, 'trim(istring):', trim(istring)
  
! BEGIN HACK OF LAREXD READ:
  PRINT*, 'reading fields from: ', cfdloc   

   ALLOCATE(data(0:nx, 0:ny))

   cfd_comm = comm
   cfd_rank = rank
   cfd_mode = MPI_MODE_RDONLY
   !CALL cfd_open_read(cfdloc)
   CALL cfd_open(cfdloc, rank, comm, MPI_MODE_RDONLY)
   
   nblocks = cfd_get_nblocks()

   !print*, 'nblocks= ', nblocks
   !print*, nx_global
    DO ix = 1, nblocks
      CALL cfd_get_next_block_info_all(name, class, type)
      !IF (rank == 0) PRINT *, ix,nblocks, name, class, type
      !print*, name
      
      IF (type == TYPE_SNAPSHOT) THEN	!if type=3
        CALL cfd_get_snapshot(time_d, snap)
        time = time_d
	ltimes(frame)=time_d
      END IF

      IF (type == TYPE_MESH) THEN	!if type=1
        ! Strangely, LARE doesn't actually read in the grid from a file
        ! This can be fixed, but for the moment, just go with the flow and
        ! Replicate the old behaviour
 
 	!cfd_get_2d_cartesian_variable_parallel(data, subtype)
        !CALL cfd_get_3d_cartesian_grid_all(myx, myy, myz)
	
        CALL cfd_skip_block()
      ELSE IF (type == TYPE_MESH_VARIABLE) THEN		!if type=2
        CALL cfd_get_common_meshtype_metadata_all(type, nd, sof)
			
	IF (nd /= DIMENSION_2D) THEN
          IF (rank == 0) PRINT *, "Non 2D Dataset found in input file, ", &
              "ignoring and continuting."
          CALL cfd_skip_block()
          CYCLE
        END IF

        IF (type /= VAR_CARTESIAN) THEN
          IF (rank == 0) PRINT *, "Non - Cartesian variable block found ", &
              "in file, ignoring and continuing"
          CALL cfd_skip_block()
          CYCLE
        END IF

        ! We now have a valid variable, let's load it up
        ! First error trapping
        CALL cfd_get_nd_cartesian_variable_metadata_all(nd, dims, extent, &
            stagger, mesh_name, mesh_class)

	IF (dims(1) /= nx_global+1 .OR. dims(2) /= ny_global+1) THEN
          IF (rank == 0) PRINT *, "Size of grid represented by one more ", &
              "variables invalid. Continuing"
          CALL cfd_skip_block
          CYCLE
        END IF

        IF (sof /= num) THEN
          IF (rank == 0) PRINT *, "Precision of data does not match ", &
              "precision of code. Continuing."
          CALL cfd_skip_block
        END IF

        ! We're not interested in the other parameters, so if we're here,
        ! load up the data

	CALL cfd_get_2d_cartesian_variable_parallel(data, subtype)

        ! Now have the data, just copy it to correct place
        	
	!IF (str_cmp(name(1:3), "Rho")) THEN
        !  rho(0:nx, 0:ny, 0:nz) = data
        !END IF
        !IF (str_cmp(name(1:6), "Energy")) THEN
        !  energy(0:nx, 0:ny, 0:nz) = data
        !END IF

        IF (str_cmp(name(1:2), "Vx")) THEN
 	 PRINT *, ix,nblocks, name, class
          vx(0:nx, 0:ny, 0,frame) = data
        END IF
        IF (str_cmp(name(1:2), "Vy")) THEN
 	 PRINT *, ix,nblocks, name, class
          vy(0:nx, 0:ny, 0,frame) = data
        END IF
        IF (str_cmp(name(1:2), "Vz")) THEN
	 PRINT *, ix,nblocks, name, class
          vz(0:nx, 0:ny, 0,frame) = data
        END IF

        IF (str_cmp(name(1:2), "Bx")) THEN
	PRINT *, ix,nblocks, name, class
	! destagger on the fly to same locations as vx, vy, vz
        !  bx(0:nx, 0:ny, 0:nz,frame) = data
	  DO ii=0,nx
	   !bx(ii,0:ny,0:nz)=stagger_right(stagger_up(data(ii,0:ny,0:nz)))
	   bx(ii,0:ny,0,frame)=stagger_bx_2d(data(ii,0:ny))
	  ENDDO  		
	END IF
        IF (str_cmp(name(1:2), "By")) THEN
	PRINT *, ix,nblocks, name, class
	! destagger on the fly to same locations as vx, vy, vz
        !  by(0:nx, 0:ny, 0:nz) = data
	!print*, 'BY started' 
	 DO ii=0,ny
	  by(0:nx,ii,0,frame)=stagger_by_2d(data(0:nx,ii))
	 ENDDO
	! print*, 'BY finished' 	 
        END IF
        IF (str_cmp(name(1:2), "Bz")) THEN
	PRINT *, ix,nblocks, name, class
	!PRINT*, 'BZ reached'
	! destagger on the fly to same locations as vx, vy, vz
        !   bz(0:nx, 0:ny, 0:nz) = data
	  bz(0:nx,0:ny,0,frame)=stagger_bz(data(0:nx,0:ny))
        END IF
	
	!IF (str_cmp(name(1:11), "Temperature")) THEN
        !  temperature(0:nx, 0:ny, 0:nz) = data
        !END IF
	
	!IF (str_cmp(name(1:8), "Pressure")) THEN
        !  pressure(0:nx, 0:ny, 0:nz) = data
        !END IF
	
	!IF (str_cmp(name(1:3), "eta")) THEN
        !  eta(0:nx, 0:ny, 0:nz) = data
        !END IF

        ! Should be at end of block, but force the point anyway
        CALL cfd_skip_block()
      ELSE
        ! Unknown block, just skip it
	
	!print*, name(1:8)
	!IF (str_cmp(name(1:8), "Pressure")) THEN
        !  pressure(0:nx, 0:ny, 0:nz) = data
        !END IF
	
        CALL cfd_skip_block()
      END IF
    END DO

    DEALLOCATE(data)
    CALL cfd_close()
    CALL MPI_BARRIER(comm, errcode)

END SUBROUTINE L2DCINIFIELDS
!----------------------------------------------------------------------;   
SUBROUTINE L2DFIELDS(R,T,E,B,DBDX,DBDY,DBDZ,DBDT,DEDX,DEDY,DEDZ,DEDT,Vf)
! wrapper module which calls a routine in lare_functions which interprets
! values at sub grid particle locations, and converts these to correctly named variables.

 REAL(num), DIMENSION(3),INTENT(OUT)	:: B,E
 REAL(num), DIMENSION(3),INTENT(OUT)	:: DBDX,DBDY,DBDZ,DBDT,DEDX,DEDY,DEDZ,DEDT
 REAL(num), DIMENSION(3),INTENT(IN)	:: R
 REAL(num), DIMENSION(3)  		:: Vf, j
 REAL(num), INTENT(IN) 			:: T
 REAL(num), DIMENSION(36)		:: iquants
 INTEGER				:: ij, jjx, jjy, jjz
 LOGICAL				:: fxflag=.FALSE., fyflag=.FALSE.!, fzflag=.FALSE.
 !LOGICAL, INTENT(OUT)			:: ERR=.FALSE.



  IF ((R(1).lt.myx(4)).and.(R(1).ge.myx(nx-3))) fxflag=.TRUE.
  IF ((R(2).lt.myy(4)).and.(R(2).ge.myy(ny-3))) fyflag=.TRUE.
  !IF ((R(3).lt.myz(4)).and.(R(3).ge.myz(nz-3))) fzflag=.TRUE.
 
 ! I don't think this is an actual error since going beyond the grid will be detected and these results ignored
  IF ((fxflag).OR.(fyflag)) THEN		
   ! ERR=.TRUE.
  !   PRINT*, 'particle not located on the grid:'
  !   PRINT*, 'R:', R
  !   PRINT*, 'x:', myx(4),'->', myx(nx-3)
  !   PRINT*, 'y:', myy(4),'->', myy(ny-3)
  !   PRINT*, 'z:', myz(4),'->', myz(nz-3)
  !   STOP
  !   PRINT*, '---'
   RETURN
  ENDIF

   iquants=T2d(R,T)
   ! compare array against itself for NANs         
   IF (ALL(iquants.eq.iquants)) THEN
   !everything is fine
   ELSE
    PRINT*, 'NAN DETECTED IN LARE OUTPUT'
    !ERR=.TRUE.
    STOP
   ENDIF
   
   B	=iquants(1:3)
   Vf	=iquants(4:6)
   E	=iquants(7:9)
   j	=iquants(10:12)
   DBDX	=iquants(13:15)
   DBDY	=iquants(16:18)
   DBDZ	=iquants(19:21)
   DEDX	=iquants(22:24)
   DEDY	=iquants(25:27)
   DEDZ	=iquants(28:30)
   DBDT	=iquants(31:33)
   DEDT	=iquants(34:36)  
   

END SUBROUTINE L2DFIELDS

END MODULE l2dc_fields

