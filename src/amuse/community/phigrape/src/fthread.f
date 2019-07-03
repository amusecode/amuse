! subroutine to start thread in fortran
! arguments: actie: action to be taken,
! 	    
!
! pt_setup is a setup routine (eg to copy data to buffer)
! viewbodies is the actual sub-program
!
! needs pt lib
	subroutine fthread(actie)
	implicit none
	character(len=5),intent(in) ::  actie
	external pipeline_init__,pipeline_done__,pipeline_execute__
	external pt_setup, viewbodies

	select case( actie)

	case('setup')
	call pipeline_init_(pt_setup, viewbodies)
	print*, 'setup thread'
	
	case('execu')
	call pipeline_execute_()
	print*, 'execute thread'
	
	case('tstop')
	call pipeline_done_()
	print*, 'stop thread'
	
	case default
	print*,'<fthread> input error'
	
	end select

	return

	end
	
	SUBROUTINE viewer
	LOGICAL initcall
	DATA initcall/.TRUE./
	SAVE initcall
	
	IF(initcall) THEN
!	CALL fthread('setup',viewbodies)
!	CALL fthread('execu',viewbodies)
	CALL fthread('setup')
	CALL fthread('execu')
	ELSE
!	call fthread('tstop',viewbodies)
	call fthread('tstop')
	ENDIF
	initcall=.NOT.initcall
	
	END
