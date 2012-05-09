!> \file mainloop.F90

!> \brief the main program loop
!! 
!! Loops through the snapshots calling the raytracing and output routines 
!<

module amuse_mainloop_mod
  use myf03_mod

  ! routines
  use gadget_general_class, only: gadget_constants_type
  use main_input_mod, only: readin_snapshot
  use oct_tree_mod, only: buildtree, setparticleorder
  use ray_mod
  use raylist_mod
  use ion_temperature_update, only: update_raylist, non_photo_update_all
  use mt19937_mod, only: genrand_real1
  use output_mod, only: output_total_snap, ion_frac_out
  use global_mod, only: set_dt_from_dtcode, set_time_elapsed_from_itime
  
  ! variables
  use global_mod, only: psys
  use global_mod, only: globalraylist
  use global_mod, only: tree
  use global_mod, only: GV
  use global_mod, only: PLAN


  
  implicit none
  
  integer(i8b), parameter :: one = 1

  integer(i4b), parameter :: rays_per_leaf = 2
  integer(i4b), parameter :: rays_per_dt = 1
  
contains
  
  !> this is the main driver of SPHRAY
  !======================================
  subroutine mainloop()
    implicit none
    
    character(clen), parameter :: myname="mainloop"
    logical, parameter :: crash=.true.
    integer, parameter :: verb=1
    character(clen) :: str,fmt

    type(raystat_type) :: raystats(raystatbuffsize)
    integer(i8b) :: raystatcnt
    type(gadget_constants_type) :: gconst

    !  local counters 
    !-----------------
    
    integer(i8b) :: snapn !< snapshot counter
    integer(i8b) :: rayn  !< ray counter
    integer(i8b) :: srcn  !< source counter
    
    real(r8b) :: rn       !< random number
    real(r8b) :: MB       !< MBs for memory consumption tracking
    
    ! work variables
    !-----------------
    logical :: srcray
    type(src_ray_type) :: ray
    real(r8b) :: outmark
    
#ifdef incHrec
    integer(i8b) :: rindx
    integer(i8b) :: pindx
    integer(i8b) :: lasthitcheck
#endif
    
    
    raystatcnt = 0
    
    ! loop over the snapshots 
    !=========================
    snaps: do snapn = GV%StartSnapNum, GV%EndSnapNum
    
       if(size(psys%src).EQ.0) then

          GV%itime = GV%itime + PLAN%snap(snapn)%SrcRays
          call set_time_elapsed_from_itime( GV )
          print*, PLAN%snap(snapn)%SrcRays, PLAN%snap(snapn)%SrcRays,GV%time_elapsed_myr
          goto 123
       endif   
         
              
       !  build oct tree.  only need to do this once per snap (for now)
       !----------------------------------------------------------------
       call buildtree(psys,tree,MB,GV%PartPerCell)
       GV%MB = GV%MB + MB
       call setparticleorder(psys, tree)             
       call prepare_raysearch(psys, globalraylist)
             
       if (GV%raystats) then
          write(GV%raystatlun) PLAN%snap(snapn)%SrcRays, raystatbuffsize 
       endif
       
       if(GV%JustInit) then
          write(str,"(A,F10.2)") "total memory allocated [MB] = ", GV%MB
          call mywrite(str,verb)
          call mywrite("just initializing", verb)
          call mywrite("",verb)
          stop
       end if

       if (GV%DoInitialOutput) then
          GV%OutputIndx = 0
          call output_total_snap(psys)      
          GV%OutputIndx = 1
       end if
              
       ! begin ray tracing 
       !------------------------- 
       src_rays: do rayn = one, PLAN%snap(snapn)%SrcRays

          GV%rayn                = GV%rayn + 1
          GV%src_rayn            = GV%src_rayn + 1
          GV%TotalSourceRaysCast = GV%TotalSourceRaysCast + 1                
          
          !  select a source randomly (weighted by their luminosity)
          rn = genrand_real1() * psys%src(size(psys%src))%Lcdf
          srcn=1
          do while(psys%src(srcn)%Lcdf.LT.rn)
             srcn=srcn+1
             if(srcn.GT.size(psys%src)) then
                write(*,*) srcn, rn, psys%src(size(psys%src))%Lcdf, size(psys%src)
                stop "src num > number of sources in mainloop.f90"
             endif
          enddo
                    
          !  create a source ray and calc the impacts
          call src_ray_make( ray, psys%src(srcn), GV%rayn, GV%dt_s, GV%Lunit, psys%box )
          
          GV%itime = GV%itime + 1
          call set_time_elapsed_from_itime( GV )
          GV%IonizingPhotonsPerSec = GV%TotalPhotonsCast / GV%time_elapsed_s

          if(mod(rayn,PLAN%snap(snapn)%SrcRays/100).EQ.0) &
              print*, rayn, PLAN%snap(snapn)%SrcRays,GV%time_elapsed_myr


          globalraylist%ray = ray
          call trace_ray(globalraylist%ray, globalraylist, psys, tree) 
          
          GV%TotalPhotonsCast = GV%TotalPhotonsCast + globalraylist%ray%pini
          
          srcray = .true.
          call update_raylist(globalraylist,psys%par,psys%box,srcray)
                    
          if (GV%raystats) then
             
             raystatcnt = raystatcnt + 1
             
             raystats(raystatcnt)%srcn  = srcn
             raystats(raystatcnt)%start = globalraylist%ray%start  
             raystats(raystatcnt)%ryd   = globalraylist%ray%freq
             
             if (raystatcnt == raystatbuffsize) then
                write(GV%raystatlun) raystats
                flush(GV%raystatlun)
                raystatcnt = 0
             end if
                          
          end if
          
                              
          
          !        output routines
          !------------------------
          
          !        check if this time step requires a small ionization frac output
          !        these are small outputs done often to monitor certain quantites
!          if( mod(GV%rayn,GV%IonFracOutRays)==0 ) then
!             call ion_frac_out(psys, tree )
!          end if
          
          
          ! check if this time step requires a full output
          if ( GV%OutputIndx <= GV%NumTotOuts ) then
             
             ! set correct time marker unit
             if ( trim(GV%OutputTiming) == "standard" ) then
                
                outmark = GV%start_time_code + GV%time_elapsed_code
                
             else if ( trim(GV%OutputTiming) == "forced" ) then
                
                if (trim(GV%ForcedUnits) == "mwionfrac") then
                   outmark = GV%mwionfrac
                else 
                   outmark = GV%time_elapsed_code
                end if
                
             else
                
                write(*,*) "output type ", trim(GV%OutputTiming), "not recognized"
                stop 
                
             end if


             
             ! check outmark against tabulated output "times"
             if ( outmark >= PLAN%OutputTimes(GV%OutputIndx) ) then
                call output_total_snap(psys)
                GV%OutputIndx = GV%OutputIndx + 1 
             end if
             
          end if
          
          
          ! if we are on the last ray and we havent gotten to the last 
          ! output, do a full output.
          if ( snapn == GV%EndSnapNum ) then
             if ( GV%OutputIndx <= GV%NumTotOuts ) then
                if ( GV%src_rayn==PLAN%snap(snapn)%SrcRays ) then
                   write(*,*) "doing an output on the last ray"
                   call output_total_snap(psys)
                end if
             end if
          end if
          
          
       end do src_rays
       
       ! free up the memory from the globalraylist.
       call kill_raylist(globalraylist)

123 call non_photo_update_all(psys%par)

    end do snaps
    
    close(GV%ionlun)
    
    
    if (GV%raystats) then
       close(GV%raystatlun)
    end if
        
  end subroutine mainloop
  
end module amuse_mainloop_mod
