!> \file amuse_mainloop.F90

!> \brief the main program loop
!! 
!! Updates one snapshot
!<

module amuse_mainloop_mod
  use myf03_mod

  ! routines
  use gadget_general_class, only: gadget_constants_type
  use main_input_mod, only: readin_snapshot
  use oct_tree_mod, only: buildtree, setparticleorder
  use ray_mod, only: src_ray_make
  use raylist_mod, only: trace_ray, kill_raylist, prepare_raysearch
  use ion_temperature_update, only: update_raylist, non_photo_update_all
  use mt19937_mod, only: genrand_real1
  use output_mod, only: output_total_snap, ion_frac_out
  use global_mod, only: set_dt_from_dtcode, set_time_elapsed_from_itime

  ! types
  use ray_mod, only: src_ray_type
  
  ! variables
  use global_mod, only: psys
  use global_mod, only: globalraylist
  use global_mod, only: tree
  use global_mod, only: GV
  use global_mod, only: PLAN

  implicit none
  

contains
  
  !> this is the main driver of SPHRAY
  !======================================
  subroutine mainloop()
    implicit none
    
    character(clen), parameter :: myname="mainloop"
    logical, parameter :: crash=.true.
    integer, parameter :: verb=1
    character(clen) :: str,fmt

    type(gadget_constants_type) :: gconst

    !  local counters 
    !-----------------
    
    integer(i8b) :: snapn !< snapshot index
    integer(i8b) :: rayn  !< ray counter
    integer(i8b) :: srcn  !< source counter
    
    real(r8b) :: rn       !< random number
    real(r8b) :: MB       !< MBs for memory consumption tracking
    
    ! work variables
    !-----------------
    logical :: srcray
    type(src_ray_type) :: ray
        
#ifdef incHrec
    integer(i8b) :: rindx
    integer(i8b) :: pindx
    integer(i8b) :: lasthitcheck
#endif
        

    snapn = GV%StartSnapNum

    
    ! if there are no sources, update the gas cooling and
    ! ionization state with zero photoionizations and using
    ! the same time step for each particle
    !----------------------------------------------------------------
    if (size(psys%src).EQ.0) then
       
       GV%itime = GV%itime + PLAN%snap(snapn)%SrcRays
       call set_time_elapsed_from_itime( GV )
       write(*,*) PLAN%snap(snapn)%SrcRays, &
                  PLAN%snap(snapn)%SrcRays, &
                  GV%time_elapsed_myr



       
    ! if there are sources, calculate the photoionization rate
    ! at each particle using ray tracing and update on individual
    ! time steps.  
    !----------------------------------------------------------------
    else   
                       
       !  build oct tree.  only need to do this once per snap (for now)
       !----------------------------------------------------------------
       call buildtree(psys,tree,MB,GV%PartPerCell)
       call setparticleorder(psys, tree)             
       call prepare_raysearch(psys, globalraylist)
       
              
       ! begin ray tracing 
       !----------------------------------------------------------------
       src_rays: do rayn = 1, PLAN%snap(snapn)%SrcRays
          
          GV%rayn                = GV%rayn + 1
          GV%src_rayn            = GV%src_rayn + 1
          GV%TotalSourceRaysCast = GV%TotalSourceRaysCast + 1                
          
          !  select a source randomly (weighted by their luminosity)
          !----------------------------------------------------------------
          rn = genrand_real1() * psys%src(size(psys%src))%Lcdf
          srcn=1
          do while( psys%src(srcn)%Lcdf < rn )
             srcn = srcn + 1
             if( srcn > size(psys%src) ) then
                stop "src num > number of sources in mainloop.f90"
             endif
          enddo
          
          !  create a source ray and calc the impacts
          !----------------------------------------------------------------
          call src_ray_make( ray, &
                             psys%src(srcn), &
                             GV%rayn, &
                             GV%dt_s, &
                             GV%Lunit, &
                             psys%box )
          
          GV%itime = GV%itime + 1
          call set_time_elapsed_from_itime( GV )
          GV%IonizingPhotonsPerSec = GV%TotalPhotonsCast / GV%time_elapsed_s
          
          if ( mod(rayn,PLAN%snap(snapn)%SrcRays/100) == 0 ) &
               write(*,*) rayn, PLAN%snap(snapn)%SrcRays, GV%time_elapsed_myr
                    
          globalraylist%ray = ray
          call trace_ray(globalraylist%ray, globalraylist, psys, tree) 
          
          GV%TotalPhotonsCast = GV%TotalPhotonsCast + globalraylist%ray%pini
          
          srcray = .true.
          call update_raylist(globalraylist, psys%par, psys%box, srcray)
          
                    
       end do src_rays
       
       ! free up the memory from the globalraylist.
       call kill_raylist(globalraylist)
       
    endif


    ! at this point particles are synchronized.  if there are sources, 
    ! particles are evolved with GHI=0 for a dt equal to the last time
    ! they were hit by a ray to the final time.  if there are no sources
    ! each particle is updated over the whole time
    call non_photo_update_all(psys%par)
    
   
    close(GV%ionlun)
    
    
    if (GV%raystats) then
       close(GV%raystatlun)
    end if
    
  end subroutine mainloop
  
end module amuse_mainloop_mod
