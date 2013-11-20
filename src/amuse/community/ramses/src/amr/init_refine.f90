!################################################################
!################################################################
!################################################################
!################################################################
subroutine init_refine
  use amr_commons
  use pm_commons
  implicit none
  !-------------------------------------------
  ! This routine builds the initial AMR grid
  !-------------------------------------------
  integer::ilevel

  if(myid==1)write(*,*)'Building initial AMR grid'
  init=.true.

  ! Base refinement
  do ilevel=1,levelmin
     call flag
     call refine
  end do

  ! Further refinements if necessary
  do ilevel=levelmin+1,nlevelmax
     if(initfile(levelmin).ne.' '.and.initfile(ilevel).eq.' ')exit
     if(hydro)call init_flow
#ifdef RT
     if(rt)call rt_init_flow
#endif
     if(ivar_refine==0)call init_refmap
     call flag
     call refine
     if(nremap>0)call load_balance
     if(numbtot(1,ilevel)==0)exit
  end do 

  ! Final pass to initialize the flow
  init=.false.
  if(hydro)call init_flow
#ifdef RT
  if(rt)call rt_init_flow
#endif

end subroutine init_refine
!################################################################
!################################################################
!################################################################
!################################################################
subroutine init_refine_2
  !--------------------------------------------------------------
  ! This routine builds additional refinements to the
  ! the initial AMR grid for filetype ne 'grafic'
  !--------------------------------------------------------------
  use amr_commons
  use hydro_commons
#ifdef RT
  use rt_hydro_commons
#endif
  use pm_commons
  use poisson_commons
  implicit none
  integer::ilevel,i,ivar

  if(filetype.eq.'grafic')return

  do i=levelmin,nlevelmax+1

     call refine_coarse
     do ilevel=1,nlevelmax
        call build_comm(ilevel)
        call make_virtual_fine_int(cpu_map(1),ilevel)
        call refine_fine(ilevel)
        if(hydro)call init_flow_fine(ilevel)
#ifdef RT
        if(rt)call rt_init_flow_fine(ilevel)
#endif
     end do

     if(nremap>0)call load_balance

     do ilevel=levelmin,nlevelmax
        if(pic)call make_tree_fine(ilevel)
        if(poisson)call rho_fine(ilevel,2)
        if(pic)then
           call kill_tree_fine(ilevel)
           call virtual_tree_fine(ilevel)
        endif
     end do

     do ilevel=nlevelmax,levelmin,-1
        if(pic)call merge_tree_fine(ilevel)
        if(hydro)then
           call upload_fine(ilevel)
#ifdef SOLVERmhd
           do ivar=1,nvar+3
#else
           do ivar=1,nvar
#endif
              call make_virtual_fine_dp(uold(1,ivar),ilevel)
#ifdef SOLVERmhd
           end do
#else
           end do
#endif
           if(simple_boundary)call make_boundary_hydro(ilevel)
        endif
#ifdef RT
        if(rt)then
           call rt_upload_fine(ilevel)
           do ivar=1,nrtvar
              call make_virtual_fine_dp(rtuold(1,ivar),ilevel)
           end do
           if(simple_boundary)call rt_make_boundary_hydro(ilevel)
        end if
#endif
     end do

     do ilevel=nlevelmax,1,-1
        call flag_fine(ilevel,2)
     end do
     call flag_coarse

  end do

#ifdef RT
  if(rt_is_init_xion .and. rt_nregion .eq. 0) then
     if(myid==1) write(*,*) 'Initializing ionization states from T profile'
     do ilevel=nlevelmax,1,-1
        call rt_init_xion(ilevel)
        call upload_fine(ilevel)
     end do
  endif
#endif  

end subroutine init_refine_2
!################################################################
!################################################################
!################################################################
!################################################################
