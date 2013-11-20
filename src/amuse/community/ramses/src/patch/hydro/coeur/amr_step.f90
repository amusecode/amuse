recursive subroutine amr_step(ilevel,icount)
  use amr_commons
  use pm_commons
  use hydro_commons
  use poisson_commons
  implicit none
  integer::ilevel,icount
  !-------------------------------------------------------------------!
  ! This routine is the adaptive-mesh/adaptive-time-step main driver. !
  ! Each routine is called using a specific order, don't change it,   !
  ! unless you check all consequences first                           !
  !-------------------------------------------------------------------!
  integer::icycle,i,idim,ivar
  if(numbtot(1,ilevel)==0)return

  if(verbose)write(*,999)icount,ilevel

  !-------------------------------------------
  ! Make new refinements and update boundaries
  !-------------------------------------------
  if(ilevel==levelmin.or.icount>1)then
     do i=ilevel,nlevelmax
        if(i>levelmin)then
           call build_comm(i) ! Build communicators
           call make_virtual_fine_int(cpu_map(1),i)
           if(hydro)then
              do ivar=1,nvar
                 call make_virtual_fine_dp(uold(1,ivar),i)
              end do
              if(simple_boundary)then
                 call make_boundary_hydro(i)
              end if
              if(poisson)then
                 do idim=1,ndim
                    call make_virtual_fine_dp(f(1,idim),i)
                 end do
              end if
           end if
        end if
        call refine_fine(i)
     end do
     ! Load balance
     if(ilevel==levelmin.and.nremap>0)then
        if(MOD(nstep_coarse,nremap)==0)then
           call load_balance
        end if
     end if
  end if

  !---------------
  ! Gravity update
  !---------------
  if(poisson)then

     ! Update ilevel particle tree to account for particle leakage
     if(pic)call make_tree_fine(ilevel)

     ! Compute Poisson source term
     call rho_fine(ilevel,icount)

     ! Synchronize hydro for gravity (first pass)
     if(hydro)call synchro_hydro_fine(ilevel,-0.5*dtnew(ilevel))

     ! Compute gravitational potential
     if(ilevel>levelmin)then
        call phi_fine_cg(ilevel)
     else
        call full_multigrid
     end if

     ! Compute gravitational acceleration
     call force_fine(ilevel)

     ! Synchronize hydro for gravity (second pass)
     if(hydro)then
        call synchro_hydro_fine(ilevel,+0.5*dtnew(ilevel))
        ! Update boundaries
        do ivar=1,nvar
           call make_virtual_fine_dp(uold(1,ivar),ilevel)
        end do
     end if

     ! Sort particles between ilevel and ilevel+1
     if(pic)then
        ! Remove particles to finer levels
        call kill_tree_fine(ilevel)
        ! Update boundary conditions for remaining particles
        call virtual_tree_fine(ilevel)
        ! Synchronize remaining particles
        call synchro_fine(ilevel)
     end if

  end if

  if(poisson) call barotrop(ilevel)
  !----------------------
  ! Compute new time step
  !----------------------
  call newdt_fine(ilevel)
  if(ilevel>levelmin)then
     dtnew(ilevel)=MIN(dtnew(ilevel-1)/real(nsubcycle(ilevel-1)),dtnew(ilevel))
  end if

  ! Set unew equal to uold
  if(hydro)call set_unew(ilevel)

  !---------------------------
  ! Recursive call to amr_step
  !---------------------------
  if(ilevel<nlevelmax)then
     if(numbtot(1,ilevel+1)>0)then
        if(nsubcycle(ilevel)==2)then
           call amr_step(ilevel+1,1)
           call amr_step(ilevel+1,2)
        else
           call amr_step(ilevel+1,1)
        endif
     else 
        ! Otherwise, update time and finer level time-step
        dtold(ilevel+1)=dtnew(ilevel)/dble(nsubcycle(ilevel))
        dtnew(ilevel+1)=dtnew(ilevel)/dble(nsubcycle(ilevel))
        call update_time(ilevel)
     end if
  else
     call update_time(ilevel)
  end if

  !-----------
  ! Hydro step
  !-----------
  if(hydro)then
     if(star)call star_formation(ilevel)                     ! Star formation
     call godunov_fine(ilevel)                            ! Hyperbolic solver
     call set_uold(ilevel)                           ! Set uold equal to unew
     if(poisson)call synchro_hydro_fine(ilevel,dtnew(ilevel))! Gravity source
     if(cooling)call cooling_fine(ilevel)                    ! Cooling source
     call upload_fine(ilevel)                                   ! Restriction
     if(poisson) call barotrop(ilevel)
     do ivar=1,nvar                                       ! Update boundaries
        call make_virtual_fine_dp(uold(1,ivar),ilevel)    
     end do
     if(simple_boundary)call make_boundary_hydro(ilevel)
  end if

  !-----------------------
  ! Compute refinement map
  !-----------------------
  call flag_fine(ilevel,icount)

  !---------------
  ! Move particles
  !---------------
  if(pic)then
     ! Move only remaining particles
     call move_fine(ilevel)
     ! Get particles from all finer levels
     call merge_tree_fine(ilevel)
  end if

  !-------------------------------
  ! Update coarser level time-step
  !-------------------------------
  if(ilevel>levelmin)then
     if(nsubcycle(ilevel-1)==1)dtnew(ilevel-1)=dtnew(ilevel)
     if(icount==2)dtnew(ilevel-1)=dtold(ilevel)+dtnew(ilevel)
  end if

999 format(' Entering amr_step',i1,' for level',i2)

end subroutine amr_step



