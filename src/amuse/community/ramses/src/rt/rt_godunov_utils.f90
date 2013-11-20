!###########################################################
!###########################################################
!###########################################################                                     
SUBROUTINE get_rt_courant_coarse(dt)

! Determine the coarse RT timestep length set by the Courant condition
!-------------------------------------------------------------------------
  use amr_parameters
  use rt_parameters
  implicit none
  integer:: nx_loc
  real(dp):: dt, scale, dx
!-------------------------------------------------------------------------
  ! Mesh spacing at coarse level
  nx_loc=icoarse_max-icoarse_min+1
  scale=boxlen/dble(nx_loc)
  dx=0.5D0**levelmin*scale
  dt = rt_courant_factor*dx/3.d0/rt_c
END SUBROUTINE get_rt_courant_coarse
!###########################################################
!###########################################################
!###########################################################
!###########################################################
subroutine rt_hydro_refine(ug,um,ud,ok,nn)
  use amr_parameters
  use const
  use rt_parameters
  implicit none
  ! dummy arguments
  integer nn
  real(dp)::ug(1:nvector,1:nrtvar)
  real(dp)::um(1:nvector,1:nrtvar)
  real(dp)::ud(1:nvector,1:nrtvar)
  logical ::ok(1:nvector)
  
  integer::k,i
  real(dp)::dg,dm,dd,pg,pm,pd,vg,vm,vd,cg,cm,cd,error
  
  if(rt .and. rt_err_grad_n >= 0.) then !---------------------------------                           
     do i=1,nGroups                                                 
        ! RT-photon density                                       
        do k=1,nn                                                 
           dg=ug(k,iGroups(i)); dm=um(k,iGroups(i)); dd=ud(k,iGroups(i))   
           error=2.0d0*MAX( &                                    
                & ABS((dd-dm)/(dd+dm+rt_floor_n)) , &            
                & ABS((dm-dg)/(dm+dg+rt_floor_n)) )              
           ok(k) = ok(k) .or. error > rt_err_grad_n              
        end do                                                   
     end do                                                      
  end if                                                         
                                                                  
end subroutine rt_hydro_refine
!###########################################################
!###########################################################
!###########################################################
!###########################################################
  


