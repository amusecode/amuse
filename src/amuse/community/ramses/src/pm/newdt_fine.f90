subroutine newdt_fine(ilevel)
  use pm_commons
  use amr_commons
  use hydro_commons
  use poisson_commons, ONLY: gravity_type
#ifdef RT
  use rt_parameters, ONLY: rt_advect
#endif
  implicit none
#ifndef WITHOUTMPI
  include 'mpif.h'
#endif
  integer::ilevel
  !-----------------------------------------------------------
  ! This routine compute the time step using 3 constraints:
  ! 1- a Courant-type condition using particle velocity
  ! 2- the gravity free-fall time
  ! 3- 10% maximum variation for aexp 
  ! 4- maximum step time for ATON
  ! This routine also compute the particle kinetic energy.
  !-----------------------------------------------------------
  integer::igrid,jgrid,ipart,jpart
  integer::npart1,ip,info
  integer,dimension(1:nvector),save::ind_part
  real(kind=8)::dt_loc,dt_all,ekin_loc,ekin_all
  real(dp)::tff,fourpi,threepi2
  real(dp)::aton_time_step,dt_aton,dt_rt

  if(numbtot(1,ilevel)==0)return
  if(verbose)write(*,111)ilevel

  ! Save old time step
  dtold(ilevel)=dtnew(ilevel)

  ! Maximum time step
  dtnew(ilevel)=boxlen/smallc
  if(poisson.and.gravity_type<=0)then
     fourpi=4.0d0*ACOS(-1.0d0)
     threepi2=3.0d0*ACOS(-1.0d0)**2
     if(cosmo)fourpi=1.5d0*omega_m*aexp
     tff=sqrt(threepi2/8./fourpi/rho_max(ilevel))
     dtnew(ilevel)=MIN(dtnew(ilevel),courant_factor*tff)
  end if
  if(cosmo)then
     dtnew(ilevel)=MIN(dtnew(ilevel),0.1/hexp)
  end if

#ifdef ATON
  ! Maximum time step for ATON
  if(aton)then
     dt_aton = aton_time_step()
     if(dt_aton>0d0)then
        dtnew(ilevel)=MIN(dtnew(ilevel),dt_aton)
     end if
  end if
#endif

#ifdef RT
  ! Maximum time step for radiative transfer
  if(rt_advect)then
     call get_rt_courant_coarse(dt_rt)
     dtnew(ilevel)=MIN(dtnew(ilevel),dt_rt/2.0**(ilevel-levelmin))
     if(static) RETURN
  endif
#endif

  if(pic) then

     dt_all=dtnew(ilevel); dt_loc=dt_all
     ekin_all=0.0; ekin_loc=0.0
     
     ! Compute maximum time step on active region
     if(numbl(myid,ilevel)>0)then
        ! Loop over grids
        ip=0
        igrid=headl(myid,ilevel)
        do jgrid=1,numbl(myid,ilevel)
           npart1=numbp(igrid)   ! Number of particles in the grid
           if(npart1>0)then
              ! Loop over particles
              ipart=headp(igrid)
              do jpart=1,npart1
                 ip=ip+1
                 ind_part(ip)=ipart
                 if(ip==nvector)then
                    call newdt2(ind_part,dt_loc,ekin_loc,ip,ilevel)
                    ip=0
                 end if
                 ipart=nextp(ipart)    ! Go to next particle
              end do
              ! End loop over particles
           end if
           igrid=next(igrid)   ! Go to next grid
        end do
        ! End loop over grids
        if(ip>0)call newdt2(ind_part,dt_loc,ekin_loc,ip,ilevel)
     end if
     
     ! Minimize time step over all cpus
#ifndef WITHOUTMPI
     call MPI_ALLREDUCE(dt_loc,dt_all,1,MPI_DOUBLE_PRECISION,MPI_MIN,&
          & MPI_COMM_WORLD,info)
     call MPI_ALLREDUCE(ekin_loc,ekin_all,1,MPI_DOUBLE_PRECISION,MPI_SUM,&
          & MPI_COMM_WORLD,info)
#endif
#ifdef WITHOUTMPI
     dt_all=dt_loc
     ekin_all=ekin_loc
#endif
     ekin_tot=ekin_tot+ekin_all
     dtnew(ilevel)=MIN(dtnew(ilevel),dt_all)
     
  end if

  if(hydro)call courant_fine(ilevel)
  
111 format('   Entering newdt_fine for level ',I2)

end subroutine newdt_fine
!#####################################################################
!#####################################################################
!#####################################################################
!#####################################################################
subroutine newdt2(ind_part,dt_loc,ekin_loc,nn,ilevel)
  use amr_commons
  use pm_commons
  use hydro_commons
  implicit none
  real(kind=8)::dt_loc,ekin_loc
  integer::nn,ilevel
  integer,dimension(1:nvector)::ind_part

  integer::i,idim,nx_loc
  real(dp)::dx,dx_loc,scale,dtpart
  real(dp),dimension(1:nvector),save::v2

  ! Compute time step
  dx=0.5D0**ilevel
  nx_loc=(icoarse_max-icoarse_min+1)
  scale=boxlen/dble(nx_loc)
  dx_loc=dx*scale

  v2(1:nn)=0.0D0
  do idim=1,ndim
     do i=1,nn
        v2(i)=MAX(v2(i),vp(ind_part(i),idim)**2)
!        v2(i)=v2(i)+vp(ind_part(i),idim)**2
     end do
  end do
  do i=1,nn
     if(v2(i)>0.0D0)then
        dtpart=courant_factor*dx_loc/sqrt(v2(i))
        dt_loc=MIN(dt_loc,dtpart)
     end if
  end do

  ! Compute kinetic energy
  do idim=1,ndim
     do i=1,nn
        ekin_loc=ekin_loc+0.5D0*mp(ind_part(i))*vp(ind_part(i),idim)**2
     end do
  end do
    
end subroutine newdt2




