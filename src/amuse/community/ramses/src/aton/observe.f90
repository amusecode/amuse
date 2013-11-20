module observe_commons
  use amr_parameters

  integer,parameter::observe_num_quantity=8
  real(dp),dimension(1:observe_num_quantity)::observe_quantity=0d0

  integer::observe_ilun_quantity=8000
  
  ! Star radiation statistics
  real(kind=8)::observe_total_star_source=0d0
  integer::observe_num_stars=0

end module observe_commons

subroutine observe_level(ilevel)
  use amr_commons
  use hydro_commons
  use cooling_module
  use observe_commons
  use radiation_commons, ONLY: Erad
  implicit none
#ifndef WITHOUTMPI
  include 'mpif.h'
#endif
  integer::ilevel
  integer::i,ivar,idim,ind,ncache,igrid,iskip
  integer::info,nleaf,ngrid,nx_loc
  integer,dimension(1:nvector),save::ind_grid,ind_cell,ind_leaf

  real(dp)::dx,cell_fraction
  real(kind=8),dimension(3)::comm_buffin,comm_buffout
  real(dp),dimension(1:nvector,1:nvar),save::uu

  real(dp),dimension(1:observe_num_quantity)::local_quantity,total_quantity
  real(dp)::density,xion,xion_mass,temperature,rad_intensity
  integer::density_bin

  real(dp)::scale_nH,scale_T2,scale_t,scale_v,scale_d,scale_l

  if(numbtot(1,ilevel)==0)return

  call units(scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2)

  ! Mesh spacing at that level
  ! TODO(tstranex): Is nx_loc actually needed and correct?
  nx_loc=icoarse_max-icoarse_min+1
  dx=0.5D0**ilevel/dble(nx_loc)
  cell_fraction=dx**ndim

  local_quantity = 0d0
  total_quantity = 0d0

  ! Loop over active grids by vector sweeps
  ncache=active(ilevel)%ngrid
  do igrid=1,ncache,nvector
     ngrid=MIN(nvector,ncache-igrid+1)
     do i=1,ngrid
        ind_grid(i)=active(ilevel)%igrid(igrid+i-1)
     end do
     
     ! Loop over cells
     do ind=1,twotondim        
        iskip=ncoarse+(ind-1)*ngridmax
        do i=1,ngrid
           ind_cell(i)=ind_grid(i)+iskip
        end do
        
        ! Gather leaf cells
        nleaf=0
        do i=1,ngrid
           if(son(ind_cell(i))==0)then
              nleaf=nleaf+1
              ind_leaf(nleaf)=ind_cell(i)
           end if
        end do

        ! Gather hydro variables
        do ivar=1,nvar
           do i=1,nleaf
              uu(i,ivar)=uold(ind_leaf(i),ivar)
           end do
        end do

        do i=1,nleaf
           density = uu(i,1)*scale_d*0.76/mH
           xion = uu(i,ixion)/uu(i,1)
           xion_mass = xion * density
           temperature = uu(i,ndim+2)
           do ivar=1,ndim
              temperature = temperature - 0.5*uu(i,ivar+1)**2/uu(i,1)
           end do
           temperature = temperature*(gamma-1.0)/uu(i,1)*scale_T2/(1 + xion)
           if (aton) then
              rad_intensity = Erad(ind_leaf(i))
           else
              rad_intensity = -1
           end if

           local_quantity(1) = local_quantity(1) + 1d0
           local_quantity(2) = local_quantity(2) + density
           local_quantity(3) = local_quantity(3) + xion
           local_quantity(4) = local_quantity(4) + xion_mass
           local_quantity(5) = local_quantity(5) + temperature
           local_quantity(6) = local_quantity(6) + rad_intensity
        end do
     end do
     ! End loop over cells
  end do
  ! End loop over grids

  local_quantity = local_quantity*cell_fraction

  local_quantity(7) = local_quantity(7) + observe_total_star_source
  local_quantity(8) = local_quantity(8) + observe_num_stars

  ! Compute global quantities
  call MPI_ALLREDUCE(local_quantity,total_quantity,observe_num_quantity, &
       & MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)

  observe_quantity = total_quantity
end subroutine observe_level

subroutine observe_init()
  use hydro_commons
  use amr_commons
  use observe_commons
  implicit none
  character(len=20)::mode
  if (myid.ne.1) return

  if (nrestart.eq.0) then
     mode = 'ASIS'
  else
     mode = 'APPEND'
  end if

  open(unit=observe_ilun_quantity,position=mode, &
       & file='averages.txt',form='formatted')
  write(observe_ilun_quantity,*) '# aexp <one> <nH> <x> <x*nH> <T> <J> <S> <num_stars>'  
end subroutine observe_init

subroutine observe_stop()
  use hydro_commons
  use amr_commons
  use observe_commons
  implicit none
  if (myid.ne.1) return
  close(unit=observe_ilun_quantity)
end subroutine observe_stop

subroutine observe_output()
  use hydro_commons
  use amr_commons
  use observe_commons
  implicit none
  integer::bin
  real(dp)::bin_value
  if (myid.ne.1) return

  write(*,*)"Output observed quantities"

  write(observe_ilun_quantity,*) aexp, &
       & observe_quantity(1), &
       & observe_quantity(2), &
       & observe_quantity(3), &
       & observe_quantity(4), &
       & observe_quantity(5), &
       & observe_quantity(6), &
       & observe_quantity(7), &
       & observe_quantity(8)

  flush(unit=observe_ilun_quantity)
end subroutine observe_output
