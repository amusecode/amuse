!################################################################
!################################################################
!################################################################
!################################################################
subroutine gas_ana
  use amr_commons
  use hydro_commons
  implicit none
#ifndef WITHOUTMPI
  include 'mpif.h'
#endif
  !----------------------------------------------------------------------
  ! Description: This module computes the PDF (probability distribution
  ! function) of log(rho) on the fly for RAMSES hydro runs. The routine
  ! loops over all leaf cells. Every cell is weighted by its volume.
  ! Andreas Bleuler 26.10.2011
  !----------------------------------------------------------------------
  ! local constants
  real(dp)::scale_nH,scale_T2,scale_l,scale_d,scale_t,scale_v
  real(kind=8),dimension(1:twotondim,1:3)::xc
  ! other variables
  integer ::ncache,ngrid
  integer ::igrid,ix,iy,iz,ind,i,iskip,nx_loc
  integer ::info
  real(kind=8),dimension(1:3)::skip_loc
  real(kind=8)::d,xx,yy,zz,vx,vy,vz,dx,dx_loc,scale,vol_loc,dx_min,vol_min,v_sound
  real(kind=8)::mini_dens,mini_dens_tot,maxi_dens,maxi_dens_tot,l_max,l_min,width
  integer ,dimension(1:nvector),save::ind_grid,ind_cell
  integer::ilevel,hist_ind
  character(LEN=5)::nchar
  real(kind=8),allocatable,dimension(:)::hist,hist_tot
  real(kind=8)::m,m_tot,v_rms,v_rms_tot

  if(.not. hydro)return


  ! Conversion factor from user units to cgs units
  call units(scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2)
  

  maxi_dens=0.d0
  mini_dens=huge(0.d0)

  xx=5.d-1; yy=5.d-1; zz=5.d-1
  vx=0.d0; vy=0.d0; vz=0.d0

  m=0.d0; m_tot=0.d0
  v_rms=0.d0; v_rms_tot=0.d0
        
  do ilevel=levelmin,nlevelmax
     if (numbtot(1,ilevel)/=0)then 
        
        ! Mesh spacing in that level
        dx=0.5D0**ilevel 
        nx_loc=(icoarse_max-icoarse_min+1)
        skip_loc=(/0.0d0,0.0d0,0.0d0/)
        if(ndim>0)skip_loc(1)=dble(icoarse_min)
        if(ndim>1)skip_loc(2)=dble(jcoarse_min)
        if(ndim>2)skip_loc(3)=dble(kcoarse_min)
        scale=boxlen/dble(nx_loc)
        dx_loc=dx*scale
        vol_loc=dx_loc**ndim
        dx_min=(0.5D0**nlevelmax)*scale
        vol_min=dx_min**ndim
        
        
        ! Cells center position relative to grid center position
        do ind=1,twotondim  
           iz=(ind-1)/4
           iy=(ind-1-4*iz)/2
           ix=(ind-1-2*iy-4*iz)
           xc(ind,1)=(dble(ix)-0.5D0)*dx
           xc(ind,2)=(dble(iy)-0.5D0)*dx
           xc(ind,3)=(dble(iz)-0.5D0)*dx
        end do
        
        !------------------------------------------------
        ! find min and max density
        !------------------------------------------------
        ! Loop over grids
        ncache=active(ilevel)%ngrid
        do igrid=1,ncache,nvector
           ngrid=MIN(nvector,ncache-igrid+1)
           do i=1,ngrid
              ind_grid(i)=active(ilevel)%igrid(igrid+i-1)
           end do
           do ind=1,twotondim
              iskip=ncoarse+(ind-1)*ngridmax
              do i=1,ngrid
                 ind_cell(i)=iskip+ind_grid(i)
              end do
              
              do i=1,ngrid

                 ! Get cell center positions 
                 xx=(xg(ind_grid(i),1)+xc(ind,1)-skip_loc(1))*scale
#if NDIM>1
                 yy=(xg(ind_grid(i),2)+xc(ind,2)-skip_loc(2))*scale
#endif
#if NDIM>2
                 zz=(xg(ind_grid(i),3)+xc(ind,3)-skip_loc(3))*scale
#endif          
                 if (son(ind_cell(i))==0)then
                    if(ana_xmi<xx .and. xx<ana_xma .and. ana_ymi<yy .and. yy<ana_yma .and. ana_zmi<zz .and. zz<ana_zma)then
                       if (uold(ind_cell(i),1) > maxi_dens)then
                          maxi_dens=uold(ind_cell(i),1)
                       endif
                       if (uold(ind_cell(i),1) < mini_dens)then
                          mini_dens=uold(ind_cell(i),1)
                       endif
                    endif
                 endif
              end do
           end do
        end do
     end if
  end do
  
  !---------------------------------
  ! communicate
  !---------------------------------

#ifndef WITHOUTMPI
  call MPI_ALLREDUCE(mini_dens,mini_dens_tot,1,MPI_DOUBLE_PRECISION,MPI_MIN,MPI_COMM_WORLD,info)
#endif
#ifdef WITHOUTMPI
  mini_dens_tot=mini_dens
#endif
#ifndef WITHOUTMPI
  call MPI_ALLREDUCE(maxi_dens,maxi_dens_tot,1,MPI_DOUBLE_PRECISION,MPI_MAX,MPI_COMM_WORLD,info)
#endif
#ifdef WITHOUTMPI
  maxi_dens_tot=maxi_dens
#endif


  l_max=log10(maxi_dens_tot*1.001)
  l_min=log10(mini_dens_tot*0.999)
  width=l_max-l_min

  allocate(hist(1:nbins))
  allocate(hist_tot(1:nbins))
  hist=0.d0
  hist_tot=0.d0


  do ilevel=levelmin,nlevelmax
     if(numbtot(1,ilevel)/=0)then

        ! Mesh spacing in that level
        dx=0.5D0**ilevel 
        nx_loc=(icoarse_max-icoarse_min+1)
        skip_loc=(/0.0d0,0.0d0,0.0d0/)
        if(ndim>0)skip_loc(1)=dble(icoarse_min)
        if(ndim>1)skip_loc(2)=dble(jcoarse_min)
        if(ndim>2)skip_loc(3)=dble(kcoarse_min)
        scale=boxlen/dble(nx_loc)
        dx_loc=dx*scale
        vol_loc=dx_loc**ndim
        dx_min=(0.5D0**nlevelmax)*scale
        vol_min=dx_min**ndim


        ! Cells center position relative to grid center position
        do ind=1,twotondim  
           iz=(ind-1)/4
           iy=(ind-1-4*iz)/2
           ix=(ind-1-2*iy-4*iz)
           xc(ind,1)=(dble(ix)-0.5D0)*dx
           xc(ind,2)=(dble(iy)-0.5D0)*dx
           xc(ind,3)=(dble(iz)-0.5D0)*dx
        end do

        ! Loop over grids
        ncache=active(ilevel)%ngrid
        do igrid=1,ncache,nvector
           ngrid=MIN(nvector,ncache-igrid+1)
           do i=1,ngrid
              ind_grid(i)=active(ilevel)%igrid(igrid+i-1)
           end do
           do ind=1,twotondim
              iskip=ncoarse+(ind-1)*ngridmax
              do i=1,ngrid
                 ind_cell(i)=iskip+ind_grid(i)
              end do

              do i=1,ngrid

                 ! Get cell center positions 
                 xx=(xg(ind_grid(i),1)+xc(ind,1)-skip_loc(1))*scale
#if NDIM>1
                 yy=(xg(ind_grid(i),2)+xc(ind,2)-skip_loc(2))*scale
#endif
#if NDIM>2
                 zz=(xg(ind_grid(i),3)+xc(ind,3)-skip_loc(3))*scale
#endif          
                 if (son(ind_cell(i))==0)then
                    if(ana_xmi<xx .and. xx<ana_xma .and. ana_ymi<yy .and. yy<ana_yma .and. ana_zmi<zz .and. zz<ana_zma)then
                       d=dble(uold(ind_cell(i),1))
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                       ! DETERMINE THE LOCAL SOUND SPEED HERE
                       ! THIS DEPENDS ON YOUR THERMODINAMICAL ASSUMPTIONS
                       ! (SEE COOLING_FINE)
                       ! MAYBE YOU WANT TO USE A PATCHED VERSION OF THIS ROUTINE
                       
                       ! assuming a polytropic equation of state
                       ! v_sound=(T2_star/scale_T2*gamma*(d*scale_nH/n_star)**(gamma-1))**0.5

                       ! isothermal
                       v_sound=dble((T2_star/scale_T2))**0.5d0
                       
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

                       vx=dble(uold(ind_cell(i),2))/(v_sound)
#if NDIM>1
                       vy=dble(uold(ind_cell(i),3))/(v_sound)
#endif
#if NDIM>2
                       vz=dble(uold(ind_cell(i),4))/(v_sound)
#endif
                       v_rms=v_rms+(vx**2.d0+vy**2.d0+vz**2.d0)*vol_loc/d 
                       m=m+vol_loc*d
                       ! log(rho) PDF
                       hist_ind=1+int((log10(d)-l_min)/width*nbins)
                       hist(hist_ind)=hist(hist_ind)+vol_loc
                    endif
                 endif
              end do
           end do
        end do
     end if         
  end do



  !---------------------------------
  ! communicate results
  !---------------------------------

#ifndef WITHOUTMPI
  call MPI_ALLREDUCE(hist,hist_tot,nbins,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
#endif
#ifdef WITHOUTMPI
  hist_tot=hist
#endif
#ifndef WITHOUTMPI
  call MPI_ALLREDUCE(m,m_tot,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
#endif
#ifdef WITHOUTMPI
  m_tot=m
#endif
#ifndef WITHOUTMPI
  call MPI_ALLREDUCE(v_rms,v_rms_tot,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
#endif
#ifdef WITHOUTMPI
  v_rms_tot=v_rms
#endif
v_rms_tot=(v_rms_tot/m_tot)**5.d-1
  !---------------------------------
  ! write textfile
  !---------------------------------
  if (myid==1)then
     write(*,*)'rms mach number in the region of interest= ',v_rms_tot
     call title(ifout-1,nchar)
     open(unit=20,file=TRIM('output_'//TRIM(nchar)//'/density_pdf.txt'),form='formatted')
     write(20,*)nbins
     do i=1,nbins
        write(20,*)(i-0.5)*width/nbins+l_min+log10(scale_d),hist_tot(i)
     end do
  end if

  deallocate(hist,hist_tot)

end subroutine gas_ana


subroutine read_gas_analytics_params()
  use amr_commons
  use hydro_commons

  implicit none
#ifndef WITHOUTMPI
  include 'mpif.h'
#endif

  !-------------------------------------------------- 
  ! Namelist definitions                              
  !-------------------------------------------------- 
  namelist/gas_analytics_params/nbins,ana_xmi,ana_xma,ana_ymi,ana_yma,ana_zmi,ana_zma
  ! Read namelist file                                                                         
  rewind(1)
  read(1,NML=gas_analytics_params,END=101)
  goto 102
101 if(myid==1)write(*,*)' You did not setup &GAS_ANALYTICS_PARAMS in parameter file. Defaults will be used...'
  ana_xmi=0.d0; ana_ymi=0.d0; ana_zmi=0.d0
  ana_xma=1.d0; ana_yma=1.d0; ana_zma=1.d0
  nbins=1000
102 rewind(1)
end subroutine read_gas_analytics_params
