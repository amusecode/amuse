!#########################################################
!#########################################################
!#########################################################
!#########################################################
subroutine gravana(x,f,dx,ncell)
  use amr_commons
  use poisson_parameters  
  implicit none
  integer ::ncell                         ! Size of input arrays
  real(dp)::dx                            ! Cell size
  real(dp),dimension(1:nvector,1:ndim)::f ! Gravitational acceleration
  real(dp),dimension(1:nvector,1:ndim)::x ! Cell center position.
  !================================================================
  ! This routine computes the acceleration using analytical models.
  ! x(i,1:ndim) are cell center position in [0,boxlen] (user units).
  ! f(i,1:ndim) is the gravitational acceleration in user units.
  !================================================================
  integer::idim,i
  real(dp)::gmass,emass,xmass,ymass,zmass,rr,rx,ry,rz
  real(dp)::rmin,rpla,qpla,epla,msun,xsun,ysun,zsun
  real(dp)::esmooth,aa,qq,phip,xa,ya,za,xb,yb,zb
  real(dp)::rxa,rya,rza,rra,rxb,ryb,rzb,rrb,ff1,ff2,ffr
  real(dp)::omega0,omega1,omega2

  ! Constant vector
  if(gravity_type==1)then 
     do idim=1,ndim
        do i=1,ncell
           f(i,idim)=gravity_params(idim)
        end do
     end do
  end if

  ! Point mass
  if(gravity_type==2)then 
     gmass=gravity_params(1) ! GM
     emass=gravity_params(2) ! Softening length
     xmass=gravity_params(3) ! Point mass coordinates
     ymass=gravity_params(4)
     zmass=gravity_params(5)
     do i=1,ncell
        rx=0.0d0; ry=0.0d0; rz=0.0d0
        rx=x(i,1)-xmass
#if NDIM>1
        ry=x(i,2)-ymass
#endif
#if NDIM>2
        rz=x(i,3)-zmass
#endif
        rr=sqrt(rx**2+ry**2+rz**2+emass**2)
        f(i,1)=-gmass*rx/rr**3
#if NDIM>1
        f(i,2)=-gmass*ry/rr**3
#endif
#if NDIM>2
        f(i,3)=-gmass*rz/rr**3
#endif
     end do
  end if

  ! Planet orbiting around the Sun
  if(gravity_type==3)then 
     msun=1.0                 ! GMsun
     rpla=1.0                 ! Planet radius
     rmin=gravity_params(1)   ! Hard sphere radius
     qpla=gravity_params(4)   ! Planet to Sun mass ratio
     epla=gravity_params(5)   ! Planet gravity softening
     omega1=1.0/rpla**(1.5)
     omega2=1.0/rmin**(1.5)
     omega0=omega1*omega2/(omega1+omega2)
     phip=omega0*t
     xsun=0.5*boxlen
     ysun=0.5*boxlen
     zsun=0.5*boxlen
     xa=rpla*(1.0-qpla)*cos(phip)+xsun
     ya=rpla*(1.0-qpla)*sin(phip)+ysun
     za=0.0                      +zsun
     xb=-rpla*qpla*cos(phip)     +xsun
     yb=-rpla*qpla*sin(phip)     +ysun
     zb= 0.0                     +zsun
     do i=1,ncell
        rxa=0.0d0; rya=0.0d0; rza=0.0d0 ! Planet
        rxb=0.0d0; ryb=0.0d0; rzb=0.0d0 ! Sun
        if(ndim>0)rxa=x(i,1)-xa
        if(ndim>1)rya=x(i,2)-ya
        if(ndim>2)rza=x(i,3)-za
        if(ndim>0)rxb=x(i,1)-xb
        if(ndim>1)ryb=x(i,2)-yb
        if(ndim>2)rzb=x(i,3)-zb
        rra=sqrt(rxa**2+rya**2+rza**2+epla**2)
        rrb=sqrt(rxb**2+ryb**2+rzb**2)
        ff1=sqrt(msun*(1.0-qpla)/rrb**3)
        ff2=sqrt(msun*(1.0-qpla)/rmin**3)
        ffr=ff1*ff2/(ff1+ff2)
        ffr=ffr**2
        if(ndim>0)f(i,1)=-rxb*ffr
        if(ndim>1)f(i,2)=-ryb*ffr
        if(ndim>2)f(i,3)=-rzb*ffr
        if(ndim>0)f(i,1)=f(i,1)-msun*qpla*rxa/rra**3
        if(ndim>1)f(i,2)=f(i,2)-msun*qpla*rya/rra**3
        if(ndim>2)f(i,3)=f(i,3)-msun*qpla*rza/rra**3
     end do
  end if

end subroutine gravana
!#########################################################
!#########################################################
!#########################################################
!#########################################################
subroutine barotrop(ilevel)
  use amr_commons
  use hydro_commons
!  use const
  implicit none
  integer::ilevel
  
  integer::igrid,ngrid,ncache,i,ind,iskip,ix,iy,iz
  integer::info,ibound,nx_loc,idim
  real(dp)::dx,dx_loc,scale,d,u,v,w,e
 ! real(kind=8)::rho_max_loc,rho_max_all,epot_loc,epot_all
 ! real(dp),dimension(1:8,1:3)::xc
  real(dp),dimension(1:3)::skip_loc

  integer ,dimension(1:nvector),save::ind_grid,ind_cell
!  real(dp),dimension(1:nvector,1:ndim),save::xx
!  real(dp),dimension(1:nvector,1:nvar),save::uu
 
  if(verbose)write(*,111)ilevel


!   write(*,*) 'barotrop'
!   write(*,*) 'temp, rho_ad ',temp, rho_ad
  rho_ad = 1.d10 * 2.2 * 1.64d-24         !in g cm^-3
  rho_ad = rho_ad / 2.d33 * (3.08d18)**3  !normalise in code units (solar mass, pc)
  rho_ad = f_rho_ad * rho_ad              !f_rho_ad is readed in the .nml

  ! Loop over all AMR linked lists (including physical boundaries).
  do ibound=1,nboundary+ncpu

     if(ibound==myid)then
        ncache=active(ilevel)%ngrid
     else if(ibound<=ncpu)then
        ncache=reception(ibound,ilevel)%ngrid
     else
        ncache=boundary(ibound-ncpu,ilevel)%ngrid
     end if

     ! Loop over grids by vector sweeps
     do igrid=1,ncache,nvector
        ngrid=MIN(nvector,ncache-igrid+1)
        if(ibound==myid)then
           do i=1,ngrid
              ind_grid(i)=active(ilevel)%igrid(igrid+i-1)
           end do
        else if(ibound<=ncpu)then
           do i=1,ngrid
              ind_grid(i)=reception(ibound,ilevel)%igrid(igrid+i-1)
           end do
        else
           do i=1,ngrid
              ind_grid(i)=boundary(ibound-ncpu,ilevel)%igrid(igrid+i-1)
           end do
        end if
     
        ! Loop over cells
        do ind=1,twotondim

           ! Gather cell indices
           iskip=ncoarse+(ind-1)*ngridmax
           do i=1,ngrid
              ind_cell(i)=iskip+ind_grid(i)
           end do

           ! Convert to specific internal energy

!              write(*,*) 'ngrid ',ngrid

           do i=1,ngrid
              d=uold(ind_cell(i),1)
 

              !! test the density if it is high enough then force
              !! barotropic equation of state
              if(d .gt. rho_ad / 10000.) then


              u=0.0; v=0.0; w=0.0
              if(ndim>0)u=uold(ind_cell(i),2)/d
              if(ndim>1)v=uold(ind_cell(i),3)/d
              if(ndim>2)w=uold(ind_cell(i),4)/d


              e = temp*( 1.D0 + (d/rho_ad)**0.6666666 ) / (gamma-1.D0) 

              uold(ind_cell(i),ndim+2) =  d*(e + 0.5*(u**2+v**2+w**2))


              endif

           end do


        end do
        ! End loop over cells

     end do
     ! End loop over grids

  end do
  ! End loop over linked lists

111 format('   Entering impose_iso for level ',I2)

end subroutine barotrop


