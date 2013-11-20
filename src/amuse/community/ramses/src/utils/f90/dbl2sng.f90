program dbl2sng
  !---------------------------------------------------------------------
  ! Ce programme lit le header d'un fichier amr.out d'une simulation RAMSES
  ! Version F90 par R. Teyssier le 28/02/00
  ! Version parallele par R. Teyssier le 13/06/03
  ! f90 head_amr.f90 -o ~/bin/header
  !---------------------------------------------------------------------
  implicit none
  integer,parameter::dp=kind(1.0D0) ! real*8
  integer,parameter::sp=kind(1.0E0) ! real*4
  integer,parameter::ndim=3
  integer,parameter :: IRandNumSize = 4, IBinarySize = 48
  integer,parameter::twotondim=2**ndim
  integer,parameter::threetondim=3**ndim
  integer,parameter::twondim=2*ndim
  integer::ibound,ilevel,i,ind,istart,igrid2
  integer::ncpu,ndim1,npart,n,icpu,ind1,idim,icell2,igrid0,iskip2,ncache
  integer::nx,ny,nz,nstep,ncell2,iskip1,igrid1,icell1
  integer::nlevelmax,noutput,iout,ifout,iback,ncoarse,ncell
  integer::ngridmax,ngridmax2,nstep_coarse,nstep_coarse_old
  real(dp)::const,mass_tot_0,rho_tot,h0,aexp_ini,epot_tot_int,epot_tot_old
  real(dp)::t,aexp,hexp,boxlen,aexp_old,mass_sph
  real(dp)::omega_m,omega_l,omega_k,omega_b
  real(dp)::scale_l,scale_d,scale_t
  logical::ok
  character*5::nchar,ncharcpu
  character*50::nomfich,nomdir
  integer::iargc,ipos,ilun
  character(len=128)::arg
  integer,dimension(IRandNumSize) :: localseed=-1
  integer::nstar_tot=0              ! Total number of star particle
  real(dp)::mstar_tot=0             ! Total star mass
  real(dp)::mstar_lost=0            ! Missing star mass

  ! Mesh parameters
  integer::geom=1             ! 1: cartesian, 2: cylindrical, 3: spherical
  integer::levelmin=1         ! Full refinement up to levelmin
  character(len=128)::ordering='hilbert'
  ! Cosmology
  real(dp)::boxlen_ini        ! Box size in h-1 Mpc
  real(dp)::xdum

  ! Output times
  integer,parameter::MAXOUT=1000
  real(dp),dimension(1:MAXOUT)::aout=1.1       ! Output expansion factors
  real(dp),dimension(1:MAXOUT)::tout=0.0       ! Output times

  integer,allocatable,dimension(:,:)::headl,headl2
  integer,allocatable,dimension(:,:)::taill,taill2
  integer,allocatable,dimension(:,:)::numbl
  integer,allocatable,dimension(:,:)::numbtot
  integer::headf,tailf,numbf,used_mem,used_mem_tot
  integer::headf2,tailf2,numbf2,used_mem2,used_mem_tot2
  real(dp),allocatable,dimension(:,:)::xg      ! grids position
  integer ,allocatable,dimension(:,:)::nbor,nbor2    ! neighboring father cells
  integer ,allocatable,dimension(:)  ::father,father2  ! father cell
  integer ,allocatable,dimension(:)  ::next,next2    ! next grid in list
  integer ,allocatable,dimension(:)  ::prev,prev2    ! previous grid in list
  integer ,allocatable,dimension(:)  ::son,son2     ! sons grids
  integer ,allocatable,dimension(:)  ::flag1,iig   ! flag for refine
  integer ,allocatable,dimension(:)  ::mapping
  real(dp),allocatable,dimension(:)  ::xdp
  real(sp),allocatable,dimension(:)  ::xsp
  real(dp),allocatable,dimension(:,:)::xp
  real(dp),allocatable,dimension(:,:)::vp
  real(dp),allocatable,dimension(:)::mp
  integer ,allocatable,dimension(:)  ::idp,lvp


  ! Global indexing
  integer ,allocatable,dimension(:)  ::cpu_map  ! domain decomposition

  ! Hilbert key
  real(kind=8),allocatable,dimension(:)::hilbert_key
  real(kind=8),allocatable,dimension(:)::bound_key,bound_key2
  real(kind=8)                         ::order_all_min,order_all_max

  n = iargc()
  if (n.NE.1) then
     print *, 'usage: header backup_dir'
     stop
  end if
  call getarg(1,arg)
  nomdir=trim(arg)//'/'
  !-----------------------------------------------
  ! Lecture des fichiers outputs au format RAMSES
  !-----------------------------------------------
  ipos=INDEX(nomdir,'backup_')
  nchar=nomdir(ipos+7:ipos+13)
  ! Lecture du fichier AMR
  nomfich=TRIM(nomdir)//'amr_'//TRIM(nchar)//'.bak00001'
  inquire(file=nomfich, exist=ok)
  if (.not. ok) then
     write(*,*)'File '//TRIM(nomfich)//' not found'
     stop
  else
     ilun=10
     open(unit=ilun,file=nomfich,status='old',form='unformatted')
     read(ilun)ncpu
     read(ilun)ndim1
     read(ilun)nx,ny,nz
     read(ilun)nlevelmax
     close(ilun)
  end if

  call system('mkdir defrag')

  allocate(headl(1:ncpu,1:nlevelmax))
  allocate(taill(1:ncpu,1:nlevelmax))
  allocate(numbl(1:ncpu,1:nlevelmax))
  allocate(numbtot(1:10,1:nlevelmax))
  allocate(bound_key(0:ncpu))

  do icpu=1,ncpu
     write(*,*)'ncpu=',icpu
     call title(icpu,ncharcpu)
     nomfich=TRIM(nomdir)//'amr_'//TRIM(nchar)//'.bak'//TRIM(ncharcpu)
     inquire(file=nomfich, exist=ok)
     if (.not. ok) then
        write(*,*)'File '//TRIM(nomfich)//' not found'
        stop
     else
        ilun=10
        open(unit=ilun,file=nomfich,status='old',form='unformatted')
        ! Write grid variables
        read(ilun)ncpu
        read(ilun)ndim1
        read(ilun)nx,ny,nz
        read(ilun)nlevelmax
        read(ilun)ngridmax
        read(ilun)boxlen
        ! Read time variables
        read(ilun)noutput,iout,ifout,iback

        ncoarse=nx*ny*nz
        ncell=ncoarse+twotondim*ngridmax

        read(ilun)tout(1:noutput)
        read(ilun)aout(1:noutput)
        read(ilun)t
        read(ilun)nstep,nstep_coarse,nstep_coarse_old
        read(ilun)const,mass_tot_0,rho_tot
        read(ilun)omega_m,omega_l,omega_k,omega_b,h0,aexp_ini,boxlen_ini
        read(ilun)aexp,hexp,aexp_old,epot_tot_int,epot_tot_old
        read(ilun)mass_sph
        ! Read levels variables        
        read(ilun)headl(1:ncpu,1:nlevelmax)
        read(ilun)taill(1:ncpu,1:nlevelmax)
        read(ilun)numbl(1:ncpu,1:nlevelmax)
        read(ilun)numbtot(1:10,1:nlevelmax)

        ngridmax2=0
        do ilevel=1,nlevelmax
           do ibound=1,ncpu
              ngridmax2=ngridmax2+numbl(ibound,ilevel)
           end do
        end do
        ngridmax2=ngridmax2+1
        ncell2=ncoarse+twotondim*ngridmax2

        write(*,*)'ngridmax =',ngridmax
        write(*,*)'ngridmax2=',ngridmax2

        allocate(prev(1:ngridmax))
        allocate(next(1:ngridmax))
        allocate(prev2(1:ngridmax2))
        allocate(next2(1:ngridmax2))
        allocate(mapping(1:ngridmax))
        mapping=0
        allocate(headl2(1:ncpu,1:nlevelmax))
        allocate(taill2(1:ncpu,1:nlevelmax))
        allocate(flag1(0:ncell2))
        allocate(cpu_map(1:ncell2))
        allocate(xg(1:ngridmax2,1:ndim))
        allocate(father(1:ngridmax2))
        allocate(nbor(1:ngridmax2,1:twondim))
        allocate(son(1:ncell2))
        allocate(father2(1:ngridmax2))
        allocate(nbor2(1:ngridmax2,1:twondim))
        allocate(son2(1:ncell2))

        ! Read linked list
        read(ilun)prev
        read(ilun)next
        ! Read free memory
        read(ilun)headf,tailf,numbf,used_mem,used_mem_tot
        ! Read cpu boundaries
        read(ilun)bound_key(0:ncpu)
        ! Read coarse level
        read(ilun)son(1:ncoarse)
        read(ilun)flag1(1:ncoarse)
        read(ilun)cpu_map(1:ncoarse)

        ! Read fine levels
        igrid2=0
        do ilevel=1,nlevelmax
           do ibound=1,ncpu

              ncache=numbl(ibound,ilevel)

              if(ncache>0)then
                 allocate(xdp(1:ncache),iig(1:ncache))
                 xdp=0d0
                 iig=0

                 ! Read grid index
                 read(ilun)iig
                 do i=1,ncache
                    mapping(iig(i))=igrid2+i
                 end do

                 ! Read grid center
                 do idim=1,ndim
                    read(ilun)xdp
                    do i=1,ncache
                      xg(igrid2+i,idim)=xdp(i)
                    end do
                 end do

                 ! Read tree structure
                 read(ilun)iig
                 do i=1,ncache
                    father(igrid2+i)=iig(i)
                 end do
                 do ind=1,twondim
                    read(ilun)iig
                    do i=1,ncache
                       nbor(igrid2+i,ind)=iig(i)
                    end do
                 end do
                 do ind=1,twotondim
                    iskip2=ncoarse+(ind-1)*ngridmax2
                    read(ilun)iig
                    do i=1,ncache
                       son(igrid2+i+iskip2)=iig(i)
                    end do
                 end do
                 do ind=1,twotondim
                    iskip2=ncoarse+(ind-1)*ngridmax2
                    read(ilun)iig
                    do i=1,ncache
                       cpu_map(igrid2+i+iskip2)=iig(i)
                    end do
                 end do
                 igrid2=igrid2+ncache
                 deallocate(xdp,iig)
              end if

           end do
        end do
        close(ilun)
     endif
     mapping(headf)=ngridmax2
     
     igrid2=0
     do ilevel=1,nlevelmax
        do ibound=1,ncpu
           ncache=numbl(ibound,ilevel)
           if(ncache>0)then
              headl2(ibound,ilevel)=igrid2+1
              taill2(ibound,ilevel)=igrid2+ncache
              prev2(igrid2+1)=0
              do i=2,ncache
                 prev2(igrid2+i)=igrid2+i-1
              end do
              do i=1,ncache-1
                 next2(igrid2+i)=igrid2+i+1
              end do
              next2(igrid2+ncache)=0
              igrid2=igrid2+ncache
           else
              headl2(ibound,ilevel)=0
              taill2(ibound,ilevel)=0
           endif
        end do
     end do
     headf2=ngridmax2
     tailf2=ngridmax2
     numbf2=1
     used_mem2=ngridmax2-1
     used_mem_tot2=used_mem2
     prev2(headf2)=0
     next2(tailf2)=0

     do i=1,ncoarse
        igrid1=son(i)
        if(igrid1>0)then
           igrid2=mapping(igrid1) 
        else
           igrid2=0
        endif
        son2(i)=igrid2 
     end do

     igrid0=0
     do ilevel=1,nlevelmax
        do ibound=1,ncpu
           
           ncache=numbl(ibound,ilevel)
           
           if(ncache>0)then

              do i=1,ncache
                 icell1=father(igrid0+i)
                 if(icell1>ncoarse)then
                    ind1=(icell1-ncoarse-1)/ngridmax+1
                    iskip1=ncoarse+(ind1-1)*ngridmax
                    igrid1=(icell1-iskip1)
                    igrid2=mapping(igrid1)
                    iskip2=ncoarse+(ind1-1)*ngridmax2
                    icell2=iskip2+igrid2
                 else
                    icell2=icell1
                 endif
                 father2(igrid0+i)=icell2
              end do

              do ind=1,twondim
                 do i=1,ncache
                    icell1=nbor(igrid0+i,ind)
                    if(icell1>ncoarse)then
                       ind1=(icell1-ncoarse-1)/ngridmax+1
                       iskip1=ncoarse+(ind1-1)*ngridmax
                       igrid1=(icell1-iskip1)
                       igrid2=mapping(igrid1)
                       iskip2=ncoarse+(ind1-1)*ngridmax2
                       icell2=iskip2+igrid2
                    else
                       icell2=icell1
                    endif
                    nbor2(igrid0+i,ind)=icell2
                 end do
              end do

              do ind=1,twotondim
                 iskip2=ncoarse+(ind-1)*ngridmax2
                 do i=1,ncache
                    igrid1=son(igrid0+i+iskip2)
                    if(igrid1>0)then
                       igrid2=mapping(igrid1)
                    else
                       igrid2=0
                    endif
                    son2(igrid0+i+iskip2)=igrid2
                 end do
              end do

              igrid0=igrid0+ncache
           end if
              
        end do
     end do

     ilun=11
     nomfich='defrag/amr_'//TRIM(nchar)//'.bak'//TRIM(ncharcpu)
     open(unit=ilun,file=nomfich,status='new',form='unformatted')
     ! Write grid variables
     write(ilun)ncpu
     write(ilun)ndim1
     write(ilun)nx,ny,nz
     write(ilun)nlevelmax
     write(ilun)ngridmax2
     write(ilun)real(boxlen)
     ! Write time variables
     write(ilun)noutput,iout,ifout,iback
     write(ilun)real(tout(1:noutput))
     write(ilun)real(aout(1:noutput))
     write(ilun)real(t)
     write(ilun)nstep,nstep_coarse,nstep_coarse_old
     write(ilun)real(const),real(mass_tot_0),real(rho_tot)
     write(ilun)real(omega_m),real(omega_l),real(omega_k),real(omega_b),real(h0),real(aexp_ini),real(boxlen_ini)
     write(ilun)real(aexp),real(hexp),real(aexp_old),real(epot_tot_int),real(epot_tot_old)
     write(ilun)real(mass_sph)
     ! Write levels variables
     write(ilun)headl2(1:ncpu,1:nlevelmax)
     write(ilun)taill2(1:ncpu,1:nlevelmax)
     write(ilun)numbl(1:ncpu,1:nlevelmax)
     write(ilun)numbtot(1:10,1:nlevelmax)
     ! Write linked list
     write(ilun)prev2
     write(ilun)next2
     ! Write free memory
     write(ilun)headf2,tailf2,numbf2,used_mem2,used_mem_tot2
     ! Write cpu boundaries
     write(ilun)bound_key(0:ncpu)
     ! Write coarse level
     write(ilun)son2(1:ncoarse)
     write(ilun)flag1(1:ncoarse)
     write(ilun)cpu_map(1:ncoarse)

     ! Write fine levels
     igrid2=0
     do ilevel=1,nlevelmax
        do ibound=1,ncpu
           ncache=numbl(ibound,ilevel)
           if(ncache>0)then
              allocate(xsp(1:ncache),iig(1:ncache))
              ! Write grid index
              do i=1,ncache
                 iig(i)=igrid2+i
              end do
              write(ilun)iig
              ! Write grid center
              do idim=1,ndim
                 do i=1,ncache
                    xsp(i)=xg(igrid2+i,idim)
                 end do
                 write(ilun)xsp
              end do
              ! Write tree structure
              do i=1,ncache
                 iig(i)=father2(igrid2+i)
              end do
              write(ilun)iig
              do ind=1,twondim
                 do i=1,ncache
                    iig(i)=nbor2(igrid2+i,ind)
                 end do
                 write(ilun)iig
              end do
              do ind=1,twotondim
                 iskip2=ncoarse+(ind-1)*ngridmax2
                 do i=1,ncache
                    iig(i)=son2(igrid2+i+iskip2)
                 end do
                 write(ilun)iig
              end do
              do ind=1,twotondim
                 iskip2=ncoarse+(ind-1)*ngridmax2
                 do i=1,ncache
                    iig(i)=cpu_map(igrid2+i+iskip2)
                 end do
                 write(ilun)iig
              end do
              igrid2=igrid2+ncache
              deallocate(xsp,iig)
           end if
        end do
     end do
     close(ilun)

     deallocate(prev)
     deallocate(next)
     deallocate(prev2)
     deallocate(next2)
     deallocate(mapping)
     deallocate(headl2)
     deallocate(taill2)
     deallocate(flag1)
     deallocate(cpu_map)
     deallocate(xg)
     deallocate(father)
     deallocate(nbor)
     deallocate(son)
     deallocate(father2)
     deallocate(nbor2)
     deallocate(son2)

     ! Now Do particles

     nomfich=TRIM(nomdir)//'part_'//TRIM(nchar)//'.bak'//TRIM(ncharcpu)
     inquire(file=nomfich, exist=ok)
     if (.not. ok) then
        write(*,*)'File '//TRIM(nomfich)//' not found'
        stop
     endif
     ilun=10
     open(unit=ilun,file=nomfich,status='old',form='unformatted')
     ! Write grid variables
     read(ilun)ndim1
     read(ilun)npart
     read(ilun)localseed
     read(ilun)nstar_tot   
     read(ilun)mstar_tot   
     read(ilun)mstar_lost
     
     allocate(xp(1:npart,1:ndim))
     allocate(vp(1:npart,1:ndim))
     allocate(mp(1:npart))
     allocate(idp(1:npart))
     allocate(lvp(1:npart))
     allocate(xsp(1:npart))
     
     do idim=1,ndim
        read(ilun)xp(1:npart,idim)
     end do
     do idim=1,ndim
        read(ilun)vp(1:npart,idim)
     end do
     read(ilun)mp(1:npart)
     read(ilun)idp(1:npart)
     read(ilun)lvp(1:npart)
     close(ilun)


     ilun=11
     nomfich='defrag/part_'//TRIM(nchar)//'.bak'//TRIM(ncharcpu)
     open(unit=ilun,file=nomfich,status='new',form='unformatted')
     ! Write grid variables
     write(ilun)ndim1
     write(ilun)npart
     write(ilun)localseed
     write(ilun)nstar_tot   
     write(ilun)real(mstar_tot)
     write(ilun)real(mstar_lost)
     
     do idim=1,ndim
        xsp=xp(1:npart,idim)
        write(ilun)xsp
     end do
     do idim=1,ndim
        xsp=vp(1:npart,idim)
        write(ilun)xsp
     end do
     xsp=mp(1:npart)
     write(ilun)xsp
     write(ilun)idp
     write(ilun)lvp
     close(ilun)

     deallocate(xp)
     deallocate(vp)
     deallocate(mp)
     deallocate(idp)
     deallocate(lvp)
     deallocate(xsp)
     
  end do

end program dbl2sng

!=======================================================================
subroutine title(n,nchar)
!=======================================================================
  implicit none
  integer::n
  character*5::nchar

  character*1::nchar1
  character*2::nchar2
  character*3::nchar3
  character*4::nchar4
  character*5::nchar5

  if(n.ge.10000)then
     write(nchar5,'(i5)') n
     nchar = nchar5
  elseif(n.ge.1000)then
     write(nchar4,'(i4)') n
     nchar = '0'//nchar4
  elseif(n.ge.100)then
     write(nchar3,'(i3)') n
     nchar = '00'//nchar3
  elseif(n.ge.10)then
     write(nchar2,'(i2)') n
     nchar = '000'//nchar2
  else
     write(nchar1,'(i1)') n
     nchar = '0000'//nchar1
  endif


end subroutine title

