!#########################################################################
!#########################################################################
!#########################################################################
!#########################################################################
subroutine dump_all
  use amr_commons
  use pm_commons
  use hydro_commons
  use cooling_module
  implicit none
#ifndef WITHOUTMPI
  include 'mpif.h'
#endif  
  character(LEN=5)::nchar
  character(LEN=80)::filename,filedir,filecmd
  integer::i,itest,info

  if(nstep_coarse==nstep_coarse_old.and.nstep_coarse>0)return
  if(nstep_coarse==0.and.nrestart>0)return
  if(verbose)write(*,*)'Entering dump_all'

  call write_screen
  call title(ifout,nchar)
  ifout=ifout+1
  if(t>=tout(iout).or.aexp>=aout(iout))iout=iout+1
  output_done=.true.

  if(ndim>1)then
     filedir='output_'//TRIM(nchar)//'/'
     filecmd='mkdir -p '//TRIM(filedir)
#ifdef NOSYSTEM
     call PXFMKDIR(TRIM(filedir),LEN(TRIM(filedir)),O'755',info)
#else
     call system(filecmd)
#endif
#ifndef WITHOUTMPI
     call MPI_BARRIER(MPI_COMM_WORLD,info)
#endif
     filename=TRIM(filedir)//'header_'//TRIM(nchar)//'.txt'
     call output_header(filename)
     if(myid==1)then
        filename=TRIM(filedir)//'info_'//TRIM(nchar)//'.txt'
        call output_info(filename)
        if(cooling)then
           filename=TRIM(filedir)//'cooling_'//TRIM(nchar)//'.out'
           call output_cool(filename)
        end if
        if(sink)then
           filename=TRIM(filedir)//'sink_'//TRIM(nchar)//'.out'
           call output_sink(filename)
           filename=TRIM(filedir)//'sink_'//TRIM(nchar)//'.csv'
           call output_sink_csv(filename)
        endif
     endif
     filename=TRIM(filedir)//'amr_'//TRIM(nchar)//'.out'
     call backup_amr(filename)
     if(hydro)then
        filename=TRIM(filedir)//'hydro_'//TRIM(nchar)//'.out'
        call backup_hydro(filename)
     end if
#ifdef RT
     if(rt.or.neq_chem)then
        filename=TRIM(filedir)//'rt_'//TRIM(nchar)//'.out'
        call rt_backup_hydro(filename)
     endif
#endif
     if(pic)then
        filename=TRIM(filedir)//'part_'//TRIM(nchar)//'.out'
        call backup_part(filename)
     end if
     if(poisson)then
        filename=TRIM(filedir)//'grav_'//TRIM(nchar)//'.out'
        call backup_poisson(filename)
     end if
#ifdef ATON
     if(aton)then
        filename=TRIM(filedir)//'rad_'//TRIM(nchar)//'.out'
        call backup_radiation(filename)
        filename=TRIM(filedir)//'radgpu_'//TRIM(nchar)//'.out'
        call store_radiation(filename)
     end if
#endif
     if (gadget_output) then
        filename=TRIM(filedir)//'gsnapshot_'//TRIM(nchar)
        call savegadget(filename)
     end if
  end if

end subroutine dump_all
!#########################################################################
!#########################################################################
!#########################################################################
!#########################################################################
subroutine backup_amr(filename)
  use amr_commons
  use hydro_commons
  use pm_commons
  implicit none
  character(LEN=80)::filename

  integer::nx_loc,ny_loc,nz_loc,ilun
  integer::ilevel,ibound,ncache,istart,i,igrid,idim,ind,iskip
  integer,allocatable,dimension(:)::ind_grid,iig
  real(dp),allocatable,dimension(:)::xdp
  real(sp),allocatable,dimension(:)::xsp
  real(dp),dimension(1:3)::skip_loc
  character(LEN=80)::fileloc
  character(LEN=5)::nchar
  real(dp)::scale_nH,scale_T2,scale_l,scale_d,scale_t,scale_v
  real(dp)::scale

  if(verbose)write(*,*)'Entering backup_amr'

  ! Conversion factor from user units to cgs units
  call units(scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2)

  ! Local constants
  nx_loc=nx; ny_loc=ny; nz_loc=nz
  if(ndim>0)nx_loc=(icoarse_max-icoarse_min+1)
  if(ndim>1)ny_loc=(jcoarse_max-jcoarse_min+1)
  if(ndim>2)nz_loc=(kcoarse_max-kcoarse_min+1)
  skip_loc=(/0.0d0,0.0d0,0.0d0/)
  if(ndim>0)skip_loc(1)=dble(icoarse_min)
  if(ndim>1)skip_loc(2)=dble(jcoarse_min)
  if(ndim>2)skip_loc(3)=dble(kcoarse_min)
  scale=boxlen/dble(nx_loc)

  !-----------------------------------
  ! Output amr grid in file
  !-----------------------------------  
  ilun=myid+10
  call title(myid,nchar)
  fileloc=TRIM(filename)//TRIM(nchar)
  open(unit=ilun,file=fileloc,form='unformatted')
  ! Write grid variables
  write(ilun)ncpu
  write(ilun)ndim
  write(ilun)nx,ny,nz
  write(ilun)nlevelmax
  write(ilun)ngridmax
  write(ilun)nboundary
  write(ilun)ngrid_current
  write(ilun)boxlen
  ! Write time variables
  write(ilun)noutput,iout,ifout
  write(ilun)tout(1:noutput)
  write(ilun)aout(1:noutput)
  write(ilun)t
  write(ilun)dtold(1:nlevelmax)
  write(ilun)dtnew(1:nlevelmax)
  write(ilun)nstep,nstep_coarse
  write(ilun)const,mass_tot_0,rho_tot
  write(ilun)omega_m,omega_l,omega_k,omega_b,h0,aexp_ini,boxlen_ini
  write(ilun)aexp,hexp,aexp_old,epot_tot_int,epot_tot_old
  write(ilun)mass_sph
  ! Write levels variables
  write(ilun)headl(1:ncpu,1:nlevelmax)
  write(ilun)taill(1:ncpu,1:nlevelmax)
  write(ilun)numbl(1:ncpu,1:nlevelmax)
  write(ilun)numbtot(1:10,1:nlevelmax)
  ! Read boundary linked list
  if(simple_boundary)then
     write(ilun)headb(1:nboundary,1:nlevelmax)
     write(ilun)tailb(1:nboundary,1:nlevelmax)
     write(ilun)numbb(1:nboundary,1:nlevelmax)
  end if
  ! Write free memory
  write(ilun)headf,tailf,numbf,used_mem,used_mem_tot
  ! Write cpu boundaries
  write(ilun)ordering
  if(ordering=='bisection') then
     write(ilun)bisec_wall(1:nbinodes)
     write(ilun)bisec_next(1:nbinodes,1:2)
     write(ilun)bisec_indx(1:nbinodes)
     write(ilun)bisec_cpubox_min(1:ncpu,1:ndim)
     write(ilun)bisec_cpubox_max(1:ncpu,1:ndim)
  else
     write(ilun)bound_key(0:ndomain)
  endif

  ! Write coarse level
  write(ilun)son(1:ncoarse)
  write(ilun)flag1(1:ncoarse)
  write(ilun)cpu_map(1:ncoarse)
  ! Write fine levels
  do ilevel=1,nlevelmax
     do ibound=1,nboundary+ncpu
        if(ibound<=ncpu)then
           ncache=numbl(ibound,ilevel)
           istart=headl(ibound,ilevel)
        else
           ncache=numbb(ibound-ncpu,ilevel)
           istart=headb(ibound-ncpu,ilevel)
        end if
        if(ncache>0)then
           allocate(ind_grid(1:ncache),xdp(1:ncache),iig(1:ncache))
           ! Write grid index
           igrid=istart
           do i=1,ncache
              ind_grid(i)=igrid
              igrid=next(igrid)
           end do
           write(ilun)ind_grid
           ! Write next index
           do i=1,ncache
              iig(i)=next(ind_grid(i))
           end do
           write(ilun)iig
           ! Write prev index
           do i=1,ncache
              iig(i)=prev(ind_grid(i))
           end do
           write(ilun)iig
           ! Write grid center
           do idim=1,ndim
              do i=1,ncache
                 xdp(i)=xg(ind_grid(i),idim)
              end do
              write(ilun)xdp
           end do
           ! Write father index
           do i=1,ncache
              iig(i)=father(ind_grid(i))
           end do
           write(ilun)iig
           ! Write nbor index
           do ind=1,twondim
              do i=1,ncache
                 iig(i)=nbor(ind_grid(i),ind)
              end do
              write(ilun)iig
           end do
           ! Write son index
           do ind=1,twotondim
              iskip=ncoarse+(ind-1)*ngridmax
              do i=1,ncache
                 iig(i)=son(ind_grid(i)+iskip)
              end do
              write(ilun)iig
           end do
           ! Write cpu map
           do ind=1,twotondim
              iskip=ncoarse+(ind-1)*ngridmax
              do i=1,ncache
                 iig(i)=cpu_map(ind_grid(i)+iskip)
              end do
              write(ilun)iig
           end do
           ! Write refinement map
           do ind=1,twotondim
              iskip=ncoarse+(ind-1)*ngridmax
              do i=1,ncache
                 iig(i)=flag1(ind_grid(i)+iskip)
              end do
              write(ilun)iig
           end do
           deallocate(xdp,iig,ind_grid)
        end if
     end do
  end do
  close(ilun)
  
end subroutine backup_amr
!#########################################################################
!#########################################################################
!#########################################################################
!#########################################################################
subroutine output_info(filename)
  use amr_commons
  use hydro_commons
  use pm_commons
  implicit none
  character(LEN=80)::filename

  integer::nx_loc,ny_loc,nz_loc,ilun,icpu,idom
  real(dp)::scale
  real(dp)::scale_nH,scale_T2,scale_l,scale_d,scale_t,scale_v
  character(LEN=80)::fileloc
  character(LEN=5)::nchar

  if(verbose)write(*,*)'Entering output_info'

  ilun=myid+10

  ! Conversion factor from user units to cgs units
  call units(scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2)

  ! Local constants
  nx_loc=nx; ny_loc=ny; nz_loc=nz
  if(ndim>0)nx_loc=(icoarse_max-icoarse_min+1)
  if(ndim>1)ny_loc=(jcoarse_max-jcoarse_min+1)
  if(ndim>2)nz_loc=(kcoarse_max-kcoarse_min+1)
  scale=boxlen/dble(nx_loc)

  ! Open file
  fileloc=TRIM(filename)
  open(unit=ilun,file=fileloc,form='formatted')
  
  ! Write run parameters
  write(ilun,'("ncpu        =",I11)')ncpu
  write(ilun,'("ndim        =",I11)')ndim
  write(ilun,'("levelmin    =",I11)')levelmin
  write(ilun,'("levelmax    =",I11)')nlevelmax
  write(ilun,'("ngridmax    =",I11)')ngridmax
  write(ilun,'("nstep_coarse=",I11)')nstep_coarse
  write(ilun,*)

  ! Write physical parameters
  write(ilun,'("boxlen      =",E23.15)')scale
  write(ilun,'("time        =",E23.15)')t
  write(ilun,'("aexp        =",E23.15)')aexp
  write(ilun,'("H0          =",E23.15)')h0
  write(ilun,'("omega_m     =",E23.15)')omega_m
  write(ilun,'("omega_l     =",E23.15)')omega_l
  write(ilun,'("omega_k     =",E23.15)')omega_k
  write(ilun,'("omega_b     =",E23.15)')omega_b
  write(ilun,'("unit_l      =",E23.15)')scale_l
  write(ilun,'("unit_d      =",E23.15)')scale_d
  write(ilun,'("unit_t      =",E23.15)')scale_t
  write(ilun,*)
  
  ! Write ordering information
  write(ilun,'("ordering type=",A80)')ordering
  if(ordering=='bisection') then
     do icpu=1,ncpu
        ! write 2*ndim floats for cpu bound box                                                         
        write(ilun,'(E23.15)')bisec_cpubox_min(icpu,:),bisec_cpubox_max(icpu,:)
        ! write 1 float for cpu load                                                                    
        write(ilun,'(E23.15)')dble(bisec_cpu_load(icpu))
     end do
  else
     write(ilun,'("   DOMAIN   ind_min                 ind_max")')
     do idom=1,ndomain
        write(ilun,'(I8,1X,E23.15,1X,E23.15)')idom,bound_key(idom-1),bound_key(idom)
     end do
  endif

  close(ilun)

end subroutine output_info
!#########################################################################
!#########################################################################
!#########################################################################
!#########################################################################
subroutine output_header(filename)
  use amr_commons
  use hydro_commons
  use pm_commons
  implicit none
#ifndef WITHOUTMPI
  include 'mpif.h'
#endif
  character(LEN=80)::filename

  integer::info,npart_tot,ilun
  character(LEN=80)::fileloc

  if(verbose)write(*,*)'Entering output_header'

  ! Compute total number of particles
#ifndef WITHOUTMPI
  call MPI_ALLREDUCE(npart,npart_tot,1,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,info)
#endif
#ifdef WITHOUTMPI
  npart_tot=npart
#endif

  if(myid==1)then

     ilun=myid+10

     ! Open file
     fileloc=TRIM(filename)
     open(unit=ilun,file=fileloc,form='formatted')
     
     ! Write header information
     write(ilun,*)'Total number of particles'
     write(ilun,*)npart_tot
     write(ilun,*)'Total number of dark matter particles'
     write(ilun,*)npart_tot-nstar_tot
     write(ilun,*)'Total number of star particles'
     write(ilun,*)nstar_tot
     write(ilun,*)'Total number of sink particles'
     write(ilun,*)nsink

     ! Keep track of what particle fields are present
     write(ilun,*)'Particle fields'
     write(ilun,'(a)',advance='no')'pos vel mass iord level '
#ifdef OUTPUT_PARTICLE_POTENTIAL
     write(ilun,'(a)',advance='no')'phi '
#endif
     if(star.or.sink) then
        write(ilun,'(a)',advance='no')'tform '
        if(metal) then
           write(ilun,'(a)',advance='no')'metal '
        endif
     endif
     close(ilun)

  endif

end subroutine output_header
!#########################################################################
!#########################################################################
!#########################################################################
!#########################################################################
subroutine savegadget(filename)
  use amr_commons
  use hydro_commons
  use pm_commons
  use gadgetreadfilemod
  implicit none
#ifndef WITHOUTMPI
  include 'mpif.h'
#endif
  character(LEN=80)::filename
  TYPE (gadgetheadertype) :: header
  real,allocatable,dimension(:,:)::pos, vel
  integer,allocatable,dimension(:)::ids
  integer::i, idim, ipart
  real:: gadgetvfact
  integer::npart_tot, info
  real, parameter:: RHOcrit = 2.7755d11

#ifndef WITHOUTMPI
  call MPI_ALLREDUCE(npart,npart_tot,1,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,info)
#endif

  allocate(pos(ndim, npart), vel(ndim, npart), ids(npart))
  gadgetvfact = 100.0 * boxlen_ini / aexp / SQRT(aexp)

  header%npart = 0
  header%npart(2) = npart
  header%mass = 0
  header%mass(2) = omega_m*RHOcrit*(boxlen_ini)**3/npart_tot/1.d10
  header%time = aexp
  header%redshift = 1.d0/aexp-1.d0
  header%flag_sfr = 0
  header%nparttotal = 0
  header%nparttotal(2) = npart_tot
  header%flag_cooling = 0
  header%numfiles = ncpu
  header%boxsize = boxlen_ini
  header%omega0 = omega_m
  header%omegalambda = omega_l
  header%hubbleparam = h0/100.0
  header%flag_stellarage = 0
  header%flag_metals = 0
  header%totalhighword = 0
  header%flag_entropy_instead_u = 0
  header%flag_doubleprecision = 0
  header%flag_ic_info = 0
  header%lpt_scalingfactor = 0
  header%unused = ' '

  do idim=1,ndim
     ipart=0
     do i=1,npartmax
        if(levelp(i)>0)then
           ipart=ipart+1
           if (ipart .gt. npart) then
                write(*,*) myid, "Ipart=",ipart, "exceeds", npart
                call clean_stop
           endif
           pos(idim, ipart)=xp(i,idim) * boxlen_ini
           vel(idim, ipart)=vp(i,idim) * gadgetvfact
           if (idim.eq.1) ids(ipart) = idp(i)
        end if
     end do
  end do

  call gadgetwritefile(filename, myid-1, header, pos, vel, ids)
  deallocate(pos, vel, ids)

end subroutine savegadget

