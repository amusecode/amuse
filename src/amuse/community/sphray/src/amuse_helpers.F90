module amuse_sphrayMod
use myf03_mod
use global_mod, only: psys, PLAN,GV
use particle_system_mod, only: particle_type,source_type
implicit none
  
type(particle_type), allocatable :: par_buffer(:)
type(source_type), allocatable :: src_buffer(:)
  
integer(i8b) :: npar_buffer=0,nsrc_buffer=0
integer(i8b) :: nmaxpar_buffer=100000
integer(i8b) :: nmaxsrc_buffer=10000

integer(i8b) :: tot_id=0
integer(i8b) :: ngas=0,nsrc=0
integer, allocatable :: gid(:),sid(:)
logical :: gas_searcheable=.FALSE.
logical :: src_searcheable=.FALSE.

character(len=clen) :: data_directory='./'

!< [Config File] file containing spectra tables
character(len=clen) :: SpectraFile='./spectra/thermal1e5.cdf'
!< [Config File] file containing b2cd tables
character(len=clen) :: b2cdFile='./column_depth/latmon_b2cd_table.txt'
!< [Config File] file containnig atomic rates tables     
character(len=clen) :: AtomicRatesFile='./atomic_rates/atomic_rates_Hui.txt'
 
  real(r8b) :: RaysperUnitTime=1.e6

real(r8b) :: default_spectral_type = -1.




 contains




subroutine sphray_set_data_directory(x)
  character(len=clen) :: x
  data_directory=x
end subroutine  

subroutine sphray_get_data_directory(x)
  character(len=clen) :: x
  x=data_directory
end subroutine  

subroutine sphray_set_output_directory(x)
  character(len=clen) :: x
  GV%OutputDir=x
end subroutine  

subroutine sphray_get_output_directory(x)
  character(len=clen) :: x
  x=GV%OutputDir
end subroutine  

function new_id()
  integer new_id
  tot_id=tot_id+1
  new_id=tot_id
end function

subroutine sphray_init
  use myf03_mod
  use mt19937_mod, only: init_mersenne_twister
  
  integer, parameter :: verb = 3    !< verbosity before config file is read

  myf03_verbosity = verb
  
  call set_default_parameters()
  
  call init_mersenne_twister(GV%IntSeed)
  
  call extend_par(par_buffer,nmaxpar_buffer)
  call extend_src(src_buffer,nmaxsrc_buffer)
    
end subroutine

function clean_gas(par) result(np)
  integer :: left,right,np
  type(particle_type), allocatable :: par(:)
  type(particle_type) :: tmp
  left=1
  if(.NOT.allocated(par)) then
    np = 0
    return
  endif
  right=size(par)
  if(right.EQ.0) then
    np=0
    return
  endif    
  do while(.TRUE.)
    do while(par(left)%mass.GE.0.AND.left.LT.right)
      left=left+1
    enddo
    do while(par(right)%mass.LT.0.AND.left.LT.right)
      right=right-1
    enddo
    if(left.LT.right) then
      tmp=par(left)
      par(left)=par(right)
      par(right)=tmp
    else
      exit
    endif  
  enddo  
  if(par(left)%mass.GE.0) left=left+1
  np=left-1
end function

function clean_src(par) result(np)
  integer :: left,right,np
  type(source_type), allocatable :: par(:)
  type(source_type) :: tmp
  left=1
  if(.NOT.allocated(par)) then
    np = 0
    return
  endif
  right=size(par)
  if(right.EQ.0) then
    np=0
    return
  endif    
  do while(.TRUE.)
    do while(par(left)%L.GE.0.AND.left.LT.right)
      left=left+1
    enddo
    do while(par(right)%L.LT.0.AND.left.LT.right)
      right=right-1
    enddo
    if(left<right) then
      tmp=par(left)
      par(left)=par(right)
      par(right)=tmp
    else
      exit
    endif  
  enddo  
  if(par(left)%L.GE.0) left=left+1
  np=left-1
end function

subroutine sphray_commit_particles
  integer :: n

 gas_searcheable=.FALSE.
 src_searcheable=.FALSE. 
 
 n=clean_gas(psys%par)
 if(npar_buffer.GT.0.OR.n.NE.size(psys%par)) then
   call extend_par(psys%par,n+npar_buffer)
   psys%par(n+1:n+npar_buffer)=par_buffer(1:npar_buffer)
   npar_buffer=0
 endif

 n=clean_src(psys%src)
 if(nsrc_buffer.GT.0.OR.n.NE.size(psys%src)) then
   call extend_src(psys%src,n+nsrc_buffer)
   psys%src(n+1:n+nsrc_buffer)=src_buffer(1:nsrc_buffer)
   nsrc_buffer=0
 endif

end subroutine

subroutine sphray_commit_parameters
  use global_mod, only: rtable,XHII_k,isoT_k
  use b2cd_mod, only: read_b2cd_file
  use spectra_mod, only: read_spectra_file
  use atomic_rates_mod, only: read_atomic_rates_file, get_atomic_rates
  use atomic_rates_mod, only: write_atomic_rates_to_log_file
  
  !< [Config File] file containing spectra tables
  GV%SpectraFile=trim(data_directory)//'/'//trim(SpectraFile) 
  !< [Config File] file containing b2cd tables
  GV%b2cdFile=trim(data_directory)//'/'//trim(b2cdFile)           
  !< [Config File] file containnig atomic rates tables
  GV%AtomicRatesFile=trim(data_directory)//'/'//trim(AtomicRatesFile) 
  
  call read_b2cd_file(GV%b2cdFile)
  call read_spectra_file(GV%SpectraFile)

  call read_atomic_rates_file(rtable, GV%AtomicRatesFile)
  call write_atomic_rates_to_log_file(rtable, GV%OutputDir)

  call get_atomic_rates(1.0d4, rtable, xHII_k)
  if (GV%IsoTemp > 0.0) then
    call get_atomic_rates(GV%IsoTemp, rtable, isoT_k)
  end if

  call amuse_planning_data
  call amuse_output_planning
  call amuse_ray_planning
  call amuse_initialize_global_variables

end subroutine

subroutine amuse_initialize_global_variables
  use gadget_General_class, only: gadget_constants_type
  use particle_system_mod, only: return_bytes_per_source
  use particle_system_mod, only: return_bytes_per_particle

  character(clen), parameter :: myname="initialize_global_variables"
  logical, parameter :: crash=.true.
  integer, parameter :: verb=2
  character(clen) :: str

  real, parameter :: zero = 0.0d0
  integer(i8b) :: nrays
  integer(i8b) :: i

  type(gadget_constants_type) :: gconst

  call mywrite("   initializing global variables",verb) 


  ! calculate total number of rays to be traced
  !----------------------------------------------------------------------
  nrays=GV%ForcedRayNumber

  ! report ionization solver being used and particle/source mem footprint
  !----------------------------------------------------------------------
  if (GV%IonTempSolver == 1) then
     call mywrite("   using Implicit Euler ionization solver",verb)
  else if (GV%IonTempSolver == 2) then
     call mywrite("   using Backwards Difference ionization solver",verb)
  end if

  GV%bytesperpar = return_bytes_per_particle()
  GV%bytespersrc = return_bytes_per_source()

  write(str,'(A,I4)') "   bytes per particle = ", GV%bytesperpar
  call mywrite(str, verb)

  write(str,'(A,I4)') "   bytes per source = ", GV%bytespersrc
  call mywrite(str, verb)

  ! open log files
  !-----------------
  GV%ionfrac_file = trim(GV%OutputDir) // "/ionfrac.log"
  call open_formatted_file_w(GV%ionfrac_file,GV%ionlun)

  GV%pardata_file = trim(GV%OutputDir) // "/particle_data.log"
  call open_formatted_file_w(GV%pardata_file,GV%pardatalun)

  GV%srcdata_file = trim(GV%OutputDir) // "/source_data.log"
  call open_formatted_file_w(GV%srcdata_file,GV%srcdatalun)


                 
  ! initialize counters and timers
  !---------------------------------
  GV%CurSnapNum = GV%StartSnapNum
  GV%rayn = 0
  GV%src_rayn = 0
  GV%MB = zero
 
  GV%itime = 0_i8b
  GV%start_time_code = 0.
  GV%start_time_s    = GV%start_time_code * GV%cgs_time / GV%LittleH
  GV%start_time_myr  = GV%start_time_s / gconst%sec_per_megayear

  GV%time_elapsed_code = zero
  GV%time_elapsed_s    = zero
  GV%time_elapsed_myr  = zero
      
  GV%TotalSourceRaysCast = zero
  GV%TotalDiffuseRaysCast = zero
  GV%IonizingPhotonsPerSec = zero
  GV%TotalPhotonsCast = zero
  GV%TotalPhotonsAbsorbed = zero
  GV%PhotonsLeavingBox = zero
  GV%TotalIonizations = zero
  GV%TotalRecombinations = zero
  
  GV%PeakUpdates = zero
  GV%AverageUpdatesPerPar = zero
  GV%ParticleCrossings = zero
  GV%TotalDerivativeCalls = zero
  
  call mywrite("",verb) 

end subroutine

subroutine amuse_planning_data()
  use gadget_general_class, only: gadget_units_type
  type(gadget_units_type) :: gunits

! one plan
  GV%Nsnaps = 1
  if (allocated(PLAN%snap)) deallocate(PLAN%snap)
  allocate( PLAN%snap(1) ) 
  
  PLAN%snap(1)%TimeAt=0.
  PLAN%snap(1)%ScaleFacAt=1.
  PLAN%snap(1)%StartTime=PLAN%snap(1)%TimeAt
  PLAN%snap(1)%SrcRays = 0 ! GV%ForcedRayNumber

! standard gadget units
  GV%cgs_len  = gunits%cgs_length
  GV%cgs_mass = gunits%cgs_mass
  GV%cgs_vel  = gunits%cgs_velocity
  GV%cgs_time = gunits%cgs_time
  GV%cgs_rho  = gunits%cgs_density
  GV%cgs_prs  = gunits%cgs_pressure
  GV%cgs_enrg = gunits%cgs_energy 

! inconsequential??
  GV%LenFac_cm=0.  !< code (input) length -> cm = cgs_len * a / h
  GV%MassFac_g=0.  !< code (input) mass -> g = cgs_mass / h
  GV%TimeFac_s=0.  !< code (input) time -> s = cgs_time / h
  GV%OmegaM=0.    !< matter / critical density z=0
  GV%OmegaB=0.    !< baryon / critical density z=0
  GV%OmegaL=0.    !< lambda / critical density z=0
  GV%LittleH=1.   !< Hubble parameter z=0 in units of 100 km/s/Mpc
! ----------------- 
 
  GV%Lunit=1.d50                  !< code source luminosity unit [photons/s]

end subroutine

subroutine amuse_output_planning
! equiv of do_output_planning
 allocate( PLAN%OutputTimes(0:1) )
 GV%NumTotOuts = GV%NumStdOuts  ! =-1
 PLAN%OutputTimes=1.d99   ! probably doesn't matter
! deferred to just before evolve 
end subroutine

subroutine amuse_ray_planning
! equiv of do_ray_planning
! I don;t think this has to do something
! deferred to just before evolve 
end subroutine

subroutine sphray_end

end subroutine 

function sphray_get_gas_particle_state(id,mass,hsml,x,y,z,rho,xe,u) result(ret)
  integer(i4b) :: id
  real(r4b) :: mass, hsml,x,y,z, &
    rho, u , xe 
  integer(i4b) :: i,index,ret

    index=find_gas(id)
    if(index.LT.0) then
      ret=index
      return
    endif  
        
    mass=psys%par(index)%mass
    x=psys%par(index)%pos(1)
    y=psys%par(index)%pos(2)
    z=psys%par(index)%pos(3)
    hsml=psys%par(index)%hsml
    rho=psys%par(index)%rho
    xe=psys%par(index)%xHII
    u=u_from_temp( real(psys%par(index)%T, r8b) , &
                   real(psys%par(index)%ye,r8b) , &
                   GV%H_mf ) 
    ret=0

end function


function sphray_set_gas_particle_state(ids,mass,hsml,x,y,z,rho,xe,u) result(ret)
  integer(i4b) ::  ids
  real(r4b) :: mass, hsml,x,y,z, &
    rho, u , xe 
  integer(i4b) :: i,index,ret

    index=find_gas(ids)
    if(index.LT.0) then
      ret=index
      return
    endif  
    psys%par(index)%mass=mass
    psys%par(index)%pos(1)=x
    psys%par(index)%pos(2)=y
    psys%par(index)%pos(3)=z
    psys%par(index)%hsml=hsml
    psys%par(index)%rho=rho
    psys%par(index)%ye=xe
    psys%par(index)%T=temp_from_u( real(u,r8b) , real(xe,r8b), GV%H_mf)  

    ret=0

end function


function sphray_set_gas_particle_pos(ids,x,y,z) result(ret)
  integer(i4b) ::  ids
  real(r4b) :: mass, hsml,x,y,z, &
    rho, u , xe 
  integer(i4b) :: i,index,ret

    index=find_gas(ids)
    if(index.LT.0) then
      ret=index
      return
    endif  
    psys%par(index)%pos(1)=x
    psys%par(index)%pos(2)=y
    psys%par(index)%pos(3)=z

    ret=0

end function


function sphray_set_gas_particle_rho(ids,rho) result(ret)
  integer(i4b) ::  ids
  real(r4b) :: mass, hsml,x,y,z, &
    rho, u , xe 
  integer(i4b) :: i,index,ret

    index=find_gas(ids)
    if(index.LT.0) then
      ret=index
      return
    endif  
    psys%par(index)%rho=rho

    ret=0

end function


function sphray_set_gas_particle_hsml(ids,hsml) result(ret)
  integer(i4b) ::  ids
  real(r4b) :: mass, hsml,x,y,z, &
    rho, u , xe 
  integer(i4b) :: i,index,ret

    index=find_gas(ids)
    if(index.LT.0) then
      ret=index
      return
    endif  
    psys%par(index)%hsml=hsml

    ret=0

end function


function sphray_set_gas_particle_u(ids,u) result(ret)
  integer(i4b) ::  ids
  real(r4b) :: mass, hsml,x,y,z, &
    rho, u , xe 
  integer(i4b) :: i,index,ret

    index=find_gas(ids)
    if(index.LT.0) then
      ret=index
      return
    endif  
    psys%par(index)%T=temp_from_u( real(u,r8b) , real(xe,r8b), GV%H_mf)  

    ret=0

end function


function sphray_set_gas_particle_vel(ids,vx,vy,vz) result(ret)
  integer(i4b) ::  ids
  real(r4b) :: mass, hsml, vx, vy, vz, rho, u , xe 
  integer(i4b) :: i,index,ret

  ret=0
#ifdef incVel
    index=find_gas(ids)
    if(index.LT.0) then
      ret=index
      return
    endif  
    psys%par(index)%vel(1)=vx
    psys%par(index)%vel(2)=vy
    psys%par(index)%vel(3)=vz
    ret=0
#endif

end function


function sphray_get_gas_particle_vel(ids,vx,vy,vz) result(ret)
  integer(i4b) ::  ids
  real(r4b) :: mass, hsml,vx,vy,vz, rho, u , xe 
  integer(i4b) :: i,index,ret

  ret=0
  vx=0;vy=0;vz=0
#ifdef incVel
    index=find_gas(ids)
    if(index.LT.0) then
      ret=index
      return
    endif  
    vx=psys%par(index)%vel(1)
    vy=psys%par(index)%vel(2)
    vz=psys%par(index)%vel(3)
    ret=0
#endif

end function



function sphray_set_gas_particle_dudt(ids,dudt) result(ret)
  integer(i4b) ::  ids
  real(r4b) :: dudt 
  integer(i4b) :: i,index,ret

  ret=0
#ifdef incHeat
    index=find_gas(ids)
    if(index.LT.0) then
      ret=index
      return
    endif  
    psys%par(index)%dinternalEdt=dudt
    ret=0
#endif

end function


function sphray_get_gas_particle_dudt(ids,dudt) result(ret)
  integer(i4b) ::  ids
  real(r4b) :: dudt 
  integer(i4b) :: i,index,ret

  ret=0
  dudt=0
#ifdef incHeat
    index=find_gas(ids)
    if(index.LT.0) then
      ret=index
      return
    endif  
    dudt=psys%par(index)%dinternalEdt
    ret=0
#endif

end function


function sphray_get_src_particle_state(id,L,x,y,z,SpcType) result(ret)
  integer(i4b) :: id
  real(r4b) :: L,x,y,z,SpcType
  integer(i4b) :: i,index,ret

    index=find_src(id)
    if(index.LT.0) then
      ret=index
      return
    endif  
    
    L=psys%src(index)%L
    x=psys%src(index)%pos(1)
    y=psys%src(index)%pos(2)
    z=psys%src(index)%pos(3)
    SpcType=psys%src(index)%SpcType
    ret=0

end function

function sphray_set_src_particle_state(id,L,x,y,z,SpcType) result(ret)
  integer(i4b) ::  id
  real(r4b) :: L,x,y,z,SpcType
  integer(i4b) :: i,index,ret

    index=find_src(id)
    if(index.LT.0) then
      ret=index
      return
    endif  
        
    psys%src(index)%L=L
    psys%src(index)%pos(1)=x
    psys%src(index)%pos(2)=y
    psys%src(index)%pos(3)=z
    if(SpcType.NE.0) then
      psys%src(index)%SpcType=SpcType
    else
      psys%src(index)%SpcType=default_spectral_type
    endif
    ret=0
end function


subroutine sphray_add_gas_particle(id,mass,hsml,x,y,z,rho,xe,u,vx,vy,vz)
  use particle_system_mod, only: particle_set_ye, particle_set_ci_eq
  integer(i4b) ::  id
  real(r4b) :: mass, hsml, x, y, z, rho, u , xe, vx, vy, vz 
  integer(i4b) :: i
  logical :: DoH
  logical :: DoHe 
  logical :: caseA(2)

  gas_searcheable=.FALSE.
  
  do while(npar_buffer+1.GT.nmaxpar_buffer) 
    nmaxpar_buffer=2*nmaxpar_buffer
    call extend_par(par_buffer,nmaxpar_buffer)
  enddo
  
  ngas=ngas+1
  id=new_id()
  npar_buffer=npar_buffer+1
  par_buffer(npar_buffer)%mass=mass
  par_buffer(npar_buffer)%id=id
  par_buffer(npar_buffer)%pos(1)=x
  par_buffer(npar_buffer)%pos(2)=y
  par_buffer(npar_buffer)%pos(3)=z
#ifdef incVel
  par_buffer(npar_buffer)%vel(1)=vx
  par_buffer(npar_buffer)%vel(2)=vy
  par_buffer(npar_buffer)%vel(3)=vz
#endif
  par_buffer(npar_buffer)%hsml=hsml
  par_buffer(npar_buffer)%rho=rho
  par_buffer(npar_buffer)%T=temp_from_u( real(u,r8b) , real(xe,r8b), GV%H_mf)   
  par_buffer(npar_buffer)%ye=xe   
  par_buffer(npar_buffer)%xHI=1-xe   
  par_buffer(npar_buffer)%xHII=xe
#ifdef incHeat  
  par_buffer(npar_buffer)%dinternalEdt=0.
#endif

!  caseA = .false.
!  if (GV%HydrogenCaseA) caseA(1) = .true.
!  if (GV%HeliumCaseA)   caseA(2) = .true.

!  DoH = .true.
!#ifdef incHe
!  DoHe = .true.
!#else
!  DoHe = .false.
!#endif

!  call particle_set_ci_eq( par_buffer(npar_buffer), caseA, DoH, DoHe, fit='hui' )
!  call particle_set_ye( par_buffer(npar_buffer), GV%H_mf, GV%He_mf, GV%NeBackground )

end subroutine

subroutine sphray_add_src_particle(id,L,x,y,z,SpcType)
  integer(i4b) ::  id
  real(r4b) :: L,  x,y,z,SpcType
  integer(i4b) :: i

  src_searcheable=.FALSE. 
  
  do while(nsrc_buffer+1.GT.nmaxsrc_buffer)
    nmaxsrc_buffer=nmaxsrc_buffer*2 
    call extend_src(src_buffer,nmaxsrc_buffer)
  enddo 
    
    nsrc=nsrc+1
    id=new_id()
    nsrc_buffer=nsrc_buffer+1
    src_buffer(nsrc_buffer)%L=L
    src_buffer(nsrc_buffer)%id=id   
    src_buffer(nsrc_buffer)%pos(1)=x
    src_buffer(nsrc_buffer)%pos(2)=y
    src_buffer(nsrc_buffer)%pos(3)=z
#ifdef incVel
    src_buffer(nsrc_buffer)%vel(1)=0.
    src_buffer(nsrc_buffer)%vel(2)=0.
    src_buffer(nsrc_buffer)%vel(3)=0.
#endif
    if(SpcType.NE.0) then
      src_buffer(nsrc_buffer)%SpcType=SpcType
    else
      src_buffer(nsrc_buffer)%SpcType=default_spectral_type
    endif
    src_buffer(nsrc_buffer)%EmisPrf=0
    src_buffer(nsrc_buffer)%lastemit=GV%rayn
    
end subroutine









function sphray_get_wall_src_state(id,L,SpcType,EmisPrf) result(ret)
  integer(i4b) :: id
  real(r4b) :: L,SpcType
  integer(i4b) :: i,index,ret,EmisPrf

    index=find_src(id)
    if(index.LT.0) then
      ret=index
      return
    endif  
    
    L=psys%src(index)%L
    SpcType=psys%src(index)%SpcType
    EmisPrf=psys%src(index)%EmisPrf
    ret=0

end function


function sphray_set_wall_src_state(id,L,SpcType,EmisPrf) result(ret)
  integer(i4b) ::  id
  real(r4b) :: L,SpcType
  integer(i4b) :: i,index,ret,EmisPrf

    index=find_src(id)
    if(index.LT.0) then
      ret=index
      return
    endif  
        
    psys%src(index)%L=L
    if(SpcType.NE.0) then
      psys%src(index)%SpcType=SpcType
    else
      psys%src(index)%SpcType=default_spectral_type
    endif

    psys%src(index)%EmisPrf=EmisPrf

    ret=0
end function



subroutine sphray_add_wall_src(id,F,SpcType,EmisPrf)
  integer(i4b) ::  id
  real(r4b) :: F,L,SpcType
  integer(i4b) :: i
  integer(i4b) :: EmisPrf

  src_searcheable=.FALSE. 
  
  do while(nsrc_buffer+1.GT.nmaxsrc_buffer)
    nmaxsrc_buffer=nmaxsrc_buffer*2 
    call extend_src(src_buffer,nmaxsrc_buffer)
  enddo 

  nsrc=nsrc+1
  id=new_id()
  nsrc_buffer=nsrc_buffer+1
    
  ! for a wall source we take F = photons / cm^2 as input
  ! and convert it to L = photons / s using the wall size
  ! this assumes the box is square

  L = F * psys%box%lens_cm(1)**2

  src_buffer(nsrc_buffer)%L=L
  src_buffer(nsrc_buffer)%id=id   

  ! position should never be used for a wall source 

  src_buffer(nsrc_buffer)%pos(1) = 0.0d0
  src_buffer(nsrc_buffer)%pos(2) = 0.0d0
  src_buffer(nsrc_buffer)%pos(3) = 0.0d0

  if(SpcType.NE.0) then
     src_buffer(nsrc_buffer)%SpcType=SpcType
  else
      src_buffer(nsrc_buffer)%SpcType=default_spectral_type
   endif
   src_buffer(nsrc_buffer)%EmisPrf=EmisPrf
   src_buffer(nsrc_buffer)%lastemit=GV%rayn
    
end subroutine



function sphray_remove_gas_particle(id) result(ret)
  integer(i4b) :: id
  integer(i4b) :: index,ret

  index=find_gas(id)
  if(index.LT.0) then
    ret=index
    return
  endif  

  psys%par(index)%mass=-1
  ngas=ngas-1
  ret=0
end function

function sphray_remove_src_particle(id) result(ret)
  integer(i4b) :: id
  integer(i4b) :: index,ret

  index=find_src(id)
  if(index.LT.0) then
    ret=index
    return
  endif  

  psys%src(index)%L=-1
  nsrc=nsrc-1
  ret=0
end function

subroutine sphray_evolve(tend)
  use particle_system_mod, only: particle_system_set_ye
  use amuse_mainloop_mod
  real(r8b) :: tend

  if(tend.LE.PLAN%snap(1)%TimeAt) return 
  PLAN%snap(1)%TimeToNext=tend
  PLAN%snap(1)%RunTime=tend-PLAN%snap(1)%TimeAt
  
  GV%ForcedRayNumber=RaysperUnitTime*(tend-PLAN%snap(1)%TimeAt)
  if(GV%ForcedRayNumber.LE.0) call myerr("no rays to trace??","sphray_evolve",.true.)
  
  PLAN%snap(1)%StartTime=PLAN%snap(1)%TimeAt
  PLAN%snap(1)%SrcRays = GV%ForcedRayNumber

  call preparemain
  call mainloop

  call particle_system_set_ye(psys, GV%H_mf, GV%He_mf, GV%NeBackground)

  PLAN%snap(1)%TimeAt=tend
  gas_searcheable=.FALSE.
  src_searcheable=.FALSE. 
    
end subroutine

function sphray_model_time() result(time)
  use amuse_mainloop_mod
  real(r8b) :: time
  time=PLAN%snap(1)%TimeAt
end function

subroutine sphray_set_model_time(time) 
  use amuse_mainloop_mod
  real(r8b) :: time
  PLAN%snap(1)%TimeAt=time   ! is this safe? 
end subroutine


subroutine preparemain()
  use gadget_General_class, only: gadget_constants_type
  use atomic_rates_mod, only: read_atomic_rates_file, get_atomic_rates
  use global_mod, only: set_dt_from_dtcode,rtable,cmbT_k
  use particle_system_mod, only: particle_system_enforce_x_and_T_minmax
  use source_input_mod, only: order_sources_lum


  type(gadget_constants_type) :: gconst
  real(r8b) :: a=1.     !< scale factor
  real(r8b) :: h=1.     !< Hubble paraemter (little H)
  integer(i4b) :: i

! most of this stuff should go into (re)commit parameters/particles

  psys%box%tops = GV%BoxUprs
  psys%box%bots = GV%BoxLwrs

  psys%box%lens    = GV%BoxUprs - GV%BoxLwrs
  psys%box%lens_cm = psys%box%lens * GV%cgs_len

  psys%box%vol    = product( psys%box%lens )
  psys%box%vol_cm = product( psys%box%lens_cm )

  psys%box%tbound = GV%BndryCond
  psys%box%bbound = GV%BndryCond

  if(size(psys%src).GT.0) call order_sources_lum(psys%src)
!  psys%src%lastemit = GV%rayn
! to recommit parameter

!  psys%src%lastemit = 0

!from readin_snaps
#ifdef outGammaHI
  psys%par(:)%gammaHI = 0.0
  psys%par(:)%time = 0.0
#endif

  ! scale the data if we need to
  !=====================================================================

  ! convert number density to flux for planar sources
  !==========================================================

  ! set EOS particles to EOS temp if you want
  !=======================================================

  ! set SFR particles to EOS temp if you want
  !=======================================================


  ! set constant temperature if we have one
  !=======================================================
  if (GV%IsoTemp > 0.0) psys%par(:)%T = GV%IsoTemp

  ! cap the ionization fractions and temperatures if we have to
  !================================================================
  call particle_system_enforce_x_and_T_minmax( &
       psys, GV%xfloor, GV%xceiling, GV%Tfloor, GV%Tceiling )

  ! and the rest of the stuff
  !===============================================================
  GV%dt_code = PLAN%snap(GV%CurSnapNum)%RunTime / PLAN%snap(GV%CurSnapNum)%SrcRays
  call set_dt_from_dtcode( GV )
  
  GV%Tcmb_cur = gconst%t_cmb0 / a
  call get_atomic_rates(GV%Tcmb_cur, rtable, cmbT_k)

  GV%total_mass = 0.0d0
  do i = 1,size(psys%par)
     GV%total_mass = GV%total_mass + psys%par(i)%mass
  end do
  
  GV%total_lum = 0.0d0
  do i = 1,size(psys%src)
     GV%total_lum = GV%total_lum + psys%src(i)%L
  end do
  
  GV%total_atoms = GV%total_mass * GV%cgs_mass * &
       (GV%H_mf  / (gconst%protonmass) + &
       GV%He_mf / (4*gconst%protonmass) )


  GV%total_photons = (GV%TotalSimTime * GV%cgs_time / GV%LittleH) * (GV%total_lum * GV%Lunit)

end subroutine

subroutine sphray_set_isothermal(flag)
  logical :: flag
  GV%FixSnapTemp=flag
end subroutine

subroutine sphray_get_isothermal(flag)
  logical :: flag
  flag=GV%FixSnapTemp
end subroutine

subroutine sphray_set_momentum_kicks(flag)
  logical :: flag
  GV%DoMomentumKicks=flag
end subroutine

subroutine sphray_get_momentum_kicks(flag)
  logical :: flag
  flag=GV%DoMomentumKicks
end subroutine

subroutine sphray_set_H_caseA(flag)
  logical :: flag
  GV%HydrogenCaseA=flag
end subroutine

subroutine sphray_get_H_caseA(flag)
  logical :: flag
  flag=GV%HydrogenCaseA
end subroutine

subroutine sphray_set_He_caseA(flag)
  logical :: flag
  GV%HeliumCaseA=flag
end subroutine

subroutine sphray_get_He_caseA(flag)
  logical :: flag
  flag=GV%HeliumCaseA
end subroutine

subroutine sphray_set_raynumber(x)
  real(r8b) :: x
  RaysperUnitTime=x
end subroutine

subroutine sphray_get_raynumber(x)
  real(r8b) :: x
  x=RaysperUnitTime
end subroutine

subroutine sphray_set_iontempsolver(N)
  integer(i4b) :: N
  GV%IonTempSolver=N
end subroutine

subroutine sphray_get_iontempsolver(N)
  integer(i4b) :: N
  N=GV%IonTempSolver
end subroutine

subroutine sphray_set_boundary(N)
  integer(i4b) :: N
  GV%BndryCond=N
end subroutine

subroutine sphray_get_boundary(N)
  integer(i4b) :: N
  N=GV%BndryCond
end subroutine

subroutine sphray_set_boxsize(x)
  real(r8b) :: x
  GV%BoxUprs=x/2
  GV%BoxLwrs=-x/2
end subroutine

subroutine sphray_get_boxsize(x)
  real(r8b) :: x
  x=maxval(GV%BoxUprs-GV%BoxLwrs)
end subroutine

subroutine sphray_set_defaultspectype(x)
  real(r8b) :: x
  default_spectral_type=x
end subroutine

subroutine sphray_get_defaultspectype(x)
  real(r8b) :: x
  x=default_spectral_type
end subroutine

function sphray_set_he_mass_frac(hemf) result(ret)
  real(r8b) :: hemf
  integer :: ret  
  if(hemf.LT.0.OR.hemf.GT.1) then
    ret=-1
    return
  endif  
#ifdef incHe
  GV%He_mf=hemf
  GV%H_mf=1.-hemf
  ret=0.
#else
  ret=-2
#endif  

end function

function sphray_get_he_mass_frac(hemf) result(ret)
  real(r8b) :: hemf
  integer :: ret
  ret=0.
#ifdef incHe
  hemf=GV%He_mf
#else
  hemf=0.
  if(GV%He_mf.NE.0) ret=-1
#endif

end function


subroutine set_default_parameters
  GV%Verbosity=3              !< [Config File] 0=silent, 1=whisper, 2=talk, 3=debug

  GV%DoTestScenario=.false.      !< [Config File] set true if performing a test problem
  GV%TestScenario='none'         !< [Config File] one of {iliev_test1, iliev_test2, iliev_test3, iliev_test4}

  GV%JustInit=.false.             !< [Config File] set true to stop after initialization
  GV%Comoving=.false.            !< [Config File] set true if values to be read are in comoving coords

  GV%IsoTemp=0                !< [Config File] if > zero all pars fixed @ IsoTemp (FixSnapTemp must be F)
  GV%FixSnapTemp=.false.         !< [Config File] if T, fix temp at snapshot values (IsoTemp must be <= 0)

  GV%EOStemp=0             !< [Config File] if non-negative, initialize EOS particles w/ T = EOStemp 
  GV%InitxHI=-1             !< [Config File] if non-negative, all xHI initialized to this value
  GV%RayDepletion=.true.        !< [Config File] remove photons from rays as they travel?

  GV%IntSeed=0123456789        !< [Config File] seed for mersenne twister

  GV%StaticFieldSimTime=1  !< [Config File] sim time for single snapshot jobs
  GV%StaticSimTimeUnit="myr"   !< [Config File] one of {codetime,myr}


  GV%InputType=1           !< [Config File] one of Gadget {1: Public 2: CosmoBH 3: OWLS/GIMIC 4: V.Bromm 5: Public HDF5}
  GV%SnapPath='./'            !< [Config File] dir where particle snapshots are
  GV%SourcePath='./'          !< [Config File] dir where source snapshots are


  GV%SpectraFile='./data/spectra/thermal1e5.cdf'        !< [Config File] file containing spectra tables
  GV%b2cdFile='./data/column_depth/latmon_b2cd_table.txt'           !< [Config File] file containing b2cd tables
  GV%AtomicRatesFile='./data/atomic_rates/atomic_rates_Hui.txt'    !< [Config File] file containnig atomic rates tables

  GV%ParFileBase='none'         !< [Config File] particle snapshot file base
  GV%SourceFileBase='none'      !< [Config File] source snapshot file base

  GV%StartSnapNum=1        !< [Config File] snapshot to start with
  GV%EndSnapNum=1          !< [Config File] snapshot to end with

  GV%ParFilesPerSnap=1     !< [Config File] files per particle snapshot
  GV%SourceFilesPerSnap=1  !< [Config File] files per source snapshot


  GV%RayScheme="raynum"           !< [Config File] one of {raynum, header}
  GV%ForcedRayNumber=30000     !< [Config File] number of rays to trace if RayScheme = raynum
  RaysperUnitTime=1.e6

  GV%RayStats=.false.            !< [Config File] T = massive output file on ray statistics in raystats.dat
  GV%BndryCond=0           !< [Config File] one of {-1:reflecting 0:vacuum 1:periodic}

  GV%RayPhotonTol=1.0d-10        !< [Config File] fractional ray depletion to stop ray
  GV%MaxRayDist=-1          !< [Config File] max ray distance in physical code units, negative=default

  GV%HydrogenCaseA=.true.       !< [Config File] T = use case A for Hydrogen Recombinations
  GV%HeliumCaseA=.true.         !< [Config File] T = use case A for Helium Recombinsations

  GV%IonTempSolver=2       !< [Config File] one of {1:euler, 2:bdf}

  GV%Tfloor=1              !< [Config File] minimum allowed temperature
  GV%Tceiling=1.d9            !< [Config File] maximum allowed temperature

  GV%xfloor=0.              !< [Config File] minimum allowed ionization fraction
  GV%xceiling=1.            !< [Config File] maximum allowed ionization fraction

  GV%NeBackground=1.e-8        !< [Config File] constant background electron number density from metals
  GV%NraysUpdateNoHits=0   !< [Config File] update all pars not hit by a ray in last NraysUpdateNoHits

  GV%H_mf=1.                !< [Config File] hydrogen mass fraction
  GV%He_mf=0.               !< [Config File] helium mass fraction

  GV%OutputDir='./'           !< [Config File] path to output directory
  GV%OutputFileBase='none'      !< [Config File] output file base

  GV%OutputType=1          !< [Config File] one of {1:Standard Binary Gadget 2:HDF5 Gadget}

  GV%OutputTiming='standard'        !< [Config File] one of {standard, forced}
  GV%NumStdOuts=-1          !< [Config File] if OutputTiming = "standard", # of outputs (maybe +1 initial)
  
  GV%DoInitialOutput=.false.     !< [Config File] produces output before any raytracing
  GV%IonFracOutRays=10000      !< [Config File] do mini output every IonFracOutRays src rays

  GV%ForcedOutFile='none'       !< [Config File] file with forced output times
  GV%ForcedUnits='myr'         !< [Config File] one of {codetime, myr, mwionfrac, vwionfrac}

  GV%PartPerCell=12         !< [Config File] minimum particles in a tree leaf

  GV%config_file='none'         !< name of the config file

  GV%BoxUprs=6.6
  GV%BoxLwrs=-6.6

  GV%DoMomentumKicks=.false.

end subroutine



function temp_from_u(u, ye, Hmf) result(T)
  use gadget_General_class, only: gadget_constants_type

  real(r8b), intent(in) :: u,ye,Hmf
  real(r8b) :: mu
  real(r8b) :: T,Tdum
  type(gadget_constants_type) :: gconst

  mu = 4.0d0 / (3.0d0 * Hmf + 1.0d0 + 4.0d0 * Hmf * ye)
  Tdum = mu * gconst%PROTONMASS / gconst%BOLTZMANN * (gconst%GAMMA - 1.0d0) * u
  T = Tdum * GV%cgs_enrg / GV%cgs_mass

end function

!> converts temperature K to internal energies / unit mass 
!=======================================================================================
function u_from_temp(T,ye,Hmf) result(U)
  use gadget_General_class, only: gadget_constants_type

  real(r8b), intent(in) :: T,ye,Hmf
  real(r8b) :: mu
  real(r8b) :: U,Udum
  type(gadget_constants_type) :: gconst

  mu = 4.0d0 / (3.0d0 * Hmf + 1.0d0 + 4.0d0 * Hmf * ye)
  Udum = gconst%BOLTZMANN * T /( (gconst%GAMMA - 1.0d0) * mu * gconst%PROTONMASS )
  U = Udum * GV%cgs_mass / GV%cgs_enrg

end function

subroutine extend_par(buf,n)
  type(particle_type), allocatable, intent (inout) :: buf(:)
  type(particle_type), allocatable :: tmpbuf(:)
  integer(i8b) :: n,m   

  m=0
  if(allocated(buf)) then
    m=min(n,size(buf))
    allocate(tmpbuf(m))
    tmpbuf(1:m)=buf(1:m)
    deallocate(buf)
  endif
  
  allocate(buf(n))
  
  
  if(m.GT.0 .and. allocated(tmpbuf)) then
    buf(1:m)=tmpbuf(1:m)
    deallocate(tmpbuf)
  endif
    
end subroutine

subroutine extend_src(buf,n)
  type(source_type), allocatable :: buf(:),tmpbuf(:)
  integer(i8b) :: n,m   

  m=0
  if(allocated(buf)) then
    m=min(n,size(buf))
    allocate(tmpbuf(m))
    tmpbuf(1:m)=buf(1:m)
    deallocate(buf)
  endif

  allocate(buf(n))
  
  if(m.GT.0) then
    buf(1:m)=tmpbuf(1:m)
    deallocate(tmpbuf)
  endif
    
end subroutine

function find_gas(id_) result(index)
  use hashMod
  type(hash_type),save ::  hash
  integer id_,index
  integer, save :: nbod=0
  
  if(.NOT.gas_searcheable) then
    nbod=size(psys%par)
    if(allocated(gid)) deallocate(gid)
    allocate(gid(nbod))
    gid(1:nbod)=psys%par(1:nbod)%id
    call initHash(nbod/2+1,nbod, gid,hash)
    gas_searcheable=.TRUE.
  endif
  
  index=find(id_,gid,hash)  

  if(index.LE.0) then
    index=-1
    return
  endif
  if(index.GT.nbod) then
    index=-2
    return
  endif
  if(gid(index).NE.id_) then
    index=-3
    return
  endif      
  
  if(psys%par(index)%mass<0) then
    index=-4
    return
  endif  
  
end function  


function find_src(id_) result(index)
  use hashMod
  type(hash_type), save :: hash
  integer id_,index
  integer, save :: nbod=0
  
  if(.NOT.src_searcheable) then
    nbod=size(psys%src)
    if(allocated(sid)) deallocate(sid)
    allocate(sid(nbod))
    sid(1:nbod)=psys%src(1:nbod)%id
    call initHash(nbod/2+1,nbod, sid,hash)
    src_searcheable=.TRUE.
  endif
  
  index=find(id_,sid,hash)  

  if(index.LE.0) then
    index=-1
    return
  endif
  if(index.GT.nbod) then
    index=-2
    return
  endif
  if(sid(index).NE.id_) then
    index=-3
    return
  endif      
   
  if(psys%src(index)%L<0) then
    index=-4
    return
  endif  

end function  


end module
