module chem_mod
  use krome_main !use krome (mandatory)
  use krome_user !use utility (for krome_idx_* constants and others)
  implicit none

  type particle_type
    integer :: id
    double precision :: density
    double precision :: temperature
    double precision :: ionrate
    double precision :: abundances(krome_nmols)
  end type

  type(particle_type), allocatable :: particles(:)
  
  double precision :: tcurrent  ! time unit = yr
  
  integer :: nparticle
  integer :: tot_id

  integer, parameter :: NMAX=1000000

  logical :: particles_searcheable=.FALSE.
  integer, allocatable :: pid(:)

contains

  function chem_initialize()  result(ret)
    integer :: ret
    tcurrent=0.
    nparticle=0
    tot_id=0
    if(.not.allocated(particles)) allocate(particles(NMAX))
    particles(:)%density=0.
    
    print*,"initialize"
    call krome_init()
    
    ret=0
  end function

  function chem_end() result(ret)
     integer :: ret 
     if(allocated(particles)) deallocate(particles)
     ret=0
  end function

  function chem_commit_parameters() result(ret)
    integer :: ret
    ret=0
  end function

  function chem_commit_particles() result (ret)
    integer :: ret
    integer :: n
    integer :: i
 
    particles_searcheable=.FALSE.
    nparticle=clean_particles(particles)
    ret=0
  end function

  function chem_model_time(outtime) result(ret)
    integer :: ret
    double precision :: outtime
    ret=0
    outtime=tcurrent
  end function

  function evolve_chem_model(tend) result(ret)
    integer :: ret
    double precision :: tend,dt
    integer :: i,iret
    ret=0
    dt=tend-tcurrent
    if(dt.LE.0) return
    do i=1,nparticle
      iret=evolve_1_particle(particles(i),dt)
      ret=min(iret,ret)
    enddo
    tcurrent=tend
  end function

  function evolve_1_particle(particle,dt) result(ret)
    integer :: ret
    type(particle_type) :: particle
    double precision :: tend,dt
    double precision :: n(krome_nmols),T,cr

    n=particle%abundances*particle%density
    T = particle%temperature
    cr = particle%ionrate
    call krome_set_user_crate(cr)
    call krome_set_user_Av(1.d0)
    call krome_set_user_Tdust(1.d1)
    call krome(n, T, dt)
    particle%temperature=T
    particle%abundances=n/particle%density
    ret=0
  end function

  function get_particle_abundance(id, aid, abundance) result(ret)
    integer :: ret
    integer :: id,index,aid
    double precision :: abundance
    index=find_particle(id)
    if(index.LT.0) then
      ret=index
      return
    endif
    if(aid.LT.1.OR.aid.GT.krome_nmols) then
      ret=-1
      return
    endif
    abundance=particles(index)%abundances(aid)
    ret=0
  end function

  function set_particle_abundance(id, aid, abundance) result(ret)
    integer :: ret
    integer :: id,index,aid
    double precision :: abundance
    index=find_particle(id)
    if(index.LT.0) then
      ret=index
      return
    endif
    if(aid.LT.1.OR.aid.GT.krome_nmols) then
      ret=-1
      return
    endif
    particles(index)%abundances(aid)=abundance
    ret=0
  end function

  function set_particle_state(id,density,temperature,ionrate) result(ret)
    integer :: ret
    integer :: id,index
    double precision :: density, temperature, ionrate

    index=find_particle(id)
    if(index.LT.0) then
      ret=index
      return
    endif

    particles(index)%density=density
    particles(index)%temperature=temperature
    particles(index)%ionrate=ionrate

  end function

  function get_particle_state(id,density,temperature,ionrate) result(ret)
    integer :: ret
    integer :: id,index
    double precision :: density, temperature, ionrate

    index=find_particle(id)

    if(index.LT.0) then
      ret=index
      return
    endif

    density=particles(index)%density
    temperature=particles(index)%temperature
    ionrate=particles(index)%ionrate
    ret=0

  end function

  function add_particle(id,density,temperature,ionrate) result(ret)
  integer :: ret
  integer :: i,id
  double precision :: density, temperature,ionrate
  particles_searcheable=.FALSE.
  id=new_id()  
  i=nparticle+1

  if(i.GT.NMAX) then
    ret=-1
    return
  endif
  particles(i)%id=id
  particles(i)%density=density
  particles(i)%temperature=temperature
  particles(i)%ionrate=ionrate
  particles(i)%abundances=0
  nparticle=nparticle+1
  ret=0
  end function

  function remove_particle(id) result(ret)
  integer :: ret
  integer :: i,id
  i=find_particle(id)
  if(i.LE.0) then
    ret=i
    return
  endif
  if(particles(i)%density.LT.0) then
    ret=-4
    return
  endif
  particles(i)%density=-1.
  ret=0
  end function

function clean_particles(par) result(np)
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
    do while(par(left)%density.GT.0.AND.left.LT.right)
      left=left+1
    enddo
    do while(par(right)%density.LE.0.AND.left.LT.right)
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
  if(par(left)%density.GT.0) left=left+1
  np=left-1
end function

function find_particle(id_) result(index)
  use hashMod
  type(hash_type),save ::  hash
  integer id_,index
  integer, save :: nbod=0
  
  if(.NOT.particles_searcheable) then
    nbod=nparticle
    if(allocated(pid)) deallocate(pid)
    allocate(pid(nbod))
    pid(1:nbod)=particles(1:nbod)%id
    call initHash(nbod/2+1,nbod, pid,hash)
    particles_searcheable=.TRUE.
  endif
  

  index=find(id_,pid,hash)

  if(index.LE.0) then   
    index=-1
    return
  endif
  if(index.GT.nbod) then
    index=-2
    return
  endif
  if(pid(index).NE.id_) then   
    index=-3
    return
  endif
  
end function

subroutine extend_particles(buf,n)
  type(particle_type), allocatable, intent (inout) :: buf(:)
  type(particle_type), allocatable :: tmpbuf(:)
  integer :: n,m

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

function new_id()
  integer new_id
  tot_id=tot_id+1
  new_id=tot_id
end function

end module
