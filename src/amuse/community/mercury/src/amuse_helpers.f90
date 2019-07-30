! what to do with close encounters?
! capture collisions tbd

! internal units? K2 scaling?

module amuse_mercuryMod
  
  private
  public :: mercury_init, mercury_end,finish_init,evolve_mercury, &
    add_particle, get_particle_state,set_particle_state,remove_particle, &
    energy_angular_momentum, get_number_of_particles, &
    get_position_src, get_velocity_src, &
    set_position_src, set_velocity_src, &
    get_density_src, set_density_src, &
    get_radius_src, set_radius_src, &
    get_spin_src, set_spin_src, &
    get_mass_src, set_mass_src, &
    energy_angular_momentum_deviation, total_energy_angular_momentum, &
    set_initial_timestep_src, &
    get_initial_timestep_src, &
    mercury_time,set_central_body,get_central_body, &
    mercury_set_begin_time, mercury_get_begin_time, &
    mercury_commit_parameters,amuse_get_gravity_at_point, &
    amuse_get_potential_at_point, &
    set_algor,get_algor, &
    mercury_set_rmax, mercury_get_rmax, &
    mercury_set_cefac, mercury_get_cefac, &
    get_outputfiles,set_outputfiles

  include 'mercury.inc'

  integer algor,nbod,nbig,stat(NMAX),lmem(NMESS)
  integer opflag,ngflag,ndump,nfun
  real*8 m(NMAX),xh(3,NMAX),vh(3,NMAX),s(3,NMAX),rho(NMAX)
  real*8 rceh(NMAX),epoch(NMAX),ngf(4,NMAX),rmax,rcen,jcen(3)
  real*8 cefac,time,tstart,tstop,dtout,h0,tol,en(3),am(3)
  character*8 id(NMAX)
  character*4096 outfile(3), dumpfile(4), mem(NMESS)
  external mdt_mvs, mdt_bs1, mdt_bs2, mdt_ra15, mdt_hy
  external mco_dh2h,mco_h2dh
  external mco_b2h,mco_h2b,mco_h2mvs,mco_mvs2h,mco_iden

  real*8,parameter :: rhocgs = AU * AU * AU * K2 / MSUN

  integer :: opt(8)=(/0,0,3,3,0,0,0,0/)

! note that opflag= -2??? (=synchronizing epoch, 0=main integration)
! is not implemented; this may be convenient/necessary at some point
! calculating actual small planets with definite epochs to their data

  real*8 :: begin_time=0.0
  integer :: tot_id=0.
  integer :: iid(NMAX)
  logical :: id_searcheable=.FALSE.

  real*8 :: kinetic_energy, potential_energy, angular_momentum(3)


 contains

function mercury_init() result(ret)
  integer :: ret

  tot_id=1
  iid(1)=1
  id(1)='centre'
  algor=10
  nbod=1
  nbig=1
  opflag=0  
  ngflag=0
  ndump=500
  nfun=100
  rmax=100
  rcen=5.d-3
  jcen=(/0.,0.,0./)
  cefac=3.
  time=0.
  tstart=0.  
  tstop=0.
  dtout=1.e30
  h0=8
  tol=1.e-12
  opt=(/0,0,3,3,0,0,0,0/)
  outfile(1)="/dev/null" !"osc_coor_vel_masses.out"
  outfile(2)="/dev/null" !"close_enc.out"
  outfile(3)="/dev/null" !"info.out"
  dumpfile(1)="/dev/null"! "bigbody_data.dmp"
  dumpfile(2)="/dev/null"!"smallbody_data.dmp"
  dumpfile(3)="/dev/null"!"int_parameters.dmp"
  dumpfile(4)="/dev/null"!"restart.dmp"

  call messages()

  m(1)=1.0*K2
  jcen(1) = jcen(1) * rcen * rcen
  jcen(2) = jcen(2) * rcen * rcen * rcen * rcen
  jcen(3) = jcen(3) * rcen * rcen * rcen * rcen * rcen * rcen

  xh(1:3,1)=0.
  vh(1:3,1)=0.
  s(1:3,1)=0.d0*K2
  begin_time=0.0
  id_searcheable=.FALSE.

  ret=0

end function

function set_initial_timestep_src(init_timestep) result(ret)
  integer :: ret
  real*8 :: init_timestep

  h0=init_timestep
  ret = 0

end function set_initial_timestep_src

function get_initial_timestep_src(init_timestep) result(ret)
  integer :: ret
  real*8 :: init_timestep

  init_timestep=h0
  ret = 0

end function get_initial_timestep_src


function get_number_of_particles(np) result(ret)
 integer :: ret,np
 np=nbod-1
 ret=0
end function

function mercury_time(timeout) result(ret)
  integer :: ret
  real*8 :: timeout
  timeout=time
  ret=0
end function

function mercury_get_begin_time(system_time) result(ret)
    implicit none
    
    integer :: ret
    real*8, intent(out) :: system_time

    system_time = begin_time

    ret = 0
end function  

function mercury_set_begin_time(system_time) result(ret)
    implicit none
    
    integer :: ret
    real*8, intent(in) :: system_time
    
    begin_time = system_time 

    ret = 0
end function  

function check_file(path) result(ret)
    implicit none
    integer :: ret, err
    character*4096, intent(in) :: path 
    character*40 :: position
    if(TRIM(path) .EQ. '/dev/null') then
        position = 'asis'
    else
        position = 'append'
    end if

    open (24,file=trim(path),status='unknown',position=position, iostat = err)
    close(24)
    if (err.GT.0)  then
       ret = -1
       PRINT *, 'error in opening file with name =', TRIM(path),', code =', err 
    else
       ret = 0
    end if
end function
    
function mercury_commit_parameters() result(ret)
    implicit none
    integer :: ret
    if(time.EQ.0.0) then
        time = begin_time
    end if
    ret = check_file(outfile(1))
    ret = check_file(outfile(2))
    ret = check_file(outfile(3))
    ret = check_file(dumpfile(1))
    ret = check_file(dumpfile(2))
    ret = check_file(dumpfile(3))
    ret = check_file(dumpfile(4))
end function


function set_central_body(mass, radius, oblateness,spin) result(ret)
  integer :: ret
  real*8, optional :: mass, radius,oblateness(3),spin(3)
  if(present(mass)) then
    m(1)=mass*K2
  endif
  if(present(radius)) then
    rcen=radius  
  endif
  if(present(oblateness)) then
    jcen(1) = oblateness(1) * rcen * rcen
    jcen(2) = oblateness(2) * rcen * rcen * rcen * rcen
    jcen(3) = oblateness(3) * rcen * rcen * rcen * rcen * rcen * rcen
  endif
  if(present(spin)) then
    s(1:3,1)=spin(1:3)*K2 
  endif
  ret=0
end function

function get_central_body(mass,radius,oblateness,spin) result(ret)
  integer :: ret
  real*8, optional :: mass, radius,oblateness(3),spin(3)
  if(present(mass)) then
    mass=m(1)/K2
  endif
  if(present(radius)) then
    radius=rcen  
  endif
  if(present(oblateness)) then
    oblateness(1) = jcen(1) / (rcen * rcen)
    oblateness(2) = jcen(2) / (rcen * rcen * rcen * rcen)
    oblateness(3) = jcen(3) / (rcen * rcen * rcen * rcen * rcen * rcen)
  endif
  if(present(spin)) then
    spin(1:3)=s(1:3,1)/K2 
  endif
  ret=0
end function

function get_particle_state(id_,mass,dens,x,y,z,vx,vy,vz,sx,sy,sz,celimit) result(ret)
  integer :: ret,id_,index
  real*8, optional :: mass,dens,x,y,z,vx,vy,vz,sx,sy,sz,celimit

  index=find_particle(id_)
  if(index.LT.0) then
    ret=index
    return
  endif  
  if (present(mass)) then
     mass=m(index)/K2
  endif
!  radius=(mass*MSUN*3/(4*PI*rho(index)*AU**3))**(1./3)
  if (present(dens)) then
     dens=rho(index)/rhocgs
  endif
  if (present(x)) then
     x=xh(1,index)
  endif
  if (present(y)) then
     y=xh(2,index)
  endif
  if (present(z)) then
     z=xh(3,index)
  endif
  if (present(vx)) then
     vx=vh(1,index)
  endif
  if (present(vy)) then
     vy=vh(2,index)
  endif
  if (present(vz)) then
     vz=vh(3,index)
  endif
  if (present(sx)) then
     sx=s(1,index)/K2
  endif
  if (present(sy)) then
     sy=s(2,index)/K2
  endif
  if (present(sz)) then
     sz=s(3,index)/K2
  endif
  if (present(celimit)) then
     celimit=rceh(index)
  endif
  ret=0
end function

function set_particle_state(id_,mass,dens,x,y,z,vx,vy,vz,sx,sy,sz,celimit) result(ret)
  integer :: ret,id_,index
  real*8, optional :: mass,dens,x,y,z,vx,vy,vz,sx,sy,sz,celimit

  index=find_particle(id_)
  if(index.LT.0) then
    ret=index
    return
  endif  
  if (present(mass)) then
     m(index)=mass*K2
  endif
!  radius=(mass*MSUN*3/(4*PI*rho(index)*AU**3))**(1./3)
  if (present(dens)) then
     rho(index)=dens*rhocgs
  endif
  if (present(x)) then
     xh(1,index)=x
  endif
  if (present(y)) then
     xh(2,index)=y
  endif
  if (present(z)) then
     xh(3,index)=z
  endif
  if (present(vx)) then
     vh(1,index)=vx
  endif
  if (present(vy)) then
     vh(2,index)=vy
  endif
  if (present(vz)) then
     vh(3,index)=vz
  endif
  if (present(sx)) then
     s(1,index)=sx*K2
  endif
  if (present(sy)) then
     s(2,index)=sy*K2
  endif
  if (present(sz)) then
     s(3,index)=sz*K2
  endif
  if (present(celimit)) then
     rceh(index)=celimit
  endif
  ret=0
end function

function mercury_end() result(ret)
  integer :: ret

! cleanup by resetting everything to starting values
  ret=mercury_init()

end function

function finish_init() result(ret)
  integer :: ret

  call mxx_en (jcen,nbod,nbig,m,xh,vh,s,en(1),am(1))
  call energy_angular_momentum(1)
  en(2)=en(1)
  am(2)=am(1)
  en(3) = 0.d0
  am(3) = 0.d0

  id_searcheable=.FALSE.
  ret=0
end function

function evolve_mercury(t_end) result(ret)
  integer :: ret
  real*8 :: t_end

  id_searcheable=.FALSE.
  tstop=t_end
  if(algor.NE.1.AND. &
     algor.NE.2.AND. &
     algor.NE.10) then
    ret=-1
    return
  endif
!  print*,time,tstart,tstop,dtout,h0
  if (algor.eq.1) call mal_hcon (time,tstart,tstop,dtout,algor,h0, &
       tol,jcen,rcen,rmax,en,am,cefac,ndump,nfun,nbod,nbig,m,xh,vh,s, &
       rho,rceh,stat,id,ngf,opt,opflag,ngflag,outfile,dumpfile,mem, &
       lmem,mdt_mvs,mco_h2mvs,mco_mvs2h)
  if (algor.eq.2) call mal_hvar (time,tstart,tstop,dtout,algor,h0, &
       tol,jcen,rcen,rmax,en,am,cefac,ndump,nfun,nbod,nbig,m,xh,vh,s, &
       rho,rceh,stat,id,ngf,opt,opflag,ngflag,outfile,dumpfile,mem, &
       lmem,mdt_bs1)
  if (algor.eq.10) call mal_hcon (time,tstart,tstop,dtout,algor,h0, &
      tol,jcen,rcen,rmax,en,am,cefac,ndump,nfun,nbod,nbig,m,xh,vh,s, &
      rho,rceh,stat,id,ngf,opt,opflag,ngflag,outfile,dumpfile,mem, &
      lmem,mdt_hy,mco_h2dh,mco_dh2h)
  call energy_angular_momentum(1)
  ret=0
  if(time.NE.t_end) then
    ret=1
  endif  
  tstart=time
end function

function new_id()
  integer new_id
  tot_id=tot_id+1
  new_id=tot_id
end function

subroutine shift_particles(first,nshift)
  integer :: first,nshift
  m(first+nshift:nbod+nshift)=m(first:nbod)
  iid(first+nshift:nbod+nshift)=iid(first:nbod)
  xh(1:3,first+nshift:nbod+nshift)=xh(1:3,first:nbod)
  vh(1:3,first+nshift:nbod+nshift)=vh(1:3,first:nbod)
  s(1:3,first+nshift:nbod+nshift)=s(1:3,first:nbod)
  ngf(1:4,first+nshift:nbod+nshift)=ngf(1:4,first:nbod)
  stat(first+nshift:nbod+nshift)=stat(first:nbod)
  rho(first+nshift:nbod+nshift)=rho(first:nbod)
  rceh(first+nshift:nbod+nshift)=rceh(first:nbod)
  epoch(first+nshift:nbod+nshift)=epoch(first:nbod)
  id_searcheable=.FALSE.
end subroutine

function add_particle(id_,mass,dens,x,y,z,vx,vy,vz,sx,sy,sz,celimit) result(ret)
  integer :: ret,id_,index
  real*8 :: mass,dens,x,y,z,vx,vy,vz,sx,sy,sz,celimit

  id_searcheable=.FALSE.  
  id_=new_id()
  nbod=nbod+1
  if(nbod.GT.NMAX) then
    ret=-1
    return
  endif  
  if(mass.GT.0) then
    call shift_particles(nbig+1,1)
    nbig=nbig+1
    index=nbig
  else
    index=nbod
  endif  
  write(id(index),'(i8)') id_
  iid(index)=id_
  m(index)=mass * K2
  xh(1,index)=x
  xh(2,index)=y
  xh(3,index)=z
  vh(1,index)=vx
  vh(2,index)=vy
  vh(3,index)=vz
  stat(index)=0
  s(1,index)=sx * K2
  s(2,index)=sy * K2
  s(3,index)=sz * K2
!  rho(index)=mass*MSUN*3/(4*PI*(AU*radius)**3)
  rho(index)=dens*rhocgs
  rceh(index)=celimit
  epoch(index)=time
  ngf(1:4,index)=0.
  ret=0
end function

function set_spin_src(id, sx, sy, sz) result(ret)
  integer :: ret, id, index
  real*8 :: sx, sy, sz
  index=find_particle(id)
  if(index.LT.0) then
     ret=index
     return
  endif
  s(1,index)=sx*K2
  s(2,index)=sy*K2
  s(3,index)=sz*K2
  ret = 0
end function

function get_spin_src(id, sx, sy, sz) result(ret)
  integer :: ret, id, index
  real*8 :: sx, sy, sz
  index=find_particle(id)
  if(index.LT.0) then
     ret=index
     return
  endif
  sx = s(1,index)/K2
  sy = s(2,index)/K2
  sz = s(3,index)/K2
  ret = 0
end function

function set_density_src(id, density) result(ret)
  integer :: ret, id, index
  real*8 :: density
  index=find_particle(id)
  if(index.LT.0) then
     ret=index
     return
  endif
  rho(index) = density*rhocgs
  ret = 0
end function

function get_density_src(id, density) result(ret)
  integer :: ret, id, index
  real*8 :: density
  index=find_particle(id)
  if(index.LT.0) then
     ret=index
     return
  endif
  density = rho(index)/rhocgs
  ret = 0
end function

function set_radius_src(id, radius) result(ret)
  integer :: ret, id, index
  real*8 :: radius
  index=find_particle(id)
  if(index.LT.0) then
     ret=index
     return
  endif
  rho(index) = radius*rhocgs
  ret = 0
end function

function get_radius_src(id, radius) result(ret)
  integer :: ret, id, index
  real*8 :: radius
  index=find_particle(id)
  if(index.LT.0) then
     ret=index
     return
  endif
  radius = rho(index)/rhocgs
  ret = 0
end function

function remove_particle(id_) result(ret)
  integer id_,ret,index
  
  index=find_particle(id_)
  if(index.LT.0) then
    ret=index
    return
  endif  
  call shift_particles(index+1,-1)
  nbod=nbod-1
  if(index.LE.nbig) nbig=nbig-1
  ret=0
  id_searcheable=.FALSE.
end function

subroutine energy_angular_momentum(mode,ek,ep,l)
  integer :: mode
  real*8,optional :: ek,ep,l(3)
  if(mode.NE.0) then
    call kin_pot_ang_mom(jcen,nbod,nbig,m,xh,vh,s)
  endif
  if(present(ek)) ek=kinetic_energy/K2
  if(present(ep)) ep=potential_energy/K2
  if(present(l)) l=angular_momentum/K2
end subroutine

subroutine total_energy_angular_momentum(e_tot,am_tot)
  real*8,optional :: e_tot,am_tot
  if(present(e_tot)) e_tot=en(2)/K2
  if(present(am_tot)) am_tot=am(2)/K2  
end subroutine

subroutine energy_angular_momentum_deviation(delta_e,delta_am)
  real*8,optional :: delta_e,delta_am
  if(present(delta_e)) delta_e=en(3)/K2
  if(present(delta_am)) delta_am=am(3)/K2  
end subroutine

! get acceleration at (heliocentric) positions x1,y1,z1 smoothed with eps1
function amuse_get_gravity_at_point(eps1, x1, y1, z1, ax,ay,az, number_of_points) result(ret)
      implicit none
      integer :: ret, i,ipart
      integer, intent(in) :: number_of_points
      real*8, intent(in) :: eps1(number_of_points)
      real*8, intent(in) :: x1(number_of_points), y1(number_of_points)
      real*8, intent(in) :: z1(number_of_points)
      real*8, intent(out) :: ax(number_of_points),ay(number_of_points),az(number_of_points)
      real*8 :: mass,x,y,z,f
      real*8 :: r2

        do ipart = 1, number_of_points
            r2=x1(ipart)**2+y1(ipart)**2+z1(ipart)**2
            mass=m(1)
            if(r2.LT.rcen**2) then
              f=mass/rcen**3
            else
              f=mass/(r2*sqrt(r2))            
            endif
            ax(ipart)=-f*x1(ipart)
            ay(ipart)=-f*y1(ipart)
            az(ipart)=-f*z1(ipart)
            do i = 1, nbig
                mass=m(i+1)
                x=xh(1,i+1)
                y=xh(2,i+1)
                z=xh(3,i+1)                
                r2 = (x-x1(ipart))**2 + (y-y1(ipart))**2 + (z-z1(ipart))**2 + eps1(ipart)**2
                if (r2.GT.0) then 
                  f=mass/(r2*sqrt(r2))
                  ax(ipart) = ax(ipart) - (x1(ipart)-x)*f
                  ay(ipart) = ay(ipart) - (y1(ipart)-y)*f
                  az(ipart) = az(ipart) - (z1(ipart)-z)*f
                endif  
            enddo
        enddo
        ret = 0
end function


! get potential at (heliocentric) positions x1,y1,z1 smoothed with eps1
function amuse_get_potential_at_point(eps1, x1, y1, z1, phi, number_of_points) result(ret)
      implicit none
      integer :: ret, i,ipart
      integer, intent(in) :: number_of_points
      real*8, intent(in) :: eps1(number_of_points)
      real*8, intent(in) :: x1(number_of_points), y1(number_of_points)
      real*8, intent(in) :: z1(number_of_points)
      real*8, intent(out) :: phi(number_of_points)
      real*8 :: mass,x,y,z
      real*8 :: r2

        do ipart = 1, number_of_points
            r2=x1(ipart)**2+y1(ipart)**2+z1(ipart)**2
            mass=m(1)
            if(r2.LT.rcen**2) then
              phi(ipart)=mass/rcen*(r2/2/rcen**2-3./2)
            else
              phi(ipart)=-mass/sqrt(r2)            
            endif
            do i = 1, nbig
                mass=m(i+1)
                x=xh(1,i+1)
                y=xh(2,i+1)
                z=xh(3,i+1)                
                r2 = (x-x1(ipart))**2 + (y-y1(ipart))**2 + (z-z1(ipart))**2 + eps1(ipart)**2
                if (r2.GT.0) phi(ipart) = phi(ipart) - mass/sqrt(r2)
            enddo
        enddo
        ret = 0
end function



function find_particle(id_) result(index)
  use hashMod
  integer id_,index
  
  if(.NOT.id_searcheable) then
    call initHash(nbod/2+1,nbod, iid)
  endif
  
  index=find(id_,iid)  

  if(index.LE.0) then
    index=-1
    return
  endif
  if(index.GT.nbod) then
    index=-2
    return
  endif
  if(iid(index).NE.id_) then
    index=-3
    return
  endif      
  
end function  

function set_algor(algor_i) result(x)
  integer x,algor_i

  algor=algor_i
  
  if(algor.NE.1.AND. &
     algor.NE.2.AND. &
     algor.NE.10) then
    x=-1
  else
    x=0
  endif

end function

function get_algor(algor_o) result(x)
  integer x,algor_o
  algor_o=algor  
  x=0
end function

function mercury_set_rmax(r_max) result(ret)
    integer :: ret
    real*8 :: r_max
    rmax = r_max
    ret = 0
end function
function mercury_get_rmax(r_max) result(ret)
    integer :: ret
    real*8 :: r_max
    r_max = rmax
    ret = 0
end function

function mercury_set_cefac(cefac_n1) result(ret)
    integer :: ret
    real*8 :: cefac_n1
    cefac = cefac_n1
    ret = 0
end function
function mercury_get_cefac(cefac_n1) result(ret)
    integer :: ret
    real*8 :: cefac_n1
    cefac_n1 = cefac
    ret = 0
end function

function set_outputfiles(f1,f2,f3,f4,f5,f6,f7) result(ret)
  integer :: ret
  character*4096, optional :: f1,f2,f3,f4,f5,f6,f7

  if(present(f1)) outfile(1)=f1
  if(present(f2)) outfile(2)=f2
  if(present(f3)) outfile(3)=f3
  if(present(f4)) dumpfile(1)=f4
  if(present(f5)) dumpfile(2)=f5
  if(present(f6)) dumpfile(3)=f6
  if(present(f7)) dumpfile(4)=f7
  ret=0
end function

function get_outputfiles(f1,f2,f3,f4,f5,f6,f7) result(ret)
  integer :: ret
  character*4096, optional :: f1,f2,f3,f4,f5,f6,f7

  if(present(f1)) f1=outfile(1)
  if(present(f2)) f2=outfile(2)
  if(present(f3)) f3=outfile(3)
  if(present(f4)) f4=dumpfile(1)
  if(present(f5)) f5=dumpfile(2)
  if(present(f6)) f6=dumpfile(3)
  if(present(f7)) f7=dumpfile(4)
  ret=0
end function

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!      MAL_HCON.FOR    (ErikSoft   28 March 2001)
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! Author: John E. Chambers
!
! Does an integration using an integrator with a constant stepsize H.
! Input and output to this routine use coordinates XH, and velocities VH,
! with respect to the central body, but the integration algorithm uses
! its own internal coordinates X, and velocities V.
!
! The programme uses the transformation routines COORD and BCOORD to change
! to and from the internal coordinates, respectively.
!
!------------------------------------------------------------------------------
!
      subroutine mal_hcon (time,tstart,tstop,dtout,algor,h0,tol,jcen, &
       rcen,rmax,en,am,cefac,ndump,nfun,nbod,nbig,m,xh,vh,s,rho,rceh, &
       stat,id,ngf,opt,opflag,ngflag,outfile,dumpfile,mem,lmem,onestep, &
       coord,bcoord)

!
      use StoppingConditions
      implicit none
      
!
! Input/Output
      integer algor,nbod,nbig,stat(nbod),opt(8),opflag,ngflag
      integer lmem(NMESS),ndump,nfun
      real*8 time,tstart,tstop,dtout,h0,tol,jcen(3),rcen,rmax
      real*8 en(3),am(3),cefac,m(nbod),xh(3,nbod),vh(3,nbod)
      real*8 s(3,nbod),rho(nbod),rceh(nbod),ngf(4,nbod)
      character*8 id(nbod)
      character*4096 outfile(3),dumpfile(4),mem(NMESS)
!
! Local
      integer i,j,k,n,itmp,nclo,nhit,jhit(CMAX),iclo(CMAX),jclo(CMAX)
      integer dtflag,ejflag,stopflag,colflag,nstored
      real*8 x(3,NMAX),v(3,NMAX),xh0(3,NMAX),vh0(3,NMAX)
      real*8 rce(NMAX),rphys(NMAX),rcrit(NMAX),epoch(NMAX)
      real*8 hby2,tout,tmp0,tdump,tfun,tlog,dtdump,dtfun
      real*8 dclo(CMAX),tclo(CMAX),dhit(CMAX),thit(CMAX)
      real*8 ixvclo(6,CMAX),jxvclo(6,CMAX),a(NMAX)
      integer clock_init, clock_current
      integer count_rate, count_max
      integer is_timeout_detection_enabled
      integer stopping_index
      integer error
      double precision timeout

      external onestep,coord,bcoord
!
!------------------------------------------------------------------------------
!
! Initialize variables. DTFLAG = 0/2: first call ever/normal
      dtout  = abs(dtout)
      dtdump = abs(h0) * ndump
      dtfun  = abs(h0) * nfun
      dtflag = 0
      nstored = 0
      hby2 = 0.500001d0 * abs(h0)
!
! Calculate close-encounter limits and physical radii
      call mce_init (tstart,algor,h0,jcen,rcen,rmax,cefac,nbod,nbig, &
       m,xh,vh,s,rho,rceh,rphys,rce,rcrit,id,opt,outfile(2),1)

!
! Set up time of next output, times of previous dump, log and periodic effect
      if (opflag.eq.-1) then
        tout = tstart
      else
        n = int (abs (time-tstart) / dtout) + 1
        tout = tstart  +  dtout * sign (dble(n), tstop - tstart)
        if ((tstop-tstart)*(tout-tstop).gt.0) tout = tstop
      end if
      tdump = time
      tfun  = time
      tlog  = time
!
! Convert to internal coordinates and velocities
      call coord (time,jcen,nbod,nbig,h0,m,xh,vh,x,v,ngf,ngflag,opt)

!
      error = reset_stopping_conditions()

      error = is_stopping_condition_enabled(&
                     TIMEOUT_DETECTION, &
                     is_timeout_detection_enabled)
      error = get_stopping_condition_timeout_parameter(timeout)

      call system_clock(clock_init, count_rate, count_max)

!
!------------------------------------------------------------------------------
!
!  MAIN  LOOP  STARTS  HERE
!
      write(6,*)"just check if stuck enter loop"
 100  continue

! timeout stopping condition

      if (is_timeout_detection_enabled.gt.0) then
         call system_clock(clock_current, count_rate, count_max)
         if ((clock_current-clock_init).gt.timeout) then
            stopping_index = next_index_for_stopping_condition()
            error = set_stopping_condition_info(stopping_index, &
                 TIMEOUT_DETECTION)
         endif
      endif
! if condition met, break
      if (is_any_condition_set().gt.0) then
	 write(6,*) "condition set"
         goto 101
      endif
!
! Is it time for output ?
      if (abs(tout-time).le.hby2.and.opflag.ge.-1) then
!
! Beware: the integration may change direction at this point!!!!
        if (opflag.eq.-1.and.dtflag.ne.0) dtflag = 1
!
! Convert to heliocentric coordinates and output data for all bodies
        call bcoord (time,jcen,nbod,nbig,h0,m,x,v,xh,vh,ngf,ngflag,opt)
        call mio_out (time,jcen,rcen,rmax,nbod,nbig,m,xh,vh,s,rho, &
         stat,id,opt,opflag,algor,outfile(1))
        call mio_ce (time,tstart,rcen,rmax,nbod,nbig,m,stat,id, &
         0,iclo,jclo,opt,stopflag,tclo,dclo,ixvclo,jxvclo,mem,lmem, &
         outfile,nstored,0)
        tmp0 = tstop - tout
        tout = tout + sign( min( abs(tmp0), abs(dtout) ), tmp0 )
!
! Update the data dump files
        do j = 2, nbod
          epoch(j) = time
        end do
        call mio_dump (time,tstart,tstop,dtout,algor,h0,tol,jcen,rcen, &
         rmax,en,am,cefac,ndump,nfun,nbod,nbig,m,xh,vh,s,rho,rceh,stat, &
         id,ngf,epoch,opt,opflag,dumpfile,mem,lmem)
        tdump = time
      end if
!
! If integration has finished, convert to heliocentric coords and return
      if (abs(tstop-time).le.hby2.and.opflag.ge.0) then
        call bcoord (time,jcen,nbod,nbig,h0,m,x,v,xh,vh,ngf,ngflag,opt)
        write(6,*)"leaving subroutine"
        return
      end if
!
! Make sure the integration is heading in the right direction
 150  continue
      tmp0 = tstop - time
      if (opflag.eq.-1) tmp0 = tstart - time
      h0 = sign (h0, tmp0)
!
! Save the current heliocentric coordinates and velocities
      if (algor.eq.1) then
        call mco_iden (time,jcen,nbod,nbig,h0,m,x,v,xh0,vh0,ngf,ngflag, &
         opt)
      else
        call bcoord(time,jcen,nbod,nbig,h0,m,x,v,xh0,vh0,ngf,ngflag,opt)
      end if
      call onestep (time,tstart,h0,tol,rmax,en,am,jcen,rcen,nbod,nbig, &
       m,x,v,s,rphys,rcrit,rce,stat,id,ngf,algor,opt,dtflag,ngflag, &
       opflag,colflag,nclo,iclo,jclo,dclo,tclo,ixvclo,jxvclo,outfile, &
       mem,lmem)
      time = time + h0
!
!------------------------------------------------------------------------------
!
!  CLOSE  ENCOUNTERS
!
! If encounter minima occurred, output details and decide whether to stop
      if (nclo.gt.0.and.opflag.ge.-1) then
        write(6,*) "encounter minima occured"
        itmp = 1
        if (colflag.ne.0) itmp = 0
        call mio_ce (time,tstart,rcen,rmax,nbod,nbig,m,stat,id,nclo, &
         iclo,jclo,opt,stopflag,tclo,dclo,ixvclo,jxvclo,mem,lmem, &
         outfile,nstored,itmp)
!rude cello        if (stopflag.eq.1) return
      end if
!
!------------------------------------------------------------------------------
!
!  COLLISIONS
!
! If collisions occurred, output details and remove lost objects
      if (colflag.ne.0) then
         write(6,*) "collision condition set"
! Reindex the surviving objects
        call bcoord (time,jcen,nbod,nbig,h0,m,x,v,xh,vh,ngf,ngflag,opt)
        call mxx_elim (nbod,nbig,m,xh,vh,s,rho,rceh,rcrit,ngf,stat, &
         id,mem,lmem,outfile(3),itmp)
!
! Reset flags, and calculate new Hill radii and physical radii
        dtflag = 1
        if (opflag.ge.0) opflag = 1
        call mce_init (tstart,algor,h0,jcen,rcen,rmax,cefac,nbod,nbig, &
         m,xh,vh,s,rho,rceh,rphys,rce,rcrit,id,opt,outfile(2),1)
        call coord (time,jcen,nbod,nbig,h0,m,xh,vh,x,v,ngf,ngflag,opt)
      end if
!
!------------------------------------------------------------------------------
!
!  COLLISIONS  WITH  CENTRAL  BODY
!
!!!!!!! this can be commented out if the collisions with central
!!!!!!! body should be completely skiped

!~ ! Check for collisions with the central body
!~       if (algor.eq.1) then
!~         call mco_iden(time,jcen,nbod,nbig,h0,m,x,v,xh,vh,ngf,ngflag,opt)
!~       else
!~         call bcoord (time,jcen,nbod,nbig,h0,m,x,v,xh,vh,ngf,ngflag,opt)
!~       end if
!~       itmp = 2
!~       if (algor.eq.11.or.algor.eq.12) itmp = 3
!~       call mce_cent (time,h0,rcen,jcen,itmp,nbod,nbig,m,xh0,vh0,xh,vh, &
!~        nhit,jhit,thit,dhit,algor,ngf,ngflag)
!~ !
!~ ! If something hit the central body, restore the coords prior to this step
!~       if (nhit.gt.0) then
!~          write(6,*) "central body HIT"
!~         call mco_iden (time,jcen,nbod,nbig,h0,m,xh0,vh0,xh,vh,ngf, &
!~          ngflag,opt)
!~         time = time - h0
!~ !
!~ ! Merge the object(s) with the central body
!~         do k = 1, nhit
!~           i = 1
!~           j = jhit(k)
!~           call mce_coll (thit(k),tstart,en(3),jcen,i,j,nbod,nbig,m,xh, &
!~            vh,s,rphys,stat,id,opt,mem,lmem,outfile(3))
!~         end do
!~ !
!~ ! Remove lost objects, reset flags and recompute Hill and physical radii
!~         call mxx_elim (nbod,nbig,m,xh,vh,s,rho,rceh,rcrit,ngf,stat, &
!~          id,mem,lmem,outfile(3),itmp)
!~         if (opflag.ge.0) opflag = 1
!~         dtflag = 1
!~         call mce_init (tstart,algor,h0,jcen,rcen,rmax,cefac,nbod,nbig, &
!~          m,xh,vh,s,rho,rceh,rphys,rce,rcrit,id,opt,outfile(2),0)
!~         if (algor.eq.1) then
!~           call mco_iden (time,jcen,nbod,nbig,h0,m,xh,vh,x,v,ngf,ngflag, &
!~            opt)
!~         else
!~           call coord (time,jcen,nbod,nbig,h0,m,xh,vh,x,v,ngf,ngflag,opt)
!~         end if
!~ !
!~ ! Redo that integration time step
!~         goto 150
!~       end if
      
!!!!!!
!!!!!!      
!
!------------------------------------------------------------------------------
!
!  DATA  DUMP  AND  PROGRESS  REPORT
!
! Convert to heliocentric coords and do the data dump
      if (abs(time-tdump).ge.abs(dtdump).and.opflag.ge.-1) then
        call bcoord (time,jcen,nbod,nbig,h0,m,x,v,xh,vh,ngf,ngflag,opt)
        do j = 2, nbod
          epoch(j) = time
        end do
        call mio_ce (time,tstart,rcen,rmax,nbod,nbig,m,stat,id, &
         0,iclo,jclo,opt,stopflag,tclo,dclo,ixvclo,jxvclo,mem,lmem, &
         outfile,nstored,0)
        call mio_dump (time,tstart,tstop,dtout,algor,h0,tol,jcen,rcen, &
         rmax,en,am,cefac,ndump,nfun,nbod,nbig,m,xh,vh,s,rho,rceh,stat, &
         id,ngf,epoch,opt,opflag,dumpfile,mem,lmem)
        tdump = time
      end if
!
! Convert to heliocentric coords and write a progress report to the log file
      if (abs(time-tlog).ge.abs(dtdump).and.opflag.ge.0) then
        call bcoord (time,jcen,nbod,nbig,h0,m,x,v,xh,vh,ngf,ngflag,opt)
        call mxx_en (jcen,nbod,nbig,m,xh,vh,s,en(2),am(2))
        call mio_log (time,tstart,en,am,opt,mem,lmem)
        tlog = time
      end if
!
!------------------------------------------------------------------------------
!
!  CHECK  FOR  EJECTIONS  AND  DO  OTHER  PERIODIC  EFFECTS
!

      if (abs(time-tfun).ge.abs(dtfun).and.opflag.ge.-1) then
          if (algor.eq.1) then
            call mco_iden (time,jcen,nbod,nbig,h0,m,x,v,xh,vh,ngf,ngflag, &
           opt)
        else
           call bcoord(time,jcen,nbod,nbig,h0,m,x,v,xh,vh,ngf,ngflag,opt)
        end if
!
! Recompute close encounter limits, to allow for changes in Hill radii
        call mce_hill (nbod,m,xh,vh,rce,a)
         do j = 2, nbod
          rce(j) = rce(j) * rceh(j)
        end do
!
! Check for ejections
        itmp = 2
        if (algor.eq.11.or.algor.eq.12) itmp = 3
        write(6,*) "3.5"
!        call mxx_ejec (time,tstart,rmax,en,am,jcen,itmp,nbod,nbig,m,xh, &
!         vh,s,stat,id,opt,ejflag,outfile(3),mem,lmem)
        write(6,*)"4"
!
! Remove ejected objects, reset flags, calculate new Hill and physical radii
        if (ejflag.ne.0) then
          write(6,*) "ejected object removal"
          call mxx_elim (nbod,nbig,m,xh,vh,s,rho,rceh,rcrit,ngf,stat, &
           id,mem,lmem,outfile(3),itmp)
          if (opflag.ge.0) opflag = 1
          dtflag = 1
          call mce_init (tstart,algor,h0,jcen,rcen,rmax,cefac,nbod,nbig, &
           m,xh,vh,s,rho,rceh,rphys,rce,rcrit,id,opt,outfile(2),0)
          if (algor.eq.1) then
            call mco_iden (time,jcen,nbod,nbig,h0,m,xh,vh,x,v,ngf, &
             ngflag,opt)
          else
            call coord (time,jcen,nbod,nbig,h0,m,xh,vh,x,v,ngf,ngflag, &
             opt)
          end if
        end if
        tfun = time
      end if
!
! Go on to the next time step
      goto 100
!
!------------------------------------------------------------------------------
!


101   continue
      write(6,*)"just check if stuck exit loop"
      end subroutine


!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!      MAL_HVAR.FOR    (ErikSoft   4 March 2001)
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! Author: John E. Chambers
!
! Does an integration using a variable-timestep integration algorithm. The
! particular integrator routine is ONESTEP and the algorithm must use
! coordinates with respect to the central body.
!
! N.B. This routine is also called by the synchronisation routine mxx_sync,
! ===  in which case OPFLAG = -2. Beware when making changes involving OPFLAG.
!
!------------------------------------------------------------------------------
!
      subroutine mal_hvar (time,tstart,tstop,dtout,algor,h0,tol,jcen, &
       rcen,rmax,en,am,cefac,ndump,nfun,nbod,nbig,m,xh,vh,s,rho,rceh, &
       stat,id,ngf,opt,opflag,ngflag,outfile,dumpfile,mem,lmem,onestep)
!
      use StoppingConditions
      implicit none
!
! Input/Output
      integer algor,nbod,nbig,stat(nbod),opt(8),opflag,ngflag,ndump,nfun
      integer lmem(NMESS)
      real*8 time,tstart,tstop,dtout,h0,tol,jcen(3),rcen,rmax
      real*8 en(3),am(3),cefac,m(nbod),xh(3,nbod),vh(3,nbod)
      real*8 s(3,nbod),rho(nbod),rceh(nbod),ngf(4,nbod)
      character*8 id(nbod)
      character*4096 outfile(3),dumpfile(4),mem(NMESS)
!
! Local
      integer i,j,k,n,itmp,nhit,ihit(CMAX),jhit(CMAX),chit(CMAX)
      integer dtflag,ejflag,nowflag,stopflag,nstored,ce(NMAX)
      integer nclo,iclo(CMAX),jclo(CMAX),nce,ice(NMAX),jce(NMAX)
      real*8 tmp0,h,hdid,tout,tdump,tfun,tlog,tsmall,dtdump,dtfun
      real*8 thit(CMAX),dhit(CMAX),thit1,x0(3,NMAX),v0(3,NMAX)
      real*8 rce(NMAX),rphys(NMAX),rcrit(NMAX),a(NMAX)
      real*8 dclo(CMAX),tclo(CMAX),epoch(NMAX)
      real*8 ixvclo(6,CMAX),jxvclo(6,CMAX)

      integer clock_init, clock_current
      integer count_rate, count_max
      integer is_timeout_detection_enabled
      integer stopping_index
      integer error
      double precision timeout

      external mfo_all,onestep
!
!------------------------------------------------------------------------------
!
! Initialize variables. DTFLAG = 0 implies first ever call to ONESTEP
      dtout  = abs(dtout)
      dtdump = abs(h0) * ndump
      dtfun  = abs(h0) * nfun
      dtflag = 0
      nstored = 0
      tsmall = h0 * 1.d-8
      h = h0
      do j = 2, nbod
        ce(j) = 0.d0
      end do
!
! Calculate close-encounter limits and physical radii for massive bodies
      call mce_init (tstart,algor,h0,jcen,rcen,rmax,cefac,nbod,nbig, &
       m,xh,vh,s,rho,rceh,rphys,rce,rcrit,id,opt,outfile(2),1)
!
! Set up time of next output, times of previous dump, log and periodic effect
      if (opflag.eq.-1) then
        tout = tstart
      else
        n = int (abs (time - tstart) / dtout) + 1
        tout = tstart  +  dtout * sign (dble(n), tstop - tstart)
        if ((tstop - tstart)*(tout - tstop).gt.0) tout = tstop
      end if
      tdump = time
      tfun  = time
      tlog  = time

!
      error = reset_stopping_conditions()

      error = is_stopping_condition_enabled(&
                     TIMEOUT_DETECTION, &
                     is_timeout_detection_enabled)
      error = get_stopping_condition_timeout_parameter(timeout)

      call system_clock(clock_init, count_rate, count_max)

!

!
!------------------------------------------------------------------------------
!
!  MAIN  LOOP  STARTS  HERE
!
 100  continue

! timeout stopping condition

      if (is_timeout_detection_enabled.gt.0) then
         call system_clock(clock_current, count_rate, count_max)
         if ((clock_current-clock_init).gt.timeout) then
            stopping_index = next_index_for_stopping_condition()
            error = set_stopping_condition_info(stopping_index, &
                 TIMEOUT_DETECTION)
         endif
      endif
! if condition met, break
      if (is_any_condition_set().gt.0) then
	 write(6,*) "condition set"
         return
      endif

!
! Is it time for output ?
      if (abs(tout-time).lt.abs(tsmall).and.opflag.ge.-1) then
!
! Beware: the integration may change direction at this point!!!!
        if (opflag.eq.-1) dtflag = 0
!
! Output data for all bodies
        call mio_out (time,jcen,rcen,rmax,nbod,nbig,m,xh,vh,s,rho, &
         stat,id,opt,opflag,algor,outfile(1))
        call mio_ce (time,tstart,rcen,rmax,nbod,nbig,m,stat,id, &
         0,iclo,jclo,opt,stopflag,tclo,dclo,ixvclo,jxvclo,mem,lmem, &
         outfile,nstored,0)
        tmp0 = tstop - tout
        tout = tout + sign( min( abs(tmp0), abs(dtout) ), tmp0 )
!
! Update the data dump files
        do j = 2, nbod
          epoch(j) = time
        end do
        call mio_dump (time,tstart,tstop,dtout,algor,h,tol,jcen,rcen, &
         rmax,en,am,cefac,ndump,nfun,nbod,nbig,m,xh,vh,s,rho,rceh,stat, &
         id,ngf,epoch,opt,opflag,dumpfile,mem,lmem)
        tdump = time
      end if
!
! If integration has finished return to the main part of programme
      if (abs(tstop-time).le.abs(tsmall).and.opflag.ne.-1) return
!
! Set the timestep
      if (opflag.eq.-1) tmp0 = tstart - time
      if (opflag.eq.-2) tmp0 = tstop  - time
      if (opflag.ge.0)  tmp0 = tout   - time
      h = sign ( max( min( abs(tmp0), abs(h) ), tsmall), tmp0 )
!
! Save the current coordinates and velocities
      call mco_iden (time,jcen,nbod,nbig,h,m,xh,vh,x0,v0,ngf,ngflag,opt)
!
! Advance one timestep
      call onestep (time,h,hdid,tol,jcen,nbod,nbig,m,xh,vh,s,rphys, &
       rcrit,ngf,stat,dtflag,ngflag,opt,nce,ice,jce,mfo_all)
      time = time + hdid
!
! Check if close encounters or collisions occurred
      nclo = 0
      call mce_stat (time,h,rcen,nbod,nbig,m,x0,v0,xh,vh,rce,rphys, &
       nclo,iclo,jclo,dclo,tclo,ixvclo,jxvclo,nhit,ihit,jhit, &
       chit,dhit,thit,thit1,nowflag,stat,outfile(3),mem,lmem)
!
!------------------------------------------------------------------------------
!
!  CLOSE  ENCOUNTERS
!
! If encounter minima occurred, output details and decide whether to stop
      if (nclo.gt.0.and.opflag.ge.-1) then
        itmp = 1
        if (nhit.ne.0) itmp = 0
        call mio_ce (time,tstart,rcen,rmax,nbod,nbig,m,stat,id,nclo, &
         iclo,jclo,opt,stopflag,tclo,dclo,ixvclo,jxvclo,mem,lmem, &
         outfile,nstored,itmp)
        if (stopflag.eq.1) return
      end if
!
!------------------------------------------------------------------------------
!
!  COLLISIONS
!
! If a collision occurred, output details and resolve the collision
      if (nhit.gt.0.and.opt(2).ne.0) then
        do k = 1, nhit
          if (chit(k).eq.1) then
            i = ihit(k)
            j = jhit(k)
            call mce_coll (thit(k),tstart,en(3),jcen,i,j,nbod,nbig,m,xh, &
             vh,s,rphys,stat,id,opt,mem,lmem,outfile(3))
          end if
        end do
!
! Remove lost objects, reset flags and recompute Hill and physical radii
        call mxx_elim (nbod,nbig,m,xh,vh,s,rho,rceh,rcrit,ngf,stat, &
         id,mem,lmem,outfile(3),itmp)
        dtflag = 1
        if (opflag.ge.0) opflag = 1
        call mce_init (tstart,algor,h0,jcen,rcen,rmax,cefac,nbod,nbig, &
         m,xh,vh,s,rho,rceh,rphys,rce,rcrit,id,opt,outfile(2),1)
      end if
!
!------------------------------------------------------------------------------
!
!  COLLISIONS  WITH  CENTRAL  BODY
!
! Check for collisions
      call mce_cent (time,hdid,rcen,jcen,2,nbod,nbig,m,x0,v0,xh,vh,nhit, &
       jhit,thit,dhit,algor,ngf,ngflag)
!
! Resolve the collisions
      if (nhit.gt.0) then
        do k = 1, nhit
          i = 1
          j = jhit(k)
          call mce_coll (thit(k),tstart,en(3),jcen,i,j,nbod,nbig,m,xh, &
           vh,s,rphys,stat,id,opt,mem,lmem,outfile(3))
        end do
!
! Remove lost objects, reset flags and recompute Hill and physical radii
        call mxx_elim (nbod,nbig,m,xh,vh,s,rho,rceh,rcrit,ngf,stat, &
         id,mem,lmem,outfile(3),itmp)
        dtflag = 1
        if (opflag.ge.0) opflag = 1
        call mce_init (tstart,algor,h0,jcen,rcen,rmax,cefac,nbod,nbig, &
         m,xh,vh,s,rho,rceh,rphys,rce,rcrit,id,opt,outfile(2),0)
      end if
!
!------------------------------------------------------------------------------
!
!  DATA  DUMP  AND  PROGRESS  REPORT
!
! Do the data dump
      if (abs(time-tdump).ge.abs(dtdump).and.opflag.ge.-1) then
        do j = 2, nbod
          epoch(j) = time
        end do
        call mio_ce (time,tstart,rcen,rmax,nbod,nbig,m,stat,id, &
         0,iclo,jclo,opt,stopflag,tclo,dclo,ixvclo,jxvclo,mem,lmem, &
         outfile,nstored,0)
        call mio_dump (time,tstart,tstop,dtout,algor,h,tol,jcen,rcen, &
         rmax,en,am,cefac,ndump,nfun,nbod,nbig,m,xh,vh,s,rho,rceh,stat,&
         id,ngf,epoch,opt,opflag,dumpfile,mem,lmem)
        tdump = time
      end if
!
! Write a progress report to the log file
      if (abs(time-tlog).ge.abs(dtdump).and.opflag.ge.0) then
        call mxx_en (jcen,nbod,nbig,m,xh,vh,s,en(2),am(2))
        call mio_log (time,tstart,en,am,opt,mem,lmem)
        tlog = time
      end if
!
!------------------------------------------------------------------------------
!
!  CHECK  FOR  EJECTIONS  AND  DO  OTHER  PERIODIC  EFFECTS
!
      if (abs(time-tfun).ge.abs(dtfun).and.opflag.ge.-1) then
!
! Recompute close encounter limits, to allow for changes in Hill radii
        call mce_hill (nbod,m,xh,vh,rce,a)
        do j = 2, nbod
          rce(j) = rce(j) * rceh(j)
        end do
!
! Check for ejections
        call mxx_ejec (time,tstart,rmax,en,am,jcen,2,nbod,nbig,m,xh,vh, &
         s,stat,id,opt,ejflag,outfile(3),mem,lmem)
!
! Remove lost objects, reset flags and recompute Hill and physical radii
        if (ejflag.ne.0) then
          call mxx_elim (nbod,nbig,m,xh,vh,s,rho,rceh,rcrit,ngf,stat, &
           id,mem,lmem,outfile(3),itmp)
          dtflag = 1
          if (opflag.ge.0) opflag = 1
          call mce_init (tstart,algor,h0,jcen,rcen,rmax,cefac,nbod,nbig, &
           m,xh,vh,s,rho,rceh,rphys,rce,rcrit,id,opt,outfile(2),0)
        end if
        tfun = time
      end if
!
! Go on to the next time step
      goto 100
!
!------------------------------------------------------------------------------
!
      end subroutine











subroutine kin_pot_ang_mom(jcen,nbod,nbig,m,xh,vh,s)
      implicit none
      integer nbod,nbig
      real*8 jcen(3),m(nbod),xh(3,nbod),vh(3,nbod),s(3,nbod)
      integer j,k,iflag,itmp(8)
      real*8 x(3,NMAX),v(3,NMAX),temp,dx,dy,dz,r2,tmp,ke,pe,l(3)
      real*8 r_1,r_2,r_4,r_6,u2,u4,u6,tmp2(4,NMAX)
!
!------------------------------------------------------------------------------
!
      ke = 0.d0
      pe = 0.d0
      l(1) = 0.d0
      l(2) = 0.d0
      l(3) = 0.d0
!
! Convert to barycentric coordinates and velocities
      call mco_h2b(temp,jcen,nbod,nbig,temp,m,xh,vh,x,v,tmp2,iflag,itmp)
!
! Do the spin angular momenta first (probably the smallest terms)
      do j = 1, nbod
        l(1) = l(1) + s(1,j)
        l(2) = l(2) + s(2,j)
        l(3) = l(3) + s(3,j)
      end do
!
! Orbital angular momentum and kinetic energy terms
      do j = 1, nbod
        l(1) = l(1)  +  m(j)*(x(2,j) * v(3,j)  -  x(3,j) * v(2,j))
        l(2) = l(2)  +  m(j)*(x(3,j) * v(1,j)  -  x(1,j) * v(3,j))
        l(3) = l(3)  +  m(j)*(x(1,j) * v(2,j)  -  x(2,j) * v(1,j))
        ke = ke + m(j)*(v(1,j)*v(1,j)+v(2,j)*v(2,j)+v(3,j)*v(3,j))
      end do
!
! Potential energy terms due to pairs of bodies
      do j = 2, nbod
        tmp = 0.d0
        do k = j + 1, nbod
          dx = x(1,k) - x(1,j)
          dy = x(2,k) - x(2,j)
          dz = x(3,k) - x(3,j)
          r2 = dx*dx + dy*dy + dz*dz
          if (r2.ne.0) tmp = tmp + m(k) / sqrt(r2)
        end do
        pe = pe  -  tmp * m(j)
      end do
!
! Potential energy terms involving the central body
      do j = 2, nbod
        dx = x(1,j) - x(1,1)
        dy = x(2,j) - x(2,1)
        dz = x(3,j) - x(3,1)
        r2 = dx*dx + dy*dy + dz*dz
        if (r2.ne.0) pe = pe  -  m(1) * m(j) / sqrt(r2)
      end do
!
! Corrections for oblateness
      if (jcen(1).ne.0.or.jcen(2).ne.0.or.jcen(3).ne.0) then
        do j = 2, nbod
          r2 = xh(1,j)*xh(1,j) + xh(2,j)*xh(2,j) + xh(3,j)*xh(3,j)
          r_1 = 1.d0 / sqrt(r2)
          r_2 = r_1 * r_1
          r_4 = r_2 * r_2
          r_6 = r_4 * r_2
          u2 = xh(3,j) * xh(3,j) * r_2
          u4 = u2 * u2
          u6 = u4 * u2
          pe = pe + m(1) * m(j) * r_1 & 
             * (jcen(1) * r_2 * (1.5d0*u2 - 0.5d0) & 
             +  jcen(2) * r_4 * (4.375d0*u4 - 3.75d0*u2 + .375d0) & 
             +  jcen(3) * r_6 & 
             *(14.4375d0*u6 - 19.6875d0*u4 + 6.5625d0*u2 - .3125d0))
        end do
      end if
!
      ke= ke*0.5

      kinetic_energy=ke
      potential_energy=pe
      angular_momentum=l

end subroutine

subroutine messages()

lmem(1)=6;mem(1)="days"
lmem(2)=6;mem(2)="years"
lmem(3)=13;mem(3)="solar masses"
lmem(4)=3;mem(4)="AU"
lmem(5)=3;mem(5)="no"
lmem(6)=3;mem(6)="yes"
lmem(7)=3;mem(7)="low"
lmem(8)=6;mem(8)="medium"
lmem(9)=4;mem(9)="high"
lmem(10)=0;mem(10)=""
lmem(11)=33;mem(11)="Integration parameters"
lmem(12)=33;mem(12)="----------------------"
lmem(13)=14;mem(13)="Algorithm:"
lmem(14)=38;mem(14)="Second-order mixed-variable symplectic"
lmem(15)=24;mem(15)="Bulirsch-Stoer (general)"
lmem(16)=37;mem(16)="Bulirsch-Stoer (conservative systems)"
lmem(17)=16;mem(17)="15th-order RADAU"
lmem(18)=0;mem(18)=""
lmem(19)=0;mem(19)=""
lmem(20)=0;mem(20)=""
lmem(21)=0;mem(21)=""
lmem(22)=5;mem(22)="Test"
lmem(23)=48;mem(23)="Hybrid symplectic integrator (mixed coordinates)"
lmem(24)=44;mem(24)="Hybrid symplectic (close binary coordinates)"
lmem(25)=43;mem(25)="Hybrid symplectic (wide binary coordinates)"
lmem(26)=32;mem(26)="Integration start epoch:"
lmem(27)=32;mem(27)="Integration stop epoch:"
lmem(28)=32;mem(28)="Output interval:"
lmem(29)=32;mem(29)="Element origin:"
lmem(30)=31;mem(30)="Initial timestep:"
lmem(31)=36;mem(31)="Accuracy parameter:"
lmem(32)=36;mem(32)="Central mass:"
lmem(33)=36;mem(33)="J_2:"
lmem(34)=36;mem(34)="J_4:"
lmem(35)=36;mem(35)="J_6:"
lmem(36)=36;mem(36)="Ejection distance:"
lmem(37)=36;mem(37)="Radius of central body:"
lmem(38)=29;mem(38)="Number of Big bodies:"
lmem(39)=29;mem(39)="Number of Small bodies:"
lmem(40)=37;mem(40)="Output precision:"
lmem(41)=40;mem(41)="Includes collisions:"
lmem(42)=40;mem(42)="Includes fragmentation:"
lmem(43)=0;mem(43)=""
lmem(44)=0;mem(44)=""
lmem(45)=40;mem(45)="Includes relativity:"
lmem(46)=40;mem(46)="Includes user-defined force routine:"
lmem(47)=10;mem(47)="barycentre"
lmem(48)=12;mem(48)="central body"
lmem(49)=0;mem(49)=""
lmem(50)=0;mem(50)=""
lmem(51)=30;mem(51)="Integration details"
lmem(52)=30;mem(52)="-------------------"
lmem(53)=29;mem(53)="Initial energy:"
lmem(54)=29;mem(54)="Initial angular momentum:"
lmem(55)=65;mem(55)="Integrating massive bodies and particles up to the same epoch."
lmem(56)=34;mem(56)="Beginning the main integration."
lmem(57)=24;mem(57)="Integration complete."
lmem(58)=48;mem(58)="Fractional energy change due to integrator:"
lmem(59)=48;mem(59)="Fractional angular momentum change:"
lmem(60)=57;mem(60)="Fractional energy change due to collisions/ejections:"
lmem(61)=57;mem(61)="Fractional angular momentum change:"
lmem(62)=47;mem(62)="Continuing integration from dump files at"
lmem(63)=6;mem(63)="Time:"
lmem(64)=6;mem(64)="Date:"
lmem(65)=9;mem(65)="dE/E:"
lmem(66)=9;mem(66)="dL/L:"
lmem(67)=35;mem(67)="collided with the central body at"
lmem(68)=12;mem(68)="ejected at"
lmem(69)=12;mem(69)="was hit by"
lmem(70)=34;mem(70)="removed due to an encounter with"
lmem(71)=4;mem(71)="at"
lmem(72)=26;mem(72)="solar masses AU^2 day^-2"
lmem(73)=26;mem(73)="solar masses AU^2 day^-1"
lmem(74)=36;mem(74)="lost mass due to rotational breakup"
lmem(75)=24;mem(75)="removed due to small a"
lmem(76)=0;mem(76)=""
lmem(77)=0;mem(77)=""
lmem(78)=0;mem(78)=""
lmem(79)=0;mem(79)=""
lmem(80)=0;mem(80)=""
lmem(81)=8;mem(81)="ERROR:"
lmem(82)=49;mem(82)="Modify mercury.inc and recompile Mercury."
lmem(83)=62;mem(83)="Check the file containing initial data for Big bodies."
lmem(84)=64;mem(84)="Check the file containing initial data for Small bodies."
lmem(85)=57;mem(85)="Check the file containing integration parameters."
lmem(86)=22;mem(86)="Check files.in"
lmem(87)=27;mem(87)="This file already exists:"
lmem(88)=34;mem(88)="This file is needed to continue:"
lmem(89)=30;mem(89)="This filename is duplicated:"
lmem(90)=40;mem(90)="The total number of bodies exceeds NMAX."
lmem(91)=68;mem(91)="Data style on first line must be Cartesian, Asteroidal or Cometary"
lmem(92)=68;mem(92)="You cannot integrate non-gravitational forces using this algorithm."
lmem(93)=64;mem(93)="You cannot integrate a user-defined force using this algorithm."
lmem(94)=64;mem(94)="You cannot integrate massive Small bodies using this algorithm."
lmem(95)=66;mem(95)="Massive Small bodies must have the same epoch as the Big bodies."
lmem(96)=49;mem(96)="Check character implies input file is corrupted."
lmem(97)=62;mem(97)="Mass, density, encounter limit must be >= 0 for this object:"
lmem(98)=46;mem(98)="This integration algorithm is not available:"
lmem(99)=50;mem(99)="A problem occurred reading the parameter on line"
lmem(100)=50;mem(100)="A problem occurred reading data for this object:"
lmem(101)=56;mem(101)="A problem occured reading the epoch for the Big bodies."
lmem(102)=67;mem(102)="You cannot use non-zero J2,J4,J6 using the close-binary algorithm."
lmem(103)=34;mem(103)="Two objects both have this name:"
lmem(104)=36;mem(104)="is corrupted at line number:"
lmem(105)=42;mem(105)="Central-body radius exceeds maximum radius."
lmem(106)=68;mem(106)="Maximum/Central radius is large. Output precision will be degraded."
lmem(107)=58;mem(107)="Coordinate origin must be Central, Barycentric or Jacobi."
lmem(108)=0;mem(108)=""
lmem(109)=0;mem(109)=""
lmem(110)=0;mem(110)=""
lmem(111)=0;mem(111)=""
lmem(112)=0;mem(112)=""
lmem(113)=0;mem(113)=""
lmem(114)=0;mem(114)=""
lmem(115)=0;mem(115)=""
lmem(116)=0;mem(116)=""
lmem(117)=0;mem(117)=""
lmem(118)=0;mem(118)=""
lmem(119)=0;mem(119)=""
lmem(120)=0;mem(120)=""
lmem(121)=10;mem(121)="WARNING:"
lmem(122)=53;mem(122)="Truncating the name of this object to 8 characters:"
lmem(123)=30;mem(123)="Main integration is backwards."
lmem(124)=26;mem(124)="No Big bodies are present."
lmem(125)=28;mem(125)="No Small bodies are present."
lmem(126)=50;mem(126)="Stopping integration due to an encounter between"
lmem(127)=45;mem(127)="Throwing this object into the central body:"
lmem(128)=42;mem(128)="Setting output threshhold DA to infinity."
lmem(129)=42;mem(129)="Setting output threshhold DE to infinity."
lmem(130)=42;mem(130)="Setting output threshhold DI to infinity."
lmem(131)=43;mem(131)="Increasing the radius of the central body."
lmem(132)=56;mem(132)="Total number of current close encounters exceeds CMAX."
lmem(133)=0;mem(133)=""
lmem(134)=0;mem(134)=""
lmem(135)=0;mem(135)=""
lmem(136)=0;mem(136)=""
lmem(137)=0;mem(137)=""
lmem(138)=0;mem(138)=""
lmem(139)=0;mem(139)=""
lmem(140)=0;mem(140)=""
lmem(141)=0;mem(141)=""
lmem(142)=0;mem(142)=""
lmem(143)=0;mem(143)=""
lmem(144)=0;mem(144)=""
lmem(145)=0;mem(145)=""
lmem(146)=0;mem(146)=""
lmem(147)=0;mem(147)=""
lmem(148)=0;mem(148)=""
lmem(149)=0;mem(149)=""
lmem(150)=0;mem(150)=""
lmem(151)=67;mem(151)=")O+_05 Integration parameters (WARNING: Do not delete this line!!)"
lmem(152)=66;mem(152)=")O+_05 Big-body initial data (WARNING: Do not delete this line!!)"
lmem(153)=68;mem(153)=")O+_05 Small-body initial data (WARNING: Do not delete this line!!)"
lmem(154)=39;mem(154)=") Lines beginning with `)' are ignored."
lmem(155)=70;mem(155)=")---------------------------------------------------------------------"
lmem(156)=43;mem(156)="style (Cartesian, Asteroidal, Cometary) ="
lmem(157)=20;mem(157)="epoch (in days) ="
lmem(158)=35;mem(158)=") Important integration parameters:"
lmem(159)=48;mem(159)="algorithm (MVS, BS, BS2, RADAU, HYBRID etc) ="
lmem(160)=21;mem(160)="start time (days) ="
lmem(161)=20;mem(161)="stop time (days) ="
lmem(162)=26;mem(162)="output interval (days) ="
lmem(163)=19;mem(163)="timestep (days) ="
lmem(164)=22;mem(164)="accuracy parameter ="
lmem(165)=22;mem(165)=") Integration options:"
lmem(166)=44;mem(166)="stop integration after a close encounter ="
lmem(167)=29;mem(167)="allow collisions to occur ="
lmem(168)=37;mem(168)="include collisional fragmentation ="
lmem(169)=33;mem(169)="express time in days or years ="
lmem(170)=51;mem(170)="express time relative to integration start time ="
lmem(171)=20;mem(171)="output precision ="
lmem(172)=24;mem(172)="< Not used at present >"
lmem(173)=37;mem(173)="include relativity in integration ="
lmem(174)=30;mem(174)="include user-defined force ="
lmem(175)=52;mem(175)=") These parameters do not need to be adjusted often:"
lmem(176)=26;mem(176)="ejection distance (AU) ="
lmem(177)=31;mem(177)="radius of central body (AU) ="
lmem(178)=31;mem(178)="central mass (solar masses) ="
lmem(179)=14;mem(179)="central J2 ="
lmem(180)=14;mem(180)="central J4 ="
lmem(181)=14;mem(181)="central J6 ="
lmem(182)=24;mem(182)="< Not used at present >"
lmem(183)=24;mem(183)="< Not used at present >"
lmem(184)=45;mem(184)="Hybrid integrator changeover (Hill radii) ="
lmem(185)=42;mem(185)="number of timesteps between data dumps ="
lmem(186)=48;mem(186)="number of timesteps between periodic effects ="
lmem(187)=41;mem(187)="origin (Central, Barycentric, Jacobi) ="
lmem(188)=0;mem(188)=""
lmem(189)=0;mem(189)=""
lmem(190)=0;mem(190)=""
lmem(191)=0;mem(191)=""
lmem(192)=0;mem(192)=""
lmem(193)=0;mem(193)=""
lmem(194)=0;mem(194)=""
lmem(195)=0;mem(195)=""
lmem(196)=0;mem(196)=""
lmem(197)=0;mem(197)=""
lmem(198)=0;mem(198)=""
lmem(199)=0;mem(199)=""
lmem(200)=0;mem(200)=""

end subroutine

end module

