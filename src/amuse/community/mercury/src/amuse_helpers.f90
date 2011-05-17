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
    mercury_time,set_central_body,get_central_body

  include 'amuse_mercury.inc'

  integer algor,nbod,nbig,stat(NMAX),lmem(NMESS)
  integer opflag,ngflag,ndump,nfun
  real*8 m(NMAX),xh(3,NMAX),vh(3,NMAX),s(3,NMAX),rho(NMAX)
  real*8 rceh(NMAX),epoch(NMAX),ngf(4,NMAX),rmax,rcen,jcen(3)
  real*8 cefac,time,tstart,tstop,dtout,h0,tol,en(3),am(3)
  character*8 id(NMAX)
  character*80 outfile(3), dumpfile(4), mem(NMESS)
  external mdt_mvs, mdt_bs1, mdt_bs2, mdt_ra15, mdt_hy
  external mco_dh2h,mco_h2dh
  external mco_b2h,mco_h2b,mco_h2mvs,mco_mvs2h,mco_iden

  real*8,parameter :: rhocgs = AU * AU * AU * K2 / MSUN

  integer :: opt(8)=(/0,1,1,2,0,1,0,0/)

! note that opflag= -2??? (=synchronizing epoch, 0=main integration)
! is not implemented; this may be convenient/necessary at some point
! calculating actual small planets with definite epochs to their data

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
  opt=(/1,1,2,3,0,0,0,0/)
  outfile(1)="osc_coor_vel_masses.out"
  outfile(2)="close_enc.out"
  outfile(3)="info.out"
  dumpfile(1)="bigbody_data.dmp"
  dumpfile(2)="smallbody_data.dmp"
  dumpfile(3)="int_parameters.dmp"
  dumpfile(4)="restart.dmp"

  call messages()

  m(1)=1.0*K2
  jcen(1) = jcen(1) * rcen * rcen
  jcen(2) = jcen(2) * rcen * rcen * rcen * rcen
  jcen(3) = jcen(3) * rcen * rcen * rcen * rcen * rcen * rcen

  xh(1:3,1)=0.
  vh(1:3,1)=0.
  s(1:3,1)=0.d0*K2

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
  if(algor.NE.10) then
    ret=-1
    return
  endif
!  print*,time,tstart,tstop,dtout,h0
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
  if(present(ek)) ek=kinetic_energy
  if(present(ep)) ep=potential_energy
  if(present(l)) l=angular_momentum
end subroutine

subroutine total_energy_angular_momentum(e_tot,am_tot)
  real*8,optional :: e_tot,am_tot
  if(present(e_tot)) e_tot=en(2)
  if(present(am_tot)) am_tot=am(2)  
end subroutine

subroutine energy_angular_momentum_deviation(delta_e,delta_am)
  real*8,optional :: delta_e,delta_am
  if(present(delta_e)) delta_e=en(3)
  if(present(delta_am)) delta_am=am(3)  
end subroutine

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
  
  if(algor.NE.10) then
    x=-1
  else
    x=0
  endif  

end function

! dummies
      subroutine mio_out (time,jcen,rcen,rmax,nbod,nbig,m,xh,vh,s,rho, &
        stat,id,opt,opflag,algor,outfile)

      implicit none
!      include 'amuse_mercury.inc'

      integer nbod, nbig, stat(nbod), opt(8), opflag, algor
      real*8 time,jcen(3),rcen,rmax,m(nbod),xh(3,nbod),vh(3,nbod)
      real*8 s(3,nbod),rho(nbod)
      character*80 outfile
      character*8 id(nbod)

      opflag = 0

      return
      end subroutine

      subroutine mio_dump (time,tstart,tstop,dtout,algor,h0,tol,jcen, &
        rcen,rmax,en,am,cefac,ndump,nfun,nbod,nbig,m,x,v,s,rho,rceh, &
        stat,id,ngf,epoch,opt,opflag,dumpfile,mem,lmem)

      implicit none
!      include 'amuse_mercury.inc'

      integer algor,nbod,nbig,stat(nbod),opt(8),opflag,ndump,nfun
      integer lmem(NMESS)
      real*8 time,tstart,tstop,dtout,h0,tol,rmax,en(3),am(3)
      real*8 jcen(3),rcen,cefac,m(nbod),x(3,nbod),v(3,nbod)
      real*8 s(3,nbod),rho(nbod),rceh(nbod),ngf(4,nbod),epoch(nbod)
      character*80 dumpfile(4),mem(NMESS)
      character*8 id(nbod)

      return
      end subroutine

      subroutine mio_log (time,tstart,en,am,opt,mem,lmem)
      implicit none
!      include 'amuse_mercury.inc'
      integer lmem(NMESS), opt(8)
      real*8 time, tstart, en(3), am(3)
      character*80 mem(NMESS)

      return
      end subroutine
!-------


!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!      MIO_CE.FOR    (ErikSoft   1 March 2001)
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! Author: John E. Chambers
!
! Writes details of close encounter minima to an output file, and decides how
! to continue the integration depending upon the close-encounter option
! chosen by the user. Close encounter details are stored until either 100
! have been accumulated, or a data dump is done, at which point the stored
! encounter details are also output.
!
! For each encounter, the routine outputs the time and distance of closest
! approach, the identities of the objects involved, and the output
! variables of the objects at this time. The output variables are:
! expressed as
!  r = the radial distance
!  theta = polar angle
!  phi = azimuthal angle
!  fv = 1 / [1 + 2(ke/be)^2], where be and ke are the object's binding and
!                             kinetic energies. (Note that 0 < fv < 1).
!  vtheta = polar angle of velocity vector
!  vphi = azimuthal angle of the velocity vector
!
!------------------------------------------------------------------------------
!
      subroutine mio_ce (time,tstart,rcen,rmax,nbod,nbig,m,stat,id, &
        nclo,iclo,jclo,opt,stopflag,tclo,dclo,ixvclo,jxvclo,mem, &
        lmem,outfile,nstored,ceflush)
!
      implicit none
!      include 'amuse_mercury.inc'
!
! Input/Output
      integer nbod,nbig,opt(8),stat(nbod),lmem(NMESS),stopflag
      integer nclo,iclo(nclo),jclo(nclo),nstored,ceflush
      real*8 time,tstart,rcen,rmax,m(nbod),tclo(nclo),dclo(nclo)
      real*8 ixvclo(6,nclo),jxvclo(6,nclo)
      character*80 outfile(3),mem(NMESS)
      character*8 id(nbod)
!
! Local
      integer k,year,month
      real*8 tmp0,t1,rfac,fr,fv,theta,phi,vtheta,vphi
      character*80 c(200)
      character*38 fstop
      character*8 mio_fl2c, mio_re2c
      character*6 tstring
!
!------------------------------------------------------------------------------
!
      save c
!
! Scaling factor (maximum possible range) for distances
      rfac = log10 (rmax / rcen)
!
! Store details of each new close-encounter minimum
      do k = 1, nclo
        nstored = nstored + 1
        c(nstored)(1:8)   = mio_fl2c(tclo(k))
        c(nstored)(9:16)  = mio_re2c(dble(iclo(k)-1),0.d0,11239423.99d0)
        c(nstored)(12:19) = mio_re2c(dble(jclo(k)-1),0.d0,11239423.99d0)
        c(nstored)(15:22) = mio_fl2c(dclo(k))
!
        call mco_x2ov (rcen,rmax,m(1),0.d0,ixvclo(1,k),ixvclo(2,k), &
          ixvclo(3,k),ixvclo(4,k),ixvclo(5,k),ixvclo(6,k),fr,theta,phi, &
          fv,vtheta,vphi)
        c(nstored)(23:30) = mio_re2c (fr    , 0.d0, rfac)
        c(nstored)(27:34) = mio_re2c (theta , 0.d0, PI)
        c(nstored)(31:38) = mio_re2c (phi   , 0.d0, TWOPI)
        c(nstored)(35:42) = mio_re2c (fv    , 0.d0, 1.d0)
        c(nstored)(39:46) = mio_re2c (vtheta, 0.d0, PI)
        c(nstored)(43:50) = mio_re2c (vphi  , 0.d0, TWOPI)
!
        call mco_x2ov (rcen,rmax,m(1),0.d0,jxvclo(1,k),jxvclo(2,k), &
          jxvclo(3,k),jxvclo(4,k),jxvclo(5,k),jxvclo(6,k),fr,theta,phi, &
          fv,vtheta,vphi)
        c(nstored)(47:54) = mio_re2c (fr    , 0.d0, rfac)
        c(nstored)(51:58) = mio_re2c (theta , 0.d0, PI)
        c(nstored)(55:62) = mio_re2c (phi   , 0.d0, TWOPI)
        c(nstored)(59:66) = mio_re2c (fv    , 0.d0, 1.d0)
        c(nstored)(63:74) = mio_re2c (vtheta, 0.d0, PI)
        c(nstored)(67:78) = mio_re2c (vphi  , 0.d0, TWOPI)
      end do
!
! If required, output the stored close encounter details
      if (nstored.ge.100.or.ceflush.eq.0) then

! no output for amuse
!  10    open (22, file=outfile(2), status='old', access='append',err=10)
!        do k = 1, nstored
!          write (22,'(a1,a2,a70)') char(12),'6b',c(k)(1:70)
!        end do
!       close (22)
 
        nstored = 0
      end if
!
! If new encounter minima have occurred, decide whether to stop integration
      stopflag = 0
      if (opt(1).eq.1.and.nclo.gt.0) then

! no output for amuse
!  20    open (23, file=outfile(3), status='old', access='append',err=20)
!! If time style is Gregorian date then...
!        tmp0 = tclo(1)
!        if (opt(3).eq.1) then
!          fstop = '(5a,/,9x,a,i10,1x,i2,1x,f4.1)'
!          call mio_jd2y (tmp0,year,month,t1)
!          write (23,fstop) mem(121)(1:lmem(121)),mem(126) &
!            (1:lmem(126)),id(iclo(1)),',',id(jclo(1)), &
!            mem(71)(1:lmem(71)),year,month,t1
!! Otherwise...
!        else
!          if (opt(3).eq.3) then
!            tstring = mem(2)
!            fstop = '(5a,/,9x,a,f14.3,a)'
!            t1 = (tmp0 - tstart) / 365.25d0
!          else
!            tstring = mem(1)
!            fstop = '(5a,/,9x,a,f14.1,a)'
!            if (opt(3).eq.0) t1 = tmp0
!            if (opt(3).eq.2) t1 = tmp0 - tstart
!          end if
!          write (23,fstop) mem(121)(1:lmem(121)),mem(126) &
!            (1:lmem(126)),id(iclo(1)),',',id(jclo(1)), &
!            mem(71)(1:lmem(71)),t1,tstring
!        end if
        stopflag = 1
!        close(23)
!

      end if
!
!------------------------------------------------------------------------------
!
      return
      end subroutine

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
      implicit none
      include "../../../../../lib/stopcond/stopcond.inc"
      include 'amuse_mercury.inc'
!
! Input/Output
      integer algor,nbod,nbig,stat(nbod),opt(8),opflag,ngflag
      integer lmem(NMESS),ndump,nfun
      real*8 time,tstart,tstop,dtout,h0,tol,jcen(3),rcen,rmax
      real*8 en(3),am(3),cefac,m(nbod),xh(3,nbod),vh(3,nbod)
      real*8 s(3,nbod),rho(nbod),rceh(nbod),ngf(4,nbod)
      character*8 id(nbod)
      character*80 outfile(3),dumpfile(4),mem(NMESS)
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
      integer is_any_condition_set
      integer is_stopping_condition_enabled
      integer is_timeout_detection_enabled
      integer get_stopping_condition_timeout_parameter 
      integer next_index_for_stopping_condition
      integer set_stopping_condition_info
      integer stopping_index
      integer reset_stopping_conditions, error
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
      if (is_any_condition_set().gt.0) goto 101

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
        itmp = 1
        if (colflag.ne.0) itmp = 0
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
! If collisions occurred, output details and remove lost objects
      if (colflag.ne.0) then
!
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
! Check for collisions with the central body
      if (algor.eq.1) then
        call mco_iden(time,jcen,nbod,nbig,h0,m,x,v,xh,vh,ngf,ngflag,opt)
      else
        call bcoord (time,jcen,nbod,nbig,h0,m,x,v,xh,vh,ngf,ngflag,opt)
      end if
      itmp = 2
      if (algor.eq.11.or.algor.eq.12) itmp = 3
      call mce_cent (time,h0,rcen,jcen,itmp,nbod,nbig,m,xh0,vh0,xh,vh, &
       nhit,jhit,thit,dhit,algor,ngf,ngflag)
!
! If something hit the central body, restore the coords prior to this step
      if (nhit.gt.0) then
        call mco_iden (time,jcen,nbod,nbig,h0,m,xh0,vh0,xh,vh,ngf, &
         ngflag,opt)
        time = time - h0
!
! Merge the object(s) with the central body
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
        if (opflag.ge.0) opflag = 1
        dtflag = 1
        call mce_init (tstart,algor,h0,jcen,rcen,rmax,cefac,nbod,nbig, &
         m,xh,vh,s,rho,rceh,rphys,rce,rcrit,id,opt,outfile(2),0)
        if (algor.eq.1) then
          call mco_iden (time,jcen,nbod,nbig,h0,m,xh,vh,x,v,ngf,ngflag, &
           opt)
        else
          call coord (time,jcen,nbod,nbig,h0,m,xh,vh,x,v,ngf,ngflag,opt)
        end if
!
! Redo that integration time step
        goto 150
      end if
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
        call mxx_ejec (time,tstart,rmax,en,am,jcen,itmp,nbod,nbig,m,xh, &
         vh,s,stat,id,opt,ejflag,outfile(3),mem,lmem)
!
! Remove ejected objects, reset flags, calculate new Hill and physical radii
        if (ejflag.ne.0) then
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


!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!      MCE_INIT.FOR    (ErikSoft   28 February 2001)
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! Author: John E. Chambers
!
! Calculates close-approach limits RCE (in AU) and physical radii RPHYS
! (in AU) for all objects, given their masses M, coordinates X, velocities
! V, densities RHO, and close-approach limits RCEH (in Hill radii).
!
! Also calculates the changeover distance RCRIT, used by the hybrid
! symplectic integrator. RCRIT is defined to be the larger of N1*HILL and
! N2*H*VMAX, where HILL is the Hill radius, H is the timestep, VMAX is the
! largest expected velocity of any body, and N1, N2 are parameters (see
! section 4.2 of Chambers 1999, Monthly Notices, vol 304, p793-799).
!
! N.B. Designed to use heliocentric coordinates, but should be adequate using
! ===  barycentric coordinates.
!
!------------------------------------------------------------------------------
!
      subroutine mce_init (tstart,algor,h,jcen,rcen,rmax,cefac,nbod, &
        nbig,m,x,v,s,rho,rceh,rphys,rce,rcrit,id,opt,outfile,rcritflag)
!
      implicit none
!      include 'mercury.inc'
!
      real*8 N2,THIRD
      parameter (N2=.4d0,THIRD=.3333333333333333d0)
!
! Input/Output
      integer nbod,nbig,algor,opt(8),rcritflag
      real*8 tstart,h,jcen(3),rcen,rmax,cefac,m(nbod),x(3,nbod)
      real*8 v(3,nbod),s(3,nbod),rho(nbod),rceh(nbod),rphys(nbod)
      real*8 rce(nbod),rcrit(nbod)
      character*8 id(nbod)
      character*80 outfile
!
! Local
      integer j
      real*8 a(NMAX),hill(NMAX),temp,amin,vmax,k_2,rhocgs,rcen_2
      character*80 header,c(NMAX)
      character*8 mio_re2c, mio_fl2c
!
!------------------------------------------------------------------------------
!
      rhocgs = AU * AU * AU * K2 / MSUN
      k_2 = 1.d0 / K2
      rcen_2 = 1.d0 / (rcen * rcen)
      amin = HUGE
!
! Calculate the Hill radii
      call mce_hill (nbod,m,x,v,hill,a)
!
! Determine the maximum close-encounter distances, and the physical radii
      temp = 2.25d0 * m(1) / PI
      do j = 2, nbod
        rce(j)   = hill(j) * rceh(j)
        rphys(j) = hill(j) / a(j) * (temp/rho(j))**THIRD
        amin = min (a(j), amin)
      end do
!
! If required, calculate the changeover distance used by hybrid algorithm
      if (rcritflag.eq.1) then
        vmax = sqrt (m(1) / amin)
        temp = N2 * h * vmax
        do j = 2, nbod
          rcrit(j) = max(hill(j) * cefac, temp)
        end do
      end if
!
! Write list of object's identities to close-encounter output file
      header(1:8)   = mio_fl2c (tstart)
      header(9:16)  = mio_re2c (dble(nbig - 1),   0.d0, 11239423.99d0)
      header(12:19) = mio_re2c (dble(nbod - nbig),0.d0, 11239423.99d0)
      header(15:22) = mio_fl2c (m(1) * k_2)
      header(23:30) = mio_fl2c (jcen(1) * rcen_2)
      header(31:38) = mio_fl2c (jcen(2) * rcen_2 * rcen_2)
      header(39:46) = mio_fl2c (jcen(3) * rcen_2 * rcen_2 * rcen_2)
      header(47:54) = mio_fl2c (rcen)
      header(55:62) = mio_fl2c (rmax)
!
      do j = 2, nbod
        c(j)(1:8) = mio_re2c (dble(j - 1), 0.d0, 11239423.99d0)
        c(j)(4:11) = id(j)
        c(j)(12:19) = mio_fl2c (m(j) * k_2)
        c(j)(20:27) = mio_fl2c (s(1,j) * k_2)
        c(j)(28:35) = mio_fl2c (s(2,j) * k_2)
        c(j)(36:43) = mio_fl2c (s(3,j) * k_2)
        c(j)(44:51) = mio_fl2c (rho(j) / rhocgs)
      end do
!
! Write compressed output to file

!  50  open (22, file=outfile, status='old', access='append', err=50)
!      write (22,'(a1,a2,i2,a62,i1)') char(12),'6a',algor,header(1:62), &
!        opt(4)
!      do j = 2, nbod
!        write (22,'(a51)') c(j)(1:51)
!      end do
!      close (22)
      
!
!------------------------------------------------------------------------------
!
      return
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

