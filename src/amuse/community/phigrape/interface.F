      FUNCTION get_index_of_first_particle(id)
      IMPLICIT NONE
#include "src/general.inc"
      INTEGER get_index_of_first_particle
      INTEGER id

      get_index_of_first_particle = 0

      END

      FUNCTION get_total_radius(id, r)
      IMPLICIT NONE
#include "src/general.inc"      
      INTEGER get_total_radius
      INTEGER id
      REAL r

      get_total_radius = 0
      END

      FUNCTION new_particle(id, mm, xx, yy, zz, vx, vy, vz, rad)
      IMPLICIT NONE
#include "src/general.inc"
      INTEGER new_particle
      INTEGER id, newid, i
      DOUBLE PRECISION mm, rad, xx, yy, zz, vx, vy, vz
    
      newid = 1
      DO i=1,NMAX
         if(myind(i).LT.0) THEN
            newid = i
            exit
         end if
      ENDDO
      N = N + 1
      id = newid
      ind(N)   = N
      myid(N)  = newid
      myind(newid)= N
c      write(*,*) "Adding", id, ind(N), myid(N), myind(newid)
      m(N)     = mm
      x(1,N)   = xx
      x(2,N)   = yy
      x(3,N)   = zz
      v(1,N)   = vx
      v(2,N)   = vy
      v(3,N)   = vz
      h(N)     = rad

      t(N)     = -1.d0
      dt(N)    = -1.d0

      new_particle = 0
      END

      FUNCTION get_total_mass(total_mass)
      IMPLICIT NONE
#include "src/general.inc"
      INTEGER get_total_mass
      INTEGER counter
      DOUBLE PRECISION total_mass

      total_mass = 0

      DO counter = 1, N, 1
          total_mass = total_mass + m(counter)
      ENDDO
      
      get_total_mass = 0
      END

      FUNCTION get_index_of_next_particle(id)
      IMPLICIT NONE
#include "src/general.inc"
      INTEGER get_index_of_next_particle
      INTEGER id

      get_index_of_next_particle = 0
      END

      FUNCTION set_state(id, mm, xx, yy, zz, vx, vy, vz, rad)
      IMPLICIT NONE
#include "src/general.inc"
   
      INTEGER          id,ip, set_state
      DOUBLE PRECISION mm,xx,yy,zz,vx,vy,vz,rad

C---  get array address
      ip = myind(id)
      IF ((ip.LE.0).OR.(ip.GT.NMAX)) THEN
         PRINT*, 'set_state: No such particle exists!'
         set_state = -1
         RETURN 
      ENDIF

      m(ip) = mm
     
      x(1,ip) = xx 
      x(2,ip) = yy
      x(3,ip) = zz
      v(1,ip) = vx
      v(2,ip) = vy
      v(3,ip) = vz
      h(ip) = rad

      set_state = 0

      END

      FUNCTION set_acceleration(id, ax, ay, az)
      IMPLICIT NONE
#include "src/general.inc"
      INTEGER set_acceleration
      INTEGER id
      DOUBLE PRECISION ax, ay, az
      
      set_acceleration = -2
      
      END

      FUNCTION get_center_of_mass_position(Cx, Cy, Cz)
      IMPLICIT NONE      
#include "src/general.inc"
      INTEGER get_center_of_mass_position
      INTEGER get_total_mass
      INTEGER counter, ret
      DOUBLE PRECISION Cx, Cy, Cz
      DOUBLE PRECISION total_mass

      ret = get_total_mass(total_mass)

      Cx=0.
      Cy=0.
      Cz=0.
      Cx = 0.0
      Cy = 0.0
      Cz = 0.0

      DO counter = 1, N, 1
          Cx = Cx +m(counter)*x(1,counter)
          Cy = Cy +m(counter)*x(2,counter)
          Cz = Cz +m(counter)*x(3,counter)
      ENDDO

      Cx = Cx/total_mass
      Cy = Cy/total_mass
      Cz = Cz/total_mass

      get_center_of_mass_position = 0
      END

      FUNCTION get_center_of_mass_velocity(Cvx, Cvy, Cvz)
      IMPLICIT NONE
#include "src/general.inc"
      INTEGER get_center_of_mass_velocity
      INTEGER get_total_mass
      INTEGER counter, ret
      DOUBLE PRECISION Cvx, Cvy, Cvz
      DOUBLE PRECISION total_mass

      ret = get_total_mass(total_mass)
      
      Cvx = 0.0
      Cvy = 0.0
      Cvz = 0.0

      DO counter = 1, N, 1
          Cvx = Cvx +m(counter)*v(1,counter)
          Cvy = Cvy +m(counter)*v(2,counter)
          Cvz = Cvz +m(counter)*v(3,counter)
      ENDDO

      Cvx = Cvx/total_mass
      Cvy = Cvy/total_mass
      Cvz = Cvz/total_mass

      get_center_of_mass_velocity = 0
      END
      
      

      FUNCTION get_gravity_at_point(eps1, x1, y1, z1, fx, fy, fz, 
     $  number_of_points)
      IMPLICIT NONE
      INTEGER, intent(IN) :: number_of_points
      DOUBLE PRECISION :: eps1(number_of_points), x1(number_of_points)
      DOUBLE PRECISION ::  y1(number_of_points), z1(number_of_points)
      DOUBLE PRECISION :: fx(number_of_points), fy(number_of_points)
      DOUBLE PRECISION :: fz(number_of_points)
      DOUBLE PRECISION p(number_of_points)
      INTEGER get_gravity_at_point
#include "src/general.inc"
         
      CALL amuse_init_grape
      CALL update_grape(1) 
      
      CALL gravity_at_point(time_cur, x1, y1, z1, fx, fy, fz, p,
     $  number_of_points)
      
      CALL amuse_close_grape
        
      get_gravity_at_point = 0
      END
      

      FUNCTION get_potential_at_point(eps1, x1, y1, z1, p,
     $ number_of_points)
      IMPLICIT NONE
      INTEGER, intent(IN) :: number_of_points
      DOUBLE PRECISION :: eps1(number_of_points), x1(number_of_points)
      DOUBLE PRECISION ::  y1(number_of_points), z1(number_of_points)
      DOUBLE PRECISION :: fx(number_of_points), fy(number_of_points)
      DOUBLE PRECISION :: fz(number_of_points)
      DOUBLE PRECISION :: p(number_of_points)
      INTEGER get_potential_at_point   
#include "src/general.inc"  
      
      
      CALL amuse_init_grape
      CALL update_grape(1) 
      
      CALL gravity_at_point(time_cur, x1, y1, z1, fx, fy, fz, p,
     $  number_of_points)
      
      CALL amuse_close_grape

      get_potential_at_point = 0
      END

      FUNCTION get_velocity(id, vxout, vyout, vzout)
      IMPLICIT NONE
#include "src/general.inc"
      INTEGER get_velocity
      INTEGER          id,ip
      DOUBLE PRECISION vxout, vyout, vzout
      ip = myind(id)
      IF ((ip.LE.0).OR.(ip.GT.NMAX)) THEN
         PRINT*, 'get_velocity: No such particle exists',
     $           ' use add_particle'
         get_velocity = -1
         GOTO 999
      ENDIF
      

      vxout = v(1,ip)
      vyout = v(2,ip)
      vzout = v(3,ip)
      get_velocity = 0
      
  999 CONTINUE
      RETURN
      END

      FUNCTION get_position(id, rx, ry, rz)
      IMPLICIT NONE
#include "src/general.inc"
      INTEGER id, ip
      INTEGER get_position
      DOUBLE PRECISION rx, ry, rz
      ip = myind(id)
      IF ((ip.LE.0).OR.(ip.GT.NMAX)) THEN
         PRINT*, 'set_particle: No such particle exists',
     $           ' use add_particle'
         get_position = -1
         GOTO 999
      ENDIF

      rx = x(1,ip)
      ry = x(2,ip)
      rz = x(3,ip)
      get_position = 0
  999 CONTINUE
      RETURN
      END

      FUNCTION set_velocity(id, vx, vy, vz)
      IMPLICIT NONE 
#include "src/general.inc"

C---  local vars
      INTEGER          id,ip, set_velocity
      DOUBLE PRECISION vx,vy,vz


C---  get array address
      ip = myind(id)
      IF ((ip.LE.0).OR.(ip.GT.NMAX)) THEN
         PRINT*, 'set_particle: No such particle exists',
     $           ' use add_particle'
         set_velocity = -1
         GOTO 999
      ENDIF

C---  setting particle
      v(1,ip) = vx
      v(2,ip) = vy
      v(3,ip) = vz

      set_velocity = 0
 999  CONTINUE
      RETURN
      END


      FUNCTION set_position(id, xx, yy, zz)
      IMPLICIT NONE 
#include "src/general.inc"

C---  local vars
      INTEGER          id,ip, set_position
      DOUBLE PRECISION xx,yy,zz


C---  get array address
      ip = myind(id)
      IF ((ip.LE.0).OR.(ip.GT.NMAX)) THEN
         PRINT*, 'set_particle: No such particle exists',
     $           ' use add_particle'
         set_position = -1
         GOTO 999
      ENDIF

C---  setting particle
      x(1,ip) = xx
      x(2,ip) = yy
      x(3,ip) = zz

      set_position = 0
 999  CONTINUE
      RETURN
      END

      FUNCTION get_acceleration()
      IMPLICIT NONE
#include "src/general.inc"
      INTEGER get_acceleration      
      get_acceleration = -2
      END
            
      FUNCTION add_particle(id,mm,rad,xx,yy,zz,vx,vy,vz)
      IMPLICIT NONE 
#include "src/general.inc"
C---  local vars
      INTEGER          add_particle,id, newid, i
      DOUBLE PRECISION mm,xx,yy,zz,vx,vy,vz,rad
      
      !IF (id.GT.NMAX) THEN
      !   PRINT*, 'id is larger than NMAX which is not allowd!'
      !   add_particle = -1
      !   RETURN 
      !ENDIF
      !IF (myind(id).GE.0) THEN
      !   PRINT*, 'a particle with this id exists already!'
      !   add_particle = -1
      !   RETURN 
      !ENDIF
      newid = 1
      DO i=1,NMAX
         if(myind(i).LT.0) THEN
            newid = i
            exit
         end if
      ENDDO
      
C---  add a particle
      N = N + 1

      ind(N)   = newid
      myid(N)  = newid
      myind(newid)= N
      m(N)     = mm
      x(1,N)   = xx
      x(2,N)   = yy
      x(3,N)   = zz
      v(1,N)   = vx
      v(2,N)   = vy
      v(3,N)   = vz
      h(N)     = rad

      t(N)     = -1.d0
      dt(N)    = -1.d0

c      PRINT*, 'adding particle',m(N),h(N)
c      PRINT*, 'xyz  =',xx,yy,zz
c      PRINT*, 'vxyz =',vx,vy,vz

C---  done
      add_particle = newid

      RETURN
      END

      FUNCTION get_number_of_particles(number_of_particles)
      IMPLICIT NONE
#include "src/general.inc"

      INTEGER get_number_of_particles
      INTEGER number_of_particles

      number_of_particles = N
      get_number_of_particles = 0

      END

      FUNCTION get_eps2(ref_eps2)
      IMPLICIT NONE
#include "src/general.inc"

      INTEGER get_eps2
      DOUBLE PRECISION ref_eps2
      
      ref_eps2 = eps * eps
      get_eps2 = 0
    
      RETURN
      END


      FUNCTION get_dynamical_time_scale()
      IMPLICIT NONE
#include "src/general.inc"

      DOUBLE PRECISION get_dynamical_time_scale

      IF (inrg.EQ.1) THEN
         get_dynamical_time_scale =  (-0.5d0*MASS*MASS/E_pot)
     $                              / sqrt(2.d0*E_kin/MASS)
           
      ELSE
c         PRINT*, 'Warning: system not initialised yet'
         get_dynamical_time_scale = -1
      ENDIF

      RETURN
      END

      FUNCTION set_dt_max(new_dt_max)
      IMPLICIT NONE
#include "src/general.inc"
      
      INTEGER set_dt_max
      DOUBLE PRECISION new_dt_max
      
      dt_max = new_dt_max

      set_dt_max = 0
      END


      FUNCTION set_eps2(e2)
      IMPLICIT NONE
#include "src/general.inc"

      INTEGER set_eps2
      DOUBLE PRECISION e2

      eps     = SQRT(e2)

      set_eps2 = 0
      RETURN
      END

      FUNCTION set_bh_eps(e2)
      IMPLICIT NONE
#include "src/general.inc"

      INTEGER set_bh_eps
      DOUBLE PRECISION e2

      bh_eps     = SQRT(e2)

      set_bh_eps = 0
      RETURN
      END


      FUNCTION set_nbh(nbbh)
      IMPLICIT NONE
#include "src/general.inc"

      INTEGER nbbh, set_nbh

      nbh     = nbbh

      set_nbh = 0
      RETURN
      END
 
      FUNCTION set_eta(es,e)
      IMPLICIT NONE
#include "src/general.inc"

      INTEGER set_eta
      DOUBLE PRECISION es,e

      eta_s   = es
      eta     = e

      set_eta = 0
      END
      
      FUNCTION set_eta_s(value)
      IMPLICIT NONE
#include "src/general.inc"

      INTEGER set_eta_s
      DOUBLE PRECISION value

      eta_s   = value

      set_eta_s = 0
      END
      
      FUNCTION set_eta1(value)
      IMPLICIT NONE
#include "src/general.inc"

      INTEGER set_eta1
      DOUBLE PRECISION value

      eta   = value

      set_eta1 = 0
      RETURN
      END

      FUNCTION get_eta_s(new_eta_s)
      IMPLICIT NONE
#include "src/general.inc"

      INTEGER get_eta_s
      DOUBLE PRECISION new_eta_s

      new_eta_s = eta_s
      get_eta_s   = 0
      END
      
      FUNCTION get_eta(eta_out)
      IMPLICIT NONE
#include "src/general.inc"

      INTEGER get_eta
      DOUBLE PRECISION eta_out

      eta_out = eta
      get_eta = 0
      RETURN
      END
      
      FUNCTION get_initialize_once(result_value)
      IMPLICIT NONE
#include "src/general.inc"

      INTEGER get_initialize_once, result_value
      
      result_value = g_init_once
      get_initialize_once = 0
      RETURN
      END
      
      FUNCTION set_initialize_once(input_value)
      IMPLICIT NONE
#include "src/general.inc"

      INTEGER set_initialize_once, input_value
      
      IF(g_init_once.EQ.1) THEN
        IF(input_value.EQ.0) THEN
           IF(g_initialized.EQ.1) THEN
                CALL g6_close(clusterid) 
                g_initialized = 0
           ENDIF
        ENDIF
      ENDIF
      
      g_init_once = input_value
      
      set_initialize_once = 0
      
      RETURN
      END
      
      
      SUBROUTINE amuse_init_grape
      IMPLICIT NONE
#include "src/general.inc"

      IF(g_init_once.EQ.1) THEN
        IF(g_initialized.EQ.0) THEN
            CALL initgrape
            g_initialized = 1
        ENDIF
      ELSE
        CALL initgrape
        g_initialized = 1
      ENDIF
      
      END
      
      SUBROUTINE amuse_close_grape
      IMPLICIT NONE
#include "src/general.inc"

      IF(g_init_once.NE.1) THEN
        CALL g6_close(clusterid) 
        g_initialized = 0
      ENDIF
      
      END
C====================================================================
C===
C===  setup_module: function to setup the neccessary things for 
C===                running the module, here
C===
C===                * initialize MPI
C===                * initialize GRAPE
C===                * initialize myid array
C====================================================================
      FUNCTION initialize_code()
      IMPLICIT NONE
#include "src/general.inc"
#include "stopcond.inc"
      INTEGER i,initialize_code
      INTEGER error
      INTEGER mpi_setup_stopping_conditions
      INTEGER set_support_for_condition
    
      error = mpi_setup_stopping_conditions()
      
C     AMUSE STOPPING CONDITION SUPPORT
      error = set_support_for_condition(COLLISION_DETECTION)
      error = set_support_for_condition(NUMBER_OF_STEPS_DETECTION)
      error = set_support_for_condition(TIMEOUT_DETECTION)

      g_init_once = 0
      g_initialized = 0
      
      eta_s   = 0.01
      eta     = 0.02
      time_cur = 0.0
      time_begin = 0.0
C--------------------------------------------------------------------
C---  initialise MPI 
C--------------------------------------------------------------------
      CALL initmpi
C          -------

C--------------------------------------------------------------------
C---  set array to find particles by id to -1
C--------------------------------------------------------------------
      DO i=1,NMAX
         myind(i) = -1
      ENDDO

      initialize_code=0
      RETURN
      END

      FUNCTION commit_particles()
      IMPLICIT NONE
#include "src/general.inc"
      INTEGER is_any_condition_set
      INTEGER commit_particles
      DOUBLE PRECISION time
      INTEGER i

      
C (re)start grape (fip)
C  update local particles on grape
      CALL amuse_init_grape
      CALL update_grape(1) 

      
C---  set some vars
      inrg       = 0

C---  more vars to be set at creating the object
!      dt_max     = .03125d0
      dt_max     = 1.0


C---  more parameters to set
      eps2   = eps*eps
      dt_min = 2.d0**ndtmin   ! set ndtmin in paras.inc


C---  set individual particle times and time steps
      DO i=1,N
         t(i)   = time_cur
         dt(i)  = dt_min
      ENDDO

      n_loc  = N/n_proc         ! set the local number of particles

C---------------------------------------------------------------
C---  distribute body data
C---------------------------------------------------------------
      CALL sendbodydata
C          ------------
      
C--------------------------------------------------------------------
C---  send all particles to GRAPE
C--------------------------------------------------------------------
      CALL sendbodies2grape
C          ----------------

C--------------------------------------------------------------------
C---  select active particles from local particles
C--------------------------------------------------------------------
      min_t = time_cur + dt_min    ! initially all particles have 
      CALL selectactive            ! same time step dt_min
C          ------------

C--------------------------------------------------------------------
C---  predictor (for local active particles)
C--------------------------------------------------------------------
      CALL predictor(1)
C          ---------         

C--------------------------------------------------------------------
C---  gather active particles on all PE's
C--------------------------------------------------------------------
      CALL allgather_ap
C          ------------

C--------------------------------------------------------------------
C---  get local forces on all PE's
C---  (as this is the first call set ifirst to 1 in this call)
C--------------------------------------------------------------------
      CALL gravity(time_cur,1)
C     IF (is_any_condition_set().GT.0) THEN
C         print*, "stopping condition met during commit particles"
C     ENDIF
C          -------

C--------------------------------------------------------------------
C---  sum up the partial forces
C--------------------------------------------------------------------
      CALL sumforces(1)    ! call with 1 bc all forces etc.  need initialisation
C          ---------

C--------------------------------------------------------------------
C---  get individual time steps
C--------------------------------------------------------------------
      CALL timestep(0)
C          --------

C--------------------------------------------------------------------
C---  update GRAPE
C--------------------------------------------------------------------
      CALL update_grape(1)
C          ------------

C--------------------------------------------------------------------
C---  get energy, center of mass and momentum
C--------------------------------------------------------------------
      CALL gatherParts
      CALL energy
C          ------
      E_totini = E_tot          ! set initial value for total energy

      commit_particles = 0

C     close grape before exiting (fip)
      CALL amuse_close_grape

      RETURN
      END

C====================================================================
C===  
C===  
C===  assumes all particles are at the same time time_cur
C====================================================================
      FUNCTION recommit_particles()
      IMPLICIT NONE
#include "src/general.inc"
      
      INTEGER recommit_particles
      INTEGER i,k

C (re)start grape (fip)
C  update local particles on grape
      CALL amuse_init_grape
      CALL update_grape(1) 

C---  set individual particle times and time steps
      DO i=1,N
         t(i)  = time_cur      ! new particles are added at current time
         dt(i) = dt_min        ! with minimum timestep (will change)
         DO k=1,3               ! set acc and jerk to zero
            a(k,i)    = 0.d0
            adot(k,i) = 0.d0
         ENDDO
      ENDDO

      n_loc  = N/n_proc         ! set the local number of particles
  
C--------------------------------------------------------------------
C---  select new particles as active particles
C--------------------------------------------------------------------
      min_t = time_cur + dt_min    ! initially all particles have 
      CALL selectactive            ! same time step dt_min
C          ------------

C--------------------------------------------------------------------
C---  send new particle (marked as active) to  GRAPE
C--------------------------------------------------------------------
      CALL update_grape(0)
C          ------------

C--------------------------------------------------------------------
C---  predictor (for local active particles)
C--------------------------------------------------------------------
      CALL predictor(1)
C          ---------         

C--------------------------------------------------------------------
C---  gather active particles on all PE's
C--------------------------------------------------------------------
      CALL allgather_ap
C          ------------

C--------------------------------------------------------------------
C---  get local forces on all PE's
C---  (as this is the first call set ifirst to 1 in this call)
C--------------------------------------------------------------------
      CALL gravity(time_cur,1)
C          -------

C--------------------------------------------------------------------
C---  sum up the partial forces
C--------------------------------------------------------------------
      CALL sumforces(1)    ! call with 1 bc all forces etc.  need initialisation
C          ---------

C--------------------------------------------------------------------
C---  get individual time steps
C--------------------------------------------------------------------
      CALL timestep(1)
C          --------

C--------------------------------------------------------------------
C---  update GRAPE
C--------------------------------------------------------------------
      CALL update_grape(0)
C          ------------

C--------------------------------------------------------------------
C---  get energy, center of mass and momentum
C--------------------------------------------------------------------
      CALL gatherParts
      CALL energy
C          ------

      recommit_particles = 0

C     close grape before exiting (fip)
      CALL amuse_close_grape

      RETURN
      END


      FUNCTION cleanup_code()      
      IMPLICIT NONE
#include "src/general.inc"

      INTEGER cleanup_code, err, delete_all_particles
#ifndef NOMPI
      INTEGER ierr
#endif

      N=0
      err = delete_all_particles()

C---  just release GRAPE
      IF (g_init_once.EQ.1) THEN
        IF (g_initialized.EQ.1) THEN
            CALL g6_close(clusterid)
            g_initialized = 0
        ENDIF
      ENDIF
C          --------

      cleanup_code = 0
      RETURN
      END
      
      

      FUNCTION delete_all_particles()
      IMPLICIT NONE 
#include "src/general.inc"

C---  local vars
      INTEGER delete_all_particles,i,k


      delete_all_particles = -1

C---  shift particles after ip to free spot
      DO i=1,N-1
         m(i) = 0.0
         myind(i)= -1
         DO k=1,3
             x(k,i)    = 0.d0
             v(k,i)    = 0.d0
             a(k,i)    = 0.d0
             adot(k,i) = 0.d0
         ENDDO
     
         pot(i)       = 0.d0
         t(i)         = 0.d0
         dt(i)        = dt_min
         ind(i)       = 0
         myid(i)      = 0
      ENDDO
      N        = 0
      n_loc    = 0
      nbh = 0   ! removed particle was BH

      delete_all_particles = 0

      END
      
      
      


      FUNCTION set_particle(id,mm,xx,yy,zz,vx,vy,vz,rad)
      IMPLICIT NONE 
#include "src/general.inc"

C---  local vars
      INTEGER          set_particle,id,ip
      DOUBLE PRECISION mm,xx,yy,zz,vx,vy,vz,rad


C---  get array address
      ip = myind(id)
      IF ((ip.LE.0).OR.(ip.GT.NMAX)) THEN
         PRINT*, 'set_particle: No such particle exists',
     $           ' use add_particle'
         set_particle = -1
         GOTO 999
      ENDIF

C---  setting particle
      m(ip)   = mm
      x(1,ip) = xx
      x(2,ip) = yy
      x(3,ip) = zz
      v(1,ip) = vx
      v(2,ip) = vy
      v(3,ip) = vz
      h(ip)   = rad

      set_particle = 0
 999  CONTINUE
      RETURN
      END


      FUNCTION set_mass(id,mm)
      IMPLICIT NONE 
#include "src/general.inc"

C---  local vars
      INTEGER          set_mass,id,ip
      DOUBLE PRECISION mm


C---  get array address
      ip = myind(id)
      IF ((ip.LE.0).OR.(ip.GT.NMAX)) THEN
         PRINT*, 'set_mass: No such particle exists!'
         set_mass = -1
         GOTO 999
      ENDIF

C---  setting particle
      m(ip)   = mm

      set_mass = 0
 999  CONTINUE
      RETURN
      END

      FUNCTION get_mass(id, value)
      IMPLICIT NONE 
#include "src/general.inc"

C---  local vars
      INTEGER          id,ip, get_mass
      DOUBLE PRECISION value


C---  get array address
      ip = myind(id)
      
      IF ((ip.LE.0).OR.(ip.GT.NMAX)) THEN
         PRINT*, 'get_mass: No such particle exists!'
         get_mass = -1
         value = 0.d0
      GOTO 999
      ENDIF
      value = m(ip)
      get_mass = 0
 999  CONTINUE
      RETURN
      END

      FUNCTION get_potential(id, value)
      IMPLICIT NONE 
#include "src/general.inc"

C---  local vars
      INTEGER          id,ip, get_potential
      DOUBLE PRECISION value


C---  get array address
      ip = myind(id)
      IF ((ip.LE.0).OR.(ip.GT.NMAX)) THEN
         PRINT*, 'get_potential: No such particle exists!'
         get_potential = -1
         RETURN
      ENDIF

      value = pot(ip)
      
      get_potential = 0
      RETURN
      END

      FUNCTION set_radius(id,rad)
      IMPLICIT NONE 
#include "src/general.inc"

C---  local vars
      INTEGER          set_radius,id,ip
      DOUBLE PRECISION rad


C---  get array address
      ip = myind(id)
      IF ((ip.LE.0).OR.(ip.GT.NMAX)) THEN
         PRINT*, 'set_mass: No such particle exists!'
         set_radius = -1
         GOTO 999
      ENDIF

C---  setting particle
      h(ip) = rad

      set_radius = 0
 999  CONTINUE
      RETURN
      END

      FUNCTION get_radius(id, value)
      IMPLICIT NONE 
#include "src/general.inc"

C---  local vars
      INTEGER          id,ip
      INTEGER          get_radius
      DOUBLE PRECISION value


C---  get array address
      ip = myind(id)
      IF ((ip.LE.0).OR.(ip.GT.NMAX)) THEN
         PRINT*, 'set_mass: No such particle exists!'
         get_radius = -1.d0
         GOTO 999
      ENDIF
      get_radius = 0
      value = h(ip)
 999  CONTINUE
      RETURN
      END


      FUNCTION set_pos(id,xx,yy,zz)
      IMPLICIT NONE 
#include "src/general.inc"

C---  local vars
      INTEGER          set_pos,id,ip
      DOUBLE PRECISION xx,yy,zz


C---  get array address
      ip = myind(id)
      IF ((ip.LE.0).OR.(ip.GT.NMAX)) THEN
         PRINT*, 'set_pos: No such particle exists!'
         set_pos = -1
         GOTO 999
      ENDIF

C---  setting particle
      x(1,ip) = xx
      x(2,ip) = yy
      x(3,ip) = zz

      set_pos = 0
 999  CONTINUE
      RETURN
      END


      FUNCTION set_vel(id,vx,vy,vz)
      IMPLICIT NONE 
#include "src/general.inc"

C---  local vars
      INTEGER          set_vel,id,ip
      DOUBLE PRECISION vx,vy,vz


C---  get array address
      ip = myind(id)
      IF ((ip.LE.0).OR.(ip.GT.NMAX)) THEN
         PRINT*, 'set_vel: No such particle exists!'
         set_vel = -1
         GOTO 999
      ENDIF

C---  setting particle
      v(1,ip) = vx
      v(2,ip) = vy
      v(3,ip) = vz

      set_vel = 0
 999  CONTINUE
      RETURN
      END


      FUNCTION get_time(system_time)
      IMPLICIT NONE 
#include "src/general.inc"

C---  local vars
      INTEGER get_time
      DOUBLE PRECISION system_time

      system_time = time_cur

      get_time = 0

      END

      FUNCTION get_begin_time(system_time)
      IMPLICIT NONE 
#include "src/general.inc"

C---  local vars
      INTEGER get_begin_time
      DOUBLE PRECISION system_time

      system_time = time_begin

      get_begin_time = 0

      END

      FUNCTION set_begin_time(system_time)
      IMPLICIT NONE 
#include "src/general.inc"

C---  local vars
      INTEGER set_begin_time
      DOUBLE PRECISION system_time

      time_begin = system_time

      set_begin_time = 0

      END
      FUNCTION get_time_step(time_step_arg)
      IMPLICIT NONE 
#include "src/general.inc"

C---  local vars
      INTEGER  get_time_step
      DOUBLE PRECISION time_step_arg

C---  get minimum time step from particles
      CALL get_min_t
C          ---------
      time_step_arg = min_t - time_cur
      
      get_time_step = 0
      
      END

      FUNCTION get_kinetic_energy(kinetic_energy)
      IMPLICIT NONE 
#include "src/general.inc"

C---  local vars
      INTEGER get_kinetic_energy
      DOUBLE PRECISION kinetic_energy

      IF (inrg.EQ.1) THEN
          kinetic_energy = E_kin
          get_kinetic_energy = 0
      ELSE
          kinetic_energy = 0.0
          get_kinetic_energy = 0
      ENDIF
        
      END

      FUNCTION get_potential_energy(potential_energy)
      IMPLICIT NONE 
#include "src/general.inc"

C---  local vars
      INTEGER get_potential_energy
      DOUBLE PRECISION potential_energy

      IF (inrg.EQ.1) THEN
          potential_energy = E_pot
          get_potential_energy = 0
      ELSE
          potential_energy = 0.0
          get_potential_energy = 0
      ENDIF

      END

      FUNCTION get_energy_error()
      IMPLICIT NONE 
#include "src/general.inc"

C---  local vars
      DOUBLE PRECISION get_energy_error

      IF (inrg.EQ.1) THEN
         get_energy_error = (E_tot-E_totini)
      ELSE
         get_energy_error = 0.d0
      ENDIF

      RETURN
      END

      FUNCTION get_state(id,mm,xx,yy,zz,vx,vy,vz,rad)
      IMPLICIT NONE 
#include "src/general.inc"
      
      INTEGER          id,ip, get_state
      DOUBLE PRECISION mm,xx,yy,zz,vx,vy,vz,rad

C---  get array address
      ip = myind(id)
      IF ((ip.LE.0).OR.(ip.GT.NMAX)) THEN
         PRINT*, 'get_state: No such particle exists!'
         get_state = -1
         GOTO 999
      ENDIF

      mm = m(ip)
      xx = x(1,ip)
      yy = x(2,ip)
      zz = x(3,ip)
      vx = v(1,ip)
      vy = v(2,ip)
      vz = v(3,ip)
      rad= h(ip)
      get_state = 0
 999  CONTINUE
      RETURN
      END

      FUNCTION delete_particle(id)
      IMPLICIT NONE 
#include "src/general.inc"

C---  local vars
      INTEGER delete_particle,id,ip,i,k,init
      DOUBLE PRECISION a2by18(NDIM),a1by6(NDIM),aby2(NDIM)


      delete_particle = -1

C---  get array address
      ip = myind(id)
      IF ((ip.LE.0).OR.(ip.GT.NMAX)) THEN
         PRINT*, 'remove_particle: No such particle exists!'
         delete_particle = -1
         RETURN
      ENDIF

c      CALL amuse_init_grape
c      CALL update_grape(1) 

cC---  synchronize system
c      CALL selectsync(time_cur)
cC          ----------
c      CALL sync
cC          ----
cC     close grape before continuing
c      CALL amuse_close_grape

C---  shift particles after ip to free spot
c      write(*,*) "Removing", id, ind(id), myid(ip), myind(id)
      DO i=ip,N-1
         m(i) = m(i+1)
         DO k=1,3
            x(k,i)    = x(k,i+1)
            v(k,i)    = v(k,i+1)
            a(k,i)    = a(k,i+1)
            adot(k,i) = adot(k,i+1)

C---  for shifting particles on GRAPE
c            a2by18(k) = 0.d0 
c            a1by6(k)  = adot(k,i) * over6 
c            aby2(k)   = a(k,i)    * over2

         ENDDO
         pot(i)       = pot(i+1)
         t(i)         = t(i+1)
         dt(i)        = dt(i+1)
         ind(i)       = i
         myid(i)      = myid(i+1)
         myind(myid(i))= i


      ENDDO
      N        = N - 1
      n_loc    = N/n_proc
      IF (ip.LE.nbh) nbh = nbh - 1   ! removed particle was BH
      myind(id) = -1


C---  also delete last particle information on GRAPE which
C---  may interfere with force calculation
      i    = N+1
      m(i) = 0.d0
      DO k=1,3
         x(k,i)    = 0.d0
         v(k,i)    = 0.d0
         a(k,i)    = 0.d0
         adot(k,i) = 0.d0

C---  for shifting particles on GRAPE
         a2by18(k) = 0.d0 
         a1by6(k)  = 0.d0
         aby2(k)   = 0.d0

      ENDDO
      pot(i)       = 0.d0
      t(i)         = 0.d0
      dt(i)        = dt_min
      ind(i)       = 0
      myid(i)      = 0

C      CALL g6_set_j_particle(clusterid, i-1, ind(i),t(i),
C     $                       dt(i),m(i),a2by18,a1by6,aby2,
C     $                       v(1,i), x(1,i))


C skip at the moment: should it be called afterwards? (fip)
C---  recompute forces and time steps for all particles
C      init = reinitialize_particles() 
C            ----------------------

      delete_particle = 0

      END



 
      FUNCTION evolve_model(tend)
C---  isync removed      
      IMPLICIT NONE 
#include "src/general.inc"
#include "stopcond.inc"
C---  local vars
      INTEGER          evolve_model
      DOUBLE PRECISION tend

      INTEGER clock_init, clock_current
      INTEGER count_rate, count_max
      INTEGER is_any_condition_set
      INTEGER is_stopping_condition_enabled
      INTEGER is_number_of_steps_detection_enabled
      INTEGER is_timeout_detection_enabled
      INTEGER mpi_distribute_stopping_conditions
      INTEGER mpi_collect_stopping_conditions
      
      INTEGER get_stopping_condition_number_of_steps_parameter 
      INTEGER get_stopping_condition_timeout_parameter 
      INTEGER max_number_of_steps
      DOUBLE PRECISION timeout
      INTEGER number_of_steps_innerloop
      INTEGER next_index_for_stopping_condition
      INTEGER set_stopping_condition_info
      INTEGER stopping_index
      INTEGER reset_stopping_conditions, error
C (re)start grape (fip)
C  update local particles on grape
      CALL amuse_init_grape
      CALL update_grape(1) 

      t_end      = tend

      error = reset_stopping_conditions()
      number_of_steps_innerloop = 0

C--------------------------------------------------------------------
C---  main loop begins here
C--------------------------------------------------------------------
      CALL SYSTEM_CLOCK(clock_init, count_rate, count_max)
    
      error = mpi_distribute_stopping_conditions()
      
      error = is_stopping_condition_enabled(
     $               NUMBER_OF_STEPS_DETECTION, 
     $               is_number_of_steps_detection_enabled)
      error = is_stopping_condition_enabled(
     $               TIMEOUT_DETECTION,
     $               is_timeout_detection_enabled)
      error = get_stopping_condition_number_of_steps_parameter(
     $               max_number_of_steps)
      error = get_stopping_condition_timeout_parameter(
     $               timeout)

      

 100  IF (time_cur.GE.t_end) GOTO 900

C--- obsolete, use stopcond library instead..
         icollision  = -1       ! no collision yet
         idprimary   = -1       ! hence no ids for
         idsecondary = -1       ! collided particles

C---------------------------------------------------------------
C---  get the global minimum timestep 
C---------------------------------------------------------------
         CALL get_min_t
C             ---------

C--------------------------------------------------------------------
C---  select active particles from local particles
C--------------------------------------------------------------------
         CALL selectactive
C             ------------

C--------------------------------------------------------------------
C---  predictor (for local active particles)
C--------------------------------------------------------------------
         CALL predictor(0)
C             ---------         

C--------------------------------------------------------------------
C---  gather active particles on all PE's
C--------------------------------------------------------------------
         CALL allgather_ap
C             ------------

C--------------------------------------------------------------------
C---  get local forces on all PE's
C---  (as this is not the first call set ifirst to 0 in this call)
C--------------------------------------------------------------------
         CALL gravity(min_t,0)
C             -------

C--------------------------------------------------------------------
C---  sum up the partial forces
C--------------------------------------------------------------------
         CALL sumforces(0)      ! call with 0; just sum forces
C             ---------

C--------------------------------------------------------------------
C---  make the corrector step
C---  this also updates the time steps of active particles
C--------------------------------------------------------------------
         CALL corrector
C             ---------

C--------------------------------------------------------------------
C---  update local particles from local active particles
C--------------------------------------------------------------------
         CALL update_loc_p
C             ------------

C--------------------------------------------------------------------
C---  also update local particles on GRAPE from local active particles
C--------------------------------------------------------------------
         CALL update_grape(0)
C             ------------
      
         time_cur   = min_t

C timeout and number of steps detection goes here i quess...

         IF (is_number_of_steps_detection_enabled.GT.0) THEN
             number_of_steps_innerloop = 
     $           number_of_steps_innerloop + 1
             IF (number_of_steps_innerloop.GT.max_number_of_steps) THEN
                 stopping_index = next_index_for_stopping_condition()
                 error = set_stopping_condition_info(stopping_index, 
     $                       NUMBER_OF_STEPS_DETECTION);
             ENDIF

         ENDIF

         IF (is_timeout_detection_enabled.GT.0) THEN
             CALL SYSTEM_CLOCK(clock_current, count_rate, count_max)
             IF ((clock_current-clock_init).GT.timeout) THEN
                 stopping_index = next_index_for_stopping_condition()
                 error = set_stopping_condition_info(stopping_index,
     $                       TIMEOUT_DETECTION)
             ENDIF
         ENDIF
    
C---         CALL gather_collisions
         error = mpi_collect_stopping_conditions()
        
C---  break if collision occured
C--      IF (icollision.EQ.1) GOTO 900
         IF (is_any_condition_set().GT.0) GOTO 900
         GOTO 100   
 900  CONTINUE                  ! end of main loop

C----------------------------------

C----------------------------------
C---  synchronize system if flag is set
C--------------------------------------------------------------------

C-- Only synchronize_model if the user asked for it!
C      IF (isync.EQ.1) THEN
C
C         CALL selectsync(time_cur)
C             ----------
C         CALL sync
C             ----
C      ENDIF

      CALL gatherParts
C--------------------------------------------------------------------
C---  get energy, center of mass and momentum
C--------------------------------------------------------------------
      CALL energy
C          ------

C cello      evolve = idprimary        ! return primary collision member or -1

      
C     close grape before exiting (fip)
      CALL amuse_close_grape

      evolve_model = 0

      END

      FUNCTION synchronize_model()
      IMPLICIT NONE
#include "src/general.inc"
      INTEGER synchronize_model
    
      CALL amuse_init_grape
      CALL update_grape(1) 
      
      CALL selectsync(time_cur)
      
      CALL sync

      CALL gatherParts
      
      CALL energy
      
      synchronize_model = 0
      
      CALL amuse_close_grape


      END


      FUNCTION find_colliding_primary()
      IMPLICIT NONE
#include "src/general.inc"
      INTEGER  find_colliding_primary

      find_colliding_primary = idprimary

      RETURN
      END

      FUNCTION find_colliding_secondary(id1)
      IMPLICIT NONE
#include "src/general.inc"
      INTEGER  find_colliding_secondary,id1,idummy
      
      idummy = id1

      find_colliding_secondary = idsecondary

      RETURN
      END


       SUBROUTINE gather_collisions
#include "src/general.inc"
#ifndef NOMPI
#include <mpif.h>
       INCLUDE "mpi.inc"       
       
       integer cflag(NPEMAX),id1(NPEMAX),id2(NPEMAX),ierr
              
       CALL MPI_AllGather(icollision,1,MPI_INTEGER, cflag, 1, 
     &       MPI_INTEGER,MPI_COMM_WORLD,ierr)
       CALL MPI_AllGather(idprimary,1,MPI_INTEGER, id1, 1, 
     &       MPI_INTEGER,MPI_COMM_WORLD,ierr)
       CALL MPI_AllGather(idsecondary,1,MPI_INTEGER, id2, 1, 
     &       MPI_INTEGER,MPI_COMM_WORLD,ierr)
       
       icollision=-1
         DO i=1,n_proc
           IF(cflag(i).EQ.1) THEN
             icollision=cflag(i)
             idprimary=id1(i)
             idsecondary=id2(i)
             exit             
           ENDIF        
         ENDDO
         
#endif              
       END
       
      FUNCTION commit_parameters()
      IMPLICIT NONE
#include "src/general.inc"
      INTEGER  commit_parameters
      commit_parameters = 0
      time_cur = time_begin
      IF(g_init_once.EQ.1) THEN
        IF(g_initialized.EQ.0) THEN
            CALL initgrape
            g_initialized = 1
        ENDIF
      ENDIF
      
      END
      
      FUNCTION recommit_parameters()
      IMPLICIT NONE
#include "src/general.inc"
      INTEGER  recommit_parameters
      recommit_parameters = 0
      END
      
      
