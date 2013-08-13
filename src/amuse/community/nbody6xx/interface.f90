*
*             N B O D Y 6++ for AMUSE
*             ***********************
*
*       Regularized AC N-body code with triple & binary collisions.
*       --------------------------------------------------------
*
*       Hermite integration scheme with block-steps (V 4.0.0 April/99).
*       ------------------------------------------------------------------
*
*       Developed by Sverre Aarseth, IOA, Cambridge.
*       ............................................
*       Message Passing Version NBODY6++ for Massively Parallel Systems
*       Developed by Rainer Spurzem, ARI, Heidelberg
*       Ported to AMUSE by Maxwell X. Tsai, NAOC, Beijing
*
      MODULE AMUSE_INTERFACE
      integer :: last_index
      integer :: nb6_init_called
      CONTAINS

      FUNCTION main(input, output)
      include "src/common6.h"
      INTEGER input, output
      INTEGER main
      print *, "calling main()"
      amusein = 1
      CALL NBODY6
      output = input
      main = 1
      END FUNCTION
      
      FUNCTION initialize_code()
        include "src/common6.h"
        INTEGER :: initialize_code
        INTEGER :: J, VAL, N_DUMMY, ioerror
        print *, "initialization..."
        OPEN (UNIT=99, FILE='inbody6xx.nput', IOSTAT=ioerror)
        
        if (ioerror.EQ.0) then
            READ (99,*,IOSTAT=ioerror)  KSTART_AMUSE, TCOMP_AMUSE, TCRITp, 
     &      isernb, iserreg
        else 
            KSTART_AMUSE = 1
            TCOMP_AMUSE = 10000.0
            TCRITp = 1.E6
            isernb = 40
            iserreg = 40  
        end if
        
        PRINT *, KSTART_AMUSE, TCOMP_AMUSE, TCRITp, isernb, iserreg
        if (ioerror.EQ.0)  then
            READ (99,*,IOSTAT=ioerror)  N_DUMMY, NFIX, NCRIT, 
     &       NRAND_AMUSE, NNBOPT, NRUN
        end if
        
        if (ioerror.NE.0)  then
            N_DUMMY = 1000
            NFIX = 1
            NCRIT = 10
            NRAND_AMUSE = 4353
            NNBOPT = 80
            NRUN = 1
        end if
        
        PRINT *,  N, NFIX, NCRIT, NRAND_AMUSE, NNBOPT, NRUN
        if (ioerror.EQ.0)  then
            READ (99,*,IOSTAT=ioerror)  ETAI, ETAR, RS0_AMUSE, DTADJ, DELTAT, TCRIT,
     &              QE, RBAR, ZMBAR
        end if
        
        if (ioerror.NE.0)  then
            ETAI = 0.05
            ETAR = 0.05
            RS0_AMUSE = 0.12
            DTADJ = 1.0
            DELTAT = 1.0
            TCRIT = 100
            QE = 2.0E-05
            RBAR = 1.0
            ZMBAR = 0.7    
        end if
        
        if (ioerror.EQ.0)  then
            READ (99,*,IOSTAT=ioerror)  (KZ(J),J=1,50)
        end if
        
        if (ioerror.NE.0)  then
            KZ(1:50) = (/ 
     &           1, 2, 1, 0, 1, 1, 4, 0, 0, 2,
     &           1, 0, 0, 0, 2, 1, 0, 0, 3, 2,
     &           1, 0, 2, 0, 0, 2, 0, 0, 0, 2,
     &           0, 0, 2, 0, 1, 0, 1, 1, 0, 1,
     &           0, 0, 0, 0, 0, 0, 0, 0, 0, 0 
     &       /)
        end if
        
        PRINT *, (KZ(J),J=1,50)
    
        if (ioerror.EQ.0)  then
            READ (99,*,IOSTAT=ioerror)  DTMIN, RMIN, ETAU, ECLOSE, 
     &       GMIN, GMAX
        end if
        
        if (ioerror.NE.0)  then
            DTMIN = 1.0E-04 
            RMIN = 0.01 
            ETAU = 0.1 
            ECLOSE = 1.0 
            GMIN = 1.0E-06 
            GMAX = 0.01
        end if
        
        if (ioerror.EQ.0)  then
            READ (99,*,IOSTAT=ioerror) ALPHAS, BODY1, BODYN, NBIN0,
     &             EPOCH0, DTPLOT
        end if
        
        if (ioerror.NE.0)  then
            ALPHAS = 2.35
            BODY1 =  20.0 
            BODYN = 0.1 
            NBIN0 = 0 
            NHI0 = 0
            ZMET =  0.0 
            EPOCH0 = 0 
            DTPLOT = 0.0
        end if
        
        if (ioerror.EQ.0)  then
            READ (99,*,IOSTAT=ioerror)  Q, VXROT, VZROT, RSPH2
        end if
        
        if (ioerror.NE.0)  then
            Q = 0.5
            VXROT =  0.0 
            VZROT = 0.0
            RSPH2 =  0.0
        end if
        
        if (ioerror.EQ.0)  then
            READ (99,*,IOSTAT=ioerror)  SEMI0, ECC0, RATIO,RANGE,
     &         NSKIP, IDORM
        end if
        
        if (ioerror.NE.0)  then
            SEMI0 = 0.005
            ECC0 =  -1.0
            RATIO =  1.0 
            RANGE = 5.0
            NSKIP =  5
            IDORM = 0
        end if
        
        CLOSE (99, IOSTAT=ioerror)
        
        N = last_index - 1
        amusein = 1
        last_index = 1
        nb6_init_called = 0
        initialize_code = 0
      END FUNCTION
      
      FUNCTION new_particle(index_of_the_particle, mass_amuse, x_amuse, 
     & y_amuse, z_amuse, vx_amuse, vy_amuse, vz_amuse, r_amuse)
        INCLUDE 'src/common6.h'
        INTEGER :: index_of_the_particle
        DOUBLE PRECISION :: mass_amuse, x_amuse, y_amuse, z_amuse
        DOUBLE PRECISION :: vx_amuse, vy_amuse, vz_amuse, r_amuse
        INTEGER :: new_particle
        index_of_the_particle = last_index
        last_index = last_index + 1
        X(1, index_of_the_particle) = x_amuse
        X(2, index_of_the_particle) = y_amuse
        X(3, index_of_the_particle) = z_amuse
        XDOT(1, index_of_the_particle) = vx_amuse
        XDOT(2, index_of_the_particle) = vy_amuse
        XDOT(3, index_of_the_particle) = vz_amuse
        BODY(index_of_the_particle) = mass_amuse
        new_particle = 0
      END FUNCTION
      
      FUNCTION delete_particle(index_of_the_particle)
        INCLUDE 'src/common6.h'
        INTEGER :: index_of_the_particle
        INTEGER :: delete_particle
        delete_particle = -2
      END FUNCTION
      
      FUNCTION get_number_of_particles(number_of_particles)
        INCLUDE 'src/common6.h'
        INTEGER :: number_of_particles
        INTEGER :: get_number_of_particles
        number_of_particles = N
        get_number_of_particles = 0
      END FUNCTION
      
      FUNCTION get_eps2(epsilon_squared)
        INCLUDE 'src/common6.h'
        DOUBLE PRECISION :: epsilon_squared
        INTEGER :: get_eps2
        epsilon_squared = 0
        get_eps2 = -1
      END FUNCTION
      
      FUNCTION set_eps2(epsilon_squared)
        INCLUDE 'src/common6.h'
        DOUBLE PRECISION :: epsilon_squared
        INTEGER :: set_eps2
        set_eps2 = -1
      END FUNCTION
      
      FUNCTION set_radius(id_amuse, r_amuse)
        INCLUDE 'src/common6.h'
        INTEGER :: id_amuse
        DOUBLE PRECISION :: r_amuse
        INTEGER :: set_radius
        set_radius = 0 
      END FUNCTION
      
      FUNCTION set_velocity(id_particle, vx_amuse, vy_amuse, vz_amuse)
        INCLUDE 'src/common6.h'
        INTEGER :: id
        DOUBLE PRECISION :: vx_amuse, vy_amuse, vz_amuse
        INTEGER :: set_velocity
        XDOT(1, id_particle) = vx_amuse
        XDOT(2, id_particle) = vy_amuse
        XDOT(3, id_particle) = vz_amuse
        set_velocity = 0
      END FUNCTION
      
      FUNCTION get_index_of_first_particle(index_of_the_particle)
        INCLUDE 'src/common6.h'
        INTEGER :: index_of_the_particle
        INTEGER :: get_index_of_first_particle
        index_of_the_particle = 1
        get_index_of_first_particle = 0
      END FUNCTION
      
      FUNCTION get_index_of_next_particle(index_of_the_particle,
     &    index_of_the_next_particle)
        INCLUDE 'src/common6.h'
        INTEGER :: index_of_the_particle, index_of_the_next_particle
        INTEGER :: get_index_of_next_particle
        IF (index_of_the_particle .lt. N) THEN
            index_of_the_next_particle = index_of_the_particle + 1
        ELSE
            index_of_the_next_particle = 1
        END IF
        get_index_of_next_particle = 0
      END FUNCTION
      
      FUNCTION get_state(index_of_the_particle, mass_amuse,  
     & x_amuse, y_amuse, z_amuse, vx_amuse, vy_amuse, vz_amuse, r_amuse)
        INCLUDE 'src/common6.h'
        INTEGER :: index_of_the_particle
        DOUBLE PRECISION :: mass_amuse, x_amuse, y_amuse, z_amuse
        DOUBLE PRECISION :: vx_amuse, vy_amuse, vz_amuse, r_amuse
        INTEGER :: get_state
        x_amuse = X(1, index_of_the_particle)
        y_amuse = X(2, index_of_the_particle)
        z_amuse = X(3, index_of_the_particle)
        vx_amuse = XDOT(1, index_of_the_particle)
        vy_amuse = XDOT(2, index_of_the_particle)
        vz_amuse = XDOT(3, index_of_the_particle)
        mass_amuse = BODY(index_of_the_particle)
        get_state = 0
      END FUNCTION
      
      FUNCTION set_state(index_of_the_particle, mass_amuse,
     & x_amuse, y_amuse, z_amuse, vx_amuse, vy_amuse, vz_amuse, r_amuse)
        INCLUDE 'src/common6.h'
        INTEGER :: index_of_the_particle
        DOUBLE PRECISION :: mass_amuse, x_amuse, y_amuse, z_amuse
        DOUBLE PRECISION :: vx_amuse, vy_amuse, vz_amuse, r_amuse
        INTEGER :: set_state
        X(1, index_of_the_particle) = x_amuse
        X(2, index_of_the_particle) = y_amuse
        X(3, index_of_the_particle) = z_amuse
        XDOT(1, index_of_the_particle) = vx_amuse
        XDOT(2, index_of_the_particle) = vy_amuse
        XDOT(3, index_of_the_particle) = vz_amuse
        BODY(index_of_the_particle) = mass_amuse
        PRINT *, "setting state for particle ", index_of_the_particle
        set_state = 0
      END FUNCTION
      
      FUNCTION get_mass(index_of_the_particle, mass_amuse)
        INCLUDE 'src/common6.h'
        INTEGER :: index_of_the_particle
        DOUBLE PRECISION :: mass_amuse
        INTEGER :: get_mass
        mass_amuse = BODY(index_of_the_particle)
        get_mass = 0
      END FUNCTION
      
      FUNCTION set_mass(index_of_the_particle, mass_amuse)
        INCLUDE 'src/common6.h'
        INTEGER :: index_of_the_particle
        DOUBLE PRECISION :: mass_amuse
        INTEGER :: set_mass
        BODY(index_of_the_particle) = mass_amuse
        set_mass = 0
      END FUNCTION
      
      FUNCTION get_position(index_of_the_particle, x_amuse, y_amuse,
     &    z_amuse)
        INCLUDE 'src/common6.h'
        INTEGER :: index_of_the_particle
        DOUBLE PRECISION :: x_amuse, y_amuse, z_amuse
        INTEGER :: get_position
        x_amuse = X(1, index_of_the_particle)
        y_amuse = X(2, index_of_the_particle)
        z_amuse = X(3, index_of_the_particle)
        get_position = 0
      END FUNCTION
      
      FUNCTION set_position(index_of_the_particle, x_amuse, y_amuse,
     &     z_amuse)
        INCLUDE 'src/common6.h'
        INTEGER :: index_of_the_particle
        DOUBLE PRECISION :: x_amuse, y_amuse, z_amuse
        INTEGER :: set_position
        X(1, index_of_the_particle) = x_amuse
        X(2, index_of_the_particle) = y_amuse
        X(3, index_of_the_particle) = z_amuse
        set_position = 0
      END FUNCTION
      
      FUNCTION set_begin_time(input)
        INCLUDE 'src/common6.h'
        DOUBLE PRECISION :: input
        INTEGER :: set_being_time
      END FUNCTION
      
      FUNCTION get_acceleration(index_of_the_particle, ax_amuse, 
     & ay_amuse, az_amuse)
        INCLUDE 'src/common6.h'
        INTEGER :: index_of_the_particle
        DOUBLE PRECISION :: ax_amuse, ay_amuse, az_amuse
        INTEGER :: get_acceleration
        ax_amuse = FI(1, index_of_the_particle) + FR(1,
     &   index_of_the_particle)
        ay_amuse = FI(2, index_of_the_particle) + FR(2,
     &   index_of_the_particle)
        az_amuse = FI(3, index_of_the_particle) + FR(3,
     &   index_of_the_particle)
        get_acceleration = 0
      END FUNCTION
      
      FUNCTION set_acceleration(index_of_the_particle, ax, ay, az)
        INTEGER :: index_of_the_particle
        DOUBLE PRECISION :: ax, ay, az
        INTEGER :: set_acceleration
        set_acceleration = -2
      END FUNCTION
      
      FUNCTION get_potential(index_of_the_particle, potential)
        INCLUDE 'src/common6.h'
        INTEGER :: index_of_the_particle
        DOUBLE PRECISION :: potential
        INTEGER :: get_potential
      END FUNCTION
      
      FUNCTION commit_particles()
        INCLUDE 'src/common6.h'
        INTEGER :: commit_particles
        N = last_index - 1
        print *, "calling commit particles"
        commit_particles = 0
      END FUNCTION
      
      FUNCTION recommit_particles()
        INCLUDE 'src/common6.h'
        INTEGER :: recommit_particles
      END FUNCTION
      
      FUNCTION commit_parameters()
        INCLUDE 'src/common6.h'
        INTEGER :: commit_parameters

        N = last_index - 1
*        CALL NBODY6
        commit_parameters = 0
      END FUNCTION
      
      FUNCTION recommit_parameters()
        INCLUDE 'src/common6.h'
        INTEGER :: recommit_parameters
      END FUNCTION
      
      FUNCTION get_time(t_amuse)
        INCLUDE 'src/common6.h'
        DOUBLE PRECISION :: t_amuse
        INTEGER :: get_time
        t_amuse = TTOT
        get_time = 0
      END FUNCTION
      
      FUNCTION get_time_step(time_step)
        INCLUDE 'src/common6.h'
        DOUBLE PRECISION :: time_step
        INTEGER :: get_time_step
      END FUNCTION
      
      FUNCTION get_kinetic_energy(kinetic_energy)
        INCLUDE 'src/common6.h'
        DOUBLE PRECISION :: kinetic_energy
        INTEGER :: get_kinetic_energy
      END FUNCTION
      
      FUNCTION get_potential_energy(potential_energy)
        INCLUDE 'src/common6.h'
        DOUBLE PRECISION :: potential_energy
        INTEGER :: get_potential_energy
        potential_energy = POT
      END FUNCTION
      
      FUNCTION get_center_of_mass_velocity(vx, vy, vz)
        INCLUDE 'src/common6.h'
        DOUBLE PRECISION :: vx, vy, vz
        INTEGER :: get_center_of_mass_velocity
      END FUNCTION
      
      FUNCTION get_velocity(id_particle, vx_amuse, vy_amuse, vz_amuse)
        INCLUDE 'src/common6.h'
        INTEGER :: id_particle
        DOUBLE PRECISION :: vx_amuse, vy_amuse, vz_amuse
        INTEGER :: get_velocity
        vx_amuse = XDOT(1, id_particle)
        vy_amuse = XDOT(2, id_particle)
        vz_amuse = XDOT(3, id_particle)
        get_velocity = 0
      END FUNCTION
      
      FUNCTION get_center_of_mass_position(x, y, z)
        DOUBLE PRECISION :: x, y, z
        INTEGER :: get_center_of_mass_position
      END FUNCTION
      
      FUNCTION get_total_mass(mass_amuse)
        INCLUDE 'src/common6.h'
        DOUBLE PRECISION :: mass_amuse
        INTEGER :: get_total_mass
        mass_amuse = ZMASS
*        mass_amuse = ZMBAR
        get_total_mass = 0
      END FUNCTION
      
      FUNCTION get_total_radius(r_amuse)
        INCLUDE 'src/common6.h'
        DOUBLE PRECISION :: r_amuse
        INTEGER :: get_total_radius
        r_amuse = RBAR
        get_total_radius = 0
      END FUNCTION
      
      FUNCTION cleanup_code()
        INCLUDE 'src/common6.h'
        INTEGER :: cleanup_code
      END FUNCTION
      
      FUNCTION evolve_model(t_end)
        INCLUDE 'src/common6.h'
        DOUBLE PRECISION :: t_end
        INTEGER :: evolve_model
        print *, "toff = ", toff
        if(nb6_init_called.eq.0) then
            nb6_init_called = 1
            CALL NBODY6
        end if
        IF(t_end.gt.0) THEN
*            TCRIT = t_end
*            amusein = 1
*            PRINT *, "CALLING integration"
            CALL INTAMUSE
        END IF
      END FUNCTION
      
      FUNCTION get_begin_time(output)
        INCLUDE 'src/common6.h'
        DOUBLE PRECISION :: output
        INTEGER :: get_begin_time
      END FUNCTION
      
      FUNCTION get_radius(id, r_amuse)
        INCLUDE 'src/common6.h'
        INTEGER :: id
        DOUBLE PRECISION :: r_amuse
        INTEGER :: get_radius
        r_amuse = 0
        get_radius = 0
      END FUNCTION
      
      FUNCTION synchronize_model() 
        INCLUDE 'src/common6.h'
        INTEGER synchronize_model
        synchronize_model = 0
      END FUNCTION

      FUNCTION body_print() 
        INCLUDE 'src/common6.h'
        INTEGER :: body_print, pid_body
*        do pid_body = 1, N
*            print *, pid_body, BODY(pid_body)
*        end do
        body_print = 0
      END FUNCTION

      END MODULE AMUSE_INTERFACE
