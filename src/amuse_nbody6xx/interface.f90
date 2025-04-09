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
*       Ported to AMUSE by Maxwell X. TSAI/CAI, NAOC/KIAA, Beijing
*       (2013, 2014)
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

            KSTART_AMUSE = 1
            TCOMP_AMUSE = 10000.0
            TCRITp = 1.E6
            isernb = 40
            iserreg = 40  
        
        PRINT *, KSTART_AMUSE, TCOMP_AMUSE, TCRITp, isernb, iserreg
        
            N_DUMMY = 1000
            NFIX = 1
            NCRIT = 10
            NRAND_AMUSE = 4353
            NNBOPT = 2
            NRUN = 1
        
        PRINT *,  N, NFIX, NCRIT, NRAND_AMUSE, NNBOPT, NRUN

            ETAI = 0.05
            ETAR = 0.05
            RS0_AMUSE = 0.12
            DTADJ = 1.0
            DELTAT = 1.0
            TCRIT = 100
            QE = 2.0E-05
            RBAR = 1.0
            ZMBAR = 0.7    
        
            KZ(1:50) = (/ 
     &           1, 0, 1, 0, 1, 1, 4, 0, 0, 2,
     &           1, 0, 0, 0, 2, 1, 0, 0, 3, 2,
     &           1, 99, 2, 0, 0, 2, 0, 0, 0, 2,
     &           0, 0, 2, 0, 1, 0, 1, 1, 0, 1,
     &           0, 0, 0, 0, 0, 0, 0, 0, 0, 0 
     &       /)
        
        PRINT *, (KZ(J),J=1,50)
    
            DTMIN = 1.0E-04 
            RMIN = 0.01 
            ETAU = 0.1 
            ECLOSE = 1.0 
            GMIN = 1.0E-06 
            GMAX = 0.01
        
            ALPHAS = 2.35
            BODY1 =  20.0 
            BODYN = 0.1 
            NBIN0 = 0 
            NHI0 = 0
            ZMET =  0.0 
            EPOCH0 = 0 
            DTPLOT = 0.0
        
            QVIR = 0.5
            VXROT =  0.0 
            VZROT = 0.0
            RSPH2 =  0.0
        
            SEMI0 = 0.005
            ECC0 =  -1.0
            RATIO =  1.0 
            RANGE = 5.0
            NSKIP =  5
            IDORM = 0
        
        N = last_index - 1
        amusein = 1
        last_index = 1 + 2*NBIN0
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
        if (N.LT.0) then
            number_of_particles = last_index - 1
        else
            number_of_particles = N
        end if
        get_number_of_particles = 0
      END FUNCTION
      
      FUNCTION get_eps2(epsilon_squared)
        INCLUDE 'src/common6.h'
        DOUBLE PRECISION :: epsilon_squared
        INTEGER :: get_eps2
        epsilon_squared = 0
        get_eps2 = 0
      END FUNCTION
      
      FUNCTION set_eps2(epsilon_squared)
        INCLUDE 'src/common6.h'
        DOUBLE PRECISION :: epsilon_squared
        INTEGER :: set_eps2
        set_eps2 = 0
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
        set_begin_time = 0
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
      
      FUNCTION get_potential(index_of_the_particle, potential)
        INCLUDE 'src/common6.h'
        INTEGER :: index_of_the_particle
        DOUBLE PRECISION :: potential
        INTEGER :: get_potential
        potential = PHIDBL(index_of_the_particle)
        get_potential = 0
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
        recommit_particles=0
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
        recommit_parameters=0
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
        time_step = ETAI
        get_time_step = 0
      END FUNCTION
      
      FUNCTION get_kinetic_energy(kinetic_energy)
        INCLUDE 'src/common6.h'
        DOUBLE PRECISION :: kinetic_energy
        INTEGER :: get_kinetic_energy
        kinetic_energy= ZKIN
        get_kinetic_energy=0
      END FUNCTION
      
      FUNCTION get_potential_energy(potential_energy)
        INCLUDE 'src/common6.h'
        DOUBLE PRECISION :: potential_energy
        INTEGER :: get_potential_energy
        potential_energy = POT
        get_potential_energy=0
      END FUNCTION
      
      FUNCTION get_center_of_mass_velocity(cm_vx, cm_vy, cm_vz)
        INCLUDE 'src/common6.h'
        DOUBLE PRECISION :: cm_vx, cm_vy, cm_vz
        INTEGER :: get_center_of_mass_velocity
        DOUBLE PRECISION :: total_mass
        total_mass = ZMASS
        cm_vx = 0
        cm_vy = 0 
        cm_vz = 0
        DO idx = IFIRST, NTOT
            cm_vx = cm_vx + XDOT(1, idx) * BODY(idx)
            cm_vy = cm_vy + XDOT(2, idx) * BODY(idx)
            cm_vz = cm_vz + XDOT(3, idx) * BODY(idx)
        ENDDO
        cm_vx = cm_vx / total_mass
        cm_vy = cm_vy / total_mass
        cm_vz = cm_vz / total_mass
        get_center_of_mass_velocity = 0
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
      
      FUNCTION get_center_of_mass_position(cm_x, cm_y, cm_z)
        INCLUDE 'src/common6.h'
        DOUBLE PRECISION :: cm_x, cm_y, cm_z
        INTEGER :: get_center_of_mass_position
        DOUBLE PRECISION :: total_mass
        total_mass = ZMASS
        cm_x = 0
        cm_y = 0
        cm_z = 0
        DO idx = IFIRST, NTOT
            cm_x = cm_x + X(1, idx) * BODY(idx)
            cm_y = cm_y + X(2, idx) * BODY(idx)
            cm_z = cm_z + X(3, idx) * BODY(idx)
        ENDDO
        cm_x = cm_x / total_mass
        cm_y = cm_y / total_mass
        cm_z = cm_z / total_mass
        get_center_of_mass_position = 0
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
        cleanup_code=0
      END FUNCTION
      
      FUNCTION evolve_model(t_end)
        INCLUDE 'src/common6.h'
        DOUBLE PRECISION :: t_end
        INTEGER :: evolve_model
        if(nb6_init_called.eq.0) then
            nb6_init_called = 1
            CALL NBODY6
        end if
        DO WHILE(TTOT.lt.t_end)
*            TCRIT = t_end
*            amusein = 1
*            PRINT *, "CALLING integration"
            CALL INTAMUSE
        END DO
        evolve_model=0
      END FUNCTION
      
      FUNCTION get_begin_time(output)
        INCLUDE 'src/common6.h'
        DOUBLE PRECISION :: output
        INTEGER :: get_begin_time
        output = 0
        get_begin_time = 0
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

      FUNCTION set_kz(kz_option, val)
        INCLUDE 'src/common6.h'
        DOUBLE PRECISION :: val
        INTEGER :: set_kz, kz_option
        IF (kz_option.GT.0 .AND. kz_option.LE.50) THEN
            KZ(kz_option) = val
        ENDIF
        set_kz = 0
      END FUNCTION

      FUNCTION get_kz(kz_option, val)
        INCLUDE 'src/common6.h'
        DOUBLE PRECISION :: val, get_kz
        INTEGER :: kz_option
        IF (kz_option.GT.0 .AND. kz_option.LE.50) THEN
            val = KZ(kz_option)
        ENDIF
        get_kz = 0
      END FUNCTION

      FUNCTION get_eta(timestep_parameter)
        INCLUDE 'src/common6.h'
        DOUBLE PRECISION :: timestep_parameter
        INTEGER :: get_eta
        timestep_parameter = ETAI
        get_eta = 0
      END FUNCTION

      FUNCTION set_eta(timestep_parameter)
        INCLUDE 'src/common6.h'
        DOUBLE PRECISION :: timestep_parameter
        INTEGER :: set_eta
        ETAI = timestep_parameter
        ETAR = timestep_parameter
        set_eta = 0
      END FUNCTION

      FUNCTION get_etai(timestep_parameter)
        INCLUDE 'src/common6.h'
        DOUBLE PRECISION :: timestep_parameter
        INTEGER :: get_etai
        timestep_parameter = ETAI
        get_etai = 0
      END FUNCTION

      FUNCTION set_etai(timestep_parameter)
        INCLUDE 'src/common6.h'
        DOUBLE PRECISION :: timestep_parameter
        INTEGER :: set_etai
        ETAI = timestep_parameter
        set_etai = 0
      END FUNCTION

      FUNCTION get_etar(timestep_parameter)
        INCLUDE 'src/common6.h'
        DOUBLE PRECISION :: timestep_parameter
        INTEGER :: get_etar
        timestep_parameter = ETAR
        get_etar = 0
      END FUNCTION

      FUNCTION set_etar(timestep_parameter)
        INCLUDE 'src/common6.h'
        DOUBLE PRECISION :: timestep_parameter
        INTEGER :: set_etar
        ETAR = timestep_parameter
        set_etar = 0
      END FUNCTION

      FUNCTION set_rbar(val)
        INCLUDE 'src/common6.h'
        DOUBLE PRECISION :: val
        INTEGER :: set_rbar
        RBAR = val
        set_rbar = 0
      END FUNCTION

      FUNCTION get_rbar(val)
        INCLUDE 'src/common6.h'
        DOUBLE PRECISION :: val
        INTEGER :: get_rbar
        val = RBAR
        get_rbar = 0
      END FUNCTION

      FUNCTION set_zmbar(val)
        INCLUDE 'src/common6.h'
        DOUBLE PRECISION :: val
        INTEGER :: set_zmbar
        ZMBAR = val
        set_zmbar = 0
      END FUNCTION

      FUNCTION get_zmbar(val)
        INCLUDE 'src/common6.h'
        DOUBLE PRECISION :: val
        INTEGER :: get_zmbar
        val = ZMBAR
        get_zmbar = 0
      END FUNCTION

      FUNCTION set_qe(val)
        INCLUDE 'src/common6.h'
        DOUBLE PRECISION :: val
        INTEGER :: set_qe
        QE = val
        set_qe = 0
      END FUNCTION

      FUNCTION get_qe(val)
        INCLUDE 'src/common6.h'
        DOUBLE PRECISION :: val
        INTEGER :: get_qe
        val = QE
        get_qe = 0
      END FUNCTION

      FUNCTION get_gravity_at_point(eps1, x1, y1, z1, fx, fy, fz,
     &      number_of_points)
        INCLUDE 'src/common6.h'
        INTEGER, intent(IN) :: number_of_points
        DOUBLE PRECISION :: eps1(number_of_points)
        DOUBLE PRECISION :: x1(number_of_points), y1(number_of_points)
        DOUBLE PRECISION :: z1(number_of_points)
        DOUBLE PRECISION :: fx(number_of_points), fy(number_of_points)
        DOUBLE PRECISION :: fz(number_of_points)
        DOUBLE PRECISION :: r2, r2i, ri, mri, mr3i
        INTEGER get_gravity_at_point 
        INTEGER idx, ipart

        r2 = 0
        r2i = 0
        ri = 0
        mri = 0
        mr3i = 0
        fx = 0
        fy = 0 
        fz = 0
        print *, 'number_of_points', number_of_points
        DO ipart = 1, number_of_points
        DO idx = IFIRST, NTOT 
           r2 = X(1,idx)**2 + X(2,idx)**2 + X(3,idx)**2
           r2i = 1.0/(r2 + eps1(ipart))
           ri = sqrt(r2i)
           mri = BODY(idx) * ri
           mr3i = mri * r2i
           fx(ipart) = fx(ipart) + mr3i * (X(1,idx)-x1(ipart))
           fy(ipart) = fy(ipart) + mr3i * (X(2,idx)-x1(ipart))
           fz(ipart) = fz(ipart) + mr3i * (X(3,idx)-x1(ipart))
        ENDDO
        PRINT *, "fx=", fx(ipart), "fy=", fy(ipart), "fz=", fz(ipart)
        ENDDO
        get_gravity_at_point = 0 
      END FUNCTION


      FUNCTION get_potential_at_point(eps1, x1, y1, z1, phi,
     &     number_of_points)
        INCLUDE 'src/common6.h'
        INTEGER, intent(IN) :: number_of_points
        DOUBLE PRECISION :: eps1(number_of_points)
        DOUBLE PRECISION :: x1(number_of_points), y1(number_of_points)
        DOUBLE PRECISION :: z1(number_of_points)
        DOUBLE PRECISION :: phi(number_of_points)
        DOUBLE PRECISION :: r2, r2i, ri
        INTEGER get_potential_at_point
        INTEGER idx, ipart

        phi = 0
        r2 = 0
        r2i = 0
        ri = 0

        DO ipart = 1, number_of_points
        DO idx = IFIRST, NTOT
           r2 = X(1,idx)**2 + X(2,idx)**2 + X(3,idx)**2
           r2i = 1.0/(r2 + eps1(ipart))
           ri = sqrt(r2i)
           phi(ipart) = phi(ipart) - BODY(idx) * ri
        ENDDO
        PRINT *, "potential phi = ", phi(ipart)
        ENDDO
        get_potential_at_point = 0
      END FUNCTION


      END MODULE AMUSE_INTERFACE
