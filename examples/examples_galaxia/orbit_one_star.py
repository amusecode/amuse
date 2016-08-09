
#---------------------------------------------------------------------------
# Script that uses galaxia and the bridge in rotating coordinates           |
# to integrate the orbit of a star.                                         |
# For a complete explanation of the possible parameters and models          |
# included in galaxia, we refer the reader to the file: user_manual_galaxia.|
# The running time may take some seconds, so be patient!                    |
#---------------------------------------------------------------------------


from amuse import datamodel
from amuse.units import quantities, constants, units
from amuse.ext.rotating_bridge import Rotating_Bridge
from amuse.community.galaxia.interface import BarAndSpirals3D
from amuse.ext.composition_methods import *
import numpy
from matplotlib import pyplot


class drift_without_gravity(object):
    """
    This class is necesary 
    to make the evolution of 
    test particles by using bridge.
    """

    def __init__(self, particles, time= 0 |units.Myr):
        self.particles=particles
        self.model_time= time
    def evolve_model(self, t_end):
        dt= t_end- self.model_time
        self.particles.position += self.particles.velocity*dt
        self.model_time= t_end
    @property
    def potential_energy(self):
        return quantities.zero
    @property 
    def kinetic_energy(self):
        return (0.5*self.particles.mass*self.particles.velocity.lengths()**2).sum()


class IntegrateOrbitStar(object):
    """
    This class makes the integration of a star in an analytical potential
    using galaxia.
    The integration methods that can be used in the Rotating Bridge are:
    LEAPFROG
    SPLIT_4TH_S_M6
    SPLIT_4TH_S_M5
    SPLIT_4TH_S_M4
    SPLIT_6TH_SS_M11
    SPLIT_6TH_SS_M13
    SPLIT_8TH_SS_M21
    SPLIT_10TH_SS_M35
    where the ordinal number stands for the order of the integrator (i.e, 4th is fourth order);
    the S for symplectic; and M corresponds to the number of times the force in computed
    (i.e., M6 means that the force is computed 6 times).
    """
    
    def __init__(self, simulation_time= 10 |units.Myr, 
                 dt_bridge=0.5 |units.Myr, 
                 method_for_rotating_bridge= SPLIT_6TH_SS_M13, 
                 initial_phase_bar= 0, 
                 initial_phase_spiral_arms= 0, 
                 pattern_speed_spiral_arms= -20 |(units.kms/units.kpc), 
                 amplitude_spiral_arms= 1100|(units.kms**2/units.kpc), 
                 number_of_arms=4,
                 separation_locus_spiral_arms= 3.12 |units.kpc,
                 tangent_pitch_angle= 0.227194425,
                 pattern_speed_bar= -50 |(units.kms/units.kpc), 
                 mass_bar= 1.1e10 |units.MSun,
                 semimajor_axis_bar= 3.12 |units.kpc,
                 axis_ratio_bar= 0.37):
        # initialization of Simulation parameters.
        self.t_end= simulation_time
        self.dt_bridge= dt_bridge
        self.method= method_for_rotating_bridge
        self.time= 0 |units.Myr
        
        #initialization of the galaxy parameters.
        self.omega= 0 | (units.kms/units.kpc)
        self.initial_phase= 0
        self.bar_phase= initial_phase_bar
        self.spiral_phase= initial_phase_spiral_arms
        self.omega_spiral= pattern_speed_spiral_arms
        self.amplitude= amplitude_spiral_arms
        self.rsp= separation_locus_spiral_arms
        self.m= number_of_arms
        self.tan_pitch_angle= tangent_pitch_angle
        self.omega_bar= pattern_speed_bar
        self.mass_bar= mass_bar
        self.aaxis_bar= semimajor_axis_bar
        self.axis_ratio_bar= axis_ratio_bar
        return
    
    def galaxy(self):
        # Model of the Galaxy.
        # In this example, the Galaxy has two-dimensional bar and spiral arms.
        # The spiral arms model is the TWA
        # The bar does not grow adiabatically
        # The axisymmetric component has its defaul values from Allen & Santillan (1990).
        galaxy= BarAndSpirals3D()
        galaxy.kinetic_energy=quantities.zero
        galaxy.potential_energy=quantities.zero
        galaxy.parameters.bar_contribution= True
        galaxy.parameters.bar_phase= self.bar_phase
        galaxy.parameters.omega_bar= self.omega_bar
        galaxy.parameters.mass_bar= self.mass_bar
        galaxy.parameters.aaxis_bar= self.aaxis_bar
        galaxy.parameters.axis_ratio_bar= self.axis_ratio_bar 
        galaxy.parameters.spiral_contribution= True
        galaxy.parameters.spiral_model=0
        galaxy.parameters.spiral_phase= self.spiral_phase
        galaxy.parameters.omega_spiral= self.omega_spiral
        galaxy.parameters.amplitude= self.amplitude
        galaxy.parameters.rsp= self.rsp
        galaxy.parameters.m= self.m
        galaxy.parameters.tan_pitch_angle= self.tan_pitch_angle
        galaxy.commit_parameters()
        self.omega= galaxy.parameters.omega_system
        self. initial_phase= galaxy.parameters.initial_phase
        return galaxy

        
    def creation_particles_noinertial(self, particles):
        # This function creates a particle set in a 
        # rotating system. This information is used by 
        # the rotating bridge
        no_inertial_system= particles.copy()
        angle= self.initial_phase + self.omega*self.time
        C1= particles.vx + self.omega*particles.y
        C2= particles.vy - self.omega*particles.x
        no_inertial_system.x = particles.x*numpy.cos(angle) + particles.y*numpy.sin(angle)
        no_inertial_system.y = -particles.x*numpy.sin(angle) + particles.y*numpy.cos(angle) 
        no_inertial_system.z = particles.z
        no_inertial_system.vx = C1*numpy.cos(angle) + C2*numpy.sin(angle) 
        no_inertial_system.vy = C2*numpy.cos(angle) - C1*numpy.sin(angle)
        no_inertial_system.vz = particles.vz
        return no_inertial_system    

    def noinertial_to_inertial(self, part_noin, part_in):
        #Transformation function from a rotating to an inertial frame.
        angle= self.initial_phase + self.omega*self.time
        C1= part_noin.vx - part_noin.y*self.omega
        C2= part_noin.vy + part_noin.x*self.omega
        part_in.x= part_noin.x*numpy.cos(angle)-part_noin.y*numpy.sin(angle)
        part_in.y= part_noin.x*numpy.sin(angle)+part_noin.y*numpy.cos(angle)
        part_in.z= part_noin.z
        part_in.vx= C1*numpy.cos(angle) - C2*numpy.sin(angle)
        part_in.vy= C1*numpy.sin(angle) + C2*numpy.cos(angle)
        part_in.vz= part_noin.vz
        return


    def get_orbit_of_the_star(self, particle_set):
        # Function that integrates  the orbit of the star.
        # input: particle_set -> defined in an inertial system
        # steps:
        # 1. The Galaxy model is created
        # 2. A copy of particle_set is created in the rotating frame
        # where the stellar motion is computed 
        # 3. The Galaxy and the particles in the rotating frame are coupled through the rotating Bridge
        # 4. The evolution of the system is made
        # 5. The positions and velocities in particle_set are updated by transforming the
        # phase-space coordinates back to the inertial frame.   
        X=[]
        Y=[]
        
        MW= self.galaxy()
        particles_in_rotating_coordinates= self.creation_particles_noinertial(particle_set)
        gravless= drift_without_gravity(particles_in_rotating_coordinates)
        
        system= Rotating_Bridge(self.omega, timestep= self.dt_bridge, verbose= False, method= self.method)
        system.add_system(gravless, (MW,), False)
        system.add_system(MW, (), False)

        
        while (self.time< self.t_end-self.dt_bridge/2):
            X.append(particle_set[0].x.value_in(units.kpc))
            Y.append(particle_set[0].y.value_in(units.kpc))
            self.time += self.dt_bridge
            system.evolve_model(self.time)
            self.noinertial_to_inertial(particles_in_rotating_coordinates, particle_set)
      
        return X, Y


def plot_orbit(position_x, position_y):

    figure= pyplot.figure(figsize=(6,6))
    ax= figure.add_subplot(111)
    ax.plot(position_x, position_y)
    ax.set_xlabel('X [kpc]')
    ax.set_ylabel('Y [kpc]')
    ax.set_title('Orbit of a star with initial position at (%g,%g) kpc \n from the Galactic center'%(position_x[0], position_y[0]))
    ax.set_xlim(-10,10)
    ax.set_ylim(-10, 10)
    pyplot.show()

    
if __name__ in('__main__', '__plot__'):

    
    star= datamodel.Particles(1)
    star[0].mass= 1 | units.MSun
    star[0].radius= 1 |units.RSun
    star[0].position= [-8.5, 0, 0] | units.kpc
    star[0].velocity= [-7.3, 100, 0] | units.kms
                                       
    
    integrator= IntegrateOrbitStar(simulation_time= 250 |units.Myr, 
                             dt_bridge=0.5 |units.Myr, 
                             method_for_rotating_bridge= SPLIT_4TH_S_M6, 
                             initial_phase_bar= -0.35, 
                             initial_phase_spiral_arms= -0.35, 
                             pattern_speed_spiral_arms= 20 |(units.kms/units.kpc), 
                             amplitude_spiral_arms= 1100|(units.kms**2/units.kpc), 
                             number_of_arms=4,
                             pattern_speed_bar= 50 |(units.kms/units.kpc), 
                             mass_bar= 1.1e10 |units.MSun)
    
    x, y= integrator.get_orbit_of_the_star(star)
    
    plot_orbit(x,y)
    
