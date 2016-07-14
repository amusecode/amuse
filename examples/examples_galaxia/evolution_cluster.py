#---------------------------------------------------------------------------
# Script that uses galaxia and the bridge in rotating coordinates           |
# to calculate the evolution of a star cluster in the Galaxy.               |    
# The physical processes involved are self gravity, stellar evolution       |
# and the external tidal field of the Galaxy.                               |
# For a complete explanation of the possible parameters and models          |
# included in galaxia, we refer the reader to the file: user_manual_galaxia.|
# The running time may take some seconds, so be patient!                    |
#---------------------------------------------------------------------------


import numpy 
from amuse.couple import bridge
from amuse.units import units, constants, quantities, nbody_system
from amuse.ic.brokenimf import new_broken_power_law_mass_distribution
from amuse.ic.plummer import new_plummer_sphere
from amuse.community.huayno.interface import Huayno
from amuse.community.seba.interface import SeBa
from amuse.ext.composition_methods import *
from amuse.ext.rotating_bridge import Rotating_Bridge
from amuse.community.galaxia.interface import BarAndSpirals3D
from matplotlib import pyplot

def create_cluster_with_IMF(N=100, radius=3 |units.parsec):
    '''
    Creation of a cluster following a Kroupa IMF
    '''
    masses=new_broken_power_law_mass_distribution(N,
                                                  mass_boundaries= [0.08, 0.5, 100] |units.MSun, 
                                                  alphas= [-1.3,-2.3] )
    convert_nbody=nbody_system.nbody_to_si(masses.sum(),radius)
    cluster=new_plummer_sphere(N, convert_nbody)
    cluster.mass= masses
    cluster.move_to_center()
    cluster.scale_to_standard(convert_nbody)
    return cluster
    

class RealisticEvolutionCluster(object):
    """
    This class makes the integration of a star cluster 
    in an analytical potential using galaxia.
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
    
    def __init__(self, simulation_time= 100 |units.Myr,
                 dt_bridge= 1 |units.Myr, 
                 method_for_rotating_bridge= LEAPFROG, 
                 initial_phase_main_spiral_arms= 0, 
                 pattern_speed_main_spiral_arms= 20 |(units.kms/units.kpc), 
                 amplitude_main_spiral_arms= 1100|(units.kms**2/units.kpc), 
                 number_of_main_spiral_arms=2,
                 tangent_pitch_angle_main_spiral_arms= 0.227194425,
                 initial_phase_secondary_spiral_arms= 200*(numpy.pi/180),
                 pattern_speed_secondary_spiral_arms= 15 |(units.kms/units.kpc),
                 amplitude_secondary_spiral_arms= 880 |(units.kms**2/units.kpc),
                 tangent_pitch_angle_secondary_spiral_arms= numpy.tan((-14*numpy.pi)/180.),
                 number_of_secondary_spiral_arms= 2,
                 separation_locus_spiral_arms= 3.12 |units.kpc,
                 initial_phase_bar= 0, 
                 pattern_speed_bar= 40 |(units.kms/units.kpc),
                 mass_bar= 1.2e10 |units.MSun,
                 semimajor_axis_bar= 3.12 |units.kpc,
                 axis_ratio_bar= 0.37):
        
        # Simulation parameters
        self.t_end= simulation_time
        self.time= 0 |units.Myr
        self.dt_bridge= dt_bridge
        self.method= method_for_rotating_bridge
        
        #galaxy parameters
        self.omega_system= 0 | (units.kms/units.kpc)
        self.initial_phase_system= 0
        self.bar_phase= initial_phase_bar
        self.omega_bar= pattern_speed_bar
        self.mass_bar= mass_bar
        self.aaxis_bar= semimajor_axis_bar
        self.axis_ratio_bar= axis_ratio_bar
        
        self.spiral_phase_main_sp= initial_phase_main_spiral_arms
        self.omega_spiral_main_sp= pattern_speed_main_spiral_arms
        self.amplitude_main_sp= amplitude_main_spiral_arms
        self.number_main_sp= number_of_main_spiral_arms
        self.tangent_pitch_angle_main_sp= tangent_pitch_angle_main_spiral_arms
        
        self.spiral_phase_second_sp= initial_phase_secondary_spiral_arms
        self.omega_spiral_second_sp= pattern_speed_secondary_spiral_arms
        self.amplitude_second_sp= amplitude_secondary_spiral_arms
        self.number_second_sp= number_of_secondary_spiral_arms
        self.tangent_pitch_angle_second_sp= tangent_pitch_angle_secondary_spiral_arms
        self.separation_sp= separation_locus_spiral_arms
       
        
        return
        
    def softening(self, particles):
        '''
        optimum softening lenght.
        '''
        N= len(particles.mass)
        U= particles.potential_energy()
        Rvir= 0.5*constants.G*particles.mass.sum()**2/abs(U)	
        epsilon = 4*Rvir/N
        return epsilon 

    def galactic_model(self):
        '''
         Model of the Galaxy.
         In this example, the Galaxy has two-dimensional bar and spiral arms.
         The spiral arms are described by the Composite model (2+2)
         The bar does not grow adiabatically
         The axisymmetric component has its defaul values from Allen & Santillan (1990).
         '''
        galaxy= BarAndSpirals3D()
        galaxy.kinetic_energy=quantities.zero
        galaxy.potential_energy=quantities.zero
        galaxy.parameters.spiral_contribution= True
        galaxy.parameters.spiral_model= 2
        galaxy.parameters.omega_spiral= self.omega_spiral_main_sp
        galaxy.parameters.spiral_phase= self.spiral_phase_main_sp 
        galaxy.parameters.amplitude= self.amplitude_main_sp
        galaxy.parameters.m=  self.number_main_sp
        galaxy.parameters.tan_pitch_angle= self.tangent_pitch_angle_main_sp
        galaxy.parameters.phi21_spiral= self.spiral_phase_second_sp
        galaxy.parameters.omega_spiral2= self.omega_spiral_second_sp
        galaxy.parameters.amplitude2= self.amplitude_second_sp
        galaxy.parameters.m2= self.number_second_sp
        galaxy.parameters.tan_pitch_angle2= self.tangent_pitch_angle_second_sp
        galaxy.parameters.rsp= self.separation_sp  
        galaxy.parameters.bar_contribution= True
        galaxy.parameters.bar_phase= self.bar_phase
        galaxy.parameters.omega_bar= self.omega_bar
        galaxy.parameters.mass_bar= self.mass_bar
        galaxy.parameters.aaxis_bar= self.aaxis_bar
        galaxy.parameters.axis_ratio_bar= self.axis_ratio_bar  
        galaxy.commit_parameters()
        self.omega_system= galaxy.parameters.omega_system
        self.initial_phase_sytem= galaxy.parameters.initial_phase

        return galaxy

    def circular_velocity(self):
        MW= self.galactic_model()
        r= numpy.arange(15)
        vc= MW.get_velcirc(r|units.kpc, 0|units.kpc, 0|units.kpc)
        pyplot.plot(r, vc.value_in(units.kms))
        pyplot.show()
        
    def creation_cluster_in_rotating_frame(self, particles):
        "forming a cluster in a rotating frame"
        
        no_inertial_system= particles.copy()
        angle= self.initial_phase_system + self.omega_system*self.time
        C1= particles.vx + self.omega_system*particles.y
        C2= particles.vy - self.omega_system*particles.x
        no_inertial_system.x = particles.x*numpy.cos(angle) + particles.y*numpy.sin(angle)
        no_inertial_system.y = -particles.x*numpy.sin(angle) + particles.y*numpy.cos(angle) 
        no_inertial_system.z = particles.z
        no_inertial_system.vx = C1*numpy.cos(angle) + C2*numpy.sin(angle) 
        no_inertial_system.vy = C2*numpy.cos(angle) - C1*numpy.sin(angle)
        no_inertial_system.vz = particles.vz
        return no_inertial_system    
    
    def from_noinertial_to_cluster_in_inertial_frame(self, part_noin, part_in):
        'makes transformation to the inertial frame'
        
        angle= self.initial_phase_system + self.omega_system*self.time
        C1= part_noin.vx - part_noin.y*self.omega_system
        C2= part_noin.vy + part_noin.x*self.omega_system
        part_in.age= part_noin.age
        part_in.mass= part_noin.mass
        part_in.radius= part_noin.radius
        part_in.luminosity= part_noin.luminosity
        part_in.temperature= part_noin.temperature
        part_in.stellar_type=part_noin.stellar_type
        part_in.x= part_noin.x*numpy.cos(angle)-part_noin.y*numpy.sin(angle)
        part_in.y= part_noin.x*numpy.sin(angle)+part_noin.y*numpy.cos(angle)
        part_in.z= part_noin.z
        part_in.vx= C1*numpy.cos(angle) - C2*numpy.sin(angle)
        part_in.vy= C1*numpy.sin(angle) + C2*numpy.cos(angle)
        part_in.vz= part_noin.vz
        return
       
   
    def evolution_of_the_cluster(self, cluster):
        '''
        Function that makes de cluster evolution.
        input: cluster -> defined in an inertial frame (centered at the Galactic center)
        steps in this function:
        1. From cluster, another cluster is defined in a rotating frame
        2. The gravity code is initialized
        3. The stellar evolution is initialized
        4. the Galaxy model is constructed
        5. The Galaxy and the cluster in the rotating frame are coupled via the Rotating Bridge
        6. The evolution of the system is made
        7. the cluster properties are transformed back to the inertial frame 
        '''
        
        cluster_in_rotating_frame= self.creation_cluster_in_rotating_frame(cluster)
       
        # N body code
        epsilon = self.softening(cluster)
        convert_nbody= nbody_system.nbody_to_si(cluster.mass.sum(), cluster.virial_radius()) 
        
        gravity=Huayno(convert_nbody)
        gravity.parameters.timestep= self.dt_bridge/3.
        gravity.particles.add_particles(cluster_in_rotating_frame)
        gravity.parameters.epsilon_squared= epsilon**2
        channel_from_gravity_to_rotating_cluster= gravity.particles.new_channel_to(cluster_in_rotating_frame)
        channel_from_rotating_cluster_to_gravity= cluster_in_rotating_frame.new_channel_to(gravity.particles)
    
        #stellar evolution code
        se= SeBa()
        se.particles.add_particles(cluster_in_rotating_frame)
        channel_from_rotating_cluster_to_se= cluster_in_rotating_frame.new_channel_to(se.particles)
        channel_from_se_to_rotating_cluster= se.particles.new_channel_to(cluster_in_rotating_frame) 
        
        # Galaxy model and Rotating bridge
        MW= self.galactic_model()
        system= Rotating_Bridge(self.omega_system, timestep=self.dt_bridge, verbose=False, method= self.method) 
        system.add_system(gravity, (MW,), False)
        system.add_system(MW, (), False)  

        X=[]
        Y=[]
        T=[]
        #Cluster evolution
        while (self.time<= self.t_end-self.dt_bridge/2):
            self.time += self.dt_bridge
            
            system.evolve_model(self.time)
            se.evolve_model(self.time)
            
            channel_from_gravity_to_rotating_cluster.copy_attributes(['x', 'y', 'z', 'vx', 'vy', 'vz'])
            channel_from_se_to_rotating_cluster.copy_attributes(['mass', 'radius', 'luminosity', 'age', 'temperature', 'stellar_type'])
            channel_from_rotating_cluster_to_gravity.copy_attributes(['mass'])
            self.from_noinertial_to_cluster_in_inertial_frame(cluster_in_rotating_frame, cluster) 

            time= self.time.value_in(units.Myr)
            cm= cluster.center_of_mass()
           
            # write data
            if ((time==2) or (time==50) or (time==100) or (time==150)):
                X.append((cluster.x-cm[0]).value_in(units.kpc))  
                Y.append((cluster.y-cm[1]).value_in(units.kpc))
                T.append(time)  
               
                
        gravity.stop()
        se.stop()
        return T, X, Y

def plot_cluster(T, X,Y):
    
    figure= pyplot.figure(figsize=(10,10))
    figure.subplots_adjust(wspace=0.3)
    figure.suptitle('Evolution of a star cluster', fontsize=14, fontweight='bold')

    for i in range(0,len(T)): 
        ax= figure.add_subplot(2,2,i+1)
        ax.scatter(X[i], Y[i], label= 't= %g Myr'%T[i])
        ax.set_xlabel('X [kpc]')
        ax.set_ylabel('Y [kpc]')
        ax.set_xlim(-0.1, 0.1)
        ax.set_ylim(-0.1, 0.1)
        ax.legend(loc="upper right", ncol=1, shadow=False, fontsize=12)
    pyplot.show()

    
if __name__ in('__main__','__plot__'):

    
    #Create a cluster
    star_cluster= create_cluster_with_IMF(N= 200)
    star_cluster.position += [-6.5,0,0] | units.kpc
    star_cluster.velocity += [0, 50, 0] | units.kms
    
   
    # Construct the galactic model and make the evolution
    evolution= RealisticEvolutionCluster(simulation_time= 150|units.Myr)
    t, x, y= evolution.evolution_of_the_cluster(star_cluster)

    # plots
    plot_cluster(t, x, y)

   
