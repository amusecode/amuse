from amuse.lab import *
from amuse.community.gadget2.interface import Gadget2 as Gadget2
from amuse.units import units, constants, nbody_system
from amuse.ext.sink import new_sink_particles
from amuse.ext.protodisk import ProtoPlanetaryDisk
from amuse.datamodel import Particles
import numpy
import time as timing

class MinimalWorkingExample(object):
    def __init__(self,     
                 total_N=16000, #total number of disc particles
                 tend=250. | units.yr, #End time of the simulation
                 Mstar=1. | units.MSun, #Mass of the accreting star
                 dt = 1. | units.yr, #timestep
                 temp = 25. | units.K, #temperature of the ISM
                 v_ism = 3.0 | units.kms, #Velocity of the ISM with respect to the star
                 n_core = 4, #Number of cores used by the community codes
                 mu = 2.3, #mean molecular weight
                 theta = 0., #angle of the disc with respect to the flow, 0 means disc rotation axis is parallel with flow direction
                 n_dens = 5e6 | (units.cm)**(-3.), #Number density of the ism
                 dirname = './'):
        
        self.disc_N = total_N 
        self.tend = tend
        self.Mstar = Mstar
        self.dt = dt
        self.mu = mu * constants.proton_mass
        self.T = temp 
        self.cs = ((temp*constants.kB)/self.mu).sqrt()
#        print "cs=", self.cs.in_(units.kms)
        self.v_ism = v_ism
        self.ism_n_dens = n_dens #number density of the ism
        self.ism_dens = self.mu*self.ism_n_dens #mass density of the ism
        self.converter = nbody_system.nbody_to_si(self.Mstar, 1. | units.AU)
       
        self.disc_angle = numpy.radians(theta)
        self.n_core = n_core 
        self.dirname = dirname
        self.discfraction = 0.01
        self.disc_Rmin = self.converter.to_generic(10.| units.AU).value_in(nbody_system.length)
        self.disc_Rmax = self.converter.to_generic(100.| units.AU).value_in(nbody_system.length)
        self.filename = "DiskWind.h5"
    
        print 'Files will be saved in ', self.dirname


    def initialize_star(self):
        
        self.star=Particles(1)
        self.star.mass=self.Mstar
        self.star.radius= 1. | units.AU
        self.star.x=0.|units.AU
        self.star.y=0.|units.AU
        self.star.z=0.|units.AU
        self.star.vx=0.|units.kms
        self.star.vy=0.|units.kms
        self.star.vz=0.|units.kms


    def initialize_data(self):
        

        self.cylinder_radius = 500. | units.AU 
        self.cylinder_length = 2.*self.cylinder_radius
        self.cylinder_vol = self.cylinder_length*numpy.pi*self.cylinder_radius**2.
        self.sph_particle_dens = self.ism_dens/(0.01*self.Mstar/self.disc_N)
        self.initialize_star()
        self.initialize_ism()


    def initialize_ism(self):

        #specific internal energy of the ism, i.e. internal energy per unit mass
        #In case of an isothermal EOS, Gadget requires the sound-speed squared as input parameter instead of the thermal energy per unit mass
        self.ism_u = self.cs**2.
        
        self.ism_slice_length = self.v_ism*self.dt
        self.ism_slice_vol = self.ism_slice_length*numpy.pi*self.cylinder_radius**2.0
    
        self.ism_slice_mass = self.ism_dens*self.ism_slice_vol
        self.ism_mass = self.ism_dens*self.cylinder_vol
        
        #The x-coordinate of the inflow
        self.offset = self.cylinder_radius 
        
        self.ism_slice_N = int(round(self.sph_particle_dens*self.ism_slice_vol))    
        self.total_N = self.ism_slice_N*(self.cylinder_length/self.ism_slice_length)
        

    def make_slice(self):
    
        new_slice = Particles(self.ism_slice_N)
    
        rho = numpy.sqrt(numpy.random.uniform(0,1,self.ism_slice_N))*self.cylinder_radius.value_in(units.AU)
        phi = numpy.random.uniform(0,2.0*numpy.pi, self.ism_slice_N)
    

        new_slice.x = numpy.random.uniform(self.ism_slice_length.value_in(units.AU), 0, self.ism_slice_N) - self.offset.value_in(units.AU) | units.AU
        new_slice.y = rho*numpy.sin(phi) | units.AU
        new_slice.z = rho*numpy.cos(phi) | units.AU
    
        new_slice.vx = 1000*self.v_ism
        new_slice.vy = 0.0 | units.kms
        new_slice.vz = 0.0 | units.kms    
    
        new_slice.mass = self.ism_slice_mass/self.ism_slice_N
        new_slice.u = self.ism_u

        return new_slice
    

    def create_disc(self):
                
        #The following line makes sure masses for ISM and disc particles are equal:
        self.discfraction = (self.disc_N*(self.ism_slice_mass/self.ism_slice_N))/self.Mstar
        
        T_disc = self.T
        cs_disc = ((T_disc*constants.kB)/self.mu).sqrt()
        densitypower = 1.5
        g2=2-densitypower
        k_out=((1+self.discfraction)/self.disc_Rmax**3)**0.5
        sigma_out=g2*self.discfraction/(2*numpy.pi*self.disc_Rmax**densitypower*(self.disc_Rmax**g2-self.disc_Rmin**g2))
        q_out = self.converter.to_generic(cs_disc).value_in(nbody_system.length/nbody_system.time)/(numpy.pi*sigma_out/k_out)
        
        print "Number of disk particles:", self.disc_N
        
        proto=ProtoPlanetaryDisk(self.disc_N,
                                 convert_nbody=self.converter,
                                 discfraction=self.discfraction,
                                 densitypower=1.5,
                                 thermalpower=0,
                                 Rmin=self.disc_Rmin,
                                 Rmax=self.disc_Rmax,
                                 q_out=q_out)
        disc=proto.result
    
        print "The mass of a disc particle = ", disc.mass[0].value_in(units.kg)
    
        #Rotate 90 degrees with respect to the z-axis and then theta degrees with respect to the y-axis
        
        temp_x = disc[:].x
        temp_y = disc[:].y
        temp_z = disc[:].z
        temp_vx = disc[:].vx
        temp_vy = disc[:].vy
        temp_vz = disc[:].vz
    
    
        disc.x = temp_z*numpy.cos(self.disc_angle) - temp_y*numpy.sin(self.disc_angle)
        disc.y = temp_z*numpy.sin(self.disc_angle) + temp_y*numpy.cos(self.disc_angle)
        disc.z = -temp_x
    
        disc.vx = temp_vz*numpy.cos(self.disc_angle) - temp_vy*numpy.sin(self.disc_angle)
        disc.vy = temp_vz*numpy.sin(self.disc_angle) + temp_vy*numpy.cos(self.disc_angle)
        disc.vz = -temp_vx
             
        return disc

 
    def evolve_model(self):
        self.initialize_data()

        self.ism_code = Gadget2(self.converter, number_of_workers=self.n_core)#, debugger='gdb')
        self.ism_code.parameters.time_max = 1024*self.dt       
        self.ism_code.parameters.n_smooth = 64 
        self.ism_code.parameters.n_smooth_tol = 2./64.   
        self.ism_code.parameters.artificial_viscosity_alpha = 0.1 
        self.ism_code.parameters.epsilon_squared = (1. | units.AU)**2.

        self.all_particles = Particles()
        write_set_to_file(self.star, self.filename, "hdf5", append_to_file=False)
        write_set_to_file(self.all_particles, self.filename, "hdf5", append_to_file=False)
        
        self.initial_disc_particles = self.create_disc()
        self.all_particles.add_particles(self.initial_disc_particles)
        self.ism_code.gas_particles.add_particles(self.initial_disc_particles)
        
        #You can only add a sink after adding gas particles
        #starinsph refers to the corresponding particle set/id in the community code
        starinsph = self.ism_code.dm_particles.add_particles(self.star)
        #Use the build-in sink particle routine from amuse.ext.sink. Sink_radius needs to be defined manually otherwise the particle radius in gadget is taken,
        #which does not corresponding to the particle radius in the framework (since 'particles' in gadget do not have a radius, it is set to 0.01 | generic_unit_system.length
        #and corresponds to the gas gravitational smoothing epsilon.           
        sink = new_sink_particles(starinsph, sink_radius= self.star.radius)
       
        self.channel_from_ismcode_to_framework = self.ism_code.gas_particles.new_channel_to(self.all_particles)
    
        time = 0. | units.yr
    
        while time <= (self.tend+self.dt/2.):
 
 
            print "Adding new slice of ISM..."
            newslice=self.make_slice()
            self.ism_code.gas_particles.add_particles(newslice)
            self.all_particles.add_particles(newslice)

                
            start = timing.time()
            print "======================================================="
            print "Evolving to time = ", time.value_in(units.yr), " of ", self.tend.value_in(units.yr)," years..."    
            self.ism_code.evolve_model(time)
            print "This took ", (timing.time() - start), " s"
            
            out_of_bounds = self.ism_code.gas_particles.select_array(lambda x,y,z:(x > self.cylinder_radius)|((z**2+y**2).sqrt() >= self.cylinder_radius), ["x","y","z"]) 
            if len(out_of_bounds)>0:
                print "Removing ", len(out_of_bounds), " particles from the code because they were out of bounds"
                
                self.ism_code.gas_particles.remove_particles(out_of_bounds)
            
            
            self.ism_code.gas_particles.synchronize_to(self.all_particles)  
            sink.accrete(self.ism_code.gas_particles)

            write_set_to_file(self.star, self.filename, "hdf5")
            write_set_to_file(self.all_particles, self.filename, "hdf5")
            
            time += self.dt
            
        print "======================================================="
        self.ism_code.stop()


def new_option_parser():
    from amuse.units.optparse import OptionParser
    result = OptionParser()
    result.add_option("-n", 
                      dest="total_N", 
                      type="int", 
                      default = 1000,
                      help="Total number of disc particles")
    result.add_option("-t", 
                      dest="t_end", 
                      type="float", 
                      default = 100,
                      unit = units.yr, 
                      help="end time [%unit]")
    result.add_option("-d", 
                      dest="dt_diag", 
                      type="float", default = 1,
                      unit = units.yr, 
                      help="diagnostic timestep [%unit]")
    result.add_option("-m", 
                      dest="mass", 
                      type="float", 
                      default = 1.0,
                      unit = units.MSun, 
                      help="mass of the star [%unit]")
    result.add_option("--temp", 
                      dest="temp", 
                      type="float", 
                      default = 25.,
                      unit = units.K, 
                      help="Temperature of the ISM [%unit]")
    result.add_option("-v", 
                      dest="v_ism", 
                      type="float", 
                      default = 3.,
                      unit = units.kms, 
                      help="Velocity of the ISM [%unit]")
    result.add_option("--mu", 
                      dest="m_ism", 
                      type="float",
                      default = 2.3,
                      help="Mean molecular weight of ISM [%unit]")
    result.add_option("-i", 
                      dest="inclination", 
                      type="float", 
                      default = 0.0,
                      help="inclination of the disc with the y-axis")
    result.add_option("--ndens", 
                      dest="n_dens", 
                      type="int", 
                      default = 5e6 | (units.cm)**(-3.),
                      unit = (units.cm)**(-3.),
                      help="Number density of the ISM")
    result.add_option("--dir", 
                      dest="directory", 
                      type="string", 
                      default = '/home/thomas/Documents/simulations/disc',
                      help="The directory where to save the plots and data")
    result.add_option("-c", 
                      dest="n_cores", 
                      type="int", 
                      default = 4,
                      help="Number of cores used by community codes")
    return result





if __name__ in ("__main__","__plot__"):
    
    o, arguments  = new_option_parser().parse_args()

        
    code = MinimalWorkingExample( total_N=o.total_N,
                       tend=o.t_end,
                       Mstar=o.mass,
                       dt = o.dt_diag,
                       temp = o.temp,
                       v_ism = o.v_ism,
                       mu = o.m_ism,
                       n_dens = o.n_dens,
                       dirname = o.directory,
                       n_core = o.n_cores,
                       theta = o.inclination)

    code.evolve_model()
    

    
