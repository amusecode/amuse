import numpy

from amuse.support.exceptions import AmuseException, AmuseWarning
from amuse.units import units, nbody_system, generic_unit_system, constants
from amuse.datamodel import Particles
from amuse.ext.sobol import i4_sobol_generate

class EnclosedMassInterpolator(object):
    """
    Interpolator used in 'get_enclosed_mass_from_tabulated' and 'get_radius_for_enclosed_mass'.
    These two functions are required for 'new_spherical_particle_distribution'.
    """
    def __init__(self, radii = None, densities = None, core_radius = None):
        self.initialized = False
        self.four_thirds_pi = numpy.pi * 4.0/3.0
        if (radii and densities):
            self.initialize(radii, densities, core_radius = core_radius)
    
    def initialize(self, radii, densities, core_radius = None):
        self.sort_density_and_radius(densities*1.0, radii*1.0, core_radius = core_radius)
        self.calculate_enclosed_mass_table()
        self.initialized = True
    
    def sort_density_and_radius(self, densities, radii, core_radius = None):
        self.radii, self.densities = radii.sorted_with(densities)
        self.radii.prepend(core_radius or 0 | units.m)
    
    def calculate_enclosed_mass_table(self):
        self.radii_cubed = self.radii**3
        self.enclosed_mass = [0.0] | units.kg
        for rho_shell, r3_in, r3_out in zip(self.densities, self.radii_cubed, self.radii_cubed[1:]):
            self.enclosed_mass.append(self.enclosed_mass[-1] + rho_shell * (r3_out - r3_in))
        self.enclosed_mass = self.four_thirds_pi * self.enclosed_mass
    
    def get_index(self, value, sorted_vector):
        out_of_bounds = numpy.logical_or(sorted_vector[0] > value, value > sorted_vector[-1])
        if out_of_bounds.any():
            value = numpy.compress(numpy.array([out_of_bounds]).flatten(), value.number) | value.unit
            raise AmuseException("Can't find a valid index. {0} is not in "
                "the range [{1}, {2}].".format(value, sorted_vector[0], sorted_vector[-1]))
        index = numpy.searchsorted(sorted_vector.number, value.value_in(sorted_vector.unit))
        return numpy.maximum(index - 1, 0)
    
    def get_enclosed_mass(self, radius):
        if not self.initialized:
            raise AmuseException("Can't calculate enclosed mass: interpolator is not initialized")
        index = self.get_index(radius, self.radii)
        return (self.enclosed_mass[index] + self.four_thirds_pi * 
            self.densities[index] * (radius**3 - self.radii_cubed[index]))
    
    def get_radius_for_enclosed_mass(self, enclosed_mass):
        if not self.initialized:
            raise AmuseException("Can't calculate radius for enclosed mass: interpolator is not initialized")
        index = self.get_index(enclosed_mass, self.enclosed_mass)
        return (((enclosed_mass - self.enclosed_mass[index]) / (self.four_thirds_pi * self.densities[index]) 
            + self.radii_cubed[index]))**(1.0/3.0)
    


class UniformSphericalDistribution(object):
    """
    Creates a uniform spherical grid of particles. Type can be:
    "cubic":  'crystal' composed of cubes with particles on each corner
    "bcc":    as cubic but with additional particles at the center of each cube
    "body_centered_cubic": same as "bcc"
    "fcc":    as cubic but with additional particles at the face of each cube
    "face_centered_cubic": same as "fcc"
    "random": particles are randomly distributed using numpy.random.uniform
    "glass":  like random, but stabilised using hydro pressure and no gravity
    "sobol":  3D sobol sequence (low discrepancy, quasi-random)
    
    "offset" is only used for the regular grids ("cubic", "bcc", "fcc"), and 
    should contain three numbers in the half-open interval [0, 1). These 
    define the offset between the origin of the grid and the corner 
    of the unit cell, normalized to the unit cell size.
    "target_rms" is only used for "glass" as the density criterion for convergence
    """
    def __init__(self, number_of_particles, type = "bcc", 
            offset = (0.82832951,  0.27237167,  0.37096327), 
            mass_cutoff = 1, target_rms = 0.01):
        if not hasattr(self, type):
            raise TypeError("Unknown grid type option: {0}".format(type))
        self.number_of_particles = number_of_particles
        self.type = type
        self.offset = offset
        self.mass_cutoff = mass_cutoff
        self.target_rms = target_rms
    
    def cubic(self):
        n1D = numpy.ceil( 0.5*(2*self.number_of_particles)**(1./3) ) * 2 + 3 # odd number
        x, y, z = numpy.mgrid[-1. : 1. : n1D*1j,
                              -1. : 1. : n1D*1j,
                              -1. : 1. : n1D*1j]
        x = x.flatten()
        y = y.flatten()
        z = z.flatten()
        for delta, vec in zip(self.offset, (x,y,z)):
            vec += delta * 2.0 / (n1D - 1)
        return self._cutout_sphere(x, y, z)
    
    def bcc(self):
        n1D = numpy.ceil( 0.5*(self.number_of_particles)**(1./3) ) * 2 + 3 # odd number
        x1,y1,z1 = numpy.mgrid[-1. : 1. : n1D*1j,
                               -1. : 1. : n1D*1j,
                               -1. : 1. : n1D*1j]
        n_2 = n1D - 1
        x2,y2,z2 = numpy.mgrid[-1.+1./n_2 : 1.-1./n_2 : n_2*1j,
                               -1.+1./n_2 : 1.-1./n_2 : n_2*1j,
                               -1.+1./n_2 : 1.-1./n_2 : n_2*1j]
        x = numpy.concatenate( (x1.flatten(),x2.flatten()) )
        y = numpy.concatenate( (y1.flatten(),y2.flatten()) )
        z = numpy.concatenate( (z1.flatten(),z2.flatten()) )
        for delta, vec in zip(self.offset, (x,y,z)):
            vec += delta * 2.0 / n_2
        return self._cutout_sphere(x, y, z)
    
    body_centered_cubic = bcc
    
    def fcc(self):
        n1D = numpy.ceil((self.number_of_particles / 2.0)**(1.0/3.0)) + 1
        delta = 1.0 / (n1D - 1.5)
        x, y, z = numpy.mgrid[-1.0-2*delta : 1-delta : n1D*1j,
                              -1.0-2*delta : 1-delta : n1D*1j,
                              -1.0-2*delta : 1-delta : n1D*1j]
        x0 = x.flatten() + self.offset[0] * 2 * delta
        y0 = y.flatten() + self.offset[1] * 2 * delta
        z0 = z.flatten() + self.offset[2] * 2 * delta
        x1 = x0 + delta
        y1 = y0 + delta
        z1 = z0 + delta
        x = numpy.concatenate((x0, x1, x0, x1))
        y = numpy.concatenate((y0, y1, y1, y0))
        z = numpy.concatenate((z0, z0, z1, z1))
        return self._cutout_sphere(x, y, z)
    
    face_centered_cubic = fcc
    
    def _random_cube(self, number_of_particles):
        x = numpy.random.uniform(-1.0, 1.0, number_of_particles)
        y = numpy.random.uniform(-1.0, 1.0, number_of_particles)
        z = numpy.random.uniform(-1.0, 1.0, number_of_particles)
        return x, y, z
    
    def random(self, number_of_particles = None, try_number_of_particles = None):
        if number_of_particles is None:
            number_of_particles = self.number_of_particles
        if try_number_of_particles is None:
            try_number_of_particles = number_of_particles
        try_number_of_particles = int(try_number_of_particles)
        x, y, z = self._random_cube(2*try_number_of_particles)
        r_squared = x*x + y*y + z*z
        select_sphere = numpy.where( r_squared < self.mass_cutoff**(2.0/3.0))
        if len(select_sphere[0]) < number_of_particles:
            return self.random(number_of_particles, numpy.ceil(try_number_of_particles*1.1) )
        else:
            return (x[select_sphere][0:number_of_particles], 
                y[select_sphere][0:number_of_particles], 
                z[select_sphere][0:number_of_particles])
    

    def sobol(self):
        x, y, z = i4_sobol_generate(3, 2*self.number_of_particles, 3) * 2.0 - 1.0
        return self._cutout_sphere(x, y, z)
    
    def glass(self):
        from amuse.community.fi.interface import Fi
        
        if self.target_rms < 0.0001:
            print("warning: target_rms may not succeed")
        if self.number_of_particles < 1000:
            print("warning: not enough particles")
        
        N = 2 * self.number_of_particles
        L = 1 | nbody_system.length
        dt = 0.01 | nbody_system.time
        
        x, y, z = self._random_cube(N)
        vx,vy,vz= self.random(N)
        
        p = Particles(N)
        p.x = L * x
        p.y = L * y
        p.z = L * z
        p.h_smooth = 0.0 | nbody_system.length
        p.vx = (0.1 | nbody_system.speed) * vx[:N]
        p.vy = (0.1 | nbody_system.speed) * vy[:N]
        p.vz = (0.1 | nbody_system.speed) * vz[:N]
        p.u = (0.1*0.1) | nbody_system.speed**2
        p.mass = (8.0/N) | nbody_system.mass
        
        sph = Fi(mode = 'periodic', redirection = 'none')
        sph.initialize_code()
        
        sph.parameters.use_hydro_flag = True
        sph.parameters.radiation_flag = False
        sph.parameters.self_gravity_flag = False
        sph.parameters.gamma = 1.0
        sph.parameters.isothermal_flag = True
        sph.parameters.integrate_entropy_flag = False
        sph.parameters.timestep = dt
        sph.parameters.verbosity = 0
        sph.parameters.periodic_box_size = 2 * L
        sph.parameters.artificial_viscosity_alpha = 1.0
        sph.parameters.beta = 2.0
        sph.commit_parameters()
        sph.gas_particles.add_particles(p)
        sph.commit_particles()
        
        t = 0.0 | nbody_system.time
        rms = 1.0
        minrms = 1.0
        i = 0
        while rms > self.target_rms:
            i += 1
            t += (0.25 | nbody_system.time)
            sph.evolve_model(t)
            rho = sph.particles.rho.value_in(nbody_system.density)
            rms = rho.std()/rho.mean()
            minrms = min(minrms, rms)
            if (rms > 2.0*minrms) or (i > 300):
                print(" RMS(rho) convergence warning:", i, rms, minrms)
            if i > 100000:
                print("i> 100k steps - not sure about this...")
                print(" rms:", rms)
                break
        
        x = sph.particles.x.value_in(nbody_system.length)
        y = sph.particles.y.value_in(nbody_system.length)
        z = sph.particles.z.value_in(nbody_system.length)
        sph.stop()
        del sph
        return self._cutout_sphere(x, y, z)
    

    def _cutout_sphere(self, x, y, z):
        r_squared = x*x + y*y + z*z
        sorted_indices = numpy.argsort(r_squared)
        massfrac_edge = r_squared[sorted_indices[self.number_of_particles-1]]**(1.5)
        massfrac_next = r_squared[sorted_indices[self.number_of_particles]]**(1.5)
        r_max = (0.5*(massfrac_edge+massfrac_next) / self.mass_cutoff)**(1.0/3.0)
        indices = sorted_indices[:self.number_of_particles]
        return x[indices]/r_max, y[indices]/r_max, z[indices]/r_max
    
    @property
    def result(self):
        return getattr(self, self.type)()
    

keyword_arguments_doc = """    :argument keyword_arguments:    Optional arguments to UniformSphericalDistribution:
        :argument type:             Type of the basegrid. Can be:
            "cubic":  'crystal' composed of cubes with particles on each corner
            "bcc":    as cubic but with additional particles at the center of each cube
            "body_centered_cubic": same as "bcc"
            "fcc":    as cubic but with additional particles at the face of each cube
            "face_centered_cubic": same as "fcc"
            "random": particles are randomly distributed using numpy.random.uniform
            "glass":  like random, but stabilised using hydro pressure and no gravity
            "sobol":  3D sobol sequence (low discrepancy, quasi-random)
        :argument offset:           only used for the regular grids ("cubic", "bcc", "fcc"), and 
            should contain three numbers in the half-open interval [0, 1). These 
            define the offset between the origin of the grid and the corner 
            of the unit cell, normalized to the unit cell size.
        :argument target_rms        only used for "glass" as the density criterion for convergence
"""

def new_uniform_spherical_particle_distribution(number_of_particles, size, total_mass, **keyword_arguments):
    """
    Returns a Particles set with positions following a uniform 
    spherical distribution. Only the positions and masses 
    (equal-mass system) are set.
    
    :argument number_of_particles:  Number of particles in the resulting model
    :argument size:                 Radius of the sphere enclosing the model
    :argument total_mass:           Total mass of the Particles set
    """
    particles = Particles(number_of_particles)
    particles.mass = total_mass * 1.0 / number_of_particles
    x, y, z = UniformSphericalDistribution(number_of_particles, **keyword_arguments).result
    particles.x = size * x
    particles.y = size * y
    particles.z = size * z
    return particles
new_uniform_spherical_particle_distribution.__doc__ += keyword_arguments_doc

def new_spherical_particle_distribution(number_of_particles, 
        radial_density_func = None,     # not yet supported, specify radii and densities tables:
        radii = None, densities = None, 
        total_mass = None, size = None, # if total_mass is not given, it will be deduced from size or max(radii)
        **keyword_arguments):           # optional arguments for UniformSphericalDistribution
    """
    Returns a Particles set with positions following a spherical 
    distribution. The radial density profile is determined from the 
    look-up table (radii, densities). Entries in the 'radii' table 
    are interpreted as the outer radius of the shell, with uniform 
    density as defined by the corresponding entry in the 'densities' 
    table:
    rho(r) = densities[i],  for ( radii[i-1] <= r <= radii[i] )
    
    Only the positions and masses (equal-mass system) are set.
    
    :argument number_of_particles:  Number of particles in the resulting model
    :argument radii:                Table with radii for the radial density profile
    :argument densities:            Table with densities for the radial density profile
    :argument total_mass:           Total mass of the Particles set (optional, will be 
                                    deduced from size or max(radii) otherwise)
    :argument size:                 Radius of the sphere enclosing the model (optional)
    """
    if (radii is None) or (densities is None):
        raise AmuseException("Using an arbitrary radial density function is not yet "
            "supported. Radius and density tables must be passed instead.")
    
    interpolator = EnclosedMassInterpolator()
    interpolator.initialize(radii, densities)
    if total_mass is None:
        total_mass = interpolator.get_enclosed_mass(size or max(radii))
    particles = Particles(number_of_particles)
    particle_mass = total_mass * 1.0 / number_of_particles
    particles.mass = particle_mass
    x, y, z = UniformSphericalDistribution(number_of_particles, **keyword_arguments).result
    # Now scale the uniformly distributed particle positions to match the radial density profile
    r_old = numpy.sqrt(x*x + y*y + z*z)
    indices = numpy.argsort(r_old)
    if r_old[indices[0]] == 0.0:
        r_old[indices[0]] = 1.0
    f_scale = interpolator.get_radius_for_enclosed_mass(
        (numpy.arange(0.5, number_of_particles + 0.5) | units.none) * particle_mass) / r_old[indices]
    particles.x = (f_scale * x[indices]).as_quantity_in(radii.unit)
    particles.y = (f_scale * y[indices]).as_quantity_in(radii.unit)
    particles.z = (f_scale * z[indices]).as_quantity_in(radii.unit)
    return particles
new_spherical_particle_distribution.__doc__ += keyword_arguments_doc

plummer_arguments_doc = """
    :argument number_of_particles:  Number of particles in the resulting model
    :argument total_mass:           Total mass of the Particles set 
    :argument virial_radius:        Virial radius of the Plummer model
""" + keyword_arguments_doc

def new_plummer_spatial_distribution(number_of_particles, 
        total_mass = 1.0|nbody_system.mass, 
        virial_radius = 1.0|nbody_system.length,
        mass_cutoff = 0.999,
        **keyword_arguments):           # optional arguments for UniformSphericalDistribution
    """
    Returns a Particles set with positions following a Plummer 
    distribution. 
    Only the positions and masses (equal-mass system) are set.
    """
    particles = Particles(number_of_particles)
    particle_mass = total_mass * 1.0 / number_of_particles
    particles.mass = particle_mass
    x, y, z = UniformSphericalDistribution(
        number_of_particles, mass_cutoff=mass_cutoff, **keyword_arguments).result
    
    # Now scale the uniformly distributed particle positions to match the radial density profile
    r_old = numpy.sqrt(x*x + y*y + z*z)
    scale_factor = (0.1875 * numpy.pi * virial_radius.number) / numpy.sqrt(1.0 - r_old**2)
    particles.x = scale_factor * x | virial_radius.unit
    particles.y = scale_factor * y | virial_radius.unit
    particles.z = scale_factor * z | virial_radius.unit
    return particles
new_plummer_spatial_distribution.__doc__ += plummer_arguments_doc

def new_gas_plummer_distribution(number_of_particles, 
        total_mass = 1.0|nbody_system.mass, 
        virial_radius = 1.0|nbody_system.length,
        G = None,
        **keyword_arguments):           # optional arguments for UniformSphericalDistribution
    """
    Create a plummer gas model with the given number of particles. Returns a set 
    of SPH particles with equal masses and positions distributed to fit a plummer 
    distribution model. Velocities are set to zero, and internal energies are set 
    to balance the gravitational forces between the gas particles.
    """
    particles = new_plummer_spatial_distribution(number_of_particles, total_mass=total_mass, 
        virial_radius=virial_radius, **keyword_arguments)
    
    if G is None:
        G = nbody_system.G if generic_unit_system.is_generic_unit(total_mass.unit) else constants.G
    velocity_unit = (G*total_mass/virial_radius).sqrt().unit.base_unit()
    particles.velocity = [0.0, 0.0, 0.0] | velocity_unit
    
    plummer_radius = 0.1875 * numpy.pi * virial_radius
    u_unit = (velocity_unit**2).base_unit()
    particles.u = (1 + particles.position.lengths_squared()/plummer_radius**2)**(-0.5) | u_unit
    particles.u *= 0.25 * (G*total_mass**2/virial_radius) / particles.thermal_energy()
    return particles
new_gas_plummer_distribution.__doc__ += plummer_arguments_doc

def sample_from_velocity_distribution(number_of_particles):
    x = numpy.random.uniform(0.0, 1.0, number_of_particles)
    y = numpy.random.uniform(0.0, 0.1, number_of_particles)
    selected = x[y <= (x**2) * (1.0 - x**2)**3.5]
    if len(selected) < number_of_particles:
        return numpy.concatenate((selected, sample_from_velocity_distribution(number_of_particles-len(selected))))
    return selected

def random_direction(number_of_particles):
    z = numpy.random.uniform(-1.0, 1.0, number_of_particles)
    sine_theta = numpy.sqrt(1-z*z)
    phi = numpy.random.uniform(0.0, 2*numpy.pi, number_of_particles)
    return numpy.array([sine_theta * numpy.cos(phi), sine_theta * numpy.sin(phi), z]).transpose()

def new_plummer_distribution(number_of_particles, 
        total_mass = 1.0|nbody_system.mass, 
        virial_radius = 1.0|nbody_system.length,
        mass_cutoff = 0.999,
        G = None,
        **keyword_arguments):           # optional arguments for UniformSphericalDistribution
    """
    Create a plummer model with the given number of particles. Returns a set 
    of particles with equal masses and positions distributed to fit a plummer 
    distribution model. Velocities are sampled using von Neumann's rejection 
    method (Aarseth et al. 1974), balancing the gravitational forces between 
    the particles.
    """
    particles = new_plummer_spatial_distribution(number_of_particles, total_mass=total_mass, 
        virial_radius=virial_radius, **keyword_arguments)
    
    if G is None:
        G = nbody_system.G if generic_unit_system.is_generic_unit(total_mass.unit) else constants.G
    velocity_unit = (G*total_mass/virial_radius).sqrt().unit.base_unit()
    plummer_radius = 0.1875 * numpy.pi * virial_radius
    
    escape_velocity = (1 + particles.position.lengths_squared()/plummer_radius**2)**(-0.25) | velocity_unit
    velocity = escape_velocity * sample_from_velocity_distribution(number_of_particles)
    velocity *= numpy.sqrt((G*total_mass*number_of_particles) / (2*virial_radius*velocity.length_squared()))
    particles.velocity = velocity.reshape((-1,1)) * random_direction(number_of_particles)
    return particles
new_plummer_distribution.__doc__ += plummer_arguments_doc
