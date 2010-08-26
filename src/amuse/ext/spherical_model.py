import numpy
from amuse.support.data.core import Particles
from amuse.support.units import units
from amuse.support.exceptions import AmuseWarning, AmuseException


class EnclosedMassInterpolator(object):
    """
    Interpolator used in 'get_enclosed_mass_from_tabulated' and 'get_radius_for_enclosed_mass'.
    These two functions are required for 'new_spherical_particle_distribution'.
    """
    def __init__(self):
        self.initialized = False
        self.four_thirds_pi = numpy.pi * 4.0/3.0
    
    def initialize(self, radii, densities):
        self.radius_unit =  radii.unit
        self.mass_unit = (densities.unit * radii.unit**3).to_simple_form()
        self.sort_density_and_radius(densities.number, radii.number)
        self.calculate_enclosed_mass_table()
        self.initialized = True
    
    def sort_density_and_radius(self, densities, radii):
        dtype = [('density', 'float64'), ('radius', 'float64')]
        sorted = numpy.sort(numpy.array(zip(densities, radii), dtype=dtype), order='radius')
        self.densities = sorted['density']
        self.radii     = numpy.concatenate( ([0], sorted['radius']) )
    
    def calculate_enclosed_mass_table(self):
        self.radii_cubed = self.radii**3
        self.enclosed_mass = [0.0]
        for rho_shell, r3_in, r3_out in zip(self.densities, self.radii_cubed, self.radii_cubed[1:]):
            self.enclosed_mass.append(self.enclosed_mass[-1] + rho_shell * (r3_out - r3_in))
        self.enclosed_mass = self.four_thirds_pi * numpy.array(self.enclosed_mass)
    
    def get_index(self, radius):
        if not self.radii[0] <= radius <= self.radii[-1]:
            raise AmuseException("Extrapolation is not supported. {0} is not in "
                "the range [{1}, {2}].".format(radius|self.radius_unit, 
                self.radii[0]|self.radius_unit, self.radii[-1]|self.radius_unit))
        index = 1
        while self.radii[index] < radius:
            index += 1
        return index - 1
    
    def get_enclosed_mass(self, radius):
        radius = radius.value_in(self.radius_unit)
        index = self.get_index(radius)
        return (self.enclosed_mass[index] + self.four_thirds_pi * 
            self.densities[index] * (radius**3 - self.radii_cubed[index])) | self.mass_unit
    
    def get_index2(self, mass):
        if not self.enclosed_mass[0] <= mass <= self.enclosed_mass[-1]:
            raise AmuseException("No valid radius for this enclosed mass. {0} is not in "
                "the range [{1}, {2}].".format(mass|self.mass_unit, 
                self.enclosed_mass[0]|self.mass_unit, self.enclosed_mass[-1]|self.mass_unit))
        index = 1
        while self.enclosed_mass[index] < mass:
            index += 1
        return index - 1
    
    def get_radius_for_mass(self, mass):
        mass = mass.value_in(self.mass_unit)
        index = self.get_index2(mass)
        return (((mass - self.enclosed_mass[index]) / (self.four_thirds_pi * self.densities[index]) 
            + self.radii_cubed[index]))**(1.0/3.0) | self.radius_unit
    

def get_enclosed_mass_from_tabulated(radius, radii = None, densities = None, interpolator = EnclosedMassInterpolator()):
    if not (radii is densities is None):
        interpolator.initialize(radii, densities)
    if not interpolator.initialized:
        raise AmuseException("Interpolator is not initialized. Radius and density tables must be passed in the first call.")
    return interpolator.get_enclosed_mass(radius)
    
def get_radius_for_enclosed_mass(mass, radii = None, densities = None, interpolator = EnclosedMassInterpolator()):
    if not (radii is densities is None):
        interpolator.initialize(radii, densities)
    if not interpolator.initialized:
        raise AmuseException("Interpolator is not initialized. Radius and density tables must be passed in the first call.")
    return interpolator.get_radius_for_mass(mass)
    




class UniformSphericalDistribution(object):
    """
    Creates a uniform spherical grid of particles. Type can be:
    "cubic":  'crystal' composed of cubes with particles on each corner
    "bcc":    as cubic but with additional particles at the center of each cube
    "body_centered_cubic": same as "bcc"
    "random": particles are randomly distributed using numpy.random.uniform
    The offset is only used for the regular grids ("cubic", "bcc"), and 
    should contain three numbers in the half-open interval [0, 1). These 
    define the offset between the origin of the grid and the corner 
    of the unit cell, normalized to the unit cell size.
    The result can optionally be rotated around rotate = (theta, phi).
    """
    def __init__(self, number_of_particles, type = "bcc", 
            offset = (0.82832951,  0.27237167,  0.37096327), 
            rotate = None, seed = None):
        if not hasattr(self, type):
            raise TypeError("Unknown grid type option: {0}".format(type))
        self.N = number_of_particles
        self.type = type
        self.offset = offset
        self.rotate = rotate
        if rotate:
            raise AmuseWarning("rotate is not yet supported")
        self.seed = seed
    
    def cubic(self, number_of_particles):
        n1D = numpy.ceil( 0.5*(2*number_of_particles)**(1./3) ) * 2 + 3 # odd number
        x, y, z = numpy.mgrid[-1. : 1. : n1D*1j,
                              -1. : 1. : n1D*1j,
                              -1. : 1. : n1D*1j]
        x = x.flatten()
        y = y.flatten()
        z = z.flatten()
        for delta, vec in zip(self.offset, (x,y,z)):
            vec += delta * 2.0 / (n1D - 1)
        r_squared = x*x + y*y + z*z
        r_squared_max = numpy.sort(r_squared)[self.N]
        r_max = numpy.sqrt(r_squared_max)
        select_sphere = numpy.where( r_squared < r_squared_max )
        return x[select_sphere]/r_max, y[select_sphere]/r_max, z[select_sphere]/r_max
    
    def bcc(self, number_of_particles):
        n1D = numpy.ceil( 0.5*(number_of_particles)**(1./3) ) * 2 + 3 # odd number
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
        r_squared = x*x + y*y + z*z
        r_squared_max = numpy.sort(r_squared)[self.N]
        r_max = numpy.sqrt(r_squared_max)
        select_sphere = numpy.where( r_squared < r_squared_max )
        return x[select_sphere]/r_max, y[select_sphere]/r_max, z[select_sphere]/r_max
    
    body_centered_cubic = bcc
    
    def random(self, number_of_particles):
        numpy.random.seed(self.seed)
        x = numpy.random.uniform(-1., 1., 2*number_of_particles)
        y = numpy.random.uniform(-1., 1., 2*number_of_particles)
        z = numpy.random.uniform(-1., 1., 2*number_of_particles)
        r_squared = x*x + y*y + z*z
        select_sphere = numpy.where( r_squared < 1.)
        if len(select_sphere[0]) < self.N:
            return self.random( numpy.ceil(number_of_particles*1.1) )
        else:
            return x[select_sphere][0:self.N], y[select_sphere][0:self.N], z[select_sphere][0:self.N]
    
    @property
    def result(self):
        return getattr(self, self.type)(self.N)
    

"""
Returns a Particles set with positions following a uniform 
spherical distribution. Only the positions and masses 
(equal-mass system) are set.

:argument number_of_particles:  Number of particles in the resulting model
:argument size:                 Radius of the sphere enclosing the model
:argument total_mass:           Total mass of the Particles set
:argument keyword_arguments:    Optional arguments to UniformSphericalDistribution
"""
def new_uniform_spherical_particle_distribution(number_of_particles, size, total_mass, **keyword_arguments):
    particles = Particles(number_of_particles)
    particles.mass = total_mass * 1.0 / number_of_particles
    x, y, z = UniformSphericalDistribution(number_of_particles, **keyword_arguments).result
    particles.x = size * x
    particles.y = size * y
    particles.z = size * z
    return particles

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
:argument keyword_arguments:    Optional arguments to UniformSphericalDistribution
"""
def new_spherical_particle_distribution(number_of_particles, 
        radial_density_func = None,     # not yet supported, specify radii and densities tables:
        radii = None, densities = None, 
        total_mass = None, size = None, # if total_mass is not given, it will be deduced from size or max(radii)
        **keyword_arguments):           # optional arguments for UniformSphericalDistribution
    if (radii is None) or (densities is None):
        raise AmuseException("Using an arbitrary radial density function is not yet "
            "supported. Radius and density tables must be passed instead.")
    if total_mass is None:
        total_mass = get_enclosed_mass_from_tabulated((size or max(radii)), radii = radii, densities = densities)
    particles = Particles(number_of_particles)
    particle_mass = total_mass * 1.0 / number_of_particles
    particles.mass = particle_mass
    get_radius_for_enclosed_mass(0 | units.kg, radii = radii, densities = densities)
    x, y, z = UniformSphericalDistribution(number_of_particles, **keyword_arguments).result
    r_old = numpy.sqrt(x*x + y*y + z*z)
    dtype = [('r_old', 'float64'), ('x', 'float64'), ('y', 'float64'), ('z', 'float64')]
    sorted = numpy.sort(numpy.array(zip(r_old, x, y, z), dtype=dtype), order='r_old')
    f_scale = radii.unit.new_quantity(numpy.empty_like(x))
    for i, r_old_i in enumerate(sorted['r_old']):
        f_scale[i] = get_radius_for_enclosed_mass((i+0.5)*particle_mass) / r_old_i
    particles.x = f_scale * sorted['x']
    particles.y = f_scale * sorted['y']
    particles.z = f_scale * sorted['z']
    return particles
    
