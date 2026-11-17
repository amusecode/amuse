"""
Make planets based on the an assumption of oligarchical planetary growth.
(the formalism of Kokubo & Ida (1998))

After:

@article{0004-637X-807-2-157,
 author={Scott Tremaine},
 title={The Statistical Mechanics of Planet Orbits},
 journal={The Astrophysical Journal},
 volume={807},
 number={2},
 pages={157},
 url={http://stacks.iop.org/0004-637X/807/i=2/a=157},
 year={2015},
 }

Refers to (with the actual method as implemented here):


@article{0004-637X-775-1-53,
  author={Brad M. S. Hansen and Norm Murray},
  title={Testing in Situ Assembly with the Kepler Planet Candidate Sample},
  journal={The Astrophysical Journal},
  volume={775},
  number={1},
  pages={53},
  url={http://stacks.iop.org/0004-637X/775/i=1/a=53},
  year={2013}
 }

Actual model:

@article{0004-637X-581-1-666,
  author={Eiichiro Kokubo and Shigeru Ida},
  title={Formation of Protoplanet Systems and Diversity of Planetary Systems},
  journal={The Astrophysical Journal},
  volume={581},
  number={1},
  pages={666},
  url={http://stacks.iop.org/0004-637X/581/i=1/a=666},
  year={2002},
}






EXAMPLE:

from amuse.lab import units, constants
from amuse.ic.make_planets_oligarch import new_system
from amuse.ext.orbital_elements import orbital_elements


def get_combined_Rhill(Mstar, mp1, mp2, a1, a2):
    total_planet_mass = mp1 + mp2
    average_semi_major_axis = (a1 + a2) / 2.0
    return average_semi_major_axis * (total_planet_mass / (3 * Mstar))**(1/3)

Mstar = 1.0 | units.MSun
Rstar = 1.0 | units.RSun
Rmin  = 0.1 | units.AU
Rmax  = 100.0 | units.AU
Mdisk = 0.01 | units.MSun

system = new_system(Mstar, Rstar, Rmin, Rmax, Mdisk, hill_spacing=7)
star = system[0]
planets = system[1:]


for i in range(len(planets)-1):
    a1 = planets[i].semimajor_axis
    a2 = planets[i+1].semimajor_axis
    M1 = planets[i].mass
    M2 = planets[i+1].mass
    
    rhill = get_combined_Rhill(Mstar, M1, M2, a1, a2)
    k = (a2 - a1) / rhill
    print(
        f"mpi = {planets[i].mass.value_in(units.MEarth):.5e} M_Earth, "
        f"mpj = {planets[i+1].mass.value_in(units.MEarth):.5e} M_Earth, "
        f"k = {k:.2f}"
    )
    
    ke = orbital_elements(star+planets[i], G=constants.G)
    assert ke[3] < 1.0, f"Planet {i} has eccentricity {ke[3]:.5f}, which is greater than 1.0"
    assert ke[2] > 0.0 | units.au, f"Planet {i} has semi-major axis {ke[2].value_in(units.AU):.5f} AU, which is not positive"


"""


import numpy
from numpy import random


from amuse.units import quantities
from amuse.units import units
from amuse.units import nbody_system
from amuse.units.core import named_unit
from amuse.units.optparse import OptionParser
from amuse.datamodel import Particle, Particles
from amuse.community.kepler.interface import Kepler
from amuse.io import write_set_to_file


MEarth=named_unit("MEarth", "MEarth", 5.97219e+24 * units.kg)


def get_mass(a, sigma_factor=1.0, mstar=1 | units.MSun):
    s = sigma_factor * (10.0 | (units.g / units.cm**2))
    z = 2 * numpy.pi * a * s * ((a / (5.0 | units.AU) ) ** -(3.0/2.0)) * 4.0 * a
    return (2.0 * z ** 3 / (3.0 * mstar)).sqrt()


def get_orbital_separation(mass, a, mstar=1 | units.MSun, hill_spacing=10.0):
    return hill_spacing * a * (2*mass/(3.0 * mstar)) ** (1.0/3.0)


def get_mass_and_orbital_separation(
    a, 
    sigma_factor=1.0, 
    mstar=1 | units.MSun, 
    hill_spacing=10.0
    ):
    mass = get_mass(a, sigma_factor=sigma_factor, mstar=mstar)
    orbital_separation = get_orbital_separation(
        mass, a, 
        mstar=mstar, 
        hill_spacing=hill_spacing
        )
    return mass, orbital_separation


def new_distribution(x0, x1, sigma_factor=1.0, mstar=1 | units.MSun, hill_spacing=10.0):
    if hill_spacing <= 0:
        raise ValueError("hill_spacing must be finite and positive")

    x = x0
    masses = []
    positions = []
    while x < x1:
        mass, orbital_separation=get_mass_and_orbital_separation(
            x, 
            sigma_factor=sigma_factor, 
            mstar=mstar, 
            hill_spacing=hill_spacing
            )

        masses.append(mass)
        positions.append(x)
        while x + orbital_separation < x1:
            next_mass=get_mass(x + orbital_separation, sigma_factor, mstar)
            mean_a = x + 0.5 * orbital_separation
            step = hill_spacing * mean_a * ((mass + next_mass) / (3.0 * mstar)) ** (1.0/3.0)
            if abs(step - orbital_separation) <= 1.e-10 * step:
                orbital_separation=step
                break

            orbital_separation=step
            
        if x + orbital_separation == x:
            raise ValueError("hill_spacing is too small to advance the orbit")

        x += orbital_separation
    
    return (
        quantities.as_vector_quantity(masses),
        quantities.as_vector_quantity(positions)
    )
    


def is_hit(mass, target_total_mass, accuracy):
    return abs((mass - target_total_mass) / target_total_mass) < accuracy


def new_planet_distribution(
    x0, x1, 
    target_total_mass, 
    accuracy, 
    max_iterations=1000, 
    mstar=1 | units.MSun, 
    hill_spacing=10.0
    ):
    """
    Create a new planet particle set consisting of masses equal to the target_total_mass
    and positions of planets in a disk between x0 and x1.
    
    Args:
        x0:                 minimum radius of the disk
        x1:                 maximum radius of the disk
        target_total_mass:  total mass of the planets in the disk
        accuracy:           accuracy with which the disk mass should be equal to
        max_iterations:     maximum number of iterations to find a solution
        mstar:              mass of central star
        hill_spacing:       spacing in mutual Hill radii

    Returns:
        masses:                vector of planet masses
        positions:             vector of planet positions
    """
    
    min_sigma_factor = 0.1
    max_sigma_factor = 10000.0
    
    current_sigma_factor = min_sigma_factor
    mass, pos = new_distribution(
        x0, x1, current_sigma_factor, 
        mstar=mstar, 
        hill_spacing=hill_spacing
        )
    total_mass = mass.sum()
    iteration = 0
    while not is_hit(total_mass,  target_total_mass, accuracy):
        #print iteration, total_mass / (1 | MEarth), current_sigma_factor
        iteration += 1
        if iteration >= max_iterations:
            break
        
        if total_mass < target_total_mass:
            min_sigma_factor = current_sigma_factor
        elif total_mass > target_total_mass:
            max_sigma_factor = current_sigma_factor
            
        current_sigma_factor = (min_sigma_factor + max_sigma_factor) / 2.0
        
        mass, pos = new_distribution(
            x0, x1, current_sigma_factor, 
            mstar=mstar, 
            hill_spacing=hill_spacing
            )
        total_mass = mass.sum()
    
    return mass, pos


def new_rotation_matrix_from_euler_angles(phi, theta, chi):
    cosp=numpy.cos(phi)
    sinp=numpy.sin(phi)
    cost=numpy.cos(theta)
    sint=numpy.sin(theta)
    cosc=numpy.cos(chi)
    sinc=numpy.sin(chi)
    #see wikipedia: http://en.wikipedia.org/wiki/Rotation_matrix
    return numpy.array(
        [[cost*cosc, -cosp*sinc + sinp*sint*cosc, sinp*sinc + cosp*sint*cosc], 
         [cost*sinc, cosp*cosc + sinp*sint*sinc, -sinp*cosc + cosp*sint*sinc],
         [-sint,  sinp*cost,  cosp*cost]]
        )


def rotate(position, velocity, phi, theta, psi): # theta and phi in radians
    Runit = position.unit
    Vunit = velocity.unit
    matrix = new_rotation_matrix_from_euler_angles(phi, theta, psi)
    return (
        numpy.dot(matrix, position.value_in(Runit)) | Runit,
        numpy.dot(matrix, velocity.value_in(Vunit)) | Vunit
        )

# select Eurler angles randomly. 

def random_Euler_angles():
    phi   = 2 * numpy.pi * random()
    theta = numpy.acos(1 - 2 * random())
    chi   = 2 * numpy.pi * random()
    return phi, theta, chi


def posvel_from_orbital_elements(Mstar, semimajor_axis, eccentricity, kepler, rng=None):
    if rng is None:
        rng=random
        
    mean_anomaly=rng.uniform(0, 2*numpy.pi, 1)
    kepler.initialize_from_elements(
        Mstar, 
        semimajor_axis, 
        eccentricity, 
        mean_anomaly=mean_anomaly
        )
    position = quantities.as_vector_quantity(kepler.get_separation_vector())
    velocity = quantities.as_vector_quantity(kepler.get_velocity_vector())
    return position, velocity
    

def make_planets(
    central_particle, 
    masses, 
    radii, 
    density=3 | units.g/units.cm**3, 
    phi=None, 
    theta=None, 
    eccentricity=0.0, 
    kepler=None, 
    rng=None
    ):
    
    volumes=masses / density
    planet_radii=(3.0 * volumes /  (4.0 * numpy.pi))**(1.0/3.0)
    n=len(masses)
    planet_particles=Particles(n)
    planet_particles.semimajor_axis=radii
    if eccentricity is None:
        eccentricity=numpy.abs(rng.normal(-0.00001,0.00001,n))
    planet_particles.eccentricity=eccentricity
    planet_particles.mass=masses
    planet_particles.radius=planet_radii


    if phi is None:
        phi = numpy.radians(rng.uniform(0.0, 90.0, 1)[0])#rotate under x
    if theta is None:
        theta0 = numpy.radians((rng.normal(-90.0,90.0,1)[0]))#rotate under y
        theta0 = 0
        theta_inclination = numpy.radians(rng.normal(0, 1.0, n )) 
        theta_inclination[0]=0
        theta = theta0 + theta_inclination

    #psi=numpy.radians(rng.uniform(0, 180, 1))[0] #0 # numpy.radians(90) # numpy.radians(rng.uniform(0, 180, 1))[0]
    psi = numpy.radians(rng.uniform(0.0, 180.0, 1))[0] #0 # numpy.radians(90) # numpy.radians(rng.uniform(0, 180, 1))[0]
    com_particle = central_particle.copy()
    for x, t in zip(iter(planet_particles), theta):
        pos,vel = posvel_from_orbital_elements(com_particle.mass + x.mass, x.semimajor_axis, x.eccentricity, kepler, rng)
        pos,vel = rotate(pos, vel, 0, 0, psi) # theta and phi in radians            
        pos,vel = rotate(pos, vel, 0, t, 0) # theta and phi in radians            
        pos,vel = rotate(pos, vel, phi, 0, 0) # theta and phi in radians            
        x.position = pos + com_particle.position
        x.velocity = vel + com_particle.velocity

    return planet_particles

def new_system(
    star_mass=1|units.MSun, 
    star_radius=1|units.RSun, 
    rdisk_min=0.05 | units.AU,
    rdisk_max=10 | units.AU,
    disk_mass=20 | MEarth,
    accuracy=0.0001, 
    planet_density= 3 | units.g/units.cm**3,
    rng=None,
    kepler=None, 
    hill_spacing=10.0
    ):
    """
    Create a particle set consisting of a central star and planets.
    
    Args:
        star_mass:       mass of the central star
        star_radius:     radius of the central star
        rdisk_min:       minimum radius of the disk
        rdisk_max:       maximum radius of the disk
        disk_mass:       total mass of the planets in the disk
        accuracy:        accuracy with which the disk mass should be equal to
        planet_density:  density of the planets (used for planet radii calculation)
        rng:             random number generator
        kepler:          kepler code instance
        hill_spacing:    spacing in mutual Hill radii
    """
            
    central_particle = Particle()
    central_particle.mass = star_mass
    central_particle.position = (0,0,0) | units.AU
    central_particle.velocity = (0,0,0) | units.kms
    central_particle.radius = star_radius
    
    if rng is None:
        rng = numpy.random
        
    converter = nbody_system.nbody_to_si(1|units.MSun, 1 | units.AU)
    
    if kepler is None:
        kepler=Kepler(converter)
        kepler.initialize_code()
        
    m, r = new_planet_distribution(
        rdisk_min, 
        rdisk_max, 
        disk_mass,
        accuracy, 
        mstar=star_mass, 
        hill_spacing=hill_spacing
    )
    
    planets = make_planets(
        central_particle, m, r, 
        density=planet_density, 
        phi=0, 
        theta=None, 
        kepler=kepler, 
        rng=rng
    )
    
    p = Particles()
    p.add_particle(central_particle)
    p.add_particles(planets)
    
    kepler.stop()
    return p
    

def new_option_parser():
    result=OptionParser()
    result.add_option("--hill-spacing", type="float", default=10.0,
                      help="spacing in mutual Hill radii [%default]")
    result.add_option("--seed", 
                      dest="seed", type="int", default=-1,
                      help="random number seed [%default]")
    result.add_option("-o", "--output",
                      dest="output_filename", default="star.h5",
                      help="name of the output filename [%default]")
    result.add_option("-m", "--star-mass", unit=units.MSun, type="float", default=1|units.MSun,
                      dest="star_mass", help='mass of the star [%default]')
    result.add_option("-r", "--star-radius", unit=units.RSun, type="float", default=1|units.RSun,
                      dest="star_radius", help='radius of the star [%default]')
    result.add_option("--disk-r0", unit=units.AU, type="float", default=1|units.AU,
                      dest="rdisk_min", help='minimum radius of the disk [%default]')
    result.add_option("--disk-r1", unit=units.AU, type="float", default=100|units.AU,
                      dest="rdisk_max", help='maximum radius of the disk [%default]')
    result.add_option("--disk-mass", unit=MEarth, type="float", default=100|MEarth,
                      dest="disk_mass", help='mass of the disk [%default]')
    result.add_option("--planet-density", unit=units.g/units.cm**3, type="float", default=3 | units.g/units.cm**3,
                      dest="planet_density", help='density of the planets  (used for planet radii calculation) [%default]')
    result.add_option("--accuracy", type="float", default=0.0001,
                      dest="accuracy", help='how accurate should the disk mass be integrated to[%default]')
    
    return result
    
def main(
    star_mass=1|units.MSun, 
    star_radius=1|units.RSun, 
    rdisk_min=0.05 | units.AU,
    rdisk_max=10 | units.AU,
    disk_mass=200 | MEarth,
    accuracy=0.0001, 
    planet_density= 3 | units.g/units.cm**3,
    output_filename='star.h5',
    seed=-1, 
    hill_spacing=10.0
    ):
            
    if seed < 0:
        rng=random.RandomState()
    else:
        rng=random.RandomState(seed)
    
    
    output=new_system(
        star_mass, 
        star_radius, 
        rdisk_min,
        rdisk_max,
        disk_mass,
        accuracy, 
        planet_density,
        rng=rng, 
        hill_spacing=hill_spacing
        )
    
    star=output[0]
    print("Number of planets generated:", len(star.planets))
    print("Total mass:", star.planets.mass.sum().as_quantity_in(MEarth))
    for i, planet  in enumerate(star.planets):
        print("Planet: {0: 3d} , mass: {1: 8.3f}  MEarth, a: {2: 8.2f} AU".format(i, planet.mass.value_in(MEarth), planet.semimajor_axis.value_in(units.AU)))
    
    write_set_to_file(output, output_filename, 'hdf5', version="2.0",append_to_file=False)
    
    

if __name__ == "__main__":
    main(**new_option_parser().parse_args()[0].__dict__)