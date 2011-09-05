import numpy

from amuse.support.data import base
from amuse.support.data.particles import ParticlesWithUnitsConverted, AbstractParticleSet
from amuse.units import nbody_system
from amuse.units import quantities
from amuse.units import constants
from amuse.units import units
from amuse.units.quantities import zero
from amuse.units.quantities import VectorQuantity
from amuse.units.quantities import Quantity
from amuse.units.quantities import new_quantity
from amuse.units.quantities import zero

from amuse.support import exceptions

def move_to_center(particles):
    """
    Move the particle positions to the center of mass and
    move the particle velocity to the center of mass velocity.

    Implemented as::

        particles.position -= particles.center_of_mass()
        particles.velocity -= particles.center_of_mass_velocity()
    """
    particles.position -= particles.center_of_mass()
    particles.velocity -= particles.center_of_mass_velocity()


def scale_to_standard(particles, convert_nbody = None,
                      smoothing_length_squared = zero):
    """
    Scale the particles to a standard NBODY model with
    **total mass=1**, **kinetic energy=0.25** and
    **potential_energy=0.5** (or **viridial_radius=1.0**)

    :argument convert_nbody: the scaling is in nbody units,
        when the particles are in si units a convert_nbody is needed
    :argument smoothing_length_squared: needed for calculating
        the potential energy correctly.
    """
    if not convert_nbody is None:
        particles = ParticlesWithUnitsConverted(particles, convert_nbody.as_converter_from_nbody_to_si())
        if not smoothing_length_squared is zero:
            smoothing_length_squared = convert_nbody.to_nbody(smoothing_length_squared)

    # Proper order is to scale mass, then length, then velocities.
    # Simple length scaling for the potential works only in the
    # unsoftened case.  In general, it may not be possible to force
    # the potential to -0.5, so perhaps best to stop after the simple
    # scaling.  We can always scale the velocities to get the correct
    # virial ratio (and hence virial equilibrium).

    total_mass = particles.mass.sum()
    scale_factor = ((1 | total_mass.unit) / total_mass)
    particles.mass *= scale_factor

    potential_energy \
        = particles.potential_energy(G=nbody_system.G,
                         smoothing_length_squared = smoothing_length_squared)
    target_energy = -0.5 | nbody_system.energy
    scale_factor = (potential_energy / target_energy)	# unsoftened only...
    particles.position *= scale_factor
    if smoothing_length_squared == zero:
        potential_energy = target_energy
    else:
        potential_energy \
            = particles.potential_energy(G=nbody_system.G,
                         smoothing_length_squared = smoothing_length_squared)

    kinetic_energy = particles.kinetic_energy()
    target_energy =  -0.5*potential_energy
    scale_factor = (target_energy / kinetic_energy).sqrt()
    particles.velocity *= scale_factor


def center_of_mass(particles):
    """
    Returns the center of mass of the particles set.
    The center of mass is defined as the average
    of the positions of the particles, weighted by their masses.

    >>> from amuse.support.data.core import Particles
    >>> particles = Particles(2)
    >>> particles.x = [-1.0, 1.0] | units.m
    >>> particles.y = [0.0, 0.0] | units.m
    >>> particles.z = [0.0, 0.0] | units.m
    >>> particles.mass = [1.0, 1.0] | units.kg
    >>> particles.center_of_mass()
    quantity<[0.0, 0.0, 0.0] m>
    """

    masses = particles.mass
    x_values = particles.x
    y_values = particles.y
    z_values = particles.z

    total_mass = masses.sum()
    massx = (masses * x_values).sum()
    massy = (masses * y_values).sum()
    massz = (masses * z_values).sum()

    return quantities.VectorQuantity.new_from_scalar_quantities(massx/total_mass,
        massy/total_mass,
        massz/total_mass)

def center_of_mass_velocity(particles):
    """
    Returns the center of mass velocity of the particles set.
    The center of mass velocity is defined as the average
    of the velocities of the particles, weighted by their masses.

    >>> from amuse.support.data.core import Particles
    >>> particles = Particles(2)
    >>> particles.vx = [-1.0, 1.0] | units.ms
    >>> particles.vy = [0.0, 0.0] | units.ms
    >>> particles.vz = [0.0, 0.0] | units.ms
    >>> particles.mass = [1.0, 1.0] | units.kg
    >>> particles.center_of_mass_velocity()
    quantity<[0.0, 0.0, 0.0] m * s**-1>
    """


    masses = particles.mass
    x_values = particles.vx
    y_values = particles.vy
    z_values = particles.vz

    total_mass = masses.sum()
    massx = (masses * x_values).sum()
    massy = (masses * y_values).sum()
    massz = (masses * z_values).sum()

    return quantities.VectorQuantity.new_from_scalar_quantities(massx/total_mass,
        massy/total_mass,
        massz/total_mass)

def kinetic_energy(particles):
    """
    Returns the total kinetic energy of the
    particles in the particles set.

    >>> from amuse.support.data.core import Particles
    >>> particles = Particles(2)
    >>> particles.vx = [-1.0, 1.0] | units.ms
    >>> particles.vy = [0.0, 0.0] | units.ms
    >>> particles.vz = [0.0, 0.0] | units.ms
    >>> particles.mass = [1.0, 1.0] | units.kg
    >>> particles.kinetic_energy()
    quantity<1.0 m**2 * kg * s**-2>
    """

    mass = particles.mass
    vx = particles.vx
    vy = particles.vy
    vz = particles.vz
    v_squared = (vx * vx) + (vy * vy) + (vz * vz)
    m_v_squared = mass * v_squared
    return 0.5 * m_v_squared.sum()


def potential_energy(particles, smoothing_length_squared = zero, G = constants.G):
    """
    Returns the total potential energy of the particles in the particles set.

    :argument smooting_length_squared: the smoothing length is added to every distance.
    :argument G: gravitational constant, need to be changed for particles in different units systems

    >>> from amuse.support.data.core import Particles
    >>> particles = Particles(2)
    >>> particles.x = [0.0, 1.0] | units.m
    >>> particles.y = [0.0, 0.0] | units.m
    >>> particles.z = [0.0, 0.0] | units.m
    >>> particles.mass = [1.0, 1.0] | units.kg
    >>> particles.potential_energy()
    quantity<-6.67428e-11 m**2 * kg * s**-2>
    """

    if len(particles) < 2:
        return G * zero

    mass = particles.mass
    x_vector = particles.x
    y_vector = particles.y
    z_vector = particles.z

    sum_of_energies = zero

    for i in range(len(particles) - 1):
        x = x_vector[i]
        y = y_vector[i]
        z = z_vector[i]
        dx = x - x_vector[i+1:]
        dy = y - y_vector[i+1:]
        dz = z - z_vector[i+1:]
        dr_squared = (dx * dx) + (dy * dy) + (dz * dz)
        dr = (dr_squared+smoothing_length_squared).sqrt()
        m_m = mass[i] * mass[i+1:]

        energy_of_this_particle = (m_m / dr).sum()
        sum_of_energies -= energy_of_this_particle

    return G * sum_of_energies



def particle_specific_kinetic_energy(set, particle):
    """
    Returns the specific  kinetic energy of a particle.

    >>> from amuse.support.data.core import Particles
    >>> particles = Particles(2)
    >>> particles.vx = [0.0, 1.0] | units.ms
    >>> particles.vy = [0.0, 0.0] | units.ms
    >>> particles.vz = [0.0, 0.0] | units.ms
    >>> particles.mass = [1.0, 1.0] | units.kg
    >>> particles[1].specific_kinetic_energy()
    quantity<0.5 m**2 * s**-2>
    """

    return 0.5*(particle.velocity**2).sum()

def specific_kinetic_energy(particles):
    """
    Returns the specific kinetic energy of a particle.

    >>> from amuse.support.data.core import Particles
    >>> particles = Particles(2)
    >>> particles.vx = [1.0, 1.0] | units.ms
    >>> particles.vy = [0.0, 0.0] | units.ms
    >>> particles.vz = [0.0, 0.0] | units.ms
    >>> particles.mass = [1.0, 1.0] | units.kg
    >>> particles.specific_kinetic_energy()
    quantity<[0.5, 0.5] m**2 * s**-2>
    """

    return 0.5*(particles.vx**2+particles.vy**2+particles.vz**2)


def particle_potential(set, particle, smoothing_length_squared = zero, gravitationalConstant = constants.G):
    """
    Returns the potential energy of a particle.

    >>> from amuse.support.data.core import Particles
    >>> particles = Particles(2)
    >>> particles.x = [0.0, 1.0] | units.m
    >>> particles.y = [0.0, 0.0] | units.m
    >>> particles.z = [0.0, 0.0] | units.m
    >>> particles.mass = [1.0, 1.0] | units.kg
    >>> particles[1].potential()
    quantity<-6.67428e-11 m**2 * s**-2>
    """

    particles = set - particle
    dx = particle.x - particles.x
    dy = particle.y - particles.y
    dz = particle.z - particles.z
    dr_squared = (dx * dx) + (dy * dy) + (dz * dz)
    dr = (dr_squared+smoothing_length_squared).sqrt()
    return - gravitationalConstant * (particles.mass / dr).sum()

def particleset_potential(particles, smoothing_length_squared = zero, G = constants.G):
    """
    Returns the potential of the particles in the particles set.

    :argument smooting_length_squared: the smoothing length is added to every distance.
    :argument G: gravitational constant, need to be changed for particles in different units systems

    >>> from amuse.support.data.core import Particles
    >>> particles = Particles(2)
    >>> particles.x = [0.0, 1.0] | units.m
    >>> particles.y = [0.0, 0.0] | units.m
    >>> particles.z = [0.0, 0.0] | units.m
    >>> particles.mass = [1.0, 1.0] | units.kg
    >>> particles.potential()
    quantity<[-6.67428e-11, -6.67428e-11] m**2 * s**-2>
    """

    mass = particles.mass
    x_vector = particles.x
    y_vector = particles.y
    z_vector = particles.z

    potentials = VectorQuantity.zeros(len(mass),mass.unit/x_vector.unit) 

    for i in range(len(particles) - 1):
        x = x_vector[i]
        y = y_vector[i]
        z = z_vector[i]
        dx = x - x_vector[i+1:]
        dy = y - y_vector[i+1:]
        dz = z - z_vector[i+1:]
        dr_squared = (dx * dx) + (dy * dy) + (dz * dz)
        dr = (dr_squared+smoothing_length_squared).sqrt()

        potentials[i]-= (mass[i+1:]/dr).sum()
        potentials[i+1:]-= mass[i]/dr

    return G * potentials


def virial_radius(particles):
    """
    Returns the virial radius of the particles set.
    The virial radius is the inverse of the average inverse
    distance between particles, weighted by their masses.

    >>> from amuse.support.data.core import Particles
    >>> particles = Particles(2)
    >>> particles.x = [-1.0, 1.0] | units.m
    >>> particles.y = [0.0, 0.0] | units.m
    >>> particles.z = [0.0, 0.0] | units.m
    >>> particles.mass = [1.0, 1.0] | units.kg
    >>> particles.virial_radius()
    quantity<4.0 m>
    """
    if len(particles) < 2:
        raise exceptions.AmuseException("Cannot calculate virial radius for a particles set with fewer than 2 particles.")
    partial_sum = zero

    mass = particles.mass
    x_vector = particles.x
    y_vector = particles.y
    z_vector = particles.z

    for i in range(len(particles) - 1):
        x = x_vector[i]
        y = y_vector[i]
        z = z_vector[i]
        dx = x - x_vector[i+1:]
        dy = y - y_vector[i+1:]
        dz = z - z_vector[i+1:]
        dr_squared = (dx * dx) + (dy * dy) + (dz * dz)
        dr = (dr_squared).sqrt()
        m_m = mass[i] * mass[i+1:]
        partial_sum += (m_m / dr).sum()
    return (mass.sum()**2) / (2 * partial_sum)

def total_mass(particles):
    """
    Returns the total mass of the particles set.

    >>> from amuse.support.data.core import Particles
    >>> particles = Particles(3)
    >>> particles.mass = [1.0, 2.0, 3.0] | units.kg
    >>> particles.total_mass()
    quantity<6.0 kg>
    """
    return particles.mass.sum()

# move_to_center??
def get_binaries(particles,hardness=10,G = constants.G):
    """
    returns the binaries in a particleset. binaries are selected according to a hardness criterion [hardness=10]
    This function returns the binaries as a list of i,j particles. Triple detection is not done.
    
    >>> from amuse.support import data
    >>> m = [1,1,1] | units.MSun
    >>> x = [-1,1,0] | units.AU
    >>> y = [0,0,1000] | units.AU
    >>> z = [0,0,0] | units.AU
    >>> vx = [0,0,0] | units.kms
    >>> vy = [1.,-1.,0] | units.kms
    >>> vz = [0,0,0] | units.kms
    >>> particles = data.create_particle_set( mass=m,x=x,y=y,z=z,vx=vx,vy=vy,vz=vz )
    >>> binaries = particles.binaries()
    >>> print len(binaries)
    1
    
    """
    n=len(particles)
    total_Ek=(0.5*particles.mass*(particles.vx**2+particles.vy**2+particles.vz**2)).sum()
    average_Ek=total_Ek/particles.mass.sum()
    max_mass=particles.mass.amax()
    limitE=hardness*average_Ek

    a=numpy.argsort(particles.x.number)

    binaries=[]

    for i in range(n-1):
        j=i+1
        while j<n and (particles.x[a[j]]-particles.x[a[i]])<2*G*max_mass/limitE:
            r2=(particles.x[a[j]]-particles.x[a[i]])**2+ \
               (particles.y[a[j]]-particles.y[a[i]])**2+ \
               (particles.z[a[j]]-particles.z[a[i]])**2 
            v2=(particles.vx[a[j]]-particles.vx[a[i]])**2+ \
               (particles.vy[a[j]]-particles.vy[a[i]])**2+ \
               (particles.vz[a[j]]-particles.vz[a[i]])**2 
            r=r2**0.5
            eb=G*(particles.mass[i]+particles.mass[j])/r-0.5*v2
            if eb > limitE:
                binary=particles[[a[i],a[j]]].copy()
                binary.hardness=eb/average_Ek
                binaries.append(binary)
            j+=1  

    return binaries

def densitycentre_coreradius_coredens(parts):
    """
    calculate position of the density centre, coreradius and coredensity

    >>> import numpy
    >>> from amuse.ext.plummer import MakePlummerModel
    >>> numpy.random.seed(1234)
    >>> parts=MakePlummerModel(100).result
    >>> pos,coreradius,coredens=parts.densitycentre_coreradius_coredens()
    >>> print coreradius
    0.286582946447 length
    """

    from amuse.community.hop.interface import Hop
    hop=Hop()
    hop.particles.add_particles(parts)
    hop.parameters.density_method=2
    hop.parameters.number_of_neighbors_for_local_density=7
    hop.calculate_densities()

    dens=hop.particles.density
    x=hop.particles.x
    y=hop.particles.y
    z=hop.particles.z
    rho=dens.amax()

    tdens=numpy.sum(dens)
    x_core=numpy.sum(dens*x)/tdens
    y_core=numpy.sum(dens*y)/tdens
    z_core=numpy.sum(dens*z)/tdens

    rc=numpy.sqrt(
        numpy.sum(dens**2*((x-x_core)**2+(y-y_core)**2+(z-z_core)**2))/numpy.sum(dens**2))
    hop.stop()
    return [x_core,y_core,z_core],rc,rho

def LagrangianRadii(stars,
                       cm=None,
                       mf=[0.01,0.02,0.05,0.1,0.2,0.5,0.75,0.9,1]
                       ):
    """
    Calculate lagrangian radii. Output is radii, mass fraction 

    >>> import numpy
    >>> from amuse.ext.plummer import MakePlummerModel
    >>> numpy.random.seed(1234)
    >>> parts=MakePlummerModel(100).result
    >>> lr,mf=parts.LagrangianRadii()
    >>> print lr[5]
    0.856966667972 length
    """
    import bisect
    if cm is None:
        cm,rcore,rhocore=stars.densitycentre_coreradius_coredens()
    cmx,cmy,cmz=cm
    r2=(stars.x-cmx)**2+(stars.y-cmy)**2+(stars.z-cmz)**2
    a=numpy.argsort(r2.number)
    rsorted=r2[a]**0.5
    msorted=stars.mass[a].number
    mcum=msorted.cumsum()
    lr=cmx.unit([])
    for f in mf:
        i=bisect.bisect(mcum,mcum[-1]*f)
        if i<=0:
            lr.append(rsorted[0]/2.)
        else:
            lr.append(rsorted[i-1])
    return lr,mf


AbstractParticleSet.add_global_function_attribute("center_of_mass", center_of_mass)
AbstractParticleSet.add_global_function_attribute("center_of_mass_velocity", center_of_mass_velocity)
AbstractParticleSet.add_global_function_attribute("kinetic_energy", kinetic_energy)
AbstractParticleSet.add_global_function_attribute("potential_energy", potential_energy)
AbstractParticleSet.add_global_function_attribute("virial_radius", virial_radius)
AbstractParticleSet.add_global_function_attribute("total_mass", total_mass)

AbstractParticleSet.add_global_vector_attribute("position", ["x","y","z"])
AbstractParticleSet.add_global_vector_attribute("velocity", ["vx","vy","vz"])
AbstractParticleSet.add_global_vector_attribute("acceleration", ["ax","ay","az"])
AbstractParticleSet.add_global_vector_attribute("angularmomentum", ["Lx","Ly","Lz"])
AbstractParticleSet.add_global_vector_attribute("oblateness", ["j2","j4","j6"])

AbstractParticleSet.add_global_function_attribute("specific_kinetic_energy", specific_kinetic_energy, particle_specific_kinetic_energy)
AbstractParticleSet.add_global_function_attribute("potential", particleset_potential, particle_potential)

AbstractParticleSet.add_global_function_attribute("move_to_center", move_to_center)
AbstractParticleSet.add_global_function_attribute("scale_to_standard", scale_to_standard)

AbstractParticleSet.add_global_function_attribute("binaries", get_binaries)

AbstractParticleSet.add_global_function_attribute("densitycentre_coreradius_coredens", densitycentre_coreradius_coredens)

AbstractParticleSet.add_global_function_attribute("LagrangianRadii", LagrangianRadii)
