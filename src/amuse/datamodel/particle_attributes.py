import numpy

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

from amuse.datamodel import base
from amuse.datamodel import rotation
from amuse.datamodel import ParticlesWithUnitsConverted, AbstractParticleSet, Particle

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
        particles = ParticlesWithUnitsConverted(particles, convert_nbody.as_converter_from_generic_to_si())
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
    scale_factor = numpy.sqrt(target_energy / kinetic_energy)
    particles.velocity *= scale_factor


def center_of_mass(particles):
    """
    Returns the center of mass of the particles set.
    The center of mass is defined as the average
    of the positions of the particles, weighted by their masses.

    >>> from amuse.datamodel import Particles
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

    >>> from amuse.datamodel import Particles
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

def total_momentum(particles):
    """
    Returns the total momentum of the particles set.

    >>> from amuse.datamodel import Particles
    >>> particles = Particles(2)
    >>> particles.vx = [-1.0, 1.0] | units.ms
    >>> particles.vy = [0.0, 0.0] | units.ms
    >>> particles.vz = [0.0, 0.0] | units.ms
    >>> particles.mass = [1.0, 1.0] | units.kg
    >>> particles.total_momentum()
    quantity<[0.0, 0.0, 0.0] m * kg * s**-1>
    """
    masses = particles.mass
    x_values = particles.vx
    y_values = particles.vy
    z_values = particles.vz

    momx = (masses * x_values).sum()
    momy = (masses * y_values).sum()
    momz = (masses * z_values).sum()

    return quantities.VectorQuantity.new_from_scalar_quantities(momx,
        momy, momz)

def total_angular_momentum(particles):
    """
    Returns the total angular momentum of the particles set.

    >>> from amuse.datamodel import Particles
    >>> particles = Particles(2)
    >>> particles.x = [-1.0, 1.0] | units.m
    >>> particles.y = [0.0, 0.0] | units.m
    >>> particles.z = [0.0, 0.0] | units.m
    >>> particles.vx = [0.0, 0.0] | units.ms
    >>> particles.vy = [-1.0, 1.0] | units.ms
    >>> particles.vz = [0.0, 0.0] | units.ms
    >>> particles.mass = [1.0, .5] | units.kg
    >>> particles.total_angular_momentum()
    quantity<[0.0, 0.0, 1.5] m**2 * kg * s**-1>
    """
    m = particles.mass
    x = particles.x
    y = particles.y
    z = particles.z
    vx = particles.vx
    vy = particles.vy
    vz = particles.vz

    lx=(m*(y*vz-z*vy)).sum()
    ly=(m*(z*vx-x*vz)).sum()
    lz=(m*(x*vy-y*vx)).sum()

    return quantities.VectorQuantity.new_from_scalar_quantities(lx,
        ly, lz)

def moment_of_inertia(particles):
    """
    Returns the total moment of inertia (about the Z axis) of the particle
    set.
    """
    m = particles.mass
    x = particles.x
    y = particles.y

    return (m * (x**2 + y**2)).sum()


def kinetic_energy(particles):
    """
    Returns the total kinetic energy of the
    particles in the particles set.

    >>> from amuse.datamodel import Particles
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

    >>> from amuse.datamodel import Particles
    >>> particles = Particles(2)
    >>> particles.x = [0.0, 1.0] | units.m
    >>> particles.y = [0.0, 0.0] | units.m
    >>> particles.z = [0.0, 0.0] | units.m
    >>> particles.mass = [1.0, 1.0] | units.kg
    >>> particles.potential_energy()
    quantity<-6.67428e-11 m**2 * kg * s**-2>
    """

    if len(particles) < 2:
        return zero * G

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

def thermal_energy(particles):
    """
    Returns the total internal energy of the (gas)
    particles in the particles set.

    >>> from amuse.datamodel import Particles
    >>> particles = Particles(2)
    >>> particles.u = [0.5, 0.5] | units.ms**2
    >>> particles.mass = [1.0, 1.0] | units.kg
    >>> particles.thermal_energy()
    quantity<1.0 m**2 * kg * s**-2>
    """
    return (particles.mass * particles.u).sum()


def particle_specific_kinetic_energy(set, particle):
    """
    Returns the specific  kinetic energy of a particle.

    >>> from amuse.datamodel import Particles
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

    >>> from amuse.datamodel import Particles
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

    >>> from amuse.datamodel import Particles
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

    >>> from amuse.datamodel import Particles
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

    >>> from amuse.datamodel import Particles
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

    >>> from amuse.datamodel import Particles
    >>> particles = Particles(3)
    >>> particles.mass = [1.0, 2.0, 3.0] | units.kg
    >>> particles.total_mass()
    quantity<6.0 kg>
    """
    return particles.mass.sum()

def total_radius(particles):
    """
    Returns the total radius (maximum distance from center) of the particles set.

    >>> from amuse.datamodel import Particles
    >>> particles = Particles(3)
    >>> particles.mass = [1.0, 2.0, 3.0] | units.kg
    >>> particles.position = [0.0, 0.0, 0.0] | units.m
    >>> particles.x = [0.0, 3.0, 6.0] | units.m
    >>> particles.total_radius()
    quantity<4.0 m>
    """
    return (particles.position - particles.center_of_mass()).lengths_squared().amax().sqrt()

# move_to_center??
def get_binaries(particles,hardness=10,G = constants.G):
    """
    returns the binaries in a particleset. binaries are selected according to a hardness criterion [hardness=10]
    This function returns the binaries as a list of i,j particles. Triple detection is not done.
    
    >>> from amuse import datamodel
    >>> m = [1,1,1] | units.MSun
    >>> x = [-1,1,0] | units.AU
    >>> y = [0,0,1000] | units.AU
    >>> z = [0,0,0] | units.AU
    >>> vx = [0,0,0] | units.kms
    >>> vy = [1.,-1.,0] | units.kms
    >>> vz = [0,0,0] | units.kms
    >>> particles = datamodel.create_particle_set( mass=m,x=x,y=y,z=z,vx=vx,vy=vy,vz=vz )
    >>> binaries = particles.get_binaries()
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

def densitycentre_coreradius_coredens(particles,unit_converter=None,number_of_neighbours=7):
    """
    calculate position of the density centre, coreradius and coredensity

    >>> import numpy
    >>> from amuse.ic.plummer import new_plummer_sphere
    >>> numpy.random.seed(1234)
    >>> particles=new_plummer_sphere(100)
    >>> pos,coreradius,coredens=particles.densitycentre_coreradius_coredens()
    >>> print coreradius
    0.404120092331 length
    """

    from amuse.community.hop.interface import Hop
    hop=Hop(unit_converter=unit_converter)
    hop.particles.add_particles(particles)
    hop.parameters.density_method=2
    hop.parameters.number_of_neighbors_for_local_density=number_of_neighbours
    hop.calculate_densities()

    density=hop.particles.density
    x=hop.particles.x
    y=hop.particles.y
    z=hop.particles.z
    rho=density.amax()

    total_density=numpy.sum(density)
    x_core=numpy.sum(density*x)/total_density
    y_core=numpy.sum(density*y)/total_density
    z_core=numpy.sum(density*z)/total_density

    rc = (density * ((x-x_core)**2+(y-y_core)**2+(z-z_core)**2).sqrt()).sum() / total_density
    hop.stop()
    
    return [x_core,y_core,z_core],rc,rho

def new_particle_from_cluster_core(particles, unit_converter=None, density_weighting_power=2,cm=None):
    """
    Uses Hop to find the density centre (core) of a particle distribution
    and stores the properties of this core on a particle:
    position, velocity, (core) radius and (core) density.
    
    Particles are assigned weights that depend on the density (as determined by 
    Hop) to a certain power.
    The default weighting power is 2, which is most commonly used. Set 
    density_weighting_power to 1 in order to get the original weighting of 
    Casertano & Hut (1985, ApJ, 298, 80).
    
    :argument unit_converter: Required if the particles are in SI units
    :argument density_weighting_power: Particle properties are weighted by density to this power
    """
    from amuse.community.hop.interface import Hop
    hop = Hop(unit_converter=unit_converter)
    in_hop = hop.particles.add_particles(particles)
    hop.parameters.density_method = 2
    hop.parameters.number_of_neighbors_for_local_density = 7
    hop.calculate_densities()
    density = in_hop.density.copy()
    hop.stop()
    
    weights = (density**density_weighting_power).reshape((-1,1))
    # Reshape makes sure that density can be multiplied with vectors, e.g. position
    
    result = Particle()
    result.density = density.amax()
    total_weight = weights.sum()
    if cm is None:
      result.position = (weights * particles.position).sum(axis=0) / total_weight
    else:
      result.position = cm
    result.velocity = (weights * particles.velocity).sum(axis=0) / total_weight
    result.radius = (weights.flatten() * (particles.position - result.position).lengths()).sum() / total_weight
    return result

def bound_subset(particles,tidal_radius=None,unit_converter=None, density_weighting_power=2,
                     smoothing_length_squared = zero, G = constants.G,core=None  ):
    """
    find the particles bound to the cluster. Returns a subset of bound particles.

    :argument tidal_radius: particles beyond this are considered not bound
    :argument unit_converter: Required if the particles are in SI units
    :argument density_weighting_power: Particle properties are weighted by density to this power
    :argument smooting_length_squared: the smoothing length for gravity.
    :argument G: gravitational constant, need to be changed for particles in different units systems
    :argument core: (optional) core of the cluster

    >>> from amuse.ic.plummer import new_plummer_model
    >>> from amuse.units import nbody_system
    >>> plum=new_plummer_model(100)
    >>> print len(plum.bound_subset(G=nbody_system.G))
    100
    >>> plum[0].velocity*=100  
    >>> plum[0].position*=100  
    >>> print len(plum.bound_subset(G=nbody_system.G))
    99
    """
    if core is None:
      core=particles.cluster_core( unit_converter, density_weighting_power)
    position=particles.position-core.position
    velocity=particles.velocity-core.velocity
    
    v2=velocity.lengths_squared()
    r2=position.lengths_squared()
    pot=particles.potential(smoothing_length_squared, G)
    
    if tidal_radius is None:
      boundary_radius2=r2.max()
    else:
      boundary_radius2=tidal_radius**2
    
    bs=numpy.where( (r2 <= boundary_radius2) & (pot+0.5*v2 < zero) )[0]
    return particles[bs]

def mass_segregation_Gini_coefficient(particles,unit_converter=None, density_weighting_power=2,
                                        core=None):
    """
    Converse & Stahler 2008 Gini coefficient for cluster.

    :argument unit_converter: Required if the particles are in SI units
    :argument density_weighting_power: Particle properties are weighted by density to this power
    :argument core: (optional) core of the cluster
    
    >>> import numpy
    >>> from amuse.ic.plummer import new_plummer_model
    >>> from amuse.units import nbody_system
    >>> plum=new_plummer_model(100)
    >>> a=numpy.argsort(plum.position.lengths_squared().number)
    >>> plum.mass=0|nbody_system.mass
    >>> plum[a[0]].mass=1|nbody_system.mass
    >>> print plum.mass_segregation_Gini_coefficient()
    1.0
    """                   
    if core is None:
      core=particles.cluster_core( unit_converter, density_weighting_power)

    position=particles.position-core.position

    r2=position.lengths_squared().number
    a=numpy.argsort(r2)
    m=particles.mass.number[a]
    
    nf=1.*numpy.array(range(len(m)))/(len(m)-1.)
    mf=m.cumsum()
    mf=mf/mf[-1]
    
    mfmnf=2*(mf-nf)
    
    return (mfmnf[1:]+mfmnf[:-1]).sum()/2/(len(mf)-1.)

def LagrangianRadii(stars,
                       cm=None,
                       mf=[0.01,0.02,0.05,0.1,0.2,0.5,0.75,0.9,1],
                       unit_converter=None,number_of_neighbours=7):
    """
    Calculate lagrangian radii. Output is radii, mass fraction 

    >>> import numpy
    >>> from amuse.ic.plummer import new_plummer_sphere
    >>> numpy.random.seed(1234)
    >>> parts=new_plummer_sphere(100)
    >>> lr,mf=parts.LagrangianRadii()
    >>> print lr[5]
    0.856966667972 length
    """
    import bisect
    if cm is None:
        cm,rcore,rhocore=stars.densitycentre_coreradius_coredens(
                           unit_converter=unit_converter,number_of_neighbours=number_of_neighbours)
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
    
def find_closest_particle_to(particles,x,y,z):
    """
    return closest particle to x,y,z position

    >>> from amuse.datamodel import Particles
    >>> particles = Particles(2)
    >>> particles.x = [0.0, 2.0] | units.m
    >>> particles.y = [0.0, 0.0] | units.m
    >>> particles.z = [0.0, 0.0] | units.m
    >>> print particles.find_closest_particle_to( -1 | units.m,0.| units.m,0.| units.m).x
    0.0 m
    """
    d2=(particles.x-x)**2+(particles.y-y)**2+(particles.z-z)**2
    return particles[d2.number.argmin()]

def potential_energy_in_field(particles, field_particles, smoothing_length_squared = zero, G = constants.G):
    """
    Returns the total potential energy of the particles in the particles set.

    :argument field_particles: the external field consists of these (i.e. potential energy is calculated relative to the field particles) 
    :argument smooting_length_squared: the smoothing length is added to every distance.
    :argument G: gravitational constant, need to be changed for particles in different units systems

    >>> from amuse.datamodel import Particles
    >>> field_particles = Particles(2)
    >>> field_particles.x = [0.0, 2.0] | units.m
    >>> field_particles.y = [0.0, 0.0] | units.m
    >>> field_particles.z = [0.0, 0.0] | units.m
    >>> field_particles.mass = [1.0, 1.0] | units.kg
    >>> particles = Particles(2)
    >>> particles.x = [1.0, 3.0] | units.m
    >>> particles.y = [0.0, 0.0] | units.m
    >>> particles.z = [0.0, 0.0] | units.m
    >>> particles.mass = [1.0, 1.0] | units.kg
    >>> particles.potential_energy_in_field(field_particles)
    quantity<-2.22476e-10 m**2 * kg * s**-2>
    """
    if len(field_particles) == 0:
        return zero * G
        
    n = len(particles)
    dimensions = particles.position.shape[-1]
    transposed_positions = particles.position.reshape([n,1,dimensions]) 
    dxdydz = transposed_positions - field_particles.position
    dr_squared = (dxdydz**2).sum(-1)
    dr = (dr_squared+smoothing_length_squared).sqrt()
    m_m = particles.mass.reshape([n,1]) * field_particles.mass
    return -G * (m_m / dr).sum()
    
def distances_squared(particles,field_particles):
    """
    Returns the total potential energy of the particles in the particles set.

    :argument field_particles: the external field consists of these (i.e. potential energy is calculated relative to the field particles) 
    
    >>> from amuse.datamodel import Particles
    >>> field_particles = Particles(2)
    >>> field_particles.x = [0.0, 2.0] | units.m
    >>> field_particles.y = [0.0, 0.0] | units.m
    >>> field_particles.z = [0.0, 0.0] | units.m
    >>> particles = Particles(3)
    >>> particles.x = [1.0, 3.0, 4] | units.m
    >>> particles.y = [0.0, 0.0] | units.m
    >>> particles.z = [0.0, 0.0] | units.m
    >>> distances_squared(particles, field_particles)
    quantity<[[ 1.  1.], [ 9.  1.], [ 16.   4.]] m**2>
    """

    n = len(particles)
    dimensions = particles.position.shape[-1]
    transposed_positions = particles.position.reshape([n,1,dimensions]) 
    dxdydz = transposed_positions - field_particles.position
    return (dxdydz**2).sum(-1)
    

def velocity_diff_squared(particles,field_particles):
    """
    Returns the total potential energy of the particles in the particles set.

    :argument field_particles: the external field consists of these (i.e. potential energy is calculated relative to the field particles) 
    
    >>> from amuse.datamodel import Particles
    >>> field_particles = Particles(2)
    >>> field_particles.vx = [0.0, 2.0] | units.m
    >>> field_particles.vy = [0.0, 0.0] | units.m
    >>> field_particles.vz = [0.0, 0.0] | units.m
    >>> particles = Particles(3)
    >>> particles.vx = [1.0, 3.0, 4] | units.m
    >>> particles.vy = [0.0, 0.0] | units.m
    >>> particles.vz = [0.0, 0.0] | units.m
    >>> velocity_diff_squared(particles, field_particles)
    quantity<[[ 1.  1.], [ 9.  1.], [ 16.   4.]] m**2>
    """

    n = len(particles)
    dimensions = particles.velocity.shape[-1]
    transposed_positions = particles.velocity.reshape([n,1,dimensions]) 
    dxdydz = transposed_positions - field_particles.velocity
    return (dxdydz**2).sum(-1)

def Qparameter(parts, distfunc=None):
    """
    Calculates the minimum spanning tree Q parameter (Cartwright & Whitworth 2004)
    for a projection of the particle set.
    
    :argument distfunc:  distfunc is the distance function which can be used to select
    the projection plane.

    """
    from amuse.ext.basicgraph import Graph, MinimumSpanningTreeFromEdges
    if distfunc is None:
      def distfunc(p,q):
        return (((p.x-q.x)**2+(p.y-q.y)**2)**0.5).value_in(p.x.unit)
    N=len(parts)
  
    graph=Graph()
  
    print "making graph"
    for p in parts:
      d=distfunc(p,parts)
      for i,q in enumerate(parts):
        if p!=q:
          graph.add_edge(p,q, d[i] ) 

    all_edges=graph.all_edges()
  
    ml=reduce(lambda x,y: x+y[0],all_edges,zero )/len(all_edges)
  
    print "constructing MST"
    mst=MinimumSpanningTreeFromEdges(all_edges)
  
    mlmst=reduce(lambda x,y: x+y[0],mst, zero )/len(mst)
# normalize  
    mlmst=mlmst/(N*numpy.pi)**0.5*(N-1)
  
    print "Q:",mlmst/ml
    return mlmst/ml

def connected_components(parts, threshold=None, distfunc=None, verbose=False):
    """
    return a list of connected component subsets of particles, connected if the distfunc
    is smaller than the threshold.
    
    :argument threshold: value of the threshold. Must have consistent units with distfunc
    :argument distfunc: distance or weight function. Must have consistent units with threshold
    """
    if threshold is None:
      threshold=1. | parts.x.unit
    
    if distfunc is None:
      def distfunc(p,q):
        return (((p.x-q.x)**2+(p.y-q.y)**2+(p.z-q.z)**2)**0.5)
  
    if verbose: print "making CC"
    tocheck=range(len(parts))
    cc=[]
    while len(tocheck)>0:
       p=tocheck.pop()
       stack=[p]
       currentcc=[p]
       
       while len(stack)>0 and len(tocheck)>0:
         p=stack.pop()
         
         d=distfunc(parts[p],parts[tocheck]).value_in(threshold.unit)
         toadd=[ tocheck.pop(i) for i in reversed(range(len(tocheck))) if d[i] < threshold.number ]
         stack.extend(toadd)
         currentcc.extend(toadd)
       cc.append(parts[currentcc])  
         
    if verbose: print "done"
    if verbose: print "number of CC:",len(cc)
    return cc


AbstractParticleSet.add_global_function_attribute("center_of_mass", center_of_mass)
AbstractParticleSet.add_global_function_attribute("center_of_mass_velocity", center_of_mass_velocity)
AbstractParticleSet.add_global_function_attribute("kinetic_energy", kinetic_energy)
AbstractParticleSet.add_global_function_attribute("potential_energy", potential_energy)
AbstractParticleSet.add_global_function_attribute("thermal_energy", thermal_energy)
AbstractParticleSet.add_global_function_attribute("virial_radius", virial_radius)
AbstractParticleSet.add_global_function_attribute("total_mass", total_mass)
AbstractParticleSet.add_global_function_attribute("total_radius", total_radius)
AbstractParticleSet.add_global_function_attribute("total_momentum", total_momentum)
AbstractParticleSet.add_global_function_attribute("total_angular_momentum", total_angular_momentum)
AbstractParticleSet.add_global_function_attribute("moment_of_inertia", moment_of_inertia)

AbstractParticleSet.add_global_function_attribute("potential_energy_in_field", potential_energy_in_field)

AbstractParticleSet.add_global_vector_attribute("position", ["x","y","z"])
AbstractParticleSet.add_global_vector_attribute("velocity", ["vx","vy","vz"])
AbstractParticleSet.add_global_vector_attribute("acceleration", ["ax","ay","az"])
AbstractParticleSet.add_global_vector_attribute("angularmomentum", ["Lx","Ly","Lz"])
AbstractParticleSet.add_global_vector_attribute("oblateness", ["j2","j4","j6"])

AbstractParticleSet.add_global_function_attribute("specific_kinetic_energy", specific_kinetic_energy, particle_specific_kinetic_energy)
AbstractParticleSet.add_global_function_attribute("potential", particleset_potential, particle_potential)

AbstractParticleSet.add_global_function_attribute("move_to_center", move_to_center)
AbstractParticleSet.add_global_function_attribute("scale_to_standard", scale_to_standard)
AbstractParticleSet.add_global_function_attribute("rotate", rotation.rotate)

AbstractParticleSet.add_global_function_attribute("binaries", get_binaries)
AbstractParticleSet.add_global_function_attribute("get_binaries", get_binaries)

AbstractParticleSet.add_global_function_attribute("densitycentre_coreradius_coredens", densitycentre_coreradius_coredens)
AbstractParticleSet.add_global_function_attribute("new_particle_from_cluster_core", new_particle_from_cluster_core)

AbstractParticleSet.add_global_function_attribute("cluster_core", new_particle_from_cluster_core)
AbstractParticleSet.add_global_function_attribute("bound_subset", bound_subset)
AbstractParticleSet.add_global_function_attribute("mass_segregation_Gini_coefficient", mass_segregation_Gini_coefficient)

AbstractParticleSet.add_global_function_attribute("LagrangianRadii", LagrangianRadii)
AbstractParticleSet.add_global_function_attribute("find_closest_particle_to", find_closest_particle_to)

AbstractParticleSet.add_global_function_attribute("Qparameter", Qparameter)
AbstractParticleSet.add_global_function_attribute("connected_components", connected_components)

