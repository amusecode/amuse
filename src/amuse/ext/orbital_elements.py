import numpy

from amuse.units import units,nbody_system,constants
from amuse.datamodel import Particles,rotation

# wrong??
def new_binary_from_orbital_elements2(
        mass1,
        mass2,
        semi_major_axis, 
        eccentricity = 0,
        true_anomaly = 0, 
        inclination = 0,
        longitude_of_the_ascending_node = 0,
        argument_of_periapsis = 0,
        G = nbody_system.G
    ):
    inclination = numpy.radians(inclination)
    argument_of_periapsis = numpy.radians(argument_of_periapsis)
    longitude_of_the_ascending_node = numpy.radians(longitude_of_the_ascending_node)
    true_anomaly = numpy.radians(true_anomaly)

    q   = mass2 / mass1
    mu  = G * (mass2 + mass1)
    
    semilatus_rectum  = semi_major_axis*(1-eccentricity**2)
    radius            = semilatus_rectum/(1+eccentricity*numpy.cos(true_anomaly))
    velocity          = numpy.sqrt(mu/semilatus_rectum)*eccentricity*numpy.sin(true_anomaly)
    
    specific_relative_angular_momentum = numpy.sqrt(mu*semilatus_rectum)
    
    r = radius.new_zeros_array(3)
    v = velocity.new_zeros_array(3)
    
    r[0] = radius * numpy.cos(true_anomaly)
    r[1] = radius * numpy.sin(true_anomaly) 
    
    v[0] = velocity * numpy.cos(true_anomaly) - specific_relative_angular_momentum / radius * numpy.sin(true_anomaly)
    v[1] = velocity * numpy.sin(true_anomaly) + specific_relative_angular_momentum / radius * numpy.cos(true_anomaly)

    result = Particles(2)
    result.position = radius.new_zeros_array(3)
    result.velocity = velocity.new_zeros_array(3)
    result[0].mass = mass1
    result[1].mass = mass2

    
    result[1].position = r
    result[1].velocity = v
    rotation.rotate(result, argument_of_periapsis, inclination, longitude_of_the_ascending_node)
    return result

def new_binary_from_orbital_elements(
        mass1,
        mass2,
        semimajor_axis, 
        eccentricity = 0,
        true_anomaly = 0, 
        inclination = 0,
        longitude_of_the_ascending_node = 0,
        argument_of_periapsis = 0,
        G=nbody_system.G
    ):
    """ 

    Function that returns two-particle Particle set, with the second 
    particle position and velocities computed from the input orbital 
    elements. angles in degrees, inclination between 0 and 180

    """
    
    inclination = numpy.radians(inclination)
    argument_of_periapsis = numpy.radians(argument_of_periapsis)
    longitude_of_the_ascending_node = numpy.radians(longitude_of_the_ascending_node)
    true_anomaly = numpy.radians(true_anomaly)

    cos_true_anomaly = numpy.cos(true_anomaly)
    sin_true_anomaly = numpy.sin(true_anomaly)    

    cos_inclination = numpy.cos(inclination)
    sin_inclination = numpy.sin(inclination)    

    cos_arg_per = numpy.cos(argument_of_periapsis)
    sin_arg_per = numpy.sin(argument_of_periapsis)

    cos_long_asc_nodes = numpy.cos(longitude_of_the_ascending_node)
    sin_long_asc_nodes = numpy.sin(longitude_of_the_ascending_node)

    ### alpha is a unit vector directed along the line of node ###
    alphax = cos_long_asc_nodes*cos_arg_per - sin_long_asc_nodes*sin_arg_per*cos_inclination
    alphay = sin_long_asc_nodes*cos_arg_per + cos_long_asc_nodes*sin_arg_per*cos_inclination
    alphaz = sin_arg_per*sin_inclination
    alpha = [alphax,alphay,alphaz]

    ### beta is a unit vector perpendicular to alpha and the orbital angular momentum vector ###
    betax = -cos_long_asc_nodes*sin_arg_per - sin_long_asc_nodes*cos_arg_per*cos_inclination
    betay = -sin_long_asc_nodes*sin_arg_per + cos_long_asc_nodes*cos_arg_per*cos_inclination
    betaz = cos_arg_per*sin_inclination
    beta = [betax,betay,betaz]

#    print 'alpha',alphax**2+alphay**2+alphaz**2 # For debugging; should be 1
#    print 'beta',betax**2+betay**2+betaz**2 # For debugging; should be 1

    ### Relative position and velocity ###
    separation = semimajor_axis*(1.0 - eccentricity**2)/(1.0 + eccentricity*cos_true_anomaly) # Compute the relative separation
    position_vector = separation*cos_true_anomaly*alpha + separation*sin_true_anomaly*beta
    velocity_tilde = (G*(mass1 + mass2)/(semimajor_axis*(1.0 - eccentricity**2))).sqrt() # Common factor
    velocity_vector = -1.0*velocity_tilde*sin_true_anomaly*alpha + velocity_tilde*(eccentricity + cos_true_anomaly)*beta

    result = Particles(2)
    result[0].mass = mass1
    result[1].mass = mass2
    
    result[1].position = position_vector
    result[1].velocity = velocity_vector
    return result
    
def orbital_elements_from_binary( binary, G=nbody_system.G):
    """ 

    Function that computes orbital elements from given two-particle set. 
    Elements are computed for the second particle in this set and the 
    return values are: mass1, mass2, semimajor axis, eccentricity, 
    cosine of true anomaly, cosine of inclination, cosine of the 
    longitude of the ascending node and the cosine of the argument of 
    pericenter. In case of a perfectly circular orbit the true anomaly 
    and argument of pericenter cannot be determined; in this case the 
    return values are 1.0 for both cosines. 

    """
    if len(binary)>2:
      raise Exception("expects binary or single part")

    if len(binary)==2:
      mass1=binary[0].mass
      mass2=binary[1].mass
      position = binary[1].position-binary[0].position
      velocity = binary[1].velocity-binary[0].velocity
      total_mass = mass1 + mass2
    else:
      mass1=binary.mass
      mass2=0.*mass1
      position = binary.position
      velocity = binary.velocity
      total_mass = mass1
      
    specific_energy = (1.0/2.0)*velocity.lengths_squared() - G*total_mass/position.lengths()
    specific_angular_momentum = position.cross(velocity)
    specific_angular_momentum_norm = specific_angular_momentum.lengths()    
    specific_angular_momentum_unit=specific_angular_momentum/specific_angular_momentum_norm

    semimajor_axis = -G*total_mass/(2.0*specific_energy)

    eccentricity_argument = 2.0*specific_angular_momentum_norm**2*specific_energy/(G**2*total_mass**2)
    if (eccentricity_argument <= -1): eccentricity = 0.0
    else: eccentricity = numpy.sqrt(1.0 + eccentricity_argument)

    ### Orbital inclination ###
    inclination = numpy.degrees(numpy.arccos(specific_angular_momentum.z/specific_angular_momentum_norm))

    ### Longitude of ascending nodes, with reference direction along x-axis ###
    z_vector = [0.,0.,1.] | units.none
    ascending_node_vector = z_vector.cross(specific_angular_momentum)
    if ascending_node_vector.lengths().number==0:
      ascending_node_vector_unit= numpy.array([1.,0.,0.]) 
    else:
      ascending_node_vector_unit = ascending_node_vector/ascending_node_vector.lengths()
    long_asc_node=numpy.degrees(numpy.arctan2(ascending_node_vector_unit[1],ascending_node_vector_unit[0]))

    ### Argument of periapsis and true anomaly, using eccentricity a.k.a. Laplace-Runge-Lenz vector ###
    mu = G*total_mass ### Argument of pericenter ###
    position_unit = position/position.lengths()
    e_vector = ( (1.0/mu)*velocity.cross(specific_angular_momentum) - position_unit ) | units.none
    if (e_vector.lengths() == 0.0): ### Argument of pericenter and true anomaly cannot be determined for e = 0, in this case return 1.0 for the cosines ###
        cos_arg_per = 1.0
        arg_per=0.
        cos_true_anomaly = 1.0
        true_anomaly=0.
    else:
        e_vector_unit = e_vector/e_vector.lengths()
        
        cos_arg_per = e_vector_unit.dot(ascending_node_vector_unit)
        e_cross_an=numpy.cross(e_vector_unit,ascending_node_vector_unit)
        ss=-numpy.sign(specific_angular_momentum_unit.dot(e_cross_an))
        sin_arg_per = ss*(e_cross_an**2).sum()**0.5
        arg_per=numpy.degrees(numpy.arctan2(sin_arg_per,cos_arg_per))


        cos_true_anomaly = e_vector_unit.dot(position_unit)
        e_cross_pos=numpy.cross(e_vector_unit,position_unit)
        ss=numpy.sign(specific_angular_momentum_unit.dot(e_cross_pos))
        sin_true_anomaly = ss*(e_cross_pos**2).sum()**0.5
        true_anomaly=numpy.degrees(numpy.arctan2(sin_true_anomaly,cos_true_anomaly))

        
    return mass1, mass2, semimajor_axis, eccentricity, true_anomaly, inclination, long_asc_node, arg_per
