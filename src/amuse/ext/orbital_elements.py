import numpy

from amuse.units import units,nbody_system,constants
from amuse.datamodel import Particles,rotation

def newton(f,x0,fprime=None,args=(),tol=1.48e-8,maxiter=50):
    if fprime is None:
        print "provide fprime"
        return x0
    i=0
    x=x0
    while (i<maxiter):
        fv=f(x,*args)
        dfv=fprime(x,*args)
        if(dfv==0):
            return x0,-2
        delta=-fv/dfv
        if(abs(delta)<tol):
            return x+delta,0
        x=x+delta
        i=i+1
    return x,-1    


def true_anomaly_from_eccentric_anomaly(E,e):
  return 2*numpy.arctan2((1+e)**0.5*numpy.sin(E/2),(1-e)**0.5*numpy.cos(E/2))

# E from M,e
# newton solver for M=E-e sin E

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
    
    result.move_to_center()
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
        
        cos_arg_per = numpy.dot(e_vector_unit,ascending_node_vector_unit)
        #cos_arg_per = e_vector_unit.dot(ascending_node_vector_unit)
        e_cross_an=numpy.cross(e_vector_unit,ascending_node_vector_unit)
        ss=-numpy.sign(numpy.dot(specific_angular_momentum_unit,e_cross_an))
        #ss=-numpy.sign(specific_angular_momentum_unit.dot(e_cross_an))
        sin_arg_per = ss*(e_cross_an**2).sum()**0.5
        arg_per=numpy.degrees(numpy.arctan2(sin_arg_per,cos_arg_per))


        cos_true_anomaly = numpy.dot(e_vector_unit,position_unit)
        #cos_true_anomaly = e_vector_unit.dot(position_unit)
        e_cross_pos=numpy.cross(e_vector_unit,position_unit)
        ss=numpy.sign(numpy.dot(specific_angular_momentum_unit,e_cross_pos))
        #ss=numpy.sign(specific_angular_momentum_unit.dot(e_cross_pos))
        sin_true_anomaly = ss*(e_cross_pos**2).sum()**0.5
        true_anomaly=numpy.degrees(numpy.arctan2(sin_true_anomaly,cos_true_anomaly))

        
    return mass1, mass2, semimajor_axis, eccentricity, true_anomaly, inclination, long_asc_node, arg_per

def orbital_elements_for_rel_posvel_arrays(rel_position_raw, rel_velocity_raw, total_masses, G=nbody_system.G):
    """
    Orbital elements from array of relative positions and velocities vectors,
    based on orbital_elements_from_binary and adapted to work for arrays (each line
    characterises a two body problem).
    
    For circular orbits (eccentricity=0): returns argument of pericenter = 0.
    
    For equatorial orbits (inclination=0): longitude of ascending node = 0,
        argument of pericenter = arctan2(e_y,e_x).
    
    :argument rel_position: array of vectors of relative positions of the two-body systems
    :argument rel_velocity: array of vectors of relative velocities of the two-body systems
    :argument total_masses: array of total masses for two-body systems
    :argument G: gravitational constant
    
    :output semimajor_axis: array of semi-major axes
    :output eccentricity: array of eccentricities
    :output period: array of orbital periods
    :output inc: array of inclinations [radians]
    :output long_asc_node: array of longitude of ascending nodes [radians]
    :output arg_per_mat: array of argument of pericenters [radians]
    """
    if len(numpy.shape(rel_position_raw))==1:
        rel_position = numpy.zeros([1,3]) * rel_position_raw[0]
        rel_position[0,0] = rel_position_raw[0]
        rel_position[0,1] = rel_position_raw[1]
        rel_position[0,2] = rel_position_raw[2]
        rel_velocity = numpy.zeros([1,3]) * rel_velocity_raw[0]
        rel_velocity[0,0] = rel_velocity_raw[0]
        rel_velocity[0,1] = rel_velocity_raw[1]
        rel_velocity[0,2] = rel_velocity_raw[2]
    else:
        rel_position = rel_position_raw
        rel_velocity = rel_velocity_raw
    
    separation = rel_position.lengths()    
    n_vec = len(rel_position)
    
    speed_squared = rel_velocity.lengths_squared()
    
    semimajor_axis = (G * total_masses * separation / \
        (2. * G * total_masses - separation * speed_squared))
    
    neg_ecc_arg = (rel_position.cross(rel_velocity)**2).sum(axis=-1)/(G * total_masses * semimajor_axis)   
    filter_ecc0 = (1. <= neg_ecc_arg)
    eccentricity = numpy.zeros(separation.shape)
    eccentricity[~filter_ecc0] = numpy.sqrt( 1.0 - neg_ecc_arg[~filter_ecc0])
    eccentricity[filter_ecc0] = 0.
        
    period = (2 * numpy.pi * (semimajor_axis**1.5) / ((G * total_masses).sqrt()))
    
    # angular momentum
    mom = rel_position.cross(rel_velocity)        
    
    # inclination
    inc = numpy.arccos(mom[:,2]/mom.lengths())
    
    # Longitude of ascending nodes, with reference direction along x-axis
    asc_node_matrix_unit = numpy.zeros(rel_position.shape)
    z_vectors = numpy.zeros([n_vec,3])
    z_vectors[:,2] = 1.
    z_vectors = z_vectors | units.none
    ascending_node_vectors = z_vectors.cross(mom)
    filter_non0_incl = (ascending_node_vectors.lengths().number>0.)
    asc_node_matrix_unit[~filter_non0_incl] = numpy.array([1.,0.,0.])
    an_vectors_len = ascending_node_vectors[filter_non0_incl].lengths()
    asc_node_matrix_unit[filter_non0_incl] = normalize_vector(ascending_node_vectors[filter_non0_incl], an_vectors_len)
    long_asc_node = numpy.arctan2(asc_node_matrix_unit[:,1],asc_node_matrix_unit[:,0])
    
    # Argument of periapsis using eccentricity a.k.a. Laplace-Runge-Lenz vector
    mu = G*total_masses
    pos_unit_vecs = normalize_vector(rel_position, separation)
    mom_len = mom.lengths()
    mom_unit_vecs = normalize_vector(mom, mom_len)
    e_vecs = (normalize_vector(rel_velocity.cross(mom), mu) - pos_unit_vecs)
    
    # Argument of pericenter cannot be determined for e = 0,
    # in this case return 0.0 and 1.0 for the cosines
    e_vecs_norm = (e_vecs**2).sum(axis=1)**0.5
    filter_non0_ecc = (e_vecs_norm > 1.e-15)
    arg_per_mat = numpy.zeros(long_asc_node.shape)
    cos_arg_per = numpy.zeros(long_asc_node.shape)
    arg_per_mat[~filter_non0_ecc] = 0.
    cos_arg_per[~filter_non0_ecc] = 1.
    
    e_vecs_unit = numpy.zeros(rel_position.shape)
    e_vecs_unit[filter_non0_ecc] = normalize_vector(e_vecs[filter_non0_ecc],
                                                    e_vecs_norm[filter_non0_ecc])
    #~ cos_arg_per = numpy.einsum('ij,ji->i', e_vecs_unit[filter_non0_ecc], asc_node_matrix_unit[filter_non0_ecc].T)
    cos_arg_per = (e_vecs_unit[filter_non0_ecc]*asc_node_matrix_unit[filter_non0_ecc]).sum(axis=-1)
    e_cross_an = numpy.zeros(e_vecs_unit.shape)
    e_cross_an[filter_non0_ecc] = numpy.cross(e_vecs_unit[filter_non0_ecc],asc_node_matrix_unit[filter_non0_ecc])
    e_cross_an_norm=(e_cross_an**2).sum(axis=1)**0.5
    filter_non0_e_cross_an = (e_cross_an_norm != 0.)
    #~ ss = -numpy.sign(numpy.einsum('ij,ji->i',mom_unit_vecs[filter_non0_e_cross_an], e_cross_an[filter_non0_e_cross_an].T))
    ss = -numpy.sign((mom_unit_vecs[filter_non0_e_cross_an]*e_cross_an[filter_non0_e_cross_an]).sum(axis=-1))
    sin_arg_per = ss*e_cross_an_norm[filter_non0_e_cross_an]
    arg_per_mat[filter_non0_e_cross_an] = numpy.arctan2(sin_arg_per,cos_arg_per)
    
    # in case longitude of ascenfing node is 0, omega=arctan2(e_y,e_x)
    arg_per_mat[~filter_non0_e_cross_an & filter_non0_ecc] = \
        numpy.arctan2(e_vecs[~filter_non0_e_cross_an & filter_non0_ecc,1], \
                      e_vecs[~filter_non0_e_cross_an & filter_non0_ecc,0])
    filter_negative_zmom = (~filter_non0_e_cross_an & filter_non0_ecc & (mom[:,2]<0.*mom[0,0]))
    arg_per_mat[filter_negative_zmom] = 2.*numpy.pi - arg_per_mat[filter_negative_zmom]
    
    return semimajor_axis, eccentricity, period, inc, long_asc_node, arg_per_mat
    
def normalize_vector(vecs, norm, one_dim = False):
    """
    normalize array of vector quantities
    """
    if one_dim:
        vecs_norm = numpy.zeros(vecs.shape)
        vecs_norm[0] = vecs[0]/norm
        vecs_norm[1] = vecs[1]/norm
        vecs_norm[2] = vecs[2]/norm
    else:
        vecs_norm = numpy.zeros(vecs.shape)
        vecs_norm[:,0] = vecs[:,0]/norm
        vecs_norm[:,1] = vecs[:,1]/norm
        vecs_norm[:,2] = vecs[:,2]/norm
    return vecs_norm
    
    
    

