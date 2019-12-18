"""
orbital element conversion and utility functions

this module provides:

generate_binaries
orbital_elements
get_orbital_elements_from_binary
get_orbital_elements_from_binaries
get_orbital_elements_from_arrays

and the following deprecated functions (assume input
or output angle floats to be degrees):

new_binary_from_orbital_elements
orbital_elements_from_binary
orbital_elements_for_rel_posvel_arrays

"""

import numpy

import warnings

from amuse.units import units, nbody_system
from amuse.units.trigo import cos, sin, arccos, arctan2
from amuse.datamodel import Particles, Particle

from amuse.units.quantities import to_quantity, VectorQuantity


def newton(f, x0, fprime=None, args=(), tol=1.48e-8, maxiter=50):
    if fprime is None:
        warnings.warn("provide fprime")
        return x0
    i = 0
    x = x0
    while (i < maxiter):
        fv = f(x, *args)
        dfv = fprime(x, *args)
        if(dfv == 0):
            return x0, -2
        delta = -fv/dfv
        if(abs(delta) < tol):
            return x+delta, 0
        x = x+delta
        i = i+1
    return x, -1


def true_anomaly_from_eccentric_anomaly(E, e):
    return 2*arctan2((1+e)**0.5*sin(E/2), (1-e)**0.5*cos(E/2))


def equal_length_array_or_scalar(
        array, length=1, mode="continue"
        ):
    """
    Returns 'array' if its length is equal to 'length'. If this is not the
    case, returns an array of length 'length' with values equal to the first
    value of the array (or if 'array' is a scalar, that value. If mode is
    "warn", issues a warning if this happens; if mode is "exception" raises an
    exception in this case.
    """
    try:
        array_length = len(array)
        if array_length == length:
            return array
        else:
            if mode == "warn":
                warnings.warn("Length of array is not equal to %i. Using only\
                        the first value." % length)
                try:
                    unit = array.unit
                    value = array[0].value_in(unit)
                except:
                    unit = units.none
                    value = array[0]
                array = VectorQuantity(
                        array=numpy.ones(length) * value,
                        unit=unit,
                        )
                return array
            elif mode == "exception":
                raise Exception("Length of array is not equal to %i. This is\
                not supported." % length)
    except:
        try:
            unit = array.unit
            value = array.value_in(unit)
        except:
            unit = units.none
            value = array
        array = VectorQuantity(
                array=numpy.ones(length) * value,
                unit=unit,
                )
        if mode == "warn":
            warnings.warn("Using single value for all cases.")
        return array


def center_of_mass_array(
        vectors,
        primary_mass,
        secondary_mass,
        ):
    """
    Returns array of center_of_mass vectors, where primaries are considered to
    be at (0,0,0) and secondaries at 'vectors'.
    """
    total_mass = (primary_mass + secondary_mass).reshape(
            (len(primary_mass), 1)
            )
    center_of_mass_array = (
            (
                vectors
                * secondary_mass.reshape(
                    (len(secondary_mass), 1)
                    )
                )
            / total_mass
            )
    return center_of_mass_array


def orbital_period_to_semimajor_axis(
        T, M1, M2, G=nbody_system.G
        ):
    mu = G * (M1 + M2)
    semi_major_axis = ((T / (2*numpy.pi))**2 * mu)**(1./3.)

    return semi_major_axis


def rel_posvel_arrays_from_orbital_elements(
        primary_mass,
        secondary_mass,
        semi_major_axis,
        eccentricity=0 | units.rad,
        true_anomaly=0 | units.rad,
        inclination=0 | units.rad,
        longitude_of_the_ascending_node=0 | units.rad,
        argument_of_periapsis=0 | units.rad,
        G=nbody_system.G
        ):
    """
    Returns relative positions/velocities for secondaries orbiting primaries.
    If primary_mass is a scalar, assumes the same primary for all secondaries.
    """
    try:
        number_of_secondaries = len(secondary_mass)
    except:
        number_of_secondaries = 1

    # arrays need to be equal to number of secondaries, or have just one value
    primary_mass = equal_length_array_or_scalar(
            primary_mass, length=number_of_secondaries)
    semi_major_axis = equal_length_array_or_scalar(
            semi_major_axis, length=number_of_secondaries)
    eccentricity = equal_length_array_or_scalar(
            eccentricity, length=number_of_secondaries)
    true_anomaly = equal_length_array_or_scalar(
            true_anomaly, length=number_of_secondaries)
    inclination = equal_length_array_or_scalar(
            inclination, length=number_of_secondaries)
    longitude_of_the_ascending_node = equal_length_array_or_scalar(
            longitude_of_the_ascending_node, length=number_of_secondaries)
    argument_of_periapsis = equal_length_array_or_scalar(
            argument_of_periapsis, length=number_of_secondaries)

    cos_true_anomaly = cos(true_anomaly)
    sin_true_anomaly = sin(true_anomaly)

    cos_inclination = cos(inclination)
    sin_inclination = sin(inclination)

    cos_arg_per = cos(argument_of_periapsis)
    sin_arg_per = sin(argument_of_periapsis)

    cos_long_asc_nodes = cos(longitude_of_the_ascending_node)
    sin_long_asc_nodes = sin(longitude_of_the_ascending_node)

    # alpha is a unit vector directed along the line of node
    alphax = (
            cos_long_asc_nodes*cos_arg_per
            - sin_long_asc_nodes*sin_arg_per*cos_inclination
            )
    alphay = (
            sin_long_asc_nodes*cos_arg_per
            + cos_long_asc_nodes*sin_arg_per*cos_inclination
            )
    alphaz = sin_arg_per*sin_inclination
    alpha = numpy.array([alphax, alphay, alphaz])

    # beta is a unit vector perpendicular to alpha and the orbital angular
    # momentum vector
    betax = (
            - cos_long_asc_nodes*sin_arg_per
            - sin_long_asc_nodes*cos_arg_per*cos_inclination
            )
    betay = (
            - sin_long_asc_nodes*sin_arg_per
            + cos_long_asc_nodes*cos_arg_per*cos_inclination
            )
    betaz = cos_arg_per*sin_inclination
    beta = numpy.array([betax, betay, betaz])

    # Relative position and velocity
    separation = (  # Compute the relative separation
            semi_major_axis*(1.0 - eccentricity**2)
            / (1.0 + eccentricity*cos_true_anomaly)
            )
    position_vector = (
            separation*cos_true_anomaly*alpha
            + separation*sin_true_anomaly*beta
            ).T
    velocity_tilde = (
            (
                G*(primary_mass + secondary_mass)
                / (semi_major_axis*(1.0 - eccentricity**2))
                )**0.5
            )  # Common factor
    velocity_vector = (
            -1.0 * velocity_tilde * sin_true_anomaly * alpha
            + velocity_tilde*(eccentricity + cos_true_anomaly)*beta
            ).T

    return position_vector, velocity_vector


def generate_binaries(
        primary_mass,
        secondary_mass,
        semi_major_axis,
        eccentricity=0 | units.rad,
        true_anomaly=0 | units.rad,
        inclination=0 | units.rad,
        longitude_of_the_ascending_node=0 | units.rad,
        argument_of_periapsis=0 | units.rad,
        G=nbody_system.G
        ):
    """
    returns two particlesets, which contain the primaries and the secondaries
    in binary pairs.
    """
    mass_unit = primary_mass.unit
    try:
        number_of_primaries = len(primary_mass)
    except:
        number_of_primaries = 1
        primary_mass = numpy.array(
                [primary_mass.value_in(mass_unit)]
                ) | mass_unit
    try:
        number_of_secondaries = len(secondary_mass)
    except:
        number_of_secondaries = 1
        secondary_mass = numpy.array(
                [secondary_mass.value_in(mass_unit)]
                ) | mass_unit

    # mass arrays need to be the same length
    if number_of_secondaries != number_of_primaries:
        raise Exception("The number of primaries is not the same as the number\
                of secondaries, this is not supported.")

    position_vector, velocity_vector = rel_posvel_arrays_from_orbital_elements(
        primary_mass,
        secondary_mass,
        semi_major_axis,
        eccentricity=eccentricity,
        true_anomaly=true_anomaly,
        inclination=inclination,
        longitude_of_the_ascending_node=longitude_of_the_ascending_node,
        argument_of_periapsis=argument_of_periapsis,
        G=G
        )
    number_of_primaries

    primaries = Particles(number_of_primaries)
    secondaries = Particles(number_of_secondaries)
    primaries.mass = primary_mass
    secondaries.mass = secondary_mass

    centers_of_mass = center_of_mass_array(
            position_vector, primary_mass, secondary_mass)
    centers_of_mass_velocity = center_of_mass_array(
            velocity_vector, primary_mass, secondary_mass)

    primaries.position = - centers_of_mass
    secondaries.position = position_vector - centers_of_mass
    primaries.velocity = - centers_of_mass_velocity
    secondaries.velocity = velocity_vector - centers_of_mass_velocity
    return primaries, secondaries


def new_binary_from_orbital_elements(
        mass1,
        mass2,
        semimajor_axis,
        eccentricity=0,
        true_anomaly=0 | units.deg,
        inclination=0 | units.deg,
        longitude_of_the_ascending_node=0 | units.deg,
        argument_of_periapsis=0 | units.deg,
        G=nbody_system.G
        ):
    """
    returns a two-particle Particle set, with the second particle's position
    and velocities computed from the input orbital elements.
    inclination is given between 0 and 180 deg.
    angles are assumed to be in deg if no unit is given.
    """
    def angle_with_unit(angle, default_unit=units.deg):
        try:
            default_unit = angle.unit
        except:
            angle = angle | default_unit
        return angle

    # If no unit is given for angles, assume they are in degrees
    true_anomaly = angle_with_unit(true_anomaly, default_unit=units.deg)
    inclination = angle_with_unit(inclination, default_unit=units.deg)
    argument_of_periapsis = angle_with_unit(
            argument_of_periapsis,
            default_unit=units.deg
            )
    longitude_of_the_ascending_node = angle_with_unit(
            longitude_of_the_ascending_node,
            default_unit=units.deg
            )
    primary, secondary = generate_binaries(
            mass1, mass2, semimajor_axis,
            eccentricity=eccentricity,
            true_anomaly=true_anomaly,
            inclination=inclination,
            longitude_of_the_ascending_node=longitude_of_the_ascending_node,
            argument_of_periapsis=argument_of_periapsis,
            G=G
            )

    result = Particles()
    result.add_particle(primary[0])
    result.add_particle(secondary[0])
    return result


def get_orbital_elements_from_binary(binary, G=nbody_system.G):
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
    primaries = Particles()
    secondaries = Particles()
    if len(binary) > 2:
        raise Exception("expects binary or single part")
    elif len(binary) == 2:
        primaries.add_particle(binary[0])
        secondaries.add_particle(binary[1])
    else:
        # FIXME: in case of one particle, what do we calculate the orbit of?
        # The method below is what was default before.
        primaries.add_particle(binary[0])
        primaries[0].position *= 0
        primaries[0].velocity *= 0
        secondaries.add_particle(Particle())
        secondaries[0].mass = 0 * primaries[0].mass
        secondaries[0].position = binary.position
        secondaries[0].velocity = binary.velocity

    (
            mass1, mass2, semimajor_axis, eccentricity, true_anomaly,
            inclination, long_asc_node, arg_per
            ) = get_orbital_elements_from_binaries(primaries, secondaries, G=G)
    return (
            mass1[0], mass2[0], semimajor_axis[0], eccentricity[0],
            true_anomaly[0], inclination[0], long_asc_node[0], arg_per[0])


def orbital_elements_from_binary(binary, G=nbody_system.G):
    (
            mass1, mass2, semimajor_axis, eccentricity, true_anomaly,
            inclination, long_asc_node, arg_per
            ) = get_orbital_elements_from_binary(binary, G=G)
    return (
            mass1, mass2, semimajor_axis, eccentricity,
            true_anomaly.value_in(units.deg),
            inclination.value_in(units.deg),
            long_asc_node.value_in(units.deg),
            arg_per.value_in(units.deg))


def get_orbital_elements_from_binaries(
        primaries, secondaries, G=nbody_system.G):
    """
    Function that computes orbital elements from given primaries and
    secondaries.
    Elements are computed for the second particle in this set and the
    return values are: mass1, mass2, semimajor axis, eccentricity,
    cosine of true anomaly, cosine of inclination, cosine of the
    longitude of the ascending node and the cosine of the argument of
    pericenter. In case of a perfectly circular orbit the true anomaly
    and argument of pericenter cannot be determined; in this case the
    return values are 1.0 for both cosines.
    """

    position = secondaries.position - primaries.position
    velocity = secondaries.velocity - primaries.velocity
    mass1 = primaries.mass
    mass2 = secondaries.mass
    total_mass = mass1 + mass2
    semimajor_axis, eccentricity, true_anomaly, inclination, long_asc_node, \
        arg_per = get_orbital_elements_from_arrays(
            position, velocity, total_mass, G=G)

    return (
            mass1, mass2, semimajor_axis, eccentricity, true_anomaly,
            inclination, long_asc_node, arg_per)


def get_orbital_elements_from_arrays(
        rel_position_raw, rel_velocity_raw,
        total_masses, G=nbody_system.G):
    """
    Orbital elements from array of relative positions and velocities vectors,
    based on orbital_elements_from_binary and adapted to work for arrays (each
    line characterises a two body problem).

    For circular orbits (eccentricity=0): returns argument of pericenter = 0.,
        true anomaly = 0.

    For equatorial orbits (inclination=0): longitude of ascending node = 0,
        argument of pericenter = arctan2(e_y,e_x).

    :argument rel_position: array of vectors of relative positions of the
    two-body systems
    :argument rel_velocity: array of vectors of relative velocities of the
    two-body systems
    :argument total_masses: array of total masses for two-body systems
    :argument G: gravitational constant

    :output semimajor_axis: array of semi-major axes
    :output eccentricity: array of eccentricities
    :output period: array of orbital periods
    :output inc: array of inclinations [radians]
    :output long_asc_node: array of longitude of ascending nodes [radians]
    :output arg_per_mat: array of argument of pericenters [radians]
    :output true_anomaly: array of true anomalies [radians]
    """
    if len(numpy.shape(rel_position_raw)) == 1:
        rel_position = numpy.zeros([1, 3]) * rel_position_raw[0]
        rel_position[0, 0] = rel_position_raw[0]
        rel_position[0, 1] = rel_position_raw[1]
        rel_position[0, 2] = rel_position_raw[2]
        rel_velocity = numpy.zeros([1, 3]) * rel_velocity_raw[0]
        rel_velocity[0, 0] = rel_velocity_raw[0]
        rel_velocity[0, 1] = rel_velocity_raw[1]
        rel_velocity[0, 2] = rel_velocity_raw[2]
    else:
        rel_position = rel_position_raw
        rel_velocity = rel_velocity_raw

    separation = (rel_position**2).sum(axis=1)**0.5
    n_vec = len(rel_position)

    speed_squared = (rel_velocity**2).sum(axis=1)

    semimajor_axis = (
            G * total_masses * separation
            / (2. * G * total_masses - separation * speed_squared)
            )

    neg_ecc_arg = (
            (
                to_quantity(rel_position).cross(rel_velocity)**2
                ).sum(axis=-1)
            / (G * total_masses * semimajor_axis)
            )
    filter_ecc0 = (1. <= neg_ecc_arg)
    eccentricity = numpy.zeros(separation.shape)
    eccentricity[~filter_ecc0] = numpy.sqrt(1.0 - neg_ecc_arg[~filter_ecc0])
    eccentricity[filter_ecc0] = 0.

    # angular momentum
    mom = to_quantity(rel_position).cross(rel_velocity)

    # inclination
    inc = arccos(mom[:, 2]/to_quantity(mom).lengths())

    # Longitude of ascending nodes, with reference direction along x-axis
    asc_node_matrix_unit = numpy.zeros(rel_position.shape)
    z_vectors = numpy.zeros([n_vec, 3])
    z_vectors[:, 2] = 1.
    z_vectors = z_vectors | units.none
    ascending_node_vectors = z_vectors.cross(mom)
    filter_non0_incl = (
            to_quantity(ascending_node_vectors).lengths().number > 0.)
    asc_node_matrix_unit[~filter_non0_incl] = numpy.array([1., 0., 0.])
    an_vectors_len = to_quantity(
            ascending_node_vectors[filter_non0_incl]).lengths()
    asc_node_matrix_unit[filter_non0_incl] = normalize_vector(
            ascending_node_vectors[filter_non0_incl],
            an_vectors_len)
    long_asc_node = arctan2(
            asc_node_matrix_unit[:, 1],
            asc_node_matrix_unit[:, 0])

    # Argument of periapsis using eccentricity a.k.a. Laplace-Runge-Lenz vector
    mu = G*total_masses
    pos_unit_vecs = normalize_vector(rel_position, separation)
    mom_len = to_quantity(mom).lengths()
    mom_unit_vecs = normalize_vector(mom, mom_len)
    e_vecs = (
            normalize_vector(
                to_quantity(rel_velocity).cross(mom), mu)
            - pos_unit_vecs
            )

    # Argument of pericenter cannot be determined for e = 0,
    # in this case return 0.0 and 1.0 for the cosines
    e_vecs_norm = (e_vecs**2).sum(axis=1)**0.5
    filter_non0_ecc = (e_vecs_norm > 1.e-15)
    arg_per_mat = VectorQuantity(
            array=numpy.zeros(long_asc_node.shape),
            unit=units.rad)
    cos_arg_per = numpy.zeros(long_asc_node.shape)
    arg_per_mat[~filter_non0_ecc] = 0. | units.rad
    cos_arg_per[~filter_non0_ecc] = 1.

    e_vecs_unit = numpy.zeros(rel_position.shape)
    e_vecs_unit[filter_non0_ecc] = normalize_vector(
            e_vecs[filter_non0_ecc],
            e_vecs_norm[filter_non0_ecc]
            )
    cos_arg_per[filter_non0_ecc] = (
            e_vecs_unit[filter_non0_ecc]
            * asc_node_matrix_unit[filter_non0_ecc]
            ).sum(axis=-1)
    e_cross_an = numpy.zeros(e_vecs_unit.shape)
    e_cross_an[filter_non0_ecc] = numpy.cross(
            e_vecs_unit[filter_non0_ecc],
            asc_node_matrix_unit[filter_non0_ecc]
            )
    e_cross_an_norm = (e_cross_an**2).sum(axis=1)**0.5
    filter_non0_e_cross_an = (e_cross_an_norm != 0.)
    ss = -numpy.sign(
            (
                mom_unit_vecs[filter_non0_e_cross_an]
                * e_cross_an[filter_non0_e_cross_an]
                ).sum(axis=-1)
            )
    # note change in size in sin_arg_per and cos_arg_per; they are not used further        
    sin_arg_per = ss*e_cross_an_norm[filter_non0_e_cross_an]
    cos_arg_per = cos_arg_per[filter_non0_e_cross_an]
    arg_per_mat[filter_non0_e_cross_an] = arctan2(sin_arg_per, cos_arg_per)

    # in case longitude of ascending node is 0, omega=arctan2(e_y,e_x)
    arg_per_mat[~filter_non0_e_cross_an & filter_non0_ecc] = (
            arctan2(
                e_vecs[~filter_non0_e_cross_an & filter_non0_ecc, 1],
                e_vecs[~filter_non0_e_cross_an & filter_non0_ecc, 0]
                )
            )
    filter_negative_zmom = (
            ~filter_non0_e_cross_an
            & filter_non0_ecc
            & (mom[:, 2] < 0.*mom[0, 0])
            )
    arg_per_mat[filter_negative_zmom] = (
            2. * numpy.pi
            - arg_per_mat[filter_negative_zmom]
            )

    # true anomaly
    cos_true_anomaly = (e_vecs_unit*pos_unit_vecs).sum(axis=-1)
    e_cross_pos = numpy.cross(e_vecs_unit, pos_unit_vecs)
    ss2 = numpy.sign((mom_unit_vecs*e_cross_pos).sum(axis=-1))
    sin_true_anomaly = ss2*(e_cross_pos**2).sum(axis=1)**0.5
    true_anomaly = arctan2(sin_true_anomaly, cos_true_anomaly)

    return (
            semimajor_axis, eccentricity, true_anomaly,
            inc, long_asc_node, arg_per_mat
            )


def orbital_elements(*args, **kwargs):
    try:
        if len(args) == 1:
            return get_orbital_elements_from_binary(*args, **kwargs)
        elif len(args) == 2:
            return get_orbital_elements_from_binaries(*args, **kwargs)
        elif len(args) == 3:
            return get_orbital_elements_from_arrays(*args, **kwargs)
        else:
            raise Exception
    except Exception as ex:
        if not ex.args: 
            ex.args=()
        ex.args = ex.args + ("""
  note: orbital elements takes as input either:
  - single two particle set,
  - two sets of primaries and secondaries
  - arrays of rel. position, rel. velocity and masses
  """,)
        raise

def orbital_elements_for_rel_posvel_arrays(
        rel_position_raw, rel_velocity_raw,
        total_masses, G=nbody_system.G):
    (semimajor_axis, eccentricity, true_anomaly, inc, long_asc_node,
        arg_per_mat) = get_orbital_elements_from_arrays(
                rel_position_raw, rel_velocity_raw, total_masses, G)

    true_anomaly = true_anomaly.value_in(units.deg)
    inc = inc.value_in(units.deg)
    long_asc_node = long_asc_node.value_in(units.deg)
    arg_per_mat = arg_per_mat.value_in(units.deg)

    return (
            semimajor_axis, eccentricity, true_anomaly,
            inc, long_asc_node, arg_per_mat
            )


def normalize_vector(vecs, norm, one_dim=False):
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
        vecs_norm[:, 0] = vecs[:, 0]/norm
        vecs_norm[:, 1] = vecs[:, 1]/norm
        vecs_norm[:, 2] = vecs[:, 2]/norm
    return vecs_norm
