import numpy
from amuse.lab import *



def get_relative_velocity(total_mass, semimajor_axis, ecc):
    return (constants.G * total_mass * ((1.0 + ecc)/(1.0 - ecc)) / semimajor_axis).sqrt()

def make_circular_binary(M, m, a):
    masses = [M.value_in(units.MSun), m.value_in(units.MSun)] | units.MSun
    orbital_period = (4 * numpy.pi**2 * a**3 / 
        (constants.G * masses.sum())).sqrt().as_quantity_in(units.day)
    
    print "   Initializing inner binary"
    print "   Orbital period inner binary:", orbital_period
    stars =  Particles(2)
    stars.mass = masses
    stars.position = [0.0, 0.0, 0.0] | units.AU
    stars.velocity = [0.0, 0.0, 0.0] | units.km / units.s
    stars[0].x = a
    stars[0].vy = get_relative_velocity(stars.total_mass(), a, 0.0)
    stars.move_to_center()
    return stars

# from Hills 1983 (Eq. 1) 1983ApJ...267..322H
def post_supernova_semimajor_axis(M0, m0, M1, m1, a0, r0):
    dM = M0-M1
    a = 0.5 * (M0+m0 - dM)/(0.5*(M0+m0) - (a0/r0)*dM)
    return a

# from Hills 1983 (Eq. 6) 1983ApJ...267..322H
def post_supernova_eccentricity(M0, m0, M1, m1, a0, r0):
    dM = M0-M1
    e0 = 0.0
    e = numpy.sqrt(1 - (1-e0**2)
                   * ( (1-(2*a0/r0)*(dM/(M0+m0)))) / ((1-dM/(m0+m0))**2))
    return e


from amuse.ext.orbital_elements import orbital_elements_from_binary
def supernova_in_binary_nbody(M0, m0, a0, tsn):
    stars = make_circular_binary(M0, m0, a0)
    M, m, a, e, ta_out, inc, lan_out, aop_out = orbital_elements_from_binary(stars, G=constants.G)
    print "Initial binary: a=", a.in_(units.AU), "e=", e, "M=", stars[0].mass, "and m=", stars[1].mass

    converter = nbody_system.nbody_to_si(M+m, a)
    gravity = Hermite(converter)
    gravity.particles.add_particles(stars)

    print "Integrate binary to t=", tsn.in_(units.day)
    gravity.evolve_model(tsn)
    M, m, a, e, ta_out, inc, lan_out, aop_out = orbital_elements_from_binary(stars, G=constants.G)
    print "Pre supernova orbit: a=", a.in_(units.AU), "e=", e
    stars[0].mass *= 0.1
    print "Reduce stellar mass to: M=", stars[0].mass, "and m=", stars[1].mass

    v_kick = (0, 0, 0) | units.kms
    stars[0].velocity += v_kick
    gravity.evolve_model(2*tsn)
    M, m, a, e, ta_out, inc, lan_out, aop_out = orbital_elements_from_binary(stars, G=constants.G)
    print "Post supernova orbit: a=", a.in_(units.AU), "e=", e

    r0 = a
    a_ana = post_supernova_semimajor_axis(M0, m0, stars[0].mass, stars[1].mass, a, r0)
    e_ana = post_supernova_eccentricity(M0, m0, stars[0].mass, stars[1].mass, a, r0)
    print "Analytic solution to post orbit orbital parameters: a=", a_ana, "e=", e_ana
    gravity.stop()
    
if __name__ in ('__main__'):
    M = 10|units.MSun
    m = 10|units.MSun
    a = 1 | units.AU
    orbital_period = (4 * numpy.pi**2 * a**3 / (constants.G * (M+m))).sqrt()
    tsn = 10*orbital_period
    supernova_in_binary_nbody(M, m, a, tsn)

