from amuse.lab import *
from amuse.units.quantities import as_vector_quantity
from amuse.couple import encounters
from amuse.units import quantities
import numpy
import logging

def new_smalln(converter):
    result = SmallN(converter)
    result.parameters.timestep_parameter = 0.1
    result.parameters.cm_index = 2001
    return result

def new_kepler(converter):
    kepler = Kepler(converter)
    kepler.initialize_code()
    kepler.set_longitudinal_unit_vector(1.0,0.0, 0.0)
    kepler.set_transverse_unit_vector(0.0, 1.0, 0)
    return kepler

def new_binary_orbit(mass1, mass2, semi_major_axis,
                     eccentricity = 0, keyoffset = 1):
    total_mass = mass1 + mass2
    mass_fraction_particle_1 = mass1 / (total_mass)
    
    binary = Particles(2)
    binary[0].mass = mass1
    binary[1].mass = mass2
    
    mu = constants.G * total_mass
    
    velocity_perihelion \
        = numpy.sqrt( mu / semi_major_axis * ((1.0 + eccentricity)
                                               /(1.0 - eccentricity)))
    radius_perihelion = semi_major_axis * (1.0 - eccentricity)
    print(velocity_perihelion)
    
    binary[0].position = ((1.0 - mass_fraction_particle_1) \
                           * radius_perihelion * [1.0,0.0,0.0])
    binary[1].position = -(mass_fraction_particle_1 \
                            * radius_perihelion * [1.0,0.0,0.0])
    
    binary[0].velocity = ((1.0 - mass_fraction_particle_1) \
                           * velocity_perihelion * [0.0,1.0,0.0])
    binary[1].velocity = -(mass_fraction_particle_1 \
                            * velocity_perihelion * [0.0,1.0,0.0])

    return binary

# see Eggleton 2006 Equation 1.6.3 (2006epbm.book.....E)
def random_semimajor_axis_PPE(Mprim, Msec, P_min=10.|units.day,
                              P_max=100.|units.yr):

    Pmax = P_max.value_in(units.day)
    Pmin = P_min.value_in(units.day)
    mpf = (Mprim.value_in(units.MSun)**2.5)/5.e+4
    rnd_max = (Pmax * mpf)**(1./3.3) / (1 + (Pmin * mpf)**(1./3.3))
    rnd_min = (Pmin * mpf)**(1./3.3) / (1 + (Pmax * mpf)**(1./3.3))
    rnd_max = min(rnd_max, 1)
    rnd =numpy.random.uniform(rnd_min, rnd_max, 1)
    Porb = ((rnd/(1.-rnd))**3.3)/mpf | units.day
    Mtot = Mprim + Msec
    a = ((constants.G*Mtot) * (Porb/(2*numpy.pi))**2)**(1./3.)
    return a

def make_secondaries(center_of_masses, Nbin):

    resulting_binaries = Particles()
    singles_in_binaries = Particles()
    binaries = center_of_masses.random_sample(Nbin)
    mmin = center_of_masses.mass.min()
    for bi in binaries:
        mp = bi.mass
        ms = numpy.random.uniform(mmin.value_in(units.MSun),
                                  mp.value_in(units.MSun)) | units.MSun
        a = random_semimajor_axis_PPE(mp, ms)
        e = numpy.sqrt(numpy.random.random())

        nb = new_binary_orbit(mp, ms, a, e) 
        nb.position += bi.position
        nb.velocity += bi.velocity
        nb = singles_in_binaries.add_particles(nb)
        nb.radius = 0.01 * a 

        bi.radius = 3*a 
        binary_particle = bi.copy()
        binary_particle.child1 = nb[0]
        binary_particle.child2 = nb[1]
        binary_particle.semi_major_axis = a
        binary_particle.eccentricity = e
        resulting_binaries.add_particle(binary_particle)

    single_stars = center_of_masses-binaries
    return single_stars, resulting_binaries, singles_in_binaries

def calculate_orbital_elementss(bi, converter):
    kep = new_kepler(converter)
    comp1 = bi.child1
    comp2 = bi.child2
    mass = (comp1.mass + comp2.mass)
    pos = (comp2.position - comp1.position)
    vel = (comp2.velocity - comp1.velocity)
    kep.initialize_from_dyn(mass, pos[0], pos[1], pos[2],
                            vel[0], vel[1], vel[2])
    a,e = kep.get_elements()
    kep.stop()
    return a, e

###BOOKLISTSTART###
def resolve_changed_binaries(stopping_condition, stellar, converter):
    new_binaries = stopping_condition.particles(0)
    for bi in new_binaries:
        print(("add new binary:", bi))
        a, e = calculate_orbital_elementss(bi, converter)
        bi.semi_major_axis = a
        bi.eccentricity = e
        stellar.binaries.add_particle(bi)
        print(("new binary parameters", a, e))
        print(bi)

    lost_binaries = stopping_condition.particles(1)
    for bi in lost_binaries:
        print(("remove old binary:", bi.key))
        stellar.binaries.remove_particle(bi)

    changed_binaries = stopping_condition.particles(2)
    for bi in changed_binaries:
        bs = bi.as_particle_in_set(stellar.binaries)
        a, e = calculate_orbital_elementss(bi, converter)
        bs.semi_major_axis = a
        bs.eccentricity = e
        print(("Modified binary parameters", a, e))
        print(bs)
###BOOKLISTSTOP###

def update_dynamical_binaries_from_stellar(stellar, multiples_code, converter):
    kep = new_kepler(converter)

    # THIS NEEDS TO BE CHECKED!
    print("++++++++++++ THIS NEEDS TO BE CHECKED ++++++++++++++++++++")

    print(("Number of binaries=", len(stellar.binaries)))
    for bi in stellar.binaries:
        bs = bi.as_particle_in_set(multiples_code.binaries)
        total_mass = bi.child1.mass+bi.child2.mass 
        kep.initialize_from_elements(total_mass, bi.semi_major_axis,
                                     bi.eccentricity)
        rel_position = as_vector_quantity(kep.get_separation_vector())
        rel_velocity = as_vector_quantity(kep.get_velocity_vector())
        mu = bi.child1.mass / total_mass 
        bs.child1.position = mu * rel_position 
        bs.child2.position = -(1-mu) * rel_position 
        bs.child1.velocity = mu * rel_velocity
        bs.child2.velocity = -(1-mu) * rel_velocity
        print(("semi_major_axis=", bi.semi_major_axis, total_mass, \
              bi.child1.mass, bi.child2.mass, bi.eccentricity))
    kep.stop()
        
def kira(tend, N, R, Nbin):
    logging.basicConfig(level=logging.ERROR)

    mass = new_salpeter_mass_distribution(N, mass_min=10|units.MSun)
    converter = nbody_system.nbody_to_si(mass.sum(), R)
    code = Hermite(converter)
    stars = new_plummer_model(N, convert_nbody=converter)
    stars.mass = mass
    stars.radius = 0.01/len(stars) | R.unit

    single_stars, binary_stars, singles_in_binaries \
        = make_secondaries(stars, Nbin)
    print(binary_stars)

    stellar = SeBa()
    stellar.particles.add_particles(single_stars)
    stellar.particles.add_particles(singles_in_binaries)
    stellar.binaries.add_particles(binary_stars)
    channel_to_stars = stellar.particles.new_channel_to(stars)

    encounter_code = encounters.HandleEncounter(
        kepler_code =  new_kepler(converter),
        resolve_collision_code = new_smalln(converter),
        interaction_over_code = None,
        G = constants.G
    )
    multiples_code = encounters.Multiples(
        gravity_code = code,
        handle_encounter_code = encounter_code,
        G = constants.G
    )
    multiples_code.particles.add_particles((stars-binary_stars).copy())
    multiples_code.singles_in_binaries.add_particles(singles_in_binaries)
    multiples_code.binaries.add_particles(binary_stars)
    multiples_code.commit_particles()
    channel_from_stars_to_particles \
        = stellar.particles.new_channel_to(multiples_code.particles)

    stopping_condition \
        = multiples_code.stopping_conditions.binaries_change_detection
    stopping_condition.enable()

    from matplotlib import pyplot
    from distinct_colours import get_distinct
    pyplot.rcParams.update({'font.size': 30})
    figure = pyplot.figure(figsize=(12, 9))
    ax = pyplot.gca()
    ax.get_yaxis().get_major_formatter().set_useOffset(False)
    ax.xaxis._autolabelpos = True
    ax.yaxis._autolabelpos = True
    
    color = get_distinct(2)
    pyplot.scatter(numpy.log10(stellar.binaries.semi_major_axis.
                               value_in(units.AU)),
                   stellar.binaries.eccentricity, c=color[0], s=200, lw=0)
    
    t = quantities.linspace(0*tend, tend, 11)
    for ti in t:
        print(("t, Energy=", ti, multiples_code.particles.mass.sum(), \
              multiples_code.get_total_energy()))
        multiples_code.evolve_model(ti)
        print(("at t=", multiples_code.model_time, \
              "Nmultiples:", len(multiples_code.multiples)))

        if stopping_condition.is_set():
            resolve_changed_binaries(stopping_condition, stellar, converter)

        stellar.evolve_model(ti)
        channel_from_stars_to_particles.copy_attributes(["mass", "radius"])
        update_dynamical_binaries_from_stellar(stellar, multiples_code,
                                               converter)

        print(("Lagrangian radii:", \
              multiples_code.all_singles.LagrangianRadii(converter)))
        print(("MC.particles", multiples_code.particles))
        print(("Lagrangian radii:", \
              multiples_code.particles.LagrangianRadii(converter)))
        print(("t, Energy=", ti, multiples_code.get_total_energy()))

    pyplot.scatter(numpy.log10(stellar.binaries.semi_major_axis
                               .value_in(units.AU)),
                   stellar.binaries.eccentricity, c=color[1], lw=0, s=50)
    pyplot.xlabel("$\log_{10}(a/R_\odot)$")
    pyplot.ylabel("eccentricity")

    save_file = 'kira_a_vs_e.pdf'
    pyplot.savefig(save_file)
    print('\nSaved figure in file', save_file,'\n')
    pyplot.show()
        
    stellar.stop()
        
def new_option_parser():
    from amuse.units.optparse import OptionParser
    result = OptionParser()
    result.add_option("-t", unit=units.Myr,
                      dest="tend",type="float", default=10.|units.Myr)
    result.add_option("-R", unit=units.parsec,
                      dest="R",type="float", default=1|units.parsec)
    result.add_option("-N", 
                      dest="N",type="float", default=100)
    result.add_option("--Nbin", 
                      dest="Nbin",type="int", default=50)
    result.add_option("--seed", 
                      dest="seed",type="int", default=2)
    return result

if __name__ == "__main__":
    set_printing_strategy("custom", 
                      preferred_units = [units.MSun, units.parsec, units.Myr], 
                      precision = 4, prefix = "", 
                      separator = " [", suffix = "]")

    options, arguments  = new_option_parser().parse_args()
    if options.seed >= 0:
        numpy.random.seed(options.seed)
        # This is only for random.sample, which apparently does not use numpy
        import random
        random.seed(options.seed)
    kira(options.tend, options.N, options.R, options.Nbin)
