import numpy
from amuse.units import nbody_system, units, constants
from amuse.ic.plummer import new_plummer_model
from amuse.community.ph4.interface import ph4
from amuse.community.smalln.interface import SmallN
from amuse.community.kepler.interface import Kepler
from amuse.couple import multiples

# Awkward syntax here because multiples needs a function that resets
# and returns a small-N integrator.

###BOOKLISTSTART1###
SMALLN = None
def init_smalln(converter):
    global SMALLN
    SMALLN = SmallN(convert_nbody=converter)
def new_smalln():
    SMALLN.reset()
    return SMALLN
###BOOKLISTSTOP1###

def stop_smalln():
    global SMALLN
    SMALLN.stop()

def print_diagnostics(grav, E0=None):

    # Simple diagnostics.

    ke = grav.kinetic_energy
    pe = grav.potential_energy
    Nmul, Nbin, Emul = grav.get_total_multiple_energy()
    print ''
    print 'Time =', grav.get_time().in_(units.Myr)
    print '    top-level kinetic energy =', ke
    print '    top-level potential energy =', pe
    print '    total top-level energy =', ke + pe
    print '   ', Nmul, 'multiples,', 'total energy =', Emul
    E = ke + pe + Emul
    print '    uncorrected total energy =', E
    
    # Apply known corrections.
    
    Etid = grav.multiples_external_tidal_correction \
            + grav.multiples_internal_tidal_correction  # tidal error
    Eerr = grav.multiples_integration_energy_error	# integration error

    E -= Etid + Eerr
    print '    corrected total energy =', E

    if E0 is not None: print '    relative energy error=', (E-E0)/E0
    
    return E

###BOOKLISTSTART2###
def integrate_system(N, t_end, seed=None):

    total_mass = N|units.MSun
    length = 1|units.parsec
    converter = nbody_system.nbody_to_si(total_mass, length)
    gravity = ph4(convert_nbody=converter)
    gravity.initialize_code()
    gravity.parameters.set_defaults()
    gravity.parameters.epsilon_squared = (0.0|units.parsec)**2

    if seed is not None: numpy.random.seed(seed)
    stars = new_plummer_model(N, convert_nbody=converter)
    stars.mass = total_mass/N
    stars.scale_to_standard(convert_nbody=converter,
                            smoothing_length_squared
                             = gravity.parameters.epsilon_squared)
    id = numpy.arange(N)
    stars.id = id+1
    stars.radius = 0.5/N | units.parsec

    gravity.particles.add_particles(stars)

    stopping_condition = gravity.stopping_conditions.collision_detection
    stopping_condition.enable()

    init_smalln(converter)
    kep = Kepler(unit_converter=converter)
    kep.initialize_code()
    multiples_code = multiples.Multiples(gravity, new_smalln, kep,
                                         constants.G)
    multiples_code.neighbor_perturbation_limit = 0.05
    multiples_code.global_debug = 1
###BOOKLISTSTOP2###

    #	global_debug = 0: no output from multiples
    #	               1: minimal output
    #	               2: debugging output
    #	               3: even more output

    print ''
    print 'multiples_code.neighbor_veto =', \
        multiples_code.neighbor_veto
    print 'multiples_code.neighbor_perturbation_limit =', \
        multiples_code.neighbor_perturbation_limit
    print 'multiples_code.retain_binary_apocenter =', \
        multiples_code.retain_binary_apocenter
    print 'multiples_code.wide_perturbation_limit =', \
        multiples_code.wide_perturbation_limit

    time = numpy.sqrt(length**3/(constants.G*total_mass))
    print '\ntime unit =', time.in_(units.Myr)    
###BOOKLISTSTART3###

    E0 = print_diagnostics(multiples_code)
    multiples_code.evolve_model(t_end)
    print_diagnostics(multiples_code, E0)
###BOOKLISTSTOP3###

    gravity.stop()
    kep.stop()
    stop_smalln()
    
if __name__ in ('__main__'):
    N = 100
    t_end = 10.0 | units.Myr
    integrate_system(N, t_end, 42)
