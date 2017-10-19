import numpy
from amuse.units import nbody_system
from amuse.ic.plummer import new_plummer_model
from amuse.community.ph4.interface import ph4
from amuse.community.smalln.interface import SmallN
from amuse.community.kepler.interface import Kepler
from amuse.couple import multiples

# Awkward syntax here because multiples needs a function that resets
# and returns a small-N integrator.

SMALLN = None
def init_smalln():
    global SMALLN
    SMALLN = SmallN()

def new_smalln():
    SMALLN.reset()
    return SMALLN

def stop_smalln():
    global SMALLN
    SMALLN.stop()

def print_diagnostics(grav, E0=None):

    # Simple diagnostics.

    ke = grav.kinetic_energy
    pe = grav.potential_energy
    Nmul, Nbin, Emul = grav.get_total_multiple_energy()
    print ''
    print 'Time =', grav.get_time()
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

def integrate_system(N, t_end, seed=None):

    gravity = ph4()
    gravity.initialize_code()
    gravity.parameters.set_defaults()
    
    if seed is not None: numpy.random.seed(seed)
    stars = new_plummer_model(N)
    stars.mass = 1./N | nbody_system.mass
    stars.scale_to_standard(smoothing_length_squared
                             = gravity.parameters.epsilon_squared)

    id = numpy.arange(N)
    stars.id = id+1                       # used in multiples bookkeeping

    # Set dynamical radii for encounters.

    stars.radius = 0.5*stars.mass.number | nbody_system.length

    gravity.particles.add_particles(stars)

    stopping_condition = gravity.stopping_conditions.collision_detection
    stopping_condition.enable()

    # Combine modules into the multiples code.    

    init_smalln()
    kep = Kepler(unit_converter=None)
    kep.initialize_code()
    multiples_code = multiples.Multiples(gravity, new_smalln, kep)
    multiples_code.neighbor_perturbation_limit = 0.05

    multiples_code.global_debug = 1

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

    # Advance the system.

    E0 = print_diagnostics(multiples_code)
    multiples_code.evolve_model(t_end)
    print_diagnostics(multiples_code, E0)

    gravity.stop()
    kep.stop()
    stop_smalln()
    
if __name__ in ('__main__'):
    N = 100
    t_end = 10.0 | nbody_system.time
    integrate_system(N, t_end) #, 42)
