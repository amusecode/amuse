import os.path
from amuse.test.amusetest import get_path_to_results
try:
    from matplotlib import pyplot
    HAS_MATPLOTLIB = True
    from amuse.plot import plot, semilogy, xlabel, ylabel, loglog
except ImportError:
    HAS_MATPLOTLIB = False

from amuse.support.data.core import Particles
from amuse.support.units import units, generic_unit_system, nbody_system, constants
from amuse.support.units.generic_unit_converter import ConvertBetweenGenericAndSiUnits
from amuse.support.exceptions import AmuseException
from amuse.community.mesa.interface import MESA
from amuse.community.gadget2.interface import Gadget2
from amuse.community.fi.interface import Fi
from amuse.ext.star_to_sph import *

def inject_supernova_energy():
    pass

def setup_stellar_evolution_model():
    out_pickle_file = os.path.join(get_path_to_results(), "super_giant_stellar_structure.pkl")
    if os.path.exists(out_pickle_file):
        return out_pickle_file
    
    stellar_evolution = MESA(redirection = "none")
    stars =  Particles(1)
    stars.mass = 10.0 | units.MSun
    stellar_evolution.initialize_module_with_default_parameters() 
    stellar_evolution.particles.add_particles(stars)
    stellar_evolution.initialize_stars()

    print "Evolving a MESA star with mass:", stellar_evolution.particles[0].mass
    try:
        while True:
            stellar_evolution.evolve_model()
    except AmuseException as ex:
        print "Evolved star to", stellar_evolution.particles[0].age
        print "Radius:", stellar_evolution.particles[0].radius
    
    pickle_stellar_model(stellar_evolution.particles[0], out_pickle_file)
    stellar_evolution.stop()
    return out_pickle_file

def run_supernova():
    # options:
    use_hydro_code = Gadget2 # Fi -or- Gadget2
    hydro_code_options = dict() # e.g. dict(use_gl = True)
    number_of_sph_particles = 1000
    t_end = 1.0e6 | units.s
    
    pickle_file = setup_stellar_evolution_model()
    
    print "Creating initial conditions from a MESA stellar evolution model..."
    core, gas_without_core, core_radius = convert_stellar_model_to_SPH(
        None, 
        number_of_sph_particles, 
        seed = 12345,
        pickle_file = pickle_file,
        with_core_particle = True
    )
    if len(core):
        print "Created", len(gas_without_core), "SPH particles and one 'core-particle':\n", core
        print "Setting gravitational smoothing to:", core_radius
    else:
        print "Warning: Only SPH particles created."
    
    inject_supernova_energy()
    
    print "Evolving (SPH) to:", t_end
    n_steps = 100
    
    unit_converter = ConvertBetweenGenericAndSiUnits(1.0 | units.RSun, constants.G, t_end)
    hydro_code = use_hydro_code(unit_converter, redirection = "none", **hydro_code_options)
    
    try:
        hydro_code.parameters.timestep = t_end / n_steps
    except Exception as exc:
        if not "parameter is read-only" in str(exc): raise
    
    hydro_code.parameters.epsilon_squared = core_radius**2
    hydro_code.gas_particles.add_particles(gas_without_core)
    hydro_code.dm_particles.add_particles(core)
    
    times = [] | units.s
    kinetic_energies =   [] | units.J
    potential_energies = [] | units.J
    thermal_energies =   [] | units.J
    for time in [i*t_end/n_steps for i in range(1, n_steps+1)]:
        hydro_code.evolve_model(time)
        times.append(time)
        kinetic_energies.append(   hydro_code.kinetic_energy)
        potential_energies.append( hydro_code.potential_energy)
        thermal_energies.append(   hydro_code.thermal_energy)
    
    energy_plot(times, kinetic_energies, potential_energies, thermal_energies, 
        os.path.join(get_path_to_results(), "supernova_energy_evolution.png"))
    
    hydro_code.stop()
    print "All done!\n"
    

def energy_plot(time, E_kin, E_pot, E_therm, figname):
    if not HAS_MATPLOTLIB:
        return
    pyplot.figure(figsize = (5, 5))
    plot(time, E_kin.as_quantity_in(units.erg), label='E_kin')
    plot(time, E_pot, label='E_pot')
    plot(time, E_therm, label='E_therm')
    plot(time, E_kin+E_pot+E_therm, label='E_total')
    xlabel('Time')
    ylabel('Energy')
    pyplot.legend(loc=3)
    pyplot.savefig(figname)
    print "\nPlot of energy evolution was saved to: ", figname
    pyplot.close()

if __name__ == "__main__":
    print "Test run to mimic a supernova in SPH"
    print
    print "Details:"
    print "First a high-mass star is evolved up to the super giant phase using MESA. " \
        "Then it is converted to SPH particles with the convert_stellar_model_to_SPH " \
        "procedure (with a non-SPH 'core' particle). Finally the internal energies of " \
        "the inner particles are increased, such that the star gains the 10^51 ergs " \
        "released in supernova explosions."
    print
    run_supernova()
