import sys
import numpy 
import os

try:
    from matplotlib import pyplot
    HAS_MATPLOTLIB = True
except ImportError:
    HAS_MATPLOTLIB = False

from amuse.support.units import units
from amuse.support.data import core
from amuse.legacy.sse.interface import SSE
from amuse.legacy.evtwin.interface import EVtwin
from amuse.legacy.support.core import is_mpd_running
import path_to_test_results

def simulate_stellar_evolution(mass_in = 5.0 | units.MSun, end_time = 250.0 | units.Myr,
        name_of_the_figure="handover_stellar_evolution.png"):
    """
    A star will be created with the given mass. The stellar evolution will be 
    simulated using EVtwin until it returns an error code. This usually means a 
    star is approaching a troublesome phase in stellar evolution, e.g. the Helium 
    flash. From that point SSE will take over as best it can.
    """
    print "The evolution of a star with mass=",str(mass_in)," will be", \
        "simulated until t=", str(end_time), "..."
    
    stellar_evolution = EVtwin()
    stellar_evolution.initialize_module_with_default_parameters()
        
    print "Initializing the particles..."
    stars = core.Stars(1)
    stars.add_calculated_attribute("temperature",calculate_effective_temperature)
    star = stars[0]
    star.mass = mass_in
    star.radius = 0.0 | units.RSun
    stellar_evolution.setup_particles(stars)
    stellar_evolution.initialize_stars()
    from_stellar_evolution_to_model = stellar_evolution.particles.new_channel_to(stars)
    from_stellar_evolution_to_model.copy()

    memory=[]
    age_at_code_switch = -1.0 | units.Myr
    current_time = 0.0 | units.Myr
    no_exception = True
    print "Start evolving with EVtwin..."
    while current_time < end_time and no_exception:
        try:
            stellar_evolution.evolve_model()
        except Exception as exception_message:
            no_exception=False
            print exception_message
        else:
            from_stellar_evolution_to_model.copy()
            memory.append((star.age, star.mass, star.radius, star.temperature, star.luminosity))
            current_time = star.age
    
    if no_exception:
        print "The simulation has been performed using EVtwin only - no exceptions occured."
    else:
        stellar_evolution.stop()
        age_at_code_switch = current_time
        print "Exception in EVtwin interface, continue using SSE from t=", \
            str(age_at_code_switch)
        stellar_evolution = SSE()
        stellar_evolution.initialize_module_with_default_parameters()
        stellar_evolution.setup_particles(stars)
        stellar_evolution.initialize_stars()
        from_stellar_evolution_to_model = stellar_evolution.particles.new_channel_to(stars)
        from_stellar_evolution_to_model.copy()
        while current_time < end_time:
            stellar_evolution.evolve_model()
            from_stellar_evolution_to_model.copy()
            memory.append((star.age, star.mass, star.radius,
                star.temperature, star.luminosity))
            current_time = star.age
        
    print "Evolved model successfully."
    stellar_evolution.stop()
    
    print "Starting the reference run (SSE only) for comparison..."
    stellar_evolution = SSE()
    stellar_evolution.initialize_module_with_default_parameters()
    stars2 = core.Stars(1)
    stars2.add_calculated_attribute("temperature",calculate_effective_temperature)
    star2 = stars2[0]
    star2.mass = mass_in
    star2.radius = 0.0 | units.RSun
    stellar_evolution.setup_particles(stars2)
    stellar_evolution.initialize_stars()
    from_stellar_evolution_to_model = stellar_evolution.particles.new_channel_to(stars2)
    from_stellar_evolution_to_model.copy()
    refmemory=[]
    current_time = 0.0 | units.Myr
    while current_time < end_time:
        current_time = current_time + (1.0 | units.Myr)
        stellar_evolution.evolve_model(current_time)
        from_stellar_evolution_to_model.copy()
        refmemory.append((star2.age, star2.mass, star2.radius,
            star2.temperature, star2.luminosity))
    
    print "Reference run (SSE only) finished."
    stellar_evolution.print_refs()
    stellar_evolution.stop()

    create_plots(memory, age_at_code_switch, end_time, name_of_the_figure, refmemory)
    print
    print "All done!"

def create_plots(memory, age_at_code_switch, end_time, name_of_the_figure, refmemory):
    if HAS_MATPLOTLIB:
        # An obscure but elegant way to transpose the list of tuples:
        transposed_memory = map(list, zip(*memory))
        n_plots = len(transposed_memory)-1
        transposed_refmemory = map(list, zip(*refmemory))

        print "Plotting the data..."
        figure = pyplot.figure(figsize = (8, 20))
        plots = map(lambda x : figure.add_subplot(n_plots,1,n_plots-x), range(n_plots))
        plots.reverse()
        if age_at_code_switch > 0.0 | units.Myr:
            pyplot.suptitle('Stellar evolution was handed over from EVtwin to SSE at t='+ \
                str(age_at_code_switch))
        else:
            pyplot.suptitle('Stellar evolution was simulated using EVtwin only.', fontsize=16)
        pyplot.title('SSE-only simulation shown for comparison (with offset, magenta).', fontsize=12, color="m")
        
        x_values = [x.value_in(units.Myr) for x in transposed_memory[0]]
        ref_x_values = [x.value_in(units.Myr) for x in transposed_refmemory[0]]
        y_labels = [r'mass (M$_\odot$)',r'radius (R$_\odot$)',r'T$_{\rm eff}$ (K)', \
            r'luminosity (L$_\odot$)']
        y_units = [units.MSun,units.RSun,units.K,units.LSun]
        for i, y_values in enumerate(transposed_memory):
            if i==0: continue
            y_values = [y.value_in(y_units[i-1]) for y in y_values]
            plots[i-1].semilogy(x_values, y_values, "b-")
            f_offset=pow((plots[i-1].axis()[3])/(plots[i-1].axis()[2]),0.1)
            f_offset = 3.0
            ref_y_values = [f_offset*(y.value_in(y_units[i-1])) for y in transposed_refmemory[i]]
            plots[i-1].semilogy(ref_x_values, ref_y_values, "m-")
            
            # Mark the time of the code switch:
            x_mark = age_at_code_switch.value_in(units.Myr)
            plots[i-1].semilogy([x_mark,x_mark], [plots[i-1].axis()[2],plots[i-1].axis()[3]], "k--")
            plots[i-1].axis(xmin=0.0, xmax=end_time.value_in(units.Myr))
            plots[i-1].set_ylabel(y_labels[i-1])
            if not i==n_plots:
                pyplot.setp(plots[i-1].get_xticklabels(), visible=False)

        plots[n_plots-1].set_xlabel('star age (Myr)')
        pyplot.savefig(name_of_the_figure)

def calculate_effective_temperature(luminosity,radius):
    Stefan_Boltzmann_constant = 5.670400e-8 | units.J * units.s**-1 * units.m**-2 * units.K**-4
    return ((luminosity/(4*numpy.pi*Stefan_Boltzmann_constant*radius**2))**.25).in_(units.K)

def test_simulate_short():
    assert is_mpd_running()
    test_results_path = path_to_test_results.get_path_to_test_results()
    output_file = os.path.join(test_results_path, "handover_stellar_evolution.png")
    simulate_stellar_evolution(end_time = 2.0 | units.Myr, \
        name_of_the_figure=output_file)
    
if __name__ == '__main__':
    if len(sys.argv) == 1:
        simulate_stellar_evolution()
    elif len(sys.argv) == 2:
        simulate_stellar_evolution(name_of_the_figure = sys.argv[1])
    elif len(sys.argv) == 3:
        simulate_stellar_evolution(mass_in=(float(sys.argv[1]) | units.MSun),
            end_time=(float(sys.argv[2]) | units.Myr))
    elif len(sys.argv) == 4:
        simulate_stellar_evolution(name_of_the_figure = sys.argv[1],
            mass_in=(float(sys.argv[2]) | units.MSun),
            end_time=(float(sys.argv[3]) | units.Myr))
    else:
        print "I'm confused: too many arguments!"
