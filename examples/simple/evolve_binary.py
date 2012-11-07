"""
Evolves a stellar binary and reports the mass of each star during the evolution.

Shows the type of each star as they change.
"""
from amuse.lab import *
import numpy
from matplotlib import pyplot

def evolve_binary(mass_of_star1, mass_of_star2, orbital_period, eccentricity):
    code = SeBa()
    code.commit_parameters()
    
    stars =  Particles(2)
    stars[0].mass = mass_of_star1
    stars[1].mass = mass_of_star2
    
    
    mu = stars.mass.sum() * constants.G
    semi_major_axis = (((orbital_period / (2.0 * numpy.pi))**2)*mu)**(1.0/3.0)
    
    code.particles.add_particles(stars)
    
    binaries =  Particles(1)
    
    binary = binaries[0]
    binary.semi_major_axis = semi_major_axis
    binary.eccentricity = eccentricity
    binary.child1 = stars[0]
    binary.child2 = stars[1]
    
    code.binaries.add_particles(binaries)
    
    from_seba_to_model = code.particles.new_channel_to(stars)
    from_seba_to_model.copy()

    from_seba_to_model_binaries = code.binaries.new_channel_to(binaries)
    from_seba_to_model_binaries.copy()
    
    previous_type_child1 = binary.child1.stellar_type
    previous_type_child2 = binary.child2.stellar_type
    
    results = []
    current_time = 0 | units.Myr
    while current_time < (1000 | units.Myr):
        code.update_time_steps()
        # The next line appears a bit weird, but saves time for this simple test.
        deltat = max(1.0*code.binaries[0].time_step, 0.1 | units.Myr)
        current_time = current_time + deltat
        code.evolve_model(current_time)
        from_seba_to_model.copy()
        from_seba_to_model_binaries.copy()
        
        if not binary.child1.stellar_type == previous_type_child1:
            print binary.age, "Child 1, change of stellar type", previous_type_child1, ' -> ',binary.child1.stellar_type
            previous_type_child1 = binary.child1.stellar_type
        if not binary.child2.stellar_type == previous_type_child2:
            print binary.age, "Child 2, change of stellar type", previous_type_child2, ' -> ',binary.child2.stellar_type
            previous_type_child2 = binary.child2.stellar_type
        results.append((binary.age, binary.child1.mass, binary.child1.stellar_type, binary.child2.mass, binary.child2.stellar_type))
        
        
    
    code.stop() 
    return results
    

def plot_masses(table):
    time = [] 
    mass_child1 = []
    mass_child2 = []
    
    for age, mass1, type1, mass2, type2 in table:
        time.append(age.value_in(units.Myr))
        mass_child1.append(mass1.value_in(units.MSun))
        mass_child2.append(mass2.value_in(units.MSun))
        
    
    figure = pyplot.figure(figsize = (8, 6))
    plot = figure.add_subplot(1,1,1)
    plot.plot(time, mass_child1)
    plot.plot(time, mass_child2)
    plot.set_xlabel('Time')
    plot.set_ylabel('Mass')
   
    pyplot.show()
    
if __name__ in ('__main__', '__plot__'):        
    table = evolve_binary(
        3.0 | units.MSun,
        0.3 | units.MSun,
        200 | units.day,
        0.5
    )
    plot_masses(table)

