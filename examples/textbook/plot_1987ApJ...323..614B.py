from amuse.lab import *
#from amuse.plot import sph_particles_plot
from amuse.plot import *

from prepare_figure import *
#from distinct_colours import get_distinct

def get_zones(stars):
    stars.move_to_center()
    from amuse.ext.LagrangianRadii import LagrangianRadii
    massf = [0.25, 0.50, 0.75, 1.0]
    lagrad = LagrangianRadii(stars,  massf = massf)
    for zi in lagrad:
        print "R=", zi.in_(units.RSun)
    center_of_mass = [0,0,0] | units.RSun #initial_stars.center_of_mass()
    zoneA = stars.select(lambda r: (center_of_mass-r).length()<lagrad[0], ["position"])
    zoneB = stars.select(lambda r: (center_of_mass-r).length()<lagrad[1], ["position"])
    zoneC = stars.select(lambda r: (center_of_mass-r).length()<lagrad[2], ["position"])
    zoneC -= zoneB
    zoneB -= zoneA
    zoneD = stars.select(lambda r: (center_of_mass-r).length()>lagrad[2], ["position"])
    print "N=", len(zoneA), len(zoneB), len(zoneC), len(zoneD)
    return (zoneA, zoneB, zoneC, zoneD)

def _get_number_of_particles_in_zone(iz, initial_stars, initial_zones, final_zones):
    range_list = [0,1,2,3]
    range_list.pop(iz)
    print range_list
    fzA = initial_stars.copy()
    for i in range_list:
        fzA -= final_zones[i]  
    nzA = [0,0,0,0]
    for izi in range(len(initial_zones)):
        for fzi in fzA:
            if fzi in initial_zones[izi]:
                nzA[izi] +=1
    print "iz=", iz, nzA
    return nzA

def get_number_of_particles_in_zone(iz, initial_stars, initial_zones, final_star, final_zones):
    range_list = [0,1,2,3]
    range_list.pop(iz)
    print range_list
    fzA = final_zones[iz]
    not_selected = final_star - final_zones[iz]
    print iz, len(initial_stars), len(not_selected)
    initial_subset = initial_stars - not_selected
    print len(initial_subset)

    nzA = [0,0,0,0]
    for iz in range(len(initial_zones)):
        for siz in initial_subset:
            if siz in initial_zones[iz]:
                nzA[iz] +=1

    """
    # Normalize
    import numpy
    nzA_sum = float(numpy.sum(nzA))
    print "sum:", nzA_sum
    for i in range(len(nzA)):
        nzA[i] /= nzA_sum
    print nzA
    """
    
    """
    print nzA
    from matplotlib import pyplot
    pyplot.scatter(final_star.x.value_in(units.RSun), final_star.y.value_in(units.RSun), s=1)
    pyplot.scatter(fzA.x.value_in(units.RSun), fzA.y.value_in(units.RSun))
    pyplot.show()
    """
    print nzA
    return nzA

def find_clumps(particles):

    unit_converter = nbody_system.nbody_to_si(1|units.MSun, 1|units.RSun)
    hop = Hop(unit_converter)
    hop.particles.add_particles(particles)
    hop.calculate_densities()

    mean_densty = hop.particles.density.mean() 
    hop.parameters.peak_density_threshold = mean_densty
    hop.parameters.saddle_density_threshold = 0.99*mean_densty
    hop.parameters.outer_density_threshold = 0.01*mean_densty

    hop.do_hop()
    result = [x.get_intersecting_subset_in(particles) for x in hop.groups()]
    hop.stop()
    return result

def plot_collide_two_stars(time, filename):

    from matplotlib import pyplot
    figure = pyplot.figure(figsize=(14, 8))

    ax1 = figure.add_subplot(111)
    ax1.axis('off')
    if filename:
        process_file(time, filename, figure)
    else:
        n = 0.8
        bars = ax1.bar(range(1, 5), [n,0,0,0], color='None', ecolor='black', hatch='//') +\
               ax1.bar(range(1, 5), [n,0,0,0], color='None', ecolor='black', hatch='\\\\') +\
               ax1.bar(range(1, 5), [0,n,0,0], color='None', ecolor='black', hatch='\\\\') +\
               ax1.bar(range(1, 5), [0,0,n,0], color='None', ecolor='black', hatch='*') +\
               ax1.bar(range(1, 5), [0,0,0,n], color='None', ecolor='black', hatch='.')
        pyplot.savefig("1987ApJ323_614B")

    #pyplot.show()
        

def process_file(time, filename, figure):
    fileA = "Hydro_AM06MSun.h5"
    pstar = read_set_from_file(fileA, format='hdf5')
    fileB = "Hydro_BM06MSun.h5"
    sstar = read_set_from_file(fileB, format='hdf5')
    initial_stars = Particles(0)
    initial_stars.add_particles(pstar)
    initial_stars.add_particles(sstar)

    final_star = Particles(0)
    merger = read_set_from_file(filename, "hdf5")
    for si in merger.history:
        ts = si.get_timestamp()
        print ts, time, len(si)
        if ts >= time: 
            final_star.add_particles(si)            
            break

    """
    clumps = find_clumps(final_star)
    for clump in clumps:
        print "N=", len(clump)
    final_star = clumps[0]
    """
        
    initial_zones = get_zones(initial_stars)
    final_zones = get_zones(final_star)

    nA0 = get_number_of_particles_in_zone(0, initial_stars, initial_zones, final_star, final_zones)
    nA1 = get_number_of_particles_in_zone(1, initial_stars, initial_zones, final_star, final_zones)
    nA2 = get_number_of_particles_in_zone(2, initial_stars, initial_zones, final_star, final_zones)
    nA3 = get_number_of_particles_in_zone(3, initial_stars, initial_zones, final_star, final_zones)


    ax2 = pyplot.gca()
    from operator import add
    bars = ax2.bar(range(1, 5), nA0, bottom=[0,0,0,0], color='None', ecolor='black', hatch='//') +\
           ax2.bar(range(1, 5), nA0, color='None', ecolor='black', hatch='\\\\') +\
           ax2.bar(range(1, 5), nA1, bottom=nA0, color='None', ecolor='black', hatch='\\\\') +\
           ax2.bar(range(1, 5), nA2, bottom=map(add, nA0, nA1), color='None', ecolor='black', hatch='/') +\
           ax2.bar(range(1, 5), nA3, bottom=map(add, nA0, map(add, nA1, nA2)), color='None', ecolor='black', hatch='.')
#    ax2.set_xticks([1, 2, 3, 4])
    
    #pyplot.show()
    pyplot.savefig("1987ApJ323_614B")

def new_option_parser():
    from amuse.units.optparse import OptionParser
    result = OptionParser()
    result.add_option("-f", 
                      dest="filename", default="1987ApJ...323..614B_headon.h5",
                      help="filename [%default]")
    result.add_option("-t",  unit=units.hour,
                      dest="time", type="float", default=7|units.hour,
                      help="plotting time [%default]")
    return result

if __name__ in ('__main__', '__plot__'):
    o, arguments  = new_option_parser().parse_args()
    print o.time
    plot_collide_two_stars(o.time, o.filename)

