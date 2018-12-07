import numpy
from amuse.lab import *
from amuse.ext.orbital_elements import orbital_elements_from_binary
from amuse.community.fractalcluster.interface import new_fractal_cluster_model
from matplotlib import pyplot

def plot_single_image(groups_of_particles, lim=10):
    left, width = 0.1, 0.4
    bottom, height = 0.1, 0.4
    bottom_h = left_h = left+width+0.05
    rect_xy = [left, bottom, width, height]
    rect_xz = [left, bottom_h, width, 0.4]
    rect_yz = [left_h, bottom, 0.4, height]

    from distinct_colours import get_distinct
    colors = get_distinct(12)

    fig = pyplot.figure(figsize=(10,10))

    xy = pyplot.axes(rect_xy)
    xz = pyplot.axes(rect_xz)
    yz = pyplot.axes(rect_yz)
    xy.set_xlabel("X [pc]")
    xy.set_ylabel("Y [pc]")
    xz.set_ylabel("Z [pc]")
    yz.set_xlabel("Z [pc]")

    i = 0
    for group in groups_of_particles:
        x = group.x.value_in(units.parsec)
        y = group.y.value_in(units.parsec)
        z = group.z.value_in(units.parsec)
        xy.scatter(x, y, lw=0, c=colors[min(11, i)], s=8)
        xz.scatter(x, z, lw=0, c=colors[min(11, i)], s=8)
        yz.scatter(z, y, lw=0, c=colors[min(11, i)], s=8)
        i += 1

    xy.set_xlim((-lim, lim))
    xy.set_ylim((-lim, lim))
    xz.set_xlim(xy.get_xlim())
    yz.set_ylim(xy.get_xlim())
    yz.set_xlim(xy.get_xlim())
    xz.set_ylim(xy.get_xlim())

    save_file = 'FractalClusterHop.png'
    pyplot.savefig(save_file)
    print '\nSaved figure in file', save_file,'\n'
    pyplot.show()

def find_clumps_with_hop(particles, unit_converter):

###BOOKLISTSTART###
    hop = Hop(unit_converter)
    hop.particles.add_particles(particles)
    hop.calculate_densities()

    mean_density = hop.particles.density.mean() 
    hop.parameters.peak_density_threshold = mean_density
    hop.parameters.saddle_density_threshold = 0.99*mean_density
    hop.parameters.outer_density_threshold = 0.01*mean_density

    hop.do_hop()
    result = [x.get_intersecting_subset_in(particles) for x in hop.groups()]
    hop.stop()
###BOOKLISTSTOP###

    return result

def main(N, Rvir, Qvir, Fd, seed):
    numpy.random.seed(seed)
    
    masses = new_kroupa_mass_distribution(N)

    converter = nbody_system.nbody_to_si(masses.sum(), Rvir)
    bodies = new_fractal_cluster_model(N=N, fractal_dimension=Fd,
                                       random_seed=seed,
                                       convert_nbody=converter)
    bodies.mass = masses
    bodies.move_to_center()
    bodies.scale_to_standard(converter, virial_ratio=Qvir)

    clumps = find_clumps_with_hop(bodies, converter)
    
    for clump in clumps:
        clump.scale_to_standard(converter, virial_ratio=0.5)

    plot_single_image(clumps, lim=50)
    
def new_option_parser():
    from amuse.units.optparse import OptionParser
    result = OptionParser()
    result.add_option("-N", dest="N", type="int",default = 2000,
                      help="number of stars [%default]")
    result.add_option("-R", dest="Rvir", type="float",
                      unit=units.parsec, default = 0.5|units.parsec,
                      help="cluser virial radius [%default]")
    result.add_option("-Q", dest="Qvir", type="float",default = 0.5,
                      help="virial ratio [%default]")
    result.add_option("-F", dest="Fd", type="float",default = 1.6,
                      help="fractal dimension [%default]")
    result.add_option("-s", dest="seed", type="int",default = 12345,
                      help="random number seed [%default]")
    return result

if __name__ in ('__main__', '__plot__'):
    o, arguments  = new_option_parser().parse_args()
    main(**o.__dict__)
