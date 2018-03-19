"""
   Visualization for simple N-body integration.
   Reads particle set from file (nbody.hdf5) and prints subsequent frames.
"""
import os
import sys
import numpy

import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot
#from matplotlib import pyplot
from amuse.plot import scatter, xlabel, ylabel
from amuse.lab import *
from amuse.io import store
from optparse import OptionParser
from amuse.units.optparse import OptionParser
from time import sleep
from amuse.ext.orbital_elements import orbital_elements_from_binary

def plot_single_image(planets, disk, lim, index):

    #centered on the Sun
    com = planets[0].position
    planets.position -= com
    disk.position -= com
    
    left, width = 0.1, 0.65
    bottom, height = 0.1, 0.65
    bottom_h = left_h = left+width+0.05
    rect_scatter = [left, bottom, width, height]
    rect_histx = [left, bottom_h, width, 0.2]
    rect_histy = [left_h, bottom, 0.2, height]

    fig = pyplot.figure(figsize=(12,12))	
    time = disk.get_timestamp()
    # pyplot.title("Cluster at t="+str(time.in_(units.Gyr)))
    xy = pyplot.axes(rect_scatter)
    xy.text(110,110, "protoplanetary disk (img#"+str(index)+")",
            ha='left', va='bottom')
    xz = pyplot.axes(rect_histx)
    yz = pyplot.axes(rect_histy)
    xy.set_xlabel("X [AU]")
    xy.set_ylabel("Y [AU]")
    xz.set_ylabel("Z [AU]")
    yz.set_xlabel("Z [AU]")

    positions = disk.position
    x, y, z = positions.x.value_in(units.AU), positions.y.value_in(units.AU), \
              positions.z.value_in(units.AU)

    alpha = 0.01
    us        = disk.u
    u_min, u_max = min(us), max(us)
    print("u=", u_min, u_max)
    log_u = numpy.log((us / u_min)) / numpy.log((u_max / u_min))
    clipped_log_u = numpy.minimum(numpy.ones_like(log_u),
                                  numpy.maximum(numpy.zeros_like(log_u), log_u))

    ps = disk.rho
    p_min, p_max = min(ps), max(ps)
    log_p = numpy.log((ps / p_min)) / numpy.log((p_max / p_min))
    clipped_log_p = numpy.minimum(numpy.ones_like(log_p),
                                  numpy.maximum(numpy.zeros_like(log_p), log_p))

    red = 1.0 - clipped_log_u
    blue = clipped_log_u
    green = clipped_log_p
    colors = numpy.transpose(numpy.array([red, green, blue]))

    sizes = 2000*disk.h_smooth/disk.h_smooth.max()
    xy.scatter(x, y, sizes, c=colors, edgecolors = "none", alpha = alpha)
    xy.set_xlim( (-lim, lim) )
    xy.set_ylim( (-lim, lim) )
    sizes = 100*disk.h_smooth/disk.h_smooth.max()
    xz.scatter(x, z, sizes, c=colors, edgecolors = "none", alpha = alpha)
    yz.scatter(z, y, sizes, c=colors, edgecolors = "none", alpha = alpha)
    xz.set_xlim( xy.get_xlim() )
    yz.set_ylim( xy.get_xlim() )
    yz.set_xlim( (-0.05*lim, 0.05*lim) )
    xz.set_ylim( (-0.05*lim, 0.05*lim) )

    from distinct_colours import get_distinct
    c = get_distinct(len(planets))
#    m = 1000 * planets.mass/planets.mass.max()
    m = 100 * numpy.log10(1+planets.mass/planets.mass.min())
    #m[0] = min(10*m[1:].max(), 30)
    xy.scatter(planets.x.value_in(units.AU), planets.y.value_in(units.AU),
               s=m, c=c, lw=0) 
    xz.scatter(planets.x.value_in(units.AU), planets.z.value_in(units.AU),
               s=m, c=c, lw=0) 
    yz.scatter(planets.z.value_in(units.AU), planets.y.value_in(units.AU),
               s=m, c=c, lw=0) 

    filename = "planetary_system_i{0:04}.png".format(index)
    fig.savefig(filename)

def XX_main(filename, lim=-1|units.AU, image_id=-1):
    if image_id < 0:
        pyplot.ion()
    star = read_set_from_file("earlysolarsystem.amuse", "amuse")
    planet = read_set_from_file("earlysolarsystem.amuse", "amuse")
    disk = read_set_from_file("earlysolarsystem.amuse", "amuse")

    mc = read_set_from_file(filename, "hdf5") # particles
    time = 0
    snapshot_id = 0
    for si, ssi in zip(mc.history, star.history):
        snapshot_id += 1
        time = si.get_timestamp()
        print("Snapshot=", snapshot_id, time)
        if image_id < 0 or image_id == snapshot_id:
            m = 1
            plot_single_image(si, ssi, lim.value_in(units.AU), snapshot_id)
        if image_id == snapshot_id:
            print("Stop plotting")
            break
        if image_id<0:
            pyplot.draw()
            pyplot.cla()

def calculate_orbital_elements(star, planet):
    p = Particles()
    p.add_particle(star)
    p.add_particle(planet)
    M, m, a, e, ta_out, inc_out, lan_out, aop_out \
        = orbital_elements_from_binary(p, G=constants.G)
    return a, e

def main(filename, lim=-1|units.AU, image_id=-1):
    if image_id < 0:
        output_multiple_images(lim)
    else:
        output_single_image(lim, image_id)

def output_single_image(lim, snapshot_id):
    filename = "planetary_system_i{0:04}.amuse".format(snapshot_id)
    bodies = read_set_from_file(filename, "amuse")
    for bi in bodies.history:
        if len(bi)<=20:
            planets = bi.copy()
            time = bi.get_timestamp()
            print("Orbits at t=", time, planets.semimajor_axis.in_(units.AU), \
                  planets.eccentricity)
        else:
            disk = bi.copy()
            time = bi.get_timestamp()
            print("Snapshot=", snapshot_id, time)
            plot_single_image(planets, disk, lim.value_in(units.AU),
                              snapshot_id)

def output_multiple_images(lim):
    snapshot_id = 0
    filename = "planetary_system_i{0:04}.amuse".format(snapshot_id)
    while os.path.exists(filename):
        bodies = read_set_from_file(filename, "amuse")
        for bi in bodies.history:
            if len(bi) <= 20:
                planets = bi.copy()
                time = bi.get_timestamp()
                print("Orbits at t=", time, \
                    planets.semimajor_axis.in_(units.AU), planets.eccentricity)
            else:
                disk = bi.copy()

                time = bi.get_timestamp()
                print("Snapshot=", snapshot_id, time)
                if image_id < 0 or image_id == snapshot_id:
                    plot_single_image(planets, disk, lim.value_in(units.AU),
                                      snapshot_id)
                if image_id == snapshot_id:
                    print("Stop plotting")
                    break
        snapshot_id += 1
        filename = "planetary_system_i{0:04}.amuse".format(snapshot_id)
            
def new_option_parser():
    result = OptionParser()
    result.add_option("-f", dest="filename", default = "hydro.hdf5",
                      help="output filename [%default]")
    result.add_option("-l", unit=units.AU, 
                      dest="lim", type="float", default = 100|units.AU,
                      help="axis length [%default]")
    result.add_option("-i", 
                      dest="image_id", type="int", default = 281,
                      help="image id [%default]")
    return result

if __name__ in ('__main__', '__plot__'):
    o, arguments  = new_option_parser().parse_args()
    main(**o.__dict__)


