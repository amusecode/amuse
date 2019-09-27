"""
   Visualization for simple N-body integration.
   Reads particle set from file (nbody.hdf5) and prints subsequent frames.
"""
import os
import sys
import numpy

import matplotlib
#matplotlib.use('Agg')
from matplotlib import pyplot
#from matplotlib import pyplot
from amuse.plot import scatter, xlabel, ylabel
from amuse.lab import *
from amuse.io import store
from optparse import OptionParser
from amuse.units.optparse import OptionParser
from time import sleep
from amuse.ext.orbital_elements import orbital_elements_from_binary

def read_planetary_system(filename): #lim, snapshot_id):
    source = Particles(0)
    gas = Particles(0)
    time = 0 | units.yr
    bodies = read_set_from_file(filename, "amuse")
    i = 0
    for bi in bodies.history:
        i+=1
        if i%2:
            source = bi
        else:
            gas = bi
            break
    print "N=", len(source), len(gas)
    return source, gas

def mu(X = None, Y = 0.25, Z = 0.02, x_ion = 0.1):
    """
    Compute the mean molecular weight in kg (the average weight of particles in a gas)
    X, Y, and Z are the mass fractions of Hydrogen, of Helium, and of metals, respectively.
    x_ion is the ionisation fraction (0 < x_ion < 1), 1 means fully ionised
    """
    if X is None:
        X = 1.0 - Y - Z
    elif abs(X + Y + Z - 1.0) > 1e-6:
        raise Exception("Error in calculating mu: mass fractions do not sum to 1.0")
    return constants.proton_mass / (X*(1.0+x_ion) + Y*(1.0+2.0*x_ion)/4.0 + Z*x_ion/2.0)

def plot_single_image(source, disk):
    lim = 150#|units.AU
    index = 4
    
    left, width = 0.1, 0.65
    bottom, height = 0.1, 0.65
    bottom_h = left_h = left+width+0.05
    rect_scatter = [left, bottom, width, height]
    rect_histx = [left, bottom_h, width, 0.2]
    rect_histy = [left_h, bottom, 0.2, height]

    fig = pyplot.figure(figsize=(16,14))	
    time = disk.get_timestamp()
    xy = pyplot.axes(rect_scatter)
    xz = pyplot.axes(rect_histx)
    yz = pyplot.axes(rect_histy)
    xy.set_xlabel("X [AU]")
    xy.set_ylabel("Y [AU]")
    xz.set_ylabel("Z [AU]")
    yz.set_xlabel("Z [AU]")

    positions = disk.position
    x, y, z = positions.x.value_in(units.AU), positions.y.value_in(units.AU), positions.z.value_in(units.AU)
    vx, vy, vz = disk.vx.value_in(units.kms), disk.vy.value_in(units.kms), disk.vz.value_in(units.kms)
    v = numpy.sqrt(vx**2+vy**2+vz**2)
    v_max = max(v)
    v_min = min(v)

    stars = source
    xs, ys, zs = stars.position.x.value_in(units.AU), stars.position.y.value_in(units.AU), stars.position.z.value_in(units.AU)

#    sizes = 1000
#    sizes = 1000*disk.rho/disk.rho.max()
    alpha = 0.2
    us = disk.u
    xion = disk.xion
    print xion
    u_min, u_max = min(us), max(us)
    xion_min, xion_max = min(xion), max(xion)
    print "u=", u_min, u_max
    Ts = mu() / constants.kB * disk.u
    T_min = max(20|units.K, mu() / constants.kB * u_min)
    T_max = min(1800|units.K, mu() / constants.kB * u_max)
    print "T=", T_min, T_max
    print "X=", xion_min, xion_max

    log_u = numpy.log((us / u_min)) / numpy.log((u_max / u_min))
    clipped_log_u = numpy.minimum(numpy.ones_like(log_u), numpy.maximum(numpy.zeros_like(log_u), log_u))

    log_v = numpy.log((v / v_min)) / numpy.log((v_max / v_min))
    clipped_log_v = numpy.minimum(numpy.ones_like(log_v), numpy.maximum(numpy.zeros_like(log_v), log_v))


    log_T = numpy.log((Ts / T_min)) / numpy.log((T_max / T_min))
    clipped_log_T = numpy.minimum(numpy.ones_like(log_T), numpy.maximum(numpy.zeros_like(log_T), log_T))
    
    ps = disk.rho
    p_min, p_max = min(ps), max(ps)
    log_p = numpy.log((ps / p_min)) / numpy.log((p_max / p_min))
    clipped_log_p = numpy.minimum(numpy.ones_like(log_p), numpy.maximum(numpy.zeros_like(log_p), log_p))

    blue = clipped_log_T
    red  = clipped_log_p
    green = clipped_log_v
    colors = numpy.transpose(numpy.array([red, green, blue]))

    sizes = 2*2000*disk.h_smooth/disk.h_smooth.max()
    xy.scatter(x, y, sizes, c=colors, edgecolors = "none", alpha = alpha)
    s = 50*stars.mass/stars.mass.min()
    xy.scatter(xs, ys, s, c='b', edgecolors = "none", alpha = 1)
    
    xy.set_xlim( (-lim, lim) )
    xy.set_ylim( (-lim, lim) )
    sizes = 20*200*disk.h_smooth/disk.h_smooth.max()
    xz.scatter(x, z, sizes, c=colors, edgecolors = "none", alpha = alpha)
    xz.scatter(xs, zs, s, c='b', edgecolors = "none", alpha = 1)

    yz.scatter(z, y, sizes, c=colors, edgecolors = "none", alpha = alpha)
    yz.scatter(zs, ys, s, c='b', edgecolors = "none", alpha = 1)

    xz.set_xlim( xy.get_xlim() )
    yz.set_ylim( xy.get_xlim() )
    yz.set_xlim( (-0.2*lim, 0.2*lim) )
    xz.set_ylim( (-0.2*lim, 0.2*lim) )
#    yz.set_xlim( (-lim, lim) )
#    xz.set_ylim( (-lim, lim) )

    filename = "fig_hydro_disk_with_bump_i{0:04}".format(index)
    fig.savefig(filename)
#    pyplot.show()

def main(filename=None): #, lim=-1|units.AU, image_id=-1):

    output_single_image(filename)
    """    
    if filename:
        output_single_image(lim, image_id)
    else:
        output_multiple_images(lim)
    """    

def calculate_orbital_elements(star, planet):
    from amuse.ext.orbital_elements import orbital_elements_from_binary

    p = Particles()
    p.add_particle(star)
    p.add_particle(planet)
    M, m, a, e, ta_out, inc, lan_out, aop_out = orbital_elements_from_binary(p, G=constants.G)
    #print "Orbital elements:", M, m, a, e, ta_out, inc, lan_out, aop_out
    return a, e, inc
    
def output_single_image(filename):
    source, disk = read_planetary_system(filename)
    print disk[0].position.in_(units.AU)
    snapshot_id = 1
    lim = 150|units.AU
    plot_single_image(source, disk)

def output_multiple_images(lim):
    snapshot_id = 0
    filename = "planetary_system_i{0:04}.amuse".format(snapshot_id)
    while os.path.exists(filename):
        bodies = read_set_from_file(filename, "amuse")
        for bi in bodies.history:
            if len(bi)<=20:
                planets = bi.copy()
                time = bi.get_timestamp()
                print "Orbits at t=", time, planets.semimajor_axis.in_(units.AU), planets.eccentricity
            else:
                disk = bi.copy()

                time = bi.get_timestamp()
                print "Snapshot=", snapshot_id, time
                if image_id<0 or image_id == snapshot_id:
                    plot_single_image(planets, disk)
                if image_id == snapshot_id:
                    print "Stop plotting"
                    break
        snapshot_id += 1
        filename = "planetary_system_i{0:04}.amuse".format(snapshot_id)
            
def new_option_parser():
    result = OptionParser()
    result.add_option("-f", dest="filename", default = "hydro_disk_with_bump_i0004.amuse",
                      help="output filename [%default]")
    return result

if __name__ in ('__main__', '__plot__'):
    o, arguments  = new_option_parser().parse_args()
    main(**o.__dict__)


