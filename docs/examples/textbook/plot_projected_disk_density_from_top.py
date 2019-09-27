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

def plot_single_image(planets, debris, disk, shell, figfile=None):

    if len(shell)>0:
        disk.add_particles(shell)
    lim = 150
    """
    #centered on the Sun
    if len(planets)>1:
        print "Center on COM"
        com = planets[0].position
        vcom = planets[0].velocity
        planets.position -= com
        planets.velocity -= vcom
        disk.position -= com
        disk.velocity -= vcom
        debris.position -= com
        debris.velocity -= vcom
    """
    
    left, width = 0.1, 0.65
    bottom, height = 0.1, 0.65
    bottom_h = left_h = left+width+0.05
    rect_scatter = [left, bottom, width, height]
    rect_histx = [left, bottom_h, width, 0.2]
    rect_histy = [left_h, bottom, 0.2, height]

    fig = pyplot.figure(figsize=(12,12))	
    time = disk.get_timestamp()
#    pyplot.title("Cluster at t="+str(time.in_(units.Gyr)))
    xy = pyplot.axes(rect_scatter)
    #xy.text(110,110, "protoplanetary disk (img#"+str(index)+")",  ha='left', va='bottom')
    xz = pyplot.axes(rect_histx)
    yz = pyplot.axes(rect_histy)
    xy.set_xlabel("X [AU]")
    xy.set_ylabel("Y [AU]")
    xz.set_ylabel("Z [AU]")
    yz.set_xlabel("Z [AU]")

    positions = disk.position
    x, y, z = positions.x.value_in(units.AU), positions.y.value_in(units.AU), positions.z.value_in(units.AU)

#    sizes = 1000
#    sizes = 1000*disk.rho/disk.rho.max()
    alpha = 0.01
    us        = disk.u
    if hasattr(disk, "xion"):
        xs        = disk.xion
    else:
        xs = numpy.zeros(len(disk))
    u_min, u_max = min(us), max(us)
    xion_min, xion_max = min(xs), max(xs)
#    if xion_min<=0:
#        xion_min = xion_max/100.
#        xs += xion_min
#        xion_max += xion_min
    print "u=", u_min, u_max
    Ts = mu() / constants.kB * disk.u
    T_min = max(20|units.K, mu() / constants.kB * u_min)
    T_max = min(1800|units.K, mu() / constants.kB * u_max)
    print "T=", T_min, T_max
    print "X=", xion_min, xion_max

    log_u = numpy.log((us / u_min)) / numpy.log((u_max / u_min))
    clipped_log_u = numpy.minimum(numpy.ones_like(log_u), numpy.maximum(numpy.zeros_like(log_u), log_u))

#    log_x = numpy.log((xs / xion_min)) / numpy.log((xion_max / xion_min))
#    clipped_log_x = numpy.minimum(numpy.ones_like(log_x), numpy.maximum(numpy.zeros_like(log_x), log_x))
#    xrange = (xs/xion_min) / (xion_max/xion_min)
#    print "xrange=", xrange.min(), xrange.max()
    clipped_log_x = xs


    log_T = numpy.log((Ts / T_min)) / numpy.log((T_max / T_min))
    clipped_log_T = numpy.minimum(numpy.ones_like(log_T), numpy.maximum(numpy.zeros_like(log_T), log_T))
    
    ps = disk.rho
    p_min, p_max = min(ps), max(ps)
    log_p = numpy.log((ps / p_min)) / numpy.log((p_max / p_min))
    clipped_log_p = numpy.minimum(numpy.ones_like(log_p), numpy.maximum(numpy.zeros_like(log_p), log_p))

    #    red   = 1 - clipped_log_u**(1./2.)
#    blue  = clipped_log_u**(1./2.)
#    green = numpy.minimum(red, blue)
#    red   = 1.0 - clipped_log_T
#    red   = 1-clipped_log_x**(1./2.)
    blue  = clipped_log_T
    red  = 1-clipped_log_T
    #green = numpy.minimum(red, blue)
    green = clipped_log_p
    colors = numpy.transpose(numpy.array([red, green, blue]))

    sizes = 2*2000*disk.h_smooth/disk.h_smooth.max()
    xy.scatter(x, y, sizes, c=colors, edgecolors = "none", alpha = alpha)
    xy.set_xlim( (-lim, lim) )
    xy.set_ylim( (-lim, lim) )
    sizes = 2*200*disk.h_smooth/disk.h_smooth.max()
    xz.scatter(x, z, sizes, c=colors, edgecolors = "none", alpha = alpha)
    yz.scatter(z, y, sizes, c=colors, edgecolors = "none", alpha = alpha)
    xz.set_xlim( xy.get_xlim() )
    yz.set_ylim( xy.get_xlim() )
    yz.set_xlim( (-0.1*lim, 0.1*lim) )
    xz.set_ylim( (-0.1*lim, 0.1*lim) )
#    yz.set_xlim( (-lim, lim) )
#    xz.set_ylim( (-lim, lim) )

    if len(debris)>0:
        c = 'k'
        m = 0.1
        xy.scatter(debris.x.value_in(units.AU), debris.y.value_in(units.AU), s=m, c=c, lw=0) 
        xz.scatter(debris.x.value_in(units.AU), debris.z.value_in(units.AU), s=m, c=c, lw=0) 
        yz.scatter(debris.z.value_in(units.AU), debris.y.value_in(units.AU), s=m, c=c, lw=0) 

    if len(planets)>0:
        from distinct_colours import get_distinct
        c = get_distinct(len(planets))
        m = 1000 * planets.mass/planets.mass.max()
        m[0] = min(10*m[1:].max(), 30)
        xy.scatter(planets.x.value_in(units.AU), planets.y.value_in(units.AU), s=m, c=c, lw=0) 
        xz.scatter(planets.x.value_in(units.AU), planets.z.value_in(units.AU), s=m, c=c, lw=0) 
        yz.scatter(planets.z.value_in(units.AU), planets.y.value_in(units.AU), s=m, c=c, lw=0) 

    filename = "planetary_system.png"
    fig.savefig(filename)

def myplot(x, y, s, bins=1000):
    from scipy.ndimage.filters import gaussian_filter
    heatmap, xedges, yedges = numpy.histogram2d(x, y, bins=bins)
    heatmap = gaussian_filter(heatmap, sigma=s)
    extent = [xedges[0], xedges[-1], yedges[0], yedges[-1]]
    return heatmap.T, extent

def plot_view(x, y, z, xlim, ylim, levels=4):
    from scipy.interpolate import griddata
    xi = numpy.linspace(-xlim,xlim,100)
    yi = numpy.linspace(-ylim,ylim,100)
#    xi = numpy.linspace(-xlim,xlim,xlim)
#    yi = numpy.linspace(-ylim,ylim,ylim)
    # grid the data.
    print len(x), len(y), len(z)
    zi = griddata((x, y), z, (xi[None,:], yi[:,None]), method='linear')
    # contour the gridded data, plotting dots at the randomly spaced data points.
#    CS = pyplot.contourf(xi,yi,zi,levels,cmap=pyplot.cm.plasma,
#                       anorm=pyplot.Normalize(vmax=0.0, vmin=0.9*abs(z).max()))
    CS = pyplot.contourf(xi,yi,zi,levels,cmap=pyplot.cm.plasma,
                         anorm=pyplot.Normalize(vmax=0.0, vmin=9000), vmin=0, vmax=9000)
    # now determine nice limits by hand:
    pyplot.xlim((-xlim, xlim))
    pyplot.ylim((-ylim, ylim))
    return CS

def Xplot_density_view(ax_xy, x, y, z, xlim, ylim, levels=5):
    import numpy as np
    import matplotlib.pyplot as pl
    import scipy.stats as st

    xmin = -xlim
    xmax = xlim
    ymin = -ylim
    ymax = ylim
    xx, yy = np.mgrid[xmin:xmax:100j, ymin:ymax:100j]
    positions = np.vstack([xx.ravel(), yy.ravel()])
    values = np.vstack([x, y])
    kernel = st.gaussian_kde(values)
    f = np.reshape(kernel(positions).T, xx.shape)

    ax = ax_xy
    ax.set_xlim(xmin, xmax)
    ax.set_ylim(ymin, ymax)
    # Contourf plot
    cfset = ax.contourf(xx, yy, f, cmap='Blues')
    ## Or kernel density estimate plot instead of the contourf plot
    #ax.imshow(np.rot90(f), cmap='Blues', extent=[xmin, xmax, ymin, ymax])
    # Contour plot
    cset = ax.contour(xx, yy, f, colors='k')
    return cset

def plot_temperature_image(planets, debris, disk, shell, figfile="fig_disk_top_view_temperature.pdf"):
    from prepare_figure import single_frame, figure_frame, set_tickmarks

    disk.temperature = mu() / constants.kB * disk.u

    x = disk.x.value_in(units.AU)
    y = disk.y.value_in(units.AU)
    z = disk.z.value_in(units.AU)
    T = disk.temperature.value_in(units.K)
    print "min temp=", T.min()

    import numpy as np
    import matplotlib.pyplot as plt
    from matplotlib.ticker import NullFormatter
    nullfmt = NullFormatter()         # no labels

    if "0000" in figfile:
        fig = plt.figure(1, figsize=(12, 10))
    else:
        fig = plt.figure(1, figsize=(10, 10))
    pyplot.axis('equal')

    view = 10
    img = plot_view(x, y, T, 100, 100, view)

    if "0000" in figfile:
        cbar = pyplot.colorbar(img, orientation='vertical')
        font = matplotlib.font_manager.FontProperties(family='times new roman', size=16)
        cbar.set_label('T [K]', rotation=270)

    pyplot.xlabel("X [au]")
    pyplot.ylabel("Y [au]")

    pyplot.savefig(figfile)
    plt.show()

def plot_density_view(ax_xy, x, y, z, xlim, ylim, levels=5):

    import numpy as np
    import matplotlib.pyplot as pl
    import scipy.stats as st

    """
    import seaborn as sns
    sns.set(color_codes=True)
    ax = sns.kdeplot(x, y, shade=True, ax=ax_xy)

    cmap = sns.cubehelix_palette(as_cmap=True, dark=0, light=1, reverse=True)
    sns.kdeplot(x, y, cmap=cmap, n_levels=60, shade=True);
    return cmap
    """

    xmin = -xlim
    xmax = xlim
    ymin = -ylim
    ymax = ylim
    xx, yy = np.mgrid[xmin:xmax:100j, ymin:ymax:100j]
    positions = np.vstack([xx.ravel(), yy.ravel()])
    values = np.vstack([x, y])
    kernel = st.gaussian_kde(values)
    f = np.reshape(kernel(positions).T, xx.shape)

    ax = ax_xy
    ax.set_xlim(xmin, xmax)
    ax.set_ylim(ymin, ymax)
    # Contourf plot
#    cfset = ax.contourf(xx, yy, f, cmap='Blues')
    cfset = ax.contourf(xx, yy, f, levels, cmap='jet')

#    import matplotlib as mpl
#    norm = mpl.colors.Normalize(vmin=0,vmax=1.)
    
    ## Or kernel density estimate plot instead of the contourf plot
#    ax.imshow(np.rot90(f), cmap='Blues', extent=[xmin, xmax, ymin, ymax],
#              vmin=0, vmax=1)
    ax.imshow(np.rot90(f), cmap='Blues', extent=[xmin, xmax, ymin, ymax],
              norm=matplotlib.colors.LogNorm())
    # Contour plot
    #cset = ax.contour(xx, yy, f, colors='k')

    return cfset


def fmt(x, pos):
    a, b = '{:.2e}'.format(x).split('e')
    b = int(b)
    return r'${} \times 10^{{{}}}$'.format(a, b)

def plot_projected_density_image(planets, debris, disk, shell, figfile):
    from prepare_figure import single_frame, figure_frame, set_tickmarks

    disk.temperature = mu() / constants.kB * disk.u

    x = disk.x.value_in(units.AU)
    y = disk.y.value_in(units.AU)
    z = disk.z.value_in(units.AU)
    T = disk.temperature.value_in(units.K)

    SLICE_DATA = False
    if SLICE_DATA:
        selection = z < 10
        x = x[selection]
        y = y[selection]
        z = z[selection]
        T = T[selection]
        selection = z > -10
        x = x[selection]
        y = y[selection]
        z = z[selection]
        T = T[selection]
        print "After slicing the data: N=", len(x)
    
    rho = numpy.log10(disk.rho.value_in(units.g/units.cm**3))
    print "Density=", rho.min(), rho.max(), "in log(g/cm^3)"
    #-12.9251289255 -9.9657991654 in log(g/cm^3)

    import numpy as np
    import matplotlib.pyplot as plt
    from matplotlib.ticker import NullFormatter
    nullfmt = NullFormatter()         # no labels

    left, width = 0.2, 0.45
    bottom, height = 0.1, 0.65
    bottom_h = left_h = left + width + 0.02

    rect_xy = [left, bottom, width, height]
    rect_xz = [left, bottom_h+0.11, width+0.114, 0.2]
    rect_zy = [left_h, bottom, 0.15, height]

    # start with a rectangular Figure
    fig = plt.figure(1, figsize=(12, 8))

    ax_xy = plt.axes(rect_xy)
    ax_xz = plt.axes(rect_xz)
    ax_zy = plt.axes(rect_zy)
    # no labels
    ax_xz.xaxis.set_major_formatter(nullfmt)
    ax_zy.yaxis.set_major_formatter(nullfmt)

#    levels = [-16, -14, -13, -11] 
    levels = 20
    plot_density_view(ax_xy, x, y, rho, 100, 100, levels)
    plot_density_view(ax_xz, x, z, rho, 100, 30, levels)
    img = plot_density_view(ax_zy, z, y, rho, 30, 100, levels)

#    import matplotlib.ticker as ticker
#    cbar = pyplot.colorbar(img, ax=ax_xz, orientation='vertical', format=ticker.FuncFormatter(fmt))
#    cbar = pyplot.colorbar(img, ax=ax_xz, orientation='vertical', format='%.0e')

    from matplotlib.ticker import LogFormatter 
    formatter = LogFormatter(10, labelOnlyBase=False)
    bounds = [-4, -3, -2, -1]
    import matplotlib as mpl
    cmap = mpl.cm.jet ##cool
    norm = mpl.colors.BoundaryNorm(bounds, cmap.N)
    cbar = pyplot.colorbar(img, ax=ax_xz, orientation='vertical', format=formatter,
                           norm=norm,
                           ticks=bounds)
#    cbar = pyplot.colorbar(img, ax=ax_xz, orientation='vertical', format=formatter,
#                           ticks=[-4, -3, -2])
#    cb = pyplot.colorbar(ticks=[1,5,10,20,50], format=formatter)

    ax_bar = cbar.ax
    text = ax_bar.yaxis.label
    font = matplotlib.font_manager.FontProperties(family='times new roman', size=16)
    text.set_font_properties(font)
    cbar.set_label(r'$log(\rho) [g/cm^3]$', rotation=270, labelpad=+20)

    ax_xz.set_xlim(ax_xy.get_xlim())
    ax_zy.set_ylim(ax_xy.get_ylim())

    ax_xy.set_xlabel("X [au]")
    ax_xy.set_ylabel("Y [au]")
    ax_xz.set_ylabel("Z [au]")
    ax_zy.set_xlabel("Z [au]")
    
    pyplot.savefig(figfile)
    plt.show()
    

def main(filename=None, pp=False):

    source, planets, debris, disk, shell = read_planetary_system(filename)
    print "N=", len(source), len(planets), len(debris), len(disk), len(shell)
    figfile = filename.split(".amuse")[0] + ".pdf"
    if pp:
        #plot_temperature_image(planets, debris, disk, shell, figfile)
        plot_projected_density_image(planets, debris, disk, shell, figfile)
        #plot_single_image_simple(planets, debris, disk, shell)
    else:
        plot_single_image(planets, debris, disk, shell)

def calculate_orbital_elements(star, planet):
    from amuse.ext.orbital_elements import orbital_elements_from_binary

    p = Particles()
    p.add_particle(star)
    p.add_particle(planet)
    M, m, a, e, ta_out, inc, lan_out, aop_out = orbital_elements_from_binary(p, G=constants.G)
    #print "Orbital elements:", M, m, a, e, ta_out, inc, lan_out, aop_out
    return a, e, inc

def read_planetary_system(filename): #lim, snapshot_id):
    planets = Particles(0)
    debris = Particles(0)
    source = Particles(0)
    gas = Particles(0)
    shell = Particles(0)
    time = 0 | units.yr
    bodies = read_set_from_file(filename, "amuse")
    for bi in bodies.history:
      print bi
      if len(bi)>0:
        if "gas" in bi.name:
            gas.add_particles(bi.copy())
        elif "debris" in bi.name:
            debris.add_particles(bi.copy())
        elif "Jupiter" in bi.name:
            planets.add_particles(bi.copy())
        elif "shell" in bi.name:
            shell.add_particles(bi.copy())
        else:
            source.add_particles(bi.copy())
    print len(gas), len(shell)
#    print gas[1]
#    print shell[1]
            
#    shell = gas.select(lambda n: "shell" in n,["name"])
#    gas -= shell
#    print shell[1]
    print "N=", len(source), len(planets), len(debris), len(gas), len(shell)
    print "Read planetary system at time", time.in_(units.yr)
    return source, planets, debris, gas, shell

def output_single_image(filename): #lim, snapshot_id):
    bodies = read_set_from_file(filename, "amuse")
    planets = Particles(0)
    debris = Particles(0)
    for bi in bodies.history:
      if len(bi)>1:
        if "debris" in bi.name:
            debris = bi.copy()
            #print "Orbits at t=", time, debris.name, debris.semimajor_axis.in_(units.AU), debris.eccentricity
        elif "gas" in bi.name:
            disk = bi.copy()
            lim = 150|units.AU
            snapshot_id = 1
#            for pi in planets[1:]:
#                a, e, inc = calculate_orbital_elements(planets[0], pi)
#                print "Planet orbits:", a.in_(units.AU), e, inc
            
        else:
            planets = bi.copy()
            time = bi.get_timestamp()
            print "Orbits at t=", time, planets.name, planets.semimajor_axis.in_(units.AU), planets.eccentricity

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
                    plot_single_image(planets, disk, lim.value_in(units.AU), snapshot_id)
                if image_id == snapshot_id:
                    print "Stop plotting"
                    break
        snapshot_id += 1
        filename = "planetary_system_i{0:04}.amuse".format(snapshot_id)
            
def new_option_parser():
    result = OptionParser()
    result.add_option("-f", dest="filename", default = "ID_Nd5Nr7Ms23Rs02pcff-2_i0009.amuse",
                      help="output filename [%default]")
    result.add_option("-p", 
                      action="store_false", dest="pp", default=True,
                      help="plot temperature")
    return result

if __name__ in ('__main__', '__plot__'):
    o, arguments  = new_option_parser().parse_args()
    main(**o.__dict__)


