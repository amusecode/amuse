try:
    import matplotlib.pyplot as native_plot
    from mpl_toolkits.mplot3d import Axes3D
except ImportError:
    class FakePlotLibrary(object):
        def stub(self, *args, **kwargs):
            raise Exception("No plot library available")
        plot = stub
        scatter = stub
        hist = stub
        xlabel = stub
        ylabel = stub

    native_plot = FakePlotLibrary()

import numpy
try:
    from pynbody.array import SimArray
    from pynbody.snapshot import SimSnap, _new
    import pynbody.plot.sph as pynbody_sph
    HAS_PYNBODY = True
except ImportError:
    HAS_PYNBODY = False

from amuse.support.exceptions import AmuseException
from amuse.units import units, constants
from amuse.units import quantities
from amuse.support import console

auto_label = "{0}"
custom_label = "{0} {1}"

class UnitlessArgs(object):
    current_plot = None

    @classmethod
    def strip(self, *args, **kwargs):
        if self.current_plot is native_plot.gca():
            args = [arg.as_quantity_in(unit) if quantities.is_quantity(arg) else arg
                for arg, unit in map(None, args, self.arg_units)]
        self.clear()
        self.current_plot = native_plot.gca()
        for arg in args:
            if quantities.is_quantity(arg):
                arg = console.current_printing_strategy.convert_quantity(arg)

                self.stripped_args.append(arg.value_in(arg.unit))
                self.arg_units.append(arg.unit)
                self.unitnames_of_args.append("["+str(arg.unit)+"]")
            else:
                self.stripped_args.append(arg)
                self.arg_units.append(None)
                self.unitnames_of_args.append("")

        return self.stripped_args

    @classmethod
    def clear(self):
        self.stripped_args = []
        self.arg_units = []
        self.unitnames_of_args = []

    @classmethod
    def x_label(self, s=None):
        unit_name = self.unitnames_of_args[0]

        return self.label(s, unit_name)

    @classmethod
    def y_label(self, s=None):
        unit_name = self.unitnames_of_args[1]

        return self.label(s, unit_name)

    @classmethod
    def z_label(self, s=None):
        unit_name = self.unitnames_of_args[2]

        return self.label(s, unit_name)

    @classmethod
    def label(self, s, unit_name):
        if s is None:
            return auto_label.format(unit_name)
        else:
            return custom_label.format(s, unit_name)

def latex_support():
    from matplotlib import rc
    #rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
    #rc('font',**{'family':'serif','serif':['Palatino']})
    rc('text', usetex=True)

def plot(*args, **kwargs):
    args = UnitlessArgs.strip(*args, **kwargs)
    result = native_plot.plot(*args, **kwargs)
    native_plot.xlabel(UnitlessArgs.x_label())
    native_plot.ylabel(UnitlessArgs.y_label())
    return result

def plot3(*args, **kwargs):
    args = UnitlessArgs.strip(*args, **kwargs)
    fig = native_plot.figure()
    ax = fig.gca(projection='3d')
    return ax.plot(*args, **kwargs)

def semilogx(*args, **kwargs):
    args = UnitlessArgs.strip(*args, **kwargs)
    result = native_plot.semilogx(*args, **kwargs)
    native_plot.xlabel(UnitlessArgs.x_label())
    native_plot.ylabel(UnitlessArgs.y_label())
    return result

def semilogy(*args, **kwargs):
    args = UnitlessArgs.strip(*args, **kwargs)
    result = native_plot.semilogy(*args, **kwargs)
    native_plot.xlabel(UnitlessArgs.x_label())
    native_plot.ylabel(UnitlessArgs.y_label())
    return result

def loglog(*args, **kwargs):
    args = UnitlessArgs.strip(*args, **kwargs)
    result = native_plot.loglog(*args, **kwargs)
    native_plot.xlabel(UnitlessArgs.x_label())
    native_plot.ylabel(UnitlessArgs.y_label())
    return result

def scatter(x, y, **kwargs):
    args = UnitlessArgs.strip(x,y)
    result = native_plot.scatter(*UnitlessArgs.stripped_args, **kwargs)
    native_plot.xlabel(UnitlessArgs.x_label())
    native_plot.ylabel(UnitlessArgs.y_label())
    return result

def hist(x, bins=10, range=None, normed=False, weights=None, cumulative=False, bottom=None, histtype='bar', align='mid', orientation='vertical', rwidth=None, log=False, hold=None, **kwargs):
    args = UnitlessArgs.strip(x)
    result = native_plot.hist(args[0], bins, range, normed, weights, cumulative, bottom, histtype, align, orientation, rwidth, log, hold, **kwargs)
    UnitlessArgs.unitnames_of_args.append("")
    return result

def errorbar(*args, **kwargs):
    for label in ['yerr', 'xerr']:
        if label in kwargs:
            args += (kwargs.pop(label),)
        else:
            args += (None,)

    yerr, xerr = args[2:4]

    args1 = UnitlessArgs.strip(*args[:2])
    args2 = UnitlessArgs.strip(xerr, yerr)
    args = args1 + args2[::-1]

    result = native_plot.errorbar(*args, **kwargs)
    native_plot.xlabel(UnitlessArgs.x_label())
    native_plot.ylabel(UnitlessArgs.y_label())
    return result

def text(*args, **kwargs):
    args = UnitlessArgs.strip(*args, **kwargs)
    return native_plot.text(*args, **kwargs)

def xlabel(s, *args, **kwargs):
    if not '[' in s:
        s = UnitlessArgs.x_label(s)
    return native_plot.xlabel(s, *args, **kwargs)

def ylabel(s, *args, **kwargs):
    if not '[' in s:
        s = UnitlessArgs.y_label(s)
    return native_plot.ylabel(s, *args, **kwargs)

def xlim(*args, **kwargs):
    if len(UnitlessArgs.arg_units) is 0:
        raise AmuseException("Cannot call xlim function before plotting")

    x_unit = UnitlessArgs.arg_units[0]

    if x_unit is not None:
        args = [arg.value_in(x_unit) for arg in args]

        for name in ("xmin", "xmax"):
            if name in kwargs:
                kwargs[name] = kwargs[name].value_in(x_unit)

    native_plot.xlim(*args, **kwargs)

def ylim(*args, **kwargs):
    if len(UnitlessArgs.arg_units) is 0:
        raise AmuseException("Cannot call ylim function before plotting")

    y_unit = UnitlessArgs.arg_units[1]

    if y_unit is not None:
        args = [arg.value_in(y_unit) for arg in args]

        for name in ("ymin", "ymax"):
            if name in kwargs:
                kwargs[name] = kwargs[name].value_in(y_unit)

    native_plot.ylim(*args, **kwargs)

def fix_xyz_axes(X, Y, Z):
    if not (X.shape == Z.shape and Y.shape == Z.shape):
        X, Y = numpy.meshgrid(X, Y)
    return X, Y, Z

def log_norm(Z):
    # for log scale, 0 is considered a missing value
    masked_Z = numpy.ma.masked_equal(Z, 0.0, copy=False)
    zmin, zmax = masked_Z.min(), masked_Z.max()

    from matplotlib.colors import LogNorm
    return LogNorm(vmin=zmin, vmax=zmax)

def has_log_scaling(array):
    diff = numpy.diff(array)
    if numpy.all(diff - diff[0] < diff[0]/10.):
        return False

    logdiff = numpy.diff(numpy.log10(array))
    if numpy.all(logdiff - logdiff[0] < logdiff[0]/10.):
        return True

    raise AmuseException("This method cannot be used for non regular arrays")

def density_plot(x, y, z, label=None, add_colorbar=True, zlog=False, **kwargs):
    """
        Plot a density matrix as a color map, similar to pyplot.contour.
        I don't understand why pyplot does not have this out of the box.
        It only works if x and y are regular (linear or logarithmic).
    """
    X, Y, Z = UnitlessArgs.strip(x, y, z)
    X, Y, Z = fix_xyz_axes(X, Y, Z)

    xlow = X[0,0]
    xhigh = X[-1,-1]
    ylow = Y[0,0]
    yhigh = Y[-1,-1]
    extent = (xlow, xhigh, ylow, yhigh)

    norm = log_norm(Z) if zlog else None

    cax = native_plot.imshow(Z, origin='lower', aspect='auto', extent=extent, norm=norm, **kwargs)

    if has_log_scaling(X[0,:]):
        native_plot.gca().set_xscale('log')
    if has_log_scaling(Y[:,0]):
        native_plot.gca().set_yscale('log')

    native_plot.xlabel(UnitlessArgs.x_label())
    native_plot.ylabel(UnitlessArgs.y_label())

    if add_colorbar:
        bar = native_plot.colorbar(cax)
        bar.set_label(UnitlessArgs.z_label(label))

        return cax, bar
    else:
        return cax

def contour(*args, **kwargs):
    if len(args)%2 == 0:
        stripped_args = UnitlessArgs.strip(*args[:-1])

        levels = args[-1]
        z_unit = UnitlessArgs.arg_units[-1]

        if quantities.is_quantity(levels):
            stripped_args.append(levels.value_in(z_unit))
    else:
        stripped_args = UnitlessArgs.strip(*args)


    if 'levels' in kwargs:
        levels = kwargs['levels']
        z_unit = UnitlessArgs.arg_units[-1]

        if quantities.is_quantity(levels):
            kwargs['levels'] = levels.value_in(z_unit)

    result = native_plot.contour(*stripped_args, **kwargs)

    native_plot.xlabel(UnitlessArgs.x_label())
    native_plot.ylabel(UnitlessArgs.y_label())

    return result

def smart_length_units_for_vector_quantity(quantity):
    length_units = [units.Mpc, units.kpc, units.parsec, units.AU, units.RSun, units.km]
    total_size = max(quantity) - min(quantity)
    for length_unit in length_units:
        if total_size > (1 | length_unit):
            return length_unit
    return units.m

def sph_particles_plot(particles, u_range = None, min_size = 100, max_size = 10000,
        alpha = 0.1, gd_particles=None, width=None, view=None):
    """
    Very simple and fast procedure to make a plot of the hydrodynamics state of
    a set of SPH particles. The particles must have the following attributes defined:
    position, u, h_smooth.
    For a more accurate plotting procedure, see for example:
    examples/applications/christmas_card_2010.py

    :argument particles: the SPH particles to be plotted
    :argument u_range: range of internal energy for color scale [umin, umax]
    :argument min_size: minimum size to use for plotting particles, in pixel**2
    :argument max_size: maximum size to use for plotting particles, in pixel**2
    :argument alpha: the opacity of each particle
    :argument gd_particles: non-SPH particles can be indicated with white circles
    :argument view: the (physical) region to plot [xmin, xmax, ymin, ymax]
    """
    positions = particles.position
    us        = particles.u
    h_smooths = particles.h_smooth
    x, y, z = positions.x, positions.y, positions.z
    z, x, y, us, h_smooths = z.sorted_with(x, y, us, h_smooths)

    if u_range:
        u_min, u_max = u_range
    else:
        u_min, u_max = min(us), max(us)
    log_u = numpy.log((us / u_min)) / numpy.log((u_max / u_min))
    clipped_log_u = numpy.minimum(numpy.ones_like(log_u), numpy.maximum(numpy.zeros_like(log_u), log_u))

    red   = 1.0 - clipped_log_u**4
    blue  = clipped_log_u**4
    green = numpy.minimum(red, blue)

    colors = numpy.transpose(numpy.array([red, green, blue]))
    n_pixels = native_plot.gcf().get_dpi() * native_plot.gcf().get_size_inches()

    current_axes = native_plot.gca()
    current_axes.set_axis_bgcolor('#101010')
    if width is not None:
        view = width * [-0.5, 0.5, -0.5, 0.5]

    if view:
        current_axes.set_aspect("equal", adjustable = "box")
        length_unit = smart_length_units_for_vector_quantity(view)
        current_axes.set_xlim(view[0].value_in(length_unit),
            view[1].value_in(length_unit), emit=True, auto=False)
        current_axes.set_ylim(view[2].value_in(length_unit),
            view[3].value_in(length_unit), emit=True, auto=False)
        phys_to_pix2 = n_pixels[0]*n_pixels[1] / ((view[1]-view[0])**2 + (view[3]-view[2])**2)
    else:
        current_axes.set_aspect("equal", adjustable = "datalim")
        length_unit = smart_length_units_for_vector_quantity(x)
        phys_to_pix2 = n_pixels[0]*n_pixels[1] / ((max(x)-min(x))**2 + (max(y)-min(y))**2)
    sizes = numpy.minimum(numpy.maximum((h_smooths**2 * phys_to_pix2), min_size), max_size)

    x = x.as_quantity_in(length_unit)
    y = y.as_quantity_in(length_unit)
    scatter(x, y, s=sizes, c=colors, edgecolors="none", alpha=alpha)
    if gd_particles:
        scatter(gd_particles.x, gd_particles.y, c='w', marker='o')
    xlabel('x')
    ylabel('y')

def convert_particles_to_pynbody_data(particles, length_unit=units.kpc, pynbody_unit="kpc"):
    if not HAS_PYNBODY:
        raise AmuseException("Couldn't find pynbody")

    if hasattr(particles, "u"):
        pynbody_data = _new(gas=len(particles))
    else:
        pynbody_data = _new(dm=len(particles))
    pynbody_data._filename = "AMUSE"
    if hasattr(particles, "mass"):
        pynbody_data['mass'] = SimArray(particles.mass.value_in(units.MSun), "Msol")
    if hasattr(particles, "position"):
        pynbody_data['x'] = SimArray(particles.x.value_in(length_unit), pynbody_unit)
        pynbody_data['y'] = SimArray(particles.y.value_in(length_unit), pynbody_unit)
        pynbody_data['z'] = SimArray(particles.z.value_in(length_unit), pynbody_unit)
    if hasattr(particles, "velocity"):
        pynbody_data['vx'] = SimArray(particles.vx.value_in(units.km / units.s), "km s^-1")
        pynbody_data['vy'] = SimArray(particles.vy.value_in(units.km / units.s), "km s^-1")
        pynbody_data['vz'] = SimArray(particles.vz.value_in(units.km / units.s), "km s^-1")
    if hasattr(particles, "h_smooth"):
        pynbody_data['smooth'] = SimArray(particles.h_smooth.value_in(length_unit), pynbody_unit)
    if hasattr(particles, "rho"):
        pynbody_data['rho'] = SimArray(particles.rho.value_in(units.g / units.cm**3),
            "g cm^-3")
    if hasattr(particles, "u"):
#        pynbody_data['u'] = SimArray(particles.u.value_in(units.km**2 / units.s**2), "km^2 s^-2")
        temp = 2.0/3.0 * particles.u * mu() / constants.kB
        pynbody_data['temp'] = SimArray(temp.value_in(units.K), "K")
    return pynbody_data

def mu(X = None, Y = 0.25, Z = 0.02, x_ion = 0.1):
    """
    Compute the mean molecular weight in kg (the average weight of particles in a gas)
    X, Y, and Z are the mass fractions of Hydrogen, of Helium, and of metals, respectively.
    x_ion is the ionisation fraction (0 < x_ion < 1), 1 means fully ionised
    """
    if X is None:
        X = 1.0 - Y - Z
    elif abs(X + Y + Z - 1.0) > 1e-6:
        raise AmuseException("Error in calculating mu: mass fractions do not sum to 1.0")
    return constants.proton_mass / (X*(1.0+x_ion) + Y*(1.0+2.0*x_ion)/4.0 + Z*x_ion/2.0)

def _smart_length_units_for_pynbody_data(length):
    length_units = [(units.Gpc, "Gpc"), (units.Mpc, "Mpc"), (units.kpc, "kpc"),
        (units.parsec, "pc"), (units.AU, "au"), (1.0e9*units.m, "1.0e9 m"),
        (1000*units.km, "1000 km"), (units.km, "km")]
    for length_unit, pynbody_unit in length_units:
        if length > (1 | length_unit):
            return length_unit, pynbody_unit
    return units.m, "m"

def pynbody_column_density_plot(particles, width=None, qty='rho', units=None,
        sideon=False, faceon=False, **kwargs):
    if not HAS_PYNBODY:
        raise AmuseException("Couldn't find pynbody")

    if width is None:
        width = 2.0 * particles.position.lengths_squared().amax().sqrt()
    length_unit, pynbody_unit = _smart_length_units_for_pynbody_data(width)
    pyndata = convert_particles_to_pynbody_data(particles, length_unit, pynbody_unit)
    UnitlessArgs.strip([1]|length_unit, [1]|length_unit)

    if sideon:
        function = pynbody_sph.sideon_image
    elif faceon:
        function = pynbody_sph.faceon_image
    else:
        function = pynbody_sph.image

    if units is None and qty == 'rho':
        units = 'm_p cm^-2'

    result = function(pyndata, width=width.value_in(length_unit), qty=qty, units=units, **kwargs)
    UnitlessArgs.current_plot = native_plot.gca()
    return result

def effective_iso_potential_plot(gravity_code,
        omega,
        center_of_rotation = [0, 0]|units.AU,
        xlim = [-1.5, 1.5] | units.AU,
        ylim = [-1.5, 1.5] | units.AU,
        resolution = [1000, 1000],
        number_of_contours = 20,
        fraction_screen_filled = 0.5,
        quadratic_contour_levels = True,
        contour_kwargs = dict(),
        omega2 = None,
        center_of_rotation2 = [0, 0]|units.AU,
        fraction_screen_filled2 = 0.2):
    """
    Create a contour plot of the effective potential of particles in a gravity code.
    The code needs to support 'get_potential_at_point' only, so it can also be an
    instance of Bridge.

    :argument gravity_code: an instance of a gravity code
    :argument omega: The angular velocity of the system
    :argument center_of_rotation: The (2D) center around which the system rotates, usually the center of mass
    :argument xlim: Range in x coordinate; width of window
    :argument ylim: Range in y coordinate; width of window
    :argument resolution: Number of points to sample potential for x and y direction
    :argument number_of_contours: How many contour lines to plot
    :argument fraction_screen_filled: Lowest contour will enclose this fraction of the screen
    :argument quadratic_contour_levels: Quadratic or linear scaling between contour levels
    :argument contour_kwargs: Optional keyword arguments for pyplot.contour
    """
    UnitlessArgs.strip(xlim, ylim)
    xlim, ylim = UnitlessArgs.stripped_args

    x_num = numpy.linspace(xlim[0], xlim[1], resolution[0])
    y_num = numpy.linspace(ylim[0], ylim[1], resolution[1])
    x_num, y_num = numpy.meshgrid(x_num, y_num)
    x = (x_num | UnitlessArgs.arg_units[0]).flatten()
    y = (y_num | UnitlessArgs.arg_units[1]).flatten()
    zeros = x.aszeros()
    potential = gravity_code.get_potential_at_point(zeros, x, y, zeros)
    potential -= omega**2 * ((x-center_of_rotation[0])**2 + (y-center_of_rotation[1])**2) / 2.0

    levels = set_contour_levels(potential, number_of_contours, fraction_screen_filled, quadratic_contour_levels)
    CS = native_plot.contour(x_num, y_num, potential.number.reshape(resolution[::-1]), levels, **contour_kwargs)
    #~native_plot.clabel(CS, inline=1, fontsize=10)

    if omega2 is None:
        return potential

    potential2 = potential - omega2**2 * ((x-center_of_rotation2[0])**2 + (y-center_of_rotation2[1])**2) / 2.0
    #~levels = set_contour_levels(potential, number_of_contours2, fraction_screen_filled2, quadratic_contour_levels2)
    levels = set_contour_levels(potential2, number_of_contours, fraction_screen_filled2, quadratic_contour_levels)
    CS = native_plot.contour(x_num, y_num, potential2.number.reshape(resolution[::-1]), levels, **contour_kwargs)
    return potential.reshape(resolution[::-1]), potential2.reshape(resolution[::-1])

def set_contour_levels(potential, number_of_contours, fraction_screen_filled, quadratic_contour_levels):
    uniform = numpy.linspace(0.0, 1.0, number_of_contours)
    V_max = potential.amax().number
    V_min = potential.sorted().number[int(len(potential)*(1-fraction_screen_filled))]
    if quadratic_contour_levels:
        levels = V_min + (V_max-V_min) * uniform * (2 - uniform)
    else:
        levels = V_min + (V_max-V_min) * uniform
    return levels


if __name__ == '__main__':
    import numpy as np

    latex_support()

    x = np.pi/20.0 * (range(-10,10) | units.m)
    y1 = units.m.new_quantity(np.sin(x.number))
    y2 = x
    native_plot.subplot(2,2,1)
    plot(x, y1, label='model')
    scatter(x, y2, label='data')
    xlabel('x')
    ylabel('$M_\odot y$')
    native_plot.legend(loc=2)

    x = range(50) | units.Myr
    y1 = quantities.new_quantity(np.sin(np.arange(0,1.5,0.03)), 1e50*units.erg)
    y2 = -(1e43 | units.J) - y1
    native_plot.subplot(2,2,2)
    plot(x, y1, label='$E_{kin}$')
    plot(x, y2, label='$E_{pot}$')
    xlabel('t')
    ylabel('E')
    native_plot.legend()

    x = range(7) | units.day
    y1 = [0, 4, 2, 3, 2, 5, 1]
    y2 = [3, 0, 2, 2, 3, 0, 4]
    native_plot.subplot(2,2,3)
    plot(x, y1, 'ks', label='coffee')
    plot(x, y2, 'yo', label='tea')
    xlabel('time')
    ylabel('consumption / day')
    native_plot.legend()

    y1 = units.N.new_quantity(np.random.normal(0.,1.,100))
    x = (units.g * units.cm * units.s**-2).new_quantity(np.arange(-3.0e5, 3.0e5, 0.1e5))
    y2 = np.exp(-np.arange(-3., 3., 0.1)**2)/np.sqrt(np.pi)
    native_plot.subplot(2,2,4)
    hist(y1, bins=12, range=(-3,3), normed=True, label='data')
    plot(x, y2, 'y--', label='model')
    xlabel('force')
    ylabel('pdf')
    native_plot.legend()
    native_plot.show()
