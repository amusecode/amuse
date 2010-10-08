try:
    import matplotlib.pyplot as native_plot
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
    
from amuse.support.units import units
from amuse.support.data import values

auto_label = "{0}"
custom_label = "{0} {1}"

class UnitlessArgs(object):
    current_plot = None

    @classmethod
    def strip(cli, *args, **kwargs):
        if cli.current_plot is native_plot.gca():
            args = [arg.as_quantity_in(unit) if isinstance(arg, values.Quantity) else arg 
                for arg, unit in map(None, args, cli.arg_units)]
        cli.clear()
        cli.current_plot = native_plot.gca()
        for i, v in enumerate(args):
            if isinstance(v, values.Quantity):
                cli.stripped_args.append(v.value_in(v.unit))
                cli.arg_units.append(v.unit)
                cli.unitnames_of_args.append("["+str(v.unit)+"]")
            else:
                cli.stripped_args.append(v)
                cli.arg_units.append(None)
                cli.unitnames_of_args.append("")
        
    @classmethod
    def clear(cli):
        cli.stripped_args = []
        cli.arg_units = []
        cli.unitnames_of_args = []

def latex_support():
    from matplotlib import rc
    #rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
    #rc('font',**{'family':'serif','serif':['Palatino']})
    rc('text', usetex=True)

def plot(*args, **kwargs):
    UnitlessArgs.strip(*args, **kwargs)
    args = UnitlessArgs.stripped_args
    native_plot.plot(*args, **kwargs)
    native_plot.xlabel(auto_label.format(UnitlessArgs.unitnames_of_args[0]))
    native_plot.ylabel(auto_label.format(UnitlessArgs.unitnames_of_args[1]))

def semilogx(*args, **kwargs):
    UnitlessArgs.strip(*args, **kwargs)
    args = UnitlessArgs.stripped_args
    native_plot.semilogx(*args, **kwargs)
    native_plot.xlabel(auto_label.format(UnitlessArgs.unitnames_of_args[0]))
    native_plot.ylabel(auto_label.format(UnitlessArgs.unitnames_of_args[1]))

def semilogy(*args, **kwargs):
    UnitlessArgs.strip(*args, **kwargs)
    args = UnitlessArgs.stripped_args
    native_plot.semilogy(*args, **kwargs)
    native_plot.xlabel(auto_label.format(UnitlessArgs.unitnames_of_args[0]))
    native_plot.ylabel(auto_label.format(UnitlessArgs.unitnames_of_args[1]))

def loglog(*args, **kwargs):
    UnitlessArgs.strip(*args, **kwargs)
    args = UnitlessArgs.stripped_args
    native_plot.loglog(*args, **kwargs)
    native_plot.xlabel(auto_label.format(UnitlessArgs.unitnames_of_args[0]))
    native_plot.ylabel(auto_label.format(UnitlessArgs.unitnames_of_args[1]))

def scatter(x, y, s=20, c='b', marker='o', cmap=None, norm=None, vmin=None, vmax=None, alpha=1.0, linewidths=None, faceted=True, verts=None, hold=None, **kwargs):
    UnitlessArgs.strip(x,y)
    args = UnitlessArgs.stripped_args
    native_plot.scatter(args[0], args[1], s, c, marker, cmap, norm, vmin, vmax, alpha, linewidths, faceted, verts, hold, **kwargs)

def hist(x, bins=10, range=None, normed=False, weights=None, cumulative=False, bottom=None, histtype='bar', align='mid', orientation='vertical', rwidth=None, log=False, hold=None, **kwargs):
    UnitlessArgs.strip(x)
    args = UnitlessArgs.stripped_args
    native_plot.hist(args[0], bins, range, normed, weights, cumulative, bottom, histtype, align, orientation, rwidth, log, hold, **kwargs)
    UnitlessArgs.unitnames_of_args.append("")

def xlabel(s, *args, **kwargs):
    if not '[' in s:
        s = custom_label.format(s, UnitlessArgs.unitnames_of_args[0])
    native_plot.xlabel(s, *args, **kwargs)

def ylabel(s, *args, **kwargs):
    if not '[' in s:
        s = custom_label.format(s, UnitlessArgs.unitnames_of_args[1])
    native_plot.ylabel(s, *args, **kwargs)

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
    y1 = values.new_quantity(np.sin(np.arange(0,1.5,0.03)), 1e50*units.erg)
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
