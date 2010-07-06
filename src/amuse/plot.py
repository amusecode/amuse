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
    
import numpy as np
from amuse.support.units import units
from amuse.support.data import values

auto_label = "[{0}]"
custom_label = "{0} [{1}]"

class UnitlessArgs(object):
    stripped_args = []
    unitnames_of_args = []

    @classmethod
    def strip(cli, *args, **kwargs):
        cli.clear()
        for i, v in enumerate(args):
            if isinstance(v, values.Quantity):
                stripped = v.value_in(v.unit)
                cli.stripped_args.append(stripped)
                cli.unitnames_of_args.append(v.unit.name)
            else:
                cli.stripped_args.append(v)
        
    @classmethod
    def clear(cli):
        stripped_args = []
        unitnames_of_args = []

def plot(*args, **kwargs):
    UnitlessArgs.strip(*args, **kwargs)
    args = UnitlessArgs.stripped_args
    native_plot.plot(*args, **kwargs)
    native_plot.xlabel(auto_label.format(UnitlessArgs.unitnames_of_args[0]))
    native_plot.ylabel(auto_label.format(UnitlessArgs.unitnames_of_args[1]))

def scatter(x, y, s=20, c='b', marker='o', cmap=None, norm=None, vmin=None, vmax=None, alpha=1.0, linewidths=None, faceted=True, verts=None, hold=None, **kwargs):
    UnitlessArgs.strip(x,y)
    args = UnitlessArgs.stripped_args
    native_plot.scatter(args[0], args[1], s, c, marker, cmap, norm, vmin, vmax, alpha, linewidths, faceted, verts, hold, **kwargs)

def hist(x, bins=10, range=None, normed=False, weights=None, cumulative=False, bottom=None, histtype='bar', align='mid', orientation='vertical', rwidth=None, log=False, hold=None, **kwargs):
    UnitlessArgs.strip(x)
    args = UnitlessArgs.stripped_args
    native_plot.scatter(args[0], bins, range, normed, weights, cumulative, bottom, histtype, align, orientation, rwidth, log, hold, **kwargs)

def xlabel(s, *args, **kwargs):
    s = custom_label.format(s, UnitlessArgs.unitnames_of_args[0])
    native_plot.xlabel(s, *args, **kwargs)

def ylabel(s, *args, **kwargs):
    s = custom_label.format(s, UnitlessArgs.unitnames_of_args[1])
    native_plot.ylabel(s, *args, **kwargs)

if __name__ == '__main__':
    x = range(-10,10)|units.m
    y = [i.number**2 for i in x]|units.m
    native_plot.subplot(2,1,1)
    plot(x,y)
    xlabel('x')
    ylabel('y')
    native_plot.subplot(2,1,2)
    scatter(x,y)
    xlabel('x')
    ylabel('y')
    native_plot.show()
