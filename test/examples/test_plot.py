from amuse.lab import *
from amuse.plot import *
import numpy as np

if __name__ == '__main__':

    latex_support()

    x = np.pi/20.0 * (range(-10,10) | units.m)
    y1 = units.MSun.new_quantity(np.sin(x.number))
    y2 = units.MSun.new_quantity(x.number)
    native_plot.subplot(2,2,1)
    plot(x, y1, label='model')
    scatter(x, y2, label='data')
    xlabel('x')
    ylabel('mass [$M_\odot$]')
    native_plot.legend(loc=2)

    x = range(50) | units.Myr
    y1 = values.new_quantity(np.sin(np.arange(0,1.5,0.03)), 1e50*units.erg)
    y2 = -(1e43 | units.J) - y1
    native_plot.subplot(2,2,2)
    plot(x, y1, label='$E_\mathrm{kin}$')
    plot(x, y2, label='$E_\mathrm{pot}$')
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

    y1 = units.N.new_quantity(np.random.normal(0.0,1.0,100))
    x = units.N.new_quantity(np.arange(-3, 3, 0.1))
    y2 = np.exp(-np.arange(-3, 3, 0.1)**2)/np.sqrt(np.pi)
    native_plot.subplot(2,2,4)
    plot(x, y2, 'y--', label='model')
    hist(y1, bins=12, range=(-3,3), normed=True, label='data')
    xlabel('force')
    ylabel('pdf')
    native_plot.legend()
    native_plot.show()
