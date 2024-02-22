import numpy as np
import matplotlib.pyplot as plt
import scipy.constants as physcons
from scipy.special import ive

from amuse.units import units
from amuse.community.vader.interface import Vader


kB = physcons.k*1e7 | units.erg/units.K
mH = 1e3*physcons.physical_constants['atomic mass constant'][0]*1.00794 | units.g
mu = 0.61
sigma = physcons.sigma*1e3 | units.erg/units.s / units.cm**2 / units.K**4
c = physcons.c*1e2 | units.cm/units.s
a = 4*sigma/c
dyne = units.g*units.cm/units.s**2


def setup_vader (ring_mass, init_teff, col_ratio, ring_loc, kinematic_visc):

    viscous = Vader(mode='ringrad', redirection='none')

    viscous.parameters.alpha_function = True
    viscous.parameters.inner_boundary_function = True
    viscous.parameters.inner_pressure_boundary_type = 3
    viscous.parameters.outer_boundary_function = True
    viscous.parameters.outer_pressure_boundary_type = 3
    viscous.parameters.equation_of_state_function = True
    viscous.parameters.minimum_timestep = 1e-20
    viscous.parameters.maximum_tolerated_change = 1.
    viscous.parameters.use_backwards_euler = True
    viscous.parameters.verbosity = 1
    viscous.parameters.number_of_user_parameters = 5
    viscous.parameters.interpolation_order = 1

    viscous.initialize_keplerian_grid(4096, True, 1.5e10|units.cm, 1.5e12|units.cm,
        5.97e33|units.g)

    idx = np.argmax(viscous.grid.r > ring_loc)
    init_col = ring_mass / viscous.grid.area[idx]

    viscous.set_parameter(0, kinematic_visc.value_in(units.cm**2/units.s))
    viscous.set_parameter(1, ring_loc.value_in(units.cm))
    viscous.set_parameter(2, ring_mass.value_in(units.g))
    viscous.set_parameter(3, init_col.value_in(units.g/units.cm**2)/col_ratio)
    viscous.set_parameter(4, (init_teff*kB/(mu*mH)).value_in((units.cm/units.s)**2))

    return viscous


def setup_initial_conditions (viscous, ring_mass, init_teff, col_ratio, ring_loc, 
        kinematic_visc):

    gamma = 5./3.

    idx = np.argmax(viscous.grid.r > ring_loc)
    init_col = ring_mass / viscous.grid.area[idx]
    col = np.ones(len(viscous.grid.r)) * init_col/col_ratio
    col[idx] = init_col
    pres = col * init_teff*kB / (mu*mH) + 1./3.*f_z0*a*init_teff**4
    eInt = col * init_teff*kB / (mu*mH*(gamma-1.)) + f_z0*a*init_teff**4

    viscous.grid.column_density = col
    viscous.grid.pressure = pres
    viscous.grid.internal_energy = eInt


def plot_results (viscous, times, pres, ring_mass, init_teff, col_ratio, ring_loc, 
        kinematic_visc):

    colors = ['k', 'r', 'g']

    fig = plt.figure()
    ax1 = fig.add_subplot(211)
    ax2 = fig.add_subplot(212)

    x = viscous.grid.r/ring_loc

    t0 = ring_loc**2/(12.*kinematic_visc)
    #col0 = ring_mass / (np.pi*ring_loc**2)

    idx = np.argmax(x > 1.)
    init_col = ring_mass / viscous.grid.area[idx]


    fig = plt.figure()
    ax = fig.add_subplot(111)

    for i in range(len(times)):
        ax.plot(x, (pres[i]/f_z0).value_in(dyne/units.cm**2), c=colors[i],
            label='$t/t_0=${a}'.format(a=times[i]/t0))

    ax.set_yscale('log')

    ax.set_xlim(0., 2.)
    ax.set_ylim(1e3, 1e9)

    ax.set_xlabel('$r/R_0$')
    ax.set_ylabel('$P/f_{z_0}$ [dyn cm$^{-2}$]')

    ax.legend(loc='upper right')

    plt.savefig('KrumholzForbes2015_Fig8.pdf')


def run_ringrad (ring_mass, init_teff, col_ratio, ring_loc, kinematic_visc, f_z0):

    viscous = setup_vader(ring_mass, init_teff, col_ratio, ring_loc, kinematic_visc)

    setup_initial_conditions(viscous, ring_mass, init_teff, col_ratio, ring_loc, 
        kinematic_visc)

    grid_copy = viscous.grid.copy()
    ch_from_code = viscous.grid.new_channel_to(grid_copy)


    t0 = ring_loc**2/(12.*kinematic_visc)
    times = np.array([0.002, 0.032, 0.128]) * t0

    pres = np.zeros((len(times),len(grid_copy))) | units.g/units.s**2


    for i in range(len(times)):

        viscous.evolve_model( times[i] )

        ch_from_code.copy()
        pres[i] = grid_copy.pressure


    plot_results(viscous, times, pres, ring_mass, init_teff, col_ratio, ring_loc, 
        kinematic_visc)


if __name__ == '__main__':

    ring_mass = 1.99e27 | units.g
    init_teff = 1e4 | units.K
    col_ratio = 1e10
    ring_loc = 7.5e11 | units.cm
    kinematic_visc = 1.483e11 | units.cm**2/units.s
    f_z0 = 7.5e9 | units.cm

    print ("(Partial) reproduction of Figure 8 in Krumholz & Forbes 2015")

    run_ringrad(ring_mass, init_teff, col_ratio, ring_loc, kinematic_visc, f_z0)
