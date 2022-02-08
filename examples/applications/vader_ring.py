import numpy as np
import matplotlib.pyplot as plt
import scipy.constants as physcons
from scipy.special import ive

from amuse.units import units
from amuse.community.vader.interface import Vader


kB = physcons.k*1e7 | units.erg/units.K
mH = 1e3*physcons.physical_constants['atomic mass constant'][0]*1.00794 | units.g
mu = 2.33


def setup_vader (ring_mass, init_temp, col_ratio, ring_loc, kinematic_visc):

    viscous = Vader(mode='ring', redirection='none')

    viscous.parameters.alpha_function = True
    viscous.parameters.inner_pressure_boundary_type = 3
    viscous.parameters.inner_boundary_function = True
    viscous.parameters.outer_pressure_boundary_type = 3
    viscous.parameters.outer_boundary_function = True
    viscous.parameters.gamma = 1.000001
    viscous.parameters.minimum_timestep = 1e-20
    viscous.parameters.number_of_user_parameters = 5
    viscous.parameters.verbosity = 1

    viscous.initialize_keplerian_grid(4096, True, 1.5e12|units.cm, 1.5e14|units.cm,
        1.99e33|units.g)

    idx = np.argmax(viscous.grid.r > ring_loc)
    init_col = ring_mass / viscous.grid.area[idx]

    viscous.set_parameter(0, kinematic_visc.value_in(units.cm**2/units.s))
    viscous.set_parameter(1, ring_loc.value_in(units.cm))
    viscous.set_parameter(2, ring_mass.value_in(units.g))
    viscous.set_parameter(3, init_col.value_in(units.g/units.cm**2)/col_ratio)
    viscous.set_parameter(4, (init_temp*kB/(mu*mH)).value_in((units.cm/units.s)**2))

    return viscous


def setup_initial_conditions (viscous, ring_mass, init_temp, col_ratio, ring_loc, 
        kinematic_visc):

    t0 = ring_loc**2/(12.*kinematic_visc)
    col0 = ring_mass / (np.pi*ring_loc**2)
    pres0 = col0 * init_temp*kB/(mu*mH)

    idx = np.argmax(viscous.grid.r > ring_loc)
    init_col = ring_mass / viscous.grid.area[idx]
    col = np.ones(len(viscous.grid.r))*init_col/col_ratio
    col[idx] = init_col
    pres = col*init_temp*kB/(mu*mH)

    viscous.grid.column_density = col
    viscous.grid.pressure = pres


def plot_results (viscous, times, col, ring_mass, init_temp, col_ratio, ring_loc, 
        kinematic_visc):

    fig = plt.figure()
    ax1 = fig.add_subplot(211)
    ax2 = fig.add_subplot(212)

    colors = ['k', 'r', 'g', 'b']

    x = viscous.grid.r/ring_loc

    t0 = ring_loc**2/(12.*kinematic_visc)
    col0 = ring_mass / (np.pi*ring_loc**2)

    idx = np.argmax(x > 1.)
    init_col = ring_mass / viscous.grid.area[idx]

    for i in range(len(times)):

        ax1.scatter(x[::64], (col[i]/col0)[::64], c=colors[i])

        col_analytic = (col0 / (times[i]/t0) * x**-0.25 * np.exp(-(1-x)**2 / \
            (times[i]/t0)) * ive(0.25, 2*x/(times[i]/t0)))

        ax1.plot(x, col_analytic/col0, c=colors[i], 
            label='$t/t_0=$'+str(times[i]/t0))

        ax2.plot(x, np.abs((col[i] - col_analytic - init_col/col_ratio) / \
            (col_analytic + init_col/col_ratio**0.5)), c=colors[i])

    ax1.set_yscale('log')
    ax2.set_yscale('log')

    ax1.set_xlim(0., 2.)
    ax1.set_ylim(1e-6, 1e7)

    ax2.set_xlim(0., 2.)
    ax2.set_ylim(1e-6, 1e-1)

    ax2.set_xlabel('r/R$_0$')
    ax1.set_ylabel('$\\Sigma/\\Sigma_0$')
    ax2.set_ylabel('Error')

    ax1.axes.get_xaxis().set_visible(False)
    ax1.legend()

    plt.subplots_adjust(hspace=0)

    plt.savefig('KrumholzForbes2015_Fig4.pdf')


def run_ring (ring_mass, init_temp, col_ratio, ring_loc, kinematic_visc):

    viscous = setup_vader(ring_mass, init_temp, col_ratio, ring_loc, kinematic_visc)

    setup_initial_conditions(viscous, ring_mass, init_temp, col_ratio, ring_loc, 
        kinematic_visc)

    grid_copy = viscous.grid.copy()
    ch_from_code = viscous.grid.new_channel_to(grid_copy)


    t0 = ring_loc**2/(12.*kinematic_visc)
    times = np.array([0.004, 0.008, 0.032, 0.128]) * t0

    col = np.zeros((len(times),len(grid_copy))) | units.g/units.cm**2


    for i in range(len(times)):

        viscous.evolve_model( times[i] )

        ch_from_code.copy()
        col[i] = grid_copy.column_density


    plot_results(viscous, times, col, ring_mass, init_temp, col_ratio, ring_loc, 
        kinematic_visc)


if __name__ == '__main__':

    ring_mass = 1.99e27 | units.g
    init_temp = 1e2 | units.K
    col_ratio = 1e10
    ring_loc = 7.5e13 | units.cm
    kinematic_visc = 5.93e12 | units.cm**2/units.s

    print ("(Partial) reproduction of Figure 4 in Krumholz & Forbes 2015")

    run_ring(ring_mass, init_temp, col_ratio, ring_loc, kinematic_visc)
