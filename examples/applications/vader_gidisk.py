import numpy as np
import matplotlib.pyplot as plt
import scipy.constants as physcons

from amuse.units import units
from amuse.community.vader.interface import Vader


G = physcons.G*1e3 | units.cm**3 / units.g / units.s**2


def setup_vader (params):

    torb = 2.*np.pi*params['R_out']/params['vphi']
    chi = G*params['mdot']/params['vphi']**3
    s = 1./2.**0.5 * (chi/params['eta'])**(1./3.)

    h_steady = params['gamma']/(params['gamma']-1.) * (s*params['vphi'])**2
    torb = 2.*np.pi*params['R_out']/params['vphi']

    viscous = Vader(mode='gidisk', redirection='none')

    viscous.parameters.alpha_function = True
    viscous.parameters.inner_pressure_boundary_type = 3
    viscous.parameters.inner_boundary_function = True
    viscous.parameters.outer_pressure_boundary_type = 1
    viscous.parameters.outer_pressure_boundary_mass_flux = -params['mdot']
    viscous.parameters.outer_enthalpy_boundary_enthalpy = h_steady * \
        params['obc_vdisp']**2 / params['init_col']
    viscous.parameters.gamma = params['gamma']
    viscous.parameters.internal_energy_source_function = True
    viscous.parameters.number_of_user_parameters = 3
    viscous.parameters.verbosity = 1
    viscous.parameters.initial_timestep = params['dt_init'] * torb

    viscous.initialize_flat_grid(512, False, 3.09e20|units.cm, params['R_out'],
        params['vphi'])

    viscous.set_parameter(0, params['eta'])
    viscous.set_parameter(1, chi)
    viscous.set_parameter(2, params['t_Q'])

    return viscous


def setup_initial_conditions (viscous, params):

    chi = G*params['mdot']/params['vphi']**3
    s = 1./2.**0.5 * (chi/params['eta'])**(1./3.)

    col1 = params['vphi']**2 * (chi/params['eta'])**(1./3.) / \
        (np.pi*G*params['R_out'])
    colSteady = col1 * (params['R_out']/viscous.grid.r)
    presSteady = colSteady * (s*params['vphi'])**2

    col = colSteady * params['init_col']
    pres = presSteady * params['init_col'] * params['init_vdisp']**2

    viscous.grid.column_density = col
    viscous.grid.pressure = pres


def run_gidisk (params):

    viscous = setup_vader(params)

    setup_initial_conditions(viscous, params)

    grid_copy = viscous.grid.copy()
    ch_from_code = viscous.grid.new_channel_to(grid_copy)


    fig = plt.figure()
    ax1 = fig.add_subplot(311)
    ax2 = fig.add_subplot(312)
    ax3 = fig.add_subplot(313)

    ax1.plot(grid_copy.r/params['R_out'], params['R_out']/grid_copy.r, c='k',
        linestyle='--', label='Steady state', lw=4)
    ax2.plot(grid_copy.r/params['R_out'], 
        np.ones(len(grid_copy.r))*2.**-0.5 * \
        (G*params['mdot']/params['vphi']**3/params['eta'])**(1./3.), c='k',
        linestyle='--', lw=4)
    ax3.plot(grid_copy.r/params['R_out'], np.ones(len(grid_copy.r)), c='k',
        linestyle='--', lw=4)

    ax1.plot(grid_copy.r/params['R_out'], 
        grid_copy.column_density/grid_copy.column_density[-1], c='b', 
        label='Simulation, $T=${a}'.format(a=0.))
    ax2.plot(grid_copy.r/params['R_out'], 
        (grid_copy.pressure/grid_copy.column_density)**0.5/params['vphi'], c='b')
    Q = 2.**0.5 * grid_copy.rotational_velocity/grid_copy.r * \
        (grid_copy.pressure/grid_copy.column_density)**0.5 / \
        (np.pi*G*grid_copy.column_density)
    ax3.plot(grid_copy.r/params['R_out'], Q, c='b')


    torb = 2.*np.pi*params['R_out']/params['vphi']

    times = np.array([0.001, 0.1, 1.]) * torb
    colors = ['g', 'r', 'c']

    for i in range(len(times)):
        viscous.evolve_model( times[i] )
        ch_from_code.copy()


        ax1.plot(grid_copy.r/params['R_out'], 
            grid_copy.column_density/grid_copy.column_density[-1], c=colors[i],
            label='Simulation, $T=${a}'.format(a=times[i]/torb))
        ax2.plot(grid_copy.r/params['R_out'], 
            (grid_copy.pressure/grid_copy.column_density)**0.5/params['vphi'], 
            c=colors[i])
        Q = 2.**0.5 * grid_copy.rotational_velocity/grid_copy.r * \
            (grid_copy.pressure/grid_copy.column_density)**0.5 / \
            (np.pi*G*grid_copy.column_density)
        ax3.plot(grid_copy.r/params['R_out'], Q, c=colors[i])


    ax1.set_xscale('log')
    ax1.set_yscale('log')

    ax2.set_xscale('log')
    ax2.set_yscale('log')

    ax3.set_xscale('log')

    ax1.set_xlim(1e-2, 1e0)
    ax1.set_ylim(1e0, 3e3)

    ax2.set_xlim(1e-2, 1e0)
    ax2.set_ylim(1e-2, 2e-1)

    ax3.set_xlim(1e-2, 1e0)
    ax3.set_ylim(0., 2.5)

    ax3.set_xlabel('$r/R$')
    ax1.set_ylabel('$\\Sigma/\\Sigma(R)$')
    ax2.set_ylabel('$\\sigma/v_\\phi$')
    ax3.set_ylabel('$Q$')

    ax1.axes.get_xaxis().set_visible(False)
    ax2.axes.get_xaxis().set_visible(False)

    ax1.legend(loc='upper right', ncol=2, frameon=False)

    plt.subplots_adjust(hspace=0)

    plt.savefig('KrumholzForbes2015_Fig6.pdf')


if __name__ == '__main__':

    params = {
        'eta': 1.5,
        'n_orbit': 4.,
        't_Q': 1.,
        'init_col': 1.,
        'init_vdisp': 0.5,
        'obc_vdisp': 0.5,
        'dt_init': 1e-5,
        'gamma': 1.666666666667,

        'vphi': 2.2e7 | units.cm/units.s,
        'mdot': 6.3e25 | units.g/units.s,
        'R_out': 3.09e22 | units.cm,
    }

    print ("Reproduction of Figure 6 in Krumholz & Forbes 2015")

    run_gidisk(params)
