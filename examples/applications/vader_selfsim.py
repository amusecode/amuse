import numpy as np
import matplotlib.pyplot as plt
import scipy.constants as physcons

from amuse.units import units
from amuse.community.vader.interface import Vader


kB = physcons.k*1e7 | units.erg/units.K
mH = 1e3*physcons.physical_constants['atomic mass constant'][0]*1.00794 | units.g
mu = 2.33


def setup_vader (R0, nu0, Mdot0, init_temp):

    viscous = Vader(mode='selfsim', redirection='none')

    viscous.parameters.alpha_function = True
    viscous.parameters.inner_pressure_boundary_type = 3
    viscous.parameters.inner_boundary_function = True
    viscous.parameters.outer_pressure_boundary_type = 3
    viscous.parameters.outer_boundary_function = True
    viscous.parameters.gamma = 1.000001
    viscous.parameters.number_of_user_parameters = 4
    viscous.parameters.begin_time = R0**2/(3.*nu0)
    viscous.parameters.verbosity = 1

    viscous.initialize_keplerian_grid(512, False, 1.5e12|units.cm, 3.e14|units.cm,
        1.99e33|units.g)

    viscous.set_parameter(0, nu0.value_in(units.cm**2/units.s))
    viscous.set_parameter(1, R0.value_in(units.cm))
    viscous.set_parameter(2, Mdot0.value_in(units.g/units.s))
    viscous.set_parameter(3, (init_temp*kB/(mu*mH)).value_in((units.cm/units.s)**2))

    return viscous


def setup_initial_conditions (viscous, R0, nu0, Mdot0, init_temp):

    x = viscous.grid.r/R0
    col1 = Mdot0/(3.*np.pi*nu0)

    col = col1 * np.exp(-x)/x
    pres = col * init_temp * kB/(mu*mH)

    viscous.grid.column_density = col
    viscous.grid.pressure = pres


def plot_results (viscous, col, R0, nu0, Mdot0, init_temp):

    t0 = R0**2/(3.*nu0)
    col1 = Mdot0/(3.*np.pi*nu0)

    fig = plt.figure()
    ax1 = fig.add_subplot(211)
    ax2 = fig.add_subplot(212)

    colors = ['r', 'g', 'b']

    x = viscous.grid.r/R0
    T = 1.

    col_ana = np.exp(-x/T)/(x*T**(3./2.))

    ax1.scatter(x[::16], col[0,::16]/col1, c='k')
    ax1.plot(x, col_ana, label='$t/t_s=$1', c='k')

    for i in range(3):
        T = i+2.

        col_ana = np.exp(-x/T)/(x*T**(3./2.))

        ax1.scatter(x[::16], col[i+1,::16]/col1, c=colors[i])
        ax1.plot(x, col_ana, label='$t/t_s=$'+str(i+2), c=colors[i])

        ax2.plot(x, np.abs((col[i+1]/col1 - col_ana)/col_ana), 
            label='$t/t_s=$'+str(i+2), c=colors[i])

    ax1.set_xlim(1e-1, 2e1)
    ax1.set_ylim(1e-10, 1e2)

    ax2.set_xlim(1e-1, 2e1)
    ax2.set_ylim(1e-7, 1e0)

    ax1.set_xlabel('$r/R_0$')
    ax1.set_ylabel('$\\Sigma/\\Sigma_0$')
    ax2.set_ylabel('|Error|')

    ax1.set_yscale('log')
    ax2.set_yscale('log')

    ax1.set_xscale('log')
    ax2.set_xscale('log')

    ax1.axes.get_xaxis().set_visible(False)
    ax1.legend()
    ax2.legend()

    plt.subplots_adjust(hspace=0)

    plt.savefig('KrumholzForbes2015_Fig1.pdf')


def run_selfsim (R0, nu0, Mdot0, init_temp):

    viscous = setup_vader(R0, nu0, Mdot0, init_temp)

    setup_initial_conditions(viscous, R0, nu0, Mdot0, init_temp)

    grid_copy = viscous.grid.copy()
    ch_from_code = viscous.grid.new_channel_to(grid_copy)

    col = np.zeros((4,len(grid_copy))) | units.g/units.cm**2
    col[0] = grid_copy.column_density

    t0 = R0**2/(3.*nu0)


    for i in range(3):

        viscous.evolve_model( (i+2)*t0 )

        ch_from_code.copy()
        col[i+1] = grid_copy.column_density


    plot_results(viscous, col, R0, nu0, Mdot0, init_temp)


if __name__ == '__main__':

    R0 = 1.5e13 | units.cm
    nu0 = 2.37e13 | units.cm**2/units.s
    Mdot0 = 6.3e19 | units.g/units.s
    init_temp = 1e2 | units.K

    print ("(Partial) reproduction of Figure 1 in Krumholz & Forbes 2015")
    print ("Slight numerical differences can arise due to differences in timestepping. The interface restarts every iteration with the initial timestep; VADER writes data during the run, slowing down if necessary, but not necessarily to the initial timestep.")

    run_selfsim(R0, nu0, Mdot0, init_temp)
