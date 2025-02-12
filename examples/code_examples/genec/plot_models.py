import sys
import time
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from amuse.units import units


class StellarModelPlot:
    def __init__(self, star):
        self.default_units = {
            'temperature': units.K,
            'luminosity': units.LSun,
            'radius': units.RSun,
            'mass': units.MSun,
            'age': units.Myr,
            'density': units.g * units.cm**-3,
        }
        self.star = star
        self.__n_zones = 0
        self.__redraw = True
        # self.__phase = 1
        self.__age = [] | self.default_units['age']
        self.__mass = [] | self.default_units['mass']
        self.__radius = [] | self.default_units['radius']
        self.__luminosity = [] | self.default_units['luminosity']
        self.__temperature = [] | self.default_units['temperature']
        self.__central_density = [] | self.default_units['density']
        self.__central_temperature = [] | self.default_units['temperature']
        self.__central_abundance = {}
        self.__species = np.array(star.get_names_of_species())
        self.__mainspecies = ["H", "He", "C12", "N14", "O16", "Ne20"]
        # , "Si", "Ni", "Fe"]
        # self.__zones = star.get_number_of_zones()
        self.__abundances = star.get_chemical_abundance_profiles()
        for species in self.__species:
            self.__central_abundance[species.capitalize()] = []

        self.__figures = {
            "1": plt.figure(num="Evolution", figsize=(8, 6)),
            "2": plt.figure(num="Structure", figsize=(8, 6)),
        }
        for figure in self.__figures.values():
            figure.set_tight_layout(True)

        self.__initialised = []

        ax_temp = self.__figures["2"].add_subplot(2, 2, 4)
        ax_pres = ax_temp.twinx()
        ax_lum = self.__figures["2"].add_subplot(2, 2, 3)
        ax_energy = ax_lum.twinx()
        self.__axes = {
            "HR": self.__figures["1"].add_subplot(2, 2, 1),
            "temp_dens": self.__figures["1"].add_subplot(2, 2, 2),
            "central abundance": self.__figures["1"].add_subplot(2, 2, 3),
            "Kippenhahn diagram": self.__figures["1"].add_subplot(2, 2, 4),
            "Abundance profile": self.__figures["2"].add_subplot(2, 2, 1),
            "gradient": self.__figures["2"].add_subplot(2, 2, 2),
            "Luminosity and energy production": (
                ax_lum,
                ax_energy,
            ),
            "Temperature and pressure": (
                ax_temp,
                ax_pres,
            ),
        }
        # self.__ax_hr = fig.add_subplot(2, 2, 1)
        self.update(star)
        self.plot_all()
        plt.pause(0.01)

    def update(self, star):
        self.star = star
        if self.__n_zones != star.n_zones:
            self.__redraw = True
            self.__n_zones = star.n_zones
        self.__age.append(star.age)
        self.__mass.append(star.mass)
        self.__radius.append(star.radius)
        self.__luminosity.append(star.luminosity)
        self.__temperature.append(star.temperature)
        self.__central_density.append(star.central_density)
        self.__central_temperature.append(star.central_temperature)
        # self.__zones = star.get_number_of_zones()
        self.__abundances = star.get_chemical_abundance_profiles()
        for i, species in enumerate(self.__species):
            self.__central_abundance[species.capitalize()].append(
                self.__abundances[i, 0]
            )

    def draw(
        self, verbose=False, speed=0 | units.Myr / units.minute, step=None
    ):
        for i, fig in self.__figures.items():
            if verbose:
                print(f"drawing figure {i}")

            fig.suptitle(
                f"age: {self.star.age.in_(units.Myr)}, "
                f"mass: {self.star.mass.in_(units.MSun)}, "
                f"radius: {self.star.radius.in_(units.RSun)}, "
                f"step: {step}, "
                # f"log(speedup) = {np.log10(speed):.2f}"
                f"speed: {speed.in_(units.Myr / units.minute)}"
            )
            fig.canvas.draw_idle()
            fig.canvas.flush_events()

    def initialise_hr(self):
        title = 'HR'
        unit_temp = self.default_units['temperature']
        unit_lum = self.default_units['luminosity']
        ax = self.__axes[title]
        ax.set_title(title)
        log_tstar = np.log10(
            self.star.temperature.value_in(unit_temp)
        )
        ax.set_xlim(
            log_tstar + 0.1,
            log_tstar - 0.1,
        )
        log_lstar = np.log10(
            self.star.luminosity.value_in(unit_lum)
        )
        ax.set_ylim(log_lstar - 0.1, log_lstar + 0.1)
        ax.set_xlabel(f'log (temperature / {unit_temp})')
        ax.set_ylabel(f'log (luminosity / {unit_lum})')
        ax.ticklabel_format(style='sci', scilimits=(-3, 4), axis='both')
        self.__hr_plot, = ax.plot(
            [np.log10(self.star.temperature.value_in(unit_temp))],
            [np.log10(self.star.luminosity.value_in(unit_lum))],
            # color='red',
            # linewidth=2,
        )
        self.__hr_star = ax.scatter(
            [np.log10(self.star.temperature.value_in(unit_temp))],
            [np.log10(self.star.luminosity.value_in(unit_lum))],
            facecolor='red',
            s=40,
            edgecolor=None,
            alpha=0.5,
        )
        self.__initialised.append(title)

    def plot_hr(self):
        title = 'HR'
        if title not in self.__initialised:
            self.initialise_hr()
        unit_temp = self.default_units['temperature']
        unit_lum = self.default_units['luminosity']
        ax = self.__axes[title]
        log_temperature = np.log10(
            self.__temperature.value_in(unit_temp)
        )
        log_luminosity = np.log10(
            self.__luminosity.value_in(unit_lum)
        )
        ax.set_xlim(
            max(log_temperature) + 0.1,
            min(log_temperature) - 0.1,
        )
        ax.set_ylim(
            min(log_luminosity) - 0.1,
            max(log_luminosity) + 0.1,
        )
        self.__hr_plot.set_xdata(log_temperature)
        self.__hr_plot.set_ydata(log_luminosity)
        offsets = self.__hr_star.get_offsets()
        offsets[:, 0] = [np.log10(self.star.temperature.value_in(unit_temp))]
        offsets[:, 1] = [np.log10(self.star.luminosity.value_in(unit_lum))]
        self.__hr_star.set_offsets(offsets)
        # fig = self.__figures
        # fig.canvas.draw_idle()

    def initialise_temp_dens(self):
        title = 'temp_dens'
        unit_temp = self.default_units['temperature']
        unit_dens = self.default_units['density']
        ax = self.__axes[title]
        ax.set_title(title)
        log_central_density_star = np.log10(
            self.star.central_density.value_in(unit_dens)
        )
        ax.set_xlim(
            log_central_density_star - 0.1,
            log_central_density_star + 0.1
        )
        log_central_temperature_star = np.log10(
            self.star.central_temperature.value_in(unit_temp)
        )
        ax.set_ylim(
            log_central_temperature_star - 0.1,
            log_central_temperature_star + 0.1
        )
        ax.set_xlabel(f'log (rho$_c$ / {unit_dens})')
        ax.set_ylabel(f'log (temp$_c$ / {unit_temp})')
        ax.ticklabel_format(style='sci', scilimits=(-3, 4), axis='both')
        self.__temp_dens_plot, = ax.plot(
            [],
            [],
        )
        self.__temp_dens_star = ax.scatter(
            [log_central_density_star],
            [log_central_temperature_star],
            facecolor='red',
            s=40,
            edgecolor=None,
            alpha=0.5,
        )
        self.__initialised.append(title)

    def plot_temp_dens(self):
        title = 'temp_dens'
        if title not in self.__initialised:
            self.initialise_temp_dens()
        unit_temp = self.default_units['temperature']
        unit_dens = self.default_units['density']
        ax = self.__axes[title]
        log_central_density = np.log10(
            self.__central_density.value_in(unit_dens)
        )
        log_central_temperature = np.log10(
            self.__central_temperature.value_in(unit_temp)
        )
        ax.set_xlim(
            min(log_central_density) - 0.1,
            max(log_central_density) + 0.1,
        )
        ax.set_ylim(
            min(log_central_temperature) - 0.1,
            max(log_central_temperature) + 0.1,
        )
        self.__temp_dens_plot.set_xdata(log_central_density)
        self.__temp_dens_plot.set_ydata(log_central_temperature)
        offsets = self.__temp_dens_star.get_offsets()
        offsets[:, 0] = [
            np.log10(self.star.central_density.value_in(unit_dens))
        ]
        offsets[:, 1] = [
            np.log10(self.star.central_temperature.value_in(unit_temp))
        ]
        self.__temp_dens_star.set_offsets(offsets)

    def initialise_central_abundance(self):
        title = 'central abundance'
        unit_age = self.default_units['age']
        ax = self.__axes[title]
        ax.set_title(title)

        ax.set_ylim(-0.05, 1.05)
        if (self.star.get_phase() > 1 and len(self.__age) > 2):
            ax.set_xlabel(f'log age/{unit_age}')
            time_xdata = np.log10(
                0.01 +
                (
                    1.0 * self.__age[-1].value_in(unit_age)
                    # - self.__age[-2].value_in(unit_age)
                )
                - self.__age.value_in(unit_age)
            )
        else:
            ax.set_xlabel(f'age ({unit_age})')
            time_xdata = self.__age.value_in(unit_age)
        xmax = max(time_xdata) + 0.05
        xmin = min(time_xdata)
        ax.set_xlim(xmin, xmax)

        ax.set_ylabel('abundance fraction')
        ax.ticklabel_format(
            style='sci',
            scilimits=(-3, 4),
            axis='both',
        )
        self.__central_abundance_plots = []
        for i, species in enumerate(self.__mainspecies):
            plot, = ax.plot(
                [-time_xdata],
                [self.__central_abundance[species]],
                label=species,
                # color='red',
                # linewidth=2,
            )
            self.__central_abundance_plots.append(plot)
        self.__initialised.append(title)

    def plot_central_abundance(self):
        title = 'central abundance'
        if title not in self.__initialised:
            self.initialise_central_abundance()
        unit_age = self.default_units['age']
        ax = self.__axes[title]
        if (self.star.get_phase() > 1 and len(self.__age) > 2):
            ax.set_xlabel(f'log age/{unit_age}')
            time_xdata = -np.log10(
                0.01 +
                (
                    1.0 * self.__age[-1].value_in(unit_age)
                    # - self.__age[-2].value_in(unit_age)
                )
                - self.__age.value_in(unit_age)
            )
        else:
            ax.set_xlabel(f'age ({unit_age})')
            time_xdata = self.__age.value_in(unit_age)

        xmax = max(time_xdata) + 0.05
        xmin = min(time_xdata)
        ax.set_xlim(xmin, xmax)
        for i, species in enumerate(self.__mainspecies):
            self.__central_abundance_plots[i].set_xdata(
                time_xdata
            )
            # if self.__phase == 1:
            #     self.__central_abundance_plots[i].set_xdata(
            #         self.__age.value_in(unit_age)
            #     )
            # else:
            #     self.__central_abundance_plots[i].set_xdata(
            #         # self.__age.value_in(unit_age)
            #         (self.__star.age - self.__age).value_in(unit_age)
            #     )
            self.__central_abundance_plots[i].set_ydata(
                self.__central_abundance[species.capitalize()]
            )
        # if self.__phase > 1:  # central_abundance[self.__species[0]] < 1e-4:
        #     ax.set_xscale('log')
        #     ax.set_xlim(1, xmax)
        # else:
        #     ax.set_xscale('linear')

    def initialise_kippenhahn(self):
        title = 'Kippenhahn diagram'
        unit_age = self.default_units['age']
        ax = self.__axes[title]
        ax.set_title(title)
        ax.set_xlim(
            -0.05,
            0.05,
        )
        ax.set_xlabel(f'age ({unit_age})')
        self.__initialised.append(title)

    def plot_kippenhahn(self):
        title = 'Kippenhahn diagram'
        if title not in self.__initialised:
            self.initialise_kippenhahn()
        unit_age = self.default_units['age']
        ax = self.__axes['Kippenhahn diagram']
        xmax = self.star.age.value_in(unit_age) * 1.05
        xmin = (self.star.age.value_in(unit_age) - xmax)
        ax.set_xlim(xmin, xmax)

    def initialise_abundance_profile(self):
        title = "Abundance profile"
        ax = self.__axes[title]
        ax.set_title(title)
        ax.set_xlim(0, 1)
        ax.set_ylim(-4.05, 0.05)
        ax.set_xlabel('M$_{\\rm r}$/M$_{\\rm tot}$')
        self.__initialised.append(title)
        self.__abundance_plots = []
        mass_profile = self.star.get_cumulative_mass_profile()
        for i, species in enumerate(self.__mainspecies):
            j = np.argmax(species.lower() == self.__species)
            plot, = ax.plot(
                mass_profile,
                np.log10(self.__abundances[j]),
                label=species,
            )
            self.__abundance_plots.append(plot)
        # ax.legend()
        self.__initialised.append(title)

    def plot_abundance_profile(self):
        title = "Abundance profile"
        if self.__redraw:
            ax = self.__axes[title]
            ax.clear()
        if (
            title not in self.__initialised
            or self.__redraw
        ):
            self.initialise_abundance_profile()
        mass_profile = self.star.get_cumulative_mass_profile()
        # Must reinitialise in some way if the number of zones has changed!
        # if len(mass_profile) != len(self.__abundance_plots[0].get_xdata()):
        #     self.initialise_abundance_profile()

        for i, species in enumerate(self.__mainspecies):
            j = np.argmax(species.lower() == self.__species)
            self.__abundance_plots[i].set_xdata(
                mass_profile
            )
            self.__abundance_plots[i].set_ydata(
                np.log10(self.__abundances[j])
            )

    def initialise_gradient(self):
        title = "gradient"
        ax = self.__axes[title]
        self.__gradient_plots = []
        ax.set_title(title)
        ax.set_xlim(0, 1)
        ax.set_xlabel('M$_{\\rm r}$/M$_{\\rm tot}$')
        ax.set_ylabel('Nabla_(rad-ad-mu)')
        mass_profile = self.star.get_cumulative_mass_profile()
        nabla_rad = self.star.get_nabla_rad_profile()
        nabla_ad = self.star.get_nabla_ad_profile()
        nabla_mu = self.star.get_nabla_mu_profile()
        min_nabla = min(
            min(nabla_ad),
            min(nabla_rad),
            min(nabla_mu),
            -0.1
        )
        max_nabla = max(
            max(nabla_ad),
            max(nabla_rad),
            # max(nabla_mu),
            0
        )
        plot_1, = ax.plot(mass_profile, nabla_ad, label='Nabla_ad')
        plot_2, = ax.plot(mass_profile, nabla_rad, label='Nabla_rad')
        plot_3, = ax.plot(mass_profile, nabla_mu, label='Nabla_mu', alpha=0.5)
        ax.plot(mass_profile, 0*mass_profile, linestyle='--', color='black')
        ax.set_ylim((min_nabla, min(max_nabla, 10)))
        self.__gradient_plots.append(plot_1)
        self.__gradient_plots.append(plot_2)
        self.__gradient_plots.append(plot_3)
        self.__initialised.append(title)

    def plot_gradient(self):
        title = "gradient"
        if self.__redraw:
            ax = self.__axes[title]
            ax.clear()
        if (
            title not in self.__initialised
            or self.__redraw
        ):
            self.initialise_gradient()
            return
        ax = self.__axes[title]
        mass_profile = self.star.get_cumulative_mass_profile()
        nabla_rad = self.star.get_nabla_rad_profile()
        nabla_ad = self.star.get_nabla_ad_profile()
        nabla_mu = self.star.get_nabla_mu_profile()
        min_nabla = min(
            min(nabla_ad),
            min(nabla_rad),
            min(nabla_mu),
            -0.05
        )
        max_nabla = max(
            max(nabla_ad),
            max(nabla_rad),
            0
        )
        ax.set_ylim((min_nabla, min(max_nabla, 10)))
        for i in range(3):
            self.__gradient_plots[i].set_xdata(
                mass_profile
            )
        self.__gradient_plots[0].set_ydata(nabla_ad)
        self.__gradient_plots[1].set_ydata(nabla_rad)
        self.__gradient_plots[2].set_ydata(nabla_mu)


    def initialise_luminosity_energy_production(self):
        title = "Luminosity and energy production"
        ax_lum = self.__axes[title][0]
        ax_energy = self.__axes[title][1]
        self.__lum_en_plots = []
        ax_lum.set_title(title)
        ax_lum.set_xlim(0, 1)
        ax_energy.set_xlim(0, 1)
        ax_lum.set_xlabel('M$_{\\rm r}$/M$_{\\rm tot}$')
        mass_profile = self.star.get_cumulative_mass_profile()
        luminosity_profile = self.star.get_luminosity_profile().value_in(
            units.LSun
        )
        luminosity_profile = luminosity_profile / max(luminosity_profile)
        plot_lum, = ax_lum.plot(
            mass_profile, luminosity_profile, label="lum", color="black"
        )
        ax_lum.set_ylim(
            (
                min(luminosity_profile)-0.1,
                max(luminosity_profile)+0.1,
            )
        )

        eps = self.star.get_eps_profile()
        eps[eps <= 0.] = 1e-32
        eps_H_log = np.log10(eps)
        epsy = self.star.get_epsy_profile()
        epsy[epsy <= 0.] = 1e-32
        eps_He_log = np.log10(epsy)
        max_eps = max(max(eps_H_log), max(eps_He_log))
        min_eps = max(min(eps_H_log), 2)

        ax_energy.set_ylim(
            (min_eps, max_eps)
        )
        plot_eps, = ax_energy.plot(
            mass_profile,
            eps_H_log, label="eps_H",
        )
        plot_epsy, = ax_energy.plot(
            mass_profile,
            eps_He_log, label="eps_He",
        )        
        self.__lum_en_plots.append(plot_lum)
        self.__lum_en_plots.append(plot_eps)
        self.__lum_en_plots.append(plot_epsy)
        self.__initialised.append(title)

    def plot_luminosity_energy_production(self):
        title = "Luminosity and energy production"
        if self.__redraw:
            ax_lum = self.__axes[title][0]
            ax_lum.clear()
            ax_energy = self.__axes[title][1]
            ax_energy.clear()
        if (
            title not in self.__initialised
            or self.__redraw
        ):
            self.initialise_luminosity_energy_production()
            return
        ax_lum = self.__axes[title][0]
        ax_energy = self.__axes[title][1]
        mass_profile = self.star.get_cumulative_mass_profile()
        luminosity_profile = self.star.get_luminosity_profile().value_in(
            units.LSun
        )
        luminosity_profile = luminosity_profile / max(luminosity_profile)

        for i, plot in enumerate(self.__lum_en_plots):
            plot.set_xdata(
                mass_profile
            )
        self.__lum_en_plots[0].set_ydata(
            luminosity_profile
        )
        ax_lum.set_ylim(
            (
                min(luminosity_profile)-0.1,
                max(luminosity_profile)+0.1,
            )
        )
        eps = self.star.get_eps_profile()
        eps[eps <= 0.] = 1e-32
        eps_H_log = np.log10(eps)
        epsy = self.star.get_epsy_profile()
        epsy[epsy <= 0.] = 1e-32
        eps_He_log = np.log10(epsy)
        max_eps = max(max(eps_H_log), max(eps_He_log))
        min_eps = max(min(eps_H_log), 2)

        ax_energy.set_ylim(
            (min_eps, max_eps)
        )
        for i, ydata in enumerate(
            [eps_H_log, eps_He_log, ]
        ):
            self.__lum_en_plots[1+i].set_ydata(
                ydata,
            )

    def initialise_temperature_pressure(self):
        title = "Temperature and pressure"
        ax = self.__axes[title][0]
        self.__temp_pres_plots = []
        ax.set_title(title)
        ax.set_xlim(0, 1)
        ax.set_xlabel('M$_{\\rm r}$/M$_{\\rm tot}$')
        ax.set_ylabel('log(T)')
        ax_pres = self.__axes[title][1]
        ax_pres.set_ylabel('log(P)')
        mass_profile = self.star.get_cumulative_mass_profile()
        temperature_profile = np.log10(
            self.star.get_temperature_profile().value_in(units.K)
        )
        pressure_profile = np.log10(
            self.star.get_pressure_profile().number
        )
        plot_temp, = ax.plot(
            mass_profile, temperature_profile, label="temp",
        )
        plot_pres, = ax_pres.plot(
            mass_profile, pressure_profile, label="pres",
            color='black'
        )
        self.__temp_pres_plots.append(
            plot_temp
        )
        self.__temp_pres_plots.append(
            plot_pres
        )
        self.__initialised.append(title)

    def plot_temperature_pressure(self):
        title = "Temperature and pressure"
        if self.__redraw:
            for i in (0, 1):
                ax = self.__axes[title][i]
                ax.clear()
        if (
            title not in self.__initialised
            or self.__redraw
        ):
            self.initialise_temperature_pressure()
            return
        ax = self.__axes[title][0]
        ax_pres = self.__axes[title][1]
        mass_profile = self.star.get_cumulative_mass_profile()
        temperature_profile = np.log10(
            self.star.get_temperature_profile().value_in(units.K)
        )
        pressure_profile = np.log10(
            self.star.get_pressure_profile().number
        )
        for i in (0, 1):
            self.__temp_pres_plots[i].set_xdata(
                mass_profile
            )
        self.__temp_pres_plots[0].set_ydata(
            temperature_profile
        )
        self.__temp_pres_plots[1].set_ydata(
            pressure_profile
        )


    def plot_figure_evolution(self):
        self.plot_hr()
        self.plot_temp_dens()
        self.plot_central_abundance()
        self.plot_kippenhahn()

    def plot_figure_structure(self):
        self.plot_abundance_profile()
        self.plot_gradient()
        self.plot_luminosity_energy_production()
        self.plot_temperature_pressure()

    def plot_all(self, speed=0 | units.Myr / units.minute, step=None):
        self.plot_figure_evolution()
        self.plot_figure_structure()

        self.draw(speed=speed, step=step)
        plt.tight_layout()

        for fig in self.__figures.values():
            fig.canvas.flush_events()
        self.__redraw = False
        plt.pause(0.01)
        # plt.show(block=False)


def main():
    # star = read_set_from_file(sys.argv[1])[0]
    # plotting = StellarModelPlot()

    return


if __name__ == "__main__":
    main()
