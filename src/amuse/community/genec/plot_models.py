import numpy
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
        self.__age = [] | self.default_units['age']
        self.__mass = [] | self.default_units['mass']
        self.__radius = [] | self.default_units['radius']
        self.__luminosity = [] | self.default_units['luminosity']
        self.__temperature = [] | self.default_units['temperature']
        self.__central_density = [] | self.default_units['density']
        self.__central_temperature = [] | self.default_units['temperature']
        self.__central_abundance = {}
        self.__species = star.get_names_of_species()
        # self.__zones = star.get_number_of_zones()
        self.__abundances = star.get_chemical_abundance_profiles()
        for species in self.__species:
            self.__central_abundance[species] = []

        self.__figures = {
            "1": plt.figure(num="Evolution", figsize=(8, 6)),
            "2": plt.figure(num="Structure", figsize=(8, 6)),
        }
        for figure in self.__figures.values():
            figure.set_tight_layout(True)

        self.__initialised = []

        self.__axes = {
            "HR": self.__figures["1"].add_subplot(2, 2, 1),
            "temp_dens": self.__figures["1"].add_subplot(2, 2, 2),
            "central abundance": self.__figures["1"].add_subplot(2, 2, 3),
            "Kippenhahn diagram": self.__figures["1"].add_subplot(2, 2, 4),
            "Abundance profile": self.__figures["2"].add_subplot(2, 2, 1),
            "gradient": self.__figures["2"].add_subplot(2, 2, 2),
            "Luminosity and energy production": self.__figures["2"].add_subplot(2, 2, 3),
            "Temperature and pressure": self.__figures["2"].add_subplot(2, 2, 4),
        }
        # self.__ax_hr = fig.add_subplot(2, 2, 1)
        self.update(star)
        self.plot_all()

    def update(self, star):
        self.star = star
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
            self.__central_abundance[species].append(
                self.__abundances[i, 0]
            )

    def draw(self, verbose=False, speed=0 | units.Myr / units.minute, step=None):
        for i, fig in self.__figures.items():
            if verbose:
                print(f"drawing figure {i}")

            fig.suptitle(
                f"age: {self.star.age.in_(units.Myr)}, "
                f"mass: {self.star.mass.in_(units.MSun)}, "
                f"radius: {self.star.radius.in_(units.RSun)}, "
                f"step: {step}, "
                # f"log(speedup) = {numpy.log10(speed):.2f}"
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
        log_tstar = numpy.log10(
            self.star.temperature.value_in(unit_temp)
        )
        ax.set_xlim(
            log_tstar + 0.1,
            log_tstar - 0.1,
        )
        log_lstar = numpy.log10(
            self.star.luminosity.value_in(unit_lum)
        )
        ax.set_ylim(log_lstar - 0.1, log_lstar + 0.1)
        ax.set_xlabel(f'log (temperature / {unit_temp})')
        ax.set_ylabel(f'log (luminosity / {unit_lum})')
        ax.ticklabel_format(style='sci',scilimits=(-3,4),axis='both')
        self.__hr_plot, = ax.plot(
            [numpy.log10(self.star.temperature.value_in(unit_temp))],
            [numpy.log10(self.star.luminosity.value_in(unit_lum))],
            # color='red',
            # linewidth=2,
        )
        self.__hr_star = ax.scatter(
            [numpy.log10(self.star.temperature.value_in(unit_temp))],
            [numpy.log10(self.star.luminosity.value_in(unit_lum))],
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
        log_temperature = numpy.log10(
            self.__temperature.value_in(unit_temp)
        )
        log_luminosity = numpy.log10(
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
        offsets[:, 0] = [numpy.log10(self.star.temperature.value_in(unit_temp))]
        offsets[:, 1] = [numpy.log10(self.star.luminosity.value_in(unit_lum))]
        self.__hr_star.set_offsets(offsets)

    def initialise_temp_dens(self):
        title = 'temp_dens'
        unit_temp = self.default_units['temperature']
        unit_dens = self.default_units['density']
        ax = self.__axes[title]
        ax.set_title(title)
        log_central_density_star = numpy.log10(
            self.star.central_density.value_in(unit_dens)
        )
        ax.set_xlim(
            log_central_density_star - 0.1,
            log_central_density_star + 0.1
        )
        log_central_temperature_star = numpy.log10(
            self.star.central_temperature.value_in(unit_temp)
        )
        ax.set_ylim(
            log_central_temperature_star - 0.1,
            log_central_temperature_star + 0.1
        )
        ax.set_xlabel(f'log (rho$_c$ / {unit_dens})')
        ax.set_ylabel(f'log (temp$_c$ / {unit_temp})')
        ax.ticklabel_format(style='sci',scilimits=(-3,4),axis='both')
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
        log_central_density = numpy.log10(
            self.__central_density.value_in(unit_dens)
        )
        log_central_temperature = numpy.log10(
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
            numpy.log10(self.star.central_density.value_in(unit_dens))
        ]
        offsets[:, 1] = [
            numpy.log10(self.star.central_temperature.value_in(unit_temp))
        ]
        self.__temp_dens_star.set_offsets(offsets)

    def initialise_central_abundance(self):
        title = 'central abundance'
        unit_age = self.default_units['age']
        ax = self.__axes[title]
        ax.set_title(title)

        ax.set_xlim(
            -0.05,
            0.05,
        )
        ax.set_ylim(-0.05, 1.05)
        ax.set_xlabel(f'age ({unit_age})')
        ax.set_ylabel('abundance fraction')
        ax.ticklabel_format(
            style='sci',
            scilimits=(-3, 4),
            axis='both',
        )
        self.__central_abundance_plots = []
        for i, species in enumerate(self.__species):
            plot, = ax.plot(
                [self.star.age.value_in(unit_age)],
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
        xmax = self.star.age.value_in(unit_age) * 1.05
        xmin = (self.star.age.value_in(unit_age) - xmax)
        ax.set_xlim(xmin, xmax)
        for i, species in enumerate(self.__species):
            self.__central_abundance_plots[i].set_xdata(
                self.__age.value_in(unit_age)
            )
            self.__central_abundance_plots[i].set_ydata(
                self.__central_abundance[species]
            )
        if self.__central_abundance[self.__species[0]] < 1e-4:
            ax.set_xscale('log')
            ax.set_xlim(1, xmax)
        else:
            ax.set_xscale('linear')

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
        ax.set_ylim(-4, 0)
        ax.set_xlabel('M$_{\\rm r}$/M$_{\\rm tot}$')
        self.__initialised.append(title)
        self.__abundance_plots = []
        mass_profile = self.star.get_cumulative_mass_profile()
        for i, species in enumerate(self.__species):
            plot, = ax.plot(
                mass_profile,
                numpy.log10(self.__abundances[i]),
                label=species,
            )
            self.__abundance_plots.append(plot)
        self.__initialised.append(title)

    def plot_abundance_profile(self):
        title = "Abundance profile"
        if title not in self.__initialised:
            self.initialise_abundance_profile()
        mass_profile = self.star.get_cumulative_mass_profile()
        for i, species in enumerate(self.__species):
            self.__abundance_plots[i].set_xdata(
                mass_profile
            )
            self.__abundance_plots[i].set_ydata(
                numpy.log10(self.__abundances[i])
            )

    def initialise_gradient(self):
        title = "gradient"
        ax = self.__axes[title]
        ax.set_title(title)
        ax.set_xlim(0, 1)
        ax.set_xlabel('M$_{\\rm r}$/M$_{\\rm tot}$')
        mass_profile = self.star.get_cumulative_mass_profile()
        self.__initialised.append(title)

    def plot_gradient(self):
        title = "gradient"
        if title not in self.__initialised:
            self.initialise_gradient()

    def initialise_luminosity_energy_production(self):
        title = "Luminosity and energy production"
        ax = self.__axes[title]
        ax.set_title(title)
        ax.set_xlim(0, 1)
        ax.set_xlabel('M$_{\\rm r}$/M$_{\\rm tot}$')
        mass_profile = self.star.get_cumulative_mass_profile()
        self.__initialised.append(title)

    def plot_luminosity_energy_production(self):
        title = "Luminosity and energy production"
        if title not in self.__initialised:
            self.initialise_luminosity_energy_production()

    def initialise_temperature_pressure(self):
        title = "Temperature and pressure"
        ax = self.__axes[title]
        ax.set_title(title)
        ax.set_xlim(0, 1)
        ax.set_xlabel('M$_{\\rm r}$/M$_{\\rm tot}$')
        mass_profile = self.star.get_cumulative_mass_profile()
        self.__initialised.append(title)

    def plot_temperature_pressure(self):
        title = "Temperature and pressure"
        if title not in self.__initialised:
            self.initialise_temperature_pressure()

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
        plt.pause(0.01)
        # plt.show(block=False)
