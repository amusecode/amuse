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
        for species in self.__species:
            self.__central_abundance[species] = []

        self.__figures = {
            "1": plt.figure(figsize=(12, 8)),
            # "2": plt.figure(),
        }

        self.__initialised_hr = False
        self.__initialised_temp_dens = False
        self.__initialised_central_abundance = False
        self.__initialised_kippenhahn = False

        self.__axes = {
            "HR": self.__figures["1"].add_subplot(2, 2, 1),
            "temp_dens": self.__figures["1"].add_subplot(2, 2, 2),
            "central abundance": self.__figures["1"].add_subplot(2, 2, 3),
            "Kippenhahn diagram": self.__figures["1"].add_subplot(2, 2, 4),
        }
        # self.__ax_hr = fig.add_subplot(2, 2, 1)
        self.update(star)
        self.plot_all(initialise=True)

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
        abundances = star.get_chemical_abundance_profiles()
        for i, species in enumerate(self.__species):
            self.__central_abundance[species].append(
                abundances[i, 0]
            )

    def draw(self, verbose=False, speed=0 | units.Myr / units.minute):
        for i, fig in self.__figures.items():
            if verbose:
                print(f"drawing figure {i}")

            fig.suptitle(
                f"star age = {self.star.age.in_(units.Myr)}, "
                f"radius = {self.star.radius.in_(units.RSun)}, "
                # f"log(speedup) = {numpy.log10(speed):.2f}"
                f"speedup = {speed.in_(units.Myr / units.minute)}"
            )
            fig.canvas.draw_idle()
            fig.canvas.flush_events()

    def initialise_hr(self):
        unit_temp = self.default_units['temperature']
        unit_lum = self.default_units['luminosity']
        ax = self.__axes['HR']
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
        ax.set_title('HR')
        self.__initialised_hr = True

    def plot_hr(self):
        if not self.__initialised_hr:
            self.initialise_hr()
        unit_temp = self.default_units['temperature']
        unit_lum = self.default_units['luminosity']
        ax = self.__axes['HR']
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
        unit_temp = self.default_units['temperature']
        unit_dens = self.default_units['density']
        ax = self.__axes['temp_dens']
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
        ax.set_title('temp/dens')
        self.__initialised_temp_dens = True

    def plot_temp_dens(self):
        if not self.__initialised_temp_dens:
            self.initialise_temp_dens()
        unit_temp = self.default_units['temperature']
        unit_dens = self.default_units['density']
        ax = self.__axes['temp_dens']
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
        offsets[:, 0] = [numpy.log10(self.star.central_density.value_in(unit_dens))]
        offsets[:, 1] = [numpy.log10(self.star.central_temperature.value_in(unit_temp))]
        self.__temp_dens_star.set_offsets(offsets)

    def initialise_central_abundance(self):
        unit_age = self.default_units['age']
        ax = self.__axes['central abundance']

        ax.set_xlim(
            -0.05,
            0.05,
        )
        ax.set_ylim(-0.05, 1.05)
        ax.set_xlabel(f'age {unit_age})')
        ax.set_ylabel(f'abundance fraction')
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
        ax.set_title('central abundance')
        self.__initialised_central_abundance = True

    def plot_central_abundance(self):
        if not self.__initialised_central_abundance:
            self.initialise_central_abundance()
        unit_age = self.default_units['age']
        ax = self.__axes['central abundance']
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

    def initialise_kippenhahn(self):
        unit_age = self.default_units['age']
        ax = self.__axes['Kippenhahn diagram']
        ax.set_title('Kippenhahn diagram')
        ax.set_xlim(
            -0.05,
            0.05,
        )
        ax.set_xlabel(f'age {unit_age})')
        self.__initialised_kippenhahn = True

    def plot_kippenhahn(self):
        if not self.__initialised_kippenhahn:
            self.initialise_kippenhahn()
        unit_age = self.default_units['age']
        ax = self.__axes['Kippenhahn diagram']
        xmax = self.star.age.value_in(unit_age) * 1.05
        xmin = (self.star.age.value_in(unit_age) - xmax)
        ax.set_xlim(xmin, xmax)

    def plot_all(self, initialise=False, speed=0 | units.Myr / units.minute):
        self.plot_hr()
        self.plot_temp_dens()
        self.plot_central_abundance()
        self.plot_kippenhahn()

        self.draw(speed=speed)
        plt.tight_layout()
        plt.pause(0.01)
