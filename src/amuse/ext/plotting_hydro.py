import os
import os.path
from amuse.units.quantities import zero
from amuse.plot import native_plot, sph_particles_plot


def new_plotting_hydrodynamics_code(hydrodynamics, timestep, plot_function=None, plot_function_arguments=dict(), plot_directory=None):
    """
    Returns a new subclass of the hydrodynamics code that will produce plots at
    regular intervals.
    
    :argument hydrodynamics: SPH code class to inherit from
    :argument timestep: interval for plotting (in )
    :argument plot_function: function that will be called after each timestep
    :argument plot_function_arguments: dict containing keyword arguments to the plot function
    """
    
    class PlottingHydrodynamics(hydrodynamics):
        
        _timestep = timestep
        _plot_function = staticmethod(plot_function)
        _plot_function_arguments = plot_function_arguments
        _plot_directory = plot_directory
        
        def __init__(self, *args, **kwargs):
            super(PlottingHydrodynamics, self).__init__(*args, **kwargs)
            self.time_last_plot = zero
            self.previous_plot_number = -1
            self.plot_directory = self._next_plot_directory()
            if self._plot_function is None:
                self._plot_function = sph_particles_plot
        
        def _next_plot_directory(self):
            if self._plot_directory is None:
                plot_directory = os.path.join(os.getcwd(), "plots")
                if not os.path.exists(plot_directory):
                    os.mkdir(plot_directory)
                i = 0
                while os.path.exists(os.path.join(plot_directory, "run_{0:=03}".format(i))):
                    i += 1
                new_directory = os.path.join(plot_directory, "run_{0:=03}".format(i))
                os.mkdir(new_directory)
                return new_directory
            else:
                if not os.path.exists(self._plot_directory):
                    os.mkdir(self._plot_directory)
                return self._plot_directory
        
        def _next_filename(self):
            self.previous_plot_number += 1
            return os.path.join(self.plot_directory, "hydroplot_{0:=04}".format(self.previous_plot_number))
        
        def evolve_model(self, end_time):
            time_next_plot = self.time_last_plot + self._timestep
            while time_next_plot < end_time:
                super(PlottingHydrodynamics, self).__getattr__("evolve_model")(time_next_plot)
                self._plot_function(self.gas_particles, **self._plot_function_arguments)
                native_plot.savefig(self._next_filename())
                native_plot.close()
                self.time_last_plot = time_next_plot
                time_next_plot = self.time_last_plot + self._timestep
            super(PlottingHydrodynamics, self).__getattr__("evolve_model")(end_time)
        
    PlottingHydrodynamics.__name__ = "Plotting" + hydrodynamics.__name__
    return PlottingHydrodynamics


