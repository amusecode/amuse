"""
plotting.py: File containing different routines to plot simulation results

@author: Stan Verhoeve, Anton Rijkers, and Jelmer Stroo
@credits: Stan Verhoeve, Anton Rijkers, and Jelmer Stroo
@maintainer: Stan Verhoeve
@email: verhoeve@strw.leidenuniv.nl
"""
from amuse.lab import *
from matplotlib import pyplot
from matplotlib import colormaps
from mpl_toolkits.axes_grid1 import make_axes_locatable
import numpy
import warnings
import sys
import importlib.util as ilu
from scipy.interpolate import RBFInterpolator as scipy_RBFI
from scipy.ndimage import gaussian_filter
from amuse.ext.molecular_cloud import new_ism_cube

# pyplot.rc('font', size=60/1.5)            # controls default text sizes
# pyplot.rc('axes', titlesize=96/1.5)     # fontsize of the axes title
# pyplot.rc('axes', labelsize=72/1.5)    # fontsize of the x and y labels
# pyplot.rc('xtick', labelsize=60/1.5)    # fontsize of the tick labels
# pyplot.rc('ytick', labelsize=60/1.5)    # fontsize of the tick labels
# pyplot.rc('legend', fontsize=72/1.5)    # legend fontsize
# pyplot.rc('figure', titlesize=96/1.5) 

# Check if torch and torchrbf are installed
if (ilu.find_spec("torch") is not None and ilu.find_spec("torchrbf") is not None):
    # If installed, import necessary functions
    from torch import Tensor
    from torchrbf import RBFInterpolator as torch_RBFI


def u_to_temperature(internal_energy):
    """
    Function that transforms internal energy to temperature
    (For now) assumes all gas is molecular hydrogen and gamma=5/3
    """
    temperature = (5/3 - 1) * (2.333|units.amu) * internal_energy / constants.kB

    return temperature

def create_subplot_figure(N, max_cols=2, base=15):
    """
    Function that generates subplot based on total number of plots needed

    Parameters:
    -----------
    N : int
        Total number of plots needed
    max_cols (optional) : int
        Maximum number of columns before adding a new row.
        The default is 2
    base (optional) : int or float
        Size of an individual subplot

    Returns:
    --------
    fig : matplotlib Figure object
        Figure object
    axs : ndarray
        Flattened ndarray of matplotlib Axes objects
    """
    # Calculate number of rows and columns
    if N==1:
        nrows = 1
        ncols = 1
    else:
        nrows = numpy.ceil(N / max_cols).astype(int)
        ncols = max(numpy.ceil(N / nrows).astype(int), max_cols)

    # Create figure and (flattened) axes
    fig, axs = pyplot.subplots(nrows, ncols, figsize=(base * ncols, base * nrows), squeeze=False)
    axs = axs.flatten()
    if N>1:
        # Delete all axes that will not be used
        if N%max_cols:
            for ax in axs[-(len(axs)-N):]:
                ax.remove()
    return fig, axs


class HydroPlotter:
    """
    Object to plot hydrodynamical systems. 

    Implemented functionality:
    - Plotting density slice in 2D plane
    - Plotting temperature slice in 2D plane
    - Plotting mass slice in 2D plane
    - Plotting ionization slice in 2D plane
    - Overplotting velocity field in the above slices
    - Plotting all of the above in a subplot
    - Custom vmin, vmax
        - Has a few caveats, but see examples (work in progress)
    - Custom supertitle
    - Ability to save plot (will not show)
    -  Interpolation is limited to 5000 particles
        - Significantly faster, and less memory overhead, but less accurate

    TODO
    - Add different plotting procedures
        - e.g. energy, turbulence etc.
    """

    def __init__(self, grid_resolution=100, use_torch=False):
        """
        Parameters
        ----------
        grid_resolution : int
            Number of grid points in each direction to project parameters onto
        """
        self.__grid_resolution = grid_resolution
        self.gas_particles = None
        self.star_particles = None
        self.__los_coordinate_lookup = {"yz":0, "xz":1, "xy":2}
        self.__grid_units = units.pc
        self.__use_torch = use_torch

        if self.__use_torch:
            if "torch" not in sys.modules:
                warnings.warn("Trying to use torchrbf, but could not import. Using SciPy instead")
                self.__use_torch = False

    @property
    def grid_resolution(self):
        return self.__grid_resolution

    @grid_resolution.setter
    def grid_resolution(new_res):
        self.__grid_resolution = new_res

    @property
    def grid_units(self):
        return self.__grid_units

    @grid_units.setter
    def grid_units(self, units):
        self.__grid_units = units
    
    def add_gas_particles(self, gas_particles):
        assert self.gas_particles == None, "Object already contains gas particles! Remove first with clear_gas_particles() or clear_all_particles()"
    
        self.gas_particles = gas_particles

    def add_star_particles(self, star_particles):
        assert self.star_particles == None, "Object already contains star particles! Remove first with clear_star_particles() or clear_all_particles()"
        self.star_particles = star_particles

    def clear_gas_particles(self):
        del self.gas_particles
        self.gas_particles = None
    
    def clear_star_particles(self):
        del self.star_particles
        self.star_particles = None

    def clear_all_particles(self):
        del self.gas_particles
        del self.star_particles
        self.gas_particles = None
        self.star_particles = None

    def _create_interpolated_map(self, slice_param, size, depth=0, los_index=2):
        """
        Function that creates a 2D density slice along a line of sight (los)

        Parameters
        ----------
        slice_param : ndarray
            ndarray of shape (N,) with known quantities to interpolate between
        size : int or float
            Total size of the plot in grid_units
        depth : int or float
            Depth along line of sight to take slice of. A depth of 
            0 corresponds to the center of `size`. 
            The default is 0
        los_index : int
            Index of the line of sight coordinate. A `los_index` of 0, 1, 2
            corresponds to a line of sight coordinate of `x`, `y`, `z`, respectively.
            The default is 2

        Returns
        -------
        interpolated_slice_param : ndarray
            ndarray of shape (grid_resolution, grid_resolution) containing interpolated
            quantities on the square grid
        """

        pos = self.gas_particles.position.value_in(self.grid_units)
        
        # For memory optimalisation, use only 5000 datapoints for interpolation.
        # This will reduce interpolation accuracy, but will significantly speed
        # up interpolation, as well as significantly improve memory utilisation
        # Since particle positions are stored randomly (i.e. not sorted), we can
        # grab every `jump`th particle and still achieve good (enough) coverage
        # of the plane for interpolation
        
        if len(pos) < 5000:
            jump = None
        else:
            jump = len(pos) // 5000
        pos = pos[::jump, :]

        # Create 2D grid for slice
        step = 1j * self.grid_resolution
        grid = numpy.mgrid[-0.5:0.5:step, -0.5:0.5:step] * size
        
        # Flatten grid
        grid = grid.reshape(2,-1)

        # Get line-of-sight coordinates @ depth
        line_of_sight = depth * numpy.ones(self.grid_resolution**2)
        
        # Insert line of sight coordinate at correct position
        grid = numpy.insert(grid, los_index, line_of_sight, axis=0).T
        
        if self.__use_torch:
            interpolator = torch_RBFI
            grid = Tensor(grid)
        else:
            interpolator = scipy_RBFI
        
        # Create grid interpolation object and interpolate parameters on the grid
        RBFI = interpolator(pos, slice_param[::jump])
        interpolated_slice_param = RBFI(grid)

        return interpolated_slice_param
        
    def plot_slice(self, size, depth=0, projection="xy", params_to_plot="density", \
                   plot_stars=False, plot_velocity=False, save_fig=None, \
                   vmin=None, vmax=None, cmap="jet", \
                   title=None, base=5, max_cols=2, species_index=0):
        """
        Parameters
        ----------
        size : int or float
            Total size of the plot in grid_units
        depth : int or float
            Depth along line of sight to take slice of. A depth of 
            0 corresponds to the center of `size`. 
            The default is 0
        projection : str
            Projection plane to plot onto. The default is "xy"
        params_to_plot : str
            particle parameters to plot. The default is "density"
        plot_stars : bool
            Whether to plot the position of the star
        plot_velocity : bool
            Whether to plot the gas velocity vector field on top
        save_fig : str
            Name of folder/file to save plot into
        vmin : float or iterable
            vmin to be passed to imshow(). If not provided, imshow() will determine vmin.
            If a float, assumes all subplots have same vmin. If iterable, vmin should have
            the same length as params_to_plot
        vmax : float or iterable
            vmax to be passed to imshow(). If not provided, imshow() will determine vmax.
            If a float, assumes all subplots have same vmax. If iterable, vmax should have
            the same length as params_to_plot
        cmap : str or Colormap
            Colormap for slice. The default is jet
        title : str
            Title for plot. Will always be a supertitle. If not provided, will not create a title
            The default is None
        base : int or float
            Base size of a subplot
            The default is 5
        max_cols : int
            Maximum number of columns to use when plotting subplots.
            The default is 2
        species_index : int
            index of the species of molecule in de array of abundances.
            The default is 0
        """

        assert self.gas_particles is not None, "Add particles to plotter objects using HydroPlotter.add_gas_particles() method!"
        if plot_stars and self.star_particles == None:
            warnings.warn("Trying to plot star positions without star particles. Add stars using add_star_particles() method")
            plot_stars = False
        
        # Make sure vmin, vmax have correct shape
        if hasattr(vmin, '__iter__'):
            assert len(vmin) == len(params_to_plot), "When vmin is iterable, its length must equal that of params_to_plot"
        else:
            vmin = [vmin] * len(params_to_plot)

        if hasattr(vmax, '__iter__'):
            assert len(vmax) == len(params_to_plot), "When vmax is iterable, its length must equal that of params_to_plot"
        else:
            vmax = [vmax] * len(params_to_plot)
        

        # Make params_to_plot iterable if it is a string
        if type(params_to_plot)==str:
            params_to_plot = [params_to_plot]

        # List to keep track of interpolation parameters and colorbar labels
        params_to_interpolate = list()
        cbar_labels = list()
        
        # Create subplot figure
        fig, axs = create_subplot_figure(len(params_to_plot))
        if title:
            fig.suptitle(title)

        for param in params_to_plot:
            match param:
                case "density":
                    slice_param = self.gas_particles.density.value_in(units.amu/units.cm**3)
                    cbar_label = "log density [$amu/cm^3$]"
                case "column_density":
                    pixel_area = ((size | self.grid_units) / self.grid_resolution)**2
                    slice_param = (self.gas_particles.mass / pixel_area).value_in(units.MSun/units.pc**2)
                    cbar_label = "log column density [$M_\odot/pc^2$]"
                case "temperature":
                    slice_param = self.gas_particles.u.value_in(units.kms**2)
                    # u = self.gas_particles.u
                    # slice_param = u_to_temperature(u).value_in(units.K)
                    cbar_label = "log temperature [$K$]"
                case "mass":
                    slice_param = self.gas_particles.mass.value_in(units.kg)
                    cbar_label = "log mass [$kg$]"
                case "ionization":
                    assert hasattr(self.gas_particles, "xion"), "Gas particles should have xion specified to plot ionization fraction!"
                    slice_param = 10**self.gas_particles.xion
                    cbar_label = "Ionization fraction"
                case "abundance":
                    slice_param = self.gas_particles.abundances[:,species_index]
                    cbar_label = "Abundance"
                case _:
                    raise ValueError(f"{slice_attr} is not a valid attribute to plot!")
            params_to_interpolate.append(slice_param)
            cbar_labels.append(cbar_label)
        
        
        if plot_velocity:
            # Grab velocity and add to params to be interpolated
            vel = self.gas_particles.velocity.value_in(units.kms).T
            params_to_interpolate.append(vel)
            
            # Create grid to plot velocity on
            grid_step = 1j * self.grid_resolution
            XX, YY = numpy.mgrid[-0.5:0.5:grid_step, -0.5:0.5:grid_step] * size
    
        
        # Create vertically stacked array
        params_to_interpolate = numpy.vstack(params_to_interpolate).T

        # Sort projection (e.g. 'zy' projection becomes 'yz')
        projection = "".join(sorted(projection))

        # Get the index of the line-of-sight coordinate
        los_index = self.__los_coordinate_lookup[projection]
        
        # Get interpolated image map @ depth along line-of-sight
        interpolated_image_maps = self._create_interpolated_map(params_to_interpolate, size, depth, los_index)
        
        if plot_velocity:
            interpolated_vel = interpolated_image_maps[:,len(params_to_plot):]
            interpolated_vel = interpolated_vel.reshape((self.grid_resolution, self.grid_resolution,3))
            
            # Do not consider axis along line of sight
            interpolated_vel = numpy.delete(interpolated_vel, los_index, axis=2)
        
        # Iterate over parameters to plot
        for i, param in enumerate(params_to_plot):
            # Grab corresponding image map
            interpolated_image = interpolated_image_maps[:,i]

            # If plotting temperature, convert
            if param == "temperature":
                print("Now converting to temp")
                interpolated_image = u_to_temperature(interpolated_image | (units.kms**2)).value_in(units.K)
           
            # Reshape and transpose image
            interpolated_image = interpolated_image.reshape((self.grid_resolution, self.grid_resolution)).T

            # Clean up interpolated data
            interpolated_image = numpy.where(interpolated_image < 1e-5, 1e-5, interpolated_image)
            

            # Plot image slice
            im = axs[i].imshow(numpy.log10(interpolated_image),
                               cmap=cmap,
                               extent=[-size/2,size/2,-size/2,size/2],
                               origin="lower",
                               vmin=vmin[i], 
                               vmax=vmax[i])
            
            # Create new ax to plot colorbar
            divider = make_axes_locatable(axs[i])
            cax = divider.append_axes("right", size="5%", pad=0.15)
            
            # Add colorbar to (sub)plot
            cbar = fig.colorbar(im, cax=cax, orientation="vertical")
            cbar.set_label(cbar_labels[i], rotation=270, labelpad=60)
            cbar.ax.tick_params(labelsize=20)
            
            axs[i].set_xlabel(f"{projection[0]} [{self.grid_units.unit}]")
            axs[i].set_ylabel(f"{projection[1]} [{self.grid_units.unit}]")
            
            axs[i].tick_params(axis='x', labelsize=20)
            axs[i].tick_params(axis='y', labelsize=20)

            if plot_stars:
                # Grab positions of star particles
                star_pos = numpy.atleast_2d(self.star_particles.position.value_in(self.grid_units)).T
            
                # Grab only stars within bounding plot range
                out_of_bound_idx = numpy.argwhere(numpy.any(abs(star_pos) > size/2, axis=0))
                star_pos = numpy.delete(star_pos, out_of_bound_idx, axis=1)
            
                # Only plot if there are actually stars inside frame
                if len(star_pos.flatten()) > 0:
                    star_pos_slice = numpy.delete(star_pos, los_index, axis=0)
                    star_pos_depth = star_pos[los_index]
                    
                    for hv, d in zip(star_pos_slice.T, star_pos_depth):
                        color = colormaps.get_cmap("Greys")(d)

                        axs[i].scatter(*hv, color=color, label=f"Star at z={d:.2f} pc", marker="*")
                    axs[i].legend(loc="upper right")

            if plot_velocity:
                # Plot every 10th vector in the field
                axs[i].quiver(XX[::10,::10], YY[::10,::10], \
                              interpolated_vel[::10,::10,0], \
                              interpolated_vel[::10,::10,1], \
                              color="white", angles="xy")

        pyplot.tight_layout()
        
        if save_fig:
            pyplot.savefig(save_fig, bbox_inches="tight", dpi=300)
            pyplot.close()
        else:
            pyplot.show()

        return

    def plot_projection(self, size, projection="xy", params_to_plot="density", \
                        plot_stars=False, use_gaussian_kernel=False,
                        save_fig=None, \
                        vmin=None, vmax=None, cmap="jet", \
                        title=None, base=5, max_cols=2, species_index=0):
        """
        Parameters
        ----------
        size : int or float
            Total size of the plot in grid_units
        projection : str
            Projection plane to plot onto. The default is "xy"
        params_to_plot : str
            particle parameters to plot. The default is "density"
        plot_stars : bool
            Whether to plot the position of the star
        use_gaussian_kernel : bool
            Whether to use a gaussian smoothing kernel when plotting
        save_fig : str
            Name of folder/file to save plot into
        vmin : float or iterable
            vmin to be passed to imshow(). If not provided, imshow() will determine vmin.
            If a float, assumes all subplots have same vmin. If iterable, vmin should have
            the same length as params_to_plot
        vmax : float or iterable
            vmax to be passed to imshow(). If not provided, imshow() will determine vmax.
            If a float, assumes all subplots have same vmax. If iterable, vmax should have
            the same length as params_to_plot
        cmap : str or Colormap
            Colormap for slice. The default is jet
        title : str
            Title for plot. Will always be a supertitle. If not provided, will not create a title
            The default is None
        base : int or float
            Base size of a subplot
            The default is 5
        max_cols : int
            Maximum number of columns to use when plotting subplots.
            The default is 2
        species_index : int
            index of the molecular species in the abundances array
            The default is 0
        """

        assert self.gas_particles is not None, "Add particles to plotter objects using HydroPlotter.add_gas_particles() method!"
        if plot_stars and self.star_particles == None:
            warnings.warn("Trying to plot star positions without star particles. Add stars using add_star_particles() method")
            plot_stars = False
        
        # Make sure vmin, vmax have correct shape
        if hasattr(vmin, '__iter__'):
            assert len(vmin) == len(params_to_plot), "When vmin is iterable, its length must equal that of params_to_plot"
        else:
            vmin = [vmin] * len(params_to_plot)

        if hasattr(vmax, '__iter__'):
            assert len(vmax) == len(params_to_plot), "When vmax is iterable, its length must equal that of params_to_plot"
        else:
            vmax = [vmax] * len(params_to_plot)
        

        # Make params_to_plot iterable if it is a string
        if type(params_to_plot)==str:
            params_to_plot = [params_to_plot]
        
        # Sort projection (e.g. 'zy' projection becomes 'yz')
        projection = "".join(sorted(projection))
        
        # Get the index of the line-of-sight coordinate
        los_index = self.__los_coordinate_lookup[projection]
        
        # Grab position and delete los_index
        pos = self.gas_particles.position.value_in(self.grid_units)
        pos = numpy.delete(pos, los_index, axis=1)

        # Create bins to plot projection on
        bins = numpy.linspace(-0.5, 0.5, self.grid_resolution) * size
        
        # Create subplot figure
        fig, axs = create_subplot_figure(len(params_to_plot))
        if title:
            fig.suptitle(title)

        for i, param in enumerate(params_to_plot):
            match param:
                case "density":
                    projection_param = self.gas_particles.density.value_in(units.amu/units.cm**3)
                    cbar_label = "log density [$amu/cm^3$]"
                case "column_density":
                    projection_param = self.gas_particles.mass.value_in(units.MSun)
                    cbar_label = "log column density [$M_\odot/pc^2$]"
                case "temperature":
                    u = self.gas_particles.u
                    projection_param = u_to_temperature(u).value_in(units.K)
                    cbar_label = "log temperature [$K$]"
                case "mass":
                    projection_param = self.gas_particles.mass.value_in(units.kg)
                    cbar_label = "log mass [$kg$]"
                case "ionization":
                    assert hasattr(self.gas_particles, "xion"), "Gas particles should have xion specified to plot ionization fraction!"
                    projection_param = self.gas_particles.xion
                    cbar_label = "Ionization fraction"
                case "abundance":
                    projection_param = self.gas_particles.abundances[:,species_index]
                    cbar_label = "abundance"
                case _:
                    raise ValueError(f"{param} is not a valid attribute to plot!")
            
            pixel_area = (size / len(bins))**2
            # Create projection
            projected_image, *__ = numpy.histogram2d(*pos.T, bins=bins, weights=projection_param/pixel_area)

            if use_gaussian_kernel:
                projected_image = gaussian_filter(projected_image, sigma=2, order=0)
            
            # Cleanly take log10 and reshape
            projected_image = numpy.where(projected_image<1e-5, 1e-5, projected_image)
            
            log10_projected_image = numpy.log10(projected_image)
            log10_projected_image = log10_projected_image.reshape(projected_image.shape).T
            
            print(numpy.max(log10_projected_image))
            # Plot image slice
            im = axs[i].imshow(log10_projected_image,
                               cmap=cmap,
                               extent=[-size/2,size/2,-size/2,size/2],
                               origin="lower",
                               vmin=vmin[i], 
                               vmax=vmax[i])
            
            # Create new ax to plot colorbar
            divider = make_axes_locatable(axs[i])
            cax = divider.append_axes("right", size="5%", pad=0.15)
            
            # Add colorbar to (sub)plot
            cbar = fig.colorbar(im, cax=cax, orientation="vertical")
            cbar.set_label(cbar_label, rotation=270, labelpad=60)
            
            axs[i].set_xlabel(f"{projection[0]} [{self.grid_units.unit}]")
            axs[i].set_ylabel(f"{projection[1]} [{self.grid_units.unit}]")

            if plot_stars:
                # Grab positions of star particles
                star_pos = numpy.atleast_2d(self.star_particles.position.value_in(self.grid_units)).T
            
                # Grab only stars within bounding plot range
                out_of_bound_idx = numpy.argwhere(numpy.any(abs(star_pos) > size/2, axis=0))
                star_pos = numpy.delete(star_pos, out_of_bound_idx, axis=1)
            
                # Only plot if there are actually stars inside frame
                if len(star_pos.flatten()) > 0:
                    star_pos_slice = numpy.delete(star_pos, los_index, axis=0)
                    star_pos_depth = star_pos[los_index]
                    
                    for hv, d in zip(star_pos_slice.T, star_pos_depth):
                        color = colormaps.get_cmap("Greys")(d)
                        axs[i].scatter(*hv, color=color, label=f"Star at depth={d:.2f} pc", marker="*")
                    axs[i].legend(loc="upper right")
        pyplot.tight_layout()
        
        if save_fig:
            pyplot.savefig(save_fig, bbox_inches="tight", dpi=300)
            pyplot.close()
        else:
            pyplot.show()

        return



def histogram(particles_list, particles0=None, bins=20, time=0, title="", 
                              cumulative=False, save=None, attribute='density', 
                              impact=0, logx=False, logy=False):
    """
    Function to create a histogram of an SPH particle set.

    Parameters:
    ----------
    particles: particle or array of particle sets set to be plotted.
    particles0: particleset at time=0 to compare with
    n_bins: number of bins to be plotted.
    attribute: attribute of the particles to be plotted.
    """     
    if type(particles_list) != list:
        particles_list = [particles_list]

    if particles0:
        attribute_list0 = getattr(particles0, attribute)
        attribute_unit0 = attribute_list0.unit
        attribute_list0 = attribute_list0.value_in(attribute_unit0)

    
    pyplot.figure(figsize=(15,15))
    impacts = [0,3,4,5,6,7]
    for i, particles in enumerate(particles_list):
        attribute_list = getattr(particles, attribute)
        attribute_unit = attribute_list.unit
        attribute_list = attribute_list.value_in(attribute_unit)

        if cumulative:
            ylabel = "Probability"
            # pyplot.ecdf(attribute_list, label=f"impact={(impacts[i])/10:.1f}R", lw=3)
            pyplot.ecdf(attribute_list, label=f"Binary: {i==0}", lw=3)
        else:
            ylabel = "PDF"
            hist, edges, *_ = pyplot.hist(attribute_list,  bins=bins, 
                                                           edgecolor='k', 
                                                           color='lightgray', 
                                                           label=f'impact={(impacts[i])/10:.1f}R', 
                                                           density=True)


    if particles0:
        if cumulative:
            pyplot.ecdf(attribute_list0, color="k", ls="dashdot", label="t=0", lw=3)
        else:
            pyplot.hist(attribute_list0, bins=bins, 
                                         edgecolor="k", 
                                         color="k", 
                                         label="t=0", 
                                         alpha=0.25,
                                         density=True)
    
    if logx:
        pyplot.xscale("log")
        pyplot.xlim(1e-20,1e-16)
    if logy:
        pyplot.yscale("log")

    pyplot.xlabel(f'{attribute} in {attribute_unit}')
    pyplot.ylabel(ylabel)
    pyplot.title(title)

    pyplot.legend()
    if save:
        pyplot.savefig(save, dpi=300)
        pyplot.close("all")
    pyplot.show()

def get_energy(particles):
    """
    Function to get the energy levels from the given
    particle set

    parameters:
    ----------
    particles : Amuse Particle(s)

    returns:
    ----------
    thermal : ndarray
        Thermal energy
    kinetic : ndarray
        Kinetic energy
    total : ndarray
        Sum of thermal and kinetic energy

    """
    thermal = particles.thermal_energy()
    kinetic = particles.kinetic_energy()
    total = thermal + kinetic

    return (thermal, kinetic, total)


def plot_energy(times, title, save, **energies):
    """
    Function to plot the flow of energy with time
    for potential, kinetic and total energy

    parameters:
    ----------
    times : Amuse units of time (units.kyr)
        list of times

    title: str
        Title for the plot
    save : str or None
        Path to file to save into
    **energies : dict
        Different keyword argument energies
    """

    pyplot.figure(figsize=(15,15))
    
    for energy in energies.keys():
        pyplot.plot(times.value_in(units.kyr), energies[energy].value_in(units.LSun * units.kyr), label=energy)
    
    pyplot.xlabel('Time [kyr]')
    pyplot.ylabel('Energy [$L_\odot\cdot$ kyr]')
    pyplot.title(title)
    pyplot.legend()
    
    if save:
        pyplot.savefig(save, dpi=300)
        pyplot.close("all")
    pyplot.show()


def plot_mass(times, mass, title="", save=None):
    """
    Function to plot the flow of mass with time
    for particles contained within a certain radius

    parameters:
    ----------
    times : list or ndarray
        Times in kyr

    mass : list or ndarray
        Mass contained within a region in Msun

    title: str
        Title for the plot

    save : str or None
        File name to save plot into

    """

    pyplot.figure(figsize=(15,15))
    pyplot.plot(times, mass)
    pyplot.xlabel('Time [kyr]')
    pyplot.ylabel('Mass [MSun]')
    pyplot.title(title)
    pyplot.legend()
    if save:
        pyplot.savefig(save, dpi=300)
        pyplot.close("all")
    pyplot.show()
