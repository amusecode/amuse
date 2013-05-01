from amuse.community.athena.interface import Athena
from amuse.units import generic_unit_converter
from amuse.units import generic_unit_system
from amuse.units import units
from amuse.units import constants
from amuse.test import amusetest
import numpy


try:
    from amuse import plot
    from matplotlib import pyplot
    IS_PLOT_AVAILABLE = True
except ImportError:
    IS_PLOT_AVAILABLE = False
    
dyn = units.g * units.cm / units.s**2

RUN_WITH_SELF_GRAVITY = True
 
class JeansInstability(object):

    def __init__(self, rho0, p0, wavenumber_factor = 4.0, pertubation_amplitude = 1e-3, gamma = 5.0/3.0, ncells = 100):
        self.rho0 = rho0
        self.p0 = p0
        self.pertubation_amplitude = pertubation_amplitude
        self.gamma = gamma
        self.ncells = ncells
        self.dimensions_of_mesh = (ncells, ncells, 1)
        
        self.sound_speed = ( self.gamma * self.p0 / self.rho0).sqrt()
        self.jeans_wavenumber = (4 * numpy.pi * constants.G * self.rho0).sqrt() / self.sound_speed
        self.wavenumber = wavenumber_factor * self.jeans_wavenumber
        self.length_scale = 0.5 * (numpy.pi * self.gamma * self.p0 / (constants.G * self.rho0**2)).sqrt()
        
        self.unit_converter = generic_unit_converter.ConvertBetweenGenericAndSiUnits(
            self.length_scale,
            1 | units.s,
            1 | units.g            
        )
        
    def new_code(self):
        if RUN_WITH_SELF_GRAVITY:
            result = Athena(self.unit_converter, mode = 'self-gravity')
            result.parameters.gamma = self.gamma
            result.parameters.courant_number=0.3
            result.parameters.four_pi_G = numpy.pi * 4 * constants.G
            result.parameters.gravity_mean_rho = self.rho0
        else:
            result = Athena(self.unit_converter)
            result.parameters.gamma = self.gamma
            result.parameters.courant_number=0.3
        
        return result
    
        
    def set_parameters(self):
        
        self.code.parameters.mesh_size = self.dimensions_of_mesh
        
        self.code.parameters.length_x = self.length_scale
        self.code.parameters.length_y = self.length_scale
        self.code.parameters.length_z = self.length_scale
        
        self.code.parameters.x_boundary_conditions = ("periodic","periodic")
        self.code.parameters.y_boundary_conditions = ("periodic","periodic")
        self.code.parameters.z_boundary_conditions = ("periodic","periodic")
        
        result = self.code.commit_parameters()
    
    def initialize_grid(self, grid):
        grid.momentum = [0.0, 0.0, 0.0]| units.g  / (units.s * units.cm**2)
        
        k_dot_x = self.wavenumber *  grid.x
        grid.rho = self.rho0 * (1 + self.pertubation_amplitude * numpy.cos(k_dot_x)) 
        pressure = self.p0 * (1 + self.pertubation_amplitude * self.gamma * numpy.cos(k_dot_x)) 
        grid.energy = pressure/ (self.gamma - 1)

        
    def setup(self):
        self.code = self.new_code()
        self.set_parameters()
        self.initialize_grid(self.code.grid)
        #self.update_potential_grid()
        
    def update_potential_grid(self):
        potential_grid = self.code.potential_grid
        potential_grid_inmem = potential_grid.copy()
    
        self.calculate_potential_for_grid(potential_grid_inmem, self.code.grid)
        channel = potential_grid_inmem.new_channel_to(potential_grid)
        channel.copy()
        
    def evolve_model(self, time):
        
        self.code.stopping_conditions.number_of_steps_detection.disable()
        
        self.code.stopping_conditions.number_of_steps_detection.enable()
        epsilon = 1e-12 | units.s
        while self.code.model_time + epsilon < time:
            self.code.evolve_model(time)
            
            if 0: # plot for debugging
                phi  =  self.code.grid.gravitational_potential[:,:,0]
                figure = pyplot.figure(figsize=(10,5))
                plot = figure.add_subplot(1,1,1)
                plot.imshow(phi.value_in(units.m**2 / units.s**2))
                pyplot.show()
                print phi[:,1][0:10]
                self.P0 = phi
            
            #self.update_potential_grid()
            
            
    def stop(self):
        self.code.stop()
    
    def get_gravity_for_grid(self, grid):
        """Calculates the gravitational force at each grid point,
        we do not use a code here (altoug possible)"""
        
        x = grid.x 
        y = grid.y
        z = grid.z                
        
        cell_size = (self.length_scale / self.ncells)
        cell_volume = cell_size * cell_size * (1 | units.cm)
        epsilon_squared = cell_size ** 2
        for cell in grid:
            dx = x - cell.x
            dy = x - cell.y
            dz = x - cell.z
     
            r_squared = (dx**2 + dy**2 + dz**2 + epsilon_squared)
            r = r_squared.sqrt()
            r_3 = r * r_squared       
            force = grid.rho * volume / r_3
            ax = (force * dx).sum()
            ay = (force * dy).sum()
            az = (force * dz).sum()
    
    def calculate_potential_for_grid(self, grid, field):
        x = field.x.flatten()
        
        x = field.x[:,0,0]
        y = field.y[0,:,0]
        
        # we assume constant grid spacing so
        # we just create a sample space using
        # the number of points in x and y
        
        nx = len(x)
        ny = len(y)
        dx = x[1] - x[0]
        dy = y[1] - y[0]
        
        #
        # wavenumbers from online examples
        # kx = (2 * numpy.pi) * numpy.fft.fftfreq(nx, d = dx ) * (1 | units.m)
        # ky = (2 * numpy.pi) * numpy.fft.fftfreq(ny, d = dy ) * (1 | units.m)
        #
        
        # wavenumbers as calculated in athena
        dkx = 2.0*numpy.pi/nx;
        dky = 2.0*numpy.pi/ny;
        kx = ((((2.0*numpy.cos(numpy.arange(nx) * dkx))-2.0)/(dx**2))).value_in(units.m**-2)
        ky = ((((2.0*numpy.cos(numpy.arange(ny) * dky))-2.0)/(dy**2))).value_in(units.m**-2)
        
        kx_grid , ky_grid = numpy.meshgrid(kx,ky)
        
        # convert rho to 2d field
        rho = field.rho[:,:,0]
        
        # use rho0 as mean, should be equal to rho.mean()
        # but we are comparing with athena and
        # want the most accurate solution
        rho_mean = self.rho0 #field.rho.mean()
        
        # to remove division by zero
        kx_grid[0][0] = 1
        ky_grid[0][0] = 1
        
        # for poisson solver the source term (gravity)
        # has to sum to zero over the field
        # we can ensure this by subtracting the mean from rho
        rho -= rho_mean
        
        gravity_function = 4 * numpy.pi * constants.G * rho
        
        gravity_function = gravity_function.value_in(units.s ** -2)
        gravity_function_fourier_space = numpy.fft.fftn(gravity_function)
        phi_fourier_space = gravity_function_fourier_space / (kx_grid + ky_grid)
        
        # 0,0 was divided by zero, should just be zero
        phi_fourier_space[0,0] = 0
        
        phi = numpy.fft.ifftn(phi_fourier_space).real

        # we removed the units for fft, replace these
        phi = phi |units.m**2 / units.s**2
        
        if 0: # plot for debugging
            phi_delta = phi# - self.P0
            figure = pyplot.figure(figsize=(10,5))
            plot = figure.add_subplot(1,2,1)
            plot.imshow(phi_delta.value_in(units.m**2 / units.s**2)[1:,1:])
            plot = figure.add_subplot(1,2,2)
            plot.imshow((phi - self.P0).value_in(units.m**2 / units.s**2)[1:,1:])
            pyplot.show()
            print phi[:,1][0:10]
            print  ( phi[:,1][0:10] -         self.P0[:,1][0:10]) /    phi[:,1][0:10] 
        
    
    def gravity_for_code(self, field, grid):
        x = field.x.flatten()
        y = field.y.flatten()
        z = field.z.flatten()
        rho = field.rho.flatten()
                
        cell_size = (self.length_scale / self.ncells)
        cell_volume = cell_size * cell_size * (1 | units.cm)
        epsilon_squared = cell_size ** 2
        cell_mass = rho * cell_volume
        for cell in grid.iter_cells():
            dx = x - cell.x
            dy = y - cell.y
            dz = z - cell.z
            r_squared = (dx**2 + dy**2 + dz**2 + epsilon_squared)
            r = r_squared.sqrt()
            
            cell.potential = constants.G * (cell_mass / r).sum()
        phi = grid.potential[:,:,0]
        
        figure = pyplot.figure(figsize=(10,5))
        plot = figure.add_subplot(1,1,1)
        plot.imshow(phi.value_in(units.m**2 / units.s**2))
        pyplot.show()
        
class TestJeansInstability(amusetest.TestCase):
    
    def test1(self):
        ncells = 5
        
        x = JeansInstability(
            1.5e7 | units.g / units.cm**3, 
            1.5e7 | dyn / units.cm**2,
            ncells = ncells
        )
        self.assertAlmostRelativeEquals(x.jeans_wavenumber, 2.747 | units.cm**-1, 3)
        x.setup()
        self.assertAlmostRelativeEquals( (x.code.grid.rho[:,0,0].sum()).as_quantity_in(units.g / units.cm**3) / ncells, 1.5e7 | units.g / units.cm**3)
        x.stop()
        
        
if __name__ == '__main__':
    ncells = 1000
        
    run = JeansInstability(
        1.5e7 | units.g / units.cm**3, 
        1.5e7 | dyn / units.cm**2,
        ncells = ncells
    )
    run.setup()
    t = 0 | units.s
    x = [] | units.s
    y = [] | units.g / units.cm**3
    
    x.append(t)
    y.append(run.code.grid[2,2,0].rho)
    
    while t < 10 | units.s:
        t += 0.01 | units.s
        run.evolve_model(t)
        x.append(t)
        y.append(run.code.grid[2,2,0].rho)
        print "evolved to", t, run.code.grid[2,2,0].rho
    
    if IS_PLOT_AVAILABLE:
        figure = pyplot.figure(figsize=(10,5))
        plot = figure.add_subplot(1,1,1)
        plot.plot(x.value_in(units.s), y.value_in(units.g / units.cm**3))
        
        if RUN_WITH_SELF_GRAVITY:
            figure.savefig('jeans_instability_grav.png')
        else:
            figure.savefig('jeans_instability_nograv.png')
        pyplot.show()

        
    
