"""
In this script we simulate Kelvin-Helmholtz Instability in 3d.
"""
import numpy

from amuse.support.core import late
from amuse.units.quantities import VectorQuantity
from amuse.units.generic_unit_system import *

from amuse.community.capreole.interface import Capreole
from amuse import io
from amuse.io import text

from amuse.datamodel import Grid
from amuse.datamodel.grids import SamplePointsOnMultipleGrids
from amuse.datamodel.grids import SamplePointWithIntepolation
from amuse.datamodel.grids import SamplePointOnCellCenter
try:
    from amuse import plot
    from matplotlib import pyplot
    IS_PLOT_AVAILABLE = True
except ImportError:
    IS_PLOT_AVAILABLE = False
    
class CalculateLinearWave1D(object):
    gamma = 5.0/3.0
    wave_flag = 0
    
    def __init__(self, 
        number_of_grid_points =  10, 
        number_of_workers = 1, 
        name_of_the_code = "athena",
        amplitude = 1e-4 | speed,
        vflow_factor = 1.0,
        grid_length = 1.0 | length,
        number_of_steps = 10,
    ):
        
        self.number_of_grid_points = number_of_grid_points
        self.number_of_workers = number_of_workers
        self.name_of_the_code = name_of_the_code
        self.amplitude = amplitude
        self.vflow_factor = vflow_factor
        self.grid_length = 1.0 | length
        self.number_of_steps = number_of_steps
        self.dimensions_of_mesh = (
            self.number_of_grid_points, 
            self.number_of_grid_points, 
            1
        )
        
    def new_instance_of_code(self):
        attribute = "new_instance_of_{0}_code".format(self.name_of_the_code.lower())
        return getattr(self,attribute)()
        
    def new_instance_of_capreole_code(self):
        result=Capreole(number_of_workers=self.number_of_workers)
        result.initialize_code()
        self.dimensions_of_mesh =  (
            self.number_of_grid_points, 
            self.number_of_grid_points, 
            3
        )
        return result
        
    def new_instance_of_athena_code(self):
        from amuse.community.athena.interface import Athena
        result=Athena(number_of_workers=self.number_of_workers)
        result.initialize_code()
        result.parameters.gamma = self.gamma
        result.parameters.courant_number=0.8
        return result
        

    def new_instance_of_mpiamrvac_code(self):
        from amuse.community.mpiamrvac.interface import MpiAmrVac
        result=MpiAmrVac(mode="2d", number_of_workers=self.number_of_workers, debugger="xterm")
        result.set_parameters_filename(result.default_parameters_filename)
        result.initialize_code()
        return result
        
    def set_parameters(self, instance):
        instance.parameters.mesh_size = self.dimensions_of_mesh
        
        instance.parameters.length_x = self.grid_length
        instance.parameters.length_y = self.grid_length
        instance.parameters.length_z = self.grid_length
        
        instance.parameters.x_boundary_conditions = ("periodic","periodic")
        instance.parameters.y_boundary_conditions = ("periodic","periodic")
        instance.parameters.z_boundary_conditions = ("periodic","periodic")
        
        result = instance.commit_parameters()
    
    def new_grid(self):
        grid = Grid.create(self.dimensions_of_mesh, [1,1,1] | length)
        self.clear_grid(grid)
        return grid
    
    def clear_grid(self, grid):
        density = mass / length**3
        momentum =  speed * density
        energy =  mass / (time**2 * length)

        grid.rho =  0.0 | density
        grid.rhovx = 0.0 | momentum
        grid.rhovy = 0.0 | momentum
        grid.rhovz = 0.0 | momentum
        grid.energy = 0.0 | energy
    
        return grid
    
    
    def new_rhoe_right_eigenmatrix(self, velocity, amplitude, enthalpy):
        right_eigenmatrix = numpy.zeros ( (5,5) ) | speed
        
        right_eigenmatrix[0][0] = 1.0 | speed;
        right_eigenmatrix[1][0] = velocity[0] - amplitude;
        right_eigenmatrix[2][0] = velocity[1];
        right_eigenmatrix[3][0] = velocity[2];
        right_eigenmatrix[4][0] = (1.0 | time/length) * (enthalpy - velocity[0]*amplitude);
        #right_eigenmatrix[0][1] = 0.0;
        #right_eigenmatrix[1][1] = 0.0;
        right_eigenmatrix[2][1] = 1.0 | speed;
        #right_eigenmatrix[3][1] = 0.0; 
        right_eigenmatrix[4][1] = velocity[1];

        #right_eigenmatrix[0][2] = 0.0; */
        #right_eigenmatrix[1][2] = 0.0; */
        #right_eigenmatrix[2][2] = 0.0; */
        right_eigenmatrix[3][2] = 1.0 | speed;
        right_eigenmatrix[4][2] = velocity[2];

        right_eigenmatrix[0][3] = 1.0 | speed;
        right_eigenmatrix[1][3] = velocity[0];
        right_eigenmatrix[2][3] = velocity[1];
        right_eigenmatrix[3][3] = velocity[2];
        right_eigenmatrix[4][3] = 0.5*velocity.length();

        right_eigenmatrix[0][4] = 1.0 | speed;
        right_eigenmatrix[1][4] = velocity[0] + amplitude;
        right_eigenmatrix[2][4] = velocity[1];
        right_eigenmatrix[3][4] = velocity[2];
        right_eigenmatrix[4][4] = (1.0 | time/length) *  (enthalpy + velocity[0]*amplitude);
        return right_eigenmatrix

    def initialize_grid(self, grid):
        density = mass / length**3
        momentum =  speed * density
        energy =  mass / (time**2 * length)
        rho =  1.0 | density
        pressure = (1.0/self.gamma) | (mass / (length * time**2))
        vx = (self.gamma * pressure / rho).sqrt()
        velocity = self.vflow_factor * vx * [1.0, 0.0, 0.0]
        velocity_squared = velocity.length()
        energy = (pressure/(self.gamma - 1.0) + (0.5 | length / time )*rho*velocity_squared)
        enthalpy = (energy + pressure)/rho;
        amplitude_squared = (self.gamma - 1.0) * max(enthalpy - (0.5 | length/time)* velocity_squared, 1e-100 | enthalpy.unit) 
        amplitude =amplitude_squared.sqrt()
        
        nwave = 5
        eigenvalues = ([0] * nwave) | speed
        eigenvalues[0] = velocity[0] - amplitude
        eigenvalues[1] = velocity[0]
        eigenvalues[2] = velocity[0]
        eigenvalues[3] = velocity[0]
        eigenvalues[4] = velocity[0] + amplitude
        
        right_eigenmatrix = self.new_rhoe_right_eigenmatrix(velocity, amplitude, enthalpy)
        
        grid.rho = rho
        grid.energy = energy
        grid.rhovx = rho*self.vflow_factor*(1.0 | speed)
        
        wave = self.amplitude*numpy.sin(grid.x * (2.0 | length**-1)*numpy.pi)
        
        grid.rho += wave*right_eigenmatrix[0][self.wave_flag] * (1.0 |mass * time**2 / length**5)
        grid.rhovx += wave*right_eigenmatrix[1][self.wave_flag] * (1.0 |mass * time / length**4)
        grid.rhovy += wave*right_eigenmatrix[2][self.wave_flag] * (1.0 |mass * time / length**4)
        grid.rhovz += wave*right_eigenmatrix[3][self.wave_flag] * (1.0 |mass * time / length**4)
        grid.energy += wave*right_eigenmatrix[4][self.wave_flag] *(1.0 | mass  / length**3)
        
    def store_grids(self, grids, step):
        if __name__ == '__plot__':
            return
        
        grids_in_memory = [x.copy_to_memory() for x in grids]
        io.write_set_to_file(
            grids_in_memory, 
            "linear_wave_{2}_{0}_{1}.vtu".format(self.number_of_grid_points, step, self.name_of_the_code),
            "vtu",
            is_multiple=True
        )
            
    def get_solution_at_time(self, time):
        instance=self.new_instance_of_code()
        
        self.set_parameters(instance)
        
        self.start_grids = []
        
        for x in instance.itergrids():
            inmem = x.copy_to_memory()
            self.clear_grid(inmem)
            self.initialize_grid(inmem)
            self.start_grids.append(inmem)
            from_model_to_code = inmem.new_channel_to(x)
            from_model_to_code.copy()
        
        print "start evolve"
        dt = time / self.number_of_steps
        t = dt
        step = 1
        while t <= time:
            if t == time:
                 instance.parameters.must_evolve_to_exact_time = True
            else:
                instance.parameters.must_evolve_to_exact_time = False
            instance.evolve_model(t)
            
            #print "time : ", t, instance.model_time#,  instance.parameters.must_evolve_to_exact_time 
            
            t += dt
            step += 1
        
        print "copying results"
        result = []
        for x in instance.itergrids():
            result.append(x.copy_to_memory())

        print "terminating code"
        instance.stop()

        return result
    
            
def main():
    number_of_grid_points = 128
    name_of_the_code = 'athena'
    model = CalculateLinearWave1D(
        number_of_grid_points = number_of_grid_points,
        number_of_workers = 1,
        name_of_the_code = name_of_the_code,
        amplitude = 1e-1 | speed,
        vflow_factor = 1.0,
        grid_length = 1.0 | length,
        number_of_steps = 5
    )
    if not IS_PLOT_AVAILABLE:
        return
        
    grids = model.get_solution_at_time(1.0 | time)
    rho0 = model.start_grids[0].rho[...,...,0].value_in(density)
    rho = grids[0].rho[...,...,0].value_in(density)
    drho = rho - rho0
    print drho.sum(), rho[0][0], rho0[0][0], drho[0][0]
    x = grids[0].x[...,...,0].value_in(length)
    y = grids[0].y[...,...,0].value_in(length)
    figure = pyplot.figure(figsize=(10,10))
    plot = figure.add_subplot(1,1,1, projection='3d')
    plot.plot_surface(x, y, drho)
    figure.savefig('kelvin_helmholtz_{0}_{1}.png'.format(name_of_the_code, number_of_grid_points))
    pyplot.show()
    
if __name__ == "__main__":
    main()
