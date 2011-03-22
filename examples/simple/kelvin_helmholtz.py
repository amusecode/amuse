"""
In this script we simulate Kelvin-Helmholtz Instability in 3d.
"""
import numpy

from amuse.support.core import late
from amuse.support.data.values import VectorQuantity
from amuse.support.data.core import Grid
from amuse.support import io
from amuse.support.io import text
from amuse.support.units.generic_unit_system import *
from amuse.support.data.grids import SamplePointsOnMultipleGrids, SamplePointWithIntepolation, SamplePointOnCellCenter


from amuse.community.capreole.interface import Capreole
from matplotlib import pyplot
    
class CalculateKelvinHelmholtzInstabilityIn3D(object):
    number_of_workers = 1
    number_of_grid_points = 10
    gamma = 1.4 # 5.0/3.0
    name_of_the_code = "capreole"
    
    def __init__(self, number_of_grid_points =  10, number_of_workers = 1, name_of_the_code = "capreole"):
        
        self.number_of_grid_points = number_of_grid_points
        self.number_of_workers = number_of_workers
        self.name_of_the_code = name_of_the_code
        
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
        result=MpiAmrVac(mode="2d", number_of_workers=self.number_of_workers, redirection="none", debugger="xterm")
        result.set_parameters_filename(result.default_parameters_filename)
        result.initialize_code()
        return result
        
    def set_parameters(self, instance):
        instance.parameters.mesh_size = self.dimensions_of_mesh
        
        instance.parameters.length_x = 1 | length
        instance.parameters.length_y = 1 | length
        instance.parameters.length_z = 1 | length
        
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
    
    def initialize_grid(self, grid):        
        vx = 0.5 | speed
        p = 2.5 | (mass / (length * time**2))
        
        halfway = self.dimensions_of_mesh[0]/2 - 1
        
        outerregion = numpy.logical_or(grid.y <= 0.25 | length, grid.y >= 0.75 | length)
        innerregion = numpy.logical_and(grid.y > 0.25 | length, grid.y < 0.75 | length)
        
        grid[outerregion].rho = 1  | density
        grid[outerregion].rhovx =  vx * grid[outerregion].rho
        
        grid[innerregion].rho = 2.0  | density
        grid[innerregion].rhovx = -vx * grid[innerregion].rho
        
        grid.energy = p / (self.gamma - 1)
        
    def pertubate_grid(self, grid):
        amplitude = 0.01 | speed
        grid.rhovx += grid.rho * amplitude * (numpy.random.rand(*grid.shape) - 0.5)
        grid.rhovy += grid.rho * amplitude * (numpy.random.rand(*grid.shape) - 0.5)
        
        grid.energy += 0.5 * (grid.rhovx ** 2  + grid.rhovy ** 2 + grid.rhovz ** 2) / grid.rho
        
    def store_grids(self, grids, step):
        if __name__ == '__plot__':
            return
            
        for i, x in enumerate(grids):
            mem = x.copy_to_memory()
            io.write_set_to_file(mem, "kelvin_helmholtz_{0}_{2}_{1}.vts".format(self.number_of_grid_points, step, i),"vts")
            
    def get_solution_at_time(self, time):
        instance=self.new_instance_of_code()
        
        self.set_parameters(instance)
        
        
        for x in instance.itergrids():
            inmem = x.copy_to_memory()
            self.clear_grid(inmem)
            self.initialize_grid(inmem)
            self.pertubate_grid(inmem)
            from_model_to_code = inmem.new_channel_to(x)
            from_model_to_code.copy()
        
        instance.initialize_grid()
        
        print "start evolve"
        dt = time / 10.0
        t = dt
        step = 0
        while t < time:
            instance.evolve(t)
            
            print "time : ", t
            
            self.store_grids(instance.itergrids(), step)
                
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
    number_of_grid_points = 400
    model = CalculateKelvinHelmholtzInstabilityIn3D(
        number_of_grid_points = number_of_grid_points,
        number_of_workers = 3,
        name_of_the_code = 'athena'
    )
    grids = model.get_solution_at_time(1.0 | time)
    
    rho = grids[0].rho[...,...,0].value_in(density)
    figure = pyplot.figure(figsize=(20,20))
    plot = figure.add_subplot(1,1,1)
    plot.imshow(rho, origin = 'lower')
    figure.savefig('kelvin_helmholtz.png')
    pyplot.show()
    
if __name__ in ["__main__", "__plot__"]:
    main()
