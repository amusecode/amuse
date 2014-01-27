"""
In this script we simulate a 1d linear wave in a 2d field, 
the periodic boundaries are not handled in the code but by amuse
(much slower at time of writing)
"""
import numpy
import math

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
    from matplotlib import animation
    IS_PLOT_AVAILABLE = True
except ImportError:
    IS_PLOT_AVAILABLE = False

USE_BOUNDARIES = True

class EvolveHydrodynamicsCodeWithAmusePeriodicBoundaries(object):
    
    def __init__(self, code, number_of_grid_points, nghost):
        self.code = code
        self.number_of_grid_points = number_of_grid_points
        self.nghost = nghost
        
    def set_parameters(self):
        self.code.parameters.x_boundary_conditions = ("interface","interface")
        self.code.parameters.y_boundary_conditions = ("interface","interface")
        self.code.stopping_conditions.number_of_steps_detection.enable()
        
    def init_channels(self):
        instance = self.code
        self.xbound1 = instance.get_boundary_grid('xbound1')
        self.xbound2 = instance.get_boundary_grid('xbound2')
        self.ybound1 = instance.get_boundary_grid('ybound1')
        self.ybound2 = instance.get_boundary_grid('ybound2')
    
        self.xbound1_channel = instance.grid[self.number_of_grid_points - self.nghost:,..., ...].new_channel_to(self.xbound1)
        self.xbound2_channel = instance.grid[0: self.nghost,..., ...].new_channel_to(self.xbound2)
        self.ybound1_channel = instance.grid[...,self.number_of_grid_points-self.nghost:, ...].new_channel_to(self.ybound1[self.nghost:self.number_of_grid_points+self.nghost,...,...])
        self.ybound2_channel = instance.grid[...,0:self.nghost, ...].new_channel_to(self.ybound2[self.nghost:self.number_of_grid_points+self.nghost,...,...])
        
        xbound1_bottom = self.xbound1[0:self.nghost,0:self.nghost,...]
        xbound1_top = self.xbound1[0:self.nghost,self.number_of_grid_points - self.nghost:,...]
        ybound1_left = self.ybound1[0:self.nghost,0:self.nghost,...]
        ybound2_left = self.ybound2[0:self.nghost,0:self.nghost,...]
        
        self.xbound1_ybound1_channel = xbound1_top.new_channel_to(ybound1_left)
        self.xbound1_ybound2_channel = xbound1_bottom.new_channel_to(ybound2_left)
         
        xbound2_bottom = self.xbound2[0:self.nghost,0:self.nghost,...]
        xbound2_top = self.xbound2[0:self.nghost,self.number_of_grid_points - self.nghost:,...]
        ybound1_right = self.ybound1[self.number_of_grid_points+self.nghost:,0:self.nghost,...]
        ybound2_right= self.ybound2[self.number_of_grid_points+self.nghost:,0:self.nghost,...]
        
        self.xbound2_ybound1_channel = xbound2_top.new_channel_to(ybound1_right)
        self.xbound2_ybound2_channel = xbound2_bottom.new_channel_to(ybound2_right)
        
        self.copy_to_boundary_cells()
    
    def copy_to_boundary_cells(self):
        self.xbound1_channel.copy()
        self.xbound2_channel.copy()
        self.ybound1_channel.copy()
        self.ybound2_channel.copy()
        
        self.xbound1_ybound1_channel.copy()
        self.xbound1_ybound2_channel.copy()
        self.xbound2_ybound1_channel.copy()
        self.xbound2_ybound2_channel.copy()
        
    def evolve_model(self, time, endtime):
        while self.code.model_time < time:
                
            if (self.code.get_timestep() + self.code.model_time) >= time and time == endtime:
                self.code.parameters.must_evolve_to_exact_time = True
            
            self.code.evolve_model(time)
            self.copy_to_boundary_cells()

class EvolveHydrodynamicsCodeWithAmusePeriodicBoundariesAndNCodes(object):
    
    def __init__(self, codes, number_of_grid_points, number_of_grid_points_per_code, nghost):
        self.codes = codes
        self.number_of_grid_points = number_of_grid_points
        self.number_of_grid_points_per_code = number_of_grid_points_per_code
        self.nghost = nghost
        self.min_timestep = None
        
    def set_parameters(self):
        for code in self.codes:
            code.parameters.x_boundary_conditions = ("interface","interface")
            code.parameters.y_boundary_conditions = ("interface","interface")
            code.stopping_conditions.number_of_steps_detection.enable()
        
    def init_channels(self):
        channels = []
        after_channels = []
        for index in range(0, len(self.codes)):
            instance = self.codes[index]
            nx,ny,nz = instance.grid.shape
            print nx,ny,nz
            xbound1 = instance.get_boundary_grid('xbound1')
            xbound2 = instance.get_boundary_grid('xbound2')
            xbound1_nghost_x,_,ybound1_nghost_z = xbound1.shape
            xbound2_nghost_x,_,ybound2_nghost_z = xbound2.shape
            print xbound1_nghost_x,xbound2_nghost_x
            ybound1 = instance.get_boundary_grid('ybound1')
            ybound2 = instance.get_boundary_grid('ybound2')
            ybound1_nghost_x,ybound1_nghost_y,ybound1_nghost_z = ybound1.shape
            ybound2_nghost_x,ybound2_nghost_y,ybound2_nghost_z = ybound2.shape
            ybound1_nghost_x = (ybound1_nghost_x - nx) / 2
            ybound2_nghost_x = (ybound2_nghost_x - nx) / 2
            
            
            xbound1_channel = instance.grid[nx - xbound1_nghost_x:,..., ...].new_channel_to(xbound1)
            xbound2_channel = instance.grid[0:xbound2_nghost_x,..., ...].new_channel_to(xbound2)
            if index == 0:
                instance11 = self.codes[-1]
            else:
                instance11 = self.codes[index-1]
            print ybound1[ybound1_nghost_x:nx+ybound1_nghost_x,...,...].shape
            print instance11.grid[...,ny-ybound1_nghost_y:, ...].shape
            ybound1_channel = instance11.grid[...,ny-ybound1_nghost_y:, ...].new_channel_to(ybound1[ybound1_nghost_x:nx+ybound1_nghost_x,...,...])
            nx11,ny11,nz11 = instance11.grid.shape
            ybound1_left_channel = instance11.grid[0:ybound1_nghost_x,ny11-ybound1_nghost_y:, ...].new_channel_to(ybound1[0:ybound1_nghost_x,...,...])
            ybound1_right_channel = instance11.grid[nx11-ybound1_nghost_x:,ny11-ybound1_nghost_y:, ...].new_channel_to(ybound1[nx+ybound1_nghost_x:,...,...])
            
            if index == len(self.codes)-1:
                instance12 = self.codes[0]
            else:
                instance12 = self.codes[index+1]
            ybound2_channel = instance12.grid[...,0:ybound2_nghost_y, ...].new_channel_to(ybound2[ybound2_nghost_x:nx+ybound2_nghost_x,...,...])
            nx12,ny12,nz12 = instance12.grid.shape
        
            ybound2_left_channel = instance12.grid[0:ybound2_nghost_x,0:ybound2_nghost_y, ...].new_channel_to(ybound2[0:ybound2_nghost_x,...,...])
            ybound2_right_channel = instance12.grid[nx12-ybound2_nghost_x:,0:ybound2_nghost_y, ...].new_channel_to(ybound2[nx+ybound2_nghost_x:,...,...])
            
            channels.append(xbound1_channel)
            channels.append(xbound2_channel)
            channels.append(ybound1_channel)
            channels.append(ybound2_channel)
            after_channels.append(ybound1_left_channel)
            after_channels.append(ybound1_right_channel)
            after_channels.append(ybound2_left_channel)
            after_channels.append(ybound2_right_channel)

        self.channels = channels
        self.channels.extend(after_channels)
        self.copy_to_boundary_cells()
        
    def copy_to_boundary_cells(self):
        for channel in self.channels:
            channel.copy()
        
    def evolve_model(self, time, endtime):
        
        code = self.codes[0]
        t_unit = code.get_timestep().unit
        while code.model_time < time:
            timesteps = [] | t_unit
            for x in self.codes:
                timesteps.append(x.get_timestep())
            min_timestep = timesteps.min()
            if code.model_time == 0.0 | t_unit:
                min_timestep = time
                
            if (min_timestep + code.model_time) >= time and time == endtime:
                for x in self.codes:
                    x.parameters.must_evolve_to_exact_time = True
            print min_timestep
            for x in self.codes:
                x.set_timestep(min_timestep)
                
            for x in self.codes:
                x.evolve_model(x.model_time + (min_timestep *2))
                print "MODEL_TIME:", x.model_time
                
            self.copy_to_boundary_cells()
            
class EvolveHydrodynamicsCodeWithAmusePeriodicBoundariesAndNCodesWithDifferentGridsSizes(object):
    
    """
        first version, grids connect on the y axis, whole system is periodic
    """
    def __init__(self, codes, number_of_grid_points, number_of_grid_points_per_code, nghost):
        self.codes = codes
        self.number_of_grid_points = number_of_grid_points
        self.number_of_grid_points_per_code = number_of_grid_points_per_code
        self.nghost = nghost
        self.min_timestep = None
        
    def set_parameters(self):
        for code in self.codes:
            code.parameters.x_boundary_conditions = ("interface","interface")
            code.parameters.y_boundary_conditions = ("interface","interface")
            code.stopping_conditions.number_of_steps_detection.enable()
        
    def init_channels(self):
        channels = []
        after_channels = []
        for index in range(0, len(self.codes)):
            instance = self.codes[index]
            xbound1 = instance.get_boundary_grid('xbound1')
            xbound2 = instance.get_boundary_grid('xbound2')
            
            ybound1 = instance.get_boundary_grid('ybound1')
            ybound2 = instance.get_boundary_grid('ybound2')
            
            xbound1_bottom = xbound1[0:self.nghost,0:self.nghost,...]
            xbound1_top = xbound1[0:self.nghost,self.number_of_grid_points - self.nghost:,...]
            
            xbound1_channel = instance.grid[self.number_of_grid_points - self.nghost:,..., ...].new_channel_to(xbound1)
            xbound2_channel = instance.grid[0: self.nghost,..., ...].new_channel_to(xbound2)

            if index == 0:
                instance11 = self.codes[-1]
            else:
                instance11 = self.codes[index-1]
            ybound1_channel = instance11.grid[...,self.number_of_grid_points_per_code-self.nghost:, ...].new_channel_to(ybound1[self.nghost:self.number_of_grid_points+self.nghost,...,...])
            ybound1_left_channel = instance11.grid[0:self.nghost,self.number_of_grid_points_per_code-self.nghost:, ...].new_channel_to(ybound1[0:self.nghost,...,...])
            ybound1_right_channel = instance11.grid[-self.nghost:,self.number_of_grid_points_per_code-self.nghost:, ...].new_channel_to(ybound1[self.number_of_grid_points+self.nghost:,...,...])
            
            if index == len(self.codes)-1:
                instance12 = self.codes[0]
            else:
                instance12 = self.codes[index+1]
            ybound2_channel = instance12.grid[...,0:self.nghost, ...].new_channel_to(ybound2[self.nghost:self.number_of_grid_points+self.nghost,...,...])
        
            xbound1next = instance12.get_boundary_grid('xbound1')
            xbound2next = instance12.get_boundary_grid('xbound2')
            ybound2_left_channel = xbound1next[...,0:self.nghost, ...].new_channel_to(ybound2[0:self.nghost,...,...])
            ybound2_right_channel = xbound2next[...,0:self.nghost, ...].new_channel_to(ybound2[self.number_of_grid_points+self.nghost:,...,...])
            
            channels.append(xbound1_channel)
            channels.append(xbound2_channel)
            channels.append(ybound1_channel)
            channels.append(ybound2_channel)
            after_channels.append(ybound1_left_channel)
            after_channels.append(ybound1_right_channel)
            after_channels.append(ybound2_left_channel)
            after_channels.append(ybound2_right_channel)

        self.channels = channels
        self.channels.extend(after_channels)
        self.copy_to_boundary_cells()
        
    def copy_to_boundary_cells(self):
        for channel in self.channels:
            channel.copy()
        
    def evolve_model(self, time, endtime):
        
        code = self.codes[0]
        t_unit = code.get_timestep().unit
        while code.model_time < time:
            timesteps = [] | t_unit
            for x in self.codes:
                timesteps.append(x.get_timestep())
            min_timestep = timesteps.min()
            if code.model_time == 0.0 | t_unit:
                min_timestep = time
                
            if (min_timestep + code.model_time) >= time and time == endtime:
                for x in self.codes:
                    x.parameters.must_evolve_to_exact_time = True
            
            for x in self.codes:
                x.set_timestep(min_timestep)
                
            for x in self.codes:
                x.evolve_model(time)
                
            self.copy_to_boundary_cells()



class EvolveHydrodynamicsCodeWithPeriodicBoundaries(object):
    def __init__(self, code):
        self.code = code
    
    def set_parameters(self):
        self.code.parameters.x_boundary_conditions = ("periodic","periodic")
        self.code.parameters.y_boundary_conditions = ("periodic","periodic")
        self.code.stopping_conditions.number_of_steps_detection.disable()
        
    def init_channels(self):
        pass 
    
    def evolve_model(self, time, endtime):
        while self.code.model_time < time:
                
            if (self.code.get_timestep() + self.code.model_time) >= time and time == endtime:
                self.code.parameters.must_evolve_to_exact_time = True
            
            self.code.evolve_model(time)


class CalculateLinearWave1D(object):
    gamma = 5.0/3.0
    wave_flag = 0
    
    def __init__(self, 
        number_of_grid_points =  10, 
        number_of_workers = 1, 
        name_of_the_code = "athena",
        amplitude = 1e-6 | speed,
        vflow_factor = 1.0,
        grid_length = 1.0 | length,
        number_of_steps = 10,
        use_boundaries = True,
        number_of_codes = 2
    ):
        
        self.number_of_grid_points = number_of_grid_points
        self.number_of_workers = number_of_workers
        self.name_of_the_code = name_of_the_code
        self.amplitude = amplitude
        self.vflow_factor = vflow_factor
        self.grid_length = 1.0 | length
        self.number_of_steps = number_of_steps
        self.number_of_codes = number_of_codes
        self.dimensions_of_mesh = (
            self.number_of_grid_points, 
            self.number_of_grid_points, 
            4
        )
        self.nghost = 4
        self.use_boundaries = use_boundaries
        
    def new_instance_of_code(self):
        if 1:
            if self.name_of_the_code == 'athena':
                self.name_of_the_code = 'capreole'
            else:
                self.name_of_the_code = 'athena'
            
        attribute = "new_instance_of_{0}_code".format(self.name_of_the_code.lower())
        return getattr(self,attribute)()
        
    def new_instance_of_capreole_code(self):
        result=Capreole(number_of_workers=self.number_of_workers)
        self.dimensions_of_mesh =  (
            self.number_of_grid_points, 
            self.number_of_grid_points, 
            4
        )
        self.nghost = 2
        return result
        
    def new_instance_of_athena_code(self):
        from amuse.community.athena.interface import Athena
        result=Athena(number_of_workers=self.number_of_workers)
        result.parameters.gamma = self.gamma
        result.parameters.courant_number=0.4
        self.nghost = 4
        return result
        

    def new_instance_of_mpiamrvac_code(self):
        raise Exception("MPIAMRVAC does not yet have support for detailed boundaries in amuse")
        from amuse.community.mpiamrvac.interface import MpiAmrVac
        result=MpiAmrVac(mode="2d", number_of_workers=self.number_of_workers, debugger="xterm")
        result.set_parameters_filename(result.default_parameters_filename)
        result.initialize_code()
        return result
        
    def set_parameters(self, codes, evolve):
        
        mesh_for_code = list(self.dimensions_of_mesh)
        mesh_for_code[1] /= self.number_of_codes
        for instance in codes:
            instance.parameters.mesh_size = list(mesh_for_code)
            
            instance.parameters.length_x = self.grid_length
            instance.parameters.length_y = self.grid_length / self.number_of_codes
            instance.parameters.length_z = self.grid_length
            
            instance.parameters.z_boundary_conditions = ("periodic","periodic")
        evolve.set_parameters()
        
        for instance in codes:
            instance.commit_parameters()
    
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
        grid.rhovy = rho*self.vflow_factor*(1.0 | speed)
        
        wave = self.amplitude*numpy.sin(grid.y * (2.0 | length**-1)*numpy.pi)
        
        grid.rho += wave*right_eigenmatrix[0][self.wave_flag] * (1.0 |mass * time**2 / length**5)
        grid.rhovx += wave*right_eigenmatrix[3][self.wave_flag] * (1.0 |mass * time / length**4)
        grid.rhovy += wave*right_eigenmatrix[1][self.wave_flag] * (1.0 |mass * time / length**4)
        grid.rhovz += wave*right_eigenmatrix[2][self.wave_flag] * (1.0 |mass * time / length**4)
        grid.energy += wave*right_eigenmatrix[4][self.wave_flag] *(1.0 | mass  / length**3)
        
    def store_grids(self, grids, step):
        if __name__ == '__plot__':
            return
        
        grids_in_memory = [x.copy() for x in grids]
        io.write_set_to_file(
            grids_in_memory, 
            "linear_wave_{2}_{0}_{1}.vtu".format(self.number_of_grid_points, step, self.name_of_the_code),
            "vtu",
            is_multiple=True
        )
    
    def new_evolve_object(self, instance):
        """Returns a special object to evolve to code in time"""
        if self.use_boundaries:
            #return EvolveHydrodynamicsCodeWithAmusePeriodicBoundaries(instance[0], self.number_of_grid_points, self.nghost)
            return EvolveHydrodynamicsCodeWithAmusePeriodicBoundariesAndNCodes(instance, self.number_of_grid_points, self.number_of_grid_points / self.number_of_codes , self.nghost)
        else:
            return EvolveHydrodynamicsCodeWithPeriodicBoundaries(instance[0])
            
    def initialize(self):
        self.codes = []
        for i in range(self.number_of_codes):
            self.codes.append(self.new_instance_of_code())
        
        self.evolve = self.new_evolve_object(self.codes)
        
        self.set_parameters(self.codes, self.evolve)
        
        self.start_grids = []
        
        offset = 0.0 * self.grid_length
        for code in self.codes:
            for x in code.itergrids():
                inmem = x.copy()
                self.clear_grid(inmem)
                inmem.y += offset
                self.initialize_grid(inmem)
                self.start_grids.append(inmem)
                from_model_to_code = inmem.new_channel_to(x)
                from_model_to_code.copy()
            offset += self.grid_length / self.number_of_codes
            
            
        for x in self.codes:
            x.initialize_grid()
        self.evolve.init_channels()
        
    def evolve_model(self,time, endtime):
        for code in self.codes:
            code.parameters.must_evolve_to_exact_time = False
        self.evolve.evolve_model(time, endtime)
            
    
    def get_grids(self):
        result = []
        
        offset = 0.0 * self.grid_length
        for code in self.codes:
            for x in code.itergrids():
                inmem = x.copy()
                inmem.y += offset
                result.append(inmem)
            offset += self.grid_length / self.number_of_codes

        return result
    
    def stop(self):
        print "terminating code"
        
        for code in self.codes:
            code.stop()

import sys
def main():
    number_of_grid_points = 60
    name_of_the_code = 'athena'
    number_of_steps = 2000
    vflow_factor = -1.0
    pertubation_amplitude = 1e-4 | speed
    grid_length = 1.0 | length
    number_of_codes = int(sys.argv[1])
    if number_of_grid_points % number_of_codes != 0:
        raise Exception("grid points should be dividable by the number of codes")
    model1 = CalculateLinearWave1D(
        number_of_grid_points = number_of_grid_points,
        number_of_workers = 1,
        name_of_the_code = name_of_the_code,
        amplitude = pertubation_amplitude,
        vflow_factor = vflow_factor,
        grid_length = grid_length,
        number_of_steps = number_of_steps,
        use_boundaries = True,
        number_of_codes = number_of_codes
    )
    if 0:
        model2 = CalculateLinearWave1D(
            number_of_grid_points = number_of_grid_points,
            number_of_workers = 1,
            name_of_the_code = name_of_the_code,
            amplitude = pertubation_amplitude,
            vflow_factor = vflow_factor,
            grid_length = grid_length,
            number_of_steps = number_of_steps,
            use_boundaries = False,
            number_of_codes = number_of_codes
        )
        
    if not IS_PLOT_AVAILABLE:
        return
    model1.initialize()
    #model2.initialize()
    
    
    grids1 = model1.get_grids()
    #grids2 = model2.get_grids()

    
    figure = pyplot.figure(figsize=(10,5))
    plot1 = figure.add_subplot(1,1,1)
    lines = []
    ys = []
    for grid in grids1:
        y = grid.y[0,...,0].value_in(length)
        ys.append(y)
        rho = grid.rho[0,...,0].value_in(density)
        line = plot1.plot(y,rho)[0]
        lines.append(line)
    
    end_time = 10.0 | time
    dt = end_time / number_of_steps
    
    t = dt
    step = 1
    
    title = figure.suptitle('{0:.3f}'.format(0.0))
    variables = [t, step]

    def update(xx):
        t, step = variables
        title.set_text('{0:.3f}'.format(t.value_in(time)))
        model1.evolve_model(t, end_time)
        #model2.evolve_model(t, end_time)
        t += dt
        grids1 = model1.get_grids()
        #grids2 = model2.get_grids()
        
        for line, grid in zip(lines, grids1):
            y = grid.y[0,...,0].value_in(length)
            rho = grid.rho[0,...,0].value_in(density)
            line.set_data(y,rho)
            #line = plot1.plot(y,rho)[0]
            #lines.append(line)
        print t
        step += 1
        variables[0] = t
        variables[1] = step
        return lines
    if 1:
        process = animation.FuncAnimation(
            figure, 
            update, 
            numpy.arange(1, 200), 
            interval=25,
            blit=False
        )
    else:
        update(0)
        update(0)
        pass
    
        
    pyplot.show()
    model1.stop()
    #model2.stop()
    
if __name__ == "__main__":
    main()