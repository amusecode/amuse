import numpy
from amuse.units import units
from amuse.support.data.core import Particles
from amuse.support.exceptions import AmuseException

class Grid2SPH(object):
    """
    Converts a (cartesian) hydrodynamics Grid into an SPH model consisting of the 
    specified number of particles. The Grid must have position, rho, momentum and
    energy defined on each cell.
    
    grid, number_of_sph_particles, base_distribution_type = "uniform", seed = None):
    :argument grid:  Star particle to be converted to an SPH model
    :argument number_of_sph_particles:  Number of gas particles in the resulting model
    :argument base_distribution_type:  Type of the base particle distribution ("random" or "uniform")
        "random":  particle positions are randomly sampled from the density weighted cells
        "uniform":  particles are uniformly sampled from the density weighted cells (note
            that they are assigned a random (!) position within that cell subsequently)
    :argument seed:  If provided, seed for the random number generator
    """
    
    def __init__(self, grid, number_of_sph_particles, base_distribution_type = "uniform", seed = None):
        if (grid.number_of_dimensions() != 3):
            raise AmuseException("Grid must be 3D")
        
        self.grid = grid
        self.shape = grid.shape
        self.number_of_sph_particles = number_of_sph_particles
        self.base_distribution_type = base_distribution_type # "random" or "uniform"
        if seed:
            numpy.random.seed(seed)
    
    def setup_lookup_tables(self):
        # Retrieve details of the grid and convert them to fast lookup tables
        shape_for_vector_multiply = list(self.grid.shape)
        shape_for_vector_multiply.append(1)
        density = self.grid.rho
        summed_density = density.sum()
        self.cumulative_weight = numpy.cumsum((density / summed_density).value_in(units.none))
        self.position_lookup_table = self.grid.position.reshape((-1,3))
        self.velocity_lookup_table = (self.grid.momentum / density.reshape(shape_for_vector_multiply)).reshape((-1,3))
        self.specific_internal_energy_lookup_table = (self.grid.energy / density).flatten()
        cellsize = self.grid.cellsize()
        self.cellsize_unit = cellsize.unit
        self.cellsize_number = cellsize.value_in(cellsize.unit)
        self.mass = summed_density * (cellsize[0] * cellsize[1] * cellsize[2])
    
    def setup_variates(self):
        # Generate (quasi-)random realisation
        variates = self.generate_variates(self.number_of_sph_particles)
        self.indices = numpy.searchsorted(self.cumulative_weight, variates)
    
    def generate_variates(self, number_of_variates):
        if self.base_distribution_type == "uniform":
            return numpy.linspace(0.0, 1.0, num=number_of_variates, endpoint=False)
        elif self.base_distribution_type == "random":
            return numpy.random.uniform(0.0, 1.0, number_of_variates)
        else:
            raise AmuseException("Unknown base_distribution_type: {0}. Possible "
                "options are: 'random' or 'uniform'.".format(self.base_distribution_type))
    
    def new_particle_positions(self):
        base_positions = self.position_lookup_table[self.indices]
        return base_positions + self.cellsize_unit.new_quantity(
            self.cellsize_number * numpy.random.uniform(-0.5, 0.5, (self.number_of_sph_particles, 3)))
    
    def new_particle_velocities(self):
        return self.velocity_lookup_table[self.indices]
    
    def new_particle_specific_internal_energies(self):
        return self.specific_internal_energy_lookup_table[self.indices]
    
    @property
    def result(self):
        self.setup_lookup_tables()
        self.setup_variates()
        
        sph_particles = Particles(self.number_of_sph_particles)
        sph_particles.position = self.new_particle_positions()
        sph_particles.velocity = self.new_particle_velocities()
        sph_particles.u        = self.new_particle_specific_internal_energies()
        
        sph_particles.mass = (self.mass.number * 1.0 / self.number_of_sph_particles) | self.mass.unit
        # Crude estimate of the smoothing length; the SPH code will calculate the true value itself.
        sph_particles.h_smooth = (self.grid.get_volume() * 50.0/self.number_of_sph_particles)**(1/3.0)
        
        return sph_particles

    

def convert_grid_to_SPH(grid, number_of_sph_particles, **keyword_arguments):
    """
    Converts a (cartesian) hydrodynamics Grid into an SPH model consisting of the 
    specified number of particles. The Grid must have position, rho, momentum and
    energy defined on each cell.
    
    grid, number_of_sph_particles, base_distribution_type = "uniform", seed = None):
    :argument grid:  Star particle to be converted to an SPH model
    :argument number_of_sph_particles:  Number of gas particles in the resulting model
    :argument base_distribution_type:  Type of the base particle distribution ("random" or "uniform")
        "random":  particle positions are randomly sampled from the density weighted cells
        "uniform":  particles are uniformly sampled from the density weighted cells (note
            that they are assigned a random (!) position within that cell subsequently)
    :argument seed:  If provided, seed for the random number generator
    """
    converter = Grid2SPH(grid, number_of_sph_particles, **keyword_arguments)
    return converter.result
