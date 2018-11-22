import numpy

from amuse.units import units
from amuse.support.exceptions import AmuseException
from amuse.community.gadget2.interface import Gadget2
from amuse.community.fi.interface import Fi
from amuse.datamodel import Grid

class _SPH2Grid(object):
    
    def __init__(self, sph_code, dimensions, do_scale = False):
        if (sph_code.mode != sph_code.MODE_PERIODIC_BOUNDARIES):
            raise AmuseException("Only periodic boundary conditions supported")
        
        if len(dimensions) != 3:
            raise AmuseException("Argument dimensions must contain exactly three numbers")
        
        if isinstance(sph_code, Fi):
            self.box_offset = sph_code.parameters.periodic_box_size / 2.0
        elif isinstance(sph_code, Gadget2):
            self.box_offset = sph_code.parameters.periodic_box_size * 0.0
        else:
            raise AmuseException("Unknown hydrodynamics code: {0} - don't know whether the "
            "box runs from 0 to L or from -0.5 L to 0.5 L.".format(sph_code.__class__.__name__))
        
        self.sph_code = sph_code
        self.dimensions = dimensions
        self.do_scale = do_scale
        self.box_size = sph_code.parameters.periodic_box_size
    
    @property
    def result(self):
        
        grid = Grid.create(self.dimensions, 
            self.box_size.as_vector_with_length(3))
        grid.add_vector_attribute("momentum", ["rhovx","rhovy","rhovz"])

        
        zero = numpy.zeros(self.dimensions).flatten() | units.m / units.s
        rho, rhovx, rhovy, rhovz, rhoe = self.sph_code.get_hydro_state_at_point(
            grid.x.flatten() - self.box_offset, 
            grid.y.flatten() - self.box_offset, 
            grid.z.flatten() - self.box_offset, zero, zero, zero)
        grid.rho = rho.reshape(self.dimensions)
        grid.rhovx = rhovx.reshape(self.dimensions)
        grid.rhovy = rhovy.reshape(self.dimensions)
        grid.rhovz = rhovz.reshape(self.dimensions)
        grid.energy = rhoe.reshape(self.dimensions)
        
        if self.do_scale:
            grid.rho *= self.sph_code.gas_particles.total_mass() / (grid.rho.sum() * grid.cellsize().prod())
            
            total_sph_momentum = (self.sph_code.gas_particles.mass.reshape((-1,1)) * self.sph_code.gas_particles.velocity).sum(axis=0)
            total_grid_momentum = grid.momentum.reshape((-1,3)).sum(axis=0)
            grid.momentum += total_sph_momentum / grid.cellsize().prod() - total_grid_momentum
            
            grid.energy *= self.sph_code.thermal_energy / (grid.energy.sum() * grid.cellsize().prod())
        
        return grid
    


def convert_SPH_to_grid(sph_code, dimensions, **keyword_arguments):
    """
    Currently works for periodic boundary conditions only...
    
    Converts an SPH realization into a (cartesian) hydrodynamics Grid.
    The SPH realization must reside in the SPH code 'sph_code', supporting the 
    get_hydro_state_at_point function. The number of grid cells in the x, y, 
    and z direction are determined by the 'dimensions' argument.
    
    :argument sph_code:  SPH code in which the gas particles reside
    :argument dimensions:  Tuple of three integers, defining the meshsize
    :argument do_scale:  If True, the grid density, momentum, and energy are 
        scaled to conserve the mass, momentum, and energy of the original model
    """
    converter = _SPH2Grid(sph_code, dimensions, **keyword_arguments)
    return converter.result
