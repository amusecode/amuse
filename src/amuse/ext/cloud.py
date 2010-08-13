from amuse.support.data import core
from amuse.support.units import generic_unit_system

import numpy



def fill_grid_with_cloud_and_medium(
        grid, 
        center = None,
        radius = None,
        rho_medium = 1.0 | generic_unit_system.mass / generic_unit_system.length**3,
        rho_cloud = 0.1 | generic_unit_system.mass / generic_unit_system.length**3,
        gamma = 5.0 / 3.0,
    ):
    pass



def fill_grid_with_spherical_cloud(
        grid, 
        center = None,
        radius = None,
        rho = 1.0 | generic_unit_system.mass / generic_unit_system.length**3,
        rhovx = 0.0 | generic_unit_system.mass / (generic_unit_system.time * generic_unit_system.length**2),
        rhovy = 0.0 | generic_unit_system.mass / (generic_unit_system.time * generic_unit_system.length**2),
        rhovz = 0.0 | generic_unit_system.mass / (generic_unit_system.time * generic_unit_system.length**2),
        energy = 1.0 | generic_unit_system.mass / (generic_unit_system.time**2 * generic_unit_system.length),
        subgridsize = 4,
    ):
    radii = (grid.position - center).lengths()
    
    if subgridsize <= 1:
        selection = radii <= radius
    else:
        dr = grid.cellsize().length()
        selection = radii < (radius - dr)
        
    grid.rho[selection] = rho
    grid.rhovx[selection] = rhovx
    grid.rhovy[selection] = rhovy
    grid.rhovz[selection] = rhovz
    grid.energy[selection] = energy
    
    if subgridsize <= 1:
        return
    
    selection = numpy.logical_and( radii >= (radius-dr) , radii <= (radius+dr))
    subgrid = core.Grid.create((subgridsize, subgridsize, subgridsize), grid.cellsize())
    subgrid.x -= grid.cellsize()[0] / 2.0
    subgrid.y -= grid.cellsize()[1] / 2.0
    subgrid.z -= grid.cellsize()[2] / 2.0
    x_indices, y_indices, z_indices = grid.indices()
    x_indices = x_indices[selection]
    y_indices = y_indices[selection]
    z_indices = z_indices[selection]
    print len(x_indices)
    
    position = subgrid.position
    centers = center - grid.position[selection]
    
    for i in range(len(x_indices)):
        x_index = x_indices[i]
        y_index = y_indices[i]
        z_index = z_indices[i]
            
        center_of_cloud_for_subgrid = centers[i]
        radii = (position - center_of_cloud_for_subgrid).lengths()
         
        cell = grid[x_index,y_index,z_index]
        
        subgrid.rho = cell.rho
        subgrid.rhovx = cell.rhovx
        subgrid.rhovy = cell.rhovy
        subgrid.rhovz = cell.rhovz
        subgrid.energy = cell.energy
        
        subgrid_selection = radii <= radius
        
        subgrid.rho[subgrid_selection] = rho
        subgrid.rhovx[subgrid_selection] = rhovx
        subgrid.rhovy[subgrid_selection] = rhovy
        subgrid.rhovz[subgrid_selection] = rhovz
        subgrid.energy[subgrid_selection] = energy
        
        cell.rho = subgrid.rho.mean()
        cell.rhovx = subgrid.rhovx.mean()
        cell.rhovy = subgrid.rhovy.mean()
        cell.rhovz = subgrid.rhovz.mean()
        cell.energy = subgrid.energy.mean()


