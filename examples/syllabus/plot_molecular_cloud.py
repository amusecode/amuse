"""
  example of molecular cloud evolution with explictly 
  split SPH and grav evolution

  Initial condition is a smooth spherical cloud with random velocities
  as in Bonnell et al. (2003)  
  
"""  

import numpy
  
from matplotlib import pyplot 

from amuse.lab import *
from amuse.ext.molecular_cloud import molecular_cloud
from amuse.ext.evrard_test import body_centered_grid_unit_cube
from amuse.plot import sph_particles_plot, native_plot

def create_molecular_cloud(N, Mcloud, Rcloud, t_end):

    converter = nbody_system.nbody_to_si(Mcloud,Rcloud)
    parts=molecular_cloud(targetN=N,convert_nbody=converter,
            base_grid=body_centered_grid_unit_cube, seed=100).result
#    parts = new_plummer_gas_model(N, convert_nbody=converter)

    sph=Fi(converter)
    sph.gas_particles.add_particle(parts)

    sph.evolve_model(t_end)
    ch = sph.gas_particles.new_channel_to(parts)
    ch.copy()
    sph.stop()
    return parts
  
if __name__ in ("__main__","__plot__"):
    sph_particles = create_molecular_cloud(10000, Mcloud=10000. | units.MSun, Rcloud=10. | units.parsec, t_end=1|units.day)
    native_plot.figure(figsize = (10, 10), dpi = 50)
    sph_particles_plot(sph_particles)
    native_plot.show()
