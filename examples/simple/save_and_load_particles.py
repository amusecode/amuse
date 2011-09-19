"""
example demonstrating the saving and loading of particle sets

AMUSE objects such as quantities with units can be stored using the standard
Python  Pickle module. For particlesets this does not work. Functions are 
available to read and write to a variety of formats.
"""

from amuse.units import nbody_system
from amuse.io import write_set_to_file

from amuse.ic.plummer import new_plummer_sphere
if __name__ in ('__main__', '__plot__'):

# generate a particle set
  plummer=new_plummer_sphere(128)

#write the set to file 'testfile'
# the third argument is the file format, 'amuse' is an hdf5 based format
# that saves all information. Other formats
# available are e.g. csv, txt, gadget, starlab
  write_set_to_file(plummer,'plummer128','amuse')
  del plummer

# reading back the file
  from amuse.io import read_set_from_file
  particles=read_set_from_file('plummer128','amuse')

#plotting
  from amuse.plot import *
  plot(particles.x,particles.y,'r.')
  native_plot.xlim(-5,5)
  native_plot.ylim(-5,5)
  native_plot.show()

#running a simulation
  from amuse.community.phiGRAPE.interface import PhiGRAPE
  phi=PhiGRAPE()
  phi.particles.add_particles(particles)

  phi.evolve_model( 1| nbody_system.time)

#saving the result
  write_set_to_file(phi.particles,'evolvedplummer','amuse')

  phi.stop()
