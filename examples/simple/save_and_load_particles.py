"""
example demonstrating the saving and loading of particle sets

AMUSE objects such as quantities with units can be stored using the standard
Python  Pickle module. For particlesets this does not work. Functions are 
available to read and write to a variety of formats.
"""

from amuse.units import nbody_system
from amuse.ic.plummer import new_plummer_model
from amuse.plot import plot, native_plot

# the imports needed to read and write a set
from amuse.io import write_set_to_file
from amuse.io import read_set_from_file

import os.path

if __name__ in ('__main__', '__plot__'):
    if os.path.exists('plummer128.hdf'):
        os.remove('plummer128.hdf')
        
    # generate a particle set
    plummer=new_plummer_model(128)

    # write the set to file 'testfile'
    # the third argument is the file format, 'amuse' is an hdf5 based format
    # that saves all information. Other formats
    # available are e.g. csv, txt, gadget, starlab
    write_set_to_file(plummer,'plummer128.hdf','amuse')
    del plummer

    # reading back the file
    particles=read_set_from_file('plummer128.hdf','amuse')

    # plotting
    plot(particles.x, particles.y,'r.')
    native_plot.xlim(-5,5)
    native_plot.ylim(-5,5)
    native_plot.show()
