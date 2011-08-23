from amuse.io.base import *
from amuse.io import text
from amuse.io.nemotsf import NemoFileFormatProcessor
from amuse.io.starlab import StarlabFileFormatProcessor
from amuse.io.store import HDF5FileFormatProcessor
from amuse.io.gadget import GadgetFileFormatProcessor
from amuse.io.vtk import VtkStructuredGrid
from amuse.io.vtk import VtkUnstructuredGrid

text.TableFormattedText.register()
text.CsvFileText.register()

NemoFileFormatProcessor.register()
StarlabFileFormatProcessor.register()
HDF5FileFormatProcessor.register()
GadgetFileFormatProcessor.register()
VtkStructuredGrid.register()
VtkUnstructuredGrid.register()
