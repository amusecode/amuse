"""
File input and output

AMUSE can read and write to/from several data formats. The main routines
are write_set_to_file and read_set_from_file, explained below.
"""

from amuse.io.base import *
from amuse.io import text
from amuse.io.nemotsf import NemoFileFormatProcessor
from amuse.io.nemobin import NemoBinaryFileFormatProcessor
from amuse.io.starlab import StarlabFileFormatProcessor
from amuse.io.store import HDF5FileFormatProcessor
from amuse.io.gadget import GadgetFileFormatProcessor
from amuse.io.vtk import VtkStructuredGrid
from amuse.io.vtk import VtkUnstructuredGrid

text.TableFormattedText.register()
text.CsvFileText.register()
text.AmuseText.register()

NemoFileFormatProcessor.register()
NemoBinaryFileFormatProcessor.register()
StarlabFileFormatProcessor.register()
HDF5FileFormatProcessor.register()
GadgetFileFormatProcessor.register()
VtkStructuredGrid.register()
VtkUnstructuredGrid.register()

__all__ = ["read_set_from_file", "write_set_to_file", "get_options_for_format"]
