from amuse.support.io.base import *
from amuse.support.io import text
from amuse.support.io.nemotsf import NemoFileFormatProcessor
from amuse.support.io.starlab import StarlabFileFormatProcessor
from amuse.support.io.store import HDF5FileFormatProcessor


text.TableFormattedText.register()
text.CsvFileText.register()

NemoFileFormatProcessor.register()
StarlabFileFormatProcessor.register()
HDF5FileFormatProcessor.register()
