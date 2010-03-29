from amuse.support.io.base import *
from amuse.support.io import text
from amuse.support.io.nemotsf import NemoFileFormatProcessor


text.TableFormattedText.register()
text.CsvFileText.register()

NemoFileFormatProcessor.register()
