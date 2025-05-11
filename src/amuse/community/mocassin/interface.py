from amuse.support.import_helper import load_code
from amuse.units.units import eV

MocassinInterface = load_code("mocassin", "MocassinInterface")
Mocassin = load_code("mocassin", "Mocassin")

mocassin_rydberg_unit = 13.6 * eV
