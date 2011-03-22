from amuse.support.units.units import *

pound=named('avoirdupois pound','lbm',0.45359237 * kg)
g=9.80665 * m/ s**2

poundforce=named('poundforce','lbf',pound*g)
print poundforce.as_quantity_in(N)

feet=named('feet','feet',1/3.2808399 * m)
print g.value_in(feet/s**2)
print poundforce.as_quantity_in(pound*feet/s**2)
print (1| N).as_quantity_in(poundforce)
