from amuse.units.si import *

Hz = named('hertz', 'Hz', 1/s)
MHz = named('megahertz', 'MHz', 1e6*Hz)
rad = named('radian','rad',m/m)
sr = named('steradian','sr',m**2/m**2)
N = named('newton', 'N', kg * m /s**2)
Pa = named('pascal', 'Pa', N / (m ** 2))
J = named('joule','J', kg * m **2  * s ** -2)
W = named('watt', 'W', J / s)
F = named('farad','F', s**4*A**2*m**(-2)*kg**(-1))
C = named('coulomb','C', A*s)
V = named('volt','V', J/C)
T = named('tesla','T', kg/A/s/s)
tesla = T
ohm = named('ohm','ohm', V/A)
S = named('siemens', 'S', A/V)
Wb = named('weber','Wb', V*s)
weber = Wb

# handy definitions
one = 1 | none
km = k(m)
