from amuse.support.units import core

k = core.k

system = core.system('S.I.')

m = core.base_unit('length', 'meter', 'm', system)
kg = core.base_unit('mass', 'kilogram', 'kg', system)
s = core.base_unit('time', 'second', 's' , system)
A = core.base_unit('electric current', 'ampere', 'A', system)
K = core.base_unit('thermodynamic temperature ', 'kelvin', 'K', system)
mol = core.base_unit('amount of substance', 'mole', 'mol', system)
cd = core.base_unit('luminous intensity', 'candela', 'cd', system)

no_unit = core.none_unit('no unit','')

named = core.named_unit
