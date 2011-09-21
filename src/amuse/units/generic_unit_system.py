"""
The generic unit system knows the seven base quantities in the
International System of Quantities, I.S.Q.

+-------------------+-----------------------------------+-----------------+
|Base quantity      |Name in generic unit               |Name in S.I. unit|
+-------------------+-----------------------------------+-----------------+
|length             |generic_system.length              |units.m          |
+-------------------+-----------------------------------+-----------------+
|time               |generic_system.time                |units.s          |
+-------------------+-----------------------------------+-----------------+
|mass               |generic_system.mass                |units.kg         |
+-------------------+-----------------------------------+-----------------+
|current            |generic_system.current             |units.A          |
+-------------------+-----------------------------------+-----------------+
|temperature        |generic_system.temperature         |units.K          |
+-------------------+-----------------------------------+-----------------+
|amount of substance|generic_system.amount_of_substance |units.mol        |
+-------------------+-----------------------------------+-----------------+
|luminous intensity |generic_system.luminous_intensity  |units.cd         |
+-------------------+-----------------------------------+-----------------+
"""
from amuse.units import units
from amuse.units import core

class generic_unit(core.base_unit):
    def __init__(self, unit_in_si, system):
        core.base_unit.__init__(self, unit_in_si.quantity, unit_in_si.name, unit_in_si.symbol, system)
        self.unit_in_si = unit_in_si

    def __str__(self):
        return self.unit_in_si.quantity


    def is_generic(self):
        return True
    
    
generic_system = core.system('generic')

length = generic_unit(units.m, generic_system)
time = generic_unit(units.s, generic_system)
mass = generic_unit(units.kg, generic_system)
current = generic_unit(units.A, generic_system)
temperature = generic_unit(units.K, generic_system)
luminous_intensity = generic_unit(units.cd, generic_system)

acceleration = length / (time ** 2)
force = mass*acceleration
potential = (length ** 2) / (time ** 2)
energy = mass * potential
specific_energy = potential
speed = length / time
volume = (length ** 3)
density = mass / volume
momentum_density = density * speed
energy_density = density * specific_energy
charge = current * time
pressure = mass / length / (time ** 2)

def is_generic_unit(unit):
    for factor, x in unit.base:
        if x.is_generic():
            return True
    return False
