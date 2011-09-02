from amuse.units import core

k = core.k

system = core.system('S.I.')

m = core.base_unit('length', 'meter', 'm', system)
kg = core.base_unit('mass', 'kilogram', 'kg', system)
s = core.base_unit('time', 'second', 's', system)
A = core.base_unit('electric current', 'ampere', 'A', system)
K = core.base_unit('thermodynamic temperature', 'kelvin', 'K', system)
mol = core.base_unit('amount of substance', 'mole', 'mol', system)
cd = core.base_unit('luminous intensity', 'candela', 'cd', system)

no_system = core.no_system
no_unit = core.none_unit('no unit','none')
none = core.none_unit('none','none')

named = core.named_unit

# SI prefixes
def deca(unit):
    return named('deca'+unit.name,'da'+unit.symbol,10.*unit)    
def hecto(unit):
    return named('hecto'+unit.name,'h'+unit.symbol,100.*unit)    
def kilo(unit):
    return named('kilo'+unit.name,'k'+unit.symbol,1000.*unit)    
def mega(unit):
    return named('mega'+unit.name,'M'+unit.symbol,1.e6*unit)    
def giga(unit):
    return named('giga'+unit.name,'G'+unit.symbol,1.e9*unit)    
def terra(unit):
    return named('terra'+unit.name,'T'+unit.symbol,1.e12*unit)    
def peta(unit):
    return named('peta'+unit.name,'P'+unit.symbol,1.e15*unit)    
def exa(unit):
    return named('exa'+unit.name,'E'+unit.symbol,1.e18*unit)    
def zetta(unit):
    return named('zetta'+unit.name,'Z'+unit.symbol,1.e21*unit)    
def yotta(unit):
    return named('yotta'+unit.name,'Y'+unit.symbol,1.e24*unit)    
def deci(unit):
    return named('deci'+unit.name,'d'+unit.symbol,0.1*unit)    
def centi(unit):
    return named('centi'+unit.name,'c'+unit.symbol,0.01*unit)    
def milli(unit):
    return named('milli'+unit.name,'m'+unit.symbol,0.001*unit)    
def micro(unit):
    return named('micro'+unit.name,'mu'+unit.symbol,1.e-6*unit)    
def nano(unit):
    return named('nano'+unit.name,'n'+unit.symbol,1.e-9*unit)    
def pico(unit):
    return named('pico'+unit.name,'p'+unit.symbol,1.e-12*unit)    
def femto(unit):
    return named('femto'+unit.name,'f'+unit.symbol,1.e-15*unit)    
def atto(unit):
    return named('atto'+unit.name,'a'+unit.symbol,1.e-18*unit)    
def zepto(unit):
    return named('zepto'+unit.name,'z'+unit.symbol,1.e-21*unit)    
def yocto(unit):
    return named('yocto'+unit.name,'y'+unit.symbol,1.e-24*unit)
