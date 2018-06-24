import numpy
from amuse.units import quantities
from amuse.units.units import rad, deg, rev


# trigonometric convenience functions which are "unit aware"
def sin(x):
    return numpy.sin(1.*x)


def cos(x):
    return numpy.cos(1.*x)


def tan(x):
    return numpy.tan(1.*x)


def arcsin(x):
    return numpy.arcsin(x) | rad


def arccos(x):
    return numpy.arccos(x) | rad


def arctan(x):
    return numpy.arctan(x) | rad


def arctan2(x, y):
    return numpy.arctan2(x, y) | rad


def to_rad(angle):
    return quantities.as_quantity_in(angle, rad)


def to_deg(angle):
    return quantities.as_quantity_in(angle, deg)


def to_rev(angle):
    return quantities.as_quantity_in(angle, rev)


def in_rad(angle):
    return quantities.value_in(angle, rad)


def in_deg(angle):
    return quantities.value_in(angle, deg)


def in_rev(angle):
    return quantities.value_in(angle, rev)
