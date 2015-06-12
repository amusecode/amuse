import numpy
from amuse.units import quantities
from amuse.units.units import rad,deg,rev,pi

#trigonometric convenience functions which are "unit aware"
sin=lambda x: numpy.sin(1.*x)
cos=lambda x: numpy.cos(1.*x)
tan=lambda x: numpy.tan(1.*x)
arcsin=lambda x: numpy.arcsin(x) | rad
arccos=lambda x: numpy.arccos(x) | rad
arctan=lambda x: numpy.arctan(x) | rad
arctan2=lambda x,y: numpy.arctan2(x,y) | rad

def to_rad(angle):
  return quantities.as_quantity_in(angle,rad)
def to_deg(angle):
  return quantities.as_quantity_in(angle,deg)
def to_rev(angle):
  return quantities.as_quantity_in(angle,rev)

def in_rad(angle):
  return quantities.value_in(angle,rad)
def in_deg(angle):
  return quantities.value_in(angle,deg)
def in_rev(angle):
  return quantities.value_in(angle,rev)
