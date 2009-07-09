#!/usr/bin/env python
# -*- mode: python; coding: utf-8 -*-

__author__       = 'Steve Mcmillan'
__author_email__ = '<steve@physics.drexel.edu>'
__date__         = '2008-10-24'
__version__      = 0.2

"""
Hermite Gravity module.
"""

# Watch out for
# . cvars
# . translation of fn names
# . missing fns (AttributeError exception)

import gravity.hermite0.muse_dynamics_interface as HI
from gravity.gravity import Gravity
from muse.utilities import build_muse_class

# This module is designed to complete the interface between python and
# the muse_dynamics interface.  The cvars variables are listed in
# parameters.h.  Here we make them accessible in the Hermite
# constructor, with reasonable defaults.

class Hermite(object):

  # Initializer defines default values for global variables.

  def __init__(self, dt_param=None, dt_dia=None, eps2=None,
               flag_collision=None):
    object.__init__(self)

    if dt_param is not None: HI.cvar.dt_param = dt_param
    if dt_dia   is not None: HI.cvar.dt_dia = dt_dia
    if eps2     is not None: HI.cvar.eps2 = eps2
    if flag_collision is not None: HI.cvar.flag_collision = flag_collision

    # List parameters for python listing:

    self.parameters = {'dt_param': dt_param, 'dt_dia': dt_dia,
                       'eps': eps2, 'flag_collision': flag_collision}

    # List other properties:
    
    self.properties = {'name':'hermite0',
                       'author':'Steve McMillan, Jun Makino, Piet Hut',
                       'interface author':'Steve McMillan',
                       'description':'C++ shared-timestep Hermite scheme',
                       'maintainer':'Steve McMillan',
                       'e-mail':'steve@physics.drexel.edu'}

  # Here's how to make parameters settable (by e.g. grav.eps2 = 0.1):

  def _get_dt_param(self):      return HI.cvar.dt_param
  def _set_dt_param(self, val): HI.cvar.dt_param = val
  dt_param = property(fget=_get_dt_param, fset=_set_dt_param)

  def _get_dt_dia(self):      return HI.cvar.dt_dia
  def _set_dt_dia(self, val): HI.cvar.dt_dia = val
  dt_dia = property(fget=_get_dt_dia, fset=_set_dt_dia)

  def _get_eps2(self):      return HI.cvar.eps2
  def _set_eps2(self, val): HI.cvar.eps2 = val
  eps2 = property(fget=_get_eps2, fset=_set_eps2)

  def _get_flag_collision(self):      return HI.cvar.flag_collision
  def _set_flag_collision(self, val): HI.cvar.flag_collision = val
  flag_collision = property(fget=_get_flag_collision, fset=_set_flag_collision)

  # Here's how we repackage, if necessary, an internal function
  # to conform with the expected definition:

  def find_colliding_primary(self):       return HI.get_colliding_primary()
  def find_colliding_secondary(self, id): return HI.get_colliding_secondary(id)

  # Here's how to access a "local" function (i.e. one not in gravity.py).

  def get_n_steps(self):		return HI.get_n_steps()

build_muse_class(Hermite, HI, Gravity)

# eof.

