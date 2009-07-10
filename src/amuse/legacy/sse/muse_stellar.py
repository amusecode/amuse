#!/usr/bin/env python
# -*- mode: python; coding: utf-8 -*-

__author__       = 'Jarrod Hurley, Steve McMillan'
__author_email__ = '<steve@physics.drexel.edu>'
__date__         = '2009-02-05'
__version__      = 1.2

"""
SSE stellar evolution module.
"""

import os
import math
import sys

import stellar.single.SSE.SSE_muse_interface as sse
from stellar.stellar_module_template import Stellar
from stellar.single.stellar_state import stellar_state
from stellar.binary.binary_state import binary_state
from muse.utilities import build_muse_class

# Locate the MUSE directory.
if None:
    if not os.environ.has_key('MUSEDIR'):
      cwd = os.getcwd()
      while not os.path.exists(cwd + '/.museroot') and cwd != '/':
        cwd = os.path.split(cwd)[0]
      if cwd != '/':
        musedir = cwd
      else:
        print 'Environment variable MUSEDIR is not set, and can\'t be determined'
        print 'automatically by looking for .museroot file.'
        sys.exit(1)
    else:
      musedir = os.environ['MUSEDIR']

# This module is designed to complete the interface between python and
# the muse_stellar interface.  Most of the interface specifications
# are implemented here.
#
# Functions actually performed by the module code (sse.* below) are:
#
#        initialize(metallicity)
#        evolve()
#	 delta_t()

class SSE(object):

  stars = {}

  # The "stars" hash contains the following:
  #
  # "zams_mass"     - ZAMS mass of the star
  # "epoch"	    - time of birth of the star (always 0 here - handled
  #		      by the star manager)
  # "SSE age"       - current age of the star in SSE (model time in Twin)
  # "age"           - latest time used in the evolve() method
  # "type"	    - current SSE type index (kw) of the star
  # "mass"          - current mass of the star
  # "core_mass"     - current core mass of the star
  # "envelope_mass" - current envelope mass of the star
  # "radius"        - current radius of the star
  # "core_radius"   - current core radius of the star
  # "envelope_radius" - current envelope radius (r-rc) of the star
  # "luminosity"    - current luminosity of the star
  # "spin"          - current spin of the star
  # "t_ms"          - main sequence lifetime of the star
  #
  # Also, we retain "previous" versions of several variables for
  # interpolation purposes.

  # Strings for SSE standard types, as defined in sse/sse.f (starts at 0)

  SSE_descriptor = ('Low Mass MS Star',
                    'Main sequence Star',
                    'Hertzsprung Gap',
                    'Giant Branch',
                    'Core Helium Burning',
                    'First AGB',
                    'Second AGB',
                    'Naked Helium MS',
                    'Naked Helium HG',
                    'Naked Helium GB',
                    'Helium WD',
                    'Carbon/Oxygen WD',
                    'Oxygen/Neon WD',
                    'Neutron Star',
                    'Black Hole',
                    'Massless Remnant')

  status = 0

  def __init__(self, metallicity = 0.02):
    object.__init__(self)

    # Rather than requring an evolve.in file in the run directory,
    # read it in here and pass the relevant parameters, possibly
    # (someday) modified in the constructor argument, from here.
    #
    # Fortran code:
    #
    # READ(22,*)mass,z,tphysf				# ignore
    # READ(22,*)neta,bwind,hewind,sigma			# real*8
    # READ(22,*)ifflag,wdflag,bhflag,nsflag,mxns	# 4 integer + real*8
    # READ(22,*)pts1,pts2,pts3				# real*8
    base = os.path.split(__file__)[0]
    input_file = os.path.join(base, 'sse', 'evolve.in')
    try:
      f = open(input_file, 'r')
    except:
      print "Can't open SSE input file", input_file
      sys.exit(1)

    lines = f.readlines()

    l = lines[1].strip().split()
    neta   = float(l[0])
    bwind  = float(l[1])
    hewind = float(l[2])
    sigma  = float(l[3])

    l = lines[2].strip().split()
    ifflag = int(l[0])
    wdflag = int(l[1])
    bhflag = int(l[2])
    nsflag = int(l[3])
    mxns   = float(l[4])

    l = lines[3].strip().split()
    pts1 = float(l[0])
    pts2 = float(l[1])
    pts3 = float(l[2])

    # List parameters:

    self.parameters = {'metallicity':metallicity}

    # List other properties:
    
    self.properties = {'name':'SSE',
                       'author':'Jarrod Hurley',
                       'interface author':'Steve McMillan',
                       'description':'SSE stellar evolution',
                       'maintainer':'Steve McMillan',
                       'e-mail':'steve@physics.drexel.edu'}

    self.status = sse.initialize(metallicity,
                                 neta, bwind, hewind, sigma,
                                 ifflag, wdflag, bhflag, nsflag, mxns,
                                 pts1, pts2, pts3)

  # Define the functions listed in stellar.py.

  def setup_module(self):	# do nothing
    return self.status

  def add_zams_star(self, id, mass):

    # Initialize stellar model and determine physical parameters.

    # Set up essential parameters, call evolve(epsilon) to set up the
    # rest.

    self.stars[id] = {"zams_mass":mass, 
                      "epoch":0.0,	# time of birth is always 0 here
                      "age":0.0,
                      "SSE age":0.0,
                      "type":1,
                      "mass":mass,
                      "core_mass":0.0,
                      "envelope_mass":0.0,
                      "radius":0.0,
                      "core_radius":0.0,
                      "envelope_radius":0.0,
                      "luminosity":0.0,
                      "spin":0.0,
                      "t_ms":0.0}

    epsilon = 1.e-6
    self.evolve(id, epsilon)
    self.stars[id]['age'] = 0.0
    self.stars[id]['previous SSE age'] = 0.0
    self.stars[id]['SSE age'] = 0.0

    return len(self.stars),0

  def remove_star(self, id):
    del(self.stars[id])
    return len(self.stars),0

  def get_number(self):
    return len(self.stars),0

  def get_age(self, id):
    return self.stars[id]['age']

  def evolve(self, id, time):

    self.stars[id]['age'] = time

    if self.stars[id]['SSE age'] <= time:

      # Save 'previous' entries for all quantities to be interpolated below.

      self.stars[id]['previous SSE age'] = self.stars[id]['SSE age']
      self.stars[id]['previous mass'] = self.stars[id]['mass']
      self.stars[id]['previous radius'] = self.stars[id]['radius']
      self.stars[id]['previous luminosity'] = self.stars[id]['luminosity']
      self.stars[id]['previous spin'] = self.stars[id]['spin']

      self.stars[id]['type'], \
      self.stars[id]['zams_mass'], \
      self.stars[id]['mass'], \
      self.stars[id]['radius'], \
      self.stars[id]['luminosity'], \
      self.stars[id]['core_mass'], \
      self.stars[id]['core_radius'], \
      self.stars[id]['envelope_mass'], \
      self.stars[id]['envelope_radius'], \
      self.stars[id]['spin'], \
      self.stars[id]['epoch'], \
      self.stars[id]['t_ms'], \
      self.stars[id]['SSE age'], \
      new_time \
	= sse.evolve(self.stars[id]['type'],
                     self.stars[id]['zams_mass'],
                     self.stars[id]['mass'],
                     self.stars[id]['radius'],
                     self.stars[id]['luminosity'],
                     self.stars[id]['core_mass'],
                     self.stars[id]['core_radius'],
                     self.stars[id]['envelope_mass'],
                     self.stars[id]['envelope_radius'],
                     self.stars[id]['spin'],
                     self.stars[id]['epoch'],
                     self.stars[id]['t_ms'],
                     self.stars[id]['SSE age'],
                     time)

      # Notes on the operation of evolve:
      #
      # 1. On return, new_time will be the new age of the star, and
      #    age will be set equal to it.
      #
      # 2. The new age will not necessarily be the requested time,
      #    i.e. new_time >= time, so we interpolate when properties at
      #    a specific time are needed (see below).
      #
      # 3. The SSE code adjusts the epoch to handle mass loss, WD
      #    formation, etc., so don't touch it!

    return 0

  def get_time_step(self, id, dm = 0.01, dr = 0.01):

    # This function returns the time step limited by the specified
    # changes in mass and/or radius.

    return sse.get_time_step(self.stars[id]['type'], self.stars[id]['zams_mass'], self.stars[id]['age'], self.stars[id]['mass'], self.stars[id]['t_ms'], self.stars[id]['epoch'])

  def has_real_stellar_model (self, id):
    return False

  def get_relative_velocity(self, id):
    return 0

  def get_initial_mass(self, id):
    return self.stars[id]['zams_mass']
  
  # Interpolate any quantity between the two stored models.

  def interpolate(self, id, time, prop):
    tprev = self.stars[id]['previous SSE age']
    tcurr = self.stars[id]['SSE age']
    prev_prop = ''.join(['previous ', prop])
    xprev = self.stars[id][prev_prop]
    xcurr = self.stars[id][prop]
    delta = tcurr-tprev
    if delta > 0:
      return xprev + (time-tprev)*(xcurr-xprev)/(tcurr-tprev)
    else:
      return xcurr
    
  def get_mass(self, id):
    return self.interpolate(id, self.stars[id]['SSE age'], 'mass')

  def get_radius(self, id):
    return self.interpolate(id, self.stars[id]['SSE age'], 'radius')

  def get_luminosity(self, id):
    return self.interpolate(id, self.stars[id]['SSE age'], 'luminosity')

  def get_spin(self, id):
    return self.interpolate(id, self.stars[id]['SSE age'], 'spin')

  def get_effective_temperature(self, id):
    r = self.get_radius(id)
    l = self.get_luminosity(id)
    return 5700*math.sqrt(math.sqrt(l)/r)

  def get_spin(self, id):
    return self.stars[id]['spin']

  # get_*_state: populate and return a stellar/binary state object.

  def get_stellar_state(self, id):
     return stellar_state(id = id,
                          age = self.get_age(id),
                          zams_mass = self.stars[id]['zams_mass'],
                          mass = self.get_mass(id),
                          radius = self.get_radius(id),
                          effective_temperature \
 				= self.get_effective_temperature(id),
                          luminosity = self.get_luminosity(id),
                          metallicity = self.parameters['metallicity'],
                          spin = self.get_spin(id),
                          module = self.properties['name'],
                          file = None,
                          data_pointer = None)

  def get_binary_state(self, id):
       return binary_state(id = id,
                           primary_id = 0,
                           secondary_id = 0,
                           period = 0.0,
                           eccentricity = 0.0,
                           module = None)


  def get_description(self, id):
    return self.SSE_descriptor[self.stars[id]['type']]

  #----------------------------------------------------------------------

  # Not needed, but included to satisfy references in the template.

  def add_binary():			Pass
  def get_model_time():			Pass
  def get_previous_model_time():	Pass
  def join_binary():			Pass
  def use_model_interpolation():	Pass

build_muse_class(SSE, sse, Stellar)

# eof
