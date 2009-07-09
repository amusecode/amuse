#!/usr/bin/env python

import os
import sys

# If the MUSEDIR environment variable is not set, try to determine it
# automatically by looking in successive parent directories for a
# .museroot file.

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

if not musedir in sys.path:
  sys.path.insert(0, musedir)

input_file = os.path.join(musedir, 'muse', 'in.20')

def print_diag():
  print "time = %.5f, dt = %.4e, E = %.6f, %d steps" % \
      (h.get_time(), h.get_time_step(),
       h.get_kinetic_energy() + h.get_potential_energy(), h.get_n_steps())

  # (Note: get_n_steps() is not a standard Gravity function.)

eps = 0.1

from muse_dynamics import Hermite
h = Hermite(dt_dia=1.e9,
            eps2=eps*eps,
            flag_collision=1)

h.setup_module()

dt_next = 0.25
t_max = 10
radius = 0.025

# Read data from the input file:

s = h.dynamics_state()

f = open(input_file,'r')
lines = f.readlines()

for line in lines:
  l = line.strip().split()
  if len(l) >= 8:
    s.id = int(l[0])
    s.mass = float(l[1])
    s.radius = radius
    s.x = float(l[2])
    s.y = float(l[3])
    s.z = float(l[4])
    s.vx = float(l[5])
    s.vy = float(l[6])
    s.vz = float(l[7])
    h.add_particle(s)

t_next = 0
h.initialize_particles(t_next)
print "n = %d, tdyn = %.5f" % (h.get_number(), h.get_dynamical_time_scale())
print_diag()

# Run the code:

while t_next < t_max:

  t_next += dt_next

  # Evolve is supposed to loop to time t_next, but this
  # secondary loop is necessary to handle collisions, which
  # cause an immediate return with the id of the primary.

  coll = 0
  while h.get_time() + h.get_time_step() <= t_next:

    id1 = h.evolve(t_next)

    # System is at time t_next unless a collision occurred.

    if id1 >= 0:

      h.evolve(h.get_time(), 1)	# (sync not actually necessary here)

      id2 = h.find_colliding_secondary(id1)
      if id2 >= 0:

        # Only flag one collision per output interval.

        if not coll:
          print "  detected collision of %d and %d at time %.5f" % \
              (id1, id2, h.get_time())
          coll = 1

  print_diag()

h.cleanup_module()
