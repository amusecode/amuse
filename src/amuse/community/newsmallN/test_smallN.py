import sys
import unittest
import numpy
import random
import collections
import getopt
import os

from amuse.support.units import nbody_system
from amuse.support.units import units
from amuse.support.data import core
from amuse.support.data import particle_attributes
from amuse.support.codes.core import is_mpd_running
from amuse.ext.plummer import MakePlummerModel
from amuse.ext.salpeter import new_salpeter_mass_distribution_nbody

from amuse.support.data import trees

from amuse.community.newsmallN.interface import smallN as grav
#from amuse.community.hermite0.interface import Hermite as grav

def print_log(time, gravity, E0 = 0.0 | nbody_system.energy):
    M = gravity.total_mass
    U = gravity.potential_energy
    T = gravity.kinetic_energy
    E = T + U
    if E0 == 0 | nbody_system.energy: E0 = E
    Rv = -0.5*M*M/U
    Q = -T/U
    print ""
    print "time =", time.number, " energy = ", E.number, \
	" dE/E0 = ", (E/E0 - 1).number
    print '%s %.4f %.6f %.6f %.6f %.6f %.6f %.6f' % \
	("%%", time.number, M.number, T.number, U.number, \
         E.number, Rv.number, Q.number)
    sys.stdout.flush()
    return E

def test_smallN(infile = None, number_of_stars = 10,
             end_time = 10 | nbody_system.time,
             delta_t = 1 | nbody_system.time,
             accuracy_parameter = 0.1):

    if infile != None: print "input file =", infile
    print "end_time =", end_time.number
    print "delta_t =", delta_t.number
    print "\ninitializing the gravity module"
    sys.stdout.flush()

    gravity = grav(redirection = "none")
    gravity.initialize_code()
    gravity.parameters.set_defaults()

    #-----------------------------------------------------------------

    if infile == None:

        print "making a Plummer model"
        stars = MakePlummerModel(number_of_stars).result

        id = numpy.arange(number_of_stars)
        stars.id = id+1 | units.none

        print "setting particle masses and radii"
	#stars.mass = (1.0 / number_of_stars) | nbody_system.mass
        scaled_mass = new_salpeter_mass_distribution_nbody(number_of_stars) 
        stars.mass = scaled_mass
        stars.radius = 0.0 | nbody_system.length

        print "centering stars"
        stars.move_to_center()
        print "scaling stars to virial equilibrium"
        stars.scale_to_standard(smoothing_length_squared = 0 \
				 | nbody_system.length*nbody_system.length)

        time = 0.0 | nbody_system.time
        sys.stdout.flush()

    else:

        # Read the input data.  Units are dynamical.

        print "reading file", infile

        id = []
        mass = []
        pos = []
        vel = []

        f = open(infile, 'r')
	count = 0
        for line in f:
            if len(line) > 0:
                count += 1
		cols = line.split()
                if count == 1: snap = int(cols[0])
		elif count == 2: number_of_stars = int(cols[0])
		elif count == 3: time = float(cols[0]) | nbody_system.time
		else:
                    if len(cols) >= 8:
                        id.append(int(cols[0]))
                        mass.append(float(cols[1]))
                        pos.append((float(cols[2]),
                                    float(cols[3]), float(cols[4])))
                        vel.append((float(cols[5]),
                                    float(cols[6]), float(cols[7])))
	f.close()

        stars = core.Particles(number_of_stars)
        stars.id = id | units.none
        stars.mass = mass | nbody_system.mass
        stars.position = pos | nbody_system.length
        stars.velocity = vel | nbody_system.speed
        stars.radius = 0. | nbody_system.length

    # print "IDs:", stars.id.number
    sys.stdout.flush()

    #-----------------------------------------------------------------

    gravity.parameters.timestep_parameter = accuracy_parameter

    print "adding particles"
    # print stars
    sys.stdout.flush()
    gravity.particles.add_particles(stars)
    print "committing particles"
    gravity.commit_particles()

    print ''
    print "number_of_stars =", number_of_stars
    print "evolving to time =", end_time.number, \
          "in steps of", delta_t.number
    sys.stdout.flush()

    E0 = print_log(time, gravity)
    
    # Channel to copy values from the code to the set in memory.
    channel = gravity.particles.new_channel_to(stars)

    while time < end_time:
        time += delta_t
        gravity.evolve_model(time)

        # Ensure that the stars list is consistent with the internal
        # data in the module.

        ls = len(stars)

	# Update the bookkeeping: synchronize stars with the module data.

        try:
            gravity.update_particle_set()
            gravity.particles.synchronize_to(stars)
        except:
            pass

        # Copy values from the module to the set in memory.

        channel.copy()
    
        # Copy the index (ID) as used in the module to the id field in
        # memory.  The index is not copied by default, as different
        # codes may have different indices for the same particle and
        # we don't want to overwrite silently.

        channel.copy_attribute("index_in_code", "id")

        if len(stars) != ls:
            if 0:
                print "stars:"
                for s in stars:
                    print " ", s.id.number, s.mass.number, \
			       s.x.number, s.y.number, s.z.number
            else:
		print "number of stars =", len(stars)
            sys.stdout.flush()

        print_log(time, gravity, E0)

        over = gravity.is_over()
        if over.number:
            print 'interaction is over'
            gravity.update_particle_tree()
            gravity.update_particle_set()
            gravity.particles.synchronize_to(stars)
            channel.copy()
            channel.copy_attribute("index_in_code", "id")
            for s in stars:
                print s.child1

            x = trees.BinaryTreesOnAParticleSet(stars, "child1", "child2")
            #roots = list(x.iter_roots())
            #for r in roots:
            #    print r.particle.id, r.particle.child1.id, r.particle.child2.id

            break
        else:
            print 'interaction is not over'
    
        sys.stdout.flush()

    gravity.stop()

if __name__ == '__main__':

    infile = None
    N = 10
    t_end = 5.0 | nbody_system.time
    delta_t = 1.0 | nbody_system.time
    accuracy_parameter = 0.1
    random_seed = -1

    try:
        opts, args = getopt.getopt(sys.argv[1:], "a:c:d:e:f:gGn:s:t:w:")
    except getopt.GetoptError, err:
        print str(err)
        sys.exit(1)

    for o, a in opts:
        if o == "-a":
            accuracy_parameter = float(a)
        elif o == "-d":
            delta_t = float(a) | nbody_system.time 
        elif o == "-f":
            infile = a
        elif o == "-n":
            N = int(a)
        elif o == "-s":
            random_seed = int(a)
        elif o == "-t":
            t_end = float(a) | nbody_system.time
        else:
            print "unexpected argument", o

    if random_seed <= 0:
        numpy.random.seed()
        random_seed = numpy.random.randint(1, pow(2,31)-1)
    numpy.random.seed(random_seed)
    print "random seed =", random_seed

    assert is_mpd_running()
    test_smallN(infile, N, t_end, delta_t, accuracy_parameter)
