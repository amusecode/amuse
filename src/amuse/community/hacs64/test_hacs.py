import collections
import getopt
import numpy
import os
import random
import sys
import unittest

from amuse.community.hacs64.interface import Hacs64 as grav
from amuse.units import nbody_system
from amuse.units import units

from amuse import datamodel
from amuse.datamodel import particle_attributes
from amuse.rfi.core import is_mpd_running
from amuse.ic.plummer import new_plummer_model
from amuse.ic.salpeter import new_salpeter_mass_distribution_nbody


def print_log(time, gravity, E0 = 0.0 | nbody_system.energy):
    M = gravity.total_mass
    U = gravity.potential_energy
    T = gravity.kinetic_energy
    Ebin = 0.0 | nbody_system.energy # gravity.get_binary_energy()
    Etop = T + U
    E = Etop + Ebin
    if E0 == 0 | nbody_system.energy: E0 = E
    Rv = -0.5*M*M/U
    Q = -T/U
    print ""
    print "time =", time.number, " energy = ", E.number, \
	" dE/E0 = ", (E/E0 - 1)
    print '%s %.4f %.6f %.6f %.6f %.6f %.6f %.6f %.6f' % \
	("%%", time.number, M.number, T.number, U.number, \
         E.number, Ebin.number, Rv.number, Q)
    sys.stdout.flush()
    return E

def test_hacs(infile = None,
              number_of_stars = 128,
              nmax = 2048,
              end_time = 0.1   | nbody_system.time,
              delta_t = 0.125 | nbody_system.time,
              dt_max  = 0.0625 | nbody_system.time,
              n_ngb   = 16,
              eta_irr = 0.6,
              eta_reg = 0.1,
              softening_length = 0.0 | nbody_system.length):

    if infile != None: print "input file =", infile
    print "end_time =",  end_time.number
    print "nstars= ",    number_of_stars,
    print "nmax= ",      nmax,
    print "delta_t= ",   delta_t.number
    print "dt_max= ",    dt_max.number
    print "n_ngb= ",     n_ngb,
    print "eta_irr= ",   eta_irr
    print "eta_reg= ",   eta_reg
    print "eps2=    ",   softening_length.number**2
    print "\ninitializing the gravity module"
    sys.stdout.flush()

#    gravity = grav(number_of_workers = 1, redirection = "none", mode='cpu')
    gravity = grav(number_of_workers = 1, redirection = "none", mode='cpu')
    gravity.initialize_code()
    
    #-----------------------------------------------------------------

    if infile == None:

        print "making a Plummer model"
        stars = new_plummer_model(number_of_stars)

        id = numpy.arange(number_of_stars)
        stars.id = id+1

        print "setting particle masses and radii"
	#stars.mass = (1.0 / number_of_stars) | nbody_system.mass
        scaled_mass = new_salpeter_mass_distribution_nbody(number_of_stars) 
        stars.mass = scaled_mass
        stars.radius = 0.0 | nbody_system.length

        print "centering stars"
        stars.move_to_center()
        print "scaling stars to virial equilibrium"
        stars.scale_to_standard(smoothing_length_squared
                                    = gravity.parameters.eps2)

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

        stars = datamodel.Particles(number_of_stars)
        stars.id = id 
        print len(mass), len(pos), len(vel), len(id)
        stars.mass = mass | nbody_system.mass
        stars.position = pos | nbody_system.length
        stars.velocity = vel | nbody_system.speed
        stars.radius = 0. | nbody_system.length
        nmax = 2*len(mass)

    # print "IDs:", stars.id.number
    sys.stdout.flush()

    #-----------------------------------------------------------------

    gravity.parameters.nmax    = nmax;
    gravity.parameters.dtmax   = dt_max;
#    gravity.parameters.n_ngb   = n_ngb;
    gravity.parameters.eta_irr = eta_irr;
    gravity.parameters.eta_reg = eta_reg;
    gravity.parameters.eps2    = softening_length**2
    gravity.commit_parameters();


    print "adding particles"
    # print stars
    sys.stdout.flush()
    gravity.particles.add_particles(stars)
    gravity.commit_particles()

    print ''
    print "number_of_stars =", number_of_stars
    print "evolving to time =", end_time.number, \
          "in steps of", delta_t.number
    sys.stdout.flush()

    E0 = print_log(time, gravity)
    
    # Channel to copy values from the code to the set in memory.
    channel = gravity.particles.new_channel_to(stars)

    stopping_condition = gravity.stopping_conditions.collision_detection
    stopping_condition.enable()

#    stopping_condition.disable()

    while time < end_time:
        if (gravity.get_time() >= time):
          time += delta_t

        gravity.evolve_model(time)


        # Ensure that the stars list is consistent with the internal
        # data in the module.

        ls = len(stars)

    	  # Update the bookkeeping: synchronize stars with the module data.

        # this breaks the code ...
        channel.copy()
  
        # Copy values from the module to the set in memory.
        
        channel.copy_attribute("index_in_code", "id")
    
        # Copy the index (ID) as used in the module to the id field in
        # memory.  The index is not copied by default, as different
        # codes may have different indices for the same particle and
        # we don't want to overwrite silently.

        if stopping_condition.is_set():
            star1 = stopping_condition.particles(0)[0]
            star2 = stopping_condition.particles(1)[0]
            gravity.synchronize_model()
            print '\nstopping condition set at time', \
                gravity.get_time().number,'for:\n'
            print star1
            print ''
            print star2
            print ''
            gravity.particles.remove_particle(star1)
            gravity.particles.remove_particle(star2)
           
            gravity.recommit_particles();
            
            print 'ls=', len(stars)
            
            gravity.update_particle_set()
            gravity.particles.synchronize_to(stars)
            
            
            print 'ls=', len(stars)
            

        if len(stars) != ls:
           if 0:
                print "stars:"
                for s in stars:
                    print " ", s.id.number, s.mass.number, s.x.number, s.y.number, s.z.number
           else:
             print "number of stars =", len(stars)
           sys.stdout.flush()

        print_log(gravity.get_time(), gravity, E0)
        sys.stdout.flush()

    print ''
    print_log(gravity.get_time(), gravity, E0)
    sys.stdout.flush()
    gravity.stop()

if __name__ == '__main__':

    infile = None
    N = 4009
    dt_max  = 0.0625 | nbody_system.time
    n_ngb = 16
    eta_irr = 0.8
    eta_reg = 0.14

    eta_irr = 0.6
    eta_reg = 0.1

    t_end = 1.0 | nbody_system.time
    delta_t = 0.125 | nbody_system.time
    random_seed = -1
#    random_seed = 123
#    random_seed = 15345

    try:
        opts, args = getopt.getopt(sys.argv[1:], "a:c:d:e:f:gGn:s:t:w:")
    except getopt.GetoptError, err:
        print str(err)
        sys.exit(1)

    for o, a in opts:
        if o == "-a":
            accuracy_parameter = float(a)
        elif o == "-c":
            manage_encounters = int(a)
        elif o == "-d":
            delta_t = float(a) | nbody_system.time 
        elif o == "-e":
            softening_length = float(a) | nbody_system.length
        elif o == "-f":
            infile = a
        elif o == "-g":
            use_gpu = 0
        elif o == "-G":
            use_gpu = 0
            gpu_worker = 0
        elif o == "-n":
            N = int(a)
        elif o == "-s":
            random_seed = int(a)
        elif o == "-t":
            t_end = float(a) | nbody_system.time
        elif o == "-w":
            n_workers = int(a)
        else:
            print "unexpected argument", o

    if random_seed <= 0:
        numpy.random.seed()
        random_seed = numpy.random.randint(1, pow(2,31)-1)
    numpy.random.seed(random_seed)
    print "random seed =", random_seed
    
    Nmax = N*2
    softening_length = 0.0/N  | nbody_system.length

    assert is_mpd_running()
    test_hacs(infile, 
              N,
              Nmax,
              t_end, 
              delta_t, 
              dt_max,
              n_ngb,
              eta_irr,
              eta_reg,
              softening_length)
