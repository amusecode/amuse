import collections
import getopt
import numpy
import os
import random
import sys
import unittest

from amuse.community.ph4.interface import ph4 as grav
from amuse.community.smalln.interface import SmallN
from amuse.community.kepler.interface import Kepler
from amuse.couple import multiples

from amuse.units import nbody_system
from amuse.units import units
from amuse.units import quantities

from amuse import datamodel
from amuse.datamodel import particle_attributes as pa
from amuse.rfi.core import is_mpd_running
from amuse.ic.plummer import new_plummer_model
from amuse.ic.salpeter import new_salpeter_mass_distribution_nbody

def print_log(pre, time, gravity, E0 = 0.0 | nbody_system.energy):
    N = len(gravity.particles)
    M = gravity.total_mass
    U = gravity.potential_energy
    T = gravity.kinetic_energy
    Ebin = gravity.multiples_energy_correction
    Etop = T + U
    E = Etop - Ebin
    if E0 == 0 | nbody_system.energy: E0 = E
    Rvir = -0.5*M*M/U
    Q = -T/U
    com = pa.center_of_mass(gravity.particles)
    comv = pa.center_of_mass_velocity(gravity.particles)
    dcen,rcore,rhocore = pa.densitycentre_coreradius_coredens(gravity.particles)
    cmx,cmy,cmz = dcen
    lagr,mf = pa.LagrangianRadii(gravity.particles, dcen)  # no units!

    print ''
    print pre+"time=", time.number
    print pre+"Ntot=", N
    print pre+"mass=", M.number
    print pre+"Etot=", E.number
    print pre+"Ebin=", Ebin.number
    print pre+"dE/E=", E/E0 - 1
    print pre+"Rvir=", Rvir.number
    print pre+"Qvir=", Q
    cmx,cmy,cmz = com
    print pre+"cmpos[3]= %.8f %.8f %.8f" % (cmx.number, cmy.number, cmz.number)
    cmx,cmy,cmz = comv
    print pre+"cmvel[3]= %.8f %.8f %.8f" % (cmx.number, cmy.number, cmz.number)
    cmx,cmy,cmz = dcen
    print pre+"dcpos[3]= %.8f %.8f %.8f" % (cmx.number, cmy.number, cmz.number)
    print pre+"Rcore=", rcore.number
    print pre+"Mlagr[9]=",
    for m in mf: print "%.4f" % (m),
    print ''
    print pre+"Rlagr[9]=",
    for r in lagr.number: print "%.8f" % (r),
    print ''
    kT = T/N
    Nmul,Nbin,Emul = gravity.pretty_print_multiples(pre, kT, dcen)
    print pre+"Nmul=", Nmul
    print pre+"Nbin=", Nbin
    print pre+"Emul= %.5f" % (Emul.number)
    print pre+"Emul/kT= %.5f" % (Emul/kT)
    print pre+"Emul/E= %.5f" % (Emul/E)
    print ''

    sys.stdout.flush()

    return E

SMALLN = None
def new_smalln():
    SMALLN.reset()
    return SMALLN
    
def init_smalln():
    global SMALLN
    SMALLN = SmallN()
    SMALLN.parameters.timestep_parameter = 0.1
    SMALLN.parameters.cm_index = 2001

def init_kepler(star1, star2):        
    try:
        star1.mass.value_in(units.kg) # see if SI units, throw exception if not
        unit_converter \
            = nbody_system.nbody_to_si(star1.mass + star2.mass,
                                       (star2.position-star1.position).length())
    except Exception as ex:
        unit_converter = None
        
    kep = Kepler(unit_converter, redirection = "none")
    kep.initialize_code()

    return kep
        
def test_ph4(infile = None, number_of_stars = 40,
             end_time = 10 | nbody_system.time,
             delta_t = 1 | nbody_system.time,
             n_workers = 1, use_gpu = 1, gpu_worker = 1,
             salpeter = 0,
             accuracy_parameter = 0.1,
             softening_length = 0.0 | nbody_system.length,
             manage_encounters = 1, random_seed = 1234):

    if random_seed <= 0:
        numpy.random.seed()
        random_seed = numpy.random.randint(1, pow(2,31)-1)
    numpy.random.seed(random_seed)
    print "random seed =", random_seed

    if infile != None: print "input file =", infile
    print "end_time =", end_time.number
    print "delta_t =", delta_t.number
    print "n_workers =", n_workers
    print "use_gpu =", use_gpu
    print "manage_encounters =", manage_encounters
    print "\ninitializing the gravity module"
    sys.stdout.flush()

    init_smalln()
    
    # Note that there are actually three GPU options:
    #
    #	1. use the GPU code and allow GPU use (default)
    #	2. use the GPU code but disable GPU use (-g)
    #	3. use the non-GPU code (-G)

    if gpu_worker == 1:
        try:
            gravity = grav(number_of_workers = n_workers,
                           redirection = "none", mode = "gpu")
        except Exception as ex:
            gravity = grav(number_of_workers = n_workers, redirection = "none")
    else:
        gravity = grav(number_of_workers = n_workers, redirection = "none")

    gravity.initialize_code()
    gravity.parameters.set_defaults()

    #-----------------------------------------------------------------

    if infile == None:

        print "making a Plummer model"
        stars = new_plummer_model(number_of_stars)

        id = numpy.arange(number_of_stars)
        stars.id = id+1

        print "setting particle masses and radii"
        if salpeter == 0:
            print 'equal masses'
            scaled_mass = (1.0 / number_of_stars) | nbody_system.mass
        else:
            print 'salpeter mass function'
            scaled_mass = new_salpeter_mass_distribution_nbody(number_of_stars) 
        stars.mass = scaled_mass

        print "centering stars"
        stars.move_to_center()
        print "scaling stars to virial equilibrium"
        stars.scale_to_standard(smoothing_length_squared
                                    = gravity.parameters.epsilon_squared)



        if 1:
            nbin = 20
            mfac = number_of_stars/(number_of_stars+float(nbin))
            scaled_mass *= mfac
            stars.mass = scaled_mass

            kep = Kepler(redirection = "none")
            kep.initialize_code()
            mbin = 2*scaled_mass
            amin = (0.1/number_of_stars) | nbody_system.length
            amax = 4.*amin
            ecc = 0.1

            id_count = number_of_stars
            for i in range(0, number_of_stars, number_of_stars/nbin):
                # star i is CM, becomes component, add other star at end.
                newstar = datamodel.Particles(1)
                a = amin + numpy.random.uniform()*(amax-amin)
                kep.initialize_from_elements(mbin, a, ecc)
                dr = quantities.AdaptingVectorQuantity()
                dr.extend(kep.get_separation_vector())
                dv = quantities.AdaptingVectorQuantity()
                dv.extend(kep.get_velocity_vector())
                cm = stars[i]
                newstar.mass = scaled_mass
                stars[i].mass = scaled_mass
                newstar.position = cm.position + 0.5*dr
                newstar.velocity = cm.velocity + 0.5*dv
                stars[i].position = cm.position - 0.5*dr
                stars[i].velocity = cm.velocity - 0.5*dv
                id_count += 1
                newstar.id = id_count
                stars.add_particles(newstar)

            print 'created', nbin, 'binaries'
            kep.stop()

            number_of_stars += nbin

        # Dynamical radii (assumes virial equilibrium):
        stars.radius = (2*stars.mass.number) | nbody_system.length

        time = 0.0 | nbody_system.time
        sys.stdout.flush()

    else:

        # Read the input data.  Units are dynamical (sorry).
        # Format:  id  mass  pos[3]  vel[3]

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
        stars.mass = mass | nbody_system.mass
        stars.position = pos | nbody_system.length
        stars.velocity = vel | nbody_system.speed
        #stars.radius = 0. | nbody_system.length

        # Dynamical radii (assumes virial equilibrium):
        stars.radius = (2*stars.mass.number) | nbody_system.length

    # print "IDs:", stars.id.number
    sys.stdout.flush()

    #-----------------------------------------------------------------

    if softening_length < 0.0 | nbody_system.length:

        # Use ~interparticle spacing.

        eps2 = 0.25*(float(number_of_stars))**(-0.666667) \
			| nbody_system.length**2
    else:
        eps2 = softening_length*softening_length
    print 'softening length =', eps2.sqrt()

    gravity.parameters.timestep_parameter = accuracy_parameter
    gravity.parameters.epsilon_squared = eps2
    gravity.parameters.use_gpu = use_gpu
    # gravity.parameters.manage_encounters = manage_encounters

    print ''
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

    # Channel to copy values from the code to the set in memory.
    channel = gravity.particles.new_channel_to(stars)

    stopping_condition = gravity.stopping_conditions.collision_detection
    stopping_condition.enable()

    # -----------------------------------------------------------------
    # Create the coupled code and integrate the system to the desired
    # time, managing interactions internally.

    pre = "%%% "
    kep = init_kepler(stars[0], stars[1])
    multiples_code = multiples.Multiples(gravity, new_smalln, kep)
    E0 = print_log(pre, time, multiples_code)

    while time < end_time:
        time += delta_t
        multiples_code.evolve_model(time)

        # Copy values from the module to the set in memory.

        channel.copy()
    
        # Copy the index (ID) as used in the module to the id field in
        # memory.  The index is not copied by default, as different
        # codes may have different indices for the same particle and
        # we don't want to overwrite silently.

        channel.copy_attribute("index_in_code", "id")

        E = print_log(pre, time, multiples_code, E0)
        sys.stdout.flush()

    #-----------------------------------------------------------------

    # Write data to a file.

    print ''
    gravity.stop()

if __name__ == '__main__':

    print 'command line:',
    for i in range(len(sys.argv)): print sys.argv[i],
    print '\n'

    infile = None
    N = 100
    t_end = 5.0 | nbody_system.time
    delta_t = 1.0 | nbody_system.time
    n_workers = 2
    use_gpu = 1
    gpu_worker = 1
    salpeter = 0
    accuracy_parameter = 0.1
    softening_length = 0  | nbody_system.length
    random_seed = -1
    manage_encounters = 1

    try:
        opts, args = getopt.getopt(sys.argv[1:], "a:c:d:e:f:gGn:s:St:w:")
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
        elif o == "-S":
            salpeter = 1
        elif o == "-t":
            t_end = float(a) | nbody_system.time
        elif o == "-w":
            n_workers = int(a)
        else:
            print "unexpected argument", o

    assert is_mpd_running()
    test_ph4(infile, N, t_end, delta_t, n_workers,
             use_gpu, gpu_worker,
             salpeter, accuracy_parameter, softening_length,
             manage_encounters, random_seed)
