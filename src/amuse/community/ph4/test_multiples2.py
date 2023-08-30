import collections
import getopt
import numpy
import os
import random
import sys
import unittest
from time import process_time as clock

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

def print_log(pre, time, gravity, E0 = 0.0 | nbody_system.energy, cpu0 = 0.0):
    cpu = clock()
    N = len(gravity.particles)
    M = gravity.total_mass
    U = gravity.potential_energy
    T = gravity.kinetic_energy
    Etop = T + U
    Nmul, Nbin, Emul = gravity.get_total_multiple_energy()
    tmp1,tmp2,Emul2 = gravity.get_total_multiple_energy2()
    Etot = Etop + Emul
    Eext = gravity.multiples_external_tidal_correction
    Eint = gravity.multiples_internal_tidal_correction
    Eerr = gravity.multiples_integration_energy_error
    Edel = gravity.multiples_external_tidal_correction \
		+ gravity.multiples_internal_tidal_correction \
		+ gravity.multiples_integration_energy_error
    Ecor = Etot - Edel
    if E0 == 0 | nbody_system.energy: E0 = Ecor
    Rvir = -0.5*M*M/U
    Q = -T/U
    com = pa.center_of_mass(gravity.particles)
    comv = pa.center_of_mass_velocity(gravity.particles)
    dcen,rcore,rhocore = pa.densitycentre_coreradius_coredens(gravity.particles)
    cmx,cmy,cmz = dcen
    lagr,mf = pa.LagrangianRadii(gravity.particles, cm=dcen)  # no units!

    print('')
    print(pre+"time=", time.number)
    print(pre+"CPU=", cpu - cpu0)
    print(pre+"Ntot=", N)
    print(pre+"mass=", M.number)
    print(pre+"Etot=", Etot.number)
    print(pre+"Etop=", Etop.number)
    print(pre+"Eext=", Eext.number)
    print(pre+"Eint=", Eint.number)
    print(pre+"Eerr=", Eerr.number)
    print(pre+"Edel=", Edel.number)
    print(pre+"Ecor=", Ecor.number)
    print(pre+"dE/E=", Ecor/E0 - 1)
    print(pre+"Rvir=", Rvir.number)
    print(pre+"Qvir=", Q)
    cmx,cmy,cmz = com
    print(pre+"cmpos[3]= %.8f %.8f %.8f" % (cmx.number, cmy.number, cmz.number))
    cmx,cmy,cmz = comv
    print(pre+"cmvel[3]= %.8f %.8f %.8f" % (cmx.number, cmy.number, cmz.number))
    cmx,cmy,cmz = dcen
    print(pre+"dcpos[3]= %.8f %.8f %.8f" % (cmx.number, cmy.number, cmz.number))
    print(pre+"Rcore=", rcore.number)
    print(pre+"Mcore=", (rhocore*rcore**3).number)	# fake...
    print(pre+"Mlagr[9]=", end=' ')
    for m in mf: print("%.4f" % (m), end=' ')
    print('')
    print(pre+"Rlagr[9]=", end=' ')
    for r in lagr.number: print("%.8f" % (r), end=' ')
    print('')
    kT = T/N
    Nmul,Nbin,Emul = gravity.print_multiples2(pre, kT, dcen)
    print(pre+"Nmul=", Nmul)
    print(pre+"Nbin=", Nbin)
    print(pre+"Emul= %.5f" % (Emul.number))
    print(pre+"Emul2= %.5f" % (Emul2.number))
    print(pre+"Emul/kT= %.5f" % (Emul.number/kT.number))
    print(pre+"Emul/E= %.5f" % (Emul.number/Etot.number))
    print('')

    sys.stdout.flush()
    return Ecor,cpu

SMALLN = None
def new_smalln():
    SMALLN.reset()
    return SMALLN
    
def init_smalln():
    global SMALLN
    sys.stdout.flush()
    SMALLN = SmallN()
    sys.stdout.flush()
    SMALLN.parameters.timestep_parameter = 0.1
    #SMALLN.parameters.cm_index = 2001		# don't set this here!!
    sys.stdout.flush()

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

def run_ph4(infile = None, outfile = None,
            number_of_stars = 100, number_of_binaries = 0,
            end_time = 10 | nbody_system.time,
            delta_t = 1 | nbody_system.time,
            n_workers = 1, use_gpu = 1, gpu_worker = 1,
            salpeter = 0,
            accuracy_parameter = 0.1,
            softening_length = 0.0 | nbody_system.length,
            manage_encounters = 1, random_seed = 1234,
            debug_level = 1):

    if random_seed <= 0:
        numpy.random.seed()
        random_seed = numpy.random.randint(1, pow(2,31)-1)
    numpy.random.seed(random_seed)
    print("random seed =", random_seed)

    if infile != None: print("input file =", infile)
    print("end_time =", end_time)
    print("delta_t =", delta_t)
    print("n_workers =", n_workers)
    print("use_gpu =", use_gpu)
    print("manage_encounters =", manage_encounters)
    print("\ninitializing the gravity module")
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
            gravity = grav(number_of_workers = n_workers,
                           redirection = "none")
    else:
        gravity = grav(number_of_workers = n_workers,
                       redirection = "none")

    gravity.initialize_code()
    gravity.parameters.set_defaults()

    #-----------------------------------------------------------------

    if infile == None:

        print("making a Plummer model")
        stars = new_plummer_model(number_of_stars)

        id = numpy.arange(number_of_stars)
        stars.id = id+1

        print("setting particle masses and radii")
        if salpeter == 0:
            print('equal masses')
            total_mass = 1.0 | nbody_system.mass
            scaled_mass = total_mass / number_of_stars
        else:
            print('salpeter mass function')
            scaled_mass = new_salpeter_mass_distribution_nbody(number_of_stars) 
        stars.mass = scaled_mass

        print("centering stars")
        stars.move_to_center()
        print("scaling stars to virial equilibrium")
        stars.scale_to_standard(smoothing_length_squared
                                    = gravity.parameters.epsilon_squared)

    else:

        # Read the input data.  Units are dynamical (sorry).
        # Format:  id  mass  pos[3]  vel[3]

        print("reading file", infile)

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

    total_mass = stars.mass.sum()
    ke = pa.kinetic_energy(stars)
    kT = ke/(1.5*number_of_stars)

    if number_of_binaries > 0:

        # Turn selected stars into binary components.
        # Only tested for equal-mass case.

        kep = Kepler(redirection = "none")
        kep.initialize_code()

        added_mass = 0.0 | nbody_system.mass

        # Work with energies rather than semimajor axes.

        Emin = 10*kT
        Emax = 20*kT
        ecc = 0.1

        nbin = 0
        companion_base_id = 100*(number_of_stars//10)

        for i in range(0, number_of_stars,
                       number_of_stars//number_of_binaries):

            # Star i is CM, becomes component, add other star at end.

            nbin += 1

            mass = stars[i].mass
            new_mass = numpy.random.uniform()*mass	# uniform q
            mbin = mass + new_mass
            fac = new_mass/mbin
            E = Emin + numpy.random.uniform()*(Emax-Emin)
            a = 0.5*nbody_system.G*mass*new_mass/E

            kep.initialize_from_elements(mbin, a, ecc)

            # Binaries should be approaching in order to be picked up
            # by multiples.
            
            kep.advance_to_apastron()
            kep.advance_to_radius(a)

            dr = quantities.AdaptingVectorQuantity()
            dr.extend(kep.get_separation_vector())
            dv = quantities.AdaptingVectorQuantity()
            dv.extend(kep.get_velocity_vector())

            newstar = datamodel.Particles(1)
            newstar.mass = new_mass
            newstar.position = stars[i].position + (1-fac)*dr
            newstar.velocity = stars[i].velocity + (1-fac)*dv

            stars[i].position = stars[i].position - fac*dr
            stars[i].velocity = stars[i].velocity - fac*dv

            newstar.id = companion_base_id + stars[i].id

            stars.add_particles(newstar)
            added_mass += new_mass

            if nbin >= number_of_binaries: break

        kep.stop()

        print('created', nbin, 'binaries')
        sys.stdout.flush()

        stars.mass = stars.mass * total_mass/(total_mass+added_mass)
        number_of_stars += nbin

    # Set dynamical radii (assuming virial equilibrium and standard
    # units).  Note that this choice should be refined, and updated
    # as the system evolves.  Probably the choice of radius should be
    # made entirely in the multiples module.  TODO.  In these units,
    # M = 1 and <v^2> = 0.5, so the mean 90-degree turnaround impact
    # parameter is
    #
    #		b_90 = G (m_1+m_2) / vrel^2
    #		     = 2 <m> / 2<v^2>
    #		     = 2 / N			for equal masses
    #
    # Taking r_i = m_i / 2<v^2> = m_i in virial equilibrium means
    # that, approximately, "contact" means a 90-degree deflection (r_1
    # + r_2 = b_90).  A more conservative choice with r_i less than
    # this value will isolate encounters better, but also place more
    # load on the large-N dynamical module.

    stars.radius = stars.mass.number | nbody_system.length

    time = 0.0 | nbody_system.time
    # print "IDs:", stars.id.number

    print("recentering stars")
    stars.move_to_center()
    sys.stdout.flush()

    #-----------------------------------------------------------------

    if softening_length < 0.0 | nbody_system.length:

        # Use ~interparticle spacing.  Assuming standard units here.  TODO

        eps2 = 0.25*(float(number_of_stars))**(-0.666667) \
			| nbody_system.length**2
    else:
        eps2 = softening_length*softening_length
    print('softening length =', eps2.sqrt())

    gravity.parameters.timestep_parameter = accuracy_parameter
    gravity.parameters.epsilon_squared = eps2
    gravity.parameters.use_gpu = use_gpu
    # gravity.parameters.manage_encounters = manage_encounters

    print('')
    print("adding particles")
    # print stars
    sys.stdout.flush()
    gravity.particles.add_particles(stars)
    gravity.commit_particles()

    print('')
    print("number_of_stars =", number_of_stars)
    sys.stdout.flush()

    # Channel to copy values from the code to the set in memory.
    channel = gravity.particles.new_channel_to(stars)

    stopping_condition = gravity.stopping_conditions.collision_detection
    stopping_condition.enable()

    # Debugging: prevent the multiples code from being called.
    if 0:
        stopping_condition.disable()
        print('stopping condition disabled')
        sys.stdout.flush()

    # -----------------------------------------------------------------
    # Create the coupled code and integrate the system to the desired
    # time, managing interactions internally.

    kep = init_kepler(stars[0], stars[1])
    multiples_code = multiples.Multiples(gravity, new_smalln, kep)

    multiples_code.neighbor_perturbation_limit = 0.1
    #multiples_code.neighbor_distance_factor = 2.0
    multiples_code.neighbor_veto = True
    multiples_code.global_debug = debug_level

    print('')
    print('multiples_code.initial_scale_factor =', \
        multiples_code.initial_scale_factor)
    print('multiples_code.neighbor_perturbation_limit =', \
        multiples_code.neighbor_perturbation_limit)
    print('multiples_code.neighbor_veto =', \
        multiples_code.neighbor_veto)
    print('multiples_code.final_scale_factor =', \
        multiples_code.final_scale_factor)
    print('multiples_code.initial_scatter_factor =', \
        multiples_code.initial_scatter_factor)
    print('multiples_code.final_scatter_factor =', \
        multiples_code.final_scatter_factor)
    print('multiples_code.retain_binary_apocenter =', \
        multiples_code.retain_binary_apocenter)
    print('multiples_code.wide_perturbation_limit =', \
        multiples_code.wide_perturbation_limit)

    # Find initial binaries.
    
    gravity.parameters.zero_step_mode = 1
    print('\nidentifying initial binaries')
    multiples_code.evolve_model(time)
    gravity.parameters.zero_step_mode = 0
    
    pre = "%%% "
    E0,cpu0 = print_log(pre, time, multiples_code)

    print("evolving to time =", end_time, \
            "in steps of", delta_t)

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

        print_log(pre, time, multiples_code, E0, cpu0)
        sys.stdout.flush()

    #-----------------------------------------------------------------

    if not outfile == None:

        # Write data to a file.

        f = open(outfile, 'w')

        #--------------------------------------------------
        # Need to save top-level stellar data and parameters.
        # Need to save multiple data and parameters.

        f.write('%.15g\n'%(time.number))
        for s in multiples_code.stars: write_star(s, f)

        #--------------------------------------------------

        f.close()
        print('wrote file', outfile)

    print('')
    gravity.stop()

def write_star(s, f):
    x,y,z = s.position.number
    vx,vy,vz = s.velocity.number
    f.write('%d %.15g %.15g %.15g %.15g %.15g %.15g %.15g\n' \
		%(s.id, s.mass.number, x, y, z, vx, vy, vz))

if __name__ == '__main__':

    print('\ncommand line:', end=' ')
    for a in sys.argv: print(a, end=' ')
    print('\n')

    infile = None
    outfile = None
    N = 100
    Nbin = 10
    t_end = 10.0 | nbody_system.time
    delta_t = 1.0 | nbody_system.time
    n_workers = 2
    use_gpu = 0
    gpu_worker = 0
    salpeter = 0
    accuracy_parameter = 0.1
    softening_length = 0  | nbody_system.length
    random_seed = 42
    manage_encounters = 1
    debug_level = 1

    try:
        opts, args = getopt.getopt(sys.argv[1:], "a:b:c:d:D:e:f:-F:gGn:s:St:w:")
    except getopt.GetoptError as err:
        print(str(err))
        sys.exit(1)

    for o, a in opts:
        if o == "-a":
            accuracy_parameter = float(a)
        elif o == "-b":
            Nbin = int(a)
        elif o == "-c":
            manage_encounters = int(a)
        elif o == "-d":
            delta_t = float(a) | nbody_system.time 
        elif o == "-D":
            debug_level = int(a)
        elif o == "-e":
            softening_length = float(a) | nbody_system.length
        elif o == "-f":
            infile = a
        elif o == "-F":
            outfile = a
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
            print("unexpected argument", o)

    assert is_mpd_running()
    run_ph4(infile, outfile,
            N, Nbin, t_end, delta_t, n_workers,
            use_gpu, gpu_worker,
            salpeter, accuracy_parameter, softening_length,
            manage_encounters, random_seed, debug_level)
