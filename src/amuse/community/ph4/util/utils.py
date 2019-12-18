import collections
import getopt
import numpy
import os
import random
import sys
import pickle
import numpy as np

import unittest
from time import clock

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

from amuse import io

SMALLN = None
def new_smalln():
    SMALLN.reset()
    return SMALLN
    
def init_smalln():
    global SMALLN
    sys.stdout.flush()
    SMALLN = SmallN()#redirection='none')
    sys.stdout.flush()
    SMALLN.parameters.timestep_parameter = 0.1
    #SMALLN.parameters.cm_index = 2001		# don't set this here!!
    sys.stdout.flush()

def stop_smalln():
    global SMALLN
    SMALLN.stop()

def init_kepler(star1, star2):        
    try:
        star1.mass.value_in(units.kg) # see if SI units, throw exception if not
        unit_converter \
            = nbody_system.nbody_to_si(star1.mass + star2.mass,
                                       (star2.position-star1.position).length())
    except Exception as ex:
        unit_converter = None
        
    kep = Kepler(unit_converter) # , redirection = "none")
    kep.initialize_code()

    return kep

def print_binaries(pre, stars, kT, limit):	# print all bound pairs with
						# binding energy >= limit*kT
    Emax = -(limit*kT).number
    m = np.array(stars.mass.number)
    pos = np.array(stars.position.number)
    vel = np.array(stars.velocity.number)
    n = len(stars)
    dx = np.zeros((n,3))
    dv = np.zeros((n,3))
    tiny = 1.e-12

    Ntop = 0
    Etop = 0.0
    print('bound top-level pairs:')
    for i in range(n):
        dx[:,:] = pos[:,:] - pos[i,:]
        dv[:,:] = vel[:,:] - vel[i,:]
        dr2 = (dx**2).sum(axis=1) + tiny
        dr2i = 1./dr2
        dr2i[i] = 0.0
        mi = m[i]
        dpot = -(m+mi)*np.sqrt(dr2i)
        dpot[:i] = 0.0
        dkin = 0.5*(dv**2).sum(axis=1)
        mu = mi*m/(mi+m)
        de = mu*(dpot + dkin)
        w = np.where(de < Emax)
        jw = np.arange(n)[w]
        if len(jw) > 0:
            for j in jw:
                print('   ', stars.id[i], stars.id[j], de[j]/kT.number, 'kT')
                Ntop += 1
                Etop += de[j]

    if Ntop == 0:
        print('(none)')

    print('')
    print(pre, 'Ntop =', Ntop)
    print(pre, 'Etop =', Etop)
    print(pre, 'Etop/kT =', Etop/kT.number)

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
    Edel = Eext + Eint + Eerr
    Ecor = Etot - Edel
    if E0 == 0 | nbody_system.energy: E0 = Ecor
    Rvir = -0.5*M*M/U
    Q = -T/U
    com = pa.center_of_mass(gravity.particles)
    comv = pa.center_of_mass_velocity(gravity.particles)

    # Hop complains if we have too few particles.

    dcen = com
    if N >= 100:
        dcen,rcore,rhocore \
            = pa.densitycentre_coreradius_coredens(gravity.particles)
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
    print(pre+"dE=", Ecor.number - E0.number)
    print(pre+"dE/E=", Ecor/E0 - 1)
    print(pre+"Rvir=", Rvir.number)
    print(pre+"Qvir=", Q)
    cmx,cmy,cmz = com
    print(pre+"cmpos[3]= %.8f %.8f %.8f" % (cmx.number, cmy.number, cmz.number))
    cmx,cmy,cmz = comv
    print(pre+"cmvel[3]= %.8f %.8f %.8f" % (cmx.number, cmy.number, cmz.number))

    if N >= 100:
        cmx,cmy,cmz = dcen
        print(pre+"dcpos[3]= %.8f %.8f %.8f" \
            		% (cmx.number, cmy.number, cmz.number))
        print(pre+"Rcore=", rcore.number)
        print(pre+"Mcore=", (rhocore*rcore**3).number)	# fake...
        print(pre+"Mlagr[9]=", end=' ')
        for m in mf: print("%.4f" % (m), end=' ')
        print('')
        print(pre+"Rlagr[9]=", end=' ')
        for r in lagr.number: print("%.8f" % (r), end=' ')
        print('')

    kT = 2*T/(3.*N)	# 3/2 N kT = total KE, by definition
    kT0 = -2*E0/(3.*N)
    print(pre+"kT= %.8f" % (kT.number))
    Nmul,Nbin,Emul = gravity.print_multiples2(pre, kT, dcen)
    print(pre+"Nmul=", Nmul)
    print(pre+"Nbin=", Nbin)
    print(pre+"Emul= %.5f" % (Emul.number))
    print(pre+"Emul2= %.5f" % (Emul2.number))
    print(pre+"Emul/kT= %.5f" % (Emul/kT))
    print(pre+"Emul/kT0= %.5f" % (Emul/kT0))
    print(pre+"Emul/E= %.5f" % (Emul/Etot))
    print('')

    stars = gravity.particles.copy()
    schannel = gravity.particles.new_channel_to(stars)
    schannel.copy_attribute("index_in_code", "id")
    print_binaries(pre, stars, kT, 0.1)

    sys.stdout.flush()
    return Ecor,cpu


# Save code from Brewer...

def write_state_to_file(time, stars_python, gravity_code, multiples_code,
                        write_file=None,
                        delta_t=1.0|nbody_system.time,
                        E0_log=0.0, cpu0_log=0.0):

    # Function for saving everything.  If you need to track more than
    # is normally saved, should be safe to tack it onto params.

    if write_file is not None:

        particles = gravity_code.particles.copy()
        write_channel = gravity_code.particles.new_channel_to(particles)
        write_channel.copy_attribute("index_in_code", "id")

        params = {'timestep_parameter': \
                      gravity_code.parameters.timestep_parameter,
                  'epsilon_squared': gravity_code.parameters.epsilon_squared,
                  'neighbor_veto': multiples_code.neighbor_veto,
                  'multiples_external_tidal_correction': \
                      multiples_code.multiples_external_tidal_correction,
                  'multiples_integration_energy_error': \
                      multiples_code.multiples_integration_energy_error,
                  'multiples_internal_tidal_correction': \
                      multiples_code.multiples_internal_tidal_correction,
                  'model_time': multiples_code.model_time,
                  'delta_t': delta_t,
                  'EZero': E0_log,
                  'CPUZero': cpu0_log,
                  'root_index': multiples.root_index
        }
        
        for root, tree in multiples_code.root_to_tree.items():
            root_in_particles = root.as_particle_in_set(particles)
            subset = tree.get_tree_subset().copy()
            if root_in_particles is not None:
                root_in_particles.components = subset

        io.write_set_to_file(particles,
                             write_file+".stars.hdf5",
                             'amuse',
                             version='2.0',
                             append_to_file=False,
                             copy_history=False,
                             close_file=True)

        io.write_set_to_file(stars_python,
                             write_file+".stars_python.hdf5",
                             'amuse',
                             version='2.0',
                             append_to_file=False,
                             copy_history=False,
                             close_file=True)

        config = {'time': time,
                  'py_seed': pickle.dumps(random.getstate()),
                  'numpy_seed': pickle.dumps(numpy.random.get_state()),
              }
        with open(write_file + ".conf", "wb") as f:
            pickle.dump(config, f)
            f.close()
            with open(write_file + ".params", "wb") as f:
                pickle.dump(params, f)
                f.close()

        print("\nwrote state to file ", write_file, 'at time', time)


# Restore code from Brewer.

def read_state_from_file(restart_file, gravity_code, kep, MT=0):

    # Function to load from file.  If you change params in
    # write_state_to_file, make sure you match the changes here.

    stars = io.read_set_from_file(restart_file+".stars.hdf5",
                                  'amuse', version='2.0',
                                  close_file=True).copy()
    stars_python = io.read_set_from_file(
        		restart_file+".stars_python.hdf5",
        		'amuse',version='2.0',
			close_file=True).copy()
    with open(restart_file + ".params", "rb") as f:
        params = pickle.load(f)
        f.close()
    
    if MT == 0:
        root_to_tree = {}
        for root in stars:
            if hasattr(root, 'components') and not root.components is None:
                root_to_tree[root] \
                    = datamodel.trees.BinaryTreeOnParticle(root.components[0])
    else:
        root_to_tree = MT

    gravity_code.parameters.timestep_parameter = params['timestep_parameter']
    gravity_code.parameters.epsilon_squared = params['epsilon_squared']

    gravity_code.particles.add_particles(stars)
    gravity_code.commit_particles()		# sets time steps

    multiples_code = multiples.Multiples(gravity_code, new_smalln, kep)
    multiples_code.neighbor_veto = params['neighbor_veto']
    multiples_code.multiples_external_tidal_correction \
	= params['multiples_external_tidal_correction']
    multiples_code.multiples_integration_energy_error \
	= params['multiples_integration_energy_error']
    multiples_code.multiples_internal_tidal_correction \
	= params['multiples_internal_tidal_correction']
    multiples.root_index = params['root_index']
    multiples_code.root_to_tree = root_to_tree

    print("\nread state from file ", restart_file, \
          'at time', params['model_time'])

    return stars_python, params['model_time'], params['delta_t'], \
		params['EZero'], params['CPUZero'], multiples_code
