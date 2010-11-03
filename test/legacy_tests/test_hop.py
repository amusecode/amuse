import os
import sys
import numpy
from random import random as rnd

from amuse.test import amusetest
from amuse.test.amusetest import get_path_to_results

from amuse.support.units import units, nbody_system

from amuse.legacy.fi import interface as interface
from amuse.support.legacy.core import is_mpd_running
from amuse.ext.evrard_test import regular_grid_unit_cube
from amuse.ext.molecular_cloud import molecular_cloud
from amuse.ext.evrard_test import body_centered_grid_unit_cube
from amuse.ext import plummer

from amuse.legacy.hop.interface import HopInterface

class TestHop(amusetest.TestCase):
    def test1(self):
        print "First test: adding particles, setting and getting."
        hop = HopInterface()
        n, err = hop.get_number_of_particles()
        for i in range(6):
            id, err = hop.new_particle(0, i*i, 0, 0)
            n, err = hop.get_number_of_particles()
            self.assertEquals(n, i+1)
        for i in range(6):
            x, y, z, err = hop.get_position(i)
            self.assertEquals(x, i*i)
            hop.set_position(i, x, i*i, 0)
#        for i in range(6):
#            x, y, z = hop.get_position(i)
#            self.assertEquals(x, y)
            
    def xtest2(self):  #test broken for unknown reason, needs author to fix
        print "Second test: calculating densities."
        hop = run_cloud(1000)
        n, err = hop.get_number_of_particles()
        self.assertEquals(n, 909)
        ds = {0:0.026384295895695686, 1:0.022254005074501038, 2:0.059249848127365112}
        for method in [0,1,2]:
            hop.set_density_method(method)
            hop.calculate_densities()
            d, err = hop.get_density(0)
            self.assertEquals(d,ds[method])
    
    def xtest3(self): #test broken for unknown reason, needs author to fix
        print "Third test: densest neighbors and groups."
        hop = run_cloud(1000)
        hop.calculate_densities()
        hop.do_hop()
        n, err = hop.get_densest_neighbor(0)
        self.assertEquals(n, 55)
        n, err = hop.get_group_id(0)
        self.assertEquals(n, 3)
        n, err = hop.get_densest_particle_in_group(3)
        self.assertEquals(n, 80)
        n, err = hop.get_number_of_particles_in_group(0)
        self.assertEquals(n, 223)
        n, err = hop.get_average_boundary_density_of_groups(0,1)
        self.assertEquals(n, 1.5265605449676514)
        n, err = hop.get_number_of_groups()
        self.assertEquals(n, 9)



def run_cloud(n):
    ''' partially taken from test_molecular_cloud.py '''
    cloud=molecular_cloud(32,-4.,n,base_grid=regular_grid_unit_cube,seed=299792)
    mass,x,y,z,vx,vy,vz,u=cloud.new_model()
    smooth=numpy.zeros_like(mass)

    nb = interface.FiInterface(redirection="none")
    nb.initialize_code()

    nb.set_stepout(99999)
    nb.set_steplog(99999)
    nb.set_use_hydro(1)
    nb.set_radiate(0)
    nb.set_dtime(0.05)
    nb.set_gdgop(1)
    nb.set_uentropy(0)
    nb.set_isotherm(1)
    nb.set_gamma(1.0)
    nb.set_verbosity(0)
    nb.set_unitl_in_kpc(0.01)
    nb.set_unitm_in_msun(10000.)
      
    ids,error = nb.new_sph_particle(mass,smooth,x,y,z,vx,vy,vz,u)
    if filter(lambda x: x != 0, error) != []: raise Exception

    nb.commit_particles()

    if hasattr(nb,"viewer"):
        nb.viewer()

    dt=0.05
    tnow=0.
    nb.synchronize_model()
    time,Ek,Ep,Eth=[],[],[],[]
    time.append(tnow)
    e,ret=nb.get_kinetic_energy()
    Ek.append(e)
    e,ret=nb.get_potential_energy()
    Ep.append(e)
    e,ret=nb.get_thermal_energy()
    Eth.append(e)

    while tnow<.8:
        tnow=tnow+dt
        nb.evolve(tnow)
        nb.synchronize_model()
        tnow,err=nb.get_time()
        time.append(tnow)
        e,ret=nb.get_kinetic_energy()
        Ek.append(e)
        e,ret=nb.get_potential_energy()
        Ep.append(e)
        e,ret=nb.get_thermal_energy()
        Eth.append(e)

    m,h,x,y,z,vx,vy,vz,err=nb.get_state(ids)
    nb.stop()
    hop = HopInterface()
    for i in range(len(x)):
        hop.new_particle(x[i], y[i], z[i])
    return hop
