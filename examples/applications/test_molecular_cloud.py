import os
import sys
import numpy

from amuse.community.fi import interface as interface
from amuse.ext.evrard_test import regular_grid_unit_cube
from amuse.ext.molecular_cloud import molecular_cloud
from amuse.ext.evrard_test import body_centered_grid_unit_cube

from amuse.test.amusetest import get_path_to_results
from amuse.support.codes.core import is_mpd_running

try:
    from matplotlib import pyplot
    HAS_MATPLOTLIB = True
except ImportError:
    HAS_MATPLOTLIB = False



def energy_plot(time,ek,ep,eth):
    if not HAS_MATPLOTLIB:
        return
      
    pyplot.figure(figsize = (5, 5))
    pyplot.xlabel(r'time')
    pyplot.ylabel(r'energy')
    pyplot.plot(time,ek)
    pyplot.plot(time,ep)
    pyplot.plot(time,eth)
    pyplot.plot(time,ek+ep+eth)
    test_results_path = get_path_to_results()
    pyplot.savefig(os.path.join(test_results_path, "test.png"))

def run_cloud(x):
    cloud=molecular_cloud(32,-4.,x,base_grid=regular_grid_unit_cube)
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
        nb.evolve_model(tnow)
        nb.synchronize_model()
        tnow,err=nb.get_time()
        time.append(tnow)
        e,ret=nb.get_kinetic_energy()
        Ek.append(e)
        e,ret=nb.get_potential_energy()
        Ep.append(e)
        e,ret=nb.get_thermal_energy()
        Eth.append(e)

    nb.stop()
   
    time=numpy.array(time)
    Ek=numpy.array(Ek)
    Ep=numpy.array(Ep)
    Eth=numpy.array(Eth)
    energy_plot(time,Ek,Ep,Eth)
  
def test_cloud():
    if not is_mpd_running():
        return
    run_cloud(1000)

if __name__=="__main__":
    run_cloud(1000)

