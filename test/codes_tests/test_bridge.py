import numpy 
from amuse.test.amusetest import TestWithMPI

from amuse.support.units import nbody_system
from amuse.support.units import units

from amuse.community.fi.interface import Fi
from amuse.community.hermite0.interface import Hermite
from amuse.community.phiGRAPE.interface import PhiGRAPE
from amuse.community.bhtree.interface import BHTree
from amuse.ext.bridge import bridge
from amuse.ext.kingmodel import new_king_model

def sys_from_parts(base_class,parts,converter,eps=None):
    interface=base_class(converter)
    interface.initialize_code()
    if eps is not None:
        interface.parameters.epsilon_squared = eps**2 
    interface.particles.add_particles(parts)
    return interface

class testBridge(TestWithMPI):

    def test1(self):
        convert = nbody_system.nbody_to_si(1.e5 | units.MSun, 1.0 | units.parsec)
        eps=2.e-4 | nbody_system.length

        test_class=PhiGRAPE
        number_of_particles = 50
        stars = new_king_model(number_of_particles, W0=7, convert_nbody=convert)
        stars.radius = 0.0 | units.RSun

        cluster=test_class(convert)
        cluster.parameters.epsilon_squared = eps**2 
        cluster.particles.add_particles(stars)
        cluster.synchronize_model()

        Ep1=convert.to_nbody(cluster.potential_energy).number
        Ek1=convert.to_nbody(cluster.kinetic_energy).number

        parts=cluster.particles.copy()
        parts1=parts.select_array(lambda x: (x > 0 | units.m), ['x'] )
        parts2=parts.select_array(lambda x: (x < 0 | units.m), ['x'] )
        cluster1=sys_from_parts(test_class, parts1, convert, eps)
        cluster2=sys_from_parts(test_class, parts2, convert, eps)

        cluster1.synchronize_model()
        cluster2.synchronize_model()

        bridgesys=bridge()
        bridgesys.add_system(cluster1, (cluster2,) )
        bridgesys.add_system(cluster2, (cluster1,) )

        Ep2=convert.to_nbody(bridgesys.potential_energy).number
        Ek2=convert.to_nbody(bridgesys.kinetic_energy).number
        self.assertAlmostEqual(Ek1,Ek2,12)
        self.assertAlmostEqual(Ep1,Ep2,12)
  
  
  
