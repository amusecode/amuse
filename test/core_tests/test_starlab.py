from amuse.support.data import core
from amuse.support.io import starlab
from amuse.support.units import units
from amuse.support.units import nbody_system

from amuse.test import amusetest
import os.path

class Test(amusetest.TestCase):

    def test1(self):
        """
                                                                                          
                +---------------------------------------------------------------------+   
                |                      Particle tree of test_subub.dyn                |   
                |                                                                     |   
                |                           0 16kg x=0                                |   
                |                           ,-  .                                     |   
                |                        ,-'     `._                                  |   
                |                     _,'           `.                                |   
                |                   ,'                `-.                             |   
                |                 1 8kg, x=-10        3 8kg, x=10                     |   
                |                ,-'                   - ._                           |   
                |              ,'                    ,'    `.                         |   
                |            ,'                    ,'        `._                      |   
                |          ,'                    ,'             `.                    |   
                |       2  4kg, x=-15          4 4kg, x=5      5 4kg, x=15            |   
                |                                                   .                 |   
                |                                                    `._              |   
                |                                                       `.            |   
                |                                                         `           |   
                |                                                       6 2kg, x=17   |   
                |                                                                     |   
                +---------------------------------------------------------------------+   

        """
        directory = os.path.dirname(__file__)
        convert_nbody = nbody_system.nbody_to_si(1|units.kg, 1|units.m)

        I = starlab.ParticlesFromDyn(os.path.join(directory, 'test_subsub.dyn'))

        All = I.Particles

        self.assertEquals(len(All), 7)
        self.assertEquals(len(All[0].descendents()),6)
        self.assertEquals(All[0].children().mass.value_in(nbody_system.mass)[0], 8.0)
        self.assertEquals(All[1].children().mass.value_in(nbody_system.mass)[0], 4.0)
        self.assertEquals(All[5].children().mass.value_in(nbody_system.mass)[0], 2.0)

    def test2(self):
        directory = os.path.dirname(__file__)
        convert_nbody = nbody_system.nbody_to_si(1|units.kg, 1|units.m)

        I = starlab.ParticlesFromDyn(os.path.join(directory, 'test_subsub.dyn'), convert_nbody)

        All = I.Particles

        self.assertEquals(len(All), 7)
        self.assertEquals(len(All[0].descendents()),6)
        self.assertEquals(All[0].children().mass.value_in(units.kg)[0], 8.0)
        self.assertEquals(All[1].children().mass.value_in(units.kg)[0], 4.0)
        self.assertEquals(All[5].children().mass.value_in(units.kg)[0], 2.0)
        
