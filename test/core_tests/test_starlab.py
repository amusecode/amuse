from amuse.support.data import core
from amuse.support.io import starlab
from amuse.support.units import units

import unittest
import pylab as pl
import os.path

class Test(unittest.TestCase):

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
        
        I = starlab.ParticlesFromDyn(os.path.join(directory, 'test_subsub.dyn'))

        All = I.Particles
        print All.mass

        self.assertEquals(len(All), 7)
        self.assertEquals(len(All[0].descendents()),6)
        self.assertEquals(All[0].children().mass.value_in(units.kg)[0], 8.0)
        self.assertEquals(All[1].children().mass.value_in(units.kg)[0], 4.0)
        self.assertEquals(All[5].children().mass.value_in(units.kg)[0], 2.0)
