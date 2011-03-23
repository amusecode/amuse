from amuse.support.data import core
from amuse.support.io import starlab
from amuse.support import io
from amuse.support.units import units
from amuse.support.units import nbody_system

from amuse.test import amusetest
import os.path
import StringIO
import numpy

plummer_scaled_content = """(Particle
  N = 5
(Log
 ===>  Wed Mar 23 10:36:19 2011
       Starlab 4.4.4 (user vanelteren) : makeplummer -n 5
       random number generator seed = 1300872979
  initial_mass = 1.93791189968789
 ===>  Wed Mar 23 10:36:19 2011
       Starlab 4.4.4 (user vanelteren) : makemass -l 0.2 -u 20 -x -2.35
       random number generator seed = 1300872979
       Power_Law mass function, total mass =     1.94 Solar
)Log
(Dynamics
  system_time  =  0
  m  =  1.93791189968789257
  r  =  0 0 0
  v  =  0 0 0
  com_time = 0
  com_pos = 0 0 0
  com_vel = 0 0 0
  total_energy = -0.25
)Dynamics
(Hydro
)Hydro
(Star
  mass_scale     =  0.516019329960796136
  size_scale     =  -2.25500000000000001e-08
  time_scale     =  -1
)Star
(Particle
  N = 1
(Log
)Log
(Dynamics
  m  =  0.266730587350341442
  r  =  -0.934040458729136547 -0.695803261872860679 0.564081767628105579
  v  =  0.293044501505944688 0.0259404966497079996 0.196831834670057881
)Dynamics
(Hydro
)Hydro
(Star
)Star
)Particle
(Particle
  N = 1
(Log
)Log
(Dynamics
  m  =  0.32415747586178062
  r  =  -0.177487175943360831 0.205223807753114523 -0.191956558941283828
  v  =  -0.178395285089674116 -0.0357730795197053753 0.376231470164175796
)Dynamics
(Hydro
)Hydro
(Star
)Star
)Particle
(Particle
  N = 1
(Log
)Log
(Dynamics
  m  =  0.567097501003086424
  r  =  0.607757256863783235 0.120278701131815768 0.338645014325028748
  v  =  -0.26638454687085783 0.291986511517820568 0.672417896548585303
)Dynamics
(Hydro
)Hydro
(Star
)Star
)Particle
(Particle
  N = 1
(Log
)Log
(Dynamics
  m  =  0.379463033750591316
  r  =  0.247202300114546969 0.086741731241469916 -0.164892802949227257
  v  =  -0.0888917818811222754 0.286288450446576914 -1.03818564484543829
)Dynamics
(Hydro
)Hydro
(Star
)Star
)Particle
(Particle
  N = 1
(Log
)Log
(Dynamics
  m  =  0.400463301722092879
  r  =  0.256568077694167229 0.283559021746460527 -0.545877420062623298
  v  =  0.240627112335709492 -0.568442379094400096 -0.207295556537380576
)Dynamics
(Hydro
)Hydro
(Star
)Star
)Particle
)Particle
"""

class Test(amusetest.TestCase):

    def test1(self):
        """test_starlab.test1
                                                                                          
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
        
    def test3(self):
        directory = os.path.dirname(__file__)
        set = io.read_set_from_file(os.path.join(directory, 'plummer.dyn'), 'starlab')
        self.assertEquals(len(set), 10)
        self.assertAlmostRelativeEquals(set.mass, 0.1 | nbody_system.mass)
        
    def test4(self):
        set = starlab.StarlabFileFormatProcessor().load_string(plummer_scaled_content)
        self.assertEquals(len(set), 5)
        self.assertTrue(numpy.all(set.mass > 0.2 |units.MSun))
        self.assertTrue(numpy.all(set.mass < 1.1 |units.MSun))
        self.assertTrue(numpy.all(set.x > -1 | units.parsec))
        self.assertTrue(numpy.all(set.x < 1 | units.parsec))
        self.assertTrue(numpy.all(set.y > -1 | units.parsec))
        self.assertTrue(numpy.all(set.y < 1 | units.parsec))
        self.assertTrue(numpy.all(set.z > -1 | units.parsec))
        self.assertTrue(numpy.all(set.z < 1 | units.parsec))
        
