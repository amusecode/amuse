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

with_stellar_structure = """(Particle
  N = 2
(Log
 ===>  Wed Mar 23 12:28:21 2011
       Starlab 4.4.4 (user vanelteren) : makeplummer -n 2
       random number generator seed = 1300879701
  initial_mass = 1
 ===>  Wed Mar 23 12:28:21 2011
       Starlab 4.4.4 (user vanelteren) : makemass -f 1 -x -2.0 -l 0.1 -u 20
       random number generator seed = 1300879701
       Power_Law mass function, total mass =     0.22 Solar
 ===>  Wed Mar 23 12:28:21 2011
       Starlab 4.4.4 (user vanelteren) : add_star -Q 0.5 -R 5
 ===>  Wed Mar 23 12:28:21 2011
       Starlab 4.4.4 (user vanelteren) : scale -s
  initial_total_energy = -0.25
  initial_rvirial = 1
)Log
(Dynamics
  system_time  =  0
  m  =  1.000000000000002
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
  mass_scale     =  4.62223018200823255
  size_scale     =  4.51000000000000035e-09
  time_scale     =  0.0028016566610805007
)Star
(Particle
  N = 1
(Log
)Log
(Dynamics
  m  =  0.483839787917132669
  r  =  -0.209118735762131774 -0.0880146969484976449 0.122429686361064466
  v  =  0.435701422371714053 -0.0578440884891995299 0.583282312300410277
)Dynamics
(Hydro
)Hydro
(Star
  Type   =  main_sequence
  T_cur  =  0
  M_rel  =  0.104676696933106134
  M_env  =  0.0946766969331061387
  M_core =  0.0100000000000000002
  T_eff  =  3011.01000587455155
  L_eff  =  0.00122847524117014736
)Star
)Particle
(Particle
  N = 1
(Log
)Log
(Dynamics
  m  =  0.516160212082869219
  r  =  0.196024339714903045 0.0825034772310478115 -0.114763501907016688
  v  =  -0.408419089384746692 0.0542220629403740648 -0.546758900962974193
)Dynamics
(Hydro
)Hydro
(Star
  Type   =  main_sequence
  T_cur  =  0
  M_rel  =  0.111669084350665276
  M_env  =  0.101669084350665281
  M_core =  0.0100000000000000002
  T_eff  =  3051.91244441666595
  L_eff  =  0.00147739782384390963
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
        print set.mass.as_quantity_in(units.MSun)
        self.assertTrue(numpy.all(set.mass > 0.2 |units.MSun))
        self.assertTrue(numpy.all(set.mass < 1.1 |units.MSun))
        self.assertTrue(numpy.all(set.x > -1 | units.parsec))
        self.assertTrue(numpy.all(set.x < 1 | units.parsec))
        self.assertTrue(numpy.all(set.y > -1 | units.parsec))
        self.assertTrue(numpy.all(set.y < 1 | units.parsec))
        self.assertTrue(numpy.all(set.z > -1 | units.parsec))
        self.assertTrue(numpy.all(set.z < 1 | units.parsec))
        
    def test5(self):
        set = starlab.StarlabFileFormatProcessor().load_string(with_stellar_structure)
        self.assertEquals(len(set), 2)
        self.assertAlmostRelativeEquals(set[0].envelope_mass, 0.0946766969331061387 | units.MSun)
        self.assertAlmostRelativeEquals(set[0].core_mass, 0.0100000000000000002 | units.MSun)
        self.assertAlmostRelativeEquals(set[0].relative_mass,  set[0].mass, 10)
        self.assertAlmostRelativeEquals(set[0].effective_temperature,  3011.01000587455155 | units.K)
        self.assertAlmostRelativeEquals(set[0].effective_luminocity,  0.00122847524117014736 | units.LSun)
        self.assertAlmostRelativeEquals(set[0].stellar_type, units.stellar_type("Main Sequence star"))
        self.assertAlmostRelativeEquals(set[1].relative_mass,  set[1].mass, 10)
        
    def test6(self):
        directory = os.path.dirname(__file__)
        set = io.read_set_from_file(os.path.join(directory, 'evolved.dyn'), 'starlab')
        self.assertEquals(len(set), 20)
        
        self.assertAlmostRelativeEquals(set.time, set.age, 4) #only accurate to 4, starlab stellar evolution time scale?
        
        self.assertAlmostRelativeEquals(set[0].velocity, [177.579717905, 38.5027308364, -35.8571344243] | units.km / units.hour, 8)
        self.assertAlmostRelativeEquals(set[0].acceleration, [-0.000648471729782, 0.000309476774701, -0.000356623346185] | units.parsec / (units.Myr ** 2), 8)
        self.assertAlmostRelativeEquals(set.unconverted_set()[0].specific_potential, -0.32735384622167929 | nbody_system.potential)
        #select the main sequence star, the dwarf masses don't match
        main_sequence_stars = set.select(lambda stellar_type : stellar_type ==  units.stellar_type("Main Sequence star"), ["stellar_type"])
        self.assertAlmostRelativeEquals(main_sequence_stars.mass, main_sequence_stars.relative_mass, 10)
        carbon_dwarfs = set.select(lambda stellar_type : stellar_type ==   units.stellar_type("Carbon/Oxygen White Dwarf") , ["stellar_type"])
        self.assertAlmostRelativeEquals(carbon_dwarfs.mass, carbon_dwarfs.core_mass, 10)\
        
        
        
