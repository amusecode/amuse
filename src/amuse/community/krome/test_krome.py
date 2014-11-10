import os.path
import numpy
from amuse.test.amusetest import TestWithMPI

from amuse.community.krome.interface import KromeInterface,Krome,solar_abundances
from amuse.units import units
from amuse.datamodel import Particles

from amuse.io import read_set_from_file

default_options = dict(redirection="none")
#default_options=dict(debugger="gdb")

class TestKromeInterface(TestWithMPI):

    def test1(self):
        print "Test 1: initialization"

        instance = KromeInterface(**default_options)
        self.assertEqual(0, instance.initialize_code())
        self.assertEqual(0, instance.commit_parameters())
        self.assertEqual(0, instance.cleanup_code())
        instance.stop()

    def test2(self):
        print "Test 1: add particle, get state"

        instance = KromeInterface(**default_options)
        self.assertEqual(0, instance.initialize_code())
        self.assertEqual(0, instance.commit_parameters())

        dens=1.e5
        t=500.
        ion=1.e-11
        id,err=instance.new_particle(dens,t,ion)

        self.assertEqual(err,0)
      
        self.assertEqual(instance.commit_particles(),0)

        dens_,t_,ion_,err=instance.get_state(id)

        self.assertEqual(err,0)

        self.assertEqual(dens_,dens)
        self.assertEqual(t_,t)
        self.assertEqual(ion_,ion)

        self.assertEqual(0, instance.cleanup_code())

        instance.stop()
        
    def test3(self):
        print "Test 1: add 2 particles, get state"

        instance = KromeInterface(**default_options)
        self.assertEqual(0, instance.initialize_code())
        self.assertEqual(0, instance.commit_parameters())

        dens=[1.e5,2.e5]
        t=[500.,550]
        ion=[1.e-11,2.e-11]
        id,err=instance.new_particle(dens,t,ion)

        self.assertEqual(err,0)
      
        self.assertEqual(instance.commit_particles(),0)

        for i in range(2):
          dens_,t_,ion_,err=instance.get_state(id[i])

          self.assertEqual(err,0)

          self.assertEqual(dens_,dens[i])
          self.assertEqual(t_,t[i])
          self.assertEqual(ion_,ion[i])

        self.assertEqual(0, instance.cleanup_code())

        instance.stop()

    def test4(self):
        print "Test 1: add 100 particles, get state"

        instance = KromeInterface(**default_options)
        self.assertEqual(0, instance.initialize_code())
        self.assertEqual(0, instance.commit_parameters())

        dens=1.e5*numpy.random.random(100)
        t=500.*numpy.random.random(100)
        ion=1.e-11*numpy.random.random(100)
        id,err=instance.new_particle(dens,t,ion)

        self.assertEqual(err,0)
      
        self.assertEqual(instance.commit_particles(),0)

        for i in range(100):
          dens_,t_,ion_,err=instance.get_state(id[i])

          self.assertEqual(err,0)

          self.assertEqual(dens_,dens[i])
          self.assertEqual(t_,t[i])
          self.assertEqual(ion_,ion[i])

        self.assertEqual(0, instance.cleanup_code())

        instance.stop()

    def test5(self):
        print "Test 5: can we get species?"

        instance = KromeInterface(**default_options)
        
        first,last,err=instance.get_firstlast_abundance()
        self.assertEqual(err,0)
        self.assertTrue(last-first > 0)
        
        for i in range(first,last+1):
          name,err=instance.get_name_of_species(i)
          print name
          self.assertEqual(err,0)
          index,err=instance.get_index_of_species(name)
          self.assertEqual(i,index)

        instance.stop()

    def test6(self):
        print "Test 6: add 100 particles, remove particles"

        instance = KromeInterface(**default_options)
        self.assertEqual(0, instance.initialize_code())
        self.assertEqual(0, instance.commit_parameters())

        dens=1.e5*numpy.random.random(100)
        t=500.*numpy.random.random(100)
        ion=1.e-11*numpy.random.random(100)
        ids,err=instance.new_particle(dens,t,ion)

        self.assertEqual(err,0)
      
        self.assertEqual(instance.commit_particles(),0)

        for i in ids[:10]:
          instance.delete_particle(i)

        instance.recommit_particles()

        for i in range(10,100):
          dens_,t_,ion_,err=instance.get_state(ids[i])

          self.assertEqual(err,0)

          self.assertEqual(dens_,dens[i])
          self.assertEqual(t_,t[i])
          self.assertEqual(ion_,ion[i])

        dens_,t_,ion_,err=instance.get_state(ids[0])

        self.assertEqual(err,-1)

        self.assertEqual(0, instance.cleanup_code())

        instance.stop()

    def test7(self):
        print "Test 1: add particle, set abundances"

        instance = KromeInterface(**default_options)
        self.assertEqual(0, instance.initialize_code())
        self.assertEqual(0, instance.commit_parameters())

        dens=1.e5
        t=500.
        ion=1.e-11
        id,err=instance.new_particle(dens,t,ion)
 
        instance.commit_particles()

        first,last,err=instance.get_firstlast_abundance()
        for i in range(first,last+1):
          x,err=instance.get_abundance(id,i)
          self.assertEqual(x,0.)
          self.assertEqual(err,0)
        x,err=instance.get_abundance(id,last+1)
        self.assertEqual(err,-1)
        
        for s in ["H","HE","C","SI","O"]:
          x=solar_abundances[s]
          aid,err=instance.get_index_of_species(s)          
          instance.set_abundance(id,aid,x)

        for s in ["H","HE","C","SI","O"]:
          x=solar_abundances[s]
          aid,err=instance.get_index_of_species(s)          
          xx,err=instance.get_abundance(id,aid)
          self.assertEqual(x,xx)
          self.assertEqual(err,0)
          
    def test8(self):
        print "evolve test"

        instance = KromeInterface(**default_options)
        self.assertEqual(0, instance.initialize_code())
        self.assertEqual(0, instance.commit_parameters())

        dens=1.e2
        t=100.
        ion=2.e-17
        id,err=instance.new_particle(dens,t,ion)
         
        first,last,err=instance.get_firstlast_abundance()
        for i in range(first,last+1):
          err=instance.set_abundance(id,i,1.e-40)         
         
        for s in ["H","HE","C","SI","O"]:
          x=solar_abundances[s]
          aid,err=instance.get_index_of_species(s)          
          instance.set_abundance(id,aid,x)

        instance.commit_particles()

        yr=365*24*3600.
        err=instance.evolve_model(10000.*yr)
        self.assertEqual(err,0)
        time,err=instance.get_time()
        self.assertEqual(err,0)
        self.assertEqual(time,10000.*yr)

        first,last,err=instance.get_firstlast_abundance()
        for i in range(first,last+1):
          x,err=instance.get_abundance(id,i)
          self.assertEqual(err,0)
          name,err=instance.get_name_of_species(i)
          print i,name,x
        
    def test9(self):
        print "evolve test 2"

        instance = KromeInterface(**default_options)
        self.assertEqual(0, instance.initialize_code())
        self.assertEqual(0, instance.commit_parameters())

        dens=1.e2
        t=100.
        ion=2.e-17
        id,err=instance.new_particle(dens,t,ion)
         
        first,last,err=instance.get_firstlast_abundance()
        for i in range(first,last+1):
          err=instance.set_abundance(id,i,1.e-40)         
         
        for s in ["H","HE","C","SI","O"]:
          x=solar_abundances[s]
          aid,err=instance.get_index_of_species(s)          
          instance.set_abundance(id,aid,x)

        aid,err=instance.get_index_of_species("C")          
        instance.set_abundance(id,aid,1.e-40)
        
        aid,err=instance.get_index_of_species("C+")          
        instance.set_abundance(id,aid,solar_abundances["C"])

        aid,err=instance.get_index_of_species("H2")          
        instance.set_abundance(id,aid,1.e-6)
        aid,err=instance.get_index_of_species("H+")          
        instance.set_abundance(id,aid,1.e-4)

        instance.commit_particles()

        yr=365*24*3600.
        err=instance.evolve_model(10000.*yr)
        self.assertEqual(err,0)
        time,err=instance.get_time()
        self.assertEqual(err,0)
        self.assertEqual(time,10000.*yr)

        first,last,err=instance.get_firstlast_abundance()
        for i in range(first,last+1):
          x,err=instance.get_abundance(id,i)
          self.assertEqual(err,0)
          name,err=instance.get_name_of_species(i)
          print i,name,x



    #~ def test5(self):
        #~ print "Test 1: simple evolve test"
#~ 
        #~ instance = TDCInterface(**default_options)
        #~ self.assertEqual(0, instance.initialize_code())
        #~ self.assertEqual(0, instance.commit_parameters())
#~ 
        #~ dens=1.e5
        #~ t=500.
        #~ ion=1.e-11
        #~ id,err=instance.new_particle(dens,t,ion)
#~ 
        #~ self.assertEqual(err,0)
      #~ 
        #~ self.assertEqual(instance.commit_particles(),0)
#~ 
        #~ dens_,t_,ion_,err=instance.get_state(id)
#~ 
        #~ self.assertEqual(err,0)
#~ 
        #~ self.assertEqual(dens_,dens)
        #~ self.assertEqual(t_,t)
        #~ self.assertEqual(ion_,ion)
#~ 
        #~ err=instance.evolve_model( 94.112609921226252)
#~ 
        #~ x,err=instance.get_abundance(id,1)
        #~ self.assertEqual(0,err)
        #~ self.assertAlmostRelativeEqual(x,0.38895405860649546,7)
        #~ x,err=instance.get_abundance(id,2)
        #~ self.assertEqual(0,err)
        #~ self.assertAlmostRelativeEqual(x,.30526706600870018,7)
        #~ x,err=instance.get_abundance(id,3)
        #~ self.assertEqual(0,err)
        #~ self.assertAlmostRelativeEqual(x,1.47184343246716795E-004,7)
        #~ x,err=instance.get_abundance(id,4)
        #~ self.assertEqual(0,err)
        #~ self.assertAlmostRelativeEqual(x,7.08149953046615757E-008,7)
        #~ x,err=instance.get_abundance(id,5)
        #~ self.assertEqual(0,err)
        #~ self.assertAlmostRelativeEqual(x, 6.52897069295877081E-005,7)
#~ 
        #~ self.assertEqual(0, err)
#~ 
        #~ self.assertEqual(0, instance.cleanup_code())
#~ 
        #~ instance.stop()
#~ 
    #~ def test6(self):
        #~ print "Test 1: simple evolve test 2"
#~ 
        #~ instance = TDCInterface(**default_options)
        #~ self.assertEqual(0, instance.initialize_code())
        #~ self.assertEqual(0, instance.commit_parameters())
#~ 
        #~ dens=1.e5
        #~ t=500.
        #~ ion=1.e-11
        #~ id,err=instance.new_particle(dens,t,ion)
#~ 
        #~ self.assertEqual(err,0)
      #~ 
        #~ self.assertEqual(instance.commit_particles(),0)
#~ 
        #~ dens_,t_,ion_,err=instance.get_state(id)
#~ 
        #~ self.assertEqual(err,0)
#~ 
        #~ self.assertEqual(dens_,dens)
        #~ self.assertEqual(t_,t)
        #~ self.assertEqual(ion_,ion)
#~ 
        #~ err=instance.evolve_model( 50.)
        #~ err=instance.evolve_model( 94.112609921226252)
#~ 
        #~ x,err=instance.get_abundance(id,1)
        #~ self.assertEqual(0,err)
        #~ self.assertAlmostRelativeEqual(x,0.38895405860649546,3)
        #~ x,err=instance.get_abundance(id,2)
        #~ self.assertEqual(0,err)
        #~ self.assertAlmostRelativeEqual(x,.30526706600870018,3)
        #~ x,err=instance.get_abundance(id,3)
        #~ self.assertEqual(0,err)
        #~ self.assertAlmostRelativeEqual(x,1.47184343246716795E-004,3)
        #~ x,err=instance.get_abundance(id,4)
        #~ self.assertEqual(0,err)
        #~ self.assertAlmostRelativeEqual(x,7.08149953046615757E-008,3)
        #~ x,err=instance.get_abundance(id,5)
        #~ self.assertEqual(0,err)
        #~ self.assertAlmostRelativeEqual(x, 6.52897069295877081E-005,3)
#~ 
        #~ self.assertEqual(0, err)
#~ 
        #~ self.assertEqual(0, instance.cleanup_code())
#~ 
        #~ instance.stop()
#~ 
    #~ def test7(self):
        #~ print "Test 1: simple evolve test 3"
#~ 
        #~ instance = TDCInterface(**default_options)
        #~ self.assertEqual(0, instance.initialize_code())
        #~ self.assertEqual(0, instance.commit_parameters())
#~ 
        #~ dens=1.e5
        #~ t=500.
        #~ ion=1.e-11
        #~ id1,err=instance.new_particle(dens,t,ion)
        #~ id2,err=instance.new_particle(dens,t,ion)
#~ 
        #~ self.assertEqual(err,0)
      #~ 
        #~ self.assertEqual(instance.commit_particles(),0)
#~ 
        #~ err=instance.evolve_model( 94.112609921226252)
#~ 
        #~ x,err=instance.get_abundance(id1,1)
        #~ self.assertEqual(0,err)
        #~ self.assertAlmostRelativeEqual(x,0.38895405860649546,7)
        #~ x,err=instance.get_abundance(id1,2)
        #~ self.assertEqual(0,err)
        #~ self.assertAlmostRelativeEqual(x,.30526706600870018,7)
        #~ x,err=instance.get_abundance(id1,3)
        #~ self.assertEqual(0,err)
        #~ self.assertAlmostRelativeEqual(x,1.47184343246716795E-004,7)
        #~ x,err=instance.get_abundance(id1,4)
        #~ self.assertEqual(0,err)
        #~ self.assertAlmostRelativeEqual(x,7.08149953046615757E-008,7)
        #~ x,err=instance.get_abundance(id1,5)
        #~ self.assertEqual(0,err)
        #~ self.assertAlmostRelativeEqual(x, 6.52897069295877081E-005,7)
#~ 
        #~ x,err=instance.get_abundance(id2,1)
        #~ self.assertEqual(0,err)
        #~ self.assertAlmostRelativeEqual(x,0.38895405860649546,7)
        #~ x,err=instance.get_abundance(id2,2)
        #~ self.assertEqual(0,err)
        #~ self.assertAlmostRelativeEqual(x,.30526706600870018,7)
        #~ x,err=instance.get_abundance(id2,3)
        #~ self.assertEqual(0,err)
        #~ self.assertAlmostRelativeEqual(x,1.47184343246716795E-004,7)
        #~ x,err=instance.get_abundance(id2,4)
        #~ self.assertEqual(0,err)
        #~ self.assertAlmostRelativeEqual(x,7.08149953046615757E-008,7)
        #~ x,err=instance.get_abundance(id2,5)
        #~ self.assertEqual(0,err)
        #~ self.assertAlmostRelativeEqual(x, 6.52897069295877081E-005,7)
#~ 
        #~ self.assertEqual(0, err)
#~ 
        #~ self.assertEqual(0, instance.cleanup_code())
#~ 
        #~ instance.stop()
#~ 
    #~ def test8(self):
        #~ print "Test 1: simple evolve test 4"
#~ 
        #~ instance = TDCInterface(**default_options)
        #~ self.assertEqual(0, instance.initialize_code())
        #~ self.assertEqual(0, instance.commit_parameters())
#~ 
        #~ dens=1.e5
        #~ t=500.
        #~ ion=1.e-11
        #~ id1,err=instance.new_particle(dens,t,ion)
        #~ id2,err=instance.new_particle(2*dens,t,ion)
#~ 
        #~ self.assertEqual(err,0)
      #~ 
        #~ self.assertEqual(instance.commit_particles(),0)
        #~ 
        #~ err=instance.evolve_model( 94.112609921226252)
#~ 
        #~ x,err=instance.get_abundance(id1,3)
        #~ self.assertEqual(0,err)
        #~ self.assertAlmostRelativeEqual(x,1.47184343246716795E-004,7)
        #~ x,err=instance.get_abundance(id2,3)
        #~ self.assertEqual(0,err)
        #~ self.assertAlmostRelativeEqual(x,6.7890303322737115e-05,7)
#~ 
        #~ instance.stop()
#~ 
#~ 
#~ class TestTDC(TestWithMPI):
    #~ 
    #~ def is_fortan_version_up_to_date(self):
        #~ try:
            #~ from amuse import config
            #~ is_configured = hasattr(config, 'compilers')
            #~ if is_configured:
                #~ is_configured = hasattr(config.compilers, 'gfortran_version')
        #~ except ImportError:
            #~ is_configured = False
    #~ 
        #~ if is_configured:
            #~ if not config.compilers.gfortran_version:
                #~ return True
            #~ 
            #~ try:
                #~ parts = [int(x) for x in config.compilers.gfortran_version.split('.')]
            #~ except:
                #~ parts = []
                #~ 
            #~ if len(parts) < 2:
                #~ return True
                #~ 
            #~ return parts[0] >= 4 and parts[1] > 1
        #~ else:
            #~ return True
            #~ 
    #~ def check_fortran_version(self):
        #~ if not self.is_fortan_version_up_to_date():
            #~ self.skip('cannot compile, fortran module names cannot be resolved correctly in this gfortran version')
            #~ 
    #~ def setUp(self):
        #~ super(TestWithMPI, self).setUp()
        #~ self.check_fortran_version()
    #~ 
    #~ def test0(self):
        #~ print "test1: basic startup and flow"
        #~ instance=TDC()
        #~ self.assertEquals(instance.get_name_of_current_state(), 'UNINITIALIZED')
        #~ instance.initialize_code()
        #~ self.assertEquals(instance.get_name_of_current_state(), 'INITIALIZED')
        #~ instance.commit_parameters()
        #~ self.assertEquals(instance.get_name_of_current_state(), 'EDIT')
        #~ instance.commit_particles()
        #~ self.assertEquals(instance.get_name_of_current_state(), 'RUN')
#~ 
        #~ instance.cleanup_code()
        #~ instance.stop()
#~ 
    #~ def test1(self):
        #~ print "test1: adding particles"
#~ 
        #~ instance=TDC()
#~ 
        #~ parts=self.makeparts(5)
                #~ 
        #~ self.assertEqual(len(instance.particles),0)
        #~ instance.particles.add_particles(parts)
        #~ self.assertEqual(len(instance.particles),len(parts))
#~ 
        #~ self.assertEquals(instance.get_name_of_current_state(), 'EDIT')
#~ 
        #~ part2=instance.particles.copy()
#~ 
        #~ self.assertAlmostRelativeEquals(parts.number_density,part2.number_density,12)
        #~ self.assertAlmostRelativeEquals(parts.temperature,part2.temperature,12)
        #~ self.assertAlmostRelativeEquals(parts.ionrate,part2.ionrate,12)
#~ 
        #~ self.assertAlmostRelativeEquals(part2.abundances, [
#~ 0.29910782e0,
#~ 0.35040982e0 ,
#~ 1.0037100e-08 ,
#~ 7.8517447e-14 ,
 #~ 1.6757061e-10 ,
 #~ 9.1403564e-15 ,
   #~ 1.90e-7 ,
   #~ 0.00019687628e0, 
  #~ 2.3368063e-14 ,
  #~ 2.3042991e-27 ,
  #~ 1.5078215e-06 ,
   #~ 1.1743749e-10 ,
  #~ 1.4642501e-15 ,
    #~ 4.6442955e-09, 
   #~ 1.2046541e-13 ,
   #~ 1.1560869e-07 ,
  #~ 2.0393864e-14 ,
   #~ 2.2347066e-10 ,
   #~ 9.9354275e-09 ,
   #~ 2.8124190e-09 ,
  #~ 3.2763012e-21 ,
   #~ 5.5126847e-15 ,
   #~ 6.2356408e-17 ,
  #~ 1.2307295e-15 ,
  #~ 1.0298243e-15 ,
   #~ 8.6555075e-13 ,
   #~ 9.9986328e-05 ,
   #~ 2.1199001e-15 ,
   #~ 9.8358943e-10 ,
   #~ 0.10000000e0 ,
  #~ 2.9180687e-11 ,
  #~ 4.6013858e-17 ,
   #~ 5.5236247e-08 ,
   #~ 1.7576375e-07 ,
  #~ 2.6748393e-20 ],12)
#~ 
        #~ instance.cleanup_code()
        #~ instance.stop()
#~ 
    #~ def test2(self):
        #~ print "test2: adding particles w abund."
#~ 
        #~ instance=TDC()
#~ 
        #~ parts=self.makeparts(5)
        #~ 
        #~ parts.abundances=numpy.zeros((5,35))        
#~ 
        #~ for i in range(5):
          #~ parts[i].abundances=(numpy.array(range(35))+1)/36.
 #~ 
        #~ instance.particles.add_particles(parts)
#~ 
        #~ channel=parts.new_channel_to(instance.particles)
        #~ channel.copy() 
#~ 
#~ 
        #~ part2=instance.particles.copy()
#~ 
        #~ self.assertAlmostRelativeEquals(parts.number_density,part2.number_density,12)
        #~ self.assertAlmostRelativeEquals(parts.temperature,part2.temperature,12)
        #~ self.assertAlmostRelativeEquals(parts.ionrate,part2.ionrate,12)
#~ 
        #~ for i in range(5):  
          #~ self.assertAlmostRelativeEquals(part2[i].abundances,(numpy.array(range(35))+1)/36. ,12)
#~ 
        #~ instance.cleanup_code()
        #~ instance.stop()
#~ 
    #~ def test3(self):
        #~ print "test3: evolve test"
#~ 
        #~ instance=TDC(**default_options)
#~ 
        #~ parts=Particles(5)
        #~ parts.number_density=1.e5 | units.cm**-3
        #~ parts.temperature=500 | units.K
        #~ parts.ionrate=1.e-11 | units.s**-1
        #~ parts.abundances=numpy.zeros((5,35))        
#~ 
        #~ instance.particles.add_particles(parts)
 #~ 
        #~ my_yr=24*3600*365. *units.s
#~ 
        #~ instance.evolve_model( 94.112609921226252 | my_yr )
#~ 
        #~ for i in range(5):  
          #~ self.assertAlmostRelativeEqual(instance.particles[i].abundances[0:5],
               #~ [0.38895405860649546,.30526706600870018,1.47184343246716795E-004,
                #~ 7.08149953046615757E-008,6.52897069295877081E-005],7)
#~ 
        #~ instance.stop()
#~ 
#~ 
    #~ def makeparts(self,N):
        #~ parts=Particles(N)
        #~ numpy.random.seed(1234567)
        #~ parts.number_density=(numpy.random.random(N)*1.e5+1.e5)| units.cm**-3
        #~ parts.temperature=(numpy.random.random(N)*500+100)| units.K
        #~ parts.ionrate=(numpy.random.random(N)*1.e-11+1.e-17)| units.s**-1
        #~ return parts
#~ 
