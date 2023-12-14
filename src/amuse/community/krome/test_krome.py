import os.path
import numpy
from amuse.test.amusetest import TestWithMPI

from amuse.community.krome.interface import KromeInterface,Krome,solar_abundances
from amuse.units import units
from amuse.datamodel import Particles

from amuse.io import read_set_from_file

#default_options={}
default_options = dict(redirection="none")
#default_options=dict(debugger="gdb")

class TestKromeInterface(TestWithMPI):

    def test1(self):
        print("Test 1: initialization")

        instance = self.new_instance_of_an_optional_code(KromeInterface, **default_options)
        self.assertEqual(0, instance.initialize_code())
        self.assertEqual(0, instance.commit_parameters())
        self.assertEqual(0, instance.cleanup_code())
        instance.stop()

    def test2(self):
        print("Test 1: add particle, get state")

        instance = self.new_instance_of_an_optional_code(KromeInterface, **default_options)
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
        print("Test 1: add 2 particles, get state")

        instance = self.new_instance_of_an_optional_code(KromeInterface, **default_options)
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
        print("Test 1: add 100 particles, get state")

        instance = self.new_instance_of_an_optional_code(KromeInterface, **default_options)
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
        print("Test 5: can we get species?")

        instance = self.new_instance_of_an_optional_code(KromeInterface, **default_options)
        
        first,last,err=instance.get_firstlast_abundance()
        self.assertEqual(err,0)
        self.assertTrue(last-first > 0)
        
        for i in range(first,last+1):
          name,err=instance.get_name_of_species(i)
          print(name)
          self.assertEqual(err,0)
          index,err=instance.get_index_of_species(name)
          self.assertEqual(i,index)

        instance.stop()

    def test6(self):
        print("Test 6: add 100 particles, remove particles")

        instance = self.new_instance_of_an_optional_code(KromeInterface, **default_options)
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
        print("Test 1: add particle, set abundances")

        instance = self.new_instance_of_an_optional_code(KromeInterface, **default_options)
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
          self.assertTrue((x>=0.) & (x<=1.))
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
        print("evolve test")

        instance = self.new_instance_of_an_optional_code(KromeInterface, **default_options)
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
          print(i,name,x)
        
    def test9(self):
        print("evolve test 2")

        instance = self.new_instance_of_an_optional_code(KromeInterface, **default_options)
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
          print(i,name,x)

    def test10(self):
        print("check initialization of abundances")

        instance = self.new_instance_of_an_optional_code(KromeInterface, **default_options)
        self.assertEqual(0, instance.initialize_code())
        self.assertEqual(0, instance.commit_parameters())

        dens=1.e2
        t=100.
        ion=2.e-17
        id,err=instance.new_particle(dens,t,ion)
        
        abundances={"E":0.000369180975425,
                     "H+":0.0001,"HE":0.0775, 
                     "C+":0.000269180975425, "SI": 3.2362683404e-05,
                     "O": 0.000489828841345}
         
        first,last,err=instance.get_firstlast_abundance()
        for i in range(first,last+1):
          x,err=instance.get_abundance(id,i)
          self.assertEqual(err,0)
          name,err=instance.get_name_of_species(i)
          if name in abundances:
            self.assertAlmostEqual(x,abundances[name],12)


    def test11(self):
        print("evolve test, comparison")

        instance = self.new_instance_of_an_optional_code(KromeInterface, **default_options)
        self.assertEqual(0, instance.initialize_code())
        self.assertEqual(0, instance.commit_parameters())

        dens=1.e2
        t=100.
        ion=2.e-17
        id,err=instance.new_particle(dens,t,ion)
        instance.evolve_model(1.e10)
        
        result1={}
        first,last,err=instance.get_firstlast_abundance()
        for i in range(first,last+1):
          x,err=instance.get_abundance(id,i)
          self.assertEqual(err,0)
          name,err=instance.get_name_of_species(i)
          result1[name]=x

        instance = self.new_instance_of_an_optional_code(KromeInterface, **default_options)
        self.assertEqual(0, instance.initialize_code())
        self.assertEqual(0, instance.commit_parameters())

        dens=1.e2
        t=100.
        ion=2.e-17
        id,err=instance.new_particle(dens,t,ion)
        instance.evolve_model(1.e9)
        instance.evolve_model(5.e9)
        instance.evolve_model(1.e10)
        
        result2={}
        first,last,err=instance.get_firstlast_abundance()
        for i in range(first,last+1):
          x,err=instance.get_abundance(id,i)
          self.assertEqual(err,0)
          name,err=instance.get_name_of_species(i)
          result2[name]=x

        for x in result1:
          self.assertAlmostEqual(result1[x],result2[x])


class TestKrome(TestWithMPI):
    def makeparts(self,N):
        parts=Particles(N)
        numpy.random.seed(1234567)
        parts.number_density=(numpy.random.random(N)*1.e5+1.e5)| units.cm**-3
        parts.temperature=(numpy.random.random(N)*500+100)| units.K
        parts.ionrate=(numpy.random.random(N)*1.e-11+1.e-17)| units.s**-1
        return parts

    def test0(self):
        print("test1: basic startup and flow")
        instance=self.new_instance_of_an_optional_code(Krome)
        self.assertEqual(instance.get_name_of_current_state(), 'UNINITIALIZED')
        instance.initialize_code()
        self.assertEqual(instance.get_name_of_current_state(), 'INITIALIZED')
        instance.commit_parameters()
        self.assertEqual(instance.get_name_of_current_state(), 'EDIT')
        instance.commit_particles()
        self.assertEqual(instance.get_name_of_current_state(), 'RUN')

        instance.cleanup_code()
        instance.stop()

    def test1(self):
        print("test1: adding particles")

        instance=self.new_instance_of_an_optional_code(Krome)

        parts=self.makeparts(5)
                
        self.assertEqual(len(instance.particles),0)
        instance.particles.add_particles(parts)
        self.assertEqual(len(instance.particles),len(parts))

        self.assertEqual(instance.get_name_of_current_state(), 'EDIT')

        part2=instance.particles.copy()

        self.assertAlmostRelativeEquals(parts.number_density,part2.number_density,12)
        self.assertAlmostRelativeEquals(parts.temperature,part2.temperature,12)
        self.assertAlmostRelativeEquals(parts.ionrate,part2.ionrate,12)

        for p in part2:
          i=instance.species["E"]
          self.assertAlmostEqual(p.abundances[i],0.000369180975425)
          i=instance.species["H+"]
          self.assertAlmostEqual(p.abundances[i],0.0001)
          i=instance.species["HE"]
          self.assertAlmostEqual(p.abundances[i],0.0775)
          i=instance.species["C+"]
          self.assertAlmostEqual(p.abundances[i],0.000269180975425)
          i=instance.species["SI"]
          self.assertAlmostEqual(p.abundances[i],3.2362683404e-05)
          i=instance.species["O"]
          self.assertAlmostEqual(p.abundances[i],0.000489828841345)
          
        instance.cleanup_code()
        instance.stop()


    def test2(self):
        print("test2: adding particles w abund.")

        instance=self.new_instance_of_an_optional_code(Krome)

        parts=self.makeparts(5)
        
        N=len(instance.species)

        parts.abundances=numpy.zeros((5,N))        

        for i in range(5):
          parts[i].abundances=(numpy.array(range(N))+1)/(N+1.)
 
        instance.particles.add_particles(parts)

        channel=parts.new_channel_to(instance.particles)
        channel.copy() 

        part2=instance.particles.copy()

        self.assertAlmostRelativeEquals(parts.number_density,part2.number_density,12)
        self.assertAlmostRelativeEquals(parts.temperature,part2.temperature,12)
        self.assertAlmostRelativeEquals(parts.ionrate,part2.ionrate,12)

        for i in range(5):  
          self.assertAlmostRelativeEquals(part2[i].abundances,parts[i].abundances ,12)

        instance.cleanup_code()
        instance.stop()

    def test3(self):
        print("test3: evolve test")

        instance=self.new_instance_of_an_optional_code(Krome,**default_options)

        parts=Particles(1)
        parts.number_density=1.e5 | units.cm**-3
        parts.temperature=50 | units.K
        parts.ionrate=2.e-17 | units.s**-1

        Ns=len(instance.species)

        parts.abundances=numpy.zeros((1,Ns))        

        instance.particles.add_particles(parts)
  
        instance.evolve_model( 1. | units.Myr )
 
        print(instance.particles.abundances)

        f=2*instance.particles[0].abundances[instance.species["H2"]]
        self.assertTrue(f> 0.95) # not much of a test..
        #~ for x,i in instance.species.items():
          #~ print x, instance.particles[0].abundances[i]
        
        instance.cleanup_code()
        instance.stop()

    def test4(self):
        print("test4: evolve test (10 part)")

        instance=self.new_instance_of_an_optional_code(Krome,**default_options)

        parts=Particles(10)
        parts.number_density=1.e5 | units.cm**-3
        parts.temperature=50 | units.K
        parts.ionrate=2.e-17 | units.s**-1

        Ns=len(instance.species)
        
        parts.abundances=numpy.zeros((10,Ns))        

        instance.particles.add_particles(parts)
  
        instance.evolve_model( 1. | units.Myr )

        f=2*instance.particles[0].abundances[instance.species["H2"]]
        self.assertTrue(f> 0.95) # not much of a test..
        #~ for x,i in instance.species.items():
          #~ print x, instance.particles[0].abundances[i]
        
        instance.cleanup_code()
        instance.stop()
