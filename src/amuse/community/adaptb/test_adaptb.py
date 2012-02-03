## mpiexec nosetests -v test_adaptb.py
import math

from amuse.community import *
from amuse.test.amusetest import TestWithMPI

from functions import *

from .interface import AdaptbInterface
from .interface import Adaptb

default_options = dict() #redirection='none')

class AdaptbInterfaceTests(TestWithMPI):
    
    def test1(self):
        print "Test AdaptbInterface initialization"
        instance = AdaptbInterface(**default_options)
        self.assertEquals(0, instance.initialize_code())
        self.assertEquals(0, instance.commit_parameters())
        self.assertEquals(0, instance.cleanup_code())
        instance.stop()

    def test2(self):
        print "Test AdaptbInterface new_particle / get_state"
        instance = AdaptbInterface(**default_options)
        self.assertEquals(0, instance.initialize_code())
        self.assertEquals(0, instance.commit_parameters())
        
        id, error = instance.new_particle(mass = 11.0, radius = 2.0, x = 0.0, y = 0.0, z = 0.0, vx = 0.0, vy = 0.0, vz = 0.0)
        self.assertEquals(0, error)
        self.assertEquals(0, id)
        id, error = instance.new_particle(mass = 21.0, radius = 5.0, x = 10.0, y = 0.0, z = 0.0, vx = 10.0, vy = 0.0, vz = 0.0)
        self.assertEquals(0, error)
        self.assertEquals(1, id)
        self.assertEquals(0, instance.commit_particles())
        
        retrieved_state1 = instance.get_state(0)
        retrieved_state2 = instance.get_state(1)
        self.assertEquals(0,  retrieved_state1['__result'])
        self.assertEquals(0,  retrieved_state2['__result'])
        self.assertEquals(11.0,  retrieved_state1['mass'])
        self.assertEquals(21.0,  retrieved_state2['mass'])
        self.assertEquals( 0.0,  retrieved_state1['x'])
        self.assertEquals(10.0,  retrieved_state2['x'])
        
        self.assertEquals(0, instance.cleanup_code())
        instance.stop()
    
    def test4(self):
        print "Test AdaptbInterface particle property getters/setters"
        instance = AdaptbInterface(**default_options)
        self.assertEquals(0, instance.initialize_code())
        self.assertEquals(0, instance.commit_parameters())
        self.assertEquals([0, 0], instance.new_particle(0.01, 0.1,  1, 0, 0,  0, 1, 0).values())
        self.assertEquals([1, 0], instance.new_particle(0.02, 0.1, -1, 0, 0,  0,-1, 0).values())
        self.assertEquals(0, instance.commit_particles())
        
        # getters
        mass, result = instance.get_mass(0)
        self.assertAlmostEquals(0.01, mass)
        self.assertEquals(0,result)
        radius, result = instance.get_radius(1)
        self.assertAlmostEquals(0.1, radius)
        self.assertEquals(0,result)
        self.assertEquals(-3, instance.get_mass(2)['__result']) # Particle not found
        self.assertEquals([ 1, 0, 0,  0], instance.get_position(0).values())
        self.assertEquals([-1, 0, 0,  0], instance.get_position(1).values())
        self.assertEquals([ 0, 1, 0,  0], instance.get_velocity(0).values())
        self.assertEquals([ 0,-1, 0,  0], instance.get_velocity(1).values())
        
        # setters
        self.assertEquals(0, instance.set_state(0, 0.01, 0.1, 1,2,3, 4,5,6))
        self.assertEquals([0.01, 0.1, 1.0,2.0,3.0, 4.0,5.0,6.0, 0], instance.get_state(0).values())
        self.assertEquals(0, instance.set_mass(0, 0.02))
        self.assertEquals([0.02, 0.1, 1.0,2.0,3.0, 4.0,5.0,6.0, 0], instance.get_state(0).values())
        self.assertEquals(0, instance.set_radius(0, 0.2))
        self.assertEquals([0.02, 0.2, 1.0,2.0,3.0, 4.0,5.0,6.0, 0], instance.get_state(0).values())
        self.assertEquals(0, instance.set_position(0, 10,20,30))
        self.assertEquals([0.02, 0.2, 10.0,20.0,30.0, 4.0,5.0,6.0, 0], instance.get_state(0).values())
        self.assertEquals(0, instance.set_velocity(0, 40,50,60))
        self.assertEquals([0.02, 0.2, 10.0,20.0,30.0, 40.0,50.0,60.0, 0], instance.get_state(0).values())

        self.assertEquals(0, instance.cleanup_code())
        instance.stop()

    def test5(self):
        print "Test AdaptbInterface parameters"
        instance = AdaptbInterface(**default_options)
        self.assertEquals(0, instance.initialize_code())
        
	# bs tolerance
        self.assertEquals(["1e-6", 0], instance.get_bs_tolerance().values())
        self.assertEquals(0, instance.set_bs_tolerance("1e-8"))
        #self.assertEquals(["1e-8", 0], instance.get_bs_tolerance().values())
        self.assertEquals(1e-8, eval(instance.get_bs_tolerance()["bs_tolerance"]))
        
	# word length
        self.assertEquals([64, 0], instance.get_word_length().values())
        self.assertEquals(0, instance.set_word_length(80))
        self.assertEquals([80, 0], instance.get_word_length().values())

	# softening
        self.assertEquals(["0", 0], instance.get_eps2().values())
        self.assertEquals(0, instance.set_eps2("2e-1"))
        #print instance.get_eps2()["eps2"]
        #print instance.get_eps2()["eps2"].__class__
        #print eval("1+2")
        self.assertEquals(2.0e-1, eval(instance.get_eps2()["eps2"]))
        
	# print intervals
        #self.assertEquals(["1e-1", 0], instance.get_dt_print().values())
	self.assertEquals(1e-1, eval(instance.get_dt_print()["dt_print"]))
        self.assertEquals(0, instance.set_dt_print("1e-2"))
        #self.assertEquals(["1e-2", 0], instance.get_dt_print().values())
        self.assertEquals(1e-2, eval(instance.get_dt_print()["dt_print"]))

	# max cpu time
        self.assertEquals(["3600", 0], instance.get_max_cpu_time().values())
        self.assertEquals(0, instance.set_max_cpu_time("120"))
        self.assertEquals(["120", 0], instance.get_max_cpu_time().values())
        
	# output dir
        self.assertEquals(["./", 0], instance.get_adaptb_output_directory().values())
        self.assertEquals(0, instance.set_adaptb_output_directory("./out"))
        self.assertEquals(["./out/", 0], instance.get_adaptb_output_directory().values())
        self.assertEquals(0, instance.set_adaptb_output_directory(instance.output_directory))
        self.assertEquals([instance.output_directory+"/", 0], instance.get_adaptb_output_directory().values())
        
        self.assertEquals(0, instance.commit_parameters())
        self.assertEquals(0, instance.cleanup_code())
        instance.stop()
    
    def test6(self):
        print "Test AdaptbInterface evolve_model, equal-mass binary"
        instance = AdaptbInterface(**default_options)
        self.assertEquals(0, instance.initialize_code())

        self.assertEquals(0, instance.set_dt_print("1e-1"))

        self.assertEquals(0, instance.set_word_length(64))
        self.assertEquals(0, instance.set_bs_tolerance("1e-8"))

        self.assertEquals(0, instance.commit_parameters())
        
        self.assertEquals([0, 0], instance.new_particle(0.5, 0.0,  0.5, 0, 0,  0, 0.5, 0).values())
        self.assertEquals([1, 0], instance.new_particle(0.5, 0.0, -0.5, 0, 0,  0,-0.5, 0).values())
        self.assertEquals(0, instance.commit_particles())
        
        self.assertEquals(0, instance.evolve_model(math.pi)) # half an orbit
        for result, expected in zip(instance.get_position(0).values(), [-0.5, 0.0, 0.0, 0]):
            self.assertAlmostEquals(result, expected, 5)
        
        self.assertEquals(0, instance.evolve_model(2 * math.pi)) # full orbit
        for result, expected in zip(instance.get_position(0).values(), [0.5, 0.0, 0.0, 0]):
            self.assertAlmostEquals(result, expected, 5)
        
        self.assertEquals(0, instance.cleanup_code())
        instance.stop()

    def xtest7(self):
        print "Test AdaptbInterface evolve_model, pythagorean problem"
        instance = AdaptbInterface(**default_options)
        self.assertEquals(0, instance.initialize_code())

        self.assertEquals(0, instance.set_dt_print("1e-1"))

        self.assertEquals(0, instance.set_bs_tolerance("1e-6"))
        self.assertEquals(0, instance.set_word_length(64))

        self.assertEquals(0, instance.commit_parameters())

        self.assertEquals([0, 0], instance.new_particle("3", "0",  "1",  "3", "0", "0", "0", "0").values())
        self.assertEquals([1, 0], instance.new_particle("4", "0", "-2", "-1", "0", "0", "0", "0").values())
        self.assertEquals([2, 0], instance.new_particle("5", "0",  "1", "-1", "0", "0", "0", "0").values())
        self.assertEquals(0, instance.commit_particles())
        
        self.assertEquals(0, instance.evolve_model(100))
        
        self.assertEquals(0, instance.cleanup_code())
        instance.stop()

    def xtest8(self):
        print "Test AdaptbInterface evolve_model, pythagorean problem, show convergence"

	tolerance = ["1e-0", "1e-2", "1e-4", "1e-6", "1e-8", "1e-10", "1e-12", "1e-14", "1e-16", "1e-18", "1e-20", "1e-22", "1e-24"]
	word_length = [64,    64,     64   ,     64,     64,     64,      80,      80,      96,      96,      112,    112,     128]
        tolerance = tolerance[8:9]
        word_length = word_length[8:9]

	tcpu = []
	dE = []
	x0 = []
	vx0 = []

	i=0
	while i<len(tolerance):

          instance = AdaptbInterface(**default_options)
          self.assertEquals(0, instance.initialize_code())

          self.assertEquals(0, instance.set_dt_print("1e-1"))

          self.assertEquals(0, instance.set_bs_tolerance(tolerance[i]))
          self.assertEquals(0, instance.set_word_length(word_length[i]))

          self.assertEquals(0, instance.commit_parameters())

          self.assertEquals([0, 0], instance.new_particle("3", "0",  "1",  "3", "0", "0", "0", "0").values())
          self.assertEquals([1, 0], instance.new_particle("4", "0", "-2", "-1", "0", "0", "0", "0").values())
          self.assertEquals([2, 0], instance.new_particle("5", "0",  "1", "-1", "0", "0", "0", "0").values())
          self.assertEquals(0, instance.commit_particles())
        
          self.assertEquals(0, instance.evolve_model(100))
        
          self.assertEquals(0, instance.cleanup_code())
          instance.stop()

	  file_log = "file.log"
	  file_out = "file.out"

	  data_log = read_log(file_log)
	  data_out = read_out(file_out)

	  tcpu.append( float(data_log[9]) )
	  dE.append( float(data_log[10]) )
	  x0.append(data_out[0])
	  vx0.append(data_out[1])

	  print tolerance[i], word_length[i], tcpu[i], dE[i], x0[i], vx0[i]

	  i += 1

    def xtest9(self):
        print "Test AdaptbInterface evolve_model, pythagorean problem, plot orbits"

	tolerance = ["1e-4", "1e-8"] #, "1e-12", "1e-16"]
	word_length = [64  ,     64] #,      80,      96]

	x1 = []
	y1 = []
	x2 = []
	y2 = []
	x3 = []
	y3 = []

	i=0
	while i<len(tolerance):

          instance = AdaptbInterface(**default_options)
          self.assertEquals(0, instance.initialize_code())

          self.assertEquals(0, instance.set_dt_print("1e-1"))

          self.assertEquals(0, instance.set_bs_tolerance(tolerance[i]))
          self.assertEquals(0, instance.set_word_length(word_length[i]))

          self.assertEquals(0, instance.commit_parameters())

          self.assertEquals([0, 0], instance.new_particle("3", "0",  "1",  "3", "0", "0", "0", "0").values())
          self.assertEquals([1, 0], instance.new_particle("4", "0", "-2", "-1", "0", "0", "0", "0").values())
          self.assertEquals([2, 0], instance.new_particle("5", "0",  "1", "-1", "0", "0", "0", "0").values())
          self.assertEquals(0, instance.commit_particles())
        
          self.assertEquals(0, instance.evolve_model(100))
        
          self.assertEquals(0, instance.cleanup_code())
          instance.stop()

	  file_out = "file.out"
	  my_x1, my_y1, my_x2, my_y2, my_x3, my_y3 = read_xy(file_out)

	  x1.append(my_x1)
	  y1.append(my_y1)
	  x2.append(my_x2)
	  y2.append(my_y2)
	  x3.append(my_x3)
	  y3.append(my_y3)

	  i += 1

	print x1


