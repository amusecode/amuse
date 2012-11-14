from amuse.test.amusetest import TestWithMPI

#cello
from amuse.support.codes import stopping_conditions
from amuse.support.interface import InCodeComponentImplementation

import subprocess
import os
import shlex

from amuse.units import nbody_system
from amuse.units import units
from amuse import datamodel
from amuse.support.exceptions import AmuseException
from amuse.rfi.tools import create_c
from amuse.rfi.tools import create_fortran
from amuse.rfi import channel
from amuse.rfi.core import *
from amuse.community import NO_UNIT
codestring = """
#include <stopcond.h>
#ifdef __cplusplus
extern "C" {
#endif
int initialize_code()
{
    
    // AMUSE STOPPING CONDITIONS SUPPORT
    supported_conditions = COLLISION_DETECTION_BITMAP | PAIR_DETECTION_BITMAP | TIMEOUT_DETECTION_BITMAP | OUT_OF_BOX_DETECTION_BITMAP;
    // -----------------------
    return 0;
}
#ifdef __cplusplus
}
#endif
"""

codestringF = """
      FUNCTION initialize_code()
      IMPLICIT NONE
      include "stopcond.inc"
      INTEGER :: initialize_code
      INTEGER :: set_support_for_condition
      INTEGER :: return
      initialize_code = 0
      return = set_support_for_condition(COLLISION_DETECTION)
      return = set_support_for_condition(PAIR_DETECTION)
      RETURN
      END FUNCTION
"""

codestringFModule = """
MODULE AmuseInterface
CONTAINS
      FUNCTION initialize_code()
          use StoppingConditions
          IMPLICIT NONE
          INTEGER :: initialize_code
          INTEGER :: return
          initialize_code = 0
          return = set_support_for_condition(COLLISION_DETECTION)
          return = set_support_for_condition(PAIR_DETECTION)
      END FUNCTION
      
      FUNCTION fire_condition(condition_to_set, particle_index1, particle_index2, rank)
          use StoppingConditions
          use mpi
          IMPLICIT NONE
          INTEGER :: fire_condition
          INTEGER :: my_rank
          INTEGER :: error, stopping_index
          INTEGER, intent(in) :: condition_to_set, particle_index1, particle_index2, rank
          fire_condition = 0
          call mpi_comm_rank(MPI_COMM_WORLD, my_rank, error)
          
          if (rank.GE.0 .AND. rank.NE.my_rank) then
            return
          end if
          stopping_index = next_index_for_stopping_condition()
          error = set_stopping_condition_info(stopping_index, condition_to_set)
          if(particle_index1 .GT.  0) then
            error = set_stopping_condition_particle_index(stopping_index, 0, particle_index1)
          end if
          if(particle_index2 .GT.  0) then
            error = set_stopping_condition_particle_index(stopping_index, 1, particle_index2)
          end if
      END FUNCTION
END MODULE
"""

class ForTestingInterface(CodeInterface, stopping_conditions.StoppingConditionInterface):
    def __init__(self, exefile, **options):
        CodeInterface.__init__(self, exefile, **options)

    include_headers = ['stopcond.h']
    @legacy_function
    def initialize_code():
        function = LegacyFunctionSpecification()  
        function.result_type = 'int32'
        function.can_handle_array = False
        return function     

    @legacy_function
    def reset_stopping_conditions():
        function = LegacyFunctionSpecification()  
        function.result_type = 'int32'
        function.can_handle_array = False
        return function     
        
    @legacy_function
    def next_index_for_stopping_condition():
        function = LegacyFunctionSpecification()
        function.result_type = 'int32'
        function.result_unit = NO_UNIT
        function.can_handle_array = False
        return function  
        
    @legacy_function
    def set_stopping_condition_info():
        function = LegacyFunctionSpecification()  
        function.addParameter('index', dtype='int32', direction=function.IN)
        function.addParameter('index_of_the_condition', dtype='int32', direction=function.IN)
        function.result_type = 'int32'
        function.can_handle_array = True
        return function  
    
    @legacy_function
    def set_stopping_condition_particle_index():
        function = LegacyFunctionSpecification()  
        function.addParameter('index', dtype='int32', direction=function.IN)
        function.addParameter('index_of_the_condition', dtype='int32', direction=function.IN)
        function.addParameter('index_of_the_particle', dtype='int32', direction=function.IN)
        function.result_type = 'int32'
        function.can_handle_array = True
        return function  


class ForTestingInterfaceFortranModule(ForTestingInterface):
    use_modules = ['StoppingConditions', 'AmuseInterface']
    
    @legacy_function
    def fire_condition():
        function = LegacyFunctionSpecification()  
        function.addParameter('condition_to_set', dtype='int32', direction=function.IN)
        function.addParameter('particle_index_1', dtype='int32', direction=function.IN)
        function.addParameter('particle_index_2', dtype='int32', direction=function.IN)
        function.addParameter('rank', dtype='int32', direction=function.IN, default = -1)
        function.result_type = 'int32'
        function.can_handle_array = True
        return function  
        
    
    @legacy_function
    def mpi_setup_stopping_conditions():
        function = LegacyFunctionSpecification()  
        function.result_type = 'int32'
        function.can_handle_array = False
        return function    
         
    @legacy_function
    def mpi_collect_stopping_conditions():
        function = LegacyFunctionSpecification()  
        function.result_type = 'int32'
        function.can_handle_array = False
        return function   
          
    @legacy_function
    def mpi_distribute_stopping_conditions():
        function = LegacyFunctionSpecification()  
        function.result_type = 'int32'
        function.can_handle_array = False
        return function   
          
class ForTesting(InCodeComponentImplementation):
    def __init__(self, exefile, **options):
        if 'community_interface' in options:
            interface = options['community_interface']
        else:
            interface = ForTestingInterface
        self.stopping_conditions = stopping_conditions.StoppingConditions(self)
        InCodeComponentImplementation.__init__(self, interface(exefile, **options), **options)
        self.my_particles = datamodel.Particles()
    
    def define_methods(self, object):
        self.stopping_conditions.define_methods(object)
    
    def new_particle(self, mass):
        particles = datamodel.Particles(len(mass))
        particles.mass = mass
        self.my_particles.add_particles(particles)
        return range(len(self.my_particles)-len(mass), len(self.my_particles))
    
    def get_mass(self, indices):
        return self.my_particles.mass[indices]
    
    def delete_particle(self, particle):
        self.my_particles.remove_particle(particle)
    
    def define_particle_sets(self, object):
        object.define_set('particles', 'index_of_the_particle')
        object.set_new('particles', 'new_particle')
        object.set_delete('particles', 'delete_particle')
        object.add_getter('particles', 'get_mass', names=("mass",))
        self.stopping_conditions.define_particle_set(object)

     
        
class TestInterface(TestWithMPI):
    
            
    def get_mpicc_name(self):
        try:
            from amuse import config
            is_configured = hasattr(config, 'mpi')
        except ImportError:
            is_configured = False
    
        if is_configured:
            return config.mpi.mpicc
        else:
            return os.environ['MPICC'] if 'MPICC' in os.environ else 'mpicc'
            
    def get_mpicxx_name(self):
        try:
            from amuse import config
            is_configured = hasattr(config, 'mpi')
        except ImportError:
            is_configured = False
    
        if is_configured:
            return config.mpi.mpicxx
        else:
            return os.environ['MPICXX'] if 'MPICXX' in os.environ else 'mpicxx'
    
    def cxx_compile(self, objectname, string):
  
        root, ext = os.path.splitext(objectname)
        sourcename = root + '.cc'
        if os.path.exists(objectname):
            os.remove(objectname)
        with open(sourcename, "w") as f:
            f.write(string)
        
        rootdir = self.get_amuse_root_dir()
        arguments = [self.get_mpicxx_name(), "-I",rootdir + "/lib/stopcond", "-c",  "-o", objectname, sourcename]
        process = subprocess.Popen(
            arguments,
            stdin = subprocess.PIPE,
            stdout = subprocess.PIPE,
            stderr = subprocess.PIPE
        )
        stdout, stderr = process.communicate()
        if not os.path.exists(objectname): # or process.poll() == 1:
            raise Exception("Could not compile {0}, error = {1} ({2})".format(objectname, stderr, ' '.join(arguments)))
    
    def c_build(self, exename, objectnames):
        rootdir = self.get_amuse_root_dir()
        
        arguments = [self.get_mpicxx_name()]
        arguments.extend(objectnames)
        arguments.append("-o")
        arguments.append(exename)
        arguments.extend(["-L"+rootdir+"/lib/stopcond","-lstopcond"])
        print ' '.join(arguments)
        process = subprocess.Popen(
            arguments,
            stdin = subprocess.PIPE,
            stdout = subprocess.PIPE,
            stderr = subprocess.PIPE
        )
        stdout, stderr = process.communicate()
        if process.returncode != 0:
            raise Exception("Could not build {0}, error = {1}".format(exename, stderr))
    
    def build_worker(self):
        path = os.path.abspath(self.get_path_to_results())
        codefile = os.path.join(path,"code.o")
        interfacefile = os.path.join(path,"interface.o")
        self.exefile = os.path.join(path,"c_worker")
        
        self.cxx_compile(codefile, codestring)
        
        uc = create_c.GenerateACHeaderStringFromASpecificationClass()
        uc.specification_class = ForTestingInterface
        uc.make_extern_c = True
        uc.ignore_functions_from_specification_classes = [stopping_conditions.StoppingConditionInterface]
        header =  uc.result

        uc = create_c.GenerateACSourcecodeStringFromASpecificationClass()

        uc.specification_class = ForTestingInterface
        code =  uc.result

        string = '\n\n'.join([header, code])
        
        #print string
        
        self.cxx_compile(interfacefile, string)
        self.c_build(self.exefile, [interfacefile, codefile] )
    
    def setUp(self):
        super(TestInterface, self).setUp()
        print "building"
        self.check_can_compile_modules()
        self.build_worker()
        
    def test1(self):
        print self.exefile
        instance = ForTestingInterface(self.exefile)
        instance.reset_stopping_conditions()
        next = instance.next_index_for_stopping_condition()
        next = instance.next_index_for_stopping_condition()
        instance.stop()
        self.assertEquals(next, 1)
        
    def test2(self):
        instance = ForTesting(self.exefile) #, debugger = "xterm")
        instance.initialize_code()
        self.assertTrue(instance.stopping_conditions.pair_detection.is_supported())
        self.assertTrue(instance.stopping_conditions.collision_detection.is_supported())
        self.assertFalse(instance.stopping_conditions.escaper_detection.is_supported())
        instance.stop()
        
    def test3(self):
        instance = ForTesting(self.exefile) #, debugger = "xterm")
        instance.initialize_code()
        self.assertFalse(instance.stopping_conditions.pair_detection.is_enabled())
        instance.stopping_conditions.pair_detection.enable()
        self.assertTrue(instance.stopping_conditions.pair_detection.is_enabled())
        instance.stopping_conditions.pair_detection.disable()
        self.assertFalse(instance.stopping_conditions.pair_detection.is_enabled())
        instance.stop()
        
    def test4(self):
        instance = ForTesting(self.exefile)
        instance.reset_stopping_conditions()
        next = instance.next_index_for_stopping_condition()
        self.assertFalse(instance.stopping_conditions.pair_detection.is_set())
        print next,instance.stopping_conditions.pair_detection.type
        instance.set_stopping_condition_info(next,instance.stopping_conditions.pair_detection.type)
         
        self.assertTrue(instance.stopping_conditions.pair_detection.is_set())
        instance.stop()
        
    def test5(self):
        instance = ForTesting(self.exefile)
        instance.reset_stopping_conditions()
        next = instance.next_index_for_stopping_condition()
        self.assertFalse(instance.stopping_conditions.pair_detection.is_set())

        instance.set_stopping_condition_info(next,instance.stopping_conditions.pair_detection.type)
        instance.set_stopping_condition_particle_index(next, 0, 11)
        instance.set_stopping_condition_particle_index(next, 1, 12)
        self.assertTrue(instance.stopping_conditions.pair_detection.is_set())
        self.assertEquals(11, instance.get_stopping_condition_particle_index(next, 0))
        self.assertEquals(12, instance.get_stopping_condition_particle_index(next, 1))
        instance.stop()

    def test6(self):
        instance = ForTesting(self.exefile)
        instance.reset_stopping_conditions()
        next = instance.next_index_for_stopping_condition()
        instance.set_stopping_condition_info(next,instance.stopping_conditions.out_of_box_detection.type)
        self.assertTrue(instance.stopping_conditions.out_of_box_detection.is_set())
        instance.stop()
    
    def test7(self):
        instance = ForTesting(self.exefile)
        particles = datamodel.Particles(20)
        particles.mass = range(1, 21) | units.kg
        instance.particles.add_particles(particles)
        instance.reset_stopping_conditions()
        
        pairs = [(11, 12), (0, 4), (3, 18), (7, 2)]
        next = instance.next_index_for_stopping_condition()
        self.assertFalse(instance.stopping_conditions.pair_detection.is_set())
        instance.set_stopping_condition_info(next,instance.stopping_conditions.pair_detection.type)
        instance.set_stopping_condition_particle_index(next, 0, pairs[0][0])
        instance.set_stopping_condition_particle_index(next, 1, pairs[0][1])
        self.assertEquals(11, instance.get_stopping_condition_particle_index(next, 0))
        self.assertEquals(12, instance.get_stopping_condition_particle_index(next, 1))
        self.assertTrue(instance.stopping_conditions.pair_detection.is_set())
        self.assertEquals(len(instance.stopping_conditions.pair_detection.particles(0)), 1)
        self.assertEquals(len(instance.stopping_conditions.pair_detection.particles(1)), 1)
        
        for index1, index2 in pairs[1:]:
            next = instance.next_index_for_stopping_condition()
            instance.set_stopping_condition_info(next,instance.stopping_conditions.pair_detection.type)
            instance.set_stopping_condition_particle_index(next, 0, index1)
            instance.set_stopping_condition_particle_index(next, 1, index2)
            self.assertEquals(index1, instance.get_stopping_condition_particle_index(next, 0))
            self.assertEquals(index2, instance.get_stopping_condition_particle_index(next, 1))
        self.assertEquals(len(instance.stopping_conditions.pair_detection.particles(0)), 4)
        self.assertEquals(len(instance.stopping_conditions.pair_detection.particles(1)), 4)
        self.assertEquals(len(instance.stopping_conditions.pair_detection.particles(2)), 0)
        
        self.assertEquals(instance.stopping_conditions.pair_detection.particles(0).mass, 
            [first + 1 for first, second in pairs] | units.kg)
        self.assertEquals(instance.stopping_conditions.pair_detection.particles(1).mass, 
            [second + 1 for first, second in pairs] | units.kg)
        instance.stop()
    
    def test8(self):
        instance = ForTesting(self.exefile)
        instance.initialize_code()
        self.assertFalse(instance.stopping_conditions.escaper_detection.is_supported())
        self.assertRaises(AmuseException, instance.stopping_conditions.escaper_detection.enable, expected_message=
            "Can't enable stopping condition 'escaper_detection', since 'ForTesting' does not support this condition.")
        instance.stop()
    


class _AbstractTestInterfaceFortran(TestWithMPI):
    
    def get_mpif90_name(self):
        try:
            from amuse import config
            is_configured = hasattr(config, 'mpi')
        except ImportError:
            is_configured = False
    
        if is_configured:
            return config.mpi.mpif95
        else:
            return os.environ['MPIFC'] if 'MPIFC' in os.environ else 'mpif90'
            
    def get_mpif90_arguments(self):
        name = self.get_mpif90_name()
        return list(shlex.split(name))

    def f90_compile(self, objectname, string):
        root, ext = os.path.splitext(objectname)
        sourcename = root + '.f90'
        mpidir = self.get_mpidir()
        if os.path.exists(objectname):
            os.remove(objectname)
        with open(sourcename, "w") as f:
            f.write(string)
        
        rootdir = self.get_amuse_root_dir()
        arguments = self.get_mpif90_arguments()
        arguments.extend(["-I","{0}/lib/stopcond{1}".format(rootdir, mpidir), "-c",  "-o", objectname, sourcename])
        process = subprocess.Popen(
            arguments,
            stdin = subprocess.PIPE,
            stdout = subprocess.PIPE,
            stderr = subprocess.PIPE
        )
        stdout, stderr = process.communicate()
        if not os.path.exists(objectname): # or process.poll() == 1:
            raise Exception("Could not compile {0}, error = {1} ({2})".format(objectname, stderr, ' '.join(arguments)))
    
    def get_libname(self):
        return 'stopcond'
    
    def get_mpidir(self):
        return ''
        
    def get_codestring(self):
        return CodeStringF
        
    def get_interface_class(self):
        return ForTestingInterface
    
    def get_number_of_workers(self):
        return 1
        
    def f90_build(self, exename, objectnames):
        rootdir = self.get_amuse_root_dir()
        
        arguments = self.get_mpif90_arguments()
        arguments.extend(objectnames)
        arguments.append("-o")
        arguments.append(exename)
        arguments.extend(["-L{0}/lib/stopcond".format(rootdir),"-l"+self.get_libname()])
        print 'build command:'
        print ' '.join(arguments)
        process = subprocess.Popen(
            arguments,
            stdin = subprocess.PIPE,
            stdout = subprocess.PIPE,
            stderr = subprocess.PIPE
        )
        stdout, stderr = process.communicate()
        if process.returncode != 0:
            raise Exception("Could not build {0}, error = {1}".format(exename, stderr))
    
    def build_worker(self):
        path = os.path.abspath(self.get_path_to_results())
        codefile = os.path.join(path,"codef90.o")
        interfacefile = os.path.join(path,"interfacef90.o")
        self.exefile = os.path.join(path,"fortran_worker")
        
        self.f90_compile(codefile, self.get_codestring())
        
        uf = create_fortran.GenerateAFortranSourcecodeStringFromASpecificationClass()
        uf.specification_class = self.get_interface_class()
        string =  uf.result
        
        #print string
        
        self.f90_compile(interfacefile, string)
        self.f90_build(self.exefile, [interfacefile, codefile] )
    
    def setUp(self):
        super(_AbstractTestInterfaceFortran, self).setUp()
        print "building"
        self.check_can_compile_modules()
        self.build_worker()


class _TestInterfaceFortranSingleProcess(_AbstractTestInterfaceFortran):
    
    def get_number_of_workers(self):
        return 1
        
    def test1(self):
        instance = ForTestingInterface(self.exefile, number_of_workers = self.get_number_of_workers())
        instance.reset_stopping_conditions()
        next = instance.next_index_for_stopping_condition()
        next = instance.next_index_for_stopping_condition()
        instance.stop()
        self.assertEquals(next, 1)

    def test2(self):
        instance = ForTesting(self.exefile, number_of_workers = self.get_number_of_workers()) #, debugger = "xterm")
        instance.initialize_code()
        self.assertTrue(instance.stopping_conditions.pair_detection.is_supported())
        self.assertTrue(instance.stopping_conditions.collision_detection.is_supported())
        self.assertFalse(instance.stopping_conditions.escaper_detection.is_supported())
        instance.stop()

    def test3(self):
        instance = ForTesting(self.exefile, number_of_workers = self.get_number_of_workers()) #, debugger = "xterm")
        instance.initialize_code()
        self.assertFalse(instance.stopping_conditions.pair_detection.is_enabled())
        instance.stopping_conditions.pair_detection.enable()
        self.assertTrue(instance.stopping_conditions.pair_detection.is_enabled())
        instance.stopping_conditions.pair_detection.disable()
        self.assertFalse(instance.stopping_conditions.pair_detection.is_enabled())
        instance.stop()

                
    def test4(self):
        instance = ForTesting(self.exefile, number_of_workers = self.get_number_of_workers())
        instance.reset_stopping_conditions()

        next = instance.next_index_for_stopping_condition()
        self.assertFalse(instance.stopping_conditions.pair_detection.is_set())
    
        instance.set_stopping_condition_info(next,instance.stopping_conditions.pair_detection.type)
         
        self.assertTrue(instance.stopping_conditions.pair_detection.is_set())
        instance.stop()
    
    def test5(self):
        instance = ForTesting(self.exefile, number_of_workers = self.get_number_of_workers())
        instance.reset_stopping_conditions()
        next = instance.next_index_for_stopping_condition()
        self.assertFalse(instance.stopping_conditions.pair_detection.is_set())
        
        instance.set_stopping_condition_info(next,instance.stopping_conditions.pair_detection.type)
        instance.set_stopping_condition_particle_index(next, 0, 11)
        instance.set_stopping_condition_particle_index(next, 1, 12)
        self.assertTrue(instance.stopping_conditions.pair_detection.is_set())
        self.assertEquals(11, instance.get_stopping_condition_particle_index(next, 0))
        self.assertEquals(12, instance.get_stopping_condition_particle_index(next, 1))
        instance.stop()
    
    def test6(self):
        instance = ForTesting(self.exefile, number_of_workers = self.get_number_of_workers())
        particles = datamodel.Particles(20)
        particles.mass = range(1, 21) | units.kg
        instance.particles.add_particles(particles)
        instance.reset_stopping_conditions()
        
        pairs = [(11, 12), (0, 4), (3, 18), (7, 2)]
        next = instance.next_index_for_stopping_condition()
        self.assertFalse(instance.stopping_conditions.pair_detection.is_set())
        instance.set_stopping_condition_info(next,instance.stopping_conditions.pair_detection.type)
        instance.set_stopping_condition_particle_index(next, 0, pairs[0][0])
        instance.set_stopping_condition_particle_index(next, 1, pairs[0][1])
        self.assertEquals(11, instance.get_stopping_condition_particle_index(next, 0))
        self.assertEquals(12, instance.get_stopping_condition_particle_index(next, 1))
        self.assertTrue(instance.stopping_conditions.pair_detection.is_set())
        self.assertEquals(len(instance.stopping_conditions.pair_detection.particles(0)), 1)
        self.assertEquals(len(instance.stopping_conditions.pair_detection.particles(1)), 1)
        
        for index1, index2 in pairs[1:]:
            next = instance.next_index_for_stopping_condition()
            instance.set_stopping_condition_info(next,instance.stopping_conditions.pair_detection.type)
            instance.set_stopping_condition_particle_index(next, 0, index1)
            instance.set_stopping_condition_particle_index(next, 1, index2)
            self.assertEquals(index1, instance.get_stopping_condition_particle_index(next, 0))
            self.assertEquals(index2, instance.get_stopping_condition_particle_index(next, 1))
        self.assertEquals(len(instance.stopping_conditions.pair_detection.particles(0)), 4)
        self.assertEquals(len(instance.stopping_conditions.pair_detection.particles(1)), 4)
        self.assertEquals(len(instance.stopping_conditions.pair_detection.particles(2)), 0)
        
        self.assertEquals(instance.stopping_conditions.pair_detection.particles(0).mass, 
            [first + 1 for first, second in pairs] | units.kg)
        self.assertEquals(instance.stopping_conditions.pair_detection.particles(1).mass, 
            [second + 1 for first, second in pairs] | units.kg)
        instance.stop()
    
    def test8(self):
        instance = ForTesting(self.exefile, number_of_workers = self.get_number_of_workers())
        instance.initialize_code()
        self.assertFalse(instance.stopping_conditions.escaper_detection.is_supported())
        self.assertRaises(AmuseException, instance.stopping_conditions.escaper_detection.enable, expected_message=
            "Can't enable stopping condition 'escaper_detection', since 'ForTesting' does not support this condition.")
        instance.stop()
    

class TestInterfaceFortran(_TestInterfaceFortranSingleProcess):
    
    def get_libname(self):
        return 'stopcond'
        
    def get_codestring(self):
        return codestringF
        
    def get_interface_class(self):
        return ForTestingInterface
    
class TestInterfaceFortranModule(_TestInterfaceFortranSingleProcess):
    
    def get_libname(self):
        return 'stopcondf'
        
    def get_codestring(self):
        return codestringFModule
        
    def get_interface_class(self):
        return ForTestingInterfaceFortranModule
    
class TestInterfaceFortranModuleMultiprocess(_AbstractTestInterfaceFortran):
    
    def get_libname(self):
        return 'stopcondfmpi'
        
    def get_codestring(self):
        return codestringFModule
        
    def get_interface_class(self):
        return ForTestingInterfaceFortranModule
    
    def get_number_of_workers(self):
        return 3
    
    def get_mpidir(self):
        return '/mpi'
        
    def test1(self):
        instance = ForTesting(
            self.exefile, 
            community_interface = ForTestingInterfaceFortranModule,
            number_of_workers = self.get_number_of_workers()
        )
        instance.initialize_code()
        instance.reset_stopping_conditions()
        instance.mpi_setup_stopping_conditions()
        
        pair_detection = instance.stopping_conditions.pair_detection
        
        particles = datamodel.Particles(20)
        particles.mass = range(1, 21) | units.kg
        instance.particles.add_particles(particles)
        
        instance.stopping_conditions.pair_detection.enable()
        
        instance.mpi_collect_stopping_conditions()
        
        print pair_detection.type
        instance.fire_condition(
            pair_detection.type,
            1, 2, -1
        )
        instance.mpi_distribute_stopping_conditions()
        self.assertTrue(pair_detection.is_set())
        self.assertEquals(len(pair_detection.particles(0)),self.get_number_of_workers()) 
        self.assertEquals(len(pair_detection.particles(1)),self.get_number_of_workers()) 
        self.assertEquals(pair_detection.particles(0).key,particles[1].key)
        self.assertEquals(pair_detection.particles(1).key,particles[2].key)
        self.assertEquals(pair_detection.particles(0).mass,[2,2,2] | units.kg) 
        self.assertEquals(pair_detection.particles(1).mass,[3,3,3] | units.kg) 
        instance.stop()
        
    
        
    def test2(self):
        instance = ForTesting(
            self.exefile, 
            community_interface = ForTestingInterfaceFortranModule,
            number_of_workers = self.get_number_of_workers()
        )
        instance.initialize_code()
        instance.reset_stopping_conditions()
        instance.mpi_setup_stopping_conditions()
        
        pair_detection = instance.stopping_conditions.pair_detection
        
        particles = datamodel.Particles(20)
        particles.mass = range(1, 21) | units.kg
        instance.particles.add_particles(particles)
        
        instance.stopping_conditions.pair_detection.enable()
        
        instance.mpi_collect_stopping_conditions()
        for rank in range(self.get_number_of_workers()):
            print pair_detection.type
            instance.fire_condition(
                pair_detection.type,
                1, 2, rank
            )
            instance.mpi_distribute_stopping_conditions()
            self.assertTrue(pair_detection.is_set())
            self.assertEquals(len(pair_detection.particles(0)),1) 
            self.assertEquals(len(pair_detection.particles(1)),1) 
            self.assertEquals(pair_detection.particles(0).key,particles[1].key)
            self.assertEquals(pair_detection.particles(1).key,particles[2].key)
            self.assertEquals(pair_detection.particles(0).mass,[2] | units.kg) 
            self.assertEquals(pair_detection.particles(1).mass,[3] | units.kg) 
            instance.reset_stopping_conditions()
            instance.stopping_conditions.pair_detection.enable()
            
        instance.stop()

    def test3(self):
        instance = ForTesting(
            self.exefile, 
            community_interface = ForTestingInterfaceFortranModule,
            number_of_workers = self.get_number_of_workers()
        )
        instance.initialize_code()
        instance.reset_stopping_conditions()
        instance.mpi_setup_stopping_conditions()
        
        pair_detection = instance.stopping_conditions.pair_detection
        
        particles = datamodel.Particles(20)
        particles.mass = range(1, 21) | units.kg
        instance.particles.add_particles(particles)
        
        instance.stopping_conditions.pair_detection.enable()
        
        instance.mpi_collect_stopping_conditions()
        
        instance.fire_condition(
            pair_detection.type,
            1, 2, 0
        )
        instance.fire_condition(
            pair_detection.type,
            3, 4, 1
        )
        instance.fire_condition(
            pair_detection.type,
            5, 6, 2
        )
        instance.mpi_distribute_stopping_conditions()
        self.assertTrue(pair_detection.is_set())
        self.assertEquals(len(pair_detection.particles(0)),3) 
        self.assertEquals(len(pair_detection.particles(1)),3) 
        self.assertEquals(pair_detection.particles(0).key[0],particles[1].key)
        self.assertEquals(pair_detection.particles(1).key[0],particles[2].key)
        self.assertEquals(pair_detection.particles(0).key[1],particles[3].key)
        self.assertEquals(pair_detection.particles(1).key[1],particles[4].key)
        self.assertEquals(pair_detection.particles(0).key[2],particles[5].key)
        self.assertEquals(pair_detection.particles(1).key[2],particles[6].key)
        instance.reset_stopping_conditions()
        instance.stopping_conditions.pair_detection.enable()
            
        instance.stop()
    
    

    def test4(self):
        instance = ForTesting(
            self.exefile, 
            community_interface = ForTestingInterfaceFortranModule,
            number_of_workers = self.get_number_of_workers()
        )
        instance.initialize_code()
        instance.reset_stopping_conditions()
        instance.mpi_setup_stopping_conditions()
        
        pair_detection = instance.stopping_conditions.pair_detection
        
        particles = datamodel.Particles(20)
        particles.mass = range(1, 21) | units.kg
        instance.particles.add_particles(particles)
        
        instance.stopping_conditions.pair_detection.enable()
        
        instance.mpi_collect_stopping_conditions()
        
        instance.fire_condition(
            pair_detection.type,
            -1, -1, -1
        )
        instance.mpi_distribute_stopping_conditions()
        self.assertTrue(pair_detection.is_set())
        self.assertEquals(len(pair_detection.particles(0)),0) 
            
        instance.stop()
