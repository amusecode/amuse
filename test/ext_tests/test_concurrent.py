from amuse.test.amusetest import TestWithMPI
from amuse.ext import concurrent

from amuse.rfi.core import *
from amuse.ic.plummer import new_plummer_sphere
from amuse.units import nbody_system

import numpy

class ConcurrentTestingInterface(PythonCodeInterface):
    
    def __init__(self, implementation_factory, **options):
        PythonCodeInterface.__init__(self, implementation_factory = implementation_factory, **options)
        
    @legacy_function
    def do_concurrent_run():
        function = LegacyFunctionSpecification()  
        function.addParameter('output', dtype='string', direction=function.OUT)
        function.result_type = 'int32'
        function.can_handle_array = True 
        function.id = 10
        return function
        
        
class Test1Implementation(object):
    
    def create_particles(self):
        self.particles = new_plummer_sphere(100)
        
    def do_concurrent_run(self, output):
        processes = concurrent.MPIConcurrentProcesses()
        processes.init()
        try:
            self.particles = None 
            processes.on_root(self.create_particles)
            shared_particles = processes.share(self.particles)
            shared_particles.distribute()
            assert len(shared_particles) == 100 , '{0} != 100'.format(len(shared_particles))
            output.value = str('success')
            result = 0
        except Exception as ex:
            output.value = str(ex)
            result = -1
        
        return self.reduce_result(processes, result)
        
    def reduce_result(self, processes, result):
        mpi_comm = processes.mpi_comm
        from mpi4py import MPI
        input = numpy.zeros(1,  dtype='int64')
        output = numpy.zeros(1,  dtype='int64')
        
        input[0] = result
            
        if mpi_comm.rank == 0:
            mpi_comm.Reduce(
                [input, MPI.INTEGER8], 
                [output, MPI.INTEGER8],
                op=MPI.SUM, 
                root=0
            )
        else:
            mpi_comm.Reduce(
                [input, MPI.INTEGER8], 
                [output, MPI.INTEGER8],
                op=MPI.SUM, 
                root=0
            )
            
        return output[0]
        
class Test3Implementation(Test1Implementation):
    
    def do_concurrent_run(self, output):
        processes = concurrent.MPIConcurrentProcesses()
        processes.init()
        test = TestMPIConcurrentProcesses('test1')
        try:
            self.particles = None 
            processes.on_root(self.create_particles)
            shared_particles = processes.share(self.particles)
            shared_particles.distribute()
            test.assertEquals(len(shared_particles), 100)
            test.assertAlmostRelativeEquals(shared_particles.mass.sum() , 1.0 | nbody_system.mass)
            output.value = str('success')
            result = 0
        except Exception as ex:
            print ex
            output.value = str(ex)
            result = -1
        
        return self.reduce_result(processes, result)
        
class Test4Implementation(Test1Implementation):
    
    def root_potential_energy(self):
        self.potential_energy = self.particles.potential_energy(G = nbody_system.G)
        
    def distributed_potential_energy_on_root(self):
        test = TestMPIConcurrentProcesses('test1')        
        self.distributed_potential_energy = self.shared_particles.potential_energy(G = nbody_system.G)
        test.assertAlmostRelativeEquals(self.distributed_potential_energy, self.potential_energy)
            
    def distributed_potential_energy_not_on_root(self):
        self.distributed_potential_energy = self.shared_particles.potential_energy(G = nbody_system.G)
        
        
    def do_concurrent_run(self, output):
        processes = concurrent.MPIConcurrentProcesses()
        processes.init()
        test = TestMPIConcurrentProcesses('test1')
        self.particles = None 
        self.potential_energy = 0
        try:
            processes.on_root(self.create_particles)
            processes.on_root(self.root_potential_energy)
            self.shared_particles = processes.share(self.particles)
            self.shared_particles.distribute()
            processes.call(
                self.distributed_potential_energy_on_root,
                self.distributed_potential_energy_not_on_root,
            )
            output.value = str('success')
            result = 0
        except Exception as ex:
            print ex
            output.value = str(ex)
            result = -1
        
        return self.reduce_result(processes, result)
        
class TestMPIConcurrentProcesses(TestWithMPI):

    def test1(self):
        x = ConcurrentTestingInterface(implementation_factory = Test1Implementation)
        output, error  = x.do_concurrent_run()
        self.assertEqual(error, 0, msg = output)

    def test2(self):
        x = ConcurrentTestingInterface(
            implementation_factory = Test1Implementation,
            number_of_workers = 4
        )
        output, error  = x.do_concurrent_run()
        self.assertEqual(error, 0, msg = output)
    
    def test3(self):
        x = ConcurrentTestingInterface(
            implementation_factory = Test3Implementation,
            number_of_workers = 4
        )
        output, error  = x.do_concurrent_run()
        self.assertEqual(error, 0, msg = output)
        
    def test4(self):
        x = ConcurrentTestingInterface(
            implementation_factory = Test4Implementation,
            number_of_workers = 4
        )
        output, error  = x.do_concurrent_run()
        self.assertEqual(error, 0, msg = output)
    
