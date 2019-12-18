import sys
import traceback
import numpy
import threading

from amuse.units import units
from amuse.datamodel import ParticlesSuperset
from amuse.support.exceptions import AmuseException
from amuse.community.interface.se import merge_colliding_in_stellar_evolution_code
from amuse.support.options import option, OptionalAttributes




class ParallelStellarEvolution(OptionalAttributes):

    def __init__(
            self, 
            stellar_evolution_class, 
            number_of_workers=1, 
            individual_options=None, 
            **options
        ):
                
        OptionalAttributes.__init__(self, **options)
        
        self.code_factory = stellar_evolution_class
        self.number_of_workers = number_of_workers
        self.model_time = 0.0 | units.Myr
        
        if individual_options is None:
            options_list = [options] * number_of_workers
        else:
            options_list = [options.copy() for i in range(number_of_workers)]
            for individual, shared in zip(individual_options, options_list):
                shared.update(individual)
        
        
        threads = [ThreadWithResult(target=stellar_evolution_class, kwargs=options_list[i]) for i in range(number_of_workers)]
        self.code_instances = self._execute_all_threads(threads)
        
        self.particles = ParallelParticlesSuperset(
            [code.particles for code in self.code_instances], 
            execute_all_threads_func=self._execute_all_threads)
        self.parameters = ParallelParameters(self.code_instances)
    
    @option(type='boolean', sections=('code'))
    def must_run_threaded(self):
        return True
    
    def _execute_all_threads(self, threads):
        if self.must_run_threaded:
            return self._execute_all_threads_parallel(threads)
        else:
            return self._execute_all_threads_serial(threads)
            
    def _execute_all_threads_parallel(self, threads):
        for x in threads:
            x.start()
        for x in threads:
            x.join()
        result = [x.get_result() for x in threads]
        if not result == [None]*len(threads):
            return result


    def _execute_all_threads_serial(self, threads):
        for x in threads:
            x.run()
        result = [x.get_result() for x in threads]
        if not result == [None]*len(threads):
            return result
    
    def _run_threaded(self, function_name, args=()):
        threads = [ThreadWithResult(target=getattr(code, function_name), args=args) for code in self.code_instances]
        return self._execute_all_threads(threads)
    
    def initialize_code(self):
        self._run_threaded("initialize_code")
    
    def commit_parameters(self):
        self._run_threaded("commit_parameters")
    
    def recommit_parameters(self):
        self._run_threaded("recommit_parameters")
    
    def commit_particles(self):
        self._run_threaded("commit_particles")
    
    def recommit_particles(self):
        self._run_threaded("recommit_particles")
    
    def evolve_model(self, end_time):
        self._run_threaded("evolve_model", args=(end_time,))
        self.model_time = end_time
    
    def cleanup_code(self):
        self._run_threaded("cleanup_code")
    
    def stop(self):
        self._run_threaded("stop")
    
    def merge_colliding(self, *args, **kwargs):
        return merge_colliding_in_stellar_evolution_code(self, *args, **kwargs)
    
    def new_particle_from_model(self, internal_structure, current_age, key=None):
        index = self.particles.next_index_of_code_instance_for_new_particle_from_model()
        return self.code_instances[index].new_particle_from_model(internal_structure, current_age, key=key)


class ThreadWithResult(threading.Thread):
    
    def __init__(self, target=None, args=(), kwargs=dict()):
        self.__target = target
        self.__args = args
        self.__kwargs = kwargs
        self.result = None
        self.caught_exception = False
        threading.Thread.__init__(self, target=target, args=args, kwargs=kwargs)
    
    def run(self):
        try:
            self.result = self.__target(*self.__args, **self.__kwargs)
        except Exception:
            self.caught_exception = True
            self.result = sys.exc_info()
    
    def get_result(self):
        if self.caught_exception:
            raise self.result[1].with_traceback(self.result[2])
        return self.result
    

class ParallelParticlesSuperset(ParticlesSuperset):
    
    def __init__(self, particle_sets, execute_all_threads_func=None):
        ParticlesSuperset.__init__(self, particle_sets)
        self._private.number_of_particles = 0
        self._private.number_of_sets = len(particle_sets)
        self._private.execute_all_threads_func = execute_all_threads_func
    
    def add_particles_to_store(self, keys, attributes = [], values = []):
        slices = [slice((i-self._private.number_of_particles) % self._private.number_of_sets, len(keys), 
            self._private.number_of_sets) for i in range(self._private.number_of_sets)]
        threads = [
            ThreadWithResult(
                target=particle_set.add_particles_to_store, 
                args=(keys[one_slice], attributes, [v[one_slice] for v in values])
            ) for particle_set, one_slice in zip(self._private.particle_sets, slices) if len(keys[one_slice])]
        self._private.execute_all_threads_func(threads)
        self._private.number_of_particles += len(keys)
    
    def next_index_of_code_instance_for_new_particle_from_model(self):
        next_index = -self._private.number_of_particles % self._private.number_of_sets
        self._private.number_of_particles += 1
        return next_index
    
class ParallelParameters(object):
    
    def __init__(self, code_instances):
        object.__setattr__(self, "code_instances", code_instances)
    
    def __getattr__(self, attribute_name):
        return getattr(object.__getattribute__(self, "code_instances")[0].parameters, attribute_name)
    
    def __setattr__(self, attribute_name, value):
        for code in object.__getattribute__(self, "code_instances"):
            setattr(code.parameters, attribute_name, value)
    

