import numpy
import threading

from amuse.units import units
from amuse.datamodel import ParticlesSuperset
from amuse.support.exceptions import AmuseException

class ParallelStellarEvolution(object):

    def __init__(self, stellar_evolution_class, number_of_workers=1, **options):
        self.code_factory = stellar_evolution_class
        self.number_of_workers = number_of_workers
        self.model_time = 0.0 | units.Myr
        self.code_instances = [stellar_evolution_class(**options) for i in range(number_of_workers)]
        self.particles = ParallelParticlesSuperset([code.particles for code in self.code_instances])
    
    @property 
    def parameters(self):
        raise AmuseException("Not implemented for parallel stellar evolution")
    
    def run_threaded(self, function_name, args=()):
        threads=[]
        for code in self.code_instances:
            threads.append(threading.Thread(target=getattr(code, function_name), args=args))
        for x in threads:
            x.start()
        for x in threads:
            x.join()
    
    def initialize_code(self):
        self.run_threaded("initialize_code")
    
    def commit_parameters(self):
        self.run_threaded("commit_parameters")
    
    def recommit_parameters(self):
        self.run_threaded("recommit_parameters")
    
    def commit_particles(self):
        self.run_threaded("commit_particles")
    
    def recommit_particles(self):
        self.run_threaded("recommit_particles")
    
    def evolve_model(self, end_time):
        self.run_threaded("evolve_model", args=(end_time,))
        self.model_time = end_time
    
    def cleanup_code(self):
        self.run_threaded("cleanup_code")
    
    def stop(self):
        self.run_threaded("stop")
    

class ParallelParticlesSuperset(ParticlesSuperset):
    
    def __init__(self, particle_sets):
        ParticlesSuperset.__init__(self, particle_sets)
        self._private.number_of_particles = 0
        self._private.number_of_sets = len(particle_sets)
    
    def add_particles_to_store(self, keys, attributes = [], values = []):
        slices = [slice((i+self._private.number_of_particles) % self._private.number_of_sets, len(keys), 
            self._private.number_of_sets) for i in range(self._private.number_of_sets)]
        for particle_set, one_slice in zip(self._private.particle_sets, slices):
            particle_set.add_particles_to_store(keys[one_slice], attributes, [v[one_slice] for v in values])
        self._private.number_of_particles += len(keys)
    

