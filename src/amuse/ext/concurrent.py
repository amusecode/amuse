import random
from amuse.rfi.channel import MpiChannel
from amuse.datamodel import Particles

from amuse.units.quantities import zero
from amuse.units import constants

import pickle

import numpy
try:
    from mpi4py import MPI
except ImportError:
    MPI = None
    
class ConcurrentProcesses(object):
    pass
    
class MPIConcurrentProcesses(object):
    
    ROOT = 0
    @classmethod
    def is_available(self):
        return not MPI is None
        
    def init(self):
        MpiChannel.ensure_mpi_initialized()
        self.mpi_comm = MPI.COMM_WORLD
        self.shared_particles_ids = set([])
    
    def share(self, particles = None):
        sendbuffer = numpy.zeros(1,  dtype='int64')
        if self.mpi_comm == self.ROOT:
            new_id = random.getrandbits(64)
            sendbuffer[0] = new_id
            
        self.mpi_comm.Bcast([sendbuffer, MPI.INTEGER8], root=self.ROOT)
        shared_id = sendbuffer[0]
        
        self.shared_particles_ids.add(shared_id)
            
        return MPISharedParticlesProxy(particles, shared_id, self)
        
    def is_on_root(self):
        return self.mpi_comm.rank == self.ROOT
        
    def on_root(self, callable):
        if self.mpi_comm.rank == self.ROOT:
            callable()
        
    def not_on_root(self, callable):
        if self.mpi_comm.rank != self.ROOT:
            callable()
    
    def call(self, on_root, not_on_root):
        if self.mpi_comm.rank == self.ROOT:
            on_root()
        else:
            not_on_root()

class MPISharedParticlesProxy(object):
    def __init__(self, particles, shared_id, concurrent_processes):
        self.shared_id = shared_id
        self.particles = particles
        self.concurrent_processes = concurrent_processes
        
    def __getattr__(self, name):
        return self.particles.__getattr__(name)
    
    def distribute(self):
        self.concurrent_processes.call(
            self.distribute_on_root,
            self.distribute_not_on_root
        )
    
    def distribute_on_root(self):
        attribute_names = self.particles.get_attribute_names_defined_in_store()
        
        values = self.particles.get_values_in_store(
            self.particles.get_all_indices_in_store(), 
            attribute_names
        )
        units = [x.unit for x in values]
        units_dump = pickle.dumps(units)
        attributes_dump = pickle.dumps(attribute_names)
        
        units_dump = numpy.fromstring(units_dump,dtype='uint8')
        attributes_dump = numpy.fromstring(attributes_dump,dtype='uint8')
        
        sendbuffer = numpy.zeros(4,  dtype='int64')
        sendbuffer[0] = self.shared_id
        sendbuffer[1] = len(self.particles)
        sendbuffer[2] = len(units_dump)
        sendbuffer[3] = len(attributes_dump)
        
        self.concurrent_processes.mpi_comm.Bcast([sendbuffer, MPI.INTEGER8], root=self.concurrent_processes.ROOT)
        
        sendbuffer = self.particles.key
        self.concurrent_processes.mpi_comm.Bcast([sendbuffer, MPI.INTEGER8], root=self.concurrent_processes.ROOT)
        
        attribute_names = self.particles.get_attribute_names_defined_in_store()
       
        self.concurrent_processes.mpi_comm.Bcast([units_dump, MPI.CHARACTER], root=self.concurrent_processes.ROOT)
        self.concurrent_processes.mpi_comm.Bcast([attributes_dump, MPI.CHARACTER], root=self.concurrent_processes.ROOT)
        
        for x, unit in zip(values, units):
            value = x.value_in(unit)
            self.concurrent_processes.mpi_comm.Bcast([value, MPI.DOUBLE], root=self.concurrent_processes.ROOT)
        
        
    def distribute_not_on_root(self):
        sendbuffer = numpy.zeros(4,  dtype='int64')
        self.concurrent_processes.mpi_comm.Bcast([sendbuffer, MPI.INTEGER8], root=self.concurrent_processes.ROOT)
        shared_id = sendbuffer[0]
        number_of_particles = sendbuffer[1]
        units_dump_len = sendbuffer[2]
        attributes_dump_len = sendbuffer[3]
        
        sendbuffer = numpy.zeros(number_of_particles,  dtype='int64')
        self.concurrent_processes.mpi_comm.Bcast([sendbuffer, MPI.INTEGER8], root=self.concurrent_processes.ROOT)
        
        units_dump = numpy.zeros(units_dump_len, dtype='uint8')
        self.concurrent_processes.mpi_comm.Bcast([units_dump, MPI.CHARACTER], root=self.concurrent_processes.ROOT)
        
        attributes_dump = numpy.zeros(attributes_dump_len, dtype='uint8')
        self.concurrent_processes.mpi_comm.Bcast([attributes_dump, MPI.CHARACTER], root=self.concurrent_processes.ROOT)
        
        units = pickle.loads(units_dump.tobytes())
        attributes = pickle.loads(attributes_dump.tobytes())
        values = []
        for x in units:
            value = numpy.zeros(number_of_particles, dtype='float64')
            self.concurrent_processes.mpi_comm.Bcast([value, MPI.DOUBLE], root=self.concurrent_processes.ROOT)
            values.append(x.new_quantity(value))
            
        self.particles = Particles(keys = sendbuffer)
        self.particles.set_values_in_store(self.particles.get_all_indices_in_store(), attributes, values)
    
    def potential_energy(self, smoothing_length_squared = zero, G = constants.G):
        mpi_comm = self.concurrent_processes.mpi_comm
        
        mass = self.mass
        x_vector = self.x
        y_vector = self.y
        z_vector = self.z

        sum_of_energies = zero
        
        number_of_particles = len(self)
        block_size = (number_of_particles - 1) // mpi_comm.size
        start = mpi_comm.rank * block_size
        if mpi_comm.rank == (mpi_comm.size - 1):
            block_size = (number_of_particles - 1) - start
            
        for i in range(start, start + block_size):
            x = x_vector[i]
            y = y_vector[i]
            z = z_vector[i]
            dx = x - x_vector[i+1:]
            dy = y - y_vector[i+1:]
            dz = z - z_vector[i+1:]
            dr_squared = (dx * dx) + (dy * dy) + (dz * dz)
            dr = (dr_squared+smoothing_length_squared).sqrt()
            m_m = mass[i] * mass[i+1:]

            energy_of_this_particle = (m_m / dr).sum()
            sum_of_energies -= energy_of_this_particle
        
        value = sum_of_energies.value_in(sum_of_energies.unit)
        # for not assume unit is the same accross processes, 
        # so units are not distributed!
        
        input = numpy.zeros(1,  dtype='float64')
        output = numpy.zeros(1,  dtype='float64')
        
        input[0] = value
            
        mpi_comm.Reduce(
            [input, MPI.DOUBLE], 
            [output, MPI.DOUBLE],
            op=MPI.SUM, 
            root=0
        )
        
        return G * sum_of_energies.unit.new_quantity(output[0])
        
    def __len__(self):
        return len(self.particles)
