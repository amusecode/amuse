import os.path
from mpi4py import MPI
import numpy

from amuse.legacy.support import core

from amuse.legacy.support.core import RemoteFunction, legacy_global
from amuse.support.units import nbody_system

class PhiGRAPE(object):            
    def __init__(self, convert_nbody = None):
        directory_of_this_module = os.path.dirname(__file__);
        full_name_of_the_worker = os.path.join(directory_of_this_module , 'muse_worker')
        
        self.intercomm = MPI.COMM_SELF.Spawn(full_name_of_the_worker, None,  1)
        self.channel = core.MpiChannel(self.intercomm)
        self.convert_nbody = convert_nbody
        
    def __del__(self):
        self.stop_worker()
        
    @core.legacy_function
    def stop_worker():
        function = RemoteFunction() 
        function.id = 0
        return function

    @core.legacy_function   
    def setup_module():
        function = RemoteFunction()  
        function.result_type = 'i'
        return function
        
    @core.legacy_function      
    def cleanup_module():
        function = RemoteFunction()  
        function.result_type = 'i'
        return function
    
    @core.legacy_function    
    def initialize_particles():
        function = RemoteFunction()  
        function.addParameter('time', dtype='d', direction=function.IN)
        function.result_type = 'i'
        return function;

    @core.legacy_function  
    def reinitialize_particles():
        function = RemoteFunction()  
        function.result_type = 'i'
        return function
                
    @core.legacy_function    
    def add_particle():
        function = RemoteFunction()  
        function.addParameter('id', dtype='i', direction=function.IN)
        for x in ['mass','radius','x','y','z','vx','vy','vz']:
            function.addParameter(x, dtype='d', direction=function.IN)
        function.result_type = 'i'
        return function

    @core.legacy_function    
    def evolve():
        function = RemoteFunction()  
        function.addParameter('time_end', dtype='d', direction=function.IN)
        function.addParameter('synchronize', dtype='i', direction=function.IN)
        function.result_type = 'i'
        return function
        
    @core.legacy_function   
    def get_number():
        function = RemoteFunction()  
        function.result_type = 'i'
        return function;
             
    @core.legacy_function    
    def get_state():
        function = RemoteFunction()  
        function.addParameter('id', dtype='i', direction=function.IN)
        for x in ['mass','radius','x','y','z','vx','vy','vz']:
            function.addParameter(x, dtype='d', direction=function.OUT)
        return function
        
    @core.legacy_function
    def set_mass():
        function = RemoteFunction()  
        function.result_type = 'i'
        function.addParameter('id', dtype='i', direction=function.IN)
        function.addParameter('mass', dtype='d', direction=function.IN)
        return function;

    @core.legacy_function      
    def get_time():
        function = RemoteFunction()
        function.result_type = 'd'
        return function

    @core.legacy_function      
    def get_time_step():
        function = RemoteFunction()  
        function.result_type = 'd'
        return function

    @core.legacy_function      
    def set_eps():
        function = RemoteFunction()  
        function.addParameter('eps2', dtype='d', direction=function.IN)
        return function

    @core.legacy_function      
    def set_eta():
        function = RemoteFunction()  
        function.addParameter('etas', dtype='d', direction=function.IN)
        function.addParameter('eta', dtype='d', direction=function.IN)
        return function

    @core.legacy_function      
    def get_kinetic_energy():
        function = RemoteFunction()  
        function.result_type = 'd'
        return function

    @core.legacy_function      
    def get_potential_energy():
        function = RemoteFunction()  
        function.result_type = 'd'
        return function

    @core.legacy_function      
    def get_energy_error():
        function = RemoteFunction()  
        function.result_type = 'd'
        return function

    @core.legacy_function      
    def find_colliding_secondary():
        function = RemoteFunction()  
        function.addParameter('id1', dtype='i', direction=function.IN)
        function.result_type = 'i'
        return function

        
    def add_star(self, star):
        id = star.id
        mass = self.convert_nbody.to_nbody(star.mass.value())
        position = self.convert_nbody.to_nbody(star.position.value())
        velocity = self.convert_nbody.to_nbody(star.velocity.value())
        
        mass = mass.number
        x = position.number[0]
        y = position.number[1]
        z = position.number[2]
        vx = velocity.number[0]
        vy = velocity.number[1]
        vz = velocity.number[2]
        radius = self.convert_nbody.to_nbody(star.radius.value()).number
        self.add_particle(id, mass, radius, x, y, z, vx, vy, vz)
        
    def update_star(self, star):
        state = self.get_state(star.id)
        time = self.convert_nbody.to_si( 0.0 | nbody_system.time)
        #star.mass.set_value_at_time(time, self.convert_nbody.to_si(nbody_system.mass(state.mass)))
        star.position.set_value_at_time(time, self.convert_nbody.to_si(nbody_system.length(numpy.array((state['x'], state['y'], state['z'])))))
        star.velocity.set_value_at_time(time, self.convert_nbody.to_si(nbody_system.speed(numpy.array((state['vx'], state['vy'], state['vz'])))))
        return star
         
    def evolve_model(self, time_end):
        result = self.evolve(self.convert_nbody.to_nbody(time_end).number, 0.0)
        return result
        
    def set_eps2(self, eps2):
        self.eps2 = self.convert_nbody.to_nbody(eps2).number
    
    
    def add_particles(self, particles):
        for x in particles:
            self.add_star(x)
            
    def update_particles(self, particles):
        for x in particles:
            self.update_star(x)
            
    def xadd_particles(self, particles):
        mass = []
        x_ = []
        y = []
        z = []
        vx = []
        vy = []
        vz = []
        radius = []
        id = []
        for x in particles:
            #self.update_star(x)
            id.append(x.id)
            mass_ = self.convert_nbody.to_nbody(x.mass.value())
            position = self.convert_nbody.to_nbody(x.position.value())
            velocity = self.convert_nbody.to_nbody(x.velocity.value())
            
            mass.append(mass_.number)
            x_.append(position.number[0])
            y.append(position.number[1])
            z.append(position.number[2])
            vx.append(velocity.number[0])
            vy.append(velocity.number[1])
            vz.append(velocity.number[2])
            radius.append(self.convert_nbody.to_nbody(x.radius.value()).number)
        self._add_particle(id, mass, radius, x_, y, z, vx, vy, vz)
            
    def update_attributes(self, attributes):
        for id, x in attributes:
            if x.name == 'mass':
                self.set_mass(id, self.convert_nbody.to_nbody(x.value()).number)
        
    
  
