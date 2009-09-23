
import numpy


from amuse.support.units import nbody_system
from amuse.legacy import *

class BHTree(LegacyInterface):
    include_headers = ['muse_dynamics.h', 'parameters.h', 'local.h']
    
    extra_content = """
    int _add_particle(int id, double mass, double radius, double x, double y, double z, double vx, double vy, double vz) {
        dynamics_state state;
        state.id = id;
        state.mass = mass;
        state.radius = radius;
        state.x = x;
        state.y = y;
        state.z = z;
        state.vx = vx;
        state.vy = vy;
        state.vz = vz;
        return add_particle(state);
    }
    
    void _get_state(int id, int * id_out,  double * mass, double * radius, double * x, double * y, double * z, double * vx, double * vy, double * vz) {
        dynamics_state state = get_state(id);
        *id_out = state.id;
        *mass = state.mass;
        *radius = state.radius;
        *x = state.x;
        *y = state.y;
        *z = state.z;
        *vx = state.vx;
        *vy = state.vy;
        *vz = state.vz;
    }
    """
    
    
    timestep = legacy_global(name='timestep',id=21,dtype='d')
    eps2_for_gravity = legacy_global(name='eps2_for_gravity',id=22,dtype='d')
    theta_for_tree = legacy_global(name='theta_for_tree',id=23,dtype='d')
    
    use_self_gravity = legacy_global(name='use_self_gravity',id=24,dtype='i')
    ncrit_for_tree = legacy_global(name='ncrit_for_tree',id=25,dtype='i')
    
    dt_dia = legacy_global(name='dt_dia',id=246,dtype='d')
    
    
            
    def __init__(self, convert_nbody = None):
        LegacyInterface.__init__(self)
        self.convert_nbody = convert_nbody

    @legacy_function   
    def setup_module():
        function = RemoteFunction() 
        function.result_type = 'i'
        return function
    
    
    @legacy_function      
    def cleanup_module():
        function = RemoteFunction()  
        function.result_type = 'i'
        return function
    
    @legacy_function    
    def initialize_particles():
        function = RemoteFunction() 
        function.addParameter('time', dtype='d', direction=function.IN)
        function.result_type = 'i'
        return function;
        
    @legacy_function    
    def add_particle():
        function = RemoteFunction()  
        function.name = '_add_particle'
        function.addParameter('id', dtype='i', direction=function.IN)
        for x in ['mass','radius','x','y','z','vx','vy','vz']:
            function.addParameter(x, dtype='d', direction=function.IN)
        function.result_type = 'i'
        return function
        
    @legacy_function    
    def get_state():
        function = RemoteFunction()  
        function.name = '_get_state'
        function.addParameter('id', dtype='i', direction=function.IN)
        function.addParameter('id_out', dtype='i', direction=function.OUT)
        for x in ['mass','radius','x','y','z','vx','vy','vz']:
            function.addParameter(x, dtype='d', direction=function.OUT)
        function.result_type = None
        return function
        
    @legacy_function    
    def evolve():
        function = RemoteFunction()  
        function.addParameter('time_end', dtype='d', direction=function.IN)
        function.addParameter('synchronize', dtype='i', direction=function.IN)
        function.result_type = 'i'
        return function
        
    @legacy_function  
    def reinitialize_particles():
        function = RemoteFunction()  
        function.result_type = 'i'
        return function
        
    @legacy_function   
    def get_number():
        function = RemoteFunction()  
        function.result_type = 'i'
        return function;
     
    @legacy_function
    def set_mass():
        function = RemoteFunction()  
        function.result_type = 'i'
        function.addParameter('id', dtype='i', direction=function.IN)
        function.addParameter('mass', dtype='d', direction=function.IN)
        return function;
    
    @legacy_function
    def get_time():
        function = RemoteFunction()  
        function.result_type = 'd'
        return function;
    
    @legacy_function      
    def get_kinetic_energy():
        function = RemoteFunction() 
        function.result_type = 'd'
        return function

    @legacy_function      
    def get_potential_energy():
        function = RemoteFunction()  
        function.result_type = 'd'
        return function
   
    def add_star(self, star):
        id = star.id
        mass = self.convert_nbody.to_nbody(star.mass.value())
        position = self.convert_nbody.to_nbody(star.position.value())
        velocity = self.convert_nbody.to_nbody(star.velocity.value())
        
        x = position.number[0]
        y = position.number[1]
        z = position.number[2]
        vx = velocity.number[0]
        vy = velocity.number[1]
        vz = velocity.number[2]
        radius = self.convert_nbody.to_nbody(star.radius.value()).number
        self.add_particle(id, mass.number, radius, x, y, z, vx, vy, vz)
        
    def update_star(self, star):
        state = self.get_state(star.id)
        time = self.convert_nbody.to_si( self.get_time() | nbody_system.time)
        #star.mass.set_value_at_time(time, self.convert_nbody.to_si(nbody_system.mass(state.mass)))
        star.position.set_value_at_time(time, self.convert_nbody.to_si(nbody_system.length(numpy.array((state['x'], state['y'], state['z'])))))
        star.velocity.set_value_at_time(time, self.convert_nbody.to_si(nbody_system.speed(numpy.array((state['vx'], state['vy'], state['vz'])))))
        return star
         
    def evolve_model(self, time_end):
        result = self.evolve(self.convert_nbody.to_nbody(time_end).number, 1)
        return result
        
    def set_eps2(self, eps2):
        self.eps2_for_gravity = self.convert_nbody.to_nbody(eps2).number
    
    
    def add_particles(self, particles):
        for x in particles:
            self.add_star(x)
            
    def update_particles(self, particles):
        for x in particles:
            self.update_star(x)
            
    def update_attributes(self, attributes):
        for id, x in attributes:
            if x.name == 'mass':
                self.set_mass(id, self.convert_nbody.to_nbody(x.value()).number)
        
    def get_energies(self):
        energy_unit = nbody_system.mass * nbody_system.length ** 2  * nbody_system.time ** -2
        kinetic_energy = self.get_kinetic_energy() | energy_unit
        potential_energy = self.get_potential_energy() | energy_unit
        return (self.convert_nbody.to_si(kinetic_energy), self.convert_nbody.to_si(potential_energy))
    
  
