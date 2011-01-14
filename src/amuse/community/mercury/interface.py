from amuse.community import *

class Mercury(LegacyInterface):
    def __init__(self, **args):
        LegacyInterface.__init__(self, name_of_the_worker = 'mercury_worker',**args)

    @legacy_function
    def initialize_code():
        function = LegacyFunctionSpecification()  
        function.result_type = 'i'
        return function

    @legacy_function   
    def cleanup_code():
        function = LegacyFunctionSpecification()  
        function.result_type = 'i'
        return function

    @legacy_function    
    def commit_particles():
        function = LegacyFunctionSpecification()  
        function.result_type = 'i'
        return function

    @legacy_function    
    def evolve():
        function = LegacyFunctionSpecification()  
        function.addParameter('time', dtype='d', direction=function.IN)
        function.result_type = 'i'
        return function

    @legacy_function    
    def get_time():
        function = LegacyFunctionSpecification()  
        function.addParameter('time', dtype='d', direction=function.OUT)
        function.result_type = 'i'
        return function

    @legacy_function          
    def delete_particle():
        function = LegacyFunctionSpecification()  
        function.addParameter('id', dtype='i', direction=function.IN)
        function.result_type = 'i'
        return function    

    @legacy_function    
    def new_orbiter():
        function = LegacyFunctionSpecification()  
        function.can_handle_array = True
        function.addParameter('id', dtype='i', direction=function.OUT)
        for x in ['mass','dens','x','y','z','vx','vy','vz','sx','sy','sz','celimit']:
            function.addParameter(x, dtype='d', direction=function.IN)
        function.result_type = 'i'
        return function

    @legacy_function    
    def get_orbiter_state():
        function = LegacyFunctionSpecification()   
        function.can_handle_array = True
        function.addParameter('id', dtype='i', direction=function.IN)
        for x in ['mass','dens','x','y','z','vx','vy','vz','sx','sy','sz','celimit']:
            function.addParameter(x, dtype='d', direction=function.OUT)
        function.result_type = 'i'
        return function

    @legacy_function    
    def get_central_mass():
        function = LegacyFunctionSpecification()   
        function.addParameter('mass', dtype='d', direction=function.OUT)
        function.result_type = 'i'
        return function

    @legacy_function    
    def set_central_mass():
        function = LegacyFunctionSpecification()   
        function.addParameter('mass', dtype='d', direction=function.IN)
        function.result_type = 'i'
        return function

    @legacy_function    
    def get_central_radius():
        function = LegacyFunctionSpecification()   
        function.addParameter('radius', dtype='d', direction=function.OUT)
        function.result_type = 'i'
        return function

    @legacy_function    
    def set_central_radius():
        function = LegacyFunctionSpecification()   
        function.addParameter('radius', dtype='d', direction=function.IN)
        function.result_type = 'i'
        return function

    @legacy_function    
    def get_central_oblateness():
        function = LegacyFunctionSpecification()   
        for x in ['j2','j4','j6']:
            function.addParameter(x, dtype='d', direction=function.OUT)
        function.result_type = 'i'
        return function

    @legacy_function    
    def set_central_oblateness():
        function = LegacyFunctionSpecification()   
        for x in ['j2','j4','j6']:
            function.addParameter(x, dtype='d', direction=function.IN)
        function.result_type = 'i'
        return function

    @legacy_function    
    def get_central_spin():
        function = LegacyFunctionSpecification()   
        for x in ['lx','ly','lz']:
            function.addParameter(x, dtype='d', direction=function.OUT)
        function.result_type = 'i'
        return function

    @legacy_function    
    def set_central_spin():
        function = LegacyFunctionSpecification()   
        for x in ['lx','ly','lz']:
            function.addParameter(x, dtype='d', direction=function.IN)
        function.result_type = 'i'
        return function

    @legacy_function    
    def get_number_of_orbiters():
        function = LegacyFunctionSpecification()   
        function.addParameter('norbiters', dtype='i', direction=function.OUT)
        function.result_type = 'i'
        return function

    @legacy_function    
    def get_energy_deviation():
        function = LegacyFunctionSpecification()   
        function.addParameter('delta_e', dtype='d', direction=function.OUT)
        function.result_type = 'i'
        return function

    @legacy_function    
    def get_kinetic_energy():
        function = LegacyFunctionSpecification()   
        function.addParameter('ek', dtype='d', direction=function.OUT)
        function.result_type = 'i'
        return function

    @legacy_function    
    def get_potential_energy():
        function = LegacyFunctionSpecification()   
        function.addParameter('ep', dtype='d', direction=function.OUT)
        function.result_type = 'i'
        return function

    @legacy_function    
    def get_total_energy():
        function = LegacyFunctionSpecification()   
        function.addParameter('e_tot', dtype='d', direction=function.OUT)
        function.result_type = 'i'
        return function

    @legacy_function    
    def get_total_angular_momentum():
        function = LegacyFunctionSpecification()   
        function.addParameter('am_tot', dtype='d', direction=function.OUT)
        function.result_type = 'i'
        return function

