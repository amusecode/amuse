from amuse.units import units, generic_unit_system
from amuse.units import nbody_system as nbody
from amuse.rfi.core import legacy_function
from amuse.rfi.core import LegacyFunctionSpecification

class StoppingConditionInterface(object):

    @legacy_function
    def has_stopping_condition():
        """
        Return 1 if the stopping condition with
        the given index is supported by the code,
        0 otherwise.
        """
        function = LegacyFunctionSpecification()  
        function.can_handle_array = True
        function.addParameter('type', dtype='int32', direction=function.IN, description = "The type index  of the stopping condition")
        function.addParameter('result', dtype='int32', direction=function.OUT, description = "1 if the stopping condition is supported")
        
        function.result_type = 'int32'
        function.result_doc = """
        0 - OK
        -1 - ERROR
        """
        return function 
        
    @legacy_function
    def enable_stopping_condition():
        """
        Will enable the stopping if it is supported
        """
        function = LegacyFunctionSpecification()   
        function.can_handle_array = True
        function.addParameter('type', dtype='int32', direction=function.IN, description = "The type index of the stopping condition")
        function.result_type = 'int32'
        function.result_doc = """
        0 - OK
        -1 - ERROR
        """
        return function 
        
    
    @legacy_function
    def disable_stopping_condition():
        """
        Will disable the stopping if it is supported
        """
        function = LegacyFunctionSpecification()  
        function.can_handle_array = True 
        function.addParameter('type', dtype='int32', direction=function.IN, description = "The index of the stopping condition")
        function.result_type = 'int32'
        function.result_doc = """
        0 - OK
        -1 - ERROR
        """
        return function 
        
    @legacy_function
    def is_stopping_condition_enabled():
        """
        Return 1 if the stopping condition with
        the given index is enabled,0 otherwise.    
        """
        function = LegacyFunctionSpecification()   
        function.can_handle_array = True
        function.addParameter('type', dtype='int32', direction=function.IN, description = "The index of the stopping condition")
        function.addParameter('result', dtype='int32', direction=function.OUT, description = "1 if the stopping condition is enabled")
        function.result_type = 'int32'
        function.result_doc = """
        0 - OK
        -1 - ERROR
        """
        return function 
        
    @legacy_function
    def is_stopping_condition_set():
        """
        Return 1 if the stopping condition with
        the given index is enabled,0 otherwise.    
        """
        function = LegacyFunctionSpecification()  
        function.can_handle_array = True 
        function.addParameter('type', dtype='int32', direction=function.IN, description = "The index of the stopping condition")
        function.addParameter('result', dtype='int32', direction=function.OUT, description = "1 if the stopping condition is enabled")
        function.result_type = 'int32'
        function.result_doc = """
        0 - OK
        -1 - ERROR
        """
        return function 
        
    
    @legacy_function
    def get_number_of_stopping_conditions_set():
        """
        Return the number of stopping conditions set, one
        condition can be set multiple times.   
        
        Stopping conditions are set when the code determines
        that the conditions are met. The objects or or information
        about the condition can be retrieved with
        the get_stopping_condition_info method.
        """
        function = LegacyFunctionSpecification()  
        function.can_handle_array = True 
        function.addParameter('result', dtype='int32', direction=function.OUT, description = "> 1 if any stopping condition is set")
        function.result_type = 'int32'
        function.result_doc = """
        0 - OK
        -1 - ERROR
        """
        return function
        
    @legacy_function
    def get_stopping_condition_info():
        """
        Generic function for getting the information connected to
        a stopping condition. Index can be between 0 and
        the result of the :method:`get_number_of_stopping_conditions_set`
        method.
        """
        function = LegacyFunctionSpecification()  
        function.can_handle_array = True 
        function.addParameter('index', dtype='int32', direction=function.IN, description = "Index in the array[0,number_of_stopping_conditions_set>")
        function.addParameter('type', dtype='int32', direction=function.OUT, description = "Kind of the condition, can be used to retrieve specific information")
        function.addParameter('number_of_particles', dtype='int32', direction=function.OUT, description = "Number of particles that met this condition")
        function.result_type = 'int32'
        function.result_doc = """
        0 - OK
        -1 - ERROR
        """
        return function
        
    @legacy_function
    def get_stopping_condition_particle_index():
        """
        For collision detection
        """
        function = LegacyFunctionSpecification()  
        function.can_handle_array = True 
        function.addParameter('index', dtype='int32', direction=function.IN, description = "Index in the array[0,number_of_stopping_conditions_set>")
        function.addParameter('index_of_the_column', dtype='int32', direction=function.IN, description = "Column index involved in the condition (for pair collisons 0 and 1 are possible)")
        function.addParameter('index_of_particle', dtype='int32', direction=function.OUT, description = "Set to the identifier of particle[index_of_the_column][index]")
        function.result_type = 'int32'
        function.result_doc = """
        0 - OK
        -1 - ERROR
        """
        return function
    
    @legacy_function
    def set_stopping_condition_timeout_parameter():
        """
        Set max computer time available (in seconds).
        """
        function = LegacyFunctionSpecification()  
        function.can_handle_array = False 
        function.addParameter('value', dtype='float64', direction=function.IN, description = "Available wallclock time in seconds")
        function.result_type = 'int32'
        function.result_doc = """
        0 - OK
        -1 - Value out of range
        """
        return function
        
    @legacy_function
    def get_stopping_condition_timeout_parameter():
        """
        Retrieve max computer time available (in seconds).
        """
        function = LegacyFunctionSpecification()  
        function.can_handle_array = False 
        function.addParameter('value', dtype='float64', direction=function.OUT, description = "Current value of avaible wallclock time in seconds")
        function.result_type = 'int32'
        function.result_doc = """
        0 - OK
        -1 - ERROR
        """
        return function

    @legacy_function
    def set_stopping_condition_number_of_steps_parameter():
        """
        Set max inner loop evaluations.
        """
        function = LegacyFunctionSpecification()
        function.can_handle_array = False
        function.addParameter('value', dtype='int32', direction=function.IN, description = "Available inner loop evaluations")
        function.result_type = 'int32'
        function.result_doc = """
        0 - OK
        -1 - Value out of range
        """
        return function

    @legacy_function
    def get_stopping_condition_number_of_steps_parameter():
        """
        Retrieve max inner loop evaluations.
        """
        function = LegacyFunctionSpecification()  
        function.can_handle_array = False 
        function.addParameter('value', dtype='int32', direction=function.OUT, description = "Current number of avaible inner loop evaluations")
        function.result_type = 'int32'
        function.result_doc = """
        0 - OK
        -1 - ERROR
        """
        return function

    @legacy_function
    def set_stopping_condition_out_of_box_parameter():
        """
        Set size of box.
        """
        function = LegacyFunctionSpecification()
        function.can_handle_array = False
        function.addParameter('value', dtype='float64', direction=function.IN, description = "Size of box")
        function.result_type = 'int32'
        function.result_doc = """
        0 - OK
        -1 - Value out of range
        """
        return function

    @legacy_function
    def get_stopping_condition_out_of_box_parameter():
        """
        Get size of box
        """
        function = LegacyFunctionSpecification()
        function.can_handle_array = False
        function.addParameter('value', dtype='float64', direction=function.OUT, description = "Size of box")
        function.result_type = 'int32'
        function.result_doc = """
        0 - OK
        -1 - Value out of range
        """
        return function
    
    @legacy_function
    def set_stopping_condition_minimum_density_parameter():
        function = LegacyFunctionSpecification()
        function.addParameter('value', dtype='float64', direction=function.IN)
        function.result_type = 'int32'
        return function
    @legacy_function
    def get_stopping_condition_minimum_density_parameter():
        function = LegacyFunctionSpecification()
        function.addParameter('value', dtype='float64', direction=function.OUT)
        function.result_type = 'int32'
        return function
    
    @legacy_function
    def set_stopping_condition_maximum_density_parameter():
        function = LegacyFunctionSpecification()
        function.addParameter('value', dtype='float64', direction=function.IN)
        function.result_type = 'int32'
        return function
    @legacy_function
    def get_stopping_condition_maximum_density_parameter():
        function = LegacyFunctionSpecification()
        function.addParameter('value', dtype='float64', direction=function.OUT)
        function.result_type = 'int32'
        return function
    
    @legacy_function
    def set_stopping_condition_minimum_internal_energy_parameter():
        function = LegacyFunctionSpecification()
        function.addParameter('value', dtype='float64', direction=function.IN)
        function.result_type = 'int32'
        return function
    @legacy_function
    def get_stopping_condition_minimum_internal_energy_parameter():
        function = LegacyFunctionSpecification()
        function.addParameter('value', dtype='float64', direction=function.OUT)
        function.result_type = 'int32'
        return function
    
    @legacy_function
    def set_stopping_condition_maximum_internal_energy_parameter():
        function = LegacyFunctionSpecification()
        function.addParameter('value', dtype='float64', direction=function.IN)
        function.result_type = 'int32'
        return function
    @legacy_function
    def get_stopping_condition_maximum_internal_energy_parameter():
        function = LegacyFunctionSpecification()
        function.addParameter('value', dtype='float64', direction=function.OUT)
        function.result_type = 'int32'
        return function
    
class StoppingCondition(object):
    
    def __init__(self, conditions, type, description):
        self.conditions = conditions
        self.type = type
        self.description = description
        self.__doc__ = description
        
    def enable(self):
        self.conditions.code.enable_stopping_condition(self.type)
        
    def disable(self):
        self.conditions.code.disable_stopping_condition(self.type)
        
    def is_enabled(self):
        return self.conditions.code.is_stopping_condition_enabled(self.type) == 1
        
    def is_supported(self):
        return self.conditions.code.has_stopping_condition(self.type) == 1

    def is_set(self):
        return self.conditions.code.is_stopping_condition_set(self.type) == 1

    def get_set_condition_indices(self, index_in_condition):
        indices = range(self.conditions.code.get_number_of_stopping_conditions_set())
        if len(indices) == 0:
            return []
        types, number_of_particles = self.conditions.code.get_stopping_condition_info(indices)
    
        result = []
        for index, type, number_of_particles_in_condition in zip(indices, types, number_of_particles):
            if type == self.type and index_in_condition < number_of_particles_in_condition:
                result.append(index)
        return result
    
    def particles(self, index_in_the_condition=0, particles_set_name="particles"):
        selected = self.get_set_condition_indices(index_in_the_condition)
        
        if len(selected) == 0:
            return self.conditions.code.__getattr__(particles_set_name)[0:0]
        else:
            return self.conditions.code.__getattr__(particles_set_name).get_stopping_condition_particle_index(
                selected, [index_in_the_condition]*len(selected))
        
class StoppingConditions(object):

    def __init__(self, code):
        self.code = code
        self.collision_detection = StoppingCondition(
            self,
            0, 
            "If enabled, the code will stop at the end of the inner loop when two stars connect"
        )
        self.pair_detection = StoppingCondition(
            self, 
            1, 
            "If enabled, the code will stop at the end of the inner loop when two stars are bound"
        )
        self.escaper_detection = StoppingCondition(
            self, 
            2, 
            "If enabled, the code will stop at the end of the inner loop when a star escapes"
        )
        self.timeout_detection = StoppingCondition(
            self, 
            3, 
            "If enabled, the code will stop at the end of the inner loop when the computer time is above a set timeout"
        )
        self.number_of_steps_detection = StoppingCondition(
            self,
            4,
            "If enabled, the code will stop at the end of the inner loop when the number of evaluations reached the set max number"
        )
        self.out_of_box_detection = StoppingCondition(
            self,
            5,
            "If enabled, the code will stop if a particle escapes the box of size out_of_box_size"
        )
        self.density_limit_detection = StoppingCondition(
            self,
            6,
            "If enabled, the code will stop if a gas particle has a density out of the range "
                "[stopping_condition_minimum_density, stopping_condition_maximum_density]"
        )
        self.internal_energy_limit_detection = StoppingCondition(
            self,
            7,
            "If enabled, the code will stop if a gas particle has an internal energy out of the range "
                "[stopping_condition_minimum_internal_energy, stopping_condition_maximum_internal_energy]"
        )

    def all_conditions(self):
        for name in dir(self):
            if name.startswith("_"):
                continue
            else:
                value = getattr(self, name)
                if isinstance(value, StoppingCondition):
                    yield name, value
                    
    def __str__(self):
        parts = []
        parts.append("Stopping conditions of a '{0}' object\n".format(type(self.code).__name__))
        supported = self.supported_conditions()
        enabled = [name for name, condition in self.all_conditions() if condition.is_enabled()]
        hit = [name for name, condition in self.all_conditions() if condition.is_set()]
        parts.append('* supported conditions: ')
        parts.append(', '.join(supported))
        parts.append('\n')
        parts.append('* enabled conditions: ')
        if enabled:
            parts.append(', '.join(enabled))
        else:
            parts.append('none')
        parts.append('\n')
        parts.append('* set conditions: ')
        if hit:
            parts.append(', '.join(hit))
        else:
            parts.append('none')
        parts.append('\n')
        return ''.join(parts)
    
    def supported_conditions(self):
        return [name for name, condition in self.all_conditions() if condition.is_supported()]
        
    def define_parameters(self, object):
        object.add_method_parameter(
            "get_stopping_condition_timeout_parameter",
            "set_stopping_condition_timeout_parameter", 
            "stopping_conditions_timeout", 
            "max wallclock time available for the evolve step", 
            default_value = 4.0 |  units.s
        )

        object.add_method_parameter(
            "get_stopping_condition_number_of_steps_parameter",
            "set_stopping_condition_number_of_steps_parameter", 
            "stopping_conditions_number_of_steps", 
            "max inner loop evals", 
            default_value = 1.0
        )

        object.add_method_parameter(
            "get_stopping_condition_out_of_box_parameter",
            "set_stopping_condition_out_of_box_parameter", 
            "stopping_conditions_out_of_box_size", 
            "size of cube", 
            default_value = 0.0 |  nbody.length
        )
        
        object.add_method_parameter(
            "get_stopping_condition_minimum_density_parameter",
            "set_stopping_condition_minimum_density_parameter", 
            "stopping_condition_minimum_density", 
            "minimum density of a gas particle", 
            default_value = -1.0 | generic_unit_system.density
        )
        
        object.add_method_parameter(
            "get_stopping_condition_maximum_density_parameter",
            "set_stopping_condition_maximum_density_parameter", 
            "stopping_condition_maximum_density", 
            "maximum density of a gas particle", 
            default_value = -1.0 | generic_unit_system.density
        )
        
        object.add_method_parameter(
            "get_stopping_condition_minimum_internal_energy_parameter",
            "set_stopping_condition_minimum_internal_energy_parameter", 
            "stopping_condition_minimum_internal_energy", 
            "minimum internal energy of a gas particle", 
            default_value = -1.0 | generic_unit_system.specific_energy
        )
        
        object.add_method_parameter(
            "get_stopping_condition_maximum_internal_energy_parameter",
            "set_stopping_condition_maximum_internal_energy_parameter", 
            "stopping_condition_maximum_internal_energy", 
            "maximum internal energy of a gas particle", 
            default_value = -1.0 | generic_unit_system.specific_energy
        )
        
    def define_methods(self, object):
        object.add_method(
            'get_stopping_condition_particle_index',
            (   
                object.NO_UNIT,
                object.NO_UNIT,
            ),
            (
                object.INDEX,
                object.ERROR_CODE,
            )
        )
        
        object.add_method(
            'has_stopping_condition',
            (   
                object.NO_UNIT,
            ),
            (
                object.NO_UNIT,
                object.ERROR_CODE,
            )
        )
        object.add_method(
            'is_stopping_condition_enabled',
            (   
                object.NO_UNIT,
            ),
            (
                object.NO_UNIT,
                object.ERROR_CODE,
            )
        )
        object.add_method(
            'is_stopping_condition_set',
            (   
                object.NO_UNIT,
            ),
            (
                object.NO_UNIT,
                object.ERROR_CODE,
            )
        )
        object.add_method(
            'get_stopping_condition_info',
            (   
                object.NO_UNIT,
            ),
            (
                object.NO_UNIT,
                object.NO_UNIT,
                object.ERROR_CODE,
            )
        )
        object.add_method(
            'get_number_of_stopping_conditions_set',
            (   
            ),
            (
                object.NO_UNIT,
                object.ERROR_CODE,
            )
        )
        object.add_method(
            'enable_stopping_condition',
            ( object.NO_UNIT,),
            (
                object.ERROR_CODE
            )
        )
        object.add_method(
            'disable_stopping_condition',
            ( object.NO_UNIT,),
            (
                object.ERROR_CODE
            )
        )
        
        object.add_method(
            "get_stopping_condition_timeout_parameter",
            (),
            (units.s, object.ERROR_CODE,)
        )
        
        object.add_method(
            "set_stopping_condition_timeout_parameter",
            (units.s, ),
            (object.ERROR_CODE,)
        )
        
        object.add_method(
            "get_stopping_condition_number_of_steps_parameter",
            (),
            (object.NO_UNIT, object.ERROR_CODE,)
        )
        
        object.add_method(
            "set_stopping_condition_number_of_steps_parameter",
            (object.NO_UNIT, ),
            (object.ERROR_CODE,)
        )
        
        object.add_method(
            "get_stopping_condition_out_of_box_parameter",
            (),
            (nbody.length, object.ERROR_CODE,)
        )
        
        object.add_method(
            "set_stopping_condition_out_of_box_parameter",
            (nbody.length, ),
            (object.ERROR_CODE,)
        )
        
        object.add_method("get_stopping_condition_minimum_density_parameter", 
            (), (generic_unit_system.density, object.ERROR_CODE,))
        object.add_method("set_stopping_condition_minimum_density_parameter", 
            (generic_unit_system.density, ), (object.ERROR_CODE,))
        
        object.add_method("get_stopping_condition_maximum_density_parameter", 
            (), (generic_unit_system.density, object.ERROR_CODE,))
        object.add_method("set_stopping_condition_maximum_density_parameter", 
            (generic_unit_system.density, ), (object.ERROR_CODE,))
        
        object.add_method("get_stopping_condition_minimum_internal_energy_parameter", 
            (), (generic_unit_system.specific_energy, object.ERROR_CODE,))
        object.add_method("set_stopping_condition_minimum_internal_energy_parameter", 
            (generic_unit_system.specific_energy, ), (object.ERROR_CODE,))
        
        object.add_method("get_stopping_condition_maximum_internal_energy_parameter", 
            (), (generic_unit_system.specific_energy, object.ERROR_CODE,))
        object.add_method("set_stopping_condition_maximum_internal_energy_parameter", 
            (generic_unit_system.specific_energy, ), (object.ERROR_CODE,))
        
    def define_particle_set(self, object, name_of_the_set = 'particles'):
        object.add_query(name_of_the_set, 'get_stopping_condition_particle_index')
