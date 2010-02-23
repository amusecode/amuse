
from amuse.support.data import parameters
from amuse.support.data import binding
from amuse.support.data import core
from amuse.support.data import values
from amuse.support.data import code_particles
from amuse.support.units import nbody_system


from amuse.support.core import late

import inspect

class OldObjectsBindingMixin(object):
                
    def setup_particles(self, particles):
        self.particles.add_particles(particles)
        
    def update_particles(self, particles):
        self.particles.copy_values_of_state_attributes_to(particles)
    
    def set_attribute(self, attribute_name, particles):
        particles.copy_values_of_attribute_to(attribute_name, self.particles)
        
    def update_attribute(self, attribute_name, particles):
        self.particles.copy_values_of_attribute_to(attribute_name, particles) 
        
        
        
class CodeInterface(OldObjectsBindingMixin):
    """Base class of next level interfaces to codes
    """
    
    parameter_definitions = []
    
    def __init__(self):
        self.parameters = parameters.Parameters(self.parameter_definitions, self)
        
class CodeInterfaceWithConvertedUnits(OldObjectsBindingMixin):
    
            
    def __init__(self, interface, converter):
        self.converter = converter
        self.interface = interface
        
    def __getattr__(self, name):
        attribute = getattr(self.interface, name)
        if inspect.ismethod(attribute):
            result = UnitsConvertionMethod(attribute, self.converter)
            setattr(self, name, result)
            return result
        elif isinstance(attribute, core.AbstractParticleSet):
            result = core.ParticlesWithUnitsConverted(attribute, self.converter)
            setattr(self, name, result)
            return result
        elif isinstance(attribute, values.Quantity):
            return self.converter.from_target_to_source(attribute)
        elif isinstance(attribute, binding.BoundCodeMethod):
            result = UnitsConvertionMethod(attribute, self.converter)
            setattr(self, name, result)
            return result
        elif isinstance(attribute, parameters.Parameters):
            result = parameters.ParametersWithUnitsConverted(attribute, self.converter)
            setattr(self, name, result)
            return result
        elif hasattr(attribute, '__iter__'):
            return self.convert_and_iterate(attribute)
        else:
            return attribute
            
    def convert_and_iterate(self, iterable):
        for x in iterable:
            yield self.converter.from_target_to_source(x)
            
    def __dir__(self):
        return dir(self.interface)
        
class CodeInterfaceWithNBodyUnitsConverted(CodeInterfaceWithConvertedUnits):
    
    def __init__(self, interface, convert_nbody):
        if convert_nbody is None:
            convert_nbody = nbody_system.nbody_to_si.get_default()
        
        CodeInterfaceWithConvertedUnits.__init__(
            self,
            interface,
            convert_nbody.as_converter_from_si_to_nbody()
        )
        

class CodeAttributeWrapper(object):
    
    def __init__(self):
        pass
        
class CodeMethodWrapper(CodeAttributeWrapper):
    
    def __init__(self, method):
        CodeAttributeWrapper.__init__(self)
        self.method = method
      
    @late
    def method_is_legacy(self):
        return hasattr(self.method, 'specification')
        
    @late
    def method_is_code(self):
        return hasattr(self.method, 'method_input_argument_names')
    
    @late
    def method_input_argument_names(self):
        if self.method_is_code:
            return self.method.method_input_argument_names
        elif self.method_is_legacy:
            return map(lambda x : x.name , self.method.specification.input_parameters)
        else:
            args = inspect.getargspec(self.method).args
            if args:
                if args[0] == 'self' or args[0] == 'cls':
                    return args[1:]
            return args
      
    @late
    def method_output_argument_names(self):
        if self.method_is_code:
            return self.method.method_output_argument_names
        elif self.method_is_legacy:
            return map(lambda x : x.name , self.method.specification.output_parameters)
        else:
            return ()


class UnitsConvertionMethod(CodeMethodWrapper):
    def __init__(self, real_method, converter):
        CodeMethodWrapper.__init__(self, real_method)
        self.converter = converter
        
    def __call__(self, *list_arguments, **keyword_arguments):
        converted_list_arguments = [self.from_source_to_target(x) for x in list_arguments]
        converted_keyword_arguments = {}
        for key, value in keyword_arguments:
            converted_keyword_arguments[key] = self.from_source_to_target(value)
            
        result = self.method(*converted_list_arguments, **converted_keyword_arguments)
        
        return self.from_target_to_source(result)
        
    def from_source_to_target(self, x):
        if isinstance(x, values.Quantity):
            return self.converter.from_source_to_target(x)
        else:
            return x
            
    def from_target_to_source(self, x):
        if isinstance(x, values.Quantity):
            return self.converter.from_target_to_source(x)
        elif hasattr(x, '__iter__'):
            return list(self.convert_and_iterate(x))
        else:
            return x
    
    def convert_and_iterate(self, iterable):
        for x in iterable:
            yield self.converter.from_target_to_source(x)

        
class HandleCodeInterfaceAttributeAccess(object):
    
    def supports(self, name, was_found):
        return False
        
    def get_attribute(self, name, result):
        return result
    
    def attribute_names(self):
        return set([])
        
        
    def setup(self, object):
        pass
    
    def has_name(self, name):
        return False
        
class LegacyInterfaceHandler(HandleCodeInterfaceAttributeAccess):
    
    def __init__(self, legacy_interface):
        self.legacy_interface = legacy_interface
        
    def supports(self, name, was_found):
        return hasattr(self.legacy_interface, name)
        
    def get_attribute(self, name, result):
        return getattr(self.legacy_interface, name)
    
    def attribute_names(self):
        return set(dir(self.legacy_interface))
        
    
    def has_name(self, name):
        return name == 'LEGACY'
    
    
        
class HandleConvertUnits(HandleCodeInterfaceAttributeAccess):
    
    def __init__(self, handler):
        self.handler = handler
        self.converter = None
        
    def supports(self, name, was_found):
        return was_found and not self.converter is None
        
    def get_attribute(self, name, attribute):
        if inspect.ismethod(attribute):
            result = attribute #UnitsConvertionMethod(attribute, self.converter)
        elif isinstance(attribute, core.AbstractParticleSet):
            result = core.ParticlesWithUnitsConverted(attribute, self.converter)
        elif isinstance(attribute, values.Quantity):
            result = self.converter.from_target_to_source(attribute)
        elif isinstance(attribute, CodeMethodWrapper):
            result = UnitsConvertionMethod(attribute, self.converter)
        elif isinstance(attribute, parameters.Parameters):
            result = parameters.ParametersWithUnitsConverted(attribute, self.converter)
        elif hasattr(attribute, '__iter__'):
            result = self.convert_and_iterate(attribute)
        else:
            result = attribute
        return result
            
    def convert_and_iterate(self, iterable):
        for x in iterable:
            yield self.converter.from_target_to_source(x)
            
    
    def set_converter(self, converter):
        self.converter = converter
    
    def set_nbody_converter(self, nbody_converter):
        self.set_converter(nbody_converter.as_converter_from_si_to_nbody())
        
    def has_name(self, name):
        return name == 'UNIT'
        
    def setup(self, object):
        object.setup_converter(self)
    


class StateMethod(CodeMethodWrapper):
    def __init__(self, definition, method):
        CodeMethodWrapper.__init__(self, method)
        self.definition = definition
    
    
    def __call__(self, *list_arguments, **keyword_arguments):
        to_state = self.definition.precall()
        
        result = self.method(*list_arguments, **keyword_arguments)
        self.definition.postcall(to_state)
        return result
    
        
        

class StateMethodDefinition(object):
    def __init__(self, handler, from_state, to_state, function_name):
        self.handler = handler
        self.transitions = []
        self.add_transition(from_state, to_state)
        self.function_name = function_name
    
    def add_transition(self, from_state, to_state):
        self.transitions.append((from_state, to_state))
    
    def new_method(self, method = None):
        if method == None:
            method = getattr(self.handler.interface, self.function_name)
        return StateMethod(self, method)
        
    def precall(self):
        stored_transitions = []
        
        for from_state, to_state in self.transitions:
            if from_state is None:
                return to_state
            elif from_state == self.handler._current_state:
                return to_state
            else:
                stored_transitions.append((from_state, to_state))
        
        for from_state, to_state  in stored_transitions:
            try:
                self.handler._do_state_transition_to(from_state)
                return to_state
            except Exception, ex:
                pass
        
        # do again to get an exception.
        self.handler._do_state_transition_to(stored_transitions[0][0])
    
    def postcall(self, to_state):
        if to_state is None:
            return
        elif to_state == self.handler._current_state:
            return
        else:
            self.handler._current_state = to_state

            
class HandleState(HandleCodeInterfaceAttributeAccess):
    class State(object):
        def __init__(self, handler, name):
            self.handler = handler
            self.name = name
            self.from_transitions = []
            self.to_transitions = []
            
        def __str__(self):
            return "state '{0}'".format(self.name)
            
    class StateTransition(object):
        def __init__(self, handler, from_state, to_state, method = None, is_default = False):
            self.method = method
            self.to_state = to_state
            self.from_state = from_state
            self.is_default = is_default
            self.handler = handler
            
        def __str__(self):
            return "transition from {0} to {1}".format(self.from_state, self.to_state)
            
        def do(self):
            if self.method is None:
                self.handler.current_state = self.to_state
            else:
                self.method.new_method()() 
                
    
        
            
    def __init__(self, interface):
        self._mapping_from_name_to_state_method = {}
        self._current_state = None
        self._do_automatic_state_transitions = True
        self.states = {}
        self.interface = interface
        
    
    def supports(self, name, was_found):
        return name in self._mapping_from_name_to_state_method
        
    def get_attribute(self, name, value):
        return self._mapping_from_name_to_state_method[name].new_method(value)
        
        
    def define_state(self, name):
        if name in self.states:
            return        
        self.states[name] = self.State(self, name)
         
    def _do_state_transition_to(self, state):
        transitions = self._get_transitions_path_from_to(self._current_state, state)
        if transitions is None:
            raise Exception("No transition from current state {0} to {1} possible".format(self._current_state, state))
        
        transitions_with_methods = filter(lambda x : not x.method is None,transitions)
        if not self._do_automatic_state_transitions and len(transitions_with_methods) > 0:
            lines = []
            lines.append("Interface is not in {0}, should transition from {1} to {0} first.\n". format(state, self._current_state))
            for x in transitions:
                if x.method is None:
                    lines.append("{0}, automatic". format(x))
                else:
                    lines.append("{0}, calling '{1}'". format(x, x.method.function_name))
            exception = Exception('\n'.join(lines))
            exception.transitions = transitions
            raise exception
            
        for transition in transitions:
            transition.do()
        
    
    def _get_transitions_path_from_to(self, from_state, to_state):
        transitions = to_state.to_transitions
        
        paths = map(lambda x : [x], transitions)
        
        def has_no_circle(path):
            seen_states = set([])
            for transition in path:
                if transition.to_state in seen_states:
                    return False
                seen_states.add(transition.to_state)
            return True
                    
            
        while paths:
            current = paths.pop()
            first = current[0]
            if first.from_state == from_state:
                return current
            else:
                transitions = first.from_state.to_transitions
                new_paths = map(lambda x : [x], transitions)

                for new_path in new_paths:
                    new_path.extend(current)
            
                new_paths = filter(has_no_circle, new_paths)
            
                paths.extend(new_paths)
                
        return None
                
        
        
        
    def _add_state_method(self, from_state, to_state, function_name):
        if not function_name in self._mapping_from_name_to_state_method:
            state_method = StateMethodDefinition(self, from_state, to_state, function_name)
            self._mapping_from_name_to_state_method[function_name] = state_method
        else:
            state_method = self._mapping_from_name_to_state_method[function_name]
            state_method.add_transition(from_state, to_state)
            


    def add_method(self, state_name, function_name):
        """
        Define a method that can run when the interface is in the
        provided state.
        """
        self.define_state(state_name)
        
        state = self.states[state_name]
        
        self._add_state_method( state, None, function_name)
            
    
    def add_transition(self, from_name, to_name, function_name):
        self.define_state(from_name)
        self.define_state(to_name)
        
        from_state = self.states[from_name]
        to_state   = self.states[to_name]
        definition = StateMethodDefinition(self, from_state, to_state, function_name)
        
        transition = self.StateTransition(self, from_state, to_state, definition)
        
        from_state.from_transitions.append(transition)
        to_state.to_transitions.append(transition)
        
        self._add_state_method(from_state, to_state, function_name)
        
        
    def add_transition_to_method(self, state_name, function_name):
        """
        Define a method that can run in any state and will transition the interface
        to the provided state.
        """
        self.define_state(state_name)
        
        state = self.states[state_name]
        
        self._add_state_method(None, state, function_name)

        
    def do_automatic_state_transitions(self, boolean):
        self._do_automatic_state_transitions = boolean
        
    def set_initial_state(self, name):
        self.define_state(name)
        self._current_state = self.states[name]
        
    
    def setup(self, object):
        object.setup_state(self)
            
    
    def has_name(self, name):
        return name == 'STATE'
        


class MethodWithUnits(CodeMethodWrapper):

    def __init__(self, definition, original_method):
        CodeMethodWrapper.__init__(self, original_method)
        self.definition = definition
        
    def __call__(self, *list_arguments, **keyword_arguments):
        
        converted_keyword_arguments = self.definition.convert_list_and_keyword_arguments(self.method_input_argument_names, list_arguments, keyword_arguments)
        return_value = self.method(**converted_keyword_arguments)
        return self.definition.handle_return_value(return_value)
      
        
            
        
class MethodWithUnitsDefinition(object):
    
    ERROR_CODE =  object()
    NO_UNIT = object()

    def __init__(self, handler, function_name, units, return_units, return_value_handler, name):
        self.function_name = function_name
        
        if hasattr(units, '__iter__'):
            self.units = units
        else:
            self.units = (units,)
            
        self.return_units = return_units
        self.handler = handler
        self.name = name
        if return_units is None:
            if return_value_handler is None:
                self.handle_return_value = self.handle_as_errorcode
            else:
                self.handle_return_value = return_value_handler
        else:
            self.handle_return_value = self.handle_as_unit
    
    def new_method(self, original_method):
        if self.has_same_name_as_original:
            return MethodWithUnits(self, original_method)
        else:
            return MethodWithUnits(self, getattr(self.handler.interface, self.function_name))

    def handle_errorcode(self, errorcode):
        if errorcode < 0:
            raise Exception("Error when calling '{0}' of a '{1}', errorcode is {2}".format(self.name, type(self.handler.interface).__name__, errorcode))
        else:
            return errorcode 
            
    def handle_as_errorcode(self, errorcode):
        if hasattr(errorcode, '__iter__'):
            for x in errorcode:
                self.handle_errorcode(x)
        else:
            self.handle_errorcode(errorcode)     
    
    def handle_as_unit(self, return_value):
        if not hasattr(self.return_units, '__iter__'):
            if self.return_units == self.NO_UNIT:
                return return_value
            else:
                return self.return_units.new_quantity(return_value)
        else:
            result = []
            for value, unit in zip(return_value, self.return_units):
                if unit == self.ERROR_CODE:
                    self.handle_as_errorcode(value)
                elif unit == self.NO_UNIT:
                    result.append(value)
                else:
                    result.append(unit.new_quantity(value))
            if len(result) == 1:
                return result[0]
            else:
                return result
    
    
    def convert_list_and_keyword_arguments(self, input_parameters, list_arguments, keyword_arguments):
        result = {}
        
        
        for index, parameter in enumerate(input_parameters):
            if parameter in keyword_arguments:
                if self.units[index] == self.NO_UNIT:
                    result[parameter] = keyword_arguments[parameter.name]
                else:
                    result[parameter] = keyword_arguments[parameter.name].value_in(self.units[index])
        
        for index, argument in enumerate(list_arguments):
            parameter = input_parameters[index]
            if self.units[index] == self.NO_UNIT:
                result[parameter] = argument
            else:
                result[parameter] = argument.value_in(self.units[index])
        
        return result
    
    @late
    def has_same_name_as_original(self):
        return self.function_name == self.name
    
    

class HandleMethodsWithUnits(object):
    ERROR_CODE = MethodWithUnitsDefinition.ERROR_CODE
    NO_UNIT = MethodWithUnitsDefinition.NO_UNIT
    
    def __init__(self, interface):
        self.method_definitions = {}
        self.interface = interface

    def supports(self, name, was_found):
        return name in self.method_definitions
        
    def get_attribute(self, name, value):
        return self.method_definitions[name].new_method(value)
        
    def add_method(self, original_name, units, return_unit = None,  public_name = None, return_value_handler = None):
        if public_name is None:
            public_name = original_name
        
        definition = MethodWithUnitsDefinition(
            self,
            original_name, 
            units, 
            return_unit, 
            return_value_handler, 
            public_name
        )
        self.method_definitions[public_name] = definition
    
    def has_name(self, name):
        return name == 'METHOD'
        
    def setup(self, object):
        object.setup_methods(self)
    


class PropertyWithUnitsDefinition(object):
    
    def __init__(self, handler, function_or_attribute_name, unit, public_name):
        self.function_or_attribute_name = function_or_attribute_name
        self.unit = unit
        self.public_name = public_name
        self.handler = handler
    
    def get_value(self, original):
        if self.has_same_name_as_original:
            function_or_attribute = original
        else:
            function_or_attribute = getattr(self.handler.interface, self.function_or_attribute_name)
            
        if hasattr(function_or_attribute, '__call__'):
            value, errorcode = function_or_attribute()
            if errorcode < 0:
                raise Exception("calling '{0}' to get the value for property '{1}' resulted in an error (errorcode {2})".format(self.function_or_attribute_name, self.public_name, errorcode))
            else:
                return self.unit.new_quantity(value)
        else:
            return self.unit.new_quantity(function_or_attribute)
            
    
    @late
    def has_same_name_as_original(self):
        return self.function_or_attribute_name == self.public_name


class HandlePropertiesWithUnits(HandleCodeInterfaceAttributeAccess):
    def __init__(self, interface):
        self.property_definitions = {}
        self.interface = interface

    def supports(self, name, was_found):
        return name in self.property_definitions
        
    def get_attribute(self, name, value):
        return self.property_definitions[name].get_value(value)
        
    def add_property(self, function_name, unit, public_name = None):
        if public_name is None:
            if function_name.startswith('get_'):
                public_name = function_name[4:]
            else:
                public_name = function_name
                        
        definition = PropertyWithUnitsDefinition(
            self,
            function_name, 
            unit,
            public_name
        )
        self.property_definitions[public_name] = definition
    
    def has_name(self, name):
        return name == 'PROPERTY'
        
    def setup(self, object):
        object.setup_properties(self)
        

class HandleParameters(HandleCodeInterfaceAttributeAccess):
    def __init__(self, interface):
        self.property_definitions = {}
        self.interface = interface
        self.definitions = []

    def supports(self, name, was_found):
        return name == 'parameters'
        
    def get_attribute(self, name, value):
        return parameters.Parameters(self.definitions, self.interface.legacy_interface)
        
    def add_attribute_parameter(self, attribute_name, name, description, unit, default_value = None):
        definition = parameters.ModuleAttributeParameterDefinition(
            attribute_name,
            name, 
            description, 
            unit, 
            default_value
        )
        self.definitions.append(definition) 
    
    def add_method_parameter(self, get_method, set_method, name, description, unit, default_value = None):
        definition = parameters.ModuleMethodParameterDefinition_Next(
            get_method, 
            set_method, 
            name, 
            description, 
            unit, 
            default_value
        )
        self.definitions.append(definition) 
        
    def add_caching_parameter(self, parameter_name, name, description, unit, default_value = None):
        definition = ModuleCachingParameterDefinition(
            parameter_name, 
            name, 
            description, 
            unit, 
            default_value = None
        )
        self.definitions.append(definition) 
    
    def has_name(self, name):
        return name == 'PARAMETER'
        
    def setup(self, object):
        object.setup_parameters(self)
    
        
class ParticleSetDefinition(object):
    
    def __init__(self, handler):
        self.handler = handler
        self.name_of_indexing_attribute = 'index_of_the_particle'
        self.new_particle_method = ('new_particle',(), None)
        self.name_of_new_particle_method = 'delete_particle'
        self.setters = []
        self.getters = []
    
    def new_storage(self, interface):
        setters = []
        for name, mapping, names in self.setters:
            x = code_particles.ParticleSetAttributesMethod(getattr(interface, name), mapping, names)
            setters.append(x)
            
        getters = []
        for name, mapping, names in self.getters:
            x = code_particles.ParticleGetAttributesMethod(getattr(interface, name), mapping, names)
            getters.append(x)
            
        name, mapping, names = self.new_particle_method
        new_particle_method = code_particles.NewParticleMethod(getattr(interface, name), mapping, names)
        
        delete_particle_method = getattr(interface, self.name_of_new_particle_method)
        number_of_particles_method = getattr(interface, self.name_of_new_particle_method)
        
        return code_particles.InCodeAttributeStorage(
            interface,
            new_particle_method, 
            delete_particle_method, 
            number_of_particles_method, 
            setters,
            getters
        )
        

        
class HandleParticles(HandleCodeInterfaceAttributeAccess):
    def __init__(self, interface):
        self.interface = interface
        self.sets = {}
        self.particle_sets = {}
        

    def supports(self, name, was_found):
        return name in self.sets
        
    def get_attribute(self, name, value):
        if name in self.particle_sets:
            return self.particle_sets[name]
        else:
            storage = self.sets[name].new_storage(self.interface)
            result = core.Particles(storage = storage)
            self.particle_sets[name] = result
            return result
        
    def has_name(self, name):
        return name == 'PARTICLES'
        
    def setup(self, object):
        object.setup_all_particles(self)
     
    def define_set(self, name, name_of_indexing_attribute = 'index_of_the_particle'):
        definition = ParticleSetDefinition(self)
        definition.name_of_indexing_attribute = name_of_indexing_attribute
        self.sets[name] = definition

    def set_new(self, name_of_the_set, name_of_new_particle_method, mapping = (), names = None):
        self.sets[name_of_the_set].new_particle_method = (name_of_new_particle_method, mapping, names)
        
    def set_delete(self, name_of_the_set, name_of_delete_particle_method):
        self.sets[name_of_the_set].name_of_delete_particle_method = name_of_delete_particle_method
        
    def add_getter(self, name_of_the_set, name_of_the_getter, mapping = (), names = None):
        
        self.sets[name_of_the_set].getters.append((name_of_the_getter, mapping, names))

    def add_setter(self, name_of_the_set, name_of_the_setter, mapping = (), names = None):
        
        self.sets[name_of_the_set].setters.append((name_of_the_setter,mapping, names))

class CodeInterface2(OldObjectsBindingMixin):
    
    def __init__(self, legacy_interface):
        object.__setattr__(self, 'legacy_interface', legacy_interface)
        object.__setattr__(self, '_handlers', [])
         
        
        self._handlers.append(LegacyInterfaceHandler(legacy_interface))
        self._handlers.append(HandleMethodsWithUnits(self))
        self._handlers.append(HandlePropertiesWithUnits(self))
        self._handlers.append(HandleParameters(self))
        self._handlers.append(HandleParticles(self))
        self._handlers.append(HandleState(self))
        self._handlers.append(HandleConvertUnits(self))
        
        self.setup()
        
    def setup(self):
        for x in self._handlers:
            x.setup(self)

    def setup_state(self, handler):
        pass
    
    def setup_methods(self, handler):
        pass
        
    def setup_properties(self, handler):
        pass
        
    def setup_converter(self, handler):
        pass
        
    def setup_parameters(self, handler):
        pass
    
    def setup_all_particles(self, handler):
        pass
        
    def get_handler(self, name):
        for x in self._handlers:
            if x.has_name(name):
                return x
        return None
        
    def __getattr__(self, name):
        result = None
        found = False
        for handler in self._handlers:
            if handler.supports(name, found):
                result = handler.get_attribute(name, result)
                found = True
        if not found:
            raise AttributeError(name)
        return result
        
            
    def __dir__(self):
        result = set(dir(self))
        for handler in self.handlers:
            result |= handler.attribute_names()
        return result
        
