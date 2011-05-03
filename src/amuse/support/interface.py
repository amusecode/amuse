from amuse.support.data import parameters
from amuse.support.data import core
from amuse.support.data import values
from amuse.support.data import incode_storage
from amuse.support.units import nbody_system
from amuse.support.units import generic_unit_system
from amuse.support.units.core import unit
from amuse.support.options import OptionalAttributes

from amuse.support.methods import CodeMethodWrapper, CodeMethodWrapperDefinition, IncorrectWrappedMethodException

from amuse.support.core import late
from amuse.support import exceptions
from amuse.support import state

import inspect

class OldObjectsBindingMixin(object):

    def setup_particles(self, particles):
        self.particles.add_particles(particles)

    def update_particles(self, particles):
        self.particles.copy_values_of_state_attributes_to(particles)

    

    





class CodeAttributeWrapper(object):

    def __init__(self):
        pass

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



class HandleConvertUnits(HandleCodeInterfaceAttributeAccess, CodeMethodWrapperDefinition):

    def __init__(self, handler):
        self.handler = handler
        self.converter = None

    def supports(self, name, was_found):
        return was_found and not self.converter is None

    def get_attribute(self, name, attribute):
        if inspect.ismethod(attribute):
            result = attribute #UnitsConvertionMethod(attribute, self.converter)
        elif isinstance(attribute, core.AbstractParticleSet):
            result = attribute #core.ParticlesWithUnitsConverted(attribute, self.converter)
        elif isinstance(attribute, values.Quantity):
            result = self.converter.from_target_to_source(attribute)
        elif isinstance(attribute, CodeMethodWrapper):
            result = CodeMethodWrapper(attribute, self)
        elif isinstance(attribute, parameters.Parameters):
            result = parameters.ParametersWithUnitsConverted(attribute, self.converter)
        elif hasattr(attribute, '__iter__'):
            result = list(self.convert_and_iterate(attribute))
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
        object.define_converter(self)

    def convert_arguments(self, method,  list_arguments, keyword_arguments):
        converted_list_arguments = [self.from_source_to_target(x) if bool else x for x,bool in zip(list_arguments,method.nbody_input_attributes)]
        converted_keyword_arguments = {}
        for key, value in keyword_arguments.iteritems():
            converted_keyword_arguments[key] = self.from_source_to_target(value)

        return converted_list_arguments, converted_keyword_arguments

    def convert_result(self, method, result):
        return self.from_target_to_source(result)

    def from_source_to_target(self, x):
        if isinstance(x, values.Quantity):
            return self.converter.from_source_to_target(x)
        else:
            return x

    def from_target_to_source(self, x):
        if isinstance(x, values.Quantity):
            if x.unit.is_non_numeric():
                return x
            else:
                return self.converter.from_target_to_source(x)
        elif hasattr(x, '__len__'):
            if len(x) > 200:
                return x
            else:
                return list(self.convert_and_iterate(x))
        elif hasattr(x, '__iter__'):
            return list(self.convert_and_iterate(x))
        else:
            return x

class StateMethodDefinition(CodeMethodWrapperDefinition):
    def __init__(self, state_machine, interface, from_state, to_state, function_name):
        self.state_machine = state_machine
        self.interface = interface
        self.transitions = []
        self.add_transition(from_state, to_state)
        self.function_name = function_name

    def add_transition(self, from_state, to_state):
        self.transitions.append((from_state, to_state))

    def new_method(self, method = None):
        if method == None:
            method = getattr(self.interface, self.function_name)
        return CodeMethodWrapper(method, self)

    def precall(self, method):
        stored_transitions = []
        for from_state, to_state in self.transitions:
            if from_state is None:
                return to_state
            elif from_state == self.state_machine._current_state:
                return to_state
            else:
                stored_transitions.append((from_state, to_state))

        for from_state, to_state  in stored_transitions:
            try:
                self.state_machine._do_state_transition_to(from_state)
                return to_state
            except Exception, ex:
                pass

        # do again to get an exception.
        self.state_machine._do_state_transition_to(stored_transitions[0][0])

    def postcall(self, method, to_state):
        if to_state is None:
            return
        elif to_state == self.state_machine._current_state:
            return
        else:
            self.state_machine._current_state = to_state


class HandleState(HandleCodeInterfaceAttributeAccess):

    def __init__(self, interface):
        self._mapping_from_name_to_state_method = {}
        self.interface = interface
        self._state_machine = state.StateMachine(interface)


    def supports(self, name, was_found):
        if name == 'state_machine':
            return True
        else:
            return self._state_machine.is_enabled() and (name in self._mapping_from_name_to_state_method)

    def get_attribute(self, name, value): 
        if name == 'state_machine':
            return self._state_machine
        else:
            return self._mapping_from_name_to_state_method[name].new_method(value)


    def attribute_names(self):
        result = set(self._mapping_from_name_to_state_method.keys())
        result.add('state_machine')
        return result

    def define_state(self, name):
        self._state_machine.new_state(name)

    def _add_state_method(self, from_state, to_state, function_name):
        if not function_name in self._mapping_from_name_to_state_method:
            state_method = StateMethodDefinition(self._state_machine, self.interface, from_state, to_state, function_name)
            self._mapping_from_name_to_state_method[function_name] = state_method
        else:
            state_method = self._mapping_from_name_to_state_method[function_name]
            state_method.add_transition(from_state, to_state)



    def add_method(self, state_name, function_name):
        """
        Define a method that can run when the interface is in the
        provided state.
        """
        self._add_state_method( self._state_machine.new_state(state_name), None, function_name)


    def add_transition(self, from_name, to_name, function_name, is_auto = True):
    
        transition = self._state_machine.new_transition(from_name, to_name, is_auto)
    
        definition = StateMethodDefinition(self._state_machine, self.interface, transition.from_state, transition.to_state, function_name)
    
        transition.method = definition
    
        self._add_state_method(transition.from_state, transition.to_state, function_name)


    def add_transition_to_method(self, state_name, function_name, is_auto = True):
        """
        Define a method that can run in any state and will transition the interface
        to the provided state.
        """
        
        transition = self._state_machine.new_transition(None, state_name, is_auto)
        
    
        definition = StateMethodDefinition(self._state_machine, self.interface, transition.from_state, transition.to_state, function_name)
        transition.method = definition
    
        self._add_state_method(None, transition.to_state, function_name)


    def do_automatic_state_transitions(self, boolean):
        self._state_machine._do_automatic_state_transitions = boolean

    def set_initial_state(self, name):
        self._state_machine.set_initial_state(name)


    def setup(self, object):
        object.define_state(self)


    def has_name(self, name):
        return name == 'STATE'




    def get_name_of_current_state(self):
        return self._state_machine.get_name_of_current_state()
    
    
class MethodWithUnits(CodeMethodWrapper):

    def __init__(self, original_method, definition):
        CodeMethodWrapper.__init__(self, original_method, definition)

    @late
    def index_input_attributes(self):
        return self.definition.index_input_attributes
    
    @late
    def nbody_input_attributes(self):
        return self.definition.nbody_input_attributes

    @late
    def index_output_attributes(self):
        return self.definition.index_output_attributes



class MethodWithUnitsDefinition(CodeMethodWrapperDefinition):

    ERROR_CODE =  object()
    NO_UNIT = object()
    INDEX = object()

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

    def check_wrapped_method(self, method):
        if method.method_is_legacy or method.method_is_code:
            self.check_outputs_of_method(method)
            self.check_inputs_of_method(method)
        
    def new_method(self, original_method):
        if self.has_same_name_as_original:
            return MethodWithUnits(original_method, self)
        else:
            return MethodWithUnits(getattr(self.handler.interface, self.function_name), self)

    def handle_errorcode(self, errorcode):
        if errorcode in self.handler.interface.errorcodes:
            raise exceptions.AmuseException("Error when calling '{0}' of a '{1}', errorcode is {2}, error is '{3}'".format(self.name, type(self.handler.interface).__name__, errorcode,  self.handler.interface.errorcodes[errorcode]))
        elif errorcode < 0:
            raise exceptions.AmuseException("Error when calling '{0}' of a '{1}', errorcode is {2}".format(self.name, type(self.handler.interface).__name__, errorcode))
        else:
            return errorcode

    def handle_as_errorcode(self, errorcode):
        if hasattr(errorcode, 'any'):
            if not errorcode.any():
                return
        if hasattr(errorcode, '__iter__'):
            for x in errorcode:
                self.handle_errorcode(x)
        else:
            self.handle_errorcode(errorcode)

    def handle_as_unit(self, return_value):
        if not hasattr(self.return_units, '__iter__'):
            if self.return_units == self.NO_UNIT or self.return_units == self.INDEX:
                return return_value
            elif self.return_units == self.ERROR_CODE:
                self.handle_as_errorcode(return_value)
            else:
                return self.return_units.new_quantity(return_value)
        else:
            if not hasattr(return_value, '__iter__'):
                return_value = [return_value]
            result = []
            for value, unit in zip(return_value, self.return_units):
                if unit == self.ERROR_CODE:
                    self.handle_as_errorcode(value)
                elif unit == self.NO_UNIT:
                    result.append(value)
                elif unit == self.INDEX:
                    result.append(value)
                else:
                    result.append(unit.new_quantity(value))
            if len(result) == 1:
                return result[0]
            else:
                return result


    def convert_arguments(self, method, list_arguments, keyword_arguments):
        result = {}
        input_parameters = method.method_input_argument_names

        for index, parameter in enumerate(input_parameters):
            if parameter in keyword_arguments:
                if self.units[index] == self.NO_UNIT or self.units[index] == self.INDEX:
                    result[parameter] = keyword_arguments[parameter]
                else:
                    result[parameter] = keyword_arguments[parameter].value_in(self.units[index])

        for index, argument in enumerate(list_arguments):
            parameter = input_parameters[index]
            if self.units[index] == self.NO_UNIT or self.units[index] == self.INDEX:
                result[parameter] = argument
            else:
                if self.units[index].is_none() and not hasattr(argument,'unit'):
                    result[parameter] = argument
                else:
                    result[parameter] = argument.value_in(self.units[index])

        return (), result

    def convert_result(self, method, result):
        return self.handle_return_value(result)

    @late
    def has_same_name_as_original(self):
        return self.function_name == self.name


    @late
    def index_input_attributes(self):
        return map(lambda x : x == self.INDEX, self.units)

    @late
    def nbody_input_attributes(self):
        return map(lambda x : isinstance(x, unit) and generic_unit_system.is_generic_unit(x), self.units)

    @late
    def index_output_attributes(self):
        if not hasattr(self.return_units, '__iter__'):
            return [self.return_units == self.INDEX]
        else:
            return map(lambda x : x == self.INDEX, self.return_units)


    def check_inputs_of_method(self, method):
    
        specification = method.legacy_specification
        if specification is None:
            return
            
        number_expected_inputs  = len(specification.input_parameters)
            
        if self.units:
            if hasattr(self.units, '__len__'):
                number_specified_inputs = len(self.units)
            else:
                number_specified_inputs = 1
        else:
            number_specified_inputs = 0
        
        if number_expected_inputs != number_specified_inputs:
            raise IncorrectMethodDefinition(self.name, type(self.handler.interface).__name__, number_expected_inputs, number_specified_inputs, 'inputs')
    
    

    def check_outputs_of_method(self, method):
    
        specification = method.legacy_specification
        if specification is None:
            return
            
        number_expected_outputs  = len(specification.output_parameters)
        if specification.result_type != None:
            number_expected_outputs  += 1
            
        if self.return_units:
            if hasattr(self.return_units, '__len__'):
                number_specified_outputs = len(self.return_units)
            else:
                number_specified_outputs = 1
        else:
            number_specified_outputs = 0
        
        if number_expected_outputs == 1 and  number_specified_outputs == 0:
            return#defualt error checks for one output
            
        if number_expected_outputs != number_specified_outputs:
            raise IncorrectMethodDefinition(self.name, type(self.handler.interface).__name__, number_expected_outputs, number_specified_outputs, 'outputs')
    
    
class HandleMethodsWithUnits(object):
    ERROR_CODE = MethodWithUnitsDefinition.ERROR_CODE
    NO_UNIT = MethodWithUnitsDefinition.NO_UNIT
    INDEX = MethodWithUnitsDefinition.INDEX

    def __init__(self, interface):
        self.method_definitions = {}
        self.interface = interface

    def supports(self, name, was_found):
        return name in self.method_definitions

    def get_attribute(self, name, value):
        return self.method_definitions[name].new_method(value)


    def attribute_names(self):
        return set(self.method_definitions.keys())

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
        object.define_methods(self)



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
            return_value = function_or_attribute()
            if hasattr(return_value, '__iter__'):
                if len(return_value) > 2:
                    return_value = list(return_value)
                    value, errorcode = return_value[:-1], return_value[-1]
                else:
                    value, errorcode = return_value
                if errorcode < 0:
                    raise exceptions.AmuseException("calling '{0}' to get the value for property '{1}' resulted in an error (errorcode {2})".format(self.function_or_attribute_name, self.public_name, errorcode))
                else:
                    return self.unit.new_quantity(value)
            else:
                return self.unit.new_quantity(return_value)
        else:
            return self.unit.new_quantity(function_or_attribute)


    @late
    def has_same_name_as_original(self):
        return self.function_or_attribute_name == self.public_name



class PropertyDefinition(object):

    def __init__(self, handler, function_or_attribute_name, public_name):
        self.function_or_attribute_name = function_or_attribute_name
        self.public_name = public_name
        self.handler = handler

    def get_value(self, original):
        if self.has_same_name_as_original:
            function_or_attribute = original
        else:
            function_or_attribute = getattr(self.handler.interface, self.function_or_attribute_name)
    
        if hasattr(function_or_attribute, '__call__'):
            return_value = function_or_attribute()
            return return_value
        else:
            return function_or_attribute


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

    def attribute_names(self):
        return set(self.property_definitions.keys())

    def add_property(self, function_name, unit = None, public_name = None):
        if public_name is None:
            if function_name.startswith('get_'):
                public_name = function_name[4:]
            else:
                public_name = function_name

        if unit is None:
            definition = PropertyDefinition(self, function_name, public_name)
        else:
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
        object.define_properties(self)


class HandleParameters(HandleCodeInterfaceAttributeAccess):
    def __init__(self, interface):
        self.property_definitions = {}
        self.interface = interface
        self.definitions = []
        self.parameters = None

    def supports(self, name, was_found):
        return name == 'parameters'

    def get_attribute(self, name, value):
        if not self.parameters:
            self.parameters =  parameters.Parameters(self.definitions, self.interface)
       
        return self.parameters

    def attribute_names(self):
        return set(['parameters'])

    def add_method_parameter(self, get_method, set_method, name, description, unit = None, default_value = None,must_set_before_get = False):
        definition = parameters.ModuleMethodParameterDefinition(
            get_method,
            set_method,
            name,
            description,
            default_value,
            must_set_before_get = must_set_before_get
        )
        self.definitions.append(definition)


    def add_caching_parameter(self, function_name, parameter_name, name, description, default_value = None):
        definition = parameters.ModuleCachingParameterDefinition(
            function_name,
            parameter_name,
            name,
            description,
            default_value
        )
        self.definitions.append(definition)

    def add_boolean_parameter(self, get_method, set_method, name, description, default_value = None):
        definition = parameters.ModuleBooleanParameterDefinition(
            get_method,
            set_method,
            name,
            description,
            default_value
        )
        self.definitions.append(definition)

    def has_name(self, name):
        return name == 'PARAMETER'

    def setup(self, object):
        object.define_parameters(self)



    def add_vector_parameter(self, name, description, parameter_names):
        definition = parameters.VectorParameterDefinition(
            name,
            description,
            parameter_names,
        )
        self.definitions.append(definition)
    
    
class HandleErrorCodes(HandleCodeInterfaceAttributeAccess):
    def __init__(self, interface):
        self.error_codes = {}
        self.interface = interface

    def supports(self, name, was_found):
        return name == 'errorcodes'

    def get_attribute(self, name, value):
        return self.error_codes

    def attribute_names(self):
        return set(['errorcodes'])

    def add_errorcode(self, number, string):
        self.error_codes[number] = string

    def has_name(self, name):
        return name == 'ERRORCODE'

    def setup(self, object):
        object.define_errorcodes(self)


class AbstractParticleSetDefinition(object):
    
    def set_new(self, name_of_new_particle_method, names = None):
        self.new_particle_method = (name_of_new_particle_method, names)
        
    def set_grid_range(self, name_of_the_get_range_method):
        self.name_of_the_get_range_method = name_of_the_get_range_method

    def set_delete(self, name_of_delete_particle_method):
        self.name_of_delete_particle_method = name_of_delete_particle_method

    def add_getter(self, name_of_the_getter, names = None):
        self.getters.append((name_of_the_getter, names))

    def add_setter(self, name_of_the_setter, names = None):
        self.setters.append((name_of_the_setter, names))

    def add_attribute(self, name_of_the_attribute, name_of_the_method, names = None):
        self.attributes.append((name_of_the_attribute,name_of_the_method, names))

    def add_query(self, name_of_the_query, names = (), public_name = None):
        if not public_name:
            public_name = name_of_the_query
        self.queries.append((name_of_the_query, names, public_name))


    def add_method(self, name_of_the_method, public_name = None):
        if not public_name:
            public_name = name_of_the_method
        self.methods.append((name_of_the_method, public_name))


    def add_select_from_particle(self, name, names = (), public_name = None):
        if not public_name:
            public_name = name
        self.selects_form_particle.append((name, names, public_name))
    
    def define_extra_keywords(self, dictionary):
        self.extra_keyword_arguments_for_getters_and_setters = dictionary
    
    

    def add_subselect_in_set(self, name, set_query_arguments_name = None, get_number_of_particles_name = None,  public_name = None):
        if not public_name:
            public_name = name
        self.subselects_in_set.append((name, set_query_arguments_name, get_number_of_particles_name, public_name))
    

    
class ParticleSetDefinition(AbstractParticleSetDefinition):

    def __init__(self, handler):
        self.handler = handler
        self.name_of_indexing_attribute = 'index_of_the_particle'
        self.new_particle_method = ('new_particle',(), None)
        self.name_of_delete_particle_method = 'delete_particle'
        self.name_of_number_of_particles_method = 'get_number_of_particles'
        self.setters = []
        self.getters = []
        self.queries = []
        self.attributes = []
    
        self.selects_form_particle = []
        self.subselects_in_set = []
        self.methods = []
        self.is_superset = False
        self.is_inmemory = False
        self.particles_factory = core.Particles

    def new_storage(self, interface):
    
        if self.is_inmemory:
            return core.get_in_memory_attribute_storage_factory()()
    
        setters = []
        for name, names in self.setters:
            x = incode_storage.ParticleSetAttributesMethod(getattr(interface, name), names)
            setters.append(x)
    
        getters = []
        for name, names in self.getters:
            x = incode_storage.ParticleGetAttributesMethod(getattr(interface, name), names)
            getters.append(x)
    
    
        name, names = self.new_particle_method
        new_particle_method = incode_storage.NewParticleMethod(getattr(interface, name), names)
    
        delete_particle_method = getattr(interface, self.name_of_delete_particle_method)
        number_of_particles_method = None#getattr(interface, self.name_of_number_of_particles_method)
    
        return incode_storage.InCodeAttributeStorage(
            interface,
            new_particle_method,
            delete_particle_method,
            number_of_particles_method,
            setters,
            getters,
            self.name_of_indexing_attribute
        )

    
    def new_set_instance(self, handler):
        storage = self.new_storage(handler.interface)
        if self.is_inmemory:
            result = self.particles_factory(handler.interface, storage = storage)
        else:
            result = self.particles_factory(storage = storage)
            
        queries = self.new_queries(handler.interface)
        for x in queries:
            result.add_function_attribute(x.public_name, x.apply)
    
        selects = self.new_selects_from_particle(handler.interface)
        for x in selects:
            result.add_function_attribute(x.public_name, x.apply_on_all)
            result.add_particle_function_attribute(x.public_name, x.apply_on_one)
            
        selects = self.new_subselects_in_set(handler.interface)
        for x in selects:
            result.add_function_attribute(x.public_name, x.apply_on_all)
    
        selects = self.new_particle_methods(handler.interface)
        for x in selects:
            result.add_function_attribute(x.public_name, x.apply_on_all, x.apply_on_one)
    
        attributes = self.attributes
        for name_of_the_attribute, name_of_the_method, names in attributes:
            result.add_calculated_attribute(name_of_the_attribute, getattr(handler.interface, name_of_the_method), names)
    
        return result

    def new_queries(self, interface):
        queries = []
        for name, names, public_name in self.queries:
            x = incode_storage.ParticleQueryMethod(getattr(interface, name), names, public_name)
            queries.append(x)

        return queries

    def new_selects_from_particle(self, interface):
        results = []
        for name, names, public_name in self.selects_form_particle:
            x = incode_storage.ParticleSpecificSelectMethod(getattr(interface, name), names, public_name)
            results.append(x)

        return results

    def new_particle_methods(self, interface):
        results = []
        for name, public_name in self.methods:
            x = incode_storage.ParticleMethod(getattr(interface, name), public_name)
            results.append(x)

        return results


    def new_subselects_in_set(self, interface):
        results = []
        for name, set_query_arguments_name, number_of_particles_name, public_name in self.subselects_in_set:
            number_of_particles_method = None if number_of_particles_name is None else getattr(interface, number_of_particles_name)
            set_query_arguments_method = None if set_query_arguments_name is None else getattr(interface, set_query_arguments_name)
            x = incode_storage.ParticleSetSelectSubsetMethod(getattr(interface, name), set_query_arguments_method, number_of_particles_method, public_name)
            results.append(x)
    
        return results
    
    
class ParticleSupersetDefinition(AbstractParticleSetDefinition):

    def __init__(self, handler, particle_subset_names, index_to_default_set=None):
        self.handler = handler
        self.particle_subset_names = particle_subset_names
        self.index_to_default_set = index_to_default_set
        self.is_superset = True
        self.particles_factory = core.ParticlesSuperset
        
    
    def new_set_instance(self, handler):
        subsets = [handler.get_attribute(subset_name, None) for subset_name in self.particle_subset_names]
        return self.particles_factory(
            subsets, 
            index_to_default_set=self.index_to_default_set
        )
    
class GridDefinition(AbstractParticleSetDefinition):

    def __init__(self, handler):
        self.handler = handler
        self.name_of_the_get_range_method = 'get_range'
        self.setters = []
        self.getters = []
        self.particles_factory = core.Grid
        self.extra_keyword_arguments_for_getters_and_setters = {}

    def new_storage(self, interface):

        setters = []
        for name, names in self.setters:
            x = incode_storage.ParticleSetAttributesMethod(getattr(interface, name), names)
            setters.append(x)

        getters = []
        for name, names in self.getters:
            x = incode_storage.ParticleGetAttributesMethod(getattr(interface, name), names)
            getters.append(x)

        range_method = getattr(interface, self.name_of_the_get_range_method)

        return incode_storage.InCodeGridAttributeStorage(
            interface,
            range_method,
            setters,
            getters,
            self.extra_keyword_arguments_for_getters_and_setters
        )

    
    
    def new_set_instance(self, handler):
        storage = self.new_storage(handler.interface)
        result = self.particles_factory(storage = storage)
        return result

class CodeInMemoryParticles(core.Particles):

    def __init__(self, code_interface = None, storage = None):
        core.Particles.__init__(self, storage = storage)
        self._private.code_interface = code_interface

class HandleParticles(HandleCodeInterfaceAttributeAccess):
    def __init__(self, interface):
        self.interface = interface
        self.mapping_from_name_to_set_definition = {}
        self.mapping_from_name_to_set_instance = {}


    def supports(self, name, was_found):
        return name in self.mapping_from_name_to_set_definition

    def get_attribute(self, name, value):
        if name in self.mapping_from_name_to_set_instance:
            return self.mapping_from_name_to_set_instance[name]
        else:
            set_definition = self.mapping_from_name_to_set_definition[name]
            result = set_definition.new_set_instance(self)
                
            self.mapping_from_name_to_set_instance[name] = result
            return result

    def attribute_names(self):
        return set(self.mapping_from_name_to_set_definition.keys())


    def has_name(self, name):
        return name == 'PARTICLES'

    def setup(self, object):
        object.define_particle_sets(self)

    def define_set(self, name, name_of_indexing_attribute = 'index_of_the_particle'):
        definition = ParticleSetDefinition(self)
        definition.name_of_indexing_attribute = name_of_indexing_attribute
        self.mapping_from_name_to_set_definition[name] = definition

    def define_super_set(self, name, particle_subsets, index_to_default_set = None):
        definition = ParticleSupersetDefinition(self, particle_subsets, index_to_default_set)
        self.mapping_from_name_to_set_definition[name] = definition

    def define_inmemory_set(self, name, particles_factory = CodeInMemoryParticles):
        definition = ParticleSetDefinition(self)
        definition.is_inmemory = True
        definition.particles_factory = particles_factory
        self.mapping_from_name_to_set_definition[name] = definition
    
    def define_grid(self, name, name_of_indexing_attribute = 'index_of_the_particle'):
        definition = GridDefinition(self)
        definition.name_of_indexing_attribute = name_of_indexing_attribute
        self.mapping_from_name_to_set_definition[name] = definition
        
    def set_new(self, name_of_the_set, name_of_new_particle_method, names = None):
        self.mapping_from_name_to_set_definition[name_of_the_set].set_new(name_of_new_particle_method, names = names)
        
    def set_grid_range(self, name_of_the_set, name_of_the_get_range_method):
        self.mapping_from_name_to_set_definition[name_of_the_set].set_grid_range(name_of_the_get_range_method)

    def set_delete(self, name_of_the_set, name_of_delete_particle_method):
        self.mapping_from_name_to_set_definition[name_of_the_set].set_delete(name_of_delete_particle_method)

    def add_getter(self, name_of_the_set, name_of_the_getter, names = None):
        self.mapping_from_name_to_set_definition[name_of_the_set].add_getter(name_of_the_getter, names = names)

    def add_setter(self, name_of_the_set, name_of_the_setter, names = None):
        self.mapping_from_name_to_set_definition[name_of_the_set].add_setter(name_of_the_setter, names = names)

    def add_attribute(self, name_of_the_set, name_of_the_attribute, name_of_the_method, names = None):
        self.mapping_from_name_to_set_definition[name_of_the_set].add_attribute(name_of_the_attribute, name_of_the_method, names = names)

    def add_query(self, name_of_the_set, name_of_the_query, names = (), public_name = None):
        self.mapping_from_name_to_set_definition[name_of_the_set].add_query(name_of_the_query, names = names, public_name = public_name)

    def add_method(self, name_of_the_set, name_of_the_method, public_name = None):
        self.mapping_from_name_to_set_definition[name_of_the_set].add_method(name_of_the_method, public_name = public_name)

    def add_select_from_particle(self, name_of_the_set, name, names = (), public_name = None):
        self.mapping_from_name_to_set_definition[name_of_the_set].add_select_from_particle(name, names = names, public_name = public_name)
        
    def define_extra_keywords(self, name_of_the_set, dictionary):
        self.mapping_from_name_to_set_definition[name_of_the_set].define_extra_keywords(dictionary)    
    
    def add_subselect_in_set(self, name_of_the_set, name, set_query_arguments_name = None, get_number_of_particles_name = None,  public_name = None):
        self.mapping_from_name_to_set_definition[name_of_the_set].add_subselect_in_set(
            name, 
            set_query_arguments_name = set_query_arguments_name, 
            get_number_of_particles_name = get_number_of_particles_name,  
            public_name = public_name
        )
        
class OverriddenCodeInterface(object):

    def __init__(self, code_interface):
        self.code_interface = code_interface

    def __getattr__(self, name):
        return self.code_interface.__getattr__(name)

class InCodeComponentImplementation(OldObjectsBindingMixin, OptionalAttributes):

    def __init__(self, legacy_interface, **options):
        OptionalAttributes.__init__(self, **options)
        self.legacy_interface = legacy_interface
        self._handlers = []


        self._handlers.append(LegacyInterfaceHandler(legacy_interface))
        self._handlers.append(HandleMethodsWithUnits(self))
        self._handlers.append(HandlePropertiesWithUnits(self))
        self._handlers.append(HandleParameters(self))
        self._handlers.append(HandleParticles(self))
        self._handlers.append(HandleState(self))
        self._handlers.append(HandleConvertUnits(self))
        self._handlers.append(HandleErrorCodes(self))

        self.setup()

    def setup(self):
        for x in self._handlers:
            x.setup(self)

    def define_state(self, handler):
        pass

    def define_methods(self, handler):
        pass

    def define_properties(self, handler):
        pass

    def define_converter(self, handler):
        pass

    def define_parameters(self, handler):
        pass

    def define_particle_sets(self, handler):
        pass

    def define_errorcodes(self, handler):
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
        result = set(dir(type(self)))
        result |= set(self.__dict__.keys())
        for handler in self._handlers:
            result |= handler.attribute_names()
        return list(result)

    def overridden(self):
        return OverriddenCodeInterface(self)

    def get_name_of_current_state(self):
        return self.state_machine.get_name_of_current_state()
        
    def _create_new_grid(self, builder_function, **extra_arguments):
        handler = self.get_handler('PARTICLES')
        definition = GridDefinition(handler)
        builder_function(definition, **extra_arguments)
        return definition.new_set_instance(handler)


class IncorrectMethodDefinition(IncorrectWrappedMethodException):
    formatstring = "Incorrect definition of method '{0}' of class '{1}', the number of {4} do not match, expected {2}, actual {3}."


class PropertyDefinition(object):

    def __init__(self, handler, functionname, publicname, keyword_arguments = {}):
        self.functionname = functionname
        self.publicname = publicname
        self.handler = handler
        self.keyword_arguments = {}

    def get_value(self, original):
        method = getattr(self.handler.interface, self.functionname)
        return method(**self.keyword_arguments)

