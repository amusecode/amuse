from amuse.units import nbody_system
from amuse.units import generic_unit_system
from amuse.units import quantities
from amuse.units.quantities import is_quantity
from amuse.units.core import unit
from amuse.units import core
from amuse.units import units
from amuse.support.options import OptionalAttributes, option

from amuse.support.methods import CodeMethodWrapper, CodeMethodWrapperDefinition, IncorrectWrappedMethodException
from amuse.support.methods import ProxyingMethodWrapper
from amuse.support.core import late
from amuse.support import exceptions
from amuse.support import state

import inspect
import numpy

from amuse import datamodel
from amuse.datamodel import base
from amuse.datamodel import parameters
from amuse.datamodel import incode_storage

import itertools
from collections import defaultdict

class ConvertArgumentsException(core.IncompatibleUnitsException):
    formatstring = "{0}"


class OldObjectsBindingMixin(object):

    def setup_particles(self, particles):
        self.particles.add_particles(particles)

    def update_particles(self, particles):
        self.particles.copy_values_of_all_attributes_to(particles)

class MethodArgumentOrResultType(object):
     _returns_result=True
     
     def append_result_value(self, method, definition, value, result):
         result.append(self.convert_result_value(method, definition, value))
         
     def convert_result_value(self, method, definition, value):
         return value
         
     def convert_argument_value(self, method, definition, value):
         return value

class NoUnitMethodArgumentOrResultType(MethodArgumentOrResultType):
    def __reduce__(self):
        return (_get_result_type, ("NO_UNIT",))
        
class UnitMethodArgumentOrResultType(MethodArgumentOrResultType):
    def __reduce__(self):
        return (_get_result_type, ("UNIT",))

class ErrorCodeMethodArgumentOrResultType(MethodArgumentOrResultType):
    _returns_result=False

    def append_result_value(self, method, definition, value, result):
        self.convert_result_value(method, definition, value)
        
    def convert_result_value(self, method, definition, errorcode):
        if hasattr(errorcode, 'any'):
            if not errorcode.any():
                return
        if hasattr(errorcode, '__iter__'):
            for x in errorcode:
                definition.handle_errorcode(x)
        else:
            definition.handle_errorcode(errorcode)
            
    def __reduce__(self):
        return (_get_result_type, ("ERROR_CODE",))
        


class IndexMethodArgumentOrResultType(MethodArgumentOrResultType):
     
    def convert_result_value(self, method, definition, value):
        return value
     
    def convert_argument_value(self, method, definition, value):
        return value
        
    def __reduce__(self):
        return (_get_result_type, ("INDEX",))

def _get_result_type(name):
    if name == "NO_UNIT":
        return MethodWithUnitsDefinition.NO_UNIT
    elif name == "UNIT":
        return MethodWithUnitsDefinition.UNIT
    elif name == "ERROR_CODE":
        return MethodWithUnitsDefinition.ERROR_CODE
    elif name == "INDEX":
        return MethodWithUnitsDefinition.INDEX
        
    
class LinkMethodArgumentOrResultType(MethodArgumentOrResultType):
    
    def __init__(self, linked_set_name):
        self.linked_set_name = linked_set_name
    
    def get_linked_set(self, method, definition):
        # method might provide a shorter path to the interface object
        # and is cleaner as definition might move to interface class later
        try:
            return getattr(definition.wrapped_object, self.linked_set_name)
        except Exception as ex:
            print ex
            raise ex
        
    def convert_result_value(self, method, definition, value):
        linked_set = self.get_linked_set(method, definition)
        storage = linked_set._private.attribute_storage
        keys = storage._get_keys_for_indices_in_the_code(value)
        result = base.LinkedArray(numpy.empty(len(keys), dtype= numpy.object))
        for index in range(len(keys)):
            key = keys[index]
            if key >= 0:
                result[index] = linked_set._get_particle(key)
        return result
     
    def convert_argument_value(self, method, definition, value):
        linked_set = self.get_linked_set(method, definition)
        storage = linked_set._private.attribute_storage
        if isinstance(value, datamodel.Particle):
            indices = storage.get_indices_of([value.key])
            return indices
        else:
            value = base.LinkedArray(numpy.asanyarray(value))
            value = value.as_set()
            valid = value.get_valid_particles_mask()
            all_keys = value.get_all_keys_in_store()
            valid_keys = all_keys[valid]
            valid_indices = numpy.asarray(storage.get_indices_of(valid_keys))
            all_indices = -1 * numpy.ones(len(value), dtype=valid_indices.dtype)
            all_indices[valid] = valid_indices
            return all_indices

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
        self.method_instances = {}

    def supports(self, name, was_found):
        return name in self.method_instances or hasattr(self.legacy_interface, name)

    def get_attribute(self, name, result):
        
        if not name in self.method_instances:
            attr = getattr(self.legacy_interface, name)      
            if hasattr(attr, '__call__'):
                self.method_instances[name] = ProxyingMethodWrapper(self.legacy_interface, name)
            else:
                return attr

        return self.method_instances[name]

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
        elif isinstance(attribute, datamodel.AbstractParticleSet):
            result = attribute #datamodel.ParticlesWithUnitsConverted(attribute, self.converter)
        elif isinstance(attribute, datamodel.AbstractGrid):
            result = attribute 
        elif isinstance(attribute, quantities.Quantity):
            result = self.converter.from_target_to_source(attribute)
        elif isinstance(attribute, CodeMethodWrapper):
            result = CodeMethodWrapper(attribute, self)
        elif isinstance(attribute, parameters.Parameters):
            result = parameters.new_parameters_with_units_converted_instance_with_docs(attribute, self.converter)
        elif isinstance(attribute, basestring):
            result = attribute
        elif isinstance(attribute, bytearray):
            result = attribute
        elif hasattr(attribute, '__iter__'):
            result = list(self.convert_and_iterate(attribute))
        else:
            result = attribute
        return result

    def convert_and_iterate(self, iterable):
        for x in iterable:
            if isinstance(x, quantities.Quantity):
                yield self.converter.from_target_to_source(x)
            else:
                yield x


    def set_converter(self, converter):
        self.converter = converter

    def set_nbody_converter(self, nbody_converter):
        self.set_converter(nbody_converter.as_converter_from_si_to_generic())

    def has_name(self, name):
        return name == 'UNIT'

    def setup(self, object):
        object.define_converter(self)

    def convert_arguments(self, method,  list_arguments, keyword_arguments):
        
        converted_list_arguments = []
        for x,is_nbody in zip(list_arguments, method.nbody_input_attributes):
            if is_nbody:
                converted_list_arguments.append(self.from_source_to_target(x))
            else:
                converted_list_arguments.append(x)
        
        converted_keyword_arguments = {}
        for key, value in keyword_arguments.iteritems():
            converted_keyword_arguments[key] = self.from_source_to_target(value)

        return converted_list_arguments, converted_keyword_arguments

    def convert_result(self, method, result):
        return self.from_target_to_source(result)

    def from_source_to_target(self, x):
        if isinstance(x, quantities.Quantity):
            return self.converter.from_source_to_target(x)
        else:
            return x

    def from_target_to_source(self, x):
        if isinstance(x, quantities.Quantity):
            if x.unit.is_non_numeric():
                return x
            else:
                return self.converter.from_target_to_source(x)
        elif isinstance(x, basestring):
            return x
        elif isinstance(x, numpy.ndarray):
            return x
        elif hasattr(x, '__len__'):
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
        self.is_determining_method = False

    def add_transition(self, from_state, to_state):
        self.transitions.append((from_state, to_state))

    def remove_transition(self, from_name, to_name):
        index = -1
        for i, transition in enumerate(self.transitions):
            from_state, to_state= transition
            if from_name == from_state.name and to_name == to_state.name:
                index = i
        if index >= 0:
            del self.transitions[index]
        
    def new_method(self, method = None):
        if method == None:
            if self.is_determining_method:
                raise Exception("A state is defined for a method with name '{0}', but the method is not implemented".format(self.function_name))
            self.is_determining_method = True
            try:
                method = getattr(self.interface, self.function_name)
            finally:
                self.is_determining_method = False
            
        return CodeMethodWrapper(method, self)

    def precall(self, method):
        stored_transitions = []
        for from_state, to_state in self.transitions:
            if from_state is None:
                return to_state
            elif from_state.matches(self.state_machine._current_state):
                return to_state
            else:
                stored_transitions.append((from_state, to_state))
        
        possible_paths = []
        for from_state, to_state  in stored_transitions:
            try:
                transition_path = self.state_machine._get_state_transition_path_to(from_state)
                possible_paths.append([transition_path, to_state])
            except Exception, ex:
                pass
        
        if len(possible_paths) == 0:            
            # do again to get an exception.
            message="While calling {0} of {1}: ".format(self.function_name, self.interface.__class__.__name__)
            try:
                self.state_machine._get_state_transition_path_to(stored_transitions[0][0])
            except exceptions.AmuseException as ex:
                args = list(ex.arguments)
                args[0] = message+str(args[0])
                ex.arguments = tuple(args)
                raise ex
            except Exception as ex:
                
                args = list(ex.args)
                args[0] = message+str(args[0])
                ex.args = tuple(args)
                raise ex
        
        for path, to_state in sorted(possible_paths, key = lambda x: len(x[0])):
            for transition in path:
                transition.do()
            return to_state


    def postcall(self, method, to_state):
        if to_state is None:
            return
        elif to_state.matches(self.state_machine._current_state):
            return
        else:
            self.state_machine._current_state = to_state
            
    def __str__(self):
        return "<StateMethod {0}>".format(self.function_name)


class HandleState(HandleCodeInterfaceAttributeAccess):

    def __init__(self, interface, **options):
        self._mapping_from_name_to_state_method = {}
        self.interface = interface
        self._state_machine = state.StateMachine(interface, **options)


    def supports(self, name, was_found):
        if name == 'state_machine':
            return True
        else:
            return self._state_machine.is_enabled and (name in self._mapping_from_name_to_state_method)

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

    def _remove_state_method(self, from_name, to_name, function_name):
        if function_name in self._mapping_from_name_to_state_method:
            state_method = self._mapping_from_name_to_state_method[function_name]
            state_method.remove_transition(from_name, to_name)

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
    
    
    def remove_transition(self, from_name, to_name, function_name):
    
        self._state_machine.remove_transition(from_name, to_name)
    
        self._remove_state_method(from_name, to_name, function_name)


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

    ERROR_CODE =  ErrorCodeMethodArgumentOrResultType()
    NO_UNIT = NoUnitMethodArgumentOrResultType()
    UNIT = UnitMethodArgumentOrResultType()
    INDEX = IndexMethodArgumentOrResultType()
    LINK = LinkMethodArgumentOrResultType

    def __init__(self, wrapped_object, function_name, units, return_units, name):
        self.function_name = function_name

        if hasattr(units, '__iter__'):
            self.units = units
        else:
            self.units = (units,)

        self.return_units = return_units
        self.is_return_units_iterable = hasattr(self.return_units, '__iter__')
        self.wrapped_object = wrapped_object
        self.name = name
        if self.return_units is None:
            self.handle_return_value = self.handle_as_errorcode
        else:
            self.handle_return_value = self.handle_as_unit
            
    def __getstate__(self):
        result = {}
        result.update(self.__dict__)
        del result['handle_return_value']
        return result
        
    def __setstate__(self, state):
        self.__dict__ = state
        if self.return_units is None:
            self.handle_return_value = self.handle_as_errorcode
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
            return MethodWithUnits(getattr(self.wrapped_object, self.function_name), self)

    def handle_errorcode(self, errorcode):
        if errorcode in self.wrapped_object.errorcodes:
            raise exceptions.AmuseException("Error when calling '{0}' of a '{1}', errorcode is {2}, error is '{3}'".format(self.name, type(self.wrapped_object).__name__, errorcode,  self.wrapped_object.errorcodes[errorcode]), errorcode)
        elif errorcode < 0:
            raise exceptions.AmuseException("Error when calling '{0}' of a '{1}', errorcode is {2}".format(self.name, type(self.wrapped_object).__name__, errorcode), errorcode)
        else:
            return errorcode

    def handle_as_errorcode(self, method, errorcode):
        if hasattr(errorcode, 'any'):
            if not errorcode.any():
                return
        if hasattr(errorcode, '__iter__'):
            for x in errorcode:
                self.handle_errorcode(x)
        else:
            self.handle_errorcode(errorcode)

    def handle_as_unit(self, method, return_value):
        if not self.is_return_units_iterable:
            return self.return_units.convert_result_value(method, self, return_value)
        else:
            if not hasattr(return_value, '__iter__'):
                return_value = [return_value]
                
            result = []
            for value, unit in itertools.izip(return_value, self.return_units):
                unit.append_result_value(method, self, value, result)
                    
            if len(result) == 1:                
                return result[0]
            else:
                return result


    def convert_arguments(self, method, list_arguments, keyword_arguments):
        result = {}
        input_parameters = method.method_input_argument_names

        for index, parameter in enumerate(input_parameters):
            if parameter in keyword_arguments:
                if self.units[index] == self.NO_UNIT:
                    arg = keyword_arguments[parameter]
                    if is_quantity(arg):
                        result[parameter] = arg.value_in(units.none)
                    else:
                        result[parameter] = arg
                elif self.units[index] == self.INDEX:
                    result[parameter] = keyword_arguments[parameter]
                elif self.units[index] == self.UNIT:
                    result[parameter] = keyword_arguments[parameter]
                else:
                    result[parameter] = quantities.value_in(keyword_arguments[parameter],self.units[index])

        for index, argument in enumerate(list_arguments):
            parameter = input_parameters[index]
            if parameter in result:
              raise ConvertArgumentsException(
              "got multiple values for argument '{0}' of method {1}".format(parameter, self.function_name))
            try:
                if self.units[index] == self.NO_UNIT:
                    if is_quantity(argument):
                        result[parameter] = argument.value_in(units.none)
                    else:
                        result[parameter] = argument
                elif self.units[index] == self.INDEX:
                    result[parameter] = argument
                elif self.units[index] == self.UNIT:
                    result[parameter] = argument
                elif type(self.units[index]) == self.LINK:
                    result[parameter] = self.units[index].convert_argument_value(method, self, argument)
                else:
                    if self.units[index].is_none() and not hasattr(argument,'unit'):
                        result[parameter] = argument
                    else:
                        result[parameter] = quantities.value_in(argument,self.units[index])
            except core.IncompatibleUnitsException as ex:
                raise ConvertArgumentsException("error while converting parameter '{0}', error: {1}".format(parameter, ex))
            except Exception as ex:
                raise exceptions.AmuseException("error while converting parameter '{0}', error: {1}".format(parameter, ex))
        
        return (), result

    def convert_result(self, method, result):
        return self.handle_return_value(method, result)
    
    # this function tries to determine the size of results from definition
    # its a bit ad-hoc. in spite of what the name suggests it determines it from scratch
    # (and not converting the index from the wrapped function's result_index)    
    def convert_result_index(self, method):
        if self.return_units is None:
            return None
        else:
            if not self.is_return_units_iterable:
                return None
            else:                
                nresult = 0
                for unit in self.return_units:
                    if not hasattr(unit,"_returns_result"):
                        nresult+=1
                    else:
                        if unit._returns_result:
                            nresult+=1
                
                return range(nresult)

    @late
    def has_same_name_as_original(self):
        return self.function_name == self.name


    @late
    def index_input_attributes(self):
        return map(lambda x : x == self.INDEX, self.units)

    @late
    def nbody_input_attributes(self):
        return map(lambda x : isinstance(x, UnitMethodArgumentOrResultType) or isinstance(x, unit) and generic_unit_system.is_generic_unit(x), self.units)

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
            raise IncorrectMethodDefinition(self.name, type(self.wrapped_object).__name__, number_expected_inputs, number_specified_inputs, 'inputs')
    
    

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
            raise IncorrectMethodDefinition(self.name, type(self.wrapped_object).__name__, number_expected_outputs, number_specified_outputs, 'outputs')
    
    
class HandleMethodsWithUnits(object):
    ERROR_CODE = MethodWithUnitsDefinition.ERROR_CODE
    NO_UNIT = MethodWithUnitsDefinition.NO_UNIT
    INDEX = MethodWithUnitsDefinition.INDEX
    LINK = MethodWithUnitsDefinition.LINK
    UNIT = MethodWithUnitsDefinition.UNIT

    def __init__(self, interface):
        self.method_definitions = {}
        self.interface = interface
        self.setup_units_from_legacy_interface()
        self.method_instances = {}
    
    def setup_units_from_legacy_interface(self):
        for name, specification in self.interface_function_specifications():
            units = [x.unit for x in specification.input_parameters]
            return_units = [x.unit for x in specification.output_parameters]
            
            if not specification.result_type is None:
                if specification.result_unit is None:
                    return_units.append(MethodWithUnitsDefinition.ERROR_CODE)
                else:
                    return_units.append(specification.result_unit)
                    
            default_to_nounit = lambda y : MethodWithUnitsDefinition.NO_UNIT if y is None else y
            
            return_units = [ default_to_nounit(x) for x in return_units]
            units = [default_to_nounit(x) for x in units]
            definition = MethodWithUnitsDefinition(
                self.interface,
                name,
                units,
                return_units,
                name
            )
            self.method_definitions[name] = definition

    def interface_function_specifications(self):
        interface_type = type(self.interface.legacy_interface)
        attribute_names = dir(interface_type)
        result = []
        for x in attribute_names:
            if x.startswith('__'):
                continue
            value = getattr(interface_type, x)
            if hasattr(value, 'specification') and hasattr(value.specification, 'input_parameters'):
                result.append( [x,value.specification])
        result.sort(key= lambda x: x[1].id)
        return result
        
    def supports(self, name, was_found):
        return name in self.method_definitions

    def get_attribute(self, name, value):
        if not name in self.method_instances:
            self.method_instances[name] = self.method_definitions[name].new_method(value)

        return self.method_instances[name]


    def attribute_names(self):
        return set(self.method_definitions.keys())

    def add_method(self, original_name, units, return_unit = None,  public_name = None):
        if public_name is None:
            public_name = original_name

        definition = MethodWithUnitsDefinition(
            self.interface,
            original_name,
            units,
            return_unit,
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
        self.definitions = defaultdict(list,parameters=[])
        self.parameters = {}

    def supports(self, name, was_found):
        return name in self.definitions.keys()

    def get_attribute(self, name, value):
        # note: parameters can be added after init, not yet removed

        name=name or 'parameters'

        if name not in self.parameters:
            d=self.definitions[name]
            self.parameters[name] = parameters.new_parameters_instance_with_docs(d, self.interface)
        else:
            self.parameters[name].update()       
        result=self.parameters[name]
        return result

    def attribute_names(self):
        return set(self.definitions.keys())

    def add_method_parameter(self, get_method, set_method, name, description, 
            default_value = None,must_set_before_get = False, is_vector = False,
            parameter_set='parameters'):
        if is_vector:
            definition = parameters.ModuleVectorMethodParameterDefinition(
                get_method,
                set_method,
                name,
                description,
                default_value,
                must_set_before_get = must_set_before_get
            )
        else:
            definition = parameters.ModuleMethodParameterDefinition(
                get_method,
                set_method,
                name,
                description,
                default_value,
                must_set_before_get = must_set_before_get
            )
        self.definitions[parameter_set].append(definition)
        
    def add_alias_parameter(self, name, aliased_name, description,parameter_set='parameters',alias_set=None):
        definition = parameters.AliasParameterDefinition(
            name,
            aliased_name,
            description,
            alias_set=alias_set
        )
        self.definitions[parameter_set].append(definition)


    def add_caching_parameter(self, function_name, parameter_name, name, description, default_value = None,parameter_set='parameters'):
        definition = parameters.ModuleCachingParameterDefinition(
            function_name,
            parameter_name,
            name,
            description,
            default_value
        )
        self.definitions[parameter_set].append(definition)

    def add_boolean_parameter(self, get_method, set_method, name, description, default_value = None,parameter_set='parameters'):
        definition = parameters.ModuleBooleanParameterDefinition(
            get_method,
            set_method,
            name,
            description,
            default_value
        )
        self.definitions[parameter_set].append(definition)

    def add_default_form_parameter(self,name,description,default,parameter_set='parameters'):
        if isinstance(default,bool):
          self.add_boolean_parameter("get_"+name,"set_"+name,name,description,default,parameter_set='parameters')
        else:
          self.add_method_parameter("get_"+name,"set_"+name,name,description,default,parameter_set='parameters')

    def add_array_parameter(self, get_method, set_method, range_method, name, description,parameter_set='parameters'):
        definition = parameters.ModuleArrayParameterDefinition(
            get_method,
            set_method,
            range_method,
            name,
            description
        )
        self.definitions[parameter_set].append(definition)

    def has_name(self, name):
        return name == 'PARAMETER'

    def setup(self, object):
        object.define_parameters(self)



    def add_vector_parameter(self, name, description, parameter_names,parameter_set='parameters'):
        default_value = None
        for parameter_name in parameter_names:
            for defined_parameter in self.definitions[parameter_set]:
                if defined_parameter.name == parameter_name:
                    if default_value is None:
                        if not is_quantity(defined_parameter.default_value):
                            default_value  = []
                        else:
                            default_value = quantities.AdaptingVectorQuantity()
                    default_value.append(defined_parameter.default_value)
        definition = parameters.VectorParameterDefinition(
            name,
            description,
            parameter_names,
            default_value
        )
        self.definitions[parameter_set].append(definition)

    def add_interface_parameter(self, name, description, default_value,state_guard=None,parameter_set='parameters'):        
        definition=parameters.InterfaceParameterDefinition(name,description,default_value,state_guard=state_guard)
        self.definitions[parameter_set].append(definition)
    
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
    
    def add_gridded_getter(self, name_of_the_getter, name_of_the_range_method, names = None):
        self.gridded_getters.append((name_of_the_getter,name_of_the_range_method, names))
        
    def add_gridded_setter(self, name_of_the_setter, name_of_the_range_method, names = None):
        self.gridded_setters.append((name_of_the_setter,name_of_the_range_method, names))
        
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
    

    
    def add_subselect_from_particle(self, name, get_number_of_particles_name = None,  public_name = None):
        if not public_name:
            public_name = name
        self.subselects_from_particle.append((name, get_number_of_particles_name, public_name))

class ParticleSetDefinition(AbstractParticleSetDefinition):

    def __init__(self, handler):
        self.handler = handler
        self.name_of_indexing_attribute = 'index_of_the_particle'
        self.new_particle_method = ('new_particle',(), None)
        self.name_of_delete_particle_method = 'delete_particle'
        self.name_of_number_of_particles_method = 'get_number_of_particles'
        self.setters = []
        self.getters = []
        self.gridded_getters = []
        self.gridded_setters = []
        self.queries = []
        self.attributes = []
    
        self.selects_form_particle = []
        self.subselects_in_set = []
        self.subselects_from_particle = []
        self.methods = []
        self.is_superset = False
        self.is_inmemory = False
        self.particles_factory = datamodel.Particles

    def new_storage(self, interface):
    
        if self.is_inmemory:
            return datamodel.get_in_memory_attribute_storage_factory()()
    
        setters = []
        for name, names in self.setters:
            x = incode_storage.ParticleSetAttributesMethod(getattr(interface, name), names)
            setters.append(x)
    
        getters = []
        for name, names in self.getters:
            x = incode_storage.ParticleGetAttributesMethod(getattr(interface, name), names)
            getters.append(x)
        
        for name, range_method_name, names in self.gridded_getters:
            x = incode_storage.ParticleGetGriddedAttributesMethod(
                    getattr(interface, name), 
                    getattr(interface, range_method_name), 
                    names
            )
            getters.append(x)
            
        for name, range_method_name, names in self.gridded_setters:
            x = incode_storage.ParticleSetGriddedAttributesMethod(
                    getattr(interface, name), 
                    getattr(interface, range_method_name), 
                    names
            )
            setters.append(x)
            
    
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
        
        
        selects = self.new_subselects_from_particle(handler.interface)
        for x in selects:
        #result.add_function_attribute(x.public_name, x.apply_on_all)
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
    
    
    def new_subselects_from_particle(self, interface):
        results = []
        for name, get_number_of_particles_name, public_name in self.subselects_from_particle:
        
            number_of_particles_method = None if get_number_of_particles_name is None else getattr(interface, get_number_of_particles_name)
            x = incode_storage.ParticleSpecificSelectSubsetMethod(getattr(interface, name), number_of_particles_method, public_name)
            results.append(x)

        return results

class ParticleSupersetDefinition(AbstractParticleSetDefinition):

    def __init__(self, handler, particle_subset_names, index_to_default_set=None):
        self.handler = handler
        self.particle_subset_names = particle_subset_names
        self.index_to_default_set = index_to_default_set
        self.is_superset = True
        self.particles_factory = datamodel.ParticlesSuperset
        self.queries = []
        
    
    def new_set_instance(self, handler):
        subsets = [handler.get_attribute(subset_name, None) for subset_name in self.particle_subset_names]
        result = self.particles_factory(
            subsets, 
            index_to_default_set=self.index_to_default_set
        )
        for one_query in self.new_queries(handler.interface):
            result.add_function_attribute(one_query.public_name, one_query.apply)
        return result
    
    def new_queries(self, interface):
        queries = []
        for name, names, public_name in self.queries:
            x = incode_storage.ParticleQueryMethod(getattr(interface, name), names, public_name, query_superset=True)
            queries.append(x)
        return queries
    

class GridDefinition(AbstractParticleSetDefinition):

    def __init__(self, handler, grid_class=datamodel.Grid):
        self.handler = handler
        self.axes_names = None
        self.name_of_the_get_range_method = 'get_range'
        self.setters = []
        self.getters = []
        self.gridded_setters=[]
        self.gridded_getters=[]
        self.particles_factory = grid_class
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

        for name, range_method_name, names in self.gridded_getters:
            x = incode_storage.ParticleGetGriddedAttributesMethod(
                    getattr(interface, name), 
                    getattr(interface, range_method_name), 
                    names
            )
            getters.append(x)
            
        for name, range_method_name, names in self.gridded_setters:
            x = incode_storage.ParticleSetGriddedAttributesMethod(
                    getattr(interface, name), 
                    getattr(interface, range_method_name), 
                    names
            )
            setters.append(x)

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
        if self.axes_names is not None:
            result.add_vector_attribute("position",self.axes_names)
        return result

class CodeInMemoryParticles(datamodel.Particles):

    def __init__(self, code_interface = None, storage = None):
        datamodel.Particles.__init__(self, storage = storage)
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
            if set_definition.state_guard:
                getattr(self.interface, set_definition.state_guard)()
            result = set_definition.new_set_instance(self)
                
            self.mapping_from_name_to_set_instance[name] = result
            return result

    def attribute_names(self):
        return set(self.mapping_from_name_to_set_definition.keys())


    def has_name(self, name):
        return name == 'PARTICLES' or name == 'DATASETS' or name == 'GRIDS'

    def setup(self, object):
        object.define_particle_sets(self)
        object.define_data_sets(self)
        object.define_grids(self)

    def define_set(self, name, name_of_indexing_attribute = 'index_of_the_particle', state_guard=None):
        definition = ParticleSetDefinition(self)
        definition.name_of_indexing_attribute = name_of_indexing_attribute
        definition.state_guard = state_guard
        self.mapping_from_name_to_set_definition[name] = definition

    def define_super_set(self, name, particle_subsets, index_to_default_set = None, state_guard=None):
        definition = ParticleSupersetDefinition(self, particle_subsets, index_to_default_set)
        definition.state_guard = state_guard
        self.mapping_from_name_to_set_definition[name] = definition

    def define_inmemory_set(self, name, particles_factory = CodeInMemoryParticles, state_guard=None):
        definition = ParticleSetDefinition(self)
        definition.is_inmemory = True
        definition.particles_factory = particles_factory
        definition.state_guard = state_guard
        self.mapping_from_name_to_set_definition[name] = definition
    
    def define_grid(self, name, name_of_indexing_attribute = 'index_of_the_particle', axes_names = None, grid_class=datamodel.Grid, state_guard=None):
        definition = GridDefinition(self, grid_class = grid_class)
        definition.name_of_indexing_attribute = name_of_indexing_attribute
        definition.axes_names = axes_names
        definition.state_guard = state_guard
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

    def add_gridded_getter(self, name_of_the_set, name_of_the_getter, name_of_the_range_method, names = None):
        self.mapping_from_name_to_set_definition[name_of_the_set].add_gridded_getter(name_of_the_getter, name_of_the_range_method, names = names)
    
    def add_gridded_setter(self, name_of_the_set, name_of_the_setter, name_of_the_range_method, names = None):
        self.mapping_from_name_to_set_definition[name_of_the_set].add_gridded_setter(name_of_the_setter, name_of_the_range_method, names = names)
    
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
        
    def add_subselect_from_particle(self, name_of_the_set, name, get_number_of_particles_name = None,  public_name = None):
        self.mapping_from_name_to_set_definition[name_of_the_set].add_subselect_from_particle(
            name, 
            get_number_of_particles_name = get_number_of_particles_name,  
            public_name = public_name
        )

    def _cleanup_instances(self):
        self.mapping_from_name_to_set_instance = {}

class OverriddenCodeInterface(object):

    def __init__(self, code_interface):
        self.code_interface = code_interface

    def __getattr__(self, name):
        return self.code_interface.__getattr__(name)

class InCodeComponentImplementation(OldObjectsBindingMixin, OptionalAttributes):

    def __init__(self, legacy_interface, **options):
        OptionalAttributes.__init__(self, **options)
        self.legacy_interface = legacy_interface
        self._options = options
        self._handlers = []
        self.__init_handlers__(legacy_interface, options)

    def __init_handlers__(self, legacy_interface, options):
        self._handlers.append(LegacyInterfaceHandler(legacy_interface))
        self._handlers.append(HandleMethodsWithUnits(self))
        self._handlers.append(HandlePropertiesWithUnits(self))
        self._handlers.append(HandleParameters(self))
        self._handlers.append(HandleParticles(self))
        if self.must_handle_state:
            self._handlers.append(HandleState(self, **options))
        self._handlers.append(HandleConvertUnits(self))
        self._handlers.append(HandleErrorCodes(self))

        self.setup()

    @option(type='boolean', sections=("code", "state",))
    def must_handle_state(self):
        return True
    
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

    def define_data_sets(self, handler):
        pass

    def define_grids(self, handler):
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
        definition = GridDefinition(handler, grid_class = extra_arguments.get("grid_class", datamodel.Grid))
        builder_function(definition, **extra_arguments)
        return definition.new_set_instance(handler)
        
    def __setstate__(self, state):
        self.__dict__ = state

    def data_store_names(self):
        self.before_get_data_store_names()
        return self.get_handler('PARTICLES').mapping_from_name_to_set_definition.keys()

    def parameter_set_names(self):
        #~ self.before_get_data_store_names()
        return self.get_handler('PARAMETER').definitions.keys()

    
class IncorrectMethodDefinition(IncorrectWrappedMethodException):
    formatstring = "Incorrect definition of method '{0}' of class '{1}', the number of {4} do not match, expected {2}, actual {3}."


class PropertyDefinition(object):

    def __init__(self, handler, functionname, publicname, keyword_arguments = {}):
        self.functionname = functionname
        self.publicname = publicname
        self.handler = handler
        self.keyword_arguments = {}
        self._method = None
        
    def get_value(self, original):
        if self._method is None:
            self._method = getattr(self.handler.interface, self.functionname)
         
        result = self._method (**self.keyword_arguments)
        if hasattr(result, "__iter__"):
            return quantities.VectorQuantity.new_from_scalar_quantities(*result)
        else:
            return result


