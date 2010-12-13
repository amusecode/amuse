import weakref

from amuse.support.units import nbody_system
from amuse.support.units import generic_unit_system
from amuse.support.data import values
from amuse.support import exceptions

import warnings

class ParameterDocumentation(object):

    def __get__(self, instance, owner):
        output = "Parameters: \n"

        for parameter_definition in instance._definitions:
            output += parameter_definition.name + "\n\n"
            output += "    " + parameter_definition.description
            output += " (default value:" + str(instance.get_default_value_for(parameter_definition.name)) + ")\n\n"

        return output

class Parameters(object):
    __doc__ = ParameterDocumentation()
    __name__ = 'Parameters'
    
    def __init__(self, definitions, instance):
        object.__setattr__(self, '_instance', weakref.ref(instance))
        object.__setattr__(self, '_definitions', definitions)
        object.__setattr__(self, '_mapping_from_name_to_definition', {})
        object.__setattr__(self, '_mapping_from_name_to_parameter', {})
    
        for x in definitions:
            self._mapping_from_name_to_definition[x.name] = x

    def __getattr__(self, name):
        #if name.startswith('__'):
        #    return object.__getattribute__(self, name)
        if not name in self._mapping_from_name_to_definition:
            raise exceptions.CoreException("tried to get unknown parameter '{0}' for a '{1}' object".format(name, type(self._instance()).__name__))
            
        return self.get_parameter(name).get_value()

    def __setattr__(self, name, value):
        if not name in self._mapping_from_name_to_definition:
            warnings.warn("tried to set unknown parameter '{0}' for a '{1}' object".format(name, type(self._instance()).__name__), exceptions.AmuseWarning)
            return
        
        if hasattr(self._instance(), "get_name_of_current_state"):
            current_state = self._instance().get_name_of_current_state()
        else:
            current_state = None
        if current_state == "UNINITIALIZED" and hasattr(self._instance(), "invoke_state_change"):
            self._instance().invoke_state_change()
        elif ((current_state == "EDIT" or current_state == "RUN") and 
                hasattr(self._instance(), "invoke_state_change2")):
            self._instance().invoke_state_change2()
            
        return self.get_parameter(name).set_value(value)

    def names(self):
        return self._mapping_from_name_to_definition.keys()

    def set_defaults(self):
        for name in self.names():
            parameter = self.get_parameter(name)
            parameter.set_default_value()

    def __dir__(self):
        result = []
        result.extend(dir(type(self)))
        result.extend(self.names())
        return result


    def get_default_value_for(self, name):
        if not name in self._mapping_from_name_to_definition:
            raise exceptions.CoreException("tried to get default value of unknown parameter '{0}' for a '{1}' object".format(name, type(self._instance()).__name__))
    
        definition = self._mapping_from_name_to_definition[name]
        return definition.default_value

    def __str__(self):
        output = ""

        for name in self.names():
            output += name + ": "
            output += str(getattr(self, name))+"\n"

        return output



    def get_parameter(self, name):
        if not name in self._mapping_from_name_to_definition:
            raise exceptions.AmuseException("{0!r} not defined as parameter".format(name))
        
        if not name in self._mapping_from_name_to_parameter:
            definition = self._mapping_from_name_to_definition[name]
            self._mapping_from_name_to_parameter[name] = Parameter(definition, self)
            
        return self._mapping_from_name_to_parameter[name]
    
    

    def iter_parameters(self):
        for name in self.names():
            yield self.get_parameter(name)
    
    

    def send_cached_parameters_to_code(self):
        cached_parameters = [x for x in self.iter_parameters() if x.definition.is_cached()]
        for x in cached_parameters:
            if not x.is_set:
                x.set_default_value()
        
        functions = {}
        for x in cached_parameters:
            definition = x.definition
            if not definition.functionname in functions:
                functions[definition.functionname] = []
            functions[definition.functionname].append(x)
            
        for functionname, parameters in functions.iteritems():
            object = self._instance()
            method = getattr(object, functionname)
            keyword_arguments = {}
            for parameter in parameters:
                keyword_arguments[parameter.definition.parameter_name] = parameter.get_cached_value()
            errorcode = method(**keyword_arguments)
    
    

    def send_not_set_parameters_to_code(self):
        parameters = [x for x in self.iter_parameters() if x.must_set_to_default()]
        for x in parameters:
            x.set_default_value()
    
    
class ParametersWithUnitsConverted(object):

    __doc__ = ParameterDocumentation()

    def __init__(self, original, converter):
        object.__setattr__(self, '_original', original)
        object.__setattr__(self, '_converter', converter)


    def __getattr__(self, name):
        original_value = getattr(self._original, name)
        if isinstance(original_value,bool) or isinstance(original_value,values.NonNumericQuantity):
            return original_value
        else:
            return self._converter.from_target_to_source(original_value)

    def __setattr__(self, name, value):
        if not name in self._original._mapping_from_name_to_definition:
            warnings.warn("tried to set unknown parameter '{0}' for a '{1}' object".format(name, type(self._instance()).__name__), exceptions.AmuseWarning)
            return
        default_value = self._original.get_default_value_for(name)
        if not isinstance(default_value,bool) and generic_unit_system.is_generic_unit(default_value.unit):
            setattr(self._original, name, self._converter.from_source_to_target(value))
        else:
            setattr(self._original, name, value)

    def names(self):
        return self._original.names()

    def set_defaults(self):
        self._original.set_defaults()

    def __dir__(self):
        return dir(self._original)

    def get_default_value_for(self, name):
        return self._converter.from_target_to_source(self._original.get_default_value_for(name))

    def __str__(self):
        output = ""

        for name in self.names():
            output += name + ": "
            output += str(getattr(self, name))+"\n"

        return output

class AbstractParameterDefinition(object):
    def __init__(self, name, description):
        self.name = name
        self.description = description

    def get_value(self, parameter, object):
        raise AmuseException("not implemented")

    def set_value(self, parameter, object, quantity):
        raise AmuseException("not implemented")

    def set_default_value(self, parameter, object):
        pass

    def is_readonly(self):
        return False
        
    def is_cached(self):
        return False



    def must_set_to_default_if_not_set(self):
        return True
    
    
class ParameterDefinition(AbstractParameterDefinition):
    def __init__(self, name, description, unit, default_value, must_set_before_get = False):
        AbstractParameterDefinition.__init__(self, name, description)
        self.unit = unit
        self.default_value = default_value
        self.must_set_before_get = must_set_before_get


    def get_value(self, parameter, object):
        if self.must_set_before_get and not parameter.is_set:
            self.set_default_value(parameter, object)
            
        result = self.unit.new_quantity(self.get_legacy_value(parameter, object))
        return result

    def set_value(self, parameter, object, quantity):
        if self.unit.is_non_numeric() or len(self.unit.base) == 0:
            if not isinstance(quantity, values.Quantity):
                quantity = quantity | self.unit
        
        self.set_legacy_value(parameter, object, quantity.value_in(self.unit))
        
        parameter.is_set = True

    def set_default_value(self, parameter, object):
        if self.default_value is None :
            return None
    
        if self.is_readonly():
            return None
    
        self.set_value(parameter, object, self.default_value)

    def is_readonly(self):
        return False


    def is_cached(self):
        return False
    
    
class ModuleAttributeParameterDefinition(ParameterDefinition):
    def __init__(self, attribute_name, name, description, unit, default_value = None, must_set_before_get = False):
        ParameterDefinition.__init__(self, name, description, unit, default_value, must_set_before_get)
        self.attribute_name = attribute_name

    def get_legacy_value(self, parameter, object):
        return getattr(object.legacy_interface, self.attribute_name)

    def set_legacy_value(self,  parameter, object, number):
        setattr(object.legacy_interface, self.attribute_name, number)

class ParameterException(AttributeError):
    template = ("Could not {0} value for parameter '{1}' of a '{2}' object, got errorcode <{3}>")

    def __init__(self, object, parameter_name, errorcode, is_get):
        AttributeError.__init__(self, self.template.format(
            "get" if is_get else "set",
            parameter_name,
            type(object).__name__,
            errorcode
        ))
        self.errorcode = errorcode
        self.parameter_name = parameter_name


class ModuleMethodParameterDefinition_Next(ParameterDefinition):
    def __init__(self, get_method, set_method, name, description, unit, default_value = None, must_set_before_get = False):
        ParameterDefinition.__init__(self, name, description, unit, default_value, must_set_before_get)
        self.get_method = get_method
        self.set_method = set_method
        self.stored_value = None


    def get_legacy_value(self, parameter, object):
        if self.get_method is None:
            return self.stored_value
        else:
            (result, error) = getattr(object, self.get_method)()
            if error < 0:
                raise ParameterException(object, self.name, error, True)
            else:
                return result

    def set_legacy_value(self, parameter, object, number):
        if self.set_method is None:
            raise exceptions.CoreException("Could not set value for parameter '{0}' of a '{1}' object, parameter is read-only".format(self.name, type(object).__name__))
    
        error = getattr(object, self.set_method)(number)
        if error < 0:
            raise ParameterException(object, self.name, error, False)
        else:
            if self.get_method is None:
                self.stored_value = number



    def is_readonly(self):
        return self.set_method is None


class ModuleBooleanParameterDefinition(ParameterDefinition):
    def __init__(self, get_method, set_method, name, description, default_value = None, must_set_before_get = False):
        ParameterDefinition.__init__(self, name, description, None, default_value, must_set_before_get)
        self.get_method = get_method
        self.set_method = set_method
        self.stored_value = None

    def get_value(self, parameter, object):
        return True if self.get_legacy_value(parameter, object) else False

    def set_value(self, parameter, object, bool):
        self.set_legacy_value(parameter, object, 1 if bool else 0)

    def get_legacy_value(self,  parameter, object):
        if self.get_method is None:
            return self.stored_value
        else:
            (result, error) = getattr(object, self.get_method)()
            if error < 0:
                raise ParameterException(object, self.name, error, True)
            else:
                return result

    def set_legacy_value(self,  parameter, object, number):
        if self.set_method is None:
            raise exceptions.CoreException("Could not set value for parameter '{0}' of a '{1}' object, parameter is read-only".format(self.name, type(object).__name__))
    
        error = getattr(object, self.set_method)(number)
        if error < 0:
            raise ParameterException(object, self.name, error, False)
        else:
            if self.get_method is None:
                self.stored_value = number

    def is_readonly(self):
        return self.set_method is None


class ModuleCachingParameterDefinition(ParameterDefinition):
    def __init__(self, functionname, parameter_name, name, description, unit, default_value = None):
        ParameterDefinition.__init__(self, name, description, unit, default_value, must_set_before_get = True)
        self.stored_value = None
        self.parameter_name = parameter_name
        self.functionname = functionname

    def get_legacy_value(self, parameter, object):
        return parameter.cached_value

    def set_legacy_value(self, parameter, object, number):
        parameter.cached_value = number


    def is_cached(self):
        return True
    
    
class Parameter(object):

    def __init__(self, definition, parameter_set):
        self.parameter_set = parameter_set
        self.definition = definition
        self.is_set = False
        
    
    def get_value(self):
        return self.definition.get_value(self, self.parameter_set._instance())
    
    def set_value(self, quantity):
        return self.definition.set_value(self, self.parameter_set._instance(), quantity)

    def set_default_value(self):
        self.definition.set_default_value(self, self.parameter_set._instance())
         
    def is_readonly(self):
        return self.definition.is_readonly()
    
    def must_set_to_default(self):
        if self.definition.is_cached():
            return False
            
        return (not self.is_set) and self.definition.must_set_to_default_if_not_set()

    def get_cached_value(self):
        return self.definition.get_legacy_value(self, self.parameter_set._instance())
    
    

class VectorParameterDefinition(AbstractParameterDefinition):
    def __init__(self, name, description, names_of_parameters):
        AbstractParameterDefinition.__init__(self, name, description)
        self.names_of_parameters = names_of_parameters

    def get_value(self, parameter, object):
        all_parameters = parameter.parameter_set
        result = []
        unit = None
        for name in self.names_of_parameters:
            parameter = all_parameters.get_parameter(name)
            element = parameter.get_value()
            if unit is None and hasattr(element, 'unit'):
                unit = element.unit
                
            result.append(element.value_in(unit))
        return unit.new_quantity(result)
        
    def set_value(self, parameter, object, quantity):
        all_parameters = parameter.parameter_set
        for index,name in enumerate(self.names_of_parameters):
            parameter = all_parameters.get_parameter(name)
            parameter.set_value(quantity[index])


    def must_set_to_default_if_not_set(self):
        return False
    
    

    def get_unit(self, parameter):
        result = None
        all_parameters = parameter.parameter_set
        for index,name in enumerate(self.names_of_parameters):
            parameter = all_parameters.get_parameter(name)
            if hasattr(parameter.definition, "unit"):
                result = parameter.definition.unit
        return result
    
    
