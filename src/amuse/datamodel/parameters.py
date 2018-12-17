import weakref
import numpy
from amuse.units import nbody_system
from amuse.units import generic_unit_system
from amuse.units import quantities
from amuse.units.core import IncompatibleUnitsException
from amuse.units.quantities import is_quantity
from amuse.support import exceptions

import warnings

from amuse.support.core import OrderedDictionary

class Parameters(object):
    __name__ = 'Parameters'
    
    def __init__(self, definitions, instance):
        object.__setattr__(self, '_instance', weakref.ref(instance))
        object.__setattr__(self, '_definitions', definitions)
        object.__setattr__(self, '_mapping_from_name_to_definition', OrderedDictionary())
        object.__setattr__(self, '_mapping_from_name_to_parameter', OrderedDictionary())

        for x in definitions:
            self._mapping_from_name_to_definition[x.name] = x
        
        

    def __getattr__(self, name):
        #if name.startswith('__'):
        #    return object.__getattribute__(self, name)
        if not name in self._mapping_from_name_to_definition:
            raise exceptions.CoreException("tried to get unknown parameter '{0}' for a '{1}' object".format(name, type(self._instance()).__name__))
        
        self._instance().before_get_parameter()
        
        return self.get_parameter(name).get_value()

    def __setattr__(self, name, value):
        if not name in self._mapping_from_name_to_definition:
            warnings.warn("tried to set unknown parameter '{0}' for a '{1}' object".format(name, type(self._instance()).__name__), exceptions.AmuseWarning)
            return
        
        self._instance().before_set_parameter()

        return self.get_parameter(name).set_value(value)

    def names(self):
        return self._mapping_from_name_to_definition.keys()

    def set_defaults(self):
        
        self._instance().before_set_parameter()
        
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
        return definition.get_default_value(self)


    def __str__(self):
        output = ""

        for name in sorted(self.names()):
            output += name + ": "
            output += str(getattr(self, name))
            if self.get_parameter(name).is_readonly():
                output += "  (read only)"
            output += "\n"

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
        
        functions = OrderedDictionary()
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
            
    def check_defaults(self):
        for x in self.iter_parameters():
            default_value = self.get_default_value_for(x.definition.name)
            try:
                value = x.get_value()
            except:
                print "could not get value for:", x.definition.name, default_value
                continue
            print x.definition.name, value, default_value
            if not value == default_value:
               print "!default value is not equal to value in code: {0}".format(x.definition.name)
            

    def copy(self):
        mapping_from_name_to_value = {}
        for name in self.names():
            mapping_from_name_to_value[name] = getattr(self, name)
            
        return ParametersMemento(mapping_from_name_to_value)
    
    def reset_from_memento(self, memento):
        for name in memento.names():
            if not name in self._mapping_from_name_to_definition:
                warnings.warn("tried to set unknown parameter '{0}' for a '{1}' object".format(name, type(self._instance()).__name__), exceptions.AmuseWarning)
                return
            
            if self.get_parameter(name).is_readonly():
                if not getattr(memento, name) == getattr(self, name):
                    warnings.warn("tried to change read-only parameter '{0}' for a '{1}' object".format(name, type(self._instance()).__name__), exceptions.AmuseWarning)
            else:
                setattr(self, name, getattr(memento, name))
            

    def has_writable_parameter(self, name):
        if not name in self._mapping_from_name_to_definition:
            return False
        return not self.get_parameter(name).is_readonly()




class ParametersMemento(object):
    __name__ = 'Parameters'
    
    def __init__(self, mapping_from_name_to_value = None):
        if mapping_from_name_to_value is None:
            mapping_from_name_to_value = {}
            
        object.__setattr__(self, '_mapping_from_name_to_value', mapping_from_name_to_value)
        
    def __getstate__(self):
        return self.__dict__

    def __setstate__(self,state):
        object.__setattr__(self, '__dict__', state)

    def __getattr__(self, name):
        if not name in self._mapping_from_name_to_value:
            raise exceptions.CoreException("tried to get unknown parameter '{0}'".format(name))
            
        
        return self._mapping_from_name_to_value[name]

    def __setattr__(self, name, value):
        if not name in self._mapping_from_name_to_value:
            warnings.warn("tried to set unknown parameter '{0}'".format(name), exceptions.AmuseWarning)
            return
            
        self._mapping_from_name_to_value[name] = value

    def names(self):
        return self._mapping_from_name_to_value.keys()

    def set_defaults(self):
        pass
        
    def __dir__(self):
        result = []
        result.extend(dir(type(self)))
        result.extend(self.names())
        return result


    def get_default_value_for(self, name):
        if not name in self._mapping_from_name_to_value:
            raise exceptions.CoreException("tried to get default value of unknown parameter '{0}'".format(name))
    
        raise exceptions.CoreException("tried to get default value, for a parameter in a parameters memento")


    def __str__(self):
        output = ""

        for name in sorted(self.names()):
            output += name + ": "
            output += str(getattr(self, name))+"\n"

        return output

    
    

    

def new_parameters_instance_with_docs(definitions, instance):
    
    class _ParametersMetaclass(type):
        def _get_doc(self):
            output = "Parameters: \n"
            for parameter_definition in definitions:
                output += parameter_definition.name + "\n\n"
                output += "    " + parameter_definition.description
                output += " (default value:" + str(parameter_definition.default_value) + ")\n\n"
            return output
        __doc__ = property(_get_doc)
    
    class ParametersWithDocs(Parameters):
        __metaclass__ = _ParametersMetaclass
        def _get_doc(self):
            output = "Parameters: \n"
            for parameter_definition in definitions:
                output += parameter_definition.name + "\n\n"
                output += "    " + parameter_definition.description
                output += " (default value:" + str(self.get_default_value_for(parameter_definition.name)) + ")\n\n"
            return output
        __doc__ = property(_get_doc)
    
    return ParametersWithDocs(definitions, instance)


def new_parameters_with_units_converted_instance_with_docs(original, converter):
    
    class _ParametersMetaclass(type):
        def _convert_from_target_to_source_if_needed(value):
            if isinstance(value, bool) or isinstance(value, quantities.NonNumericQuantity):
                return value
            else:
                return converter.from_target_to_source(value)
        def _get_doc(self):
            output = "Parameters: \n"
            for parameter_definition in original._definitions:
                value = parameter_definition.default_value
                if not isinstance(value, bool) and not isinstance(value, quantities.NonNumericQuantity):
                    value = converter.from_target_to_source(value)
                output += parameter_definition.name + "\n\n"
                output += "    " + parameter_definition.description
                output += " (default value:" + str(value) + ")\n\n"
            return output
        __doc__ = property(_get_doc)
    
    class ParametersWithDocs(ParametersWithUnitsConverted):
        __metaclass__ = _ParametersMetaclass
        def _get_doc(self):
            output = "Parameters: \n"
            for parameter_definition in original._definitions:
                output += parameter_definition.name + "\n\n"
                output += "    " + parameter_definition.description
                output += " (default value:" + str(self.get_default_value_for(parameter_definition.name)) + ")\n\n"
            return output
        __doc__ = property(_get_doc)
    
    return ParametersWithDocs(original, converter)


class ParametersWithUnitsConverted(object):

    def __init__(self, original, converter):
        object.__setattr__(self, '_original', original)
        object.__setattr__(self, '_converter', converter)


    def __getattr__(self, name):
        return self.convert_from_target_to_source_if_needed(getattr(self._original, name))

    def __setattr__(self, name, value):
        if not name in self._original._mapping_from_name_to_definition:
            warnings.warn("tried to set unknown parameter '{0}' for a '{1}' object".format(name, type(self._instance()).__name__), exceptions.AmuseWarning)
            return
        try:
            setattr(self._original, name, self._converter.from_source_to_target(value))
        except IncompatibleUnitsException as ex:
            setattr(self._original, name, value)

    def names(self):
        return self._original.names()

    def set_defaults(self):
        self._original.set_defaults()

    def __dir__(self):
        return dir(self._original)

    def convert_from_target_to_source_if_needed(self, value):
        if isinstance(value, bool) or isinstance(value, quantities.NonNumericQuantity):
            return value
        else:
            return self._converter.from_target_to_source(value)
        
    def get_default_value_for(self, name):
        return self.convert_from_target_to_source_if_needed(self._original.get_default_value_for(name))

    def __str__(self):
        output = ""

        for name in sorted(self.names()):
            output += name + ": "
            output += str(getattr(self, name))
            output += " default: " + str(self.get_default_value_for(name))
            output +="\n"

        return output
        
    
    def check_defaults(self):
        for x in self.iter_parameters():
            default_value = self.get_default_value_for(x.definition.name)
            try:
                value = x.get_value()
            except:
                print "could not get value for:", x.definition.name, default_value
                continue
            print x.definition.name, value, default_value
            if not value == default_value:
                print "default value is not equal to value in code: {0}".format(x.definition.name)
            

class AbstractParameterDefinition(object):
    def __init__(self, name, description):
        self.name = name
        self.description = description
    
    def get_default_value(self, parameterset):
        return self.default_value
        
    def get_value(self, parameter, object):
        raise exceptions.AmuseException("not implemented")

    def set_value(self, parameter, object, quantity):
        raise exceptions.AmuseException("not implemented")

    def set_default_value(self, parameter, object):
        pass

    def is_readonly(self):
        return False
        
    def is_cached(self):
        return False

    def must_set_to_default_if_not_set(self):
        return True

class AliasParameterDefinition(AbstractParameterDefinition):
    
    def __init__(self, name, aliased_name, description, alias_set=None):
        AbstractParameterDefinition.__init__(self, name, description)
        self.aliased_name = aliased_name
        self.alias_set = alias_set
        self.default_value = None
    
    def get_default_value(self, parameter_set):
        return parameter_set.get_parameter(self.aliased_name).definition.get_default_value(parameter_set)
        
    def get_value(self, parameter, object):
        if self.alias_set:
            parameter_set=getattr(object, self.alias_set)
        else:
            parameter_set=parameter.parameter_set
        return getattr(parameter_set, self.aliased_name)

    def set_value(self, parameter, object, quantity):
        if self.alias_set:
            parameter_set=getattr(object, self.alias_set)
        else:
            parameter_set=parameter.parameter_set
        return setattr(parameter_set, self.aliased_name, quantity)

    def set_default_value(self, parameter, object):
        pass

    def is_readonly(self):
        return False
        
    def is_cached(self):
        return False

    def must_set_to_default_if_not_set(self):
        return False
    
    
class ParameterDefinition(AbstractParameterDefinition):
    def __init__(self, name, description, default_value, must_set_before_get = False):
        AbstractParameterDefinition.__init__(self, name, description)
        self.default_value = default_value
        self.must_set_before_get = must_set_before_get


    def get_value(self, parameter, object):
        raise NotImplementedError()

    def set_value(self, parameter, object, quantity):
        raise NotImplementedError()
        
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
    
class InterfaceParameterDefinition(ParameterDefinition):
    def __init__(self, name, description, default_value,state_guard=None):
        AbstractParameterDefinition.__init__(self, name, description)
        self.default_value = default_value
        self.must_set_before_get = False
        self.value=default_value
        self.state_guard=state_guard
        
    def get_value(self, parameter, object):
        try:
          x=self.value.copy()
        except:
          x=self.value
        return x
        
    def set_value(self, parameter, object, quantity):
        if self.state_guard:
          getattr(object, self.state_guard)()
        try:
          self.value=quantity.copy()
        except:
          self.value=quantity

    def must_set_to_default_if_not_set(self):
        return False


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


class ModuleMethodParameterDefinition(ParameterDefinition):
    def __init__(self, get_method, set_method, name, description, default_value = None, must_set_before_get = False):
        ParameterDefinition.__init__(self, name, description, default_value, must_set_before_get)
        self.get_method = get_method
        self.set_method = set_method
        self.stored_value = None

    def get_value(self, parameter, object):
        if self.must_set_before_get and not parameter.is_set:
            self.set_default_value(parameter, object)
        
        if self.get_method is None:
            return self.stored_value
        else:
            return getattr(object, self.get_method)()

    def set_value(self, parameter, object, quantity):
        #if self.unit.is_non_numeric() or len(self.unit.base) == 0:
        #    if not isinstance(quantity, quantities.Quantity):
        #        quantity = quantity | self.unit
        
        if self.set_method is None:
            raise exceptions.CoreException("Could not set value for parameter '{0}' of a '{1}' object, parameter is read-only".format(self.name, type(object).__name__))
        
        getattr(object, self.set_method)(quantity)
        
        if self.get_method is None:
            self.stored_value = quantity
        
        parameter.is_set = True

    def is_readonly(self):
        return self.set_method is None


class ModuleBooleanParameterDefinition(ModuleMethodParameterDefinition):
    
    def __init__(self, *args, **kwargs):
        ModuleMethodParameterDefinition.__init__(self, *args, **kwargs)
    
    def get_value(self, parameter, object):
        return True if ModuleMethodParameterDefinition.get_value(self, parameter, object) else False
    
    def set_value(self, parameter, object, bool):
        return ModuleMethodParameterDefinition.set_value(self, parameter, object, 1 if bool else 0)


class ModuleCachingParameterDefinition(ParameterDefinition):
    def __init__(self, functionname, parameter_name, name, description, default_value = None):
        ParameterDefinition.__init__(self, name, description, default_value, must_set_before_get = True)
        self.stored_value = None
        self.parameter_name = parameter_name
        self.functionname = functionname

    def get_value(self, parameter, object):
        if self.must_set_before_get and not parameter.is_set:
            self.set_default_value(parameter, object)
            
        return parameter.cached_value

    def set_value(self, parameter, object, quantity):
        
        if is_quantity(self.default_value):
            unit = self.default_value.unit
            if unit.is_non_numeric() or len(unit.base) == 0:
                if not is_quantity(quantity):
                    quantity = quantity | unit
            
        parameter.cached_value = quantity
        
        parameter.is_set = True
        
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
        return self.definition.get_value(self, self.parameter_set._instance())
    
    
class ModuleVectorMethodParameterDefinition(ModuleMethodParameterDefinition):
    
    def get_value(self, parameter, object):
        if self.must_set_before_get and not parameter.is_set:
            self.set_default_value(parameter, object)
        
        if self.get_method is None:
            return self.stored_value
        else:
            list_of_scalars = getattr(object, self.get_method)()
            result = quantities.AdaptingVectorQuantity()
            result.extend(list_of_scalars)
            return result.copy()
            
            
    def set_value(self, parameter, object, vector_quantity):
        if self.set_method is None:
            raise exceptions.CoreException("Could not set value for parameter '{0}' of a '{1}' object, parameter is read-only".format(self.name, type(object).__name__))
        
        getattr(object, self.set_method)(*vector_quantity)
        
        if self.get_method is None:
            self.stored_value = vector_quantity
        
        parameter.is_set = True

    def is_readonly(self):
        return self.set_method is None
    

class VectorParameterDefinition(AbstractParameterDefinition):
    def __init__(self, name, description, names_of_parameters, default_value):
        AbstractParameterDefinition.__init__(self, name, description)
        self.names_of_parameters = names_of_parameters
        self.default_value = default_value

    def get_value(self, parameter, object):
        all_parameters = parameter.parameter_set
        result = []
        unit = None
        for name in self.names_of_parameters:
            parameter = all_parameters.get_parameter(name)
            element = parameter.get_value()
            if unit is None:
                if is_quantity(element):
                    unit = element.unit
            
            if not unit is None:
                result.append(element.value_in(unit))
            else:
                result.append(element)
               
        if not unit is None: 
            return unit.new_quantity(result)
        else:
            return numpy.asarray(result)
        
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
    
# to do: higher dimensional array parameters
class ModuleArrayParameterDefinition(ParameterDefinition):
    def __init__(self, get_method, set_method, range_method, name, description):
        ParameterDefinition.__init__(self, name, description, None, False)
        self.get_method = get_method
        self.set_method = set_method
        self.range_method = range_method
        self.stored_value = None

    def get_value(self, parameter, object):        
        if self.get_method is None:
            return self.stored_value
        else:
            irange=getattr(object, self.range_method)()
            index=numpy.arange(irange[0],irange[1]+1)
            return getattr(object, self.get_method)(index)

    def set_value(self, parameter, object, quantity):
        
        if self.set_method is None:
            raise exceptions.CoreException("Could not set value for parameter '{0}' of a '{1}' object, parameter is read-only".format(self.name, type(object).__name__))
        
        irange=getattr(object, self.range_method)()
        index=numpy.arange(irange[0],irange[1]+1)
        getattr(object, self.set_method)(index,quantity)
        
        if self.get_method is None:
            self.stored_value = quantity
        
        parameter.is_set = True

    def is_readonly(self):
        return self.set_method is None
