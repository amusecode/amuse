import weakref

from amuse.support.units import nbody_system
from amuse.support.data import values



class Parameters(object):
    def __init__(self, definitions, instance):
        object.__setattr__(self, '_instance', weakref.ref(instance))
        object.__setattr__(self, '_definitions', definitions)
        object.__setattr__(self, '_mapping_from_name_to_definition', {})
        
        for x in definitions:
            self._mapping_from_name_to_definition[x.name] = x
    
    def __getattr__(self, name):
        if not name in self._mapping_from_name_to_definition:
            raise AttributeError("tried to get unknown parameter '{0}' for a '{1}' object".format(name, type(self._instance()).__name__))
        
        definition = self._mapping_from_name_to_definition[name]
        return definition.get_value(self._instance())
    
    def __setattr__(self, name, value):
        if not name in self._mapping_from_name_to_definition:
            raise AttributeError("tried to set unknown parameter '{0}' for a '{1}' object".format(name, type(self._instance()).__name__))
        
        definition = self._mapping_from_name_to_definition[name]
        definition.set_value(self._instance(), value)
        
    def names(self):
        return self._mapping_from_name_to_definition.keys()
        
    def set_defaults(self):
        for parameter_definition in self._definitions:
            parameter_definition.set_default_value(self._instance())
            
        
class ParameterDefinition(object):
    def __init__(self, name, description, unit, default_value = None):
        self.name = name
        self.description = description
        self.unit = unit
        self.default_value = default_value
    
    
    def get_value(self, object):
        result = self.unit.new_quantity(self.get_legacy_value(object))
        if nbody_system.is_nbody_unit(self.unit):
            return object.convert_nbody.to_si(result)
        return result
        
    def set_value(self, object, quantity):
        if self.unit.is_non_numeric():
            if not isinstance(quantity, values.Quantity):
                quantity = quantity | self.unit          
                  
        if nbody_system.is_nbody_unit(self.unit):
            quantity = object.convert_nbody.to_nbody(quantity)
        self.set_legacy_value(object, quantity.value_in(self.unit))
        
    def set_default_value(self, object):
        if self.default_value is None:
            pass
        else:
            self.set_value(object, self.default_value)
        
class ModuleAttributeParameterDefinition(ParameterDefinition):
    def __init__(self, attribute_name, name, description, unit, default_value = None):
        ParameterDefinition.__init__(self, name, description, unit, default_value)
        self.attribute_name = attribute_name
        
    def get_legacy_value(self, object):
        return getattr(object, self.attribute_name)
        
    def set_legacy_value(self, object, number):
        setattr(object, self.attribute_name, number)
        
class ModuleMethodParameterDefinition(ParameterDefinition):
    def __init__(self, get_method, set_method, name, description, unit, default_value = None):
        ParameterDefinition.__init__(self, name, description, unit, default_value)
        self.get_method = get_method
        self.set_method = set_method
        self.stored_value = None
        
        
    def get_legacy_value(self, object):
        if self.get_method is None:
            return self.stored_value
        else:
            return getattr(object, self.get_method)()
        
    def set_legacy_value(self, object, number):
        getattr(object, self.set_method)(number)
        if self.get_method is None:
            self.stored_value = number
            
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
    def __init__(self, get_method, set_method, name, description, unit, default_value = None):
        ParameterDefinition.__init__(self, name, description, unit, default_value)
        self.get_method = get_method
        self.set_method = set_method
        self.stored_value = None
        
        
    def get_legacy_value(self, object):
        if self.get_method is None:
            return self.stored_value
        else:
            (result, error) = getattr(object, self.get_method)()
            if error < 0:
                raise ParameterException(object, self.name, error, True)
            else:
                return result
        
    def set_legacy_value(self, object, number):
        error = getattr(object, self.set_method)(number)
        if error < 0:
            raise ParameterException(object, self.name, error, False)
        else:
            if self.get_method is None:
                self.stored_value = number
                
                
        
class ModuleCachingParameterDefinition(ParameterDefinition):
    def __init__(self, parameter_name, name, description, unit, default_value = None):
        ParameterDefinition.__init__(self, name, description, unit, default_value)
        self.stored_value = None
        self.parameter_name = parameter_name
        
        
    def get_legacy_value(self, object):
        return self.stored_value
        
    def set_legacy_value(self, object, number):
        self.stored_value = number
        

    
