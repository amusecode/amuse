import weakref

from amuse.support.units import nbody_system



class Parameters(object):
    def __init__(self, definitions, instance):
        object.__setattr__(self, 'instance', weakref.ref(instance))
        object.__setattr__(self, 'definitions', definitions)
        object.__setattr__(self, 'mapping_from_name_to_definition', {})
        
        for x in definitions:
            self.mapping_from_name_to_definition[x.name] = x
    
    def __getattr__(self, name):
        if not name in self.mapping_from_name_to_definition:
            raise Exception("tried to get unknown parameter %s for module %s".format((name, type(instance).__name__)))
        return self.mapping_from_name_to_definition[name].get_value(self.instance())
    
    def __setattr__(self, name, value):
        if not name in self.mapping_from_name_to_definition:
            raise Exception("tried to set unknown parameter %s for module %s".format((name, type(instance).__name__)))
        self.mapping_from_name_to_definition[name].set_value(self.instance(), value)
        
class ParameterDefinition(object):
    def __init__(self, name, description, unit, default_value = None):
        self.name = name
        self.description = description
        self.unit = unit
        self.default_value = default_value
    
        
    def get_value(self, object):
        result = self.get_legacy_value(object) | self.unit
        if nbody_system.is_nbody_unit(self.unit):
            return object.convert_nbody.to_si(result)
        return result
        
    def set_value(self, object, quantity):
        if nbody_system.is_nbody_unit(self.unit):
            quantity = object.convert_nbody.to_nbody(quantity)
        self.set_legacy_value(object, quantity.value_in(self.unit))
        
    def set_default_value(self, object):
        if self.default_value:
            pass
        else:
            self.set_value(self, object, self.default_value)
        
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
        
    def get_legacy_value(self, object):
        return getattr(object, self.get_method)()
        
    def set_legacy_value(self, object, number):
        getattr(object, self.set_method)(number)
        

    
