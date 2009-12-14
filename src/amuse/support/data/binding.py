from amuse.support.data import parameters
from amuse.support.data.core import Particles, ParticleInformationChannel, Particle
import numpy

from amuse.support.units import nbody_system




class CodeMethod(object):
    
    def __init__(self, name_of_the_function, attribute_specification):
        self.name_of_the_function = name_of_the_function
        self.attribute_specification = attribute_specification
        self.mapping_from_parameter_name_to_attribute_and_unit = {}
        
        for attribute_name, parameter_name, unit in self.attribute_specification:
            self.mapping_from_parameter_name_to_attribute_and_unit[parameter_name] = (attribute_name, unit)
            
        self.mapping_from_attribute_name_to_parameter_and_unit = {}
        for attribute_name, parameter_name, unit in self.attribute_specification:
            self.mapping_from_attribute_name_to_parameter_and_unit[attribute_name] = (parameter_name, unit)
            
    def attributes(self):
        return self.mapping_from_attribute_name_to_parameter_and_unit.keys()
            
class ParticleAttributesModifier(CodeMethod):
    
    def __init__(self, name_of_the_function, attribute_specification):
        CodeMethod.__init__(self, name_of_the_function, attribute_specification)
            
    def _run(self, code_interface, particles):
            
        function = getattr(code_interface, self.name_of_the_function)
        legacy_function_specification = function.specification
        keyword_arguments = {}
        for parameter in legacy_function_specification.parameters:  
            attribute_name, unit = self.mapping_from_parameter_name_to_attribute_and_unit[parameter.name]
            if parameter.is_input():
                keyword_arguments[parameter.name] = getattr(particles, attribute_name).value_in(unit)
        
        result_dictionary = function(**keyword_arguments)
        
        for parameter in legacy_function_specification.parameters:  
            attribute_name, unit = self.mapping_from_parameter_name_to_attribute_and_unit[parameter.name]
            if parameter.is_output():
                setattr(instance.particles, attribute_name, unit.new_quantity(result_dictionary[parameter.name]))
        
            
class ParticleGetAttributesMethod(CodeMethod):
    
    def __init__(self, name_of_the_function, attribute_specification):
        CodeMethod.__init__(self, name_of_the_function, attribute_specification)
        
    
    def intersection(self, attributes):
        result = set([])
        for x in attributes:
            if x in self.mapping_from_attribute_name_to_parameter_and_unit:
                result.add(x)
        return result
    
    def apply(self, code_interface, indices, attributes):
        method = getattr(code_interface, self.name_of_the_function)
        
        if len(indices) > 1: 
            if hasattr(method, "specification") and not method.specification.can_handle_array:
                raise Exception(
                    "getter method <{0}> of a '{1}' object, cannot handle arrays".format(method_name, type(self).__name__)
                )
            
        keyword_results = method(indices)
          
        mapping_from_attribute_to_result = {}  
        for attribute in self.mapping_from_attribute_name_to_parameter_and_unit.keys():
            keyword, unit = self.mapping_from_attribute_name_to_parameter_and_unit[attribute]
            mapping_from_attribute_to_result[attribute] = unit.new_quantity(keyword_results[keyword])
        
        return mapping_from_attribute_to_result
    
class ParticleSetAttributesMethod(CodeMethod):
    
    def __init__(self, name_of_the_function, attribute_specification):
        CodeMethod.__init__(self, name_of_the_function, attribute_specification)
            
    def intersection(self, attributes):
        result = set([])
        for x in attributes:
            if x in self.mapping_from_attribute_name_to_parameter_and_unit:
                result.add(x)
        if len(result) != len(self.mapping_from_attribute_name_to_parameter_and_unit):
            return set([])
        return result
        
    
    def apply(self, code_interface,  indices, attributes = [], values = []):

        mapping = self.mapping_from_attribute_name_to_parameter_and_unit
        keyword_arguments = {}
        for attribute, quantity in zip(attributes, values):
            if attribute in mapping:
                keyword, unit = mapping[attribute]
                keyword_arguments[keyword] = quantity.value_in(unit)
        keyword_arguments['id'] = indices
        keyword_arguments['index_of_the_particle'] = indices
             
        method = getattr(code_interface, self.name_of_the_function)
        
        errors = method(**keyword_arguments)
        
        for errorcode in errors:
            if errorcode < 0:
                raise Exception("Could not add a particle")


class NewParticleMethod(CodeMethod):
    
    def __init__(self, name_of_the_function, attribute_specification):
        CodeMethod.__init__(self, name_of_the_function, attribute_specification)

    def apply(self, code_interface,  attributes = [], values = []):
        
        mapping = self.mapping_from_attribute_name_to_parameter_and_unit
        keyword_arguments = {}
        for attribute, quantity in zip(attributes, values):
            if attribute in mapping:
                keyword, unit = mapping[attribute]
                keyword_arguments[keyword] = quantity.value_in(unit)
                
        method = getattr(code_interface, self.name_of_the_function)
        
        indices, errors = method(**keyword_arguments)
        
        for errorcode in errors:
            if errorcode < 0:
                raise Exception("Could not add a particle")
        
        return indices
                
                
class InterfaceWithParametersBinding(object):
    parameter_definitions = []
    
    def __init__(self, convert_nbody = None):
               
        self.parameters = parameters.Parameters(self.parameter_definitions, self)
        

class InCodeAttributeStorage(object):
    
    
    def __init__(self, code_interface):
        self.code_interface = code_interface
        self.mapping_from_particle_key_to_index_in_the_code = {}
        self.particle_keys = []
        self.attributes = set([])
            
    def __len__(self):
        length, error = self.code_interface.get_number_of_particles()
        return length
        
    def complies_with_state(self, attributes):
        attributes_in_getstate = self._mapping_from_attribute_names_to_set_keyword_and_unit().keys()
        print attributes_in_getstate
        if not attributes_in_getstate:
            return False
        
        if len(attributes_in_getstate) == len(attributes):
            check = set(attributes)
            for x in attributes_in_getstate:
                if x in check:
                    check.remove(x)
                else:
                    return False
            return len(check) == 0
        else:
            return False
            
    def select_getters_for(self, attributes):
        attributes_in_getstate = self._mapping_from_attribute_names_to_set_keyword_and_unit().keys()
        if self.complies_with_state(attributes):
            return [("get_state", self._mapping_from_attribute_names_to_set_keyword_and_unit())]
        set_of_attributes = set(attributes)
        
        result = []
        for x in self.code_interface.attribute_definitions:
            y = x.attributes_in_getter(set_of_attributes)
            if not y is None:
                result.append(y)
                
        if set_of_attributes:
            raise Exception("Do not have attributes {0}".format(sorted(set_of_attributes)))
        return result
        
    
    def select_setters_for(self, attributes):
        attributes_in_getstate = self._mapping_from_attribute_names_to_set_keyword_and_unit().keys()
        if self.complies_with_state(attributes):
            return [("set_state", self._mapping_from_attribute_names_to_set_keyword_and_unit())]
        set_of_attributes = set(attributes)
        
        result = []
        for x in self.code_interface.attribute_definitions:
            y = x.attributes_in_setter(set_of_attributes)
            if not y is None:
                result.append(y)
                
        if set_of_attributes:
            raise Exception("Cannot set attributes {0}".format(sorted(set_of_attributes)))
        return result
        
    def _mapping_from_attribute_names_to_set_keyword_and_unit(self):
        result = {}
        for attribute_definition in self.code_interface.attribute_definitions:
            if not attribute_definition.is_required_for_setup():
                continue
            result.update(attribute_definition.state_mapping_from_name_to_keyword_and_unit())
        return result
            
    def _set_particles(self, keys, attributes = [], values = []):
        
        mapping = self._mapping_from_attribute_names_to_set_keyword_and_unit()
        keyword_arguments = {}
        for attribute, quantity in zip(attributes, values):
            if attribute in mapping:
                keyword, unit = mapping[attribute]
                unit = self.code_interface.convert_to_nbody(unit)
                keyword_arguments[keyword] = quantity.value_in(unit)
    
        indices, errors = self.code_interface.new_particle(**keyword_arguments)
        
        
        for errorcode in errors:
            if errorcode < 0:
                raise Exception("Could not setup a particle")
            
        if len(self.particle_keys) > 0:
            self.particle_keys = numpy.concatenate((self.particle_keys, keys))
        else:
            self.particle_keys = numpy.array(keys)

        d = self.mapping_from_particle_key_to_index_in_the_code
        index = 0
        for key in keys:
            d[key] = indices[index]
            index = index + 1
        

            
            

    def _get_values(self, keys, attributes):
        indices_in_the_code = []
        if keys is None:
            keys = self.particle_keys
        for particle_key in keys:
            indices_in_the_code.append(self.mapping_from_particle_key_to_index_in_the_code[particle_key])
        
        mapping_from_attribute_to_result = {}
        
        for method_name, mapping in self.select_getters_for(attributes):
            keyword_to_attribute_and_unit_mapping = {}
            
            
            method = getattr(self.code_interface, method_name)
            
            if len(indices_in_the_code) > 1: 
                if hasattr(method, "specification") and not method.specification.can_handle_array:
                    raise Exception(
                        "getter method <{0}> of a '{1}' object, cannot handle arrays".format(method_name, type(self).__name__)
                    )
            
            keyword_results = method(indices_in_the_code)
            
            for attribute in mapping.keys():
                keyword, unit = mapping[attribute]
                unit = self.code_interface.convert_to_nbody(unit)
                mapping_from_attribute_to_result[attribute] = unit.new_quantity(keyword_results[keyword])
        
        results = []
        for attribute in attributes:
            results.append(mapping_from_attribute_to_result[attribute])
        return results
        
                
        
    def _set_values(self, keys, attributes, values):
        indices_in_the_code = []
        if keys is None:
            keys = self.particle_keys
            
        for particle_key in keys:
            indices_in_the_code.append(self.mapping_from_particle_key_to_index_in_the_code[particle_key])
            
        for method_name, mapping in self.select_setters_for(attributes):
            keyword_arguments = {}
            
            for attribute, quantity in zip(attributes, values):
                keyword, unit = mapping[attribute]
                unit = self.code_interface.convert_to_nbody(unit)
                keyword_arguments[keyword] = quantity.value_in(unit)
        
            
            method = getattr(self.code_interface, method_name)
            
            if len(indices_in_the_code) > 1: 
                if hasattr(method, "specification") and not method.specification.can_handle_array:
                    raise Exception(
                        "setter method <{0}> of a '{1}' object, cannot handle arrays".format(method_name, type(self.code_interface).__name__)
                    )
                    
            keyword_arguments['id'] = indices_in_the_code
            keyword_arguments['index_of_the_particle'] = indices_in_the_code
            
            errors = method(**keyword_arguments)
            
    def _remove_particles(self, keys):
        indices_in_the_code = []
        if keys is None:
            keys = self.particle_keys
            
        for particle_key in keys:
            indices_in_the_code.append(self.mapping_from_particle_key_to_index_in_the_code[particle_key])
            
        errors = self.code_interface.delete_particle(indices_in_the_code)
        for errorcode in errors:
            if errorcode < 0:
                raise Exception("Could not remove a particle")
        
        d = self.mapping_from_particle_key_to_index_in_the_code
        index = 0
        for key in keys:
            del d[key]
            
        self.particle_keys =  numpy.delete(self.particle_keys, indices_in_the_code)
            
    def _get_attributes(self):
        return self._mapping_from_attribute_names_to_set_keyword_and_unit().keys()
    
    
    def _get_keys(self):
        return self.particle_keys

    def _has_key(self, key):
        return key in self.mapping_from_particle_key_to_index_in_the_code
    
    def _state_attributes(self):
        return self._mapping_from_attribute_names_to_set_keyword_and_unit().keys()


class InterfaceWithObjectsBinding(object):
    def __init__(self):
        self.mapping_from_particleid_to_index = {}
        self.particles = Particles()
        self.particles.attributelist = InCodeAttributeStorage(self)
    
    def convert_to_nbody(self, x):
        if nbody_system.is_nbody_unit(x):
            return self.convert_nbody.unit_to_unit_in_si(x)
        else:
            return x
                
    def setup_particles(self, particles):
        self.particles.add_particles(particles)
        
    def update_particles(self, particles):
        self.particles.copy_values_of_state_attributes_to(particles)
    
    def set_attribute(self, attribute_name, particles):
        particles.copy_values_of_attribute_to(attribute_name, self.particles)
        
    def update_attribute(self, attribute_name, particles):
        self.particles.copy_values_of_attribute_to(attribute_name, particles) 
    
    def get_attribute_definition(self, attribute_name):
        for attribute_definition in self.attribute_definitions:
            if attribute_definition.name == attribute_name:
                return attribute_definition
        return None
            
    def current_model_time(self):
        raise AttributeError("Must implement current_model_time method")
        
        
        
        

class InCodeAttributeStorage2(object):
    name_of_number_of_particles_getter = "get_number_of_particles"
    new_particle_method = NewParticleMethod("new_particle", ())
    name_of_delete_particle = "delete_particle"
    getters = ()
    setters = ()
    
    
    def __init__(self, code_interface):
        self.code_interface = code_interface
        self.mapping_from_particle_key_to_index_in_the_code = {}
        self.particle_keys = []
        
        self._get_number_of_particles = getattr(self.code_interface, self.name_of_number_of_particles_getter)
        self._delete_particle = getattr(self.code_interface, self.name_of_delete_particle)
        
        self.attributes = set([])
        for x in self.getters:
            self.attributes = set(x.attributes())
        for x in self.setters:
            self.attributes = set(x.attributes())
        
        
    def __len__(self):
        length, error = self._get_number_of_particles()
        return length
            
    def select_getters_for(self, attributes):
        set_of_attributes = set(attributes)
        result = []
        for particle_method in self.getters:
            provided_attributes = particle_method.intersection(set_of_attributes)
            if provided_attributes:
                result.append(particle_method)
                set_of_attributes -= provided_attributes
            
        if set_of_attributes:
            raise Exception("Do not have attributes {0}".format(sorted(set_of_attributes)))
        return result
        
    
    def select_setters_for(self, attributes):
        set_of_attributes = set(attributes)
        result = []
        for particle_method in self.setters:
            provided_attributes = particle_method.intersection(set_of_attributes)
            if provided_attributes:
                result.append(particle_method)
                set_of_attributes -= provided_attributes
                
        if set_of_attributes:
            raise Exception("Cannot set attributes {0}".format(sorted(set_of_attributes)))
            
        return result
        
    def _mapping_from_attribute_names_to_set_keyword_and_unit(self):
        result = {}
        for attribute_definition in self.code_interface.attribute_definitions:
            if not attribute_definition.is_required_for_setup():
                continue
            result.update(attribute_definition.state_mapping_from_name_to_keyword_and_unit())
        return result
            
    def _set_particles(self, keys, attributes = [], values = []):
        
        indices = self.new_particle_method.apply(self.code_interface, attributes, values)
            
        if len(self.particle_keys) > 0:
            self.particle_keys = numpy.concatenate((self.particle_keys, keys))
        else:
            self.particle_keys = numpy.array(keys)

        d = self.mapping_from_particle_key_to_index_in_the_code
        index = 0
        for key in keys:
            d[key] = indices[index]
            index = index + 1
        
    def get_indices_of(self, keys):
        indices_in_the_code = []
        if keys is None:
            keys = self.particle_keys
            
        for particle_key in keys:
            indices_in_the_code.append(self.mapping_from_particle_key_to_index_in_the_code[particle_key])
        return indices_in_the_code
        
    def _get_values(self, keys, attributes):
        indices_in_the_code = self.get_indices_of(keys)
        
        mapping_from_attribute_to_result = {}
        
        for getter in self.select_getters_for(attributes):
            result = getter.apply(self.code_interface, indices_in_the_code, attributes)
            mapping_from_attribute_to_result.update(result)
            
        results = []
        for attribute in attributes:
            results.append(mapping_from_attribute_to_result[attribute])
        return results
        
    def _set_values(self, keys, attributes, values):
        indices_in_the_code = self.get_indices_of(keys)
    
        for setter in self.select_setters_for(attributes):
            setter.apply(self.code_interface, indices_in_the_code, attributes, values)
    
    def _remove_particles(self, keys):
        indices_in_the_code = []
        if keys is None:
            keys = self.particle_keys
            
        for particle_key in keys:
            indices_in_the_code.append(self.mapping_from_particle_key_to_index_in_the_code[particle_key])
        
        method = getattr(self.code_interface, self.name_of_delete_particle)
        errors = method(indices_in_the_code)
        for errorcode in errors:
            if errorcode < 0:
                raise Exception("Could not remove a particle")
        
        d = self.mapping_from_particle_key_to_index_in_the_code
        index = 0
        for key in keys:
            del d[key]
            
        self.particle_keys =  numpy.delete(self.particle_keys, indices_in_the_code)
            
    def _get_attributes(self):
        return self.attributes
    
    def _state_attributes(self):
        return self._get_attributes()
        
    def _get_keys(self):
        return self.particle_keys

    def _has_key(self, key):
        return key in self.mapping_from_particle_key_to_index_in_the_code
    
    


        
        
        
class CodeProperty(object):
    
    def __init__(self, function_or_attribute_name, unit):
        self.function_or_attribute_name = function_or_attribute_name
        self.unit = unit
    
    def _name(self, instance):
        class_of_the_instance = type(instance)
        attribute_names = dir(class_of_the_instance)
        for name in attribute_names:
            if not name.startswith('_'):
                if getattr(class_of_the_instance, name) == self:
                    return name
        
        return '<unknown name>'

    def __get__(self, instance, owner):
        if instance is None:
            return self
        
        function_or_attribute = getattr(instance, self.function_or_attribute_name)
        if hasattr(function_or_attribute, '__call__'):
            value, errorcode = function_or_attribute()
            if errorcode < 0:
                raise Exception("calling '{0}' to get the value for property '{1}' resulted in an error (errorcode {2})".format(self.function_or_attribute_name, self._name(instance), errorcode))
            else:
                return self.unit.new_quantity(value)
        else:
            return self.unit.new_quantity(function_or_attribute)
            
    def __set__(self, instance, value):
        if instance is None:
            return self
            
        raise Exception("property '{0}' of a '{1}' object is read-only, you cannot change it's value".format(self._name(instance), type(instance).__name__))
