from amuse.support.data import parameters
from amuse.support.data.core import Particles, ParticleInformationChannel, Particle
import numpy

from amuse.support.units import nbody_system




class CodeMethod(object):
    """
    CodeMethod objects define an interface between the codes and
    the AMUSE data model. Abstract base class, subclasses implement
    actual implementation.
    
    The code provides functions, each function has input and output 
    parameters. These parameters are defines as scalar or vector values 
    (no units). Codes may or may not define objects.
    
    The AMUSE data model works with objects (particles) 
    and attributes (each attribute has a value and a unit).
    
    CodeMethod instances map the functions to the objects in the 
    datamodel.
    
    :argument name_of_the_function: name of the function 
        provided by the code. for example ""get_mass""
    :argument attribute_specification: list of tuples. each tuple
        consists of the name of the attribute of a particle, the name
        of the parameter in the function, and the unit of the parameter 
        in the function. Example specification::
            (
                ("mass", "m", unit.kg),
                ("radius", "r", unit.m),
            )
    
 
    """
    def __init__(self, name_of_the_function, attribute_specification):
        self.name_of_the_function = name_of_the_function
        self.attribute_specification = attribute_specification
        
        self.mapping_from_parameter_name_to_attribute_and_unit = self.to_mapping(1, self.attribute_specification)
        self.mapping_from_attribute_name_to_parameter_and_unit = self.to_mapping(0, self.attribute_specification)
            
    def attributes(self):
        """
        The list of attributes supported by this method. The method will
        transfer values of these attributes in one call.
        """
        return self.mapping_from_attribute_name_to_parameter_and_unit.keys()
    
    def to_mapping(self, index, list_of_tuples):
        result = {}
        for tuple in list_of_tuples:
            mutable = list(tuple)
            del mutable[index]
            result[tuple[index]] = mutable
        return result
            
class ParticleAttributesModifier(CodeMethod):
    """
    Objects of this class update the values of attributes in a particle set
    using a function provided by the code.
    
    Objects of this class are used when a code does not have a concept
    of a particle as something that is stored in memory, 
    for example in some stellar evolution codes. These codes instead
    provide function with a set of input and output parameters. Objects
    of this class connect one function with a set of input an output
    attributes. Objects of this class only work on functions specified
    as a *legacy_function* (this specification is used to determine
    the input and output parameters).
    """
    def __init__(self, name_of_the_function, attribute_specification, parameter_specification = ()):
        CodeMethod.__init__(self, name_of_the_function, attribute_specification)
        self.parameter_specification = parameter_specification
        
        self.mapping_from_input_name_to_parameter_name_and_unit = self.to_mapping(0, self.parameter_specification)

            
    def _run(self, code_interface, particles, *list_arguments, **input_keyword_arguments):
            
        function = getattr(code_interface, self.name_of_the_function)
        
        legacy_function_specification = function.specification
        keyword_arguments = {}
        for parameter in legacy_function_specification.parameters:
            if not parameter.name in self.mapping_from_parameter_name_to_attribute_and_unit:
                continue
            attribute_name, unit = self.mapping_from_parameter_name_to_attribute_and_unit[parameter.name]
            if parameter.is_input():
                keyword_arguments[parameter.name] = getattr(particles, attribute_name).value_in(unit)
        
        for index, argument in enumerate(list_arguments):
            input_name, parameter_name, unit = self.parameter_specification[index]
            keyword_arguments[parameter_name] = argument.as_vector_with_length(len(particles)).value_in(unit)
        
        for input_name, argument in input_keyword_arguments.iteritems():
            parameter_name, unit = self.mapping_from_input_name_to_parameter_name_and_unit[input_name]
            keyword_arguments[parameter_name] = argument.value_in(unit)
        
        result_dictionary = function(**keyword_arguments)
        
        for parameter in legacy_function_specification.parameters:
            if not parameter.name in self.mapping_from_parameter_name_to_attribute_and_unit:
                continue  
            attribute_name, unit = self.mapping_from_parameter_name_to_attribute_and_unit[parameter.name]
            if parameter.is_output():
                if not hasattr(result_dictionary, "keys"):
                    setattr(particles, attribute_name, unit.new_quantity(result_dictionary))
                else:
                    setattr(particles, attribute_name, unit.new_quantity(result_dictionary[parameter.name]))
    
        
            
class ParticleGetAttributesMethod(CodeMethod):
    """
    Objects of this class retrieve the values of the specified attributes
    from the code.
    
    Objects of this class are used to provide the back-bone of the
    the InCodeAttributeStorage objects.
    
    The function provided by the code must comply with a specific
    interface.
    For python interface functions, these must take one input
    parameter (indices of the particles)
    and return a dictionary containing a result and
    one or more ouput parameter containing the values.
    
    * one an input parameter specifying the particle indices.
    
    * one list of returncodes specificying error conditions 
      (negative for errors, zero for normal operation). This
      value must be stored under "__result" in the result dictionary.
    * one or more output parameters with the values of the parameters
      This value must be stored under the name specified
      in the given specification.
    
    An example python function is::
    
        def get_mass(self, indices):
            masses = []
            errors = []
            for i in indices:
                masses.append(i * 10.0)
                errors.append(0)
            return {
                "mass": masses,
                "__result": errors,
            }
    
    
    
    For legacy interface functions, the
    specification must be like (note: can_handle_array must be True)::
    
        @legacy_function
        def get_mass():
            function = LegacyFunctionSpecification()  
            function.addParameter(
                'index_of_the_particle',
                dtype='int32', 
                direction=function.IN,
                description = "Index of the particle to get the value of the attribute from. This index must have been returned by an earlier call to :meth:`new_particle`")
            function.addParameter(
                'parameter_name', 
                dtype='float64',
                direction=function.OUT, 
                description = "The current value of the parameter of the particle")
            function.result_type = 'int32'
            function.can_handle_array = True 
            return function
        
    
    
    """
    def __init__(self, name_of_the_function, attribute_specification):
        CodeMethod.__init__(self, name_of_the_function, attribute_specification)
        
    
    def intersection(self, attributes):
        """
        Determines the intersection between the attributes specified for this
        object and the given attributes. A get method object has an intersection
        when any of the specified attributes are in given in the parameter.
        
        Intersections are used to determine if a method can be used to get values
        of some attributes. Get functions can be used when some of the parameters are
        given. If you define a function to get the position of a particle returning
        x, y and z then you may query only the x (the y and z values will be
        sent from the code to python but these will be ignored).
        """
        result = set([])
        for x in attributes:
            if x in self.mapping_from_attribute_name_to_parameter_and_unit:
                result.add(x)
        return result
    
    def apply(self, code_interface, indices, attributes):
        """
        Implements the retrieval from the code of values of the specified attributes 
        for particles with the given indices.
        
        >>> from amuse.support.units import si
        >>> method = ParticleGetAttributesMethod("get_position", (
        ...     ("x", "X", si.m),
        ...     ("y", "Y", si.m),
        ...     ("z", "Z", si.m)
        ...     )
        ... )
        >>> class Example(object):
        ...     def get_position(self, indices):
        ...         x = []
        ...         y = []
        ...         z = []
        ...         errors = []
        ...         for i in indices:
        ...             x.append(i)
        ...             y.append(i + 1.0)
        ...             z.append(i * 2.0)
        ...             errors.append(0)
        ...         return {
        ...             "X": x,
        ...             "Y": y,
        ...             "Z": z,
        ...             "__result": errors,
        ...         }
        ...
        >>> method.apply(Example(), [1.0,2.0,3.0], ["x","z"])
        {'x': quantity<[1.0, 2.0, 3.0] m>, 'z': quantity<[2.0, 4.0, 6.0] m>}
        """
        method = getattr(code_interface, self.name_of_the_function)
        
        if len(indices) > 1: 
            if hasattr(method, "specification") and not method.specification.can_handle_array:
                raise Exception(
                    "getter method <{0}> of a '{1}' object, cannot handle arrays".format(method_name, type(self).__name__)
                )
            
        keyword_results = method(indices)
        
        set_of_attributes = set(attributes)
        mapping_from_attribute_to_result = {}  
        for attribute in self.mapping_from_attribute_name_to_parameter_and_unit.keys():
            keyword, unit = self.mapping_from_attribute_name_to_parameter_and_unit[attribute]
            if attribute in set_of_attributes:
                mapping_from_attribute_to_result[attribute] = unit.new_quantity(keyword_results[keyword])
        
        return mapping_from_attribute_to_result
    
class ParticleSetAttributesMethod(CodeMethod):
    """
    Objects of this class set the values of the specified attributes
    in the code.
    
    Objects of this class are used to provide the back-bone of the
    the InCodeAttributeStorage objects.
    
    The function provided by the code must comply with a specific
    interface.
    
    For python interface functions, these must take a input
    parameter for the indices of the particles and one or more
    input parameters for the values of the attributes. It must resturn
    and return a ordered dictionary containing a list of resultcodes
    
    
    An example python function is::
        
        def set_mass(self, indices, masses):
            for i, mass in zip(indices,masses):
                self.masses[i] = mass
                
            return {
                "__result": [0] * len(indices)
            }
    
    
    For legacy interface functions, the
    specification must be like (note: can_handle_array must be True)::
    
        @legacy_function
        def set_mass():
            function = LegacyFunctionSpecification()  
            function.addParameter(
                'index_of_the_particle',
                dtype='int32',
                direction=function.IN,
                description = "Index of the particle for which the state is to be updated. This index must have been returned by an earlier call to :meth:`new_particle`"
            )
            function.addParameter(
                'mass', 
                dtype='float64', 
                direction=function.IN, 
                description = "The new mass of the particle"
            )
            function.result_type = 'int32'
            function.can_handle_array = True
            return function    
    """
    def __init__(self, name_of_the_function, attribute_specification):
        CodeMethod.__init__(self, name_of_the_function, attribute_specification)
            
    def intersection(self, attributes):
        """
        Determines the intersection between the attributes specified for this
        object and the given attributes. A set method object only has an intersection
        when all the specified attributes are in given in the parameter.
        
        Intersections are used to determine if a method can be used to set values
        of some attributes. Set function can only be used when all parameters are
        given. If you define a function to set the position of a particle using
        x, y and z then you need to set the values in one go and you need to
        prove values for the x y and z.
        """
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
        

class InterfaceWithObjectsBinding(object):
    def __init__(self):
        self.mapping_from_particleid_to_index = {}
        
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
        
        
        
        

class InCodeAttributeStorage(object):
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
            self.attributes |= set(x.attributes())
        for x in self.setters:
            self.attributes |= set(x.attributes())
        
        
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
            self.particle_keys = numpy.concatenate((self.particle_keys, numpy.array(list(keys))))
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
