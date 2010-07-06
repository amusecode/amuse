import weakref

from amuse.support.units import nbody_system
from amuse.support.units import generic_unit_system as generic_system
from amuse.support.data import values
from amuse.support import exception

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

    def __init__(self, definitions, instance):
        object.__setattr__(self, '_instance', weakref.ref(instance))
        object.__setattr__(self, '_definitions', definitions)
        object.__setattr__(self, '_mapping_from_name_to_definition', {})

        for x in definitions:
            self._mapping_from_name_to_definition[x.name] = x

    def __getattr__(self, name):
        if not name in self._mapping_from_name_to_definition:
            raise exception.CoreException("tried to get unknown parameter '{0}' for a '{1}' object".format(name, type(self._instance()).__name__))


        definition = self._mapping_from_name_to_definition[name]
        return definition.get_value(self._instance())

    def __setattr__(self, name, value):
        if not name in self._mapping_from_name_to_definition:
            warnings.warn("tried to set unknown parameter '{0}' for a '{1}' object".format(name, type(self._instance()).__name__), exception.AmuseWarning)
            return

        definition = self._mapping_from_name_to_definition[name]
        definition.set_value(self._instance(), value)
        if hasattr(self._instance(), "invoke_state_change"):
            self._instance().invoke_state_change()

    def names(self):
        return self._mapping_from_name_to_definition.keys()

    def set_defaults(self):
        for parameter_definition in self._definitions:
            parameter_definition.set_default_value(self._instance())

    def __dir__(self):
        result = []
        result.extend(dir(type(self)))
        result.extend(self.names())
        return result


    def get_default_value_for(self, name):
        if not name in self._mapping_from_name_to_definition:
            raise exception.CoreException("tried to get default value of unknown parameter '{0}' for a '{1}' object".format(name, type(self._instance()).__name__))

        definition = self._mapping_from_name_to_definition[name]
        return definition.default_value

    def __str__(self):
        output = ""

        for name in self.names():
            output += name + ": "
            output += str(getattr(self, name))+"\n"

        return output


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
        default_value = self._original.get_default_value_for(name)
        if not isinstance(default_value,bool) and (nbody_system.is_nbody_unit(default_value.unit) 
                or generic_system.is_generic_unit(default_value.unit)):
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

class ParameterDefinition(object):
    def __init__(self, name, description, unit, default_value = None):
        self.name = name
        self.description = description
        self.unit = unit
        self.default_value = default_value


    def get_value(self, object):
        result = self.unit.new_quantity(self.get_legacy_value(object))
        #if nbody_system.is_nbody_unit(self.unit):
        #    return object.convert_nbody.to_si(result)
        return result

    def set_value(self, object, quantity):
        if self.unit.is_non_numeric() or len(self.unit.base) == 0:
            if not isinstance(quantity, values.Quantity):
                quantity = quantity | self.unit

        #if nbody_system.is_nbody_unit(self.unit):
        #    quantity = object.convert_nbody.to_nbody(quantity)
        self.set_legacy_value(object, quantity.value_in(self.unit))

    def set_default_value(self, object):
        if self.default_value is None :
            return None

        if self.is_readonly():
            return None

        self.set_value(object, self.default_value)

    def is_readonly(self):
        return False

class ModuleAttributeParameterDefinition(ParameterDefinition):
    def __init__(self, attribute_name, name, description, unit, default_value = None):
        ParameterDefinition.__init__(self, name, description, unit, default_value)
        self.attribute_name = attribute_name

    def get_legacy_value(self, object):
        return getattr(object.legacy_interface, self.attribute_name)

    def set_legacy_value(self, object, number):
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
        if self.set_method is None:
            raise exception.CoreException("Could not set value for parameter '{0}' of a '{1}' object, parameter is read-only".format(self.name, type(object).__name__))

        error = getattr(object, self.set_method)(number)
        if error < 0:
            raise ParameterException(object, self.name, error, False)
        else:
            if self.get_method is None:
                self.stored_value = number



    def is_readonly(self):
        return self.set_method is None


class ModuleBooleanParameterDefinition(ParameterDefinition):
    def __init__(self, get_method, set_method, name, description, default_value = None):
        ParameterDefinition.__init__(self, name, description, None, default_value)
        self.get_method = get_method
        self.set_method = set_method
        self.stored_value = None

    def get_value(self, object):
        return True if self.get_legacy_value(object) else False

    def set_value(self, object, bool):
        self.set_legacy_value(object, 1 if bool else 0)

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
        if self.set_method is None:
            raise exception.CoreException("Could not set value for parameter '{0}' of a '{1}' object, parameter is read-only".format(self.name, type(object).__name__))

        error = getattr(object, self.set_method)(number)
        if error < 0:
            raise ParameterException(object, self.name, error, False)
        else:
            if self.get_method is None:
                self.stored_value = number

    def is_readonly(self):
        return self.set_method is None


class ModuleCachingParameterDefinition(ParameterDefinition):
    def __init__(self, parameter_name, name, description, unit, default_value = None):
        ParameterDefinition.__init__(self, name, description, unit, default_value)
        self.stored_value = None
        self.parameter_name = parameter_name

    def get_legacy_value(self, object):
        if hasattr(object, "_parameters"):
            if self.name in object._parameters:
                return object._parameters[self.name]
            else:
                return self.stored_value
        else:
            return self.stored_value

    def set_legacy_value(self, object, number):
        if hasattr(object, "_parameters"):
            object._parameters[self.name] = number
        else:
            object._parameters = {self.name: number}
