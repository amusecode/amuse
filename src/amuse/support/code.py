from amuse.support import interface
from amuse.support.options import OptionalAttributes
from amuse.support.core import late
from amuse import datamodel

import types
import inspect
import re

def _getdoc(object):
    """Get the doc string or comments for an object."""
    result = inspect.getdoc(object) or inspect.getcomments(object)
    return result and re.sub('^ *\n', '', result.rstrip()) or ''
    
class StateMethodDescriptor(object):
    
    def __init__(self, states, function):
        self.states = states
        self.function = function
        
        state_doc = "STATE: can run in states {0}".format(self.states)
        doc_string = _getdoc(function)
        print(doc_string)
        if doc_string:
            doc_string += '\n' 
        
        doc_string += state_doc
        self.function.__doc__ =   doc_string
            
    def __get__(self, instance, owner):
        if instance is None:
            return self.function
        else:
            method = types.MethodType(self.function, instance)
            name = self.get_name()
            function_with_state = instance.state_handler.get_attribute(name , method)
            setattr(instance, name, function_with_state)
            return function_with_state
        
    def get_name(self):
        if isinstance(self.function,StateMethodDescriptor):
            return self.function.get_name()
        else:
            return self.function.__name__
    
    
    def define_on(self, handler):
        for state in self.states:
            handler.add_method(state, self.function.__name__)
        
class StateTransitionDescriptor(object):
    
    def __init__(self, from_state, to_state, function):
        self.from_state = from_state
        self.to_state = to_state
        self.function = function
        
    def __get__(self, instance, owner):
        if instance is None:
            return self.get_function()
        else:
            function = self.get_function()
                
            method = types.MethodType(function, instance)
            name = self.get_name()
            function_with_state = instance.state_handler.get_attribute(name , method)
            setattr(instance, name, function_with_state)
            return function_with_state
    
    def get_name(self):
        if isinstance(self.function,StateTransitionDescriptor):
            return self.function.get_name()
        else:
            return self.function.__name__
            
    def get_function(self):
        if isinstance(self.function,StateTransitionDescriptor):
            return self.function.get_function()
        else:
            return self.function
    
        
    def define_on(self, handler):
        handler.add_transition(self.from_state, self.to_state, self.get_name())
        if isinstance(self.function,StateTransitionDescriptor):
            self.function.define_on(handler)
            
def state_method(*arguments):
    def f(func) :
        return StateMethodDescriptor(arguments, func)
    return f
    
def state_transition(from_state, to_state):
    def f(func) :
        return StateTransitionDescriptor(from_state, to_state, func)
    return f    
    
class MetaObjectWithState(type):
    def __new__(mcs, name, bases, dict):
        definitions = []
        for x in bases:
            if hasattr(x, '__state_definitions__'):
                definitions.extend(getattr(x, '__state_definitions__'))
                
        for key, value in dict.items():
            if isinstance(value,StateTransitionDescriptor):
                definitions.append(value)
            elif isinstance(value,StateMethodDescriptor):
                definitions.append(value)
        
            
        dict['__state_definitions__'] = definitions
        result = type.__new__(mcs, name, bases, dict)
        return result
    
    
class ObjectWithState(object, metaclass=MetaObjectWithState):
    INITIAL_STATE = 'START'
    
    def __init__(self):
        self.state_handler = interface.HandleState(self)
        for x in self.__state_definitions__:
            x.define_on(self.state_handler)
        self.state_handler.set_initial_state(self.INITIAL_STATE)

class StoppingCondition(object):
    
    def __init__(self, name):
        self.name = name
        self._is_enabled = False
        self._is_set = False
        self._particles = [datamodel.Particles() for x in range(256)]
        
    def is_supported(self):
        return True
        
    def is_enabled(self):
        return self._is_enabled
        
    def is_set(self):
        return self._is_enabled and self._is_set
    
    def enable(self):
        self._is_enabled = True
        
    def disable(self):
        self._is_enabled = False
    
    def set(self, *particles):
        self._is_set = True
        for i, particle_set in enumerate(particles):
            if len(self._particles) == i:
                self._particles.append(datamodel.Particles())
            
            self._particles[i].add_particles(particle_set)
    
    def unset(self):
        if self._is_set:
            self._is_set = False
            self._particles = [datamodel.Particles() for x in range(10)]
        

    def particles(self, index):
        if index >= len(self._particles):
            return datamodel.Particles()
            
        return self._particles[index]
        
    
    
