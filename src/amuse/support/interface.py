
from amuse.support.data import parameters
from amuse.support.data import binding
from amuse.support.data import core
from amuse.support.data import values
from amuse.support.units import nbody_system


import inspect

class OldObjectsBindingMixin(object):
                
    def setup_particles(self, particles):
        self.particles.add_particles(particles)
        
    def update_particles(self, particles):
        self.particles.copy_values_of_state_attributes_to(particles)
    
    def set_attribute(self, attribute_name, particles):
        particles.copy_values_of_attribute_to(attribute_name, self.particles)
        
    def update_attribute(self, attribute_name, particles):
        self.particles.copy_values_of_attribute_to(attribute_name, particles) 
        
        
        
class CodeInterface(OldObjectsBindingMixin):
    """Base class of next level interfaces to codes
    """
    
    parameter_definitions = []
    
    def __init__(self):
        self.parameters = parameters.Parameters(self.parameter_definitions, self)
        

class CodeInterfaceWithConvertedUnits(OldObjectsBindingMixin):
    class UnitsConvertionMethod(object):
        def __init__(self, real_method, converter):
            self.real_method = real_method
            self.converter = converter
            
        def __call__(self, *list_arguments, **keyword_arguments):
            converted_list_arguments = [self.from_source_to_target(x) for x in list_arguments]
            converted_keyword_arguments = {}
            for key, value in keyword_arguments:
                converted_keyword_arguments[key] = self.from_source_to_target(value)
            
            result = self.real_method(*converted_list_arguments, **converted_keyword_arguments)
            print result
            return self.from_target_to_source(result)
            
        def from_source_to_target(self, x):
            if isinstance(x, values.Quantity):
                return self.converter.from_source_to_target(x)
            else:
                return x
                
        def from_target_to_source(self, x):
            if isinstance(x, values.Quantity):
                return self.converter.from_target_to_source(x)
            elif hasattr(x, '__iter__'):
                return self.convert_and_iterate(x)
            else:
                return x
        
        def convert_and_iterate(self, iterable):
            for x in iterable:
                yield self.converter.from_target_to_source(x)
            
    def __init__(self, interface, converter):
        self.converter = converter
        self.interface = interface
        
    def __getattr__(self, name):
        attribute = getattr(self.interface, name)
        if inspect.ismethod(attribute):
            result = self.UnitsConvertionMethod(attribute, self.converter)
            setattr(self, name, result)
            return result
        elif isinstance(attribute, core.AbstractParticleSet):
            result = core.ParticlesWithUnitsConverted(attribute, self.converter)
            setattr(self, name, result)
            return result
        elif isinstance(attribute, values.Quantity):
            return self.converter.from_target_to_source(attribute)
        elif isinstance(attribute, binding.BoundCodeMethod):
            result = self.UnitsConvertionMethod(attribute, self.converter)
            setattr(self, name, result)
            return result
        elif isinstance(attribute, parameters.Parameters):
            result = parameters.ParametersWithUnitsConverted(attribute, self.converter)
            setattr(self, name, result)
            return result
        elif hasattr(attribute, '__iter__'):
            return self.convert_and_iterate(attribute)
        else:
            return attribute
            
    def convert_and_iterate(self, iterable):
        for x in iterable:
            yield self.converter.from_target_to_source(x)
            
    def __dir__(self):
        return dir(self.interface)
        
class CodeInterfaceWithNBodyUnitsConverted(CodeInterfaceWithConvertedUnits):
    
    def __init__(self, interface, convert_nbody):
        if convert_nbody is None:
            convert_nbody = nbody_system.nbody_to_si.get_default()
        
        CodeInterfaceWithConvertedUnits.__init__(
            self,
            interface,
            convert_nbody.as_converter_from_si_to_nbody()
        )
        
        
        


class CodeInterfaceWithStateEngine(OldObjectsBindingMixin):
    class State(object):
        def __init__(self, interface, name):
            self.interface = interface
            self.name = name
            self.from_transitions = []
            self.to_transitions = []
            
        def __str__(self):
            return "state '{0}'".format(self.name)
            
    class StateTransition(object):
        def __init__(self, interface, from_state, to_state, method = None, is_default = False):
            self.method = method
            self.to_state = to_state
            self.from_state = from_state
            self.is_default = is_default
            self.interface = interface
            
        def __str__(self):
            return "transition from {0} to {1}".format(self.from_state, self.to_state)
            
        def do(self):
            if self.method is None:
                self.interface.current_state = self.to_state
            else:
                self.method()
                
    class StateMethod(object):
        def __init__(self, interface, from_state, to_state, function):
            self.interface = interface
            self.transitions = []
            self.add_transition(from_state, to_state)
            self.function = function
        
        def add_transition(self, from_state, to_state):
            self.transitions.append((from_state, to_state))
        
        def precall(self):
            stored_transitions = []
            
            for from_state, to_state in self.transitions:
                if from_state is None:
                    return to_state
                elif from_state == self.interface._current_state:
                    return to_state
                else:
                    stored_transitions.append((from_state, to_state))
            
            for from_state, to_state  in stored_transitions:
                try:
                    self.interface._do_state_transition_to(from_state)
                    return to_state
                except Exception, ex:
                    pass
            
            # do again to get an exception.
            self.interface._do_state_transition_to(stored_transitions[0][0])
                
        def postcall(self, to_state):
            if to_state is None:
                return
            elif to_state == self.interface._current_state:
                return
            else:
                self.interface._current_state = to_state

        def __call__(self, *list_arguments, **keyword_arguments):
            to_state = self.precall()
            real_method = getattr(self.interface.interface, self.function)
            result = real_method(*list_arguments, **keyword_arguments)
            self.postcall(to_state)
            return result
            
        
            
    def __init__(self, interface):
        self.interface = interface
        self._mapping_from_name_to_state_method = {}
        self._current_state = None
        self._do_automatic_state_transitions = True
        self.states = {}
        self.setup_states()
        
        
    def setup_states(self):
        """
        Abstract method, needs to be overridded by subclasses to
        actually define the states.
        """
        pass
        
    def __getattr__(self, name):
        if name in self._mapping_from_name_to_state_method:
            return self._handle_get_attribute(name)
        else:
            return getattr(self.interface, name)
        
    def define_state(self, name):
        if name in self.states:
            return
        
        self.states[name] = self.State(self, name)
    
    
    def _handle_get_attribute(self, name):
        state_method = self._mapping_from_name_to_state_method[name]
        return state_method 
         
    def _do_state_transition_to(self, state):
        transitions = self._get_transitions_path_from_to(self._current_state, state)
        if transitions is None:
            raise Exception("No transition from current state {0} to {1} possible".format(self._current_state, state))
        
        transitions_with_methods = filter(lambda x : not x.method is None,transitions)
        if not self._do_automatic_state_transitions and len(transitions_with_methods) > 0:
            lines = []
            lines.append("Interface is not in {0}, should transition from {1} to {0} first.\n". format(state, self._current_state))
            for x in transitions:
                if x.method is None:
                    lines.append("{0}, automatic". format(x))
                else:
                    lines.append("{0}, calling '{1}'". format(x, x.method.function))
            exception = Exception('\n'.join(lines))
            exception.transitions = transitions
            raise exception
            
        for transition in transitions:
            transition.do()
        
    
    def _get_transitions_path_from_to(self, from_state, to_state):
        transitions = to_state.to_transitions
        
        paths = map(lambda x : [x], transitions)
        
        def has_no_circle(path):
            seen_states = set([])
            for transition in path:
                if transition.to_state in seen_states:
                    return False
                seen_states.add(transition.to_state)
            return True
                    
            
        while paths:
            current = paths.pop()
            first = current[0]
            if first.from_state == from_state:
                return current
            else:
                transitions = first.from_state.to_transitions
                new_paths = map(lambda x : [x], transitions)

                for new_path in new_paths:
                    new_path.extend(current)
            
                new_paths = filter(has_no_circle, new_paths)
            
                paths.extend(new_paths)
                
        return None
                
        
    def add_transition(self, from_name, to_name, function_name):
        self.define_state(from_name)
        self.define_state(to_name)
        
        from_state = self.states[from_name]
        to_state   = self.states[to_name]
        method = self.StateMethod(self, from_state, to_state, function_name)
        transition = self.StateTransition(self, from_state, to_state, method)
        
        from_state.from_transitions.append(transition)
        to_state.to_transitions.append(transition)
        
        self._add_state_method(from_state, to_state, function_name)
        
    def add_method(self, state_name, function_name):
        """
        Define a method that can run when the interface is in the
        provided state.
        """
        self.define_state(state_name)
        
        state = self.states[state_name]
        
        self._add_state_method( state, None, function_name)
        
    def _add_state_method(self, from_state, to_state, function_name):
        if not function_name in self._mapping_from_name_to_state_method:
            state_method = self.StateMethod(self, from_state, to_state, function_name)
            self._mapping_from_name_to_state_method[function_name] = state_method
        else:
            state_method = self._mapping_from_name_to_state_method[function_name]
            state_method.add_transition(from_state, to_state)
            
        
        
        
    def add_transition_to_method(self, state_name, function_name):
        """
        Define a method that can run in any state and will transition the interface
        to the provided state.
        """
        self.define_state(state_name)
        
        state = self.states[state_name]
        
        self._add_state_method(None, state, function_name)

    def set_initial_state(self, name):
        self.define_state(name)
        self._current_state = self.states[name]
        
    def do_automatic_state_transitions(self, boolean):
        self._do_automatic_state_transitions = boolean
