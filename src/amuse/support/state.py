from amuse.support import exceptions
from amuse.support.options import OptionalAttributes
from amuse.support.options import option

import logging

class State(object):
    def __init__(self, handler, name):
        self.handler = handler
        self.name = name
        self.from_transitions = []
        self.to_transitions = []

    def __str__(self):
        return "state '{0}'".format(self.name)
       
    def matches(self, other):
        return other == self
    
    def add_from_transition(self, transition):
        """add a transition starting at this state"""
        self.from_transitions.append(transition)
        
    def add_to_transition(self, transition):
        """add a transition to this state"""
        self.to_transitions.append(transition)
        
        
    def get_to_transitions(self):
        return list(self.to_transitions)
        
class AllExcept(object):
    
    def __init__(self, states):
        self.states = states
        self.from_transitions = []
        self.to_transitions = []
    
    def __str__(self):
        return "all except {0}".format(', '.join([str(x) for x in self.states]))
        
    def matches(self, other):
        for x in self.states:
            if other == x:
                return False
        return True
    
    def add_from_transition(self, transition):
        """add a transition starting at this state"""
        pass # ignored
        
    def add_to_transition(self, transition):
        """add a transition to this state"""
        raise Exception('you cannot define a transition to any except one state')


    def get_to_transitions(self):
        """to transitions are the to transitions of all states except this one"""
        result = []
        state_machine = self.states[0].handler
        for state in state_machine.iter_states():
            if state == self:
                continue
            if not self.matches(state):
                continue
            else:
                for transition in state.get_to_transitions():
                    if not transition.from_state == self:
                        result.append(transition)
        return result
        
        
class StateTransition(object):
    def __init__(self, handler, from_state, to_state, method = None, is_auto = True):
        self.method = method
        self.to_state = to_state
        self.from_state = from_state
        self.is_auto = is_auto
        self.handler = handler

    def __str__(self):
        return "transition from {0} to {1}".format(self.from_state, self.to_state)

    def do(self):
        if self.handler.log_transitions:
            logging.getLogger("state").info(str(self))
            
        if self.method is None:
            self.handler.current_state = self.to_state
        else:
            self.method.new_method()()


class StateMachine(OptionalAttributes):
    
    def __init__(self, interface, **options):
        OptionalAttributes.__init__(self, **options)
        
        self.states = {}
        
        self._do_automatic_state_transitions = True
        self._current_state = State(self, None)
        self.interface = interface

    @option(type='boolean', sections=['state',])
    def is_enabled(self):
        return True
        
    @option(type='boolean', sections=['state',])
    def log_transitions(self):
        return False
        
    def enable(self):
        self.is_enabled = True

    def disable(self):
        self.is_enabled = False


    def new_transition(self, from_name, to_name, is_auto = True):
        from_state = self.new_state(from_name)
        to_state   = self.new_state(to_name)
    
        transition = StateTransition(self, from_state, to_state, None, is_auto)
    
        if not from_state is None:
            from_state.add_from_transition(transition)
        
        if not to_state is None:
            to_state.add_to_transition(transition)
    
        return transition
    
    def iter_states(self):
        return iter(self.states.values())

    def new_state(self, name):
        if name is None:
            return None
        
        if name.startswith('!'):
            return AllExcept([self.new_state(x) for x in name[1:].split('!')])
            
        if name in self.states:
            return self.states[name]
            
        self.states[name] = State(self, name)
        return self.states[name]
    
    

    def set_initial_state(self, name):
        self._current_state = self.new_state(name)
    
    

    def _get_transitions_path_from_to(self, from_state, to_state):
        transitions = filter(lambda x : x.is_auto, to_state.get_to_transitions())
    
        paths = map(lambda x : [x], transitions)
    
        def has_no_circle(path):
            seen_states = set([])
            seen_states.add(path[0].from_state)
            for transition in path:
                if transition.to_state in seen_states:
                    return False
                seen_states.add(transition.to_state)
            return True
    
    
        while paths:
            current = paths.pop()
            first = current[0]
            if first.from_state is None:
                yield current
            elif first.from_state.matches(from_state):
                yield current
            else:
                transitions = filter(lambda x : x.is_auto, first.from_state.get_to_transitions())
                new_paths = map(lambda x : [x], transitions)
    
                for new_path in new_paths:
                    new_path.extend(current)
                new_paths = filter(has_no_circle, new_paths)
    
                paths.extend(new_paths)
    
    
        return
    
    

    def _do_state_transition_to(self, state):
        all_transitions = list(self._get_transitions_path_from_to(self._current_state, state))
        transitions = []
        for x in all_transitions:
            if len(transitions) == 0 or len(x) < len(transitions):
                transitions = x
                
        if len(transitions) == 0:
            raise Exception("No transition from current state {0} to {1} possible".format(self._current_state, state))
        
        transitions_with_methods = filter(lambda x : not x.method is None,transitions)
        if not self._do_automatic_state_transitions and len(transitions_with_methods) > 0:
            lines = []
            lines.append("Interface is not in {0}, should transition from {1} to {0} first.\n". format(state, self._current_state))
            for x in transitions:
                if x.method is None:
                    lines.append("{0}, automatic". format(x))
                else:
                    lines.append("{0}, calling '{1}'". format(x, x.method.function_name))
            exception = exceptions.AmuseException('\n'.join(lines))
            exception.transitions = transitions
            raise exception
    
        for transition in transitions:
            transition.do()
    
    

    def get_name_of_current_state(self):
        return self._current_state.name
    
    
