from amuse.support import exceptions


class State(object):
    def __init__(self, handler, name):
        self.handler = handler
        self.name = name
        self.from_transitions = []
        self.to_transitions = []

    def __str__(self):
        return "state '{0}'".format(self.name)


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
        if self.method is None:
            self.handler.current_state = self.to_state
        else:
            self.method.new_method()()


class StateMachine(object):
    def __init__(self, interface):
        self.states = {}
        
        self._do_automatic_state_transitions = True
        self._current_state = None
        self.interface = interface
        self._is_enabled = True

    def is_enabled(self):
        return self._is_enabled
        
    def enable(self):
        self._is_enabled = True

    def disable(self):
        self._is_enabled = False


    def new_transition(self, from_name, to_name, is_auto = True):
        from_state = self.new_state(from_name)
        to_state   = self.new_state(to_name)
    
        transition = StateTransition(self, from_state, to_state, None, is_auto)
    
        if from_state:
            from_state.from_transitions.append(transition)
        
        if to_state:
            to_state.to_transitions.append(transition)
    
        return transition
    
    

    def new_state(self, name):
        if name is None:
            return None
            
        if name in self.states:
            return self.states[name]
            
        self.states[name] = State(self, name)
        return self.states[name]
    
    

    def set_initial_state(self, name):
        print name
        self._current_state = self.new_state(name)
    
    

    def _get_transitions_path_from_to(self, from_state, to_state):
        transitions = filter(lambda x : x.is_auto, to_state.to_transitions)
    
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
            elif first.from_state is None:
                return current
            else:
                transitions = filter(lambda x : x.is_auto, first.from_state.to_transitions)
                new_paths = map(lambda x : [x], transitions)
    
                for new_path in new_paths:
                    new_path.extend(current)
    
                new_paths = filter(has_no_circle, new_paths)
    
                paths.extend(new_paths)
    
    
        return None
    
    

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
                    lines.append("{0}, calling '{1}'". format(x, x.method.function_name))
            exception = exceptions.AmuseException('\n'.join(lines))
            exception.transitions = transitions
            raise exception
    
        for transition in transitions:
            transition.do()
    
    

    def get_name_of_current_state(self):
        return self._current_state.name
    
    
