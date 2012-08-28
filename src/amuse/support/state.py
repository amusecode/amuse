from amuse.support import exceptions
from amuse.support.options import OptionalAttributes
from amuse.support.options import option
from amuse.support.thirdparty.texttable import Texttable

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
    
    def remove_from_transition(self, to_state):
        index = -1
        for i, transition in enumerate(self.from_transitions):
            if transition.to_state is to_state:
                index = i
        if index >= 0:
            del self.from_transitions[index] 
        
    def remove_to_transition(self, from_state):
        index = -1
        for i, transition in enumerate(self.to_transitions):
            if transition.from_state is from_state:
                index = i
        if index >= 0:
            del self.to_transitions[index] 
            
    def get_to_transitions(self):
        return list(self.to_transitions)
        
    def is_named(self):
        return True
        
class AllExcept(object):
    
    def __init__(self, states):
        self.states = states
        self.from_transitions = []
        self.to_transitions = []
    
    def is_named(self):
        return False
        
    def __str__(self):
        return "all except {0}".format(', '.join([str(x) for x in self.states]))
        
    def matches(self, other):
        for x in self.states:
            if other == x:
                return False
        return True
    
    def add_from_transition(self, transition):
        """add a transition starting at this state"""
        self.from_transitions.append(transition)
        
    def add_to_transition(self, transition):
        """add a transition to this state"""
        raise Exception('you cannot define a transition to any except one state')

    def remove_from_transition(self, to_state):
        index = -1
        for i, transition in enumerate(self.from_transitions):
            if transition.to_state is to_state:
                index = i
        if index >= 0:
            del self.from_transitions[index] 
        
    def remove_to_transition(self, from_state):
        pass
            

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
        self._initial_state = None

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
        
    def remove_transition(self, from_name, to_name):
        from_state = self.new_state(from_name)
        to_state   = self.new_state(to_name)
        
        if not from_state is None:
            from_state.remove_from_transition(to_state)
    
        if not to_state is None:
            to_state.remove_to_transition(from_state)
            
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
        self._initial_state = self._current_state 
    
    

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
    
    def to_plantuml_string(self):
        lines = []
        lines.append('@startuml')
        initial_state = self._initial_state
        lines.append('[*] --> {0}'.format(initial_state.name))
        statenames = list(sorted(self.states.keys()))
        merged_transitions = {}
        for name in statenames:
            state = self.states[name]
            transitions = state.get_to_transitions()
            for transition in transitions:
                if transition.from_state.is_named():
                    if not transition.method is None:
                        transitionname = '{0}+{1}'.format(
                            transition.from_state.name,
                            transition.to_state.name
                        )
                        if transitionname in merged_transitions:
                            merged_transitions[transitionname][2].add(transition.method.function_name)
                        else:
                            merged_transitions[transitionname] = [
                                transition.from_state.name,
                                transition.to_state.name,
                                set([transition.method.function_name])
                            ]
                                
                                
                        #lines.append('{0} --> {1} : {2}'.format(
                        #        transition.from_state.name,
                        #        transition.to_state.name,
                        #        transition.method.function_name
                        #    )
                        #)
                    else:
                        
                        lines.append('{0} -> {1}'.format(
                                transition.from_state.name,
                                transition.to_state.name
                            )
                        )
                else:
                     for x in self.iter_states():
                        if x == transition.from_state:
                            continue
                        if not transition.from_state.matches(x):
                            continue
                        if not transition.method is None:
                            lines.append('{0} --> {1} : {2}'.format(
                                    x.name,
                                    transition.to_state.name,
                                    transition.method.function_name
                                )
                            )
                        else:
                            lines.append('{0} -> {1}'.format(
                                    x.name,
                                    transition.to_state.name,
                                )
                            )
        for fromname, toname, methodnames in merged_transitions.values():
            lines.append('{0} --> {1} : {2}'.format(
                    fromname,
                    toname,
                    '\\n'.join(methodnames)
                )
            )
        lines.append('@enduml')
        return '\n'.join(lines)

    def to_table_string(self, ignore_states = [], split = True):
        lines = []
        ignore_states = set(ignore_states)
        initial_state = self._initial_state
        lines.append('Initial state: {0}'.format(initial_state.name))
        statenames = list(sorted(self.states.keys()))
        merged_transitions = {}
        for name in statenames:
            state = self.states[name]
            transitions = state.get_to_transitions()
            for transition in transitions:
                if transition.from_state.is_named():
                    if not transition.method is None:
                        functionname = transition.method.function_name 
                    else:
                        functionname = '*'
                        
                    transitionname = '{0}+{1}'.format(
                        transition.from_state.name,
                        transition.to_state.name
                    )
                    if transitionname in merged_transitions:
                        merged_transitions[transitionname][2].add(functionname)
                    else:
                        merged_transitions[transitionname] = [
                            transition.from_state.name,
                            transition.to_state.name,
                            set([functionname])
                        ]
                else:
                     for x in self.iter_states():
                        if x == transition.from_state:
                            continue
                        if not transition.from_state.matches(x):
                            continue
                            
                        if not transition.method is None:
                            functionname = transition.method.function_name 
                        else:
                            functionname = '*'
                        
                    
                        transitionname = '{0}+{1}'.format(
                            x.name,
                            transition.to_state.name
                        )
                        if transitionname in merged_transitions:
                            merged_transitions[transitionname][2].add(functionname)
                        else:
                            merged_transitions[transitionname] = [
                                x.name,
                                transition.to_state.name,
                                set([functionname])
                            ]
        selectedstates = [x for x in statenames if not x in ignore_states]
        tostates = [x for x in selectedstates if not x == initial_state.name]
        
        endstates = []
        for fromstate in selectedstates:
            found = False
            for tostate in tostates:
                transitionname = '{0}+{1}'.format(
                    fromstate,
                    tostate
                )
                if transitionname in merged_transitions:
                    found = True
                    break
            if not found:
                endstates.append(fromstate)
                
        if len(endstates) == 1:
            lines.append('End state: {0}'.format(endstates[0]))
        else:
            lines.append('End states: {0}'.format(', '.join(endstates)))
        state_to_distance_from_start = {}
        state_to_distance_from_start[initial_state.name] = 0
            
        stack = [[initial_state.name, 0]]
        while len(stack) > 0:
            current, distance = stack.pop()
            for tostate in selectedstates:
                transitionname = '{0}+{1}'.format(
                    current,
                    tostate
                )
                if transitionname in merged_transitions:
                    if  tostate not in state_to_distance_from_start:
                        stack.append([tostate, distance+1])
                        state_to_distance_from_start[tostate] = distance + 1
                        
        
        state_to_distance_from_end = {}
        stack = []
        for x in endstates:
            state_to_distance_from_end[x] = 0
            stack.append([x, 0])
            
        while len(stack) > 0:
            tostate, distance = stack.pop()
            for fromstate in selectedstates:
                transitionname = '{0}+{1}'.format(
                    fromstate,
                    tostate
                )
                if transitionname in merged_transitions:
                    if  fromstate not in state_to_distance_from_end:
                        stack.append([fromstate, distance+1])
                        state_to_distance_from_end[fromstate] = distance + 1
        
        
        fromstates = [x for x in selectedstates if not x in endstates]
        fromstates = list(sorted(fromstates, key = lambda x :  -(state_to_distance_from_end[x] * len(state_to_distance_from_start)) + state_to_distance_from_start[x]))
        tostates = list(sorted(tostates, key = lambda x :  -(state_to_distance_from_end[x] * len(state_to_distance_from_start)) + state_to_distance_from_start[x]))
            
        if split:
            table = Texttable(max_width = -1)
            header = ['          to\nfrom']
            header.extend(tostates)
            rows = []
            rows.append(header)
        
            for fromstate in fromstates:
                
                row = []
                row.append(fromstate)
                for tostate in tostates:
                    transitionname = '{0}+{1}'.format(
                        fromstate,
                        tostate
                    )
                    if transitionname in merged_transitions:
                        _, _, functionnames = merged_transitions[transitionname]
                        row.append('\n'.join(functionnames))
                    else:
                        row.append('-')
                rows.append(row)
                
            table.add_rows(rows)
            table._compute_cols_width()
            widths = table._width
            splittostates = []
            currentostates = []
            w0 = widths[0]
            currentwidth = w0
            for x,state in zip(widths[1:], tostates):
                currentwidth += x
                if currentwidth > 80:
                    splittostates.append(currentostates)
                    currentostates = [state]
                    currentwidth = w0 + x
                else:
                    currentostates.append(state)
            
            splittostates.append(currentostates)
        else:
            splittostates = [tostates]
        for tostates in splittostates:
            header = ['          to\nfrom']
            header.extend(tostates)
            rows = []
            rows.append(header)
        
            for fromstate in fromstates:
                
                row = []
                row.append(fromstate)
                for tostate in tostates:
                    transitionname = '{0}+{1}'.format(
                        fromstate,
                        tostate
                    )
                    if transitionname in merged_transitions:
                        _, _, functionnames = merged_transitions[transitionname]
                        row.append('\n'.join(functionnames))
                    else:
                        row.append('-')
                rows.append(row)
       
            table = Texttable(max_width = -1)
            align = ["l",]
            align.extend("l" * len(tostates))
            table.set_cols_align(align)
            table.add_rows(rows)
            lines.append(table.draw())

        return '\n'.join(lines)

    def get_name_of_current_state(self):
        return self._current_state.name
    
    
