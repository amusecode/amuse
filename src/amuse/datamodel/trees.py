
class BinaryTreesOnAParticleSet(object):

    def __init__(self, particles_set, name_of_firstchild_attribute, name_of_secondchild_attribute):
        self.particles_set = particles_set
        self.name_of_firstchild_attribute = name_of_firstchild_attribute
        self.name_of_secondchild_attribute = name_of_secondchild_attribute
        
        
    def iter_roots(self):
        return self.iter_binary_trees()
        
    def iter_binary_trees(self):
        binaries = self._binaries()
        if len(binaries) == 0:
            return
        binaries_children1 = self._get_inner_nodes(binaries, self.name_of_firstchild_attribute)
        binaries_children2 = self._get_inner_nodes(binaries, self.name_of_secondchild_attribute)
            
        roots = (binaries - (binaries_children1 + binaries_children2))
        
        for particle in roots:
            yield BinaryTreeOnParticle(particle, self.name_of_firstchild_attribute, self.name_of_secondchild_attribute)
            
    def particles_not_in_a_multiple(self):
        binaries = self._binaries()
        binaries_children1 = self._get_descendant_nodes(self.particles_set, self.name_of_firstchild_attribute)
        binaries_children2 = self._get_descendant_nodes(self.particles_set, self.name_of_secondchild_attribute)
        if len(binaries) == 0:
            return self.particles_set
        else:
            singles = (self.particles_set - (self.roots() + binaries_children1 + binaries_children2))
            return singles
    
    def roots(self):
        binaries = self._binaries()
        if len(binaries) == 0:
            return binaries
        binaries_children1 = self._get_inner_nodes(binaries, self.name_of_firstchild_attribute)
        binaries_children2 = self._get_inner_nodes(binaries, self.name_of_secondchild_attribute)
        return (binaries - (binaries_children1 + binaries_children2))
        
    def _binaries(self):
        return self.particles_set.select_array(lambda x : x != [None], [self.name_of_firstchild_attribute,])

    def _get_inner_nodes(self, set, name_of_attribute):
        descendants = self._get_descendant_nodes(set, name_of_attribute)
        return descendants.select_array(lambda x : x != [None], [name_of_attribute,])

    def _get_descendant_nodes(self, set, name_of_attribute):
        return getattr(set, name_of_attribute).as_set().compressed()
        
        
class BinaryTreeOnParticle(object):

    def __init__(self, particle, name_of_firstchild_attribute = "child1" , name_of_secondchild_attribute = "child2"):
        self.particle = particle
        self.name_of_firstchild_attribute = name_of_firstchild_attribute
        self.name_of_secondchild_attribute = name_of_secondchild_attribute
        
        
    def iter_descendants(self):
        stack = [self.particle]
        while len(stack) > 0:
            current = stack.pop()
            children = []
            child1 = getattr(current, self.name_of_firstchild_attribute)
            if not child1 is None:
                yield child1
                children.append(child1)
                
            child2 = getattr(current, self.name_of_secondchild_attribute)
            if not child2 is None:
                yield child2
                children.append(child2)
            
            stack.extend(reversed(children))
            
    def iter_leafs(self):
        stack = [self.particle]
        while len(stack) > 0:
            current = stack.pop()
        
            children = []
            child1 = getattr(current, self.name_of_firstchild_attribute)
            if not child1 is None:
                children.append(child1)
            
            child2 = getattr(current, self.name_of_secondchild_attribute)
            if not child2 is None:
                children.append(child2)
            
            stack.extend(reversed(children))
        
            if len(children) == 0:
                yield current


    def iter_inner_nodes(self):
        stack = [self.particle]
        while len(stack) > 0:
            current = stack.pop()
        
            children = []
            child1 = getattr(current, self.name_of_firstchild_attribute)
            if not child1 is None:
                children.append(child1)
            
            child2 = getattr(current, self.name_of_secondchild_attribute)
            if not child2 is None:
                children.append(child2)
            
            stack.extend(reversed(children))
        
            if len(children) > 0:
                yield current

    def get_inner_nodes_subset(self):    
        keys = [x.key for x in self.iter_inner_nodes()]
        return self.particle.particles_set._subset(keys)

    def get_descendants_subset(self):    
        keys = [x.key for x in self.iter_descendants()]
        return self.particle.particles_set._subset(keys)

    def get_leafs_subset(self):    
        keys = [x.key for x in self.iter_leafs()]
        return self.particle.particles_set._subset(keys)
    
    def copy(self):
        copy_of_set = self.get_tree_subset().copy()
        root = copy_of_set[0]
        return BinaryTreeOnParticle(
            root,
            name_of_firstchild_attribute = self.name_of_firstchild_attribute,
            name_of_secondchild_attribute = self.name_of_secondchild_attribute
        )
        
    def get_tree_subset(self):    
        keys = [x.key for x in iter(self)]
        return self.particle.particles_set._subset(keys)

    def iter_events(self):
        stack = [('start', self.particle)]
        while len(stack) > 0:
            event, current = stack.pop()
            yield event,current
            if event == 'end':
                continue
            stack.append( ('end', current, ) ) 
        
            children = []
            child1 = getattr(current, self.name_of_firstchild_attribute)
            if not child1 is None:
                children.append( ('start', child1, ) )
                
            child2 = getattr(current, self.name_of_secondchild_attribute)
            if not child2 is None:
                children.append( ('start', child2, ) )
            
            stack.extend(reversed(children))
            
    def iter_levels(self):
        level = -1
        
        for event, particle in self.iter_events():
            if event == 'start':
                level += 1
                yield level, particle
            else:
                level -= 1
    
    
    def __iter__(self):
        stack = [self.particle]
        while len(stack) > 0:
            current = stack.pop()
            yield current
            
            children = []
            child1 = getattr(current, self.name_of_firstchild_attribute)
            if not child1 is None:
                children.append(child1)
                
            child2 = getattr(current, self.name_of_secondchild_attribute)
            if not child2 is None:
                children.append(child2)
            
            stack.extend(reversed(children))
            

class ChildTreeOnParticleSet(object):

    def __init__(self, particles_set, names_of_child_attributes = ["child1" , "child2"]):
        self.particles_set = particles_set
        self.names_of_child_attributes = names_of_child_attributes
        
    @property
    def particle(self):
        return self.particle_set #a particle set associated, maybe upgrade to some aggregate particle
        
            
    def is_leaf(self):
        return False
            
    def iter_leafs(self):
        for node in self.iter_children():
            if node.is_leaf():
                yield node
            
    def iter_branches(self):
        for node in self.iter_children():
            if not node.is_leaf():
                yield node
                
            
    def iter_descendants(self):
        stack = list(reversed(list(self.iter_children())))
        while len(stack) > 0:
            current = stack.pop()
            
            yield current
            
            children = list(current.iter_children())
            
            stack.extend(reversed(children))
            
    
    def iter_events(self):
        stack = [('start', x) for x in (reversed(list(self.iter_children())))]
        while len(stack) > 0:
            event, current = stack.pop()
            yield event,current
            if event == 'end':
                continue
            stack.append( ('end', current, ) ) 
            
            children = list(current.iter_children())
            
            stack.extend([('start', x) for x in (reversed(list(self.iter_children())))])
            
            
    def iter_descendant_leafs(self):
        stack = list(reversed(list(self.iter_children())))
        while len(stack) > 0:
            current = stack.pop()

            if current.is_leaf():
                yield current
            else:
                children = list(current.iter_children())
                stack.extend(reversed(children))


    def iter_descendant_branches(self):
        
        stack = list(reversed(list(self.iter_children())))
        while len(stack) > 0:
            current = stack.pop()
            
            if not current.is_leaf():
                yield current
                children = list(current.iter_children())
                stack.extend(reversed(children))

    def iter_children(self):
        for particle in self._children():
            yield ChildTreeOnParticle(particle, self.names_of_child_attributes)
            
        
    def iter_levels(self):
        level = -1
        
        for event, particle in self.iter_events():
            if event == 'start':
                level += 1
                yield level, particle
            else:
                level -= 1
    
    
    
    def get_children(self):
        return list(self.iter_children())
        
    def __iter__(self):
        return self.iter_children()

    def _branches(self):
        binaries = self._binaries()
        result = binaries
        
        for name in self.names_of_child_attributes:
            result -= self._get_inner_nodes(binaries, name)
        
        return result
    
    def _inner_particles(self):
        result = None
        
        for name in self.names_of_child_attributes:
            descendants = self._get_descendant_nodes(self.particles_set, name)
            if result is None:
                result = descendants
            else:
                result += descendants
                
        return result
        
    def _children(self):
        return self.particles_set - self._inner_particles()
        
    def _get_descendant_nodes(self, set, name_of_attribute):
        return getattr(set, name_of_attribute).as_set().compressed()
        


class ChildTreeOnParticle(object):

    def __init__(self, particle, names_of_child_attributes = ["child1" , "child2"]):
        self.particle = particle
        self.names_of_child_attributes = names_of_child_attributes
        
        
    def iter_descendants(self):
        stack = [self.particle]
        while len(stack) > 0:
            current = stack.pop()
            children = []
            for name in self.names_of_child_attributes:
                child  = getattr(current, name)
                if not child is None:
                    yield ChildTreeOnParticle(child1, self.names_of_child_attributes)
                    children.append(child)
            
            stack.extend(reversed(children))
            
    def is_leaf(self):
        for name in self.names_of_child_attributes:
            child  = getattr(self.particle, name)
            if not child is None:
                return False
        return True
            
    def iter_leafs(self):
        for node in self.iter_children():
            if node.is_leaf():
                yield node
                
    def iter_branches(self):
        for node in self.iter_children():
            if not node.is_leaf():
                yield node
            
    def iter_descendant_leafs(self):
        stack = list(reversed(list(self.iter_children())))
        while len(stack) > 0:
            current = stack.pop()

            if current.is_leaf():
                yield current
            else:
                children = list(current.iter_children())
                stack.extend(reversed(children))

    def iter_descendant_branches(self):
        
        stack = list(reversed(list(self.iter_children())))
        while len(stack) > 0:
            current = stack.pop()
            
            if not current.is_leaf():
                yield current
                children = list(current.iter_children())
                stack.extend(reversed(children))

    def iter_children(self):
        current = self.particle
        
        for name in self.names_of_child_attributes:
            child  = getattr(current, name)
            if not child is None:
                yield ChildTreeOnParticle(child, self.names_of_child_attributes)
    
    def get_children(self):
        return list(self.iter_children())
        
    def __iter__(self):
        return self.iter_children()
        
    def iter_levels(self):
        level = -1
        
        for event, particle in self.iter_events():
            if event == 'start':
                level += 1
                yield level, particle
            else:
                level -= 1
    
    def iter_events(self):
        stack = [('start', self,), ]
        while len(stack) > 0:
            event, current = stack.pop()
            yield event,current
            if event == 'end':
                continue
            stack.append( ('end', current, ) ) 
            
            children = list(current.iter_children())
            
            stack.extend([('start', x) for x in reversed(children)])
            
      
