
class BinaryTreesOnAParticleSet(object):

    def __init__(self, particles_set, name_of_firstchild_attribute, name_of_secondchild_attribute):
        self.particles_set = particles_set
        self.name_of_firstchild_attribute = name_of_firstchild_attribute
        self.name_of_secondchild_attribute = name_of_secondchild_attribute
        
        
    def iter_roots(self):
        return self.iter_binary_trees()
        
    def iter_binary_trees(self):
        binaries = self._binaries()
        binaries_children1 = self._get_inner_nodes(binaries, self.name_of_firstchild_attribute)
        binaries_children2 = self._get_inner_nodes(binaries, self.name_of_secondchild_attribute)
    
        roots = (binaries - (binaries_children1 + binaries_children2))
        
        for particle in roots:
            yield BinaryTreeOnParticle(particle, self.name_of_firstchild_attribute, self.name_of_secondchild_attribute)
            
    def particles_not_in_a_multiple(self):
        binaries = self._binaries()
        binaries_children1 = self._get_descendant_nodes(self.particles_set, self.name_of_firstchild_attribute)
        binaries_children2 = self._get_descendant_nodes(self.particles_set, self.name_of_secondchild_attribute)
    
        singles = (self.particles_set - (self.roots() + binaries_children1 + binaries_children2))
        return singles
    
    def roots(self):
        binaries = self._binaries()
        binaries_children1 = self._get_inner_nodes(binaries, self.name_of_firstchild_attribute)
        binaries_children2 = self._get_inner_nodes(binaries, self.name_of_secondchild_attribute)
        
        return (binaries - (binaries_children1 + binaries_children2))
        
    def _binaries(self):
        return self.particles_set.select_array(lambda x : x.get_valid_particles_mask(), [self.name_of_firstchild_attribute,])

    def _get_inner_nodes(self, set, name_of_attribute):
        descendants = self._get_descendant_nodes(set, name_of_attribute)
        return descendants.select_array(lambda x : x.get_valid_particles_mask(), [name_of_attribute,])

    def _get_descendant_nodes(self, set, name_of_attribute):
        return getattr(set, name_of_attribute).compress()
        
        
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
        copy_of_set = self.get_tree_subset().copy_to_memory()
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
            

