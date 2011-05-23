from amuse.support.units import units
from amuse.support.data import core
from amuse.support.data import values
from amuse.support import interface
from amuse.support import exceptions

from amuse.test import amusetest

class ExampleParticlesInterface(interface.InCodeComponentImplementation):
    """This is an example class to demonstrate how to work with incode particle sets
    using the object mapping interface.
    
    The particle set mapping has two functions:
    
    1) Convert function calls to attribute access. Setting or getting
       an attribute is converted to a call to the code.
          
          particles[0].mass = 10 | units.kg
          -->
          code.set_mass(index_of_the_particle = 0, mass = 10)
    
    2) Convert ids from the code into particle keys. All particles
       in AMUSE have an unique key, this key is mapped to the corresponding
       id of that particle in a code by the particle set. As every particle
       set mapping is unique for a code, every code has a mapping table (and
       different codes can have differen ids for the same particle, the key
       to id mapping will ensure that attribute values are assigned correctly)
       
    This example class is not connected to a real code, instead it 
    provides all necessary methods to create, delete and update particles. 
    
    For simplicity all attributes are stored in quantities (values with units),
    no unit conversion is provided.
    
    In this example code particles have an id (id),  position (x,y and z) and a
    mass attribute (mass).
    
    See the ExampleParticlesInterfaceTests class for examples / tests of use.
    """
    
    def log(self, message, *arguments):
        print "IN CODE >>", message.format(*arguments)
    
    def __init__(self):
        super(type(self), self).__init__(None)
        self.mapping_from_id_to_particle = {}
        self.highest_id = 0
        self.log("initialized the code")


    def define_particle_sets(self, builder):
        """In this method we define the particle set and all attribute.
        We specify wich methods to call for creating and deleting particles.
        We specify the methods to access attributes of the particles
        We specify the methods to create subsets or links
        
        Redefined from CodeInterface.
        """
        
        
        # first we define the name of the set and wich parameter will be used
        # to specify the id of a particle in the code
        builder.define_set('particles', 'index_of_the_particle')
        
        # we define new and delete methods
        builder.set_new('particles', 'new_particle')
        builder.set_delete('particles', 'delete_particle')
        
        # we define the methods to get and set
        # attributes, note these attributes may overlap
        # we need to specify the returned attribute names for the
        # get functions as the framework cannot derive these from
        # python methods (python methods do not have named return values
        builder.add_setter('particles', 'set_state')
        builder.add_getter('particles', 'get_state', names = ('mass', 'x', 'y', 'z'))
        builder.add_setter('particles', 'set_mass')
        builder.add_getter('particles', 'get_mass', names = ('mass',))
        builder.add_setter('particles', 'set_position')
        builder.add_getter('particles', 'get_position', names = ('x', 'y', 'z'))
        
        builder.add_method('particles',  'get_list', 'element_list' )
        
        self.log("defined the particle set with name {0!r}", 'particles')

    def new_particle(self, mass, x, y, z):
        """Creates new particles.
        
        Note: the parameternames are directly coupled to the attribute names
        in the particle set
        
        Note: the arguments are all arrays
        """
        
        result = []
        for mass_element, x_element, y_element, z_element in zip(mass, x, y, z):
            particle = [self.highest_id, mass_element, x_element, y_element, z_element]
            self.mapping_from_id_to_particle[self.highest_id] = particle
            result.append(self.highest_id)
            
            self.log("created new particle with id {0}, (mass={1}, x={2}, y={3}, z={4})", self.highest_id, mass_element, x_element, y_element, z_element)
            
            self.highest_id += 1
        
        return result

    def delete_particle(self, index_of_the_particle):
        """Delete particles in array index_of_the_particle
        
        Note: the parametername index_of_the_particle was specified in the "define_particle_sets" method
        """
        
        for x in index_of_the_particle:
            del self.mapping_from_id_to_particle[x]
            
            self.log("deleted the particle with id {0}", x)

    def set_state(self, index_of_the_particle, mass, x, y, z):
        """Sets the mass and the position of a particle.
        
        Note: the arguments are arrays
        """
        
        for index_element, mass_element, x_element, y_element, z_element in zip(index_of_the_particle, mass, x, y, z):
            particle = self.mapping_from_id_to_particle[index_element]
            particle[1] = mass
            particle[2] = x
            particle[3] = y
            particle[4] = z
            
            
            self.log("updated state of particle with id {0}", index_element)

    def get_state(self, index_of_the_particle):
        """Returns arrays for the mass, x, y and z values
        """
        
        massresult = values.AdaptingVectorQuantity()
        xresult = values.AdaptingVectorQuantity()
        yresult = values.AdaptingVectorQuantity()
        zresult = values.AdaptingVectorQuantity()
        
        for index_element in index_of_the_particle:
            particle = self.mapping_from_id_to_particle[index_element]
            
            massresult.append(particle[1])
            xresult.append(particle[2])
            yresult.append(particle[3])
            zresult.append(particle[4])
            
            self.log("retrieved state of particle with id {0}", index_element)
            
        return massresult, xresult, yresult, zresult

    def get_mass(self, index_of_the_particle):
        """Returns an array for the masses of the indices int the index_of_the_particle array
        """
        
        massresult = values.AdaptingVectorQuantity()
        
        for index_element in index_of_the_particle:
            particle = self.mapping_from_id_to_particle[index_element]
            massresult.append(particle[1])
            
            self.log("retrieved mass of particle with id {0} (mass = {1})", index_element, particle[1])
            
        return massresult
    

    def set_mass(self, index_of_the_particle, mass):
        """Sets the mass and the position of a particle.
        
        Note: the arguments are arrays
        """
        
        for index_element, mass_element in zip(index_of_the_particle, mass):
            particle = self.mapping_from_id_to_particle[index_element]
            
            particle[1] = mass
            
            self.log("updated mass of particle with id {0} (mass = {1})", index_element, particle[1])
            

    def get_position(self, index_of_the_particle):
        """Returns an array of the positions for the indices in index_of_the_particle
        """
        
        xresult = values.AdaptingVectorQuantity()
        yresult = values.AdaptingVectorQuantity()
        zresult = values.AdaptingVectorQuantity()
        
        for index_element in index_of_the_particle:
            particle = self.mapping_from_id_to_particle[index_element]
            
            xresult.append(particle[2])
            yresult.append(particle[3])
            zresult.append(particle[4])
            
            self.log("retrieved position of particle with id {0} (x = {1}, y = {2}, z = {3})", index_element, particle[2], particle[3], particle[4])
        return xresult, yresult, zresult

    def set_position(self, index_of_the_particle, x, y, z):
        """Sets the mass and the position of a particle.
        
        Note: the arguments are arrays
        """
        
        for index_element, x_element, y_element, z_element in zip(index_of_the_particle, x, y, z):
            particle = self.mapping_from_id_to_particle[index_element]
            
            particle[2] = x
            particle[3] = y
            particle[4] = z
            
            self.log("updated position of particle with id {0} (x = {1}, y = {2}, z = {3})", index_element, particle[2], particle[3], particle[4])
    
    def get_list_size(self, index_of_the_particle):
        """Returns the inclusive range of indices in the
        list of element coupled to a particle.
        """
        return (0,9)
    
    def get_list_element(self, index_in_the_list, index_of_the_particle ):
        """Returns an array of the positions for the indices in index_of_the_particle
        """
        if not hasattr(index_in_the_list, '__iter__'):
            index_in_the_list = [index_in_the_list,]
        if not hasattr(index_of_the_particle, '__iter__'):
            index_of_the_particle = [index_of_the_particle,]
        value1 = values.AdaptingVectorQuantity()
        value2 = values.AdaptingVectorQuantity()
        
        for index_of_one_particle, index_of_one_element in zip(index_of_the_particle, index_in_the_list):
            
            value1.append(index_of_one_particle | units.none)
            value2.append(index_of_one_element | units.none)
            
        return value1, value2
    
    def get_list(self, index_of_the_particle):
        if hasattr(index_of_the_particle, '__iter__'):
            return [self._create_new_grid(self.specify_list, index_of_the_particle = x) for x in index_of_the_particle]
        else:
            return self._create_new_grid(self.specify_list, index_of_the_particle = index_of_the_particle)
    
    def specify_list(self, definition, index_of_the_particle = 0):
        definition.set_grid_range('get_list_size')
        definition.add_getter('get_list_element', names=('value1', 'value2'))
        definition.define_extra_keywords({'index_of_the_particle':index_of_the_particle})
    

class ExampleParticlesInterfaceTests(amusetest.TestCase):
    """This class runs tests on the example particles interface
    class.
    
    """
    def test1(self):
        """
        In this test we will add and remove a particle from
        the particle set of the code.
        
        Adding a particle to the set will result in creating
        a new particle in the code.
        
        Removing a particle form the set will result in deleting
        the corresponding particle in the code.
        """
        
        self.log("adding and removing of a particle")
        
        instance = ExampleParticlesInterface()
        self.assertEquals(len(instance.particles), 0)
    
        # we create a particle in our script
        # all attributes of this particle are stored in the python space
        # when creating a particle you can set it's key
        # or let the system determine a unique key
        # to set the key in the script do: core.Particle(1000) where 1000 is the key
        
        theParticle = core.Particle()
        theParticle.mass = 10 | units.kg
        theParticle.x = 0.1 | units.m
        theParticle.y = 0.2 | units.m
        theParticle.z = 0.5 | units.m
        
        self.log("Adding particle with key {0}", theParticle.key)
        instance.particles.add_particle(theParticle)
        
        print instance.particles.index_in_code
        self.assertEquals(len(instance.particles), 1)
        
        self.log("Removing particle with key {0}", theParticle.key)
        instance.particles.remove_particle(theParticle)
        
        self.assertEquals(len(instance.particles), 0)
        

    def test2(self):
        """
        In this test we will set and get different properties
        of the particle.
        
        To limit overhead, the system will use the set_* or get_* calls 
        that are the closests match to the attributes queries.
        """
        self.log("accessing attributes of a particle")
        instance = ExampleParticlesInterface()
        self.assertEquals(len(instance.particles), 0)
        
        theParticle = core.Particle()
        theParticle.mass = 10 | units.kg
        theParticle.x = 0.1 | units.m
        theParticle.y = 0.2 | units.m
        theParticle.z = 0.5 | units.m
        
        instance.particles.add_particle(theParticle)
        
        self.log("Getting the mass of particle with key {0}, get_mass should be called", theParticle.key)
        self.assertEquals(instance.particles[0].mass, 10 | units.kg)
        
        self.log("Getting the position of particle with key {0}, get_position should be called", theParticle.key)
        self.assertEquals(instance.particles[0].position, [0.1, 0.2, 0.5] | units.m)
        
        self.log("Getting the only the x attribute of particle with key {0}, get_position should be called (y and z are discarded)", theParticle.key)
        self.assertEquals(instance.particles[0].x, 0.1 | units.m)
        
        
        self.log("Setting the position of particle with key {0}, set_position should be called", theParticle.key)
        instance.particles[0].position =  [0.2, 0.3, 0.6] | units.m
        
        self.log("Setting the x of particle with key {0}, should fail as no function can set x and no others", theParticle.key)
        def block():
            instance.particles[0].x =  0.1 | units.m
        self.assertRaises(Exception, block)
        
        
        
    def test3(self):
        """
        In this test we will get a list from a particle
        """
        
        instance = ExampleParticlesInterface()
        self.assertEquals(len(instance.particles), 0)
        
        theParticle = core.Particle()
        theParticle.mass = 10 | units.kg
        theParticle.x = 0.1 | units.m
        theParticle.y = 0.2 | units.m
        theParticle.z = 0.5 | units.m
        
        instance.particles.add_particle(theParticle)
        
        theParticle = core.Particle()
        theParticle.mass = 11 | units.kg
        theParticle.x = 0.1 | units.m
        theParticle.y = 0.2 | units.m
        theParticle.z = 0.5 | units.m
        
        instance.particles.add_particle(theParticle)
        
        self.assertEquals(len(instance.particles[0].element_list()), 10)
        list = instance.particles[0].element_list()
        self.assertEquals(list[0].value1, 0 | units.none)
        self.assertEquals(list[0].value2, 0 | units.none)
        self.assertEquals(list[1].value1, 0 | units.none)
        self.assertEquals(list[1].value2, 1 | units.none)
        for x in range(len(list)):
            self.assertEquals(list[x].value1, 0 | units.none)
            self.assertEquals(list[x].value2, x | units.none)
            
        list = instance.particles[1].element_list()
        for x in range(len(list)):
            self.assertEquals(list[x].value1, 1 | units.none)
            self.assertEquals(list[x].value2, x | units.none)
            
        #print instance.particles.element_list()
        #print instance.particles.element_list()[0][1].value2
        
        
        
    def test4(self):
        """
        In this test we will get a list from a particle
        """
        
        instance = ExampleParticlesInterface()
        self.assertEquals(len(instance.particles), 0)
        
        theParticle = core.Particle()
        theParticle.x = 0.1 | units.m
        theParticle.y = 0.2 | units.m
        theParticle.z = 0.5 | units.m
        
        self.assertRaises(exceptions.AmuseException, instance.particles.add_particle, theParticle)        
        

        
    def test5(self):
        """
        In this test we will get subsets from the incode set
        """
        
        instance = ExampleParticlesInterface()
        self.assertEquals(len(instance.particles), 0)
        
        particles = core.Particles(10)
        particles.mass = 0.1 | units.kg
        particles.x = 0.1 | units.m
        particles.y = 0.2 | units.m
        particles.z = 0.5 | units.m
        
        instance.particles.add_particle(particles)
        
        self.assertEquals(len(instance.particles), 10)
        
        subset = instance.particles[0:2]
        self.assertEquals(len(subset), 2)
        self.assertTrue(str(subset).find('key')> 0)
               
    def log(self, message, *arguments):
        print "IN TEST >>", message.format(*arguments)
    
    
