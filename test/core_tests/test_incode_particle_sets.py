from amuse.support.units import units
from amuse.support.data import core
from amuse.support.data import values
from amuse.support import interface

from amuse.test import amusetest

class ExampleParticlesInterface(interface.CodeInterface):
    """This is an example class to demonstrate how to work with incode particle sets
    using the next-level mapping interface.
    
    This class is not connected to a code, instead it provides all necessary
    methods to create, delete and update particles. For simplicity
    all attributes are stored in quantities (values with units),
    no unit conversion is provided.
    
    In this example code particles have an id (id),  position (x,y and z) and a
    mass attribute (mass).
    
    See the ExampleParticlesInterfaceTests class for examples of use.
    """
    
    def log(self, message, *arguments):
        print message.format(*arguments)
    
    def __init__(self):
        self.particles = {}
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
    
    

    def new_particle(self, mass, x, y, z):
        """Creates new particles.
        
        Note: the parameternames are directly coupled to the attribute names
        in the particle set
        
        Note: the arguments are all arrays
        """
        
        result = [] 
        for mass_element, x_element, y_element, z_element in zip(mass, x, y, z):
            particle = [self.highest_id, mass_element, x_element, y_element, z_element]
            self.particles[self.highest_id] = particle
            result.append(self.highest_id)
            self.highest_id += 1
    
    

    def delete_particle(self, index_of_the_particle):
        """Delete particles in array index_of_the_particle
        
        Note: the parametername index_of_the_particle was specified in the "define_particle_sets" method
        """
        
        for x in index_of_the_particle:
            del self.particles[x]
    
    

    def set_state(self, index_of_the_particle, mass, x, y, z):
        """Sets the mass and the position of a particle.
        
        Note: the arguments are arrays
        """
        
        for index_element, mass_element, x_element, y_element, z_element in zip(index_of_the_particle, mass, x, y, z):
            particle = self.particles[index_element]
            particle[1] = mass
            particle[2] = x
            particle[3] = y
            particle[4] = z
    
    

    def get_state(self, index_of_the_particle, mass, x, y, z):
        """Returns arrays for the mass, x, y and z values
        """
        
        massresult = values.AdaptingVectorQuantity()
        xresult = values.AdaptingVectorQuantity()
        yresult = values.AdaptingVectorQuantity()
        zresult = values.AdaptingVectorQuantity()
        
        for index_element in index_of_the_particle:
            particle = self.particles[index_element]
            
            massresult.append(particle[1])
            xresult.append(particle[2])
            yresult.append(particle[3])
            zresult.append(particle[4])
            
        return massresult, xresult, yresult, zresult
    
    

    def get_mass(self, index_of_the_particle):
        """Returns an array for the masses of the indices int the index_of_the_particle array
        """
        
        massresult = values.AdaptingVectorQuantity()
        
        for index_element in index_of_the_particle:
            particle = self.particles[index_element]
            massresult.append(particle[1])
            
        return massresult
    
    

    def get_position(self, index_of_the_particle, mass, x, y, z):
        """Returns an array of the positions for the indices in index_of_the_particle
        """
        
        xresult = values.AdaptingVectorQuantity()
        yresult = values.AdaptingVectorQuantity()
        zresult = values.AdaptingVectorQuantity()
        
        for index_element in index_of_the_particle:
            particle = self.particles[index_element]
            
            xresult.append(particle[2])
            yresult.append(particle[3])
            zresult.append(particle[4])
            
        return xresult, yresult, zresult
    
    

    def set_mass(self, index_of_the_particle, mass):
        """Sets the mass and the position of a particle.
        
        Note: the arguments are arrays
        """
        
        for index_element, mass_element in zip(index_of_the_particle, mass):
            particle = self.particles[index_element]
            
            particle[1] = mass
    
    

    def set_position(self, index_of_the_particle, x, y, z):
        """Sets the mass and the position of a particle.
        
        Note: the arguments are arrays
        """
        
        for index_element, x_element, y_element, z_element in zip(index_of_the_particle, x, y, z):
            particle = self.particles[index_element]
            
            particle[2] = x
            particle[3] = y
            particle[4] = z
    
    

class ExampleParticlesInterfaceTests(amusetest.TestCase):
    """This class runs example on the particle example particles interface
    class.
    """
    
    
    def test1(self):
        instance = ExampleParticlesInterface()
        self.assertEquals(len(instance.particles), 0)

