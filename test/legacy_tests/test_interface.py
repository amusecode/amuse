import unittest

from amuse.legacy.interface import gd
from amuse.legacy.interface.gd import GravitationalDynamics
from amuse.legacy.interface import create_definition
from amuse.legacy.support.core import RemoteFunction

class TestGravitationalDynamics(unittest.TestCase):
    def test1(self):
        x = GravitationalDynamics()
        
        function = GravitationalDynamics.new_particle
        specification = function.specification
        self.assertTrue(specification.description.startswith("Define a new particle"))
        
    def test2(self):
        specification = RemoteFunction()
        specification.name ='test'
        specification.addParameter('one','d',specification.IN, 'first parameter')
        specification.description = 'Example function'
        
        x = create_definition.CreateDescriptionOfALegacyFunctionDefinition()
        x.specification = specification
        x.start()
        
        self.assertTrue(x.out.string.find('void test(float64 one)') > 0)
        self.assertTrue(x.out.string.find('Example function') >= 0)
        self.assertTrue(x.out.string.find(':param one:') > 0)


        
    def test3(self):
        specification = RemoteFunction()
        specification.name ='test'
        specification.addParameter('one','d',specification.IN, 'first parameter')
        specification.result_type = 'i'
        specification.result_doc = 'an integer'
        specification.description = 'Example function'
        
        x = create_definition.CreateDescriptionOfALegacyFunctionDefinition()
        x.specification = specification
        x.start()
        
        self.assertTrue(x.out.string.find('int32 test(float64 one)') > 0)
        self.assertTrue(x.out.string.find(':returns:') > 0)
        
    
        
    def test4(self):
        specification = RemoteFunction()
        specification.name ='test'
        specification.addParameter('one','d',specification.IN, 'first parameter')
        specification.result_type = 'i'
        specification.result_doc = 'an integer'
        specification.description = 'Example function'
        
        x = create_definition.CreateFortranStub()
        x.specification = specification
        x.start()
        
        self.assertTrue(x.out.string.find('FUNCTION test(one)') >= 0)
        self.assertTrue(x.out.string.find('END FUNCTION') > 0)
        self.assertTrue(x.out.string.find('DOUBLE PRECISION :: one') > 0)
        
