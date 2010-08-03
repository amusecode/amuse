
from amuse.test.amusetest import TestWithMPI

from amuse.legacy.interface import gd
from amuse.legacy.interface.gd import GravitationalDynamicsInterface
from amuse.support.legacy.create_definition import LegacyDocStringProperty

from amuse.support.legacy import create_definition
from amuse.support.legacy.core import LegacyFunctionSpecification

class TestGravitationalDynamics(TestWithMPI):
    def test1(self):
        x = GravitationalDynamicsInterface()
        
        function = GravitationalDynamicsInterface.new_particle
        specification = function.specification
        self.assertTrue(specification.description.startswith("Define a new particle"))
        
    def test2(self):
        specification = LegacyFunctionSpecification()
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
        specification = LegacyFunctionSpecification()
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
        specification = LegacyFunctionSpecification()
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
        
    def test5(self):
        class WithLegacyDocStringProperty(object):
            
            def __init__(self):
                "orignal doc"
                pass
            
            @property
            def specification(self):
                specification = LegacyFunctionSpecification()
                specification.name ='test'
                specification.addParameter('one','d',specification.IN, 'first parameter')
                specification.result_type = 'i'
                specification.result_doc = 'an integer'
                specification.description = 'Example function'
                return specification
                
            __doc__ = LegacyDocStringProperty()
            
            
        self.assertEquals("orignal doc", WithLegacyDocStringProperty.__doc__)
        instance = WithLegacyDocStringProperty()
        instance_documentation =  WithLegacyDocStringProperty().__doc__
        self.assertTrue(instance_documentation.find('FUNCTION test(one)') >= 0)
        
    def test6(self):
        print "Testing description of Legacy Function with output parameter"
        specification = LegacyFunctionSpecification()
        specification.name ='test'
        specification.addParameter('one','d',specification.OUT, 'first parameter')
        specification.result_type = 'i'
        specification.result_doc = 'an integer'
        specification.description = 'Example function'
        
        x = create_definition.CreateDescriptionOfALegacyFunctionDefinition()
        x.specification = specification
        x.start()
        self.assertTrue(x.out.string.find('int32 test(float64 * one)') > 0)
        self.assertTrue(x.out.string.find(':returns:') > 0)
        
    def test7(self):
        print "Testing __str__ of Legacy Function"
        specification = LegacyFunctionSpecification()
        specification.name ='test'
        specification.addParameter('one','f',specification.IN, 'first parameter, type: float')
        specification.addParameter('two','d',specification.OUT, 'second parameter, type double')
        specification.result_type = 'i'
        specification.result_doc = 'an integer'
        specification.description = 'Example function'
        self.assertEquals(str(specification),"function: int test(float one)\noutput: double two, int __result")
        
