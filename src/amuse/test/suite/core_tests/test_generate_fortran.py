from amuse.test import amusetest
import numpy
import inspect
import collections


from amuse.rfi.tools import create_fortran
from amuse.rfi.core import *
class ForTestingInterface(CodeInterface):
    
    def __init__(self, exefile, **options):
        CodeInterface.__init__(self, exefile, **options)

    @legacy_function
    def echo_int():
        function = LegacyFunctionSpecification()  
        function.addParameter('int_in', dtype='int32', direction=function.IN)
        function.addParameter('int_out', dtype='int32', direction=function.OUT)
        function.result_type = 'int32'
        function.can_handle_array = True
        return function     
        
    @legacy_function
    def echo_double():
        function = LegacyFunctionSpecification()  
        function.addParameter('double_in', dtype='float64', direction=function.IN)
        function.addParameter('double_out', dtype='float64', direction=function.OUT)
        function.result_type = 'int32'
        function.can_handle_array = True
        return function  
        
        
class TestGenerateAFortranStubStringFromASpecificationClass(amusetest.TestCase):
    
    def assertContainsString(self, string, substring):
        index = string.find(substring)
        if index < 0:
            self.fail("{0} not found in {1}".format(substring, string))
            
    def assertNotContainsString(self, string, substring):
        index = string.find(substring)
        if index >= 0:
            self.fail("{0} found in {1}".format(substring, string))
        
    def test1(self):
        x = create_fortran.GenerateAFortranStubStringFromASpecificationClass()
        x.specification_class = ForTestingInterface
        x.start()
        outputstring = x.result        
        print(x.result)
        self.assertContainsString(outputstring, "FUNCTION echo_int(int_in, int_out)")
        self.assertContainsString(outputstring, "INTEGER :: echo_int")
        
    def test2(self):
        x = create_fortran.GenerateAFortranStubStringFromASpecificationClass()
        x.specification_class = ForTestingInterface
        x.start()
        outputstring = x.result        
        self.assertNotContainsString(outputstring, "internal__")
    
        
class TestGenerateAFortranSourcecodeStringFromASpecificationClass(amusetest.TestCase):
    
    def assertContainsString(self, string, substring):
        index = string.find(substring)
        if index < 0:
            self.fail("{0} not found in {1}".format(substring, string))
            
    def assertNotContainsString(self, string, substring):
        index = string.find(substring)
        if index >= 0:
            self.fail("{0} found in {1}".format(substring, string))
        
    def test1(self):
        x = create_fortran.GenerateAFortranSourcecodeStringFromASpecificationClass()
        x.specification_class = ForTestingInterface
        x.start()
        outputstring = x.result        
        self.assertContainsString(outputstring, "integer :: echo_int")
        self.assertContainsString(outputstring, "integer :: echo_double")
        self.assertContainsString(outputstring, "echo_int(")
        self.assertContainsString(outputstring, "echo_double(")
        
    def test2(self):
        x = create_fortran.GenerateAFortranSourcecodeStringFromASpecificationClass()
        x.specification_class = ForTestingInterface
        x.underscore_functions_from_specification_classes = [ForTestingInterface]
        x.start()
        outputstring = x.result        
        self.assertContainsString(outputstring, "integer :: echo_int_")
        self.assertContainsString(outputstring, "integer :: echo_double_")
        self.assertContainsString(outputstring, "echo_int_(")
        self.assertContainsString(outputstring, "echo_double_(")
    
    
        

