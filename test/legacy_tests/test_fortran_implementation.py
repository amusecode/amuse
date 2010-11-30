from amuse.support.legacy.core import *

from amuse.support.data import core
from amuse.support.units import nbody_system
from amuse.support.units import units
from amuse.support.legacy import channel

from amuse.test.amusetest import TestWithMPI
from amuse.support.legacy import create_fortran

import subprocess
import os
import time

codestring = """
function echo_int(int_in, int_out)
    implicit none
    integer :: int_in, int_out
    integer :: echo_int
    
    int_out = int_in
    
    echo_int = 0
end function

function echo_double(double_in, double_out)
    implicit none
    DOUBLE PRECISION :: double_in, double_out
    integer :: echo_double
    
    double_out = double_in
    
    echo_double = 0
end function

function echo_float(float_in, float_out)
    implicit none
    REAL(kind=4) :: float_in, float_out
    integer :: echo_float
    
    float_out = float_in
    
    echo_float = 0
end function

function echo_string(string_in, string_out)
    implicit none
    character(len=*) :: string_in, string_out
    integer :: echo_string
    
    string_out = string_in
    
    echo_string = 0
end function

function echo_strings(string_inout1, string_inout2)
    implicit none
    character(len=*) :: string_inout1, string_inout2
    integer :: echo_strings
    
    string_inout1(1:1) = 'A'
    string_inout2(1:1) = 'B'
    
    
    echo_strings = 0
end function


function return_string(string_in)
    implicit none
    character(len=*) :: string_in, return_string
    
    return_string = string_in
end function


function hello_string(string_out)
    implicit none
    character(len=30) :: string_out
    integer :: hello_string
    
    string_out = 'hello'
    
    hello_string = 0
end function

function print_string(string_in)
    implicit none
    character(len=*) :: string_in
    integer :: print_string
    
    write (*,*) string_in
    
    print_string = 0
end function

function print_error_string(string_in)
    implicit none
    character(len=*) :: string_in
    integer :: print_error_string
    
    write (0,*) string_in
    
    print_error_string = 0
end function

function echo_string_fixed_len(string_in, string_out)
    implicit none
    character(len=30) :: string_in, string_out
    integer :: echo_string_fixed_len
    
    string_out = string_in
    
    echo_string_fixed_len = 0
end function

function echo_array_with_result(int_in, int_out, N)
    implicit none
    integer, intent(in) :: N
    integer :: int_in(N), int_out(N)
    integer :: echo_array_with_result,  i
    
    do i = 1, N
     int_out(i) = int_in(i)
    end do
    
    echo_array_with_result = -1
end function



function echo_inout_array_with_result(inout, N) 
    implicit none
    integer, intent(in) :: N
    integer :: inout(N)
    integer :: echo_inout_array_with_result,  i
    
    do i = 1, N
     inout(i) = inout(i) + 10
    end do
    
    echo_inout_array_with_result = 11;
end function


function echo_logical(input, output)
    implicit none
    logical :: input, output
    integer :: echo_logical
    
    output = input
    print *, "INPUT=", input
    
    echo_logical = 0
end function

"""

class ForTestingInterface(LegacyInterface):
    
    def __init__(self, exefile, **options):
        LegacyInterface.__init__(self, exefile, **options)

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
        
    @legacy_function
    def echo_float():
        function = LegacyFunctionSpecification()  
        function.addParameter('float_in', dtype='float32', direction=function.IN)
        function.addParameter('float_out', dtype='float32', direction=function.OUT)
        function.result_type = 'int32'
        function.can_handle_array = True
        return function  
          
    @legacy_function
    def echo_string():
        function = LegacyFunctionSpecification()  
        function.addParameter('string_in', dtype='string', direction=function.IN)
        function.addParameter('string_out', dtype='string', direction=function.OUT)
        function.result_type = 'int32'
        function.can_handle_array = True
        return function  
          
    @legacy_function
    def echo_strings():
        function = LegacyFunctionSpecification()  
        function.addParameter('string_inout1', dtype='string', direction=function.INOUT)
        function.addParameter('string_inout2', dtype='string', direction=function.INOUT)
        function.result_type = 'int32'
        function.can_handle_array = True
        return function    
    
    @legacy_function
    def return_string():
        function = LegacyFunctionSpecification()  
        function.addParameter('string_in', dtype='string', direction=function.IN)
        function.result_type = 'string'
        function.can_handle_array = True
        return function  
        
    
    @legacy_function
    def hello_string():
        function = LegacyFunctionSpecification()  
        function.addParameter('string_out', dtype='string', direction=function.OUT)
        function.result_type = 'int32'
        function.can_handle_array = True
        return function  
    
    @legacy_function
    def echo_string_fixed_len():
        function = LegacyFunctionSpecification()  
        function.addParameter('string_in', dtype='string', direction=function.IN)
        function.addParameter('string_out', dtype='string', direction=function.OUT)
        function.result_type = 'int32'
        function.can_handle_array = True
        return function  
        
    @legacy_function
    def echo_array_with_result():
        function = LegacyFunctionSpecification()  
        function.addParameter('int_in', dtype='int32', direction=function.IN)
        function.addParameter('int_out', dtype='int32', direction=function.OUT)
        function.addParameter('len', dtype='int32', direction=function.LENGTH)
        function.result_type = 'int32'
        function.must_handle_array = True
        return function
            
    @legacy_function
    def echo_inout_array_with_result():
        function = LegacyFunctionSpecification()  
        function.addParameter('in_out', dtype='int32', direction=function.INOUT)
        function.addParameter('len', dtype='int32', direction=function.LENGTH)
        function.result_type = 'int32'
        function.must_handle_array = True
        return function   
    




    @legacy_function
    def echo_logical():
        function = LegacyFunctionSpecification()
        function.addParameter('input', dtype='bool', direction=function.IN)
        function.addParameter('output', dtype='bool', direction=function.OUT)
        function.result_type = 'int32'
        function.can_handle_array = True
        return function
        
    @legacy_function
    def print_string():
        function = LegacyFunctionSpecification()  
        function.addParameter('string_in', dtype='string', direction=function.IN)
        function.result_type = 'int32'
        function.can_handle_array = True
        return function  
        
    @legacy_function
    def print_error_string():
        function = LegacyFunctionSpecification()  
        function.addParameter('string_in', dtype='string', direction=function.IN)
        function.result_type = 'int32'
        function.can_handle_array = True
        return function  
    
    
    
class TestInterface(TestWithMPI):
    
    def get_mpif90_name(self):
        return os.environ['MPIF90'] if 'MPIF90' in os.environ else 'mpif90'
    
    def fortran_compile(self, objectname, string):
        if os.path.exists(objectname):
            os.remove(objectname)
        
        root, ext = os.path.splitext(objectname)
        sourcename = root + '.f90'
        
        with open(sourcename, "w") as f:
            f.write(string)
            
        process = subprocess.Popen(
            [self.get_mpif90_name(), "-g", "-c",  "-o", objectname, sourcename],
            stdin = subprocess.PIPE,
            stdout = subprocess.PIPE,
            stderr = subprocess.PIPE
        )
        stdout, stderr = process.communicate()
        
        if not os.path.exists(objectname) or process.returncode != 0:
            print "Could not compile {0}, error = {1}".format(objectname, stderr)
            raise Exception("Could not compile {0}, error = {1}".format(objectname, stderr))
            
    
    def fortran_build(self, exename, objectnames):
        arguments = [self.get_mpif90_name()]
        arguments.extend(objectnames)
        arguments.append("-o")
        arguments.append(exename)
        
        process = subprocess.Popen(
            arguments,
            stdin = subprocess.PIPE,
            stdout = subprocess.PIPE,
            stderr = subprocess.PIPE
        )
        stdout, stderr = process.communicate()
        if not os.path.exists(exename) or process.returncode != 0:
            print "Could not build {0}, error = {1}".format(objectname, stderr)
            raise Exception("Could not build {0}, error = {1}".format(exename, stderr))
    
    def build_worker(self):
        
        path = os.path.abspath(self.get_path_to_results())
        codefile = os.path.join(path,"code.o")
        interfacefile = os.path.join(path,"interface.o")
        self.exefile = os.path.join(path,"fortran_worker")
        
        self.fortran_compile(codefile, codestring)
        
        uc = create_fortran.MakeAFortranStringOfAClassWithLegacyFunctions()
        uc.class_with_legacy_functions = ForTestingInterface
        string =  uc.result
        self.fortran_compile(interfacefile, string)
        self.fortran_build(self.exefile, [interfacefile, codefile] )
    
    def setUp(self):
        super(TestInterface, self).setUp()
        print "building"
        self.build_worker()
        
    def test1(self):
        instance = ForTestingInterface(self.exefile)
        int_out, error = instance.echo_int(10)
        del instance
        self.assertEquals(int_out, 10)
        self.assertEquals(error, 0)
        
    def test2(self):
        instance = ForTestingInterface(self.exefile)
        out, error = instance.echo_double(4.0)
        del instance
        self.assertEquals(out, 4.0)
        self.assertEquals(error, 0)
        
        
    def test3(self):
        instance = ForTestingInterface(self.exefile)
        input = [1,2,3,4]
        output, errors = instance.echo_int(input)
        del instance
        self.assertEquals(len(errors),4)
        for actual, expected in zip(output, input):
            self.assertEquals(actual, expected)
            
    def test4(self):
        instance = ForTestingInterface(self.exefile)
        input = [1.0,2.1,3.3,4.2]
        output, errors = instance.echo_double(input)
        del instance
        self.assertEquals(len(errors),4)
        for actual, expected in zip(output, input):
            self.assertEquals(actual, expected)
            
        
    def test5(self):
        instance = ForTestingInterface(self.exefile)
        out, error = instance.echo_float(4.0)
        del instance
        self.assertEquals(out, 4.0)
        self.assertEquals(error, 0)
        
    def test6(self):
        
        instance = ForTestingInterface(self.exefile)
        out, error = instance.echo_string("abc")
        del instance
        self.assertEquals(error, 0)
        self.assertEquals(out[0], "abc")

    def test7(self):
        instance = ForTestingInterface(self.exefile)
        out, error = instance.echo_string(["abc","def"])
        del instance
        
        self.assertEquals(error[0], 0)
        self.assertEquals(error[1], 0)
        self.assertEquals(out[0], "abc")
        self.assertEquals(out[1], "def")

    def test8(self):
        instance = ForTestingInterface(self.exefile)
        out1, out2, error = instance.echo_strings("abc","def")
        del instance
        
        self.assertEquals(error, 0)
        self.assertEquals(out1[0], "Abc")
        self.assertEquals(out2[0], "Bef")
      
    def test9(self):
        instance = ForTestingInterface(self.exefile)
        str1_out, str2_out, error = instance.echo_strings(["abc", "def"], ["ghi", "jkl"])
        del instance
        
        self.assertEquals(error[0], 0)
        self.assertEquals(error[1], 0)
        self.assertEquals(str1_out[0], "Abc")
        self.assertEquals(str1_out[1], "Aef")
        self.assertEquals(str2_out[0], "Bhi")
        self.assertEquals(str2_out[1], "Bkl")
      
    def test10(self):
        instance = ForTestingInterface(self.exefile)
        out = instance.return_string("abc")
        del instance
        
        self.assertEquals(out[0], "abc")
        
    def test11(self):
        instance = ForTestingInterface(self.exefile)
        out, error = instance.hello_string()
        del instance
        
        self.assertEquals(out, "hello")
        
    def test12(self):
        
        instance = ForTestingInterface(self.exefile)
        out, error = instance.echo_string_fixed_len("abc")
        del instance
        self.assertEquals(error, 0)
        self.assertEquals(out[0], "abc")
        
        
    def test13(self):
        instance = ForTestingInterface(self.exefile, debugger="none")
        (output_ints, error) = instance.echo_array_with_result([4,5,6])
        instance.stop()
        print output_ints, error
        self.assertEquals(output_ints[0], 4)
        self.assertEquals(output_ints[1], 5)
        self.assertEquals(output_ints[2], 6)
        
        self.assertEquals(error[0], -1)
        self.assertEquals(error[1], -1)
        self.assertEquals(error[2], -1)
        
    
    def test14(self):
        
        for x in range(4):
            instance = ForTestingInterface(self.exefile)
            int_out, error = instance.echo_int(10)
            instance.stop()
            self.assertEquals(int_out, 10)
            self.assertEquals(error, 0)
            
    def test15(self):
        instance = ForTestingInterface(self.exefile)
        (output_ints, error) = instance.echo_inout_array_with_result([4,5,6])
        instance.stop()
        self.assertEquals(output_ints[0], 14)
        self.assertEquals(output_ints[1], 15)
        self.assertEquals(output_ints[2], 16)
        
        self.assertEquals(error[0], 11)
        self.assertEquals(error[1], 11)
        self.assertEquals(error[2], 11)


    


    def test16(self):
        instance = ForTestingInterface(self.exefile)
        (output1, error1) = instance.echo_logical(True)
        (output2, error2) = instance.echo_logical(False)
        instance.stop()
        self.assertEquals(error1, 0)
        self.assertEquals(error2, 0)
        self.assertTrue(output1)
        self.assertFalse(output2)
        
    def test20(self):
        if os.path.exists("pout.000"):
            os.remove("pout.000")
        if os.path.exists("perr.000"):
            os.remove("perr.000")
        
        x = ForTestingInterface(self.exefile, redirect_stderr_file = 'perr', redirect_stdout_file = 'pout', redirection="file")
        x.print_string("abc")
        x.print_error_string("exex")
        x.stop()
        
        time.sleep(0.2)
        
        self.assertTrue(os.path.exists("pout.000"))
        with open("pout.000","r") as f:
            content = f.read()
        self.assertEquals(content.strip(), "abc")
        
        self.assertTrue(os.path.exists("perr.000"))
        with open("perr.000","r") as f:
            content = f.read()
        self.assertEquals(content.strip(), "exex")
        
        x = ForTestingInterface(self.exefile, redirect_stderr_file = 'perr', redirect_stdout_file = 'pout', redirection="file")
        x.print_string("def")
        x.print_error_string("exex")
        x.stop()
        
        time.sleep(0.2)
        
        self.assertTrue(os.path.exists("pout.000"))
        with open("pout.000","r") as f:
            content = f.read()
        self.assertEquals(content.strip(), "abc\n def") 
        
        self.assertTrue(os.path.exists("perr.000"))
        with open("perr.000","r") as f:
            content = f.read()
        self.assertEquals(content.strip(), "exex\n exex")
        
        
    

    
    
