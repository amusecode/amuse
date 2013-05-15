from amuse.support.interface import InCodeComponentImplementation

from amuse.test.amusetest import TestWithMPI
from amuse.support import exceptions

import subprocess
import os
import time
import shutil
from amuse.units import nbody_system
from amuse.units import units
from amuse import datamodel
from amuse import config
from amuse.rfi.tools import create_java
from amuse.rfi import channel
from amuse.rfi.core import *
codestring = """
public class Code implements CodeInterface {
  
  // parameter "int_in" is an input parameter
  // parameter "int_out" is an output parameter
  // parameter "len" is a length parameter
  public int echo_array_with_result(int[] int_in, int[] int_out, int len) {
	for(int i = 0; i < len; i++) {
	    int_out[i] = int_in[i];
	}
	return -1;
  }
  
  
  // parameter "string_in" is an input parameter
  public int print_error_string(String[] string_in) {
	for(int i = 0; i < string_in.length; i++) {
    	    System.err.println(string_in[i]);
	}
	return 0;
  }
  
  
  // parameter "string_inout1" is an inout parameter
  // parameter "string_inout2" is an inout parameter
  public int echo_strings(String[] string_inout1, String[] string_inout2) {
	for(int i = 0; i < string_inout1.length; i++) {
	    String tmp = string_inout1[i];
	    string_inout1[i] = string_inout2[i];
	    string_inout2[i] = tmp;
	}
	return 0;
  }
  
  
  // parameter "double_in" is an input parameter
  // parameter "double_out" is an output parameter
  public int echo_double(double[] double_in, double[] double_out) {
	for(int i = 0; i < double_in.length; i++) {
	    double_out[i] = double_in[i];
	}
	return 0;
  }
  
  
  // parameter "float_in" is an input parameter
  // parameter "float_out" is an output parameter
  public int echo_float(float[] float_in, float[] float_out) {
	for(int i = 0; i < float_in.length; i++) {
	    float_out[i] = float_in[i];
	}
	return 0;
  }
  
  
  // parameter "int_in" is an input parameter
  // parameter "int_out" is an output parameter
  public int echo_int(int int_in, int[] int_out) {
	int_out[0] = int_in;
	if (int_in < 0) {
	    return -1;
	} else {
   	    return 0;
	}
  }

  // parameter "int_in" is an input parameter
  // parameter "int_out" is an output parameter
  public int echo_int_array(int[] int_in, int[] int_out) {
	for(int i = 0; i < int_in.length; i++) {
	    int_out[i] = int_in[i];
	}
	return 0;
  }
   
  
  // parameter "in_out" is an inout parameter
  // parameter "len" is a length parameter
  public int echo_inout_array_with_result(int[] in_out, int len) {
      for(int i = 0; i < len; i++) {
          in_out[i] = in_out[i] + 10;
      }
      return 11;
  }
  
  
  // parameter "int_in" is an input parameter
  // parameter "int_out" is an output parameter
  // parameter "len" is a length parameter
  public void echo_array(int[] int_in, int[] int_out, int len) {
	for(int i = 0; i < len; i++) {
	    int_out[i] = int_in[i];
	}
  }
  
  
  // parameter "in" is an input parameter
  // parameter "out" is an output parameter
  public int echo_long_long_int(long[] in, long[] out) {
	for(int i = 0; i < in.length; i++) {
 	    out[i] = in[i];
	}
	if (in[0] < 0) {
	    return -1;
	} else {
   	    return 0;
	}
 }
  
  
  // parameter "string_in" is an input parameter
  // parameter "string_out" is an output parameter
  public int echo_string(String[] string_in, String[] string_out) {
	for(int i = 0; i < string_in.length; i++) {
	    string_out[i] = string_in[i];
	}
	return 0;
  }
  
  
  // parameter "input" is an input parameter
  // parameter "output" is an output parameter
  public int echo_logical(boolean[] input, boolean[] output) {
	for(int i = 0; i < input.length; i++) {
	    output[i] = input[i];
	}
	return 0;
  }
  
  
  // parameter "string_in" is an input parameter
  public int print_string(String[] string_in) {
	for(int i = 0; i < string_in.length; i++) {
    	    System.out.println(string_in[i]);
	}
	return 0;
  }
  
  
  // parameter "int_in1" is an input parameter
  // parameter "int_in2" is an input parameter
  // parameter "int_out1" is an output parameter
  // parameter "int_out2" is an output parameter
  // parameter "len" is a length parameter
  public int echo_2_int(int[] int_in1, int[] int_in2, int[] int_out1, int[] int_out2, int len) {
	for(int i = 0; i < len; i++) {
            int_out1[i] = int_in1[i];
            int_out2[i] = int_in2[i];
    	}

	return len;
  }
  
}
"""

class ForTestingInterface(CodeInterface):

    classpath = ['worker.jar']
    cwd = ''
    
    def __init__(self, exefile, **options):
        CodeInterface.__init__(self, exefile, **options)

    #override default channel_type. We don't support the mpi channel...
    @option(choices=['mpi','remote','ibis', 'sockets'], sections=("channel",))
    def channel_type(self):
        return 'sockets'

    @legacy_function
    def echo_int():
        function = LegacyFunctionSpecification()
        function.addParameter('int_in', dtype='int32', direction=function.IN)
        function.addParameter('int_out', dtype='int32', direction=function.OUT)
        function.result_type = 'int32'
        return function 
            
    @legacy_function
    def echo_int_array():
        function = LegacyFunctionSpecification()
        function.addParameter('int_in', dtype='int32', direction=function.IN)
        function.addParameter('int_out', dtype='int32', direction=function.OUT)
        function.result_type = 'int32'
        function.can_handle_array = True
        return function 

    @legacy_function
    def echo_long_long_int():
        function = LegacyFunctionSpecification()
        function.addParameter('in', dtype='int64', direction=function.IN)
        function.addParameter('out', dtype='int64', direction=function.OUT)
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
    def echo_array():
        function = LegacyFunctionSpecification()  
        function.addParameter('int_in', dtype='int32', direction=function.IN)
        function.addParameter('int_out', dtype='int32', direction=function.OUT)
        function.addParameter('len', dtype='int32', direction=function.LENGTH)
        function.result_type = None
        function.must_handle_array = True
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
    
    #@legacy_function
    def return_string():
        function = LegacyFunctionSpecification()  
        function.addParameter('string_in', dtype='string', direction=function.IN)
        function.result_type = 'string'
        function.can_handle_array = True
        return function  
        
    @legacy_function
    def echo_2_int():
        function = LegacyFunctionSpecification()
        function.addParameter('int_in1', dtype='int32', direction=function.IN)
        function.addParameter('int_in2', dtype='int32', direction=function.IN, default = 1)
        function.addParameter('int_out1', dtype='int32', direction=function.OUT)
        function.addParameter('int_out2', dtype='int32', direction=function.OUT)
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
    
    
class ForTesting(InCodeComponentImplementation):
    
    def __init__(self, exefile, **options):
        InCodeComponentImplementation.__init__(self, ForTestingInterface(exefile, **options), **options)
    
    def define_methods(self, object):
        object.add_method(
            'echo_int',
            (units.m,)
            ,
            (
     		units.m,
                object.ERROR_CODE,
            )
        )



class TestInterface(TestWithMPI):
    
    def build_worker(self):
        path = os.path.abspath(self.get_path_to_results())
        self.exefile = os.path.join(path,"java_worker")
       
	#generate interface class
        interfacefile = os.path.join(path,"CodeInterface.java")

        uc = create_java.GenerateAJavaInterfaceStringFromASpecificationClass()
        uc.specification_class = ForTestingInterface
        interface =  uc.result

        with open(interfacefile, "w") as f:
            f.write(interface)
        
	#generate worker class
        workerfile = os.path.join(path,"Worker.java")

        uc = create_java.GenerateAJavaSourcecodeStringFromASpecificationClass()
        uc.specification_class = ForTestingInterface
        worker =  uc.result

        with open(workerfile, "w") as f:
            f.write(worker)

	#write code imlementation to a file

        codefile = os.path.join(path,"Code.java")

        with open(codefile, "w") as f:
            f.write(codestring)

	#compile all code

	if not config.java.is_enabled:
	    self.skip("java not enabled")

	javac = config.java.javac
	jar = config.java.jar

	jarfile = 'worker.jar'

	tmpdir = os.path.join(path, "tmp")

	shutil.rmtree(tmpdir, ignore_errors=True)
	os.mkdir(tmpdir)

        returncode = subprocess.call(
            [javac, "-d", tmpdir, "Worker.java", "CodeInterface.java", "Code.java"],
	    cwd = path,
        )
            
        if returncode != 0:
            print "Could not compile worker"

	#make jar file

	returncode = subprocess.call(
            [jar, '-cf', 'worker.jar', '-C', 'tmp', '.'],
            cwd = path,
        )
            
        if returncode != 0:
            print "Could not compile worker"

	#generate worker script

        uc = create_java.GenerateAJavaWorkerScript()
        uc.specification_class = ForTestingInterface
	uc.code_dir = path
        worker_script =  uc.result

        with open(self.exefile, "w") as f:
            f.write(worker_script)

	os.chmod(self.exefile, 0777)

    
    def setUp(self):
        super(TestInterface, self).setUp()
        print "building...",
        self.check_can_compile_modules()
        try:
            self.build_worker()
        except Exception as ex:
            print ex
            raise
        print "done"
        
    def test1(self):
        instance = ForTestingInterface(self.exefile)
        int_out, error = instance.echo_int(10)
        instance.stop()
        self.assertEquals(int_out, 10)
        self.assertEquals(error, 0)
        
    def test2(self):
        instance = ForTestingInterface(self.exefile)
        out, error = instance.echo_double(4.0)
        instance.stop()
        self.assertEquals(out, 4.0)
        self.assertEquals(error, 0)
        
        
    def test3(self):
        instance = ForTestingInterface(self.exefile)
        input = [1,2,3,4]
        output, errors = instance.echo_int_array(input)
        instance.stop()
        self.assertEquals(len(errors),4)
        for actual, expected in zip(output, input):
            self.assertEquals(actual, expected)
            
    def test4(self):
        instance = ForTestingInterface(self.exefile)
        input = [1.0,2.1,3.3,4.2]
        output, errors = instance.echo_double(input)
        instance.stop()
        self.assertEquals(len(errors),4)
        for actual, expected in zip(output, input):
            self.assertEquals(actual, expected)
            
        
    def test5(self):
        instance = ForTestingInterface(self.exefile)
        out, error = instance.echo_float(4.0)
        instance.stop()
        self.assertEquals(out, 4.0)
        self.assertEquals(error, 0)
        
    def test6(self):
        instance = ForTestingInterface(self.exefile)
        out, error = instance.echo_string("abc")
        instance.stop()
        self.assertEquals(error, 0)
        self.assertEquals(out[0], "abc")

    def test7(self):
        instance = ForTestingInterface(self.exefile)
        out, error = instance.echo_string(["abc","def"])
        instance.stop()
        
        self.assertEquals(error[0], 0)
        self.assertEquals(error[1], 0)
        self.assertEquals(out[0], "abc")
        self.assertEquals(out[1], "def")

    def test8(self):
        instance = ForTestingInterface(self.exefile)
        out1, out2, error = instance.echo_strings("abc","def")
        instance.stop()
        
        self.assertEquals(error, 0)
        self.assertEquals(out1[0], "def")
        self.assertEquals(out2[0], "abc")
      
    def test9(self):
        instance = ForTestingInterface(self.exefile)
        str1_out, str2_out, error = instance.echo_strings(["abc", "def"], ["ghi", "jkl"])
        instance.stop()
        
        self.assertEquals(error[0], 0)
        self.assertEquals(error[1], 0)
        self.assertEquals(str1_out[0], "ghi")
        self.assertEquals(str1_out[1], "jkl")
        self.assertEquals(str2_out[0], "abc")
        self.assertEquals(str2_out[1], "def")
      
    def xtest10(self):
        """test for ticket #74, 'running out of os file descriptors'
        
        
        Note: this test takes a very long time, to enable it
        remove the 'X' infront of the test name, to disable it
        add an 'X'.
        Also note: to test this, you best set the ulimit
        to a low number (but not too low), for example
        ulimit -n 400
        """
        for x in range(400):
            instance = ForTestingInterface(self.exefile)
            out, error = instance.echo_float(4.0)
            if x % 100 == 0:
                print "x:", x
            instance.stop()

    
    def test11(self):
        instance = ForTestingInterface(self.exefile)
        (output_ints,) = instance.echo_array([4,5,6])
        instance.stop()
        print output_ints
        self.assertEquals(output_ints[0], 4)
        self.assertEquals(output_ints[1], 5)
        self.assertEquals(output_ints[2], 6)
        
    def test12(self):
        instance = ForTestingInterface(self.exefile)
        (output_ints, error) = instance.echo_array_with_result([4,5,6])
        instance.stop()
        self.assertEquals(output_ints[0], 4)
        self.assertEquals(output_ints[1], 5)
        self.assertEquals(output_ints[2], 6)
        
        self.assertEquals(error[0], -1)
        self.assertEquals(error[1], -1)
        self.assertEquals(error[2], -1)
        
    def test13(self):
        instance = ForTesting(self.exefile)
        self.assertRaises(exceptions.AmuseException, instance.echo_int, [-1, -2]| units.m, 
            expected_message = "Error when calling 'echo_int' of a 'ForTesting', errorcode is -1")
        instance.stop()

    def test14(self):
        instance = ForTesting(self.exefile)
        self.assertRaises(exceptions.CodeException, lambda : instance.echo_int())
        instance.legacy_interface.echo_int.specification.id = -9
        self.assertRaises(exceptions.CodeException, lambda : instance.echo_int(1 | units.m))
        instance.stop()

    def test15(self):
        instance = ForTesting(self.exefile)
        output_ints1, output_ints2 = instance.echo_2_int([1,2],[3,4])
        output_ints3, output_ints4 = instance.echo_2_int([1,2,3])
        output_ints5, output_ints6 = instance.echo_2_int([5],[0])
        output_ints7, output_ints8 = instance.echo_2_int([5])
        instance.stop()
        self.assertEquals(output_ints1[1], 2)
        self.assertEquals(output_ints2[0], 3)
        self.assertEquals(output_ints2[1], 4)
        for i in range(3):
            self.assertEquals(output_ints3[i], i + 1)
            self.assertEquals(output_ints4[i], 1)
        self.assertEquals(output_ints5[0], 5)
        self.assertEquals(output_ints6[0], 0)
        self.assertEquals(output_ints7[0], 5)
        self.assertEquals(output_ints8[0], 1)

    def test16(self):
        instance = ForTesting(self.exefile)
        
        #self.assertRaises(exceptions.AmuseException, lambda : instance.echo_int([]))
        instance.stop()
        
    def test17(self):
        instance = ForTestingInterface(self.exefile)
        (output_ints, error) = instance.echo_inout_array_with_result([4,5,6])
        instance.stop()
        self.assertEquals(output_ints[0], 14)
        self.assertEquals(output_ints[1], 15)
        self.assertEquals(output_ints[2], 16)
        
        self.assertEquals(error[0], 11)
        self.assertEquals(error[1], 11)
        self.assertEquals(error[2], 11)

    def test18(self):
        instance = ForTestingInterface(self.exefile)
        out, error = instance.echo_logical([True, False, True])
        instance.stop()
        self.assertEquals(out, [True, False, True])
        self.assertEquals(error, 0)
        
    def test19(self):
        instance = ForTestingInterface(self.exefile)
        print 3935559000370003845
        int_out, error = instance.echo_long_long_int(3935559000370003845)
        instance.stop()
        self.assertEquals(int_out, 3935559000370003845)
        self.assertEquals(error, 0)
       
 
    def test20(self):
        path = os.path.abspath(self.get_path_to_results())
        output = os.path.join(path,"output.txt")
        error = os.path.join(path,"error.txt")
        
        if os.path.exists(output):
            os.remove(output)
        if os.path.exists(error):
            os.remove(error)

        instance = ForTesting(self.exefile, redirect_stderr_file = error, redirect_stdout_file = output, redirection="file")
        instance.print_string("abc")
        instance.print_error_string("exex")
        instance.stop()
        
        time.sleep(0.2)
        
        self.assertTrue(os.path.exists(output))
        with open(output,"r") as f:
            content = f.read()
        self.assertEquals(content.strip(), "abc")
        
        self.assertTrue(os.path.exists(error))
        with open(error,"r") as f:
            content = f.read()
        self.assertEquals(content.strip(), "exex")

    def test21(self):
        path = os.path.abspath(self.get_path_to_results())
        output = os.path.join(path,"output.txt")
        error = os.path.join(path,"error.txt")

        if os.path.exists(output):
            os.remove(output)
        if os.path.exists(error):
            os.remove(error)

        instance = ForTesting(self.exefile, redirect_stderr_file = output, redirect_stdout_file = output, redirection="file")
        instance.print_string("def")
        instance.print_error_string("exex")
        instance.stop()
        
        time.sleep(0.2)
        
        self.assertTrue(os.path.exists(output))
        with open(output,"r") as f:
            content = f.read()
        self.assertEquals(content.strip(), "def\nexex")
