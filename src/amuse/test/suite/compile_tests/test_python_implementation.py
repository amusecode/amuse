from amuse.support.interface import InCodeComponentImplementation

from amuse.test.amusetest import TestWithMPI
from amuse.test.amusetest import TestCase


import numpy
import parser
import sys
import os
import time
import pickle
from amuse.units import nbody_system
from amuse.units import units
from amuse import datamodel
from amuse.rfi import python_code
from amuse.rfi.core import (
    PythonCodeInterface,
    LegacyFunctionSpecification,
    legacy_function,
)
from amuse.rfi.async_request import AsyncRequestsPool
from amuse.rfi.async_request import ASyncRequestSequence
from amuse.rfi.tools.create_python_worker import CreateAPythonWorker

from amuse.support import exceptions

class ForTestingInterface(PythonCodeInterface):

    def __init__(self, **options):
        PythonCodeInterface.__init__(self, implementation_factory=ForTestingImplementation, **options)

    @legacy_function
    def get_mass():
        function = LegacyFunctionSpecification()
        function.addParameter('index_of_the_particle', dtype='int32', direction=function.IN)
        function.addParameter('mass', dtype='float64', direction=function.OUT)
        function.result_type = 'int32'
        function.can_handle_array = True
        function.id = 10
        return function

    @legacy_function
    def set_mass():
        function = LegacyFunctionSpecification()
        function.addParameter('index_of_the_particle', dtype='int32', direction=function.IN)
        function.addParameter('mass', dtype='float64', direction=function.IN, description="The new mass of the particle")
        function.result_type = 'int32'
        function.can_handle_array = True
        function.id = 11
        return function

    @legacy_function
    def new_particle():
        function = LegacyFunctionSpecification()
        function.addParameter('mass', dtype='float64', direction=function.IN, description="The new mass of the particle")
        function.addParameter('other', dtype='int32', direction=function.IN)
        function.addParameter('index_of_the_particle', dtype='int32', direction=function.OUT)
        function.result_type = 'int32'
        function.can_handle_array = True
        return function

    @legacy_function
    def delete_particle():
        function = LegacyFunctionSpecification()
        function.addParameter('index_of_the_particle', dtype='int32', direction=function.IN)
        function.result_type = 'int32'
        function.can_handle_array = True
        return function

    @legacy_function
    def echo_int():
        function = LegacyFunctionSpecification()
        function.addParameter('int_in', dtype='int32', direction=function.IN)
        function.addParameter('int_out', dtype='int32', direction=function.OUT)
        function.result_type = 'int32'
        function.can_handle_array = True
        function.id = 12
        return function

    @legacy_function
    def echo_double():
        function = LegacyFunctionSpecification()
        function.addParameter('double_in', dtype='float64', direction=function.IN)
        function.addParameter('double_out', dtype='float64', direction=function.OUT)
        function.result_type = 'int32'
        function.can_handle_array = True
        function.id = 13
        return function

    @legacy_function
    def echo_string():
        function = LegacyFunctionSpecification()
        function.addParameter('string_in', dtype='string', direction=function.IN)
        function.addParameter('string_out', dtype='string', direction=function.OUT)
        function.result_type = 'int32'
        function.can_handle_array = True
        function.id = 14
        return function

    @legacy_function
    def echo_strings():
        function = LegacyFunctionSpecification()
        function.addParameter('string_inout1', dtype='string', direction=function.INOUT)
        function.addParameter('string_inout2', dtype='string', direction=function.INOUT)
        function.result_type = 'int32'
        function.can_handle_array = True
        function.id = 15
        return function

    @legacy_function
    def echo_bool():
        function = LegacyFunctionSpecification()
        function.addParameter('string_in', dtype='bool', direction=function.IN)
        function.addParameter('string_out', dtype='bool', direction=function.OUT)
        function.result_type = 'int32'
        function.can_handle_array = True
        function.id = 16
        return function

    @legacy_function
    def sleep():
        function = LegacyFunctionSpecification()
        function.addParameter('number_of_seconds', dtype='float64', direction=function.IN)
        function.result_type = 'int32'
        function.can_handle_array = False
        return function

    @legacy_function
    def sum_doubles():
        function = LegacyFunctionSpecification()
        function.addParameter('double_in1', dtype='float64', direction=function.IN)
        function.addParameter('double_in2', dtype='float64', direction=function.IN, default=1.0)
        function.addParameter('double_out', dtype='float64', direction=function.OUT)
        function.result_type = 'int32'
        function.can_handle_array = True
        return function

    @legacy_function
    def multiply_ints():
        function = LegacyFunctionSpecification()
        function.addParameter('int_in1', dtype='int32', direction=function.IN)
        function.addParameter('int_in2', dtype='int32', direction=function.IN)
        function.addParameter('int_out', dtype='int32', direction=function.OUT)
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

    @legacy_function
    def get_position():
        function = LegacyFunctionSpecification()
        function.addParameter('index_of_the_particle', dtype='int32', direction=function.IN)
        function.addParameter('x', dtype='float64', direction=function.OUT)
        function.addParameter('y', dtype='float64', direction=function.OUT)
        function.addParameter('z', dtype='float64', direction=function.OUT)
        function.addParameter('length', 'int32', function.LENGTH)
        function.result_type = 'int32'
        function.must_handle_array = True
        return function

    @legacy_function
    def set_position():
        function = LegacyFunctionSpecification()
        function.addParameter('index_of_the_particle', dtype='int32', direction=function.IN)
        function.addParameter('x', dtype='float64', direction=function.IN)
        function.addParameter('y', dtype='float64', direction=function.IN)
        function.addParameter('z', dtype='float64', direction=function.IN)
        function.addParameter('length', 'int32', function.LENGTH)
        function.result_type = 'int32'
        function.must_handle_array = True
        return function

    @legacy_function
    def copy_over_interface():
        function = LegacyFunctionSpecification()
        function.addParameter('comm_identifier', dtype='int32', direction=function.IN)
        function.addParameter('encoded_interface', dtype='string', direction=function.IN)
        function.result_type = 'int32'
        function.can_handle_array = False
        return function

    @legacy_function
    def deep_echo_string():
        function = LegacyFunctionSpecification()
        function.addParameter('string_in', dtype='string', direction=function.IN)
        function.addParameter('string_out', dtype='string', direction=function.OUT)
        function.result_type = 'int32'
        function.can_handle_array = True
        return function

    @legacy_function
    def return_control():
        function = LegacyFunctionSpecification()
        function.result_type = 'int32'
        function.can_handle_array = False
        return function

    @legacy_function
    def echo_quantity():
        function = LegacyFunctionSpecification()
        function.addParameter('quantity_in', dtype='float64', direction=function.IN)
        function.addParameter('quantity_out', dtype='float64', direction=function.OUT)
        function.result_type = 'int32'
        function.can_handle_array = True
        function.has_units = True
        function.id = 23
        return function

    @legacy_function
    def echo_quantities():
        function = LegacyFunctionSpecification()
        function.addParameter('quantity_in', dtype='float64', direction=function.IN)
        function.addParameter('quantity_out', dtype='float64', direction=function.OUT)
        function.result_type = 'int32'
        function.must_handle_array = True
        function.has_units = True
        function.id = 23
        return function

    @legacy_function
    def echo_quantities_error():
        function = LegacyFunctionSpecification()
        function.addParameter('quantity_in', dtype='float64', direction=function.IN)
        function.addParameter('quantity_out', dtype='float64', direction=function.OUT)
        function.result_type = 'int32'
        function.must_handle_array = True
        function.has_units = True
        function.id = 25
        return function


basic_python_exe = """#!{executable}
import sys
import os
from subprocess import call

if __name__ == '__main__':
    command = sys.argv[1:]
    
    dirname=os.path.dirname(__file__)
    dirname=os.path.abspath(dirname)
    
    logfile=os.path.join(dirname, 'pythonexe.log')
    
    with open(logfile, 'a') as stream:
        stream.write('start {{0}}\\n'.format(command[0]))
        stream.flush()
    
    command.insert(0, sys.executable)
    
    returncode = call(command, close_fds=False)
    
    with open(logfile, 'a') as stream:
        stream.write('end {{0}} {{1}}\\n'.format(command[0], returncode))
        stream.flush()
        
    sys.exit(returncode)
"""


class ForTestingImplementation(object):

    def __init__(self):
        self.masses = [0.0] * 100
        self._particle_data = numpy.reshape(numpy.arange(300.0), (-1, 3))
        self.maxindex  = 0

    def new_particle(self,  mass, other, index_of_the_particle):
        try:
            self.masses[self.maxindex] = mass
            index_of_the_particle.value = self.maxindex
            self.maxindex += 1
            return 0
        except:
            return -1

    def get_mass(self, index_of_the_particle,  mass):
        try:
            mass.value = self.masses[index_of_the_particle]
            return 0
        except:
            return -1

    def set_mass(self, index_of_the_particle,  mass):
        try:
            self.masses[index_of_the_particle] = mass
            return 0
        except:
            return -1

    def echo_int(self, int_in, int_out):
        int_out.value = int_in
        return 0

    def echo_double(self, double_in, double_out):
        double_out.value = double_in
        return 0

    def echo_string(self, string_in, string_out):
        string_out.value = string_in
        return 0

    def echo_bool(self, in_, out_):
        out_.value = in_
        return 0

    def print_string(self, string_in):
        print(string_in)
        return 0

    def print_error_string(self, string_in):
        print(string_in, file=sys.stderr)
        return 0

    def echo_strings(self, string_inout1, string_inout2):
        string_inout1.value = string_inout1.value[::-1]
        string_inout2.value = string_inout2.value[::-1]
        return 0

    def sleep(self, number_of_seconds):
        import time
        time.sleep(number_of_seconds)
        return 0

    def sum_doubles(self, double_in1, double_in2, double_out):
        double_out.value = double_in1 + double_in2
        return 0

    def multiply_ints(self, int_in1, int_in2, int_out):
        int_out.value = int_in1 * int_in2
        return 0

    def get_position(self, index_of_the_particle, x, y, z, length):
        try:
            x.value = self._particle_data[index_of_the_particle, 0]
            y.value = self._particle_data[index_of_the_particle, 1]
            z.value = self._particle_data[index_of_the_particle, 2]
            return 0
        except:
            return -1

    def set_position(self, index_of_the_particle, x, y, z, length):
        try:
            self._particle_data[index_of_the_particle, 0] = x
            self._particle_data[index_of_the_particle, 1] = y
            self._particle_data[index_of_the_particle, 2] = z
            return 0
        except:
            return -1

    def copy_over_interface(self, comm_identifier, encoded_interface):
        self._other = pickle.loads(encoded_interface.encode("latin-1"))
        self._other.channel.intercomm = self._interface.communicators[comm_identifier]
        return 0

    def deep_echo_string(self, string_in, string_out):
        result, errorcode = self._other.echo_string(string_in)
        string_out.value = result[::-1]
        return errorcode

    def return_control(self):
        self._other.internal__activate_communicator(0)
        return 0

    def echo_quantity(self, quantity_in, quantity_out):
        quantity_out.value = quantity_in * (10 | (1.0/units.s))
        return 0

    def echo_quantities(self, quantity_in, quantity_out):
        quantity_out.value = quantity_in * (10 | (1.0/units.s))
        return 0

    def echo_quantities_error(self, quantity_in, quantity_out):
        raise Exception("an unexpected event")
        return 0


class ForTesting(InCodeComponentImplementation):

    def __init__(self, **options):
        InCodeComponentImplementation.__init__(self, ForTestingInterface(**options), **options)

    def define_methods(self, object):
        object.add_method("sleep", (units.s,), (object.ERROR_CODE,))
        object.add_method("new_particle", (units.kg, object.LINK("particles")), (object.INDEX, object.ERROR_CODE,))

    def define_particle_sets(self, object):
        object.define_set('particles', 'index_of_the_particle')
        object.set_new('particles', 'new_particle')

        object.set_delete('particles', 'delete_particle')
        object.add_setter('particles', 'set_mass')
        object.add_getter('particles', 'get_mass')


class TestCreatePythonWorker(TestCase):

    def test1(self):
        x = CreateAPythonWorker()
        x.implementation_factory = ForTestingImplementation
        x.channel_type = 'mpi'
        x.interface_class = ForTestingInterface

        script_string = x.new_executable_script_string()
        self.assertTrue(script_string.find('syspath = (') > 0)
        self.assertTrue(script_string.find('ForTestingInterface') > 0)
        self.assertTrue(script_string.find('ForTestingImplementation') > 0)
        self.assertTrue(script_string.find('test_python_implementation') > 0)
        self.assertTrue(script_string.find('PythonImplementation(instance, ForTestingInterface)') > 0)
        try:
            st = compile(script_string, 'test.py', 'exec')
        except SyntaxError as ex:
            self.fail("Compilation error {0}".format(ex))

    def test2(self):
        x = CreateAPythonWorker()
        x.implementation_factory = ForTestingImplementation
        x.channel_type = 'sockets'
        x.interface_class = ForTestingInterface

        script_string = x.new_executable_script_string()
        self.assertTrue(script_string.find('syspath = (') > 0)
        self.assertTrue(script_string.find('ForTestingInterface') > 0)
        self.assertTrue(script_string.find('ForTestingImplementation') > 0)
        self.assertTrue(script_string.find('test_python_implementation') > 0)
        self.assertTrue(script_string.find('PythonImplementation(instance, ForTestingInterface)') > 0)
        self.assertTrue(script_string.find('start_socket') > 0)
        try:
            st = compile(script_string, 'test.py', 'exec')
        except SyntaxError as ex:
            self.fail("Compilation error {0}".format(ex))


class TestInterface(TestWithMPI):

    def ForTesting(self, **options):
        options["worker_dir"] = self.get_path_to_results()
        return ForTesting(**options)

    def ForTestingInterface(self, **options):
        options["worker_dir"] = self.get_path_to_results()
        return ForTestingInterface(**options)

    def test02(self):
        implementation = ForTestingImplementation()
        x = python_code.PythonImplementation(implementation, ForTestingInterface)

        input_message = python_code.ClientSideMPIMessage(0, 10, 1)
        input_message.ints = [1]

        output_message = python_code.ClientSideMPIMessage(0, 10, 1)

        x.handle_message(input_message, output_message)

        self.assertEqual(len(output_message.ints), 1)
        self.assertEqual(len(output_message.doubles), 1)
        self.assertEqual(output_message.ints[0], 0)
        self.assertEqual(output_message.doubles[0], 0.0)

    def test03(self):
        implementation = ForTestingImplementation()
        x = python_code.PythonImplementation(implementation, ForTestingInterface)

        input_message = python_code.ClientSideMPIMessage(0, 11, 1)
        input_message.ints = [1]
        input_message.doubles = [12.0]

        output_message = python_code.ClientSideMPIMessage(0, 10, 1)

        x.handle_message(input_message, output_message)

        self.assertEqual(len(output_message.ints), 1)
        self.assertEqual(len(output_message.doubles), 0)
        self.assertEqual(output_message.ints[0], 0)
        self.assertEqual(implementation.masses[1], 12.0)

    def test04(self):
        implementation = ForTestingImplementation()
        x = python_code.PythonImplementation(implementation, ForTestingInterface)

        input_message = python_code.ClientSideMPIMessage(0, 11, 4)
        input_message.ints = [1, 2, 3, 4]
        input_message.doubles = [12.0, 13.0, 14.0, 15.0]

        output_message = python_code.ClientSideMPIMessage(0, 10, 4)

        x.handle_message(input_message, output_message)

        self.assertEqual(len(output_message.ints), 4)
        self.assertEqual(len(output_message.doubles), 0)
        self.assertEqual(output_message.ints[0], 0)
        self.assertEqual(output_message.ints[3], 0)
        self.assertEqual(implementation.masses[1], 12.0)
        self.assertEqual(implementation.masses[2], 13.0)
        self.assertEqual(implementation.masses[3], 14.0)
        self.assertEqual(implementation.masses[4], 15.0)

    def test05(self):
        x = self.ForTestingInterface()

        error = x.set_mass(1, 10.0)
        self.assertEqual(error, 0)

        answer, error = x.get_mass(1)
        self.assertEqual(error, 0)
        self.assertEqual(answer, 10.0)
        x.stop()

    def test06(self):
        x = self.ForTestingInterface()

        errors = x.set_mass([1, 2], [10.0, 11.0])
        self.assertEqual(errors[0], 0)
        self.assertEqual(errors[1], 0)

        answer, errors = x.get_mass([1, 2])
        self.assertEqual(errors[0], 0)
        self.assertEqual(answer[0], 10.0)
        self.assertEqual(answer[1], 11.0)
        x.stop()

        x.stop()

    def test07(self):
        x = self.ForTestingInterface()

        int_out, error = x.echo_int(20)
        self.assertEqual(error, 0)
        self.assertEqual(int_out, 20)
        x.stop()

        x.stop()

    def test08(self):
        implementation = ForTestingImplementation()
        x = python_code.PythonImplementation(implementation, ForTestingInterface)

        input_message = python_code.ClientSideMPIMessage(0, 12, 1)
        input_message.ints = [20]

        output_message = python_code.ClientSideMPIMessage(0, 10, 1)

        x.handle_message(input_message, output_message)

        self.assertEqual(len(output_message.ints), 2)
        self.assertEqual(output_message.ints[0], 0)
        self.assertEqual(output_message.ints[1], 20)

    def test09(self):
        x = self.ForTestingInterface()
        string_out, error = x.echo_string("1234567")
        self.assertEqual(error, 0)
        self.assertEqual(string_out, "1234567")
        x.stop()

    def test10(self):
        x = self.ForTestingInterface()
        string_out, error = x.echo_string(["aaaaa", "bbbb"])
        self.assertEqual(error[0], 0)
        self.assertEqual(len(string_out), 2)
        self.assertEqual(string_out[0], "aaaaa")
        self.assertEqual(string_out[1], "bbbb")
        x.stop()

    def test11(self):
        x = self.ForTestingInterface()
        string_out, error = x.echo_string(["", "bbbb"])
        self.assertEqual(error[0], 0)
        self.assertEqual(len(string_out), 2)
        self.assertEqual(string_out[0], "")
        self.assertEqual(string_out[1], "bbbb")
        x.stop()

    def test12(self):
        x = self.ForTestingInterface()
        str1_out, str2_out, error = x.echo_strings("abc", "def")
        self.assertEqual(error, 0)
        self.assertEqual(str1_out, "cba")
        self.assertEqual(str2_out, "fed")
        x.stop()

    def test13(self):
        x = self.ForTestingInterface()
        str1_out, str2_out, error = x.echo_strings(["abc", "def"], ["ghi", "jkl"])
        self.assertEqual(error[0], 0)
        self.assertEqual(error[1], 0)
        self.assertEqual(str1_out[0], "cba")
        self.assertEqual(str1_out[1], "fed")
        self.assertEqual(str2_out[0], "ihg")
        self.assertEqual(str2_out[1], "lkj")
        x.stop()

    def test14(self):
        x = self.ForTestingInterface()
        result = x.sleep(2)
        self.assertEqual(result, 0)
        request = x.sleep.asynchronous(0.01)
        request.wait()
        result = request.result()
        self.assertEqual(result, 0)
        x.stop()

    def test15(self):
        x = self.ForTestingInterface()
        y = self.ForTestingInterface()
        request1 = x.sleep.asynchronous(0.5)
        self.assertFalse(request1.is_result_available())
        request2 = y.sleep.asynchronous(1.5)
        self.assertFalse(request2.is_result_available())
        request2.wait()
        self.assertTrue(request1.is_result_available())
        self.assertTrue(request2.is_result_available())

        self.assertEqual(request1.result(), 0)
        self.assertEqual(request2.result(), 0)

        y.stop()
        x.stop()

    def test16(self):
        x = self.ForTestingInterface()
        request1 = x.sleep.asynchronous(0.4)
        x.sleep(0.01)
        self.assertTrue(request1.is_result_available(), True)
        request1.wait()
        self.assertTrue(request1.is_result_available())
        x.sleep(0.01)
        self.assertTrue(request1.is_result_available())
        x.stop()

    def test17(self):
        x = self.ForTesting()
        self.assertTrue(x.sleep.is_async_supported)
        request = x.sleep.asynchronous(0.2 | units.s)
        request.wait()
        result = request.result()
        self.assertEqual(result, [])
        x.stop()

    def test18(self):
        print("Testing the splitting of very long MPI messages into blocks")
        x = self.ForTesting(max_message_length=10)
        N = 100
        doubles = x.echo_double([1.0*i for i in range(N)])
        self.assertTrue(list(doubles) == [1.0*i for i in range(N)])
        sums = x.sum_doubles([1.0*i for i in range(N)], [1.0*i for i in range(N)])
        self.assertTrue(list(sums) == [2.0*i for i in range(N)])
        products = x.multiply_ints(list(range(N)), list(range(N)))
        self.assertTrue(list(products) == [i*i for i in range(N)])
        N = 101
        doubles = x.echo_double([1.0*i for i in range(N)])
        self.assertTrue(list(doubles) == [1.0*i for i in range(N)])
        sums = x.sum_doubles([1.0*i for i in range(N)], [1.0*i for i in range(N)])
        self.assertTrue(list(sums) == [2.0*i for i in range(N)])
        products = x.multiply_ints(list(range(N)), list(range(N)))
        self.assertTrue(list(products) == [i*i for i in range(N)])
        x.stop()

    def test19(self):
        print("Testing the splitting of very long MPI messages into blocks II: strings")
        x = self.ForTesting(max_message_length=10)
        N = 100
        strings1, strings2 = x.echo_strings(['REDRUM' for i in range(N)], ['stressed' for i in range(N)])
        self.assertTrue(list(strings1) == ['MURDER' for i in range(N)])
        self.assertTrue(list(strings2) == ['desserts' for i in range(N)])
        N = 101
        strings1, strings2 = x.echo_strings(['REDRUM' for i in range(N)], ['stressed' for i in range(N)])
        self.assertTrue(list(strings1) == ['MURDER' for i in range(N)])
        self.assertTrue(list(strings2) == ['desserts' for i in range(N)])
        x.stop()

    def xtest20(self):
        #
        # TURNED OFF support for redirection,
        # by default output is redirected to /dev/null
        # if you need file, use the support from your mpi implementation
        #
        if os.path.exists("pout.000"):
            os.remove("pout.000")
        if os.path.exists("perr.000"):
            os.remove("perr.000")

        x = self.ForTesting(redirect_stderr_file='perr', redirect_stdout_file='pout', redirection="file")
        x.print_string("abc")
        x.print_error_string("exex")
        x.stop()

        time.sleep(0.2)

        self.assertTrue(os.path.exists("pout.000"))
        with open("pout.000", "r") as f:
            content = f.read()
        self.assertEqual(content.strip(), "abc")

        self.assertTrue(os.path.exists("perr.000"))
        with open("perr.000", "r") as f:
            content = f.read()
        self.assertEqual(content.strip(), "exex")

        x = self.ForTesting(redirect_stderr_file='pout', redirect_stdout_file='pout', redirection="file")
        x.print_string("def")
        x.print_error_string("exex")
        x.stop()

        time.sleep(0.2)

        self.assertTrue(os.path.exists("pout.000"))
        with open("pout.000", "r") as f:
            content = f.read()
        self.assertEqual(content.strip(), "abc\ndef\nexex")

    def test21(self):
        print("Testing must_handle_array for Python codes")
        instance = self.ForTestingInterface()

        x, y, z, err = instance.get_position(list(range(100)))
        self.assertEqual(err, 0)
        self.assertEqual(x, numpy.arange(0.0, 300.0, 3.0))
        self.assertEqual(y, numpy.arange(1.0, 300.0, 3.0))
        self.assertEqual(z, numpy.arange(2.0, 300.0, 3.0))
        x, y, z, err = instance.get_position(list(range(101)))
        self.assertEqual(err, -1)
        self.assertEqual(list(instance.get_position(1).values()), [3.0, 4.0, 5.0, 0])

        err = instance.set_position(list(range(100)), numpy.arange(100.0), numpy.arange(100.0, 200.0), numpy.arange(200.0, 300.0))
        self.assertEqual(err, 0)
        err = instance.set_position(list(range(101)), numpy.arange(101.0), numpy.arange(101.0, 202.0), numpy.arange(202.0, 303.0))
        self.assertEqual(err, -1)
        err = instance.set_position(0, -1.0, -2.0, -3.0)

        x, y, z, err = instance.get_position(list(range(100)))
        self.assertEqual(err, 0)
        self.assertEqual(x, numpy.concatenate(([-1.0], numpy.arange(1.0, 100.0))))
        self.assertEqual(y, numpy.concatenate(([-2.0], numpy.arange(101.0, 200.0))))
        self.assertEqual(z, numpy.concatenate(([-3.0], numpy.arange(201.0, 300.0))))

        instance.stop()

    def test22(self):

        pool = AsyncRequestsPool()

        x = self.ForTestingInterface()
        y = self.ForTestingInterface()
        request1 = x.sleep.asynchronous(0.5)
        request2 = y.sleep.asynchronous(1.5)
        finished_requests = []

        def handle_result(request, index):
            self.assertTrue(request.is_result_available())
            finished_requests.append(index)

        pool.add_request(request1, handle_result, [1])
        pool.add_request(request2, handle_result, [2])

        pool.wait()
        self.assertEqual(len(finished_requests), 1)
        self.assertEqual(len(pool), 1)

        pool.wait()
        self.assertEqual(len(finished_requests), 2)
        self.assertEqual(len(pool), 0)

        self.assertTrue(request1.is_result_available())
        self.assertTrue(request2.is_result_available())

        self.assertEqual(request1.result(), 0)
        self.assertEqual(request2.result(), 0)

        y.stop()
        x.stop()

    def test22b(self):

        pool = AsyncRequestsPool()

        x = self.ForTestingInterface()
        y = self.ForTestingInterface()
        request1 = x.sleep.asynchronous(0.2)
        request2 = y.sleep.asynchronous(0.2)
        finished_requests = []

        def handle_result(request, index):
            self.assertTrue(request.is_result_available())
            finished_requests.append(index)

        pool.add_request(request1, handle_result, [1])
        pool.add_request(request2, handle_result, [2])

        time.sleep(1.0)

        pool.wait()
        pool.wait()

        self.assertEqual(len(finished_requests), 2)
        self.assertEqual(len(pool), 0)

        self.assertTrue(request1.is_result_available())
        self.assertTrue(request2.is_result_available())

        self.assertEqual(request1.result(), 0)
        self.assertEqual(request2.result(), 0)

        pool.wait()
        self.assertEqual(len(pool), 0)

        y.stop()
        x.stop()

    def test23(self):

        path = self.get_path_to_results()

        exe = os.path.join(path, "pythonexe")
        log = os.path.join(path, "pythonexe.log")

        if os.path.exists(exe):
            os.remove(exe)

        if os.path.exists(log):
            os.remove(log)

        string = basic_python_exe.format(executable=sys.executable)

        with open(exe, 'w') as f:
            f.write(string)

        os.chmod(exe, 0o777)

        instance = self.ForTestingInterface(
            use_python_interpreter=True,
            python_interpreter=exe
        )
        x, y, z, err = instance.get_position(list(range(100)))
        self.assertEqual(err, 0)
        self.assertEqual(x, numpy.arange(0.0, 300.0, 3.0))
        self.assertEqual(y, numpy.arange(1.0, 300.0, 3.0))
        self.assertEqual(z, numpy.arange(2.0, 300.0, 3.0))

        instance.stop()
        time.sleep(0.3)

        self.assertTrue(os.path.exists(log))

        with open(log, 'r') as f:
            loglines = f.read().splitlines()

        self.assertEqual(len(loglines), 2)

        self.assertTrue(loglines[0].startswith('start '))
        self.assertTrue(loglines[1].startswith('end '))

    def test24(self):

        # same as test23 but now with redirection is none
        path = self.get_path_to_results()

        exe = os.path.join(path, "pythonexe")
        log = os.path.join(path, "pythonexe.log")

        if os.path.exists(exe):
            os.remove(exe)

        if os.path.exists(log):
            os.remove(log)

        string = basic_python_exe.format(executable=sys.executable)

        with open(exe, 'w') as f:
            f.write(string)

        os.chmod(exe, 0o777)

        instance = self.ForTestingInterface(
            use_python_interpreter=True,
            python_interpreter=exe,
            redirection="none"
        )
        x, y, z, err = instance.get_position(list(range(100)))
        self.assertEqual(err, 0)
        self.assertEqual(x, numpy.arange(0.0, 300.0, 3.0))
        self.assertEqual(y, numpy.arange(1.0, 300.0, 3.0))
        self.assertEqual(z, numpy.arange(2.0, 300.0, 3.0))

        instance.stop()
        time.sleep(0.3)

        self.assertTrue(os.path.exists(log))

        with open(log, 'r') as f:
            loglines = f.read().splitlines()

        self.assertEqual(len(loglines), 2)

        self.assertTrue(loglines[0].startswith('start '))
        self.assertTrue(loglines[1].startswith('end '))

    def test25(self):
        self.check_for_mpi()
        instance = self.ForTestingInterface(polling_interval_in_milliseconds=100)
        (output1, error1) = instance.internal__get_message_polling_interval()
        instance.stop()
        self.assertEqual(error1, 0)
        self.assertEqual(output1, 100000)

    def test25(self):
        instance = self.ForTestingInterface(polling_interval_in_milliseconds=100)
        if instance.channel.is_polling_supported():
            (output1, error1) = instance.internal__get_message_polling_interval()
            self.assertEqual(error1, 0)
            self.assertEqual(output1, 100000)
        instance.stop()

    def test26(self):
        self.check_for_mpi()
        instance1 = self.ForTestingInterface()
        instance2 = self.ForTestingInterface()
        portname, error = instance1.internal__open_port()
        self.assertTrue(len(portname) > 0)
        self.assertEqual(error, 0)
        request1 = instance1.internal__accept_on_port.asynchronous(portname)
        request2 = instance2.internal__connect_to_port.asynchronous(portname)
        request1.wait()
        request2.wait()
        port_id1, error1 = request1.result()
        port_id2, error2 = request2.result()
        self.assertTrue(port_id1 >= 0)
        self.assertTrue(port_id2 >= 0)
        self.assertEqual(error1, 0)
        self.assertEqual(error2, 0)

    def test27(self):
        self.check_for_mpi()
        instance1 = self.ForTestingInterface(redirection="none")
        instance2 = self.ForTestingInterface(redirection="none")
        encoded_interface = pickle.dumps(instance1, 0)
        decoded_interface = pickle.loads(encoded_interface)
        #pickle.loads(pickle.dumps(instance1,0))
        portname, error = instance2.internal__open_port()
        request1 = instance2.internal__accept_on_port.asynchronous(portname)
        request2 = instance1.internal__connect_to_port.asynchronous(portname)
        request1.wait()
        request2.wait()
        port_id1, error1 = request1.result()
        port_id2, error2 = request2.result()
        instance2.copy_over_interface(port_id2, pickle.dumps(instance1, 0).decode('latin-1'))
        instance1.internal__activate_communicator(port_id1)
        result, errorcode = instance2.deep_echo_string("hello")
        self.assertEqual(errorcode, 0)
        self.assertEqual(result, "olleh")
        result, errorcode = instance2.deep_echo_string("world")
        self.assertEqual(errorcode, 0)
        self.assertEqual(result, "dlrow")
        instance2.return_control()
        result, errorcode = instance1.echo_string("world")
        self.assertEqual(errorcode, 0)
        self.assertEqual(result, "world")

    def test28(self):
        x = self.ForTestingInterface()
        def next_request(index):
            if index < 3:
                return x.sleep.asynchronous(0.1)
            else:
                return None

        sequence = ASyncRequestSequence(next_request)
        self.assertFalse(sequence.is_finished)
        sequence.wait()
        self.assertTrue(sequence.is_finished)
        result = sequence.result()
        self.assertEqual(len(result), 3)
        x.stop()

    def test29(self):

        pool = AsyncRequestsPool()

        x = self.ForTestingInterface()
        y = self.ForTestingInterface()
        sequenced_requests_indices = []
        def next_request(index):
            if index < 4:
                sequenced_requests_indices.append(index)
                return x.sleep.asynchronous(0.5)
            else:
                return None

        request1 = ASyncRequestSequence(next_request)
        request2 = y.sleep.asynchronous(1.0)
        finished_requests = []

        def handle_result(request, index):
            self.assertTrue(request.is_result_available())
            self.assertTrue(request.is_finished)
            finished_requests.append(index)

        pool.add_request(request1, handle_result, [1])
        pool.add_request(request2, handle_result, [2])

        pool.wait()
        self.assertEqual(len(finished_requests), 1)
        self.assertEqual(len(pool), 1)
        self.assertEqual(finished_requests, [2])
        self.assertTrue(len(sequenced_requests_indices) > 0)

        pool.wait()
        self.assertEqual(len(finished_requests), 2)
        self.assertEqual(len(pool), 0)
        x.sleep(0.1)
        self.assertEqual(sequenced_requests_indices, [0, 1, 2, 3])

        self.assertTrue(request1.is_result_available())
        self.assertTrue(request2.is_result_available())

        self.assertEqual(request1.result(), [0, 0, 0, 0])
        self.assertEqual(request2.result(), 0)

        y.stop()
        x.stop()

    def test30(self):
        instance = self.ForTesting()
        input = [1.0, 2.0, 3.0]
        output =  instance.sum_doubles(input, 5)
        self.assertAlmostRelativeEquals(output, [6.0, 7.0, 8.0])
        output =  instance.sum_doubles(5, input)
        self.assertAlmostRelativeEquals(output, [6.0, 7.0, 8.0])

    def test31(self):
        x = self.ForTesting()
        p = datamodel.Particles(5)
        p.mass = [1, 2, 3, 4, 5] | units.kg
        p.other = None
        for pi in p:
            x.particles.add_particle(pi)
        self.assertAlmostRelativeEquals(x.particles.mass, [1, 2, 3, 4, 5])
        x.stop()

    def test32(self):
        x = self.ForTestingInterface()
        quantity_out, error = x.echo_quantity(20.0 | units.m)
        self.assertEqual(error, 0)
        self.assertEqual(quantity_out, 200 | (units.m/units.s))
        quantity_out, error = x.echo_quantity(30)
        self.assertEqual(error, 0)
        self.assertEqual(quantity_out, 300 | (1.0/units.s))
        x.stop()

    def test33(self):
        x = self.ForTestingInterface()
        quantity_out, error = x.echo_quantity([20, 30, 40] | units.m)
        self.assertEqual(error, 0)
        self.assertEqual(quantity_out, [200, 300, 400] | (units.m/units.s))
        x.stop()

    def test34(self):
        x = self.ForTestingInterface()
        #self.assertException(x.echo_quantities_error, [20, 30, 40] | units.m)
        quantity_out, error = x.echo_quantities([20, 30, 40] | units.m)
        self.assertEqual(error, 0)
        self.assertEqual(quantity_out, [200, 300, 400] | (units.m/units.s))
        x.stop()

    def test35(self):
        x = self.ForTesting(max_message_length=10)
        N = 10
        doubles = x.echo_double([1.0*i for i in range(N)])
        self.assertTrue(list(doubles) == [1.0*i for i in range(N)])
        sums = x.sum_doubles([3.0*i for i in range(N)])
        self.assertTrue(list(sums) == [3.0*i + 1 for i in range(N)])
        N = 11
        doubles = x.echo_double([1.0*i for i in range(N)])
        self.assertTrue(list(doubles) == [1.0*i for i in range(N)])
        sums = x.sum_doubles([3.0*i for i in range(N)])
        self.assertTrue(list(sums) == [3.0*i + 1 for i in range(N)])
        x.stop()

    def test36(self):
        x = self.ForTestingInterface()
        self.assertRaises(exceptions.CodeException, x.echo_quantities_error, ([20, 30, 40] | units.m), expected_message="Exception when calling function 'echo_quantities_error', of code 'ForTestingInterface', exception was 'Error in code: an unexpected event'")
        x.stop()

    def test37(self):
        x = self.ForTestingInterface()
        request = x.echo_quantity.asynchronous([20, 30, 40] | units.m)
        quantity_out, error = request.result()
        self.assertEqual(error, 0)
        self.assertEqual(quantity_out, [200, 300, 400] | (units.m/units.s))
        x.stop()

    def test40(self):
        x = self.ForTesting()
        out = x.echo_bool([True, False, True])
        self.assertEqual(out, [True, False, True])
        x.stop()

