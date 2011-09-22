from amuse.support.core import late, OrderedDictionary

from mpi4py import MPI

import numpy
import sys
import os


from amuse.rfi.channel import ClientSideMPIMessage
from amuse.rfi.channel import pack_array
from amuse.rfi.channel import unpack_array
from amuse.rfi.core import legacy_function
from amuse.rfi.core import legacy_global
from amuse.rfi.core import LegacyFunctionSpecification
class ValueHolder(object):
    
    def __init__(self, value = None):
        self.value = value
        
    def __repr__(self):
        return "V({0!r})".format(self.value)

    def __str__(self):
        return "V({0!s})".format(self.value)

    
class PythonImplementation(object):
    dtype_to_message_attribute = { 
        'int32' : 'ints',
        'float64' : 'doubles',
        'float32' : 'floats',
        'string' : 'strings', 
        'bool' : 'booleans',
        'int64' : 'longs',
    }
    
    def __init__(self, implementation, interface):
        self.implementation = implementation
        self.interface = interface
        self.must_run = False
        
    def start(self):
        parent = MPI.Comm.Get_parent()
        parent = MPI.Comm.Get_parent()
        
        rank = MPI.COMM_WORLD.rank
        
        self.must_run = True
        while self.must_run:
            
            message = ClientSideMPIMessage()
            message.receive(parent)
                
            result_message = ClientSideMPIMessage(message.tag, message.length)
            
            if message.tag == 0:
                self.must_run = False
            else:
                if message.tag in self.mapping_from_tag_to_legacy_function:
                    try:
                        self.handle_message(message, result_message)
                    except Exception as ex:
                        print ex
                        result_message.tag = -1
                else:
                    result_message.tag = -1
            
            if rank == 0:
                result_message.send(parent)
        
        parent.Disconnect()
        
    def handle_message(self, input_message, output_message):
        legacy_function = self.mapping_from_tag_to_legacy_function[input_message.tag]
        specification = legacy_function.specification
        
        dtype_to_count = self.get_dtype_to_count(specification)
        
        
        for type, attribute in self.dtype_to_message_attribute.iteritems():
            count = dtype_to_count.get(type,0)
            for x in range(count):
                if type == 'string':
                    getattr(output_message, attribute).append([""] * output_message.length)
                else:
                    getattr(output_message, attribute).append(numpy.empty(output_message.length, dtype=type))
        
        if specification.name.startswith('internal__'):
            method = getattr(self, specification.name)
        else:
            method = getattr(self.implementation, specification.name)
        
        for type, attribute in self.dtype_to_message_attribute.iteritems():
            array = getattr(input_message, attribute)
            unpacked = unpack_array(array, input_message.length, type)
            setattr(input_message,attribute, unpacked)
            
        if specification.must_handle_array:
            keyword_arguments = self.new_keyword_arguments_from_message(input_message, None,  specification)
            result = method(**keyword_arguments)
            self.fill_output_message(output_message, None, result, keyword_arguments, specification)
        else:
            for index in range(input_message.length):
                keyword_arguments = self.new_keyword_arguments_from_message(input_message, index,  specification)
                try:
                    result = method(**keyword_arguments)
                except TypeError:
                    result = method(*list(keyword_arguments))
                self.fill_output_message(output_message, index, result, keyword_arguments, specification)
        
            
        for type, attribute in self.dtype_to_message_attribute.iteritems():
            array = getattr(output_message, attribute)
            packed = pack_array(array, input_message.length, type)
            setattr(output_message, attribute, packed)
    
    def new_keyword_arguments_from_message(self, input_message, index, specification):
        keyword_arguments = OrderedDictionary()
        for parameter in specification.parameters:
            attribute = self.dtype_to_message_attribute[parameter.datatype]
            argument_value = None
            if parameter.direction == LegacyFunctionSpecification.IN:
                if specification.must_handle_array:
                    argument_value = getattr(input_message, attribute)[parameter.input_index]
                else:
                    argument_value = getattr(input_message, attribute)[parameter.input_index][index]
            elif parameter.direction == LegacyFunctionSpecification.INOUT:
                if specification.must_handle_array:
                    argument_value = ValueHolder(getattr(input_message, attribute)[parameter.input_index])
                else:
                    argument_value = ValueHolder(getattr(input_message, attribute)[parameter.input_index][index])
            elif parameter.direction == LegacyFunctionSpecification.OUT:
                argument_value = ValueHolder(None)
            elif parameter.direction == LegacyFunctionSpecification.LENGTH:
                argument_value = input_message.length
            keyword_arguments[parameter.name] = argument_value
        return keyword_arguments
        
    def fill_output_message(self, output_message, index, result, keyword_arguments, specification):
        
        if not specification.result_type is None:
            attribute = self.dtype_to_message_attribute[specification.result_type]
            if specification.must_handle_array:
                getattr(output_message, attribute)[0] = result
            else:
                getattr(output_message, attribute)[0][index] = result
                
        for parameter in specification.parameters:
            attribute = self.dtype_to_message_attribute[parameter.datatype]
            if parameter.direction == LegacyFunctionSpecification.OUT or \
               parameter.direction == LegacyFunctionSpecification.INOUT:
                argument_value = keyword_arguments[parameter.name]
                if specification.must_handle_array:
                    getattr(output_message, attribute)[parameter.output_index] = argument_value.value
                else:
                    getattr(output_message, attribute)[parameter.output_index][index] = argument_value.value
    
    def get_dtype_to_count(self, specification):
        dtype_to_count = {}
        
        for parameter in specification.output_parameters:
            count = dtype_to_count.get(parameter.datatype, 0)
            dtype_to_count[parameter.datatype] = count + 1
                
        if not specification.result_type is None:
            count = dtype_to_count.get(specification.result_type, 0)
            dtype_to_count[specification.result_type] = count + 1
        
        return dtype_to_count
        
    @late
    def mapping_from_tag_to_legacy_function(self):
        result = {}
        for x in self.interface_functions:
            result[x.specification.id] = x
        return result
        
    @late
    def interface_functions(self):
        attribute_names = dir(self.interface)
        interface_functions = []
        for x in attribute_names:
            if x.startswith('__'):
                continue
            value = getattr(self.interface, x)
            if isinstance(value, legacy_function):
                interface_functions.append(value)
        
        interface_functions.sort(key= lambda x: x.specification.id)
        
        for x in interface_functions:
            x.specification.prepare_output_parameters()
            
        return interface_functions
        
    def internal__redirect_outputs(self, stdoutfile, stderrfile):
        mpi_rank = MPI.COMM_WORLD.rank
        sys.stdin.close()
        try:
            os.close(0)
        except Exception as ex:
            print ex
            
        if stdoutfile != "none":
            if stdoutfile != "/dev/null":
                fullname = "{0:s}.{1:03d}".format(stdoutfile, mpi_rank)
            else:
                fullname = stdoutfile
                
        
            sys.stdout.close()
            sys.stdout = open(fullname, "a+")
            
        if stderrfile != "none":
            if stderrfile != "/dev/null": 
                fullname = "{0:s}.{1:03d}".format(stderrfile, mpi_rank)
            else:
                fullname = stderrfile
                
            sys.stderr.close()
            sys.stderr = open(fullname, "a+") 
    
        return 0
        
