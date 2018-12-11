from amuse.support.core import late, OrderedDictionary

from mpi4py import MPI

import numpy
import sys
import os
import socket
import traceback
import types
import warnings

from amuse.rfi.channel import ClientSideMPIMessage
from amuse.rfi.channel import SocketMessage

from amuse.rfi.channel import pack_array
from amuse.rfi.channel import unpack_array
from amuse.rfi.core import legacy_function
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
        self.polling_interval = 0
        self.communicators = []
        self.lastid = -1
        self.activeid = -1
        self.id_to_activate = -1
        if not self.implementation is None:
            self.implementation._interface = self
        

    def start(self, mpi_port = None):
        if mpi_port is None:
            parent = self.intercomm
            self.communicators.append(parent)
        else:
            parent = MPI.COMM_WORLD.Connect(mpi_port, MPI.INFO_NULL, 0)
            self.communicators.append(parent)
        self.activeid = 0
        self.lastid += 1
        
        rank = parent.Get_rank()
        
        self.must_run = True
        while self.must_run:
            if self.id_to_activate >= 0 and self.id_to_activate != self.activeid:
                warnings.warn("activating: "+str(self.id_to_activate))
                self.activeid = self.id_to_activate
                self.id_to_activate = -1
                parent = self.communicators[self.activeid]
                rank = parent.Get_rank()
            message = ClientSideMPIMessage(polling_interval = self.polling_interval)
            message.receive(parent)

            result_message = ClientSideMPIMessage(message.call_id, message.function_id, message.call_count)
            
            if message.function_id == 0:
                self.must_run = False
            else:
                if message.function_id in self.mapping_from_tag_to_legacy_function:
                    try:
                        self.handle_message(message, result_message)
                    except Exception as ex:
                        warnings.warn(str(ex))
                        traceback.print_exc()
                        result_message.set_error(str(ex))
                        #for type, attribute in self.dtype_to_message_attribute.iteritems():
                        #    setattr(result_message, attribute, [])
                            
                        for type, attribute in self.dtype_to_message_attribute.iteritems():
                            array = getattr(result_message, attribute)
                            packed = pack_array(array, result_message.call_count, type)
                            setattr(result_message, attribute, packed)
                        
                else:
                    result_message.set_error("unknown function id " + str(message.function_id))
            
            if rank == 0:
                result_message.send(parent)

        if self.must_disconnect:
            for x in self.communicators:
                x.Disconnect()
            
        
    



    def start_socket(self, port, host):
        client_socket = socket.create_connection((host, port))
        
        self.must_run = True
        while self.must_run:
            
            message = SocketMessage()
            message.receive(client_socket)
                
            result_message = SocketMessage(message.call_id, message.function_id, message.call_count)
            
            if message.function_id == 0:
                self.must_run = False
            else:
                if message.function_id in self.mapping_from_tag_to_legacy_function:
                    try:
                        self.handle_message(message, result_message)
                    except  BaseException as ex:
                        traceback.print_exc()
                        result_message.set_error(ex.__str__())
                        for type, attribute in self.dtype_to_message_attribute.iteritems():
                            array = getattr(result_message, attribute)
                            packed = pack_array(array, result_message.call_count, type)
                            setattr(result_message, attribute, packed)
    
                else:
                    result_message.set_error("unknown function id " + message.function_id)
            
            result_message.send(client_socket)
        
        client_socket.close()
        

    def handle_message(self, input_message, output_message):
        legacy_function = self.mapping_from_tag_to_legacy_function[input_message.function_id]
        specification = legacy_function.specification
        dtype_to_count = self.get_dtype_to_count(specification)
        
        
        if specification.name.startswith('internal__'):
            method = getattr(self, specification.name)
        else:
            method = getattr(self.implementation, specification.name)
        
        if specification.has_units:
            input_units = self.convert_floats_to_units(input_message.encoded_units)
        else:
            input_units = ()
        
        for type, attribute in self.dtype_to_message_attribute.iteritems():
            count = dtype_to_count.get(type,0)
            for x in range(count):
                if type == 'string':
                    getattr(output_message, attribute).append([""] * output_message.call_count)
                else:
                    getattr(output_message, attribute).append(numpy.zeros(output_message.call_count, dtype=type))

        for type, attribute in self.dtype_to_message_attribute.iteritems():
            array = getattr(input_message, attribute)
            unpacked = unpack_array(array, input_message.call_count, type)
            setattr(input_message,attribute, unpacked)
        
        units = [False] * len(specification.output_parameters)
        if specification.must_handle_array:
            keyword_arguments = self.new_keyword_arguments_from_message(input_message, None,  specification, input_units)
            try: 
                result = method(**keyword_arguments)
            except TypeError, ex:
                warnings.warn("mismatch in python function specification(?): "+str(ex))
                result = method(*list(keyword_arguments))
            self.fill_output_message(output_message, None, result, keyword_arguments, specification, units)
        else:
            for index in range(input_message.call_count):
                keyword_arguments = self.new_keyword_arguments_from_message(input_message, index,  specification, input_units)
                try:
                    result = method(**keyword_arguments)
                    if result < 0:
                        warnings.warn("result <0 detected: "+str( (result, keyword_arguments) ))
                except TypeError, ex:
                    warnings.warn("mismatch in python function specification(?): "+str(ex))
                    result = method(*list(keyword_arguments))
                    if result < 0:
                        warnings.warn("result <0 detected: list "+str( (result, keyword_arguments) ))
                self.fill_output_message(output_message, index, result, keyword_arguments, specification, units)
        
            
        for type, attribute in self.dtype_to_message_attribute.iteritems():
            array = getattr(output_message, attribute)
            packed = pack_array(array, input_message.call_count, type)
            setattr(output_message, attribute, packed)
        
        if specification.has_units:
            output_message.encoded_units = self.convert_output_units_to_floats(units)
    


    def new_keyword_arguments_from_message(self, input_message, index, specification, units = []):
        keyword_arguments = OrderedDictionary()
        for parameter in specification.parameters:
            attribute = self.dtype_to_message_attribute[parameter.datatype]
            argument_value = None
            if parameter.direction == LegacyFunctionSpecification.IN:
                if specification.must_handle_array:
                    argument_value = getattr(input_message, attribute)[parameter.input_index]
                else:
                    argument_value = getattr(input_message, attribute)[parameter.input_index][index]
                if specification.has_units:
                    unit = units[parameter.index_in_input]
                    if not unit is None:
                        argument_value = argument_value | unit
            elif parameter.direction == LegacyFunctionSpecification.INOUT:
                if specification.must_handle_array:
                    argument_value = ValueHolder(getattr(input_message, attribute)[parameter.input_index])
                else:
                    argument_value = ValueHolder(getattr(input_message, attribute)[parameter.input_index][index])
                
                if specification.has_units:
                    unit = units[parameter.index_in_input]
                    if not unit is None:
                        argument_value.value = argument_value.value | unit
            elif parameter.direction == LegacyFunctionSpecification.OUT:
                argument_value = ValueHolder(None)
            elif parameter.direction == LegacyFunctionSpecification.LENGTH:
                argument_value = input_message.call_count
            name = 'in_' if parameter.name == 'in' else parameter.name
            keyword_arguments[name] = argument_value
        return keyword_arguments
        

    def fill_output_message(self, output_message, index, result, keyword_arguments, specification, units):
        from amuse.units import quantities
        
        if not specification.result_type is None:
            attribute = self.dtype_to_message_attribute[specification.result_type]
            if specification.must_handle_array:
                getattr(output_message, attribute)[0] = result
            else:
                getattr(output_message, attribute)[0][index] = result
                
        for parameter in specification.parameters:
            attribute = self.dtype_to_message_attribute[parameter.datatype]
            if (parameter.direction == LegacyFunctionSpecification.OUT or 
                parameter.direction == LegacyFunctionSpecification.INOUT):
                argument_value = keyword_arguments[parameter.name]
                output = argument_value.value
                if specification.has_units:
                    unit = output.unit if quantities.is_quantity(output) else None
                    if specification.must_handle_array or index == 0:
                        units[parameter.index_in_output] = unit
                    else:
                        unit = units[parameter.index_in_output]
                    if not unit is None:
                        output = output.value_in(unit)
                if specification.must_handle_array:
                    getattr(output_message, attribute)[parameter.output_index] = output
                else:
                    getattr(output_message, attribute)[parameter.output_index][index] = output
    
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
    
    

    def internal__set_message_polling_interval(self, inval):
        self.polling_interval = inval
        return 0
    
    def internal__get_message_polling_interval(self, outval):
        outval.value = self.polling_interval 
        return 0
        
    def get_null_info(self):
        return getattr(MPI, 'INFO_NULL') if hasattr(MPI, 'INFO_NULL') else None
        
    def internal__open_port(self, port_identifier):
        port_identifier.value = MPI.Open_port(self.get_null_info())
        return 0
        
    def internal__accept_on_port(self, port_identifier, comm_identifier):
        new_communicator = None
        rank = MPI.COMM_WORLD.Get_rank()
        if rank == 0:
            communicator = MPI.COMM_SELF.Accept(port_identifier, self.get_null_info(), 0)
            merged = communicator.Merge(False)
            new_communicator = MPI.COMM_WORLD.Create_intercomm(0, merged, 1, 65)
            merged.Disconnect()
            communicator.Disconnect()
        else:
            new_communicator = MPI.COMM_WORLD.Create_intercomm(0, MPI.COMM_WORLD, 1, 65)
        
        self.communicators.append(new_communicator)
        self.lastid += 1
        comm_identifier.value = self.lastid
        return 0
    
    
    def internal__connect_to_port(self, port_identifier, comm_identifier):
        new_communicator = None
        rank = MPI.COMM_WORLD.Get_rank()
        if rank == 0:
            communicator = MPI.COMM_SELF.Connect(port_identifier, self.get_null_info(), 0)
            merged = communicator.Merge(True)
            new_communicator = MPI.COMM_WORLD.Create_intercomm(0, merged, 0, 65)
            merged.Disconnect()
            communicator.Disconnect()
        else:
            new_communicator = MPI.COMM_WORLD.Create_intercomm(0, MPI.COMM_WORLD, 0, 65)
        
        self.communicators.append(new_communicator)
        self.lastid += 1
        comm_identifier.value = self.lastid
        return 0
        
    def internal__activate_communicator(self, comm_identifier):
        if comm_identifier > self.lastid or comm_identifier < 0:
            return -1
        self.id_to_activate = comm_identifier
        return 0
        
    
    
    def internal__redirect_outputs(self, stdoutfile, stderrfile):
        mpi_rank = MPI.COMM_WORLD.rank
        sys.stdin.close()
        try:
            os.close(0)
        except Exception as ex:
            warnings.warn( str(ex))
            
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
        
    def convert_to_unit(self, units_as_floats, index):
        return None

    def convert_unit_to_floats(self, unit):
        if unit is None:
            return numpy.zeros(9, dtype=numpy.float64)
        else:
            return unit.to_array_of_floats()

    def convert_output_units_to_floats(self, units):
        result = numpy.zeros(len(units) * 9, dtype = numpy.float64)
        for index, unit in enumerate(units):
            offset = index*9
            result[offset:offset+9] = self.convert_unit_to_floats(unit)
        return result
        

    def convert_float_to_unit(self, floats):
        from amuse.units import core
        from amuse.units import units
        if numpy.all(floats == 0):
            return None
        factor = floats[0]
        result = factor
        system_index = floats[1]
        unit_system = None
        for x in core.system.ALL.values():
            if x.index == system_index:
                unit_system = x
                break
        for x in unit_system.bases:
            power = floats[x.index + 2]
            if not power == 0.0:
                result = result * (x ** power)
        return result
        


    def convert_floats_to_units(self, floats):
        result = []
        for index in range(len(floats) // 9):
            offset = index*9
            unit_floats = floats[offset:offset+9]
            unit = self.convert_float_to_unit(unit_floats)
            result.append(unit)
        return result
        


    @late
    def intercomm(self):
        return MPI.Comm.Get_parent()
    @late
    def must_disconnect(self):
        return True

    def internal__become_code(self, number_of_workers, modulename, classname):
        warnings.warn(" possible experimental code path?")
        #~ print number_of_workers, modulename, classname
        world = self.freeworld
        color = 0 if world.rank < number_of_workers else 1
        key = world.rank if world.rank < number_of_workers else world.rank - number_of_workers
        #~ print "CC,", color, key, world.rank, world.size
        newcomm = world.Split(color, key)
        #~ print ("nc:", newcomm.size, newcomm.rank)
        #~ print ("AA", self.world, color, self.world.rank, self.world.size)
        try:
            new_intercomm = newcomm.Create_intercomm(0, self.world, 0, color)
        except Exception as ex:
            warnings.warn(str(ex))
            raise ex
        #~ print ("nccc:", new_intercomm.Get_remote_size(), new_intercomm.rank)
        
        self.communicators.append(new_intercomm)
        self.id_to_activate = len(self.communicators) - 1
        self.freeworld = newcomm
        return 0
        



class CythonImplementation(PythonImplementation):
    
        
    def handle_message(self, input_message, output_message):
        legacy_function = self.mapping_from_tag_to_legacy_function[input_message.function_id]
        specification = legacy_function.specification
        
        dtype_to_count = self.get_dtype_to_count(specification)
        
        if specification.name == '_stop_worker':
            method = lambda : None
        elif specification.name.startswith('internal__'):
            method = getattr(self, specification.name)
        else:
            method = getattr(self.implementation, specification.name)
        
        if specification.has_units:
            input_units = self.convert_floats_to_units(input_message.encoded_units)
        else:
            input_units = ()
        
        for type, attribute in self.dtype_to_message_attribute.iteritems():
            count = dtype_to_count.get(type,0)
            for x in range(count):
                if type == 'string':
                    getattr(output_message, attribute).append([""] * output_message.call_count)
                else:
                    getattr(output_message, attribute).append(numpy.zeros(output_message.call_count, dtype=type))

        for type, attribute in self.dtype_to_message_attribute.iteritems():
            array = getattr(input_message, attribute)
            unpacked = unpack_array(array, input_message.call_count, type)
            setattr(input_message,attribute, unpacked)
        
        units = [False] * len(specification.output_parameters)
        if specification.must_handle_array:
            keyword_arguments = self.new_keyword_arguments_from_message(input_message, None,  specification, input_units)
            result = method(**keyword_arguments)
            self.fill_output_message(output_message, None, result, keyword_arguments, specification, units)
        else:
            for index in range(input_message.call_count):
                #print "INDEX:", index
                keyword_arguments = self.new_keyword_arguments_from_message(input_message, index,  specification, input_units)
                try:
                    result = method(**keyword_arguments)
                except TypeError as ex:
                    result = method(*list(keyword_arguments))
                self.fill_output_message(output_message, index, result, keyword_arguments, specification, units)
        
            
        for type, attribute in self.dtype_to_message_attribute.iteritems():
            array = getattr(output_message, attribute)
            packed = pack_array(array, input_message.call_count, type)
            setattr(output_message, attribute, packed)
        
        if specification.has_units:
            output_message.encoded_units = self.convert_output_units_to_floats(units)
    
        




