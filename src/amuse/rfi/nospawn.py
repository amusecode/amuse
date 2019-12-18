
from amuse.rfi import core
from amuse.rfi.python_code import CythonImplementation
from mpi4py import MPI
from amuse.rfi import channel
from collections import namedtuple
import sys
import importlib

Code = namedtuple("Code", ['cls', 'number_of_workers', 'args', 'kwargs'])
PythonCode = namedtuple("Code", ['cls', 'number_of_workers', 'args', 'kwargs', 'implementation_factory'])

def get_number_of_workers_needed(codes):
    result = 1
    for x in codes:
        result += x.number_of_workers
    return result

def get_color(rank, codes):
    if rank == 0:
        return 0
    else:
        index = 1
        for color, x in enumerate(codes):
            if rank >= index and rank < index + x.number_of_workers:
                return color + 1
            index += x.number_of_workers
    return len(codes) + 1 #left over ranks
            
def get_key(rank, codes):
    if rank == 0:
        return 0
    else:
        index = 1
        for color, x in enumerate(codes):
            if rank >= index and rank < index + x.number_of_workers:
                return rank - index
            index += x.number_of_workers
    return rank - (len(codes) + 1) #left over ranks

def get_code_class(rank, codes):
    if rank == 0:
        return None
    else:
        index = 1
        for color, x in enumerate(codes):
            if rank >= index and rank < index + x.number_of_workers:
                return x.cls
            index += x.number_of_workers
    return None
            
            

def start_all(codes):
    
    channel.MpiChannel.ensure_mpi_initialized()
    number_of_workers_needed = get_number_of_workers_needed(codes)
    
    world = MPI.COMM_WORLD
    rank = world.rank
    if world.size < number_of_workers_needed:
        if rank == 0:
            raise Exception("cannot start all codes, the world size ({0}) is smaller than the number of requested codes ({1}) (which is always 1 + the sum of the all the number_of_worker fields)".format(world.size, number_of_workers_needed))
        else:
            return None
    
    color = get_color(world.rank, codes)
    key = get_key(world.rank, codes)
    
    newcomm = world.Split(color, key)
    
    
    localdup = world.Dup()
    if world.rank == 0:
        result = []
        remote_leader = 1
        tag = 1
        for x in codes:
            new_intercomm = newcomm.Create_intercomm(0, localdup, remote_leader, tag)
            remote_leader += x.number_of_workers
            tag += 1
            instance = x.cls(*x.args, check_mpi = False, must_start_worker = False, **x.kwargs)
            instance.legacy_interface.channel = channel.MpiChannel('_',None)
            instance.legacy_interface.channel.intercomm = new_intercomm
            result.append(instance)
            
        world.Barrier()
        
        return result    
    else:
        code_cls = get_code_class(world.rank, codes)
        if code_cls is None:
            world.Barrier()
            return None
        
        new_intercomm = newcomm.Create_intercomm(0, localdup, 0, color)
        x = get_code(world.rank, codes)
        instance =  code_cls(*x.args, check_mpi = False, must_start_worker = False, **x.kwargs)
        interface = instance.legacy_interface
       
        if hasattr(interface, '__so_module__'):
            package, _ =  code_cls.__module__.rsplit('.',1)
            modulename = package + '.' + interface.__so_module__
            module = importlib.import_module(modulename)
            module.set_comm_world(newcomm)
        else:
            module = x.implementation_factory()
        
        
        instance = CythonImplementation(module, interface.__class__)
        instance.intercomm = new_intercomm
        instance.must_disconnect = False
        world.Barrier()
        instance.start()
        
        return None
        

def stop_all(instances):
    for x in instances:
        x.stop()
def start_empty():
    
    channel.MpiChannel.ensure_mpi_initialized()
    
    world = MPI.COMM_WORLD
    rank = world.rank
    
    color = 0 if world.rank == 0 else 1
    key = 0 if world.rank == 0 else world.rank -1 
    newcomm = world.Split(color, key)
    
    localdup = world.Dup()
    if world.rank == 0:
        result = []
        remote_leader = 1
        tag = 1
        
        new_intercomm = newcomm.Create_intercomm(0, localdup, remote_leader, tag)
        
        instance = core.CodeInterface(check_mpi = False, must_start_worker = False)
        instance.channel = channel.MpiChannel('_',None)
        instance.channel.intercomm = new_intercomm
        instance.world = localdup  
        instance.remote_leader = 1
        world.Barrier()
        
        return instance    
    else:
        
        new_intercomm = newcomm.Create_intercomm(0, localdup, 0, color)
        
        
        instance = CythonImplementation(None, core.CodeInterface)
        instance.intercomm = new_intercomm
        instance.world = localdup
        instance.freeworld = newcomm
        instance.localworld = newcomm
        instance.must_disconnect = False
        world.Barrier()
        instance.start()
        print("STOP...", world.rank)
        return None
        


def get_code(rank, codes):
    if rank == 0:
        return None
    else:
        index = 1
        for color, x in enumerate(codes):
            if rank >= index and rank < index + x.number_of_workers:
                return x
            index += x.number_of_workers
    return None
            
            


