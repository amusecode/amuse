from amuse.support.interface import InCodeComponentImplementation

from amuse.test.amusetest import TestWithMPI
from amuse.support import exceptions

import subprocess
import os
import time
from amuse.units import nbody_system
from amuse.units import units
from amuse import datamodel
from amuse.rfi.tools import create_c
from amuse.rfi import channel
from amuse.rfi.core import *
codestring = """
#include <stdio.h>
#include <stdlib.h>

double grid0;
double grid1[10];
double grid2[10][10];

int set0(double in) {
    grid0=in;
    return 0;
}
int get0(double *out) {
    *out=grid0;
    return 0;
}
int get_grid0_range() {
    return 0;
}

"""

class ForTestingInterface(CodeInterface):
    include_headers = ['worker_code.h']
 
    def __init__(self, exefile, **options):
        CodeInterface.__init__(self, exefile, **options)

    @remote_function
    def set0(a=0. | units.m):
        returns ()

    @remote_function
    def get0():
        returns (a=0. | units.m)

    @remote_function
    def get_grid0_range():
        returns ()


class ForTesting(InCodeComponentImplementation):
    
    def __init__(self, exefile, **options):
        InCodeComponentImplementation.__init__(self, ForTestingInterface(exefile, **options), **options)
    
    def define_grids(self,object):
        object.define_grid('grid0', grid_class=datamodel.RectilinearGrid,state_guard="before_new_set_instance")
        object.set_grid_range('grid0', 'get_grid0_range')
        object.add_getter('grid0', 'get0',  names=['a'])
        object.add_setter('grid0', 'set0',  names=['a'])


class TestCImplementationInterface(TestWithMPI):
    
    def get_mpicc_name(self):
        try:
            from amuse import config
            is_configured = hasattr(config, 'mpi')
        except ImportError:
            is_configured = False
    
        if is_configured:
            return config.mpi.mpicc
        else:
            return os.environ['MPICC'] if 'MPICC' in os.environ else 'mpicc'
            
    def get_mpicxx_name(self):
        try:
            from amuse import config
            is_configured = hasattr(config, 'mpi')
        except ImportError:
            is_configured = False
    
        if is_configured:
            return config.mpi.mpicxx
        else:
            return os.environ['MPICXX'] if 'MPICXX' in os.environ else 'mpicxx'
    
    def get_mpicc_flags(self):
        try:
            from amuse import config
            is_configured = hasattr(config, 'compilers')
        except ImportError:
            is_configured = False
        
        if is_configured:
            return config.compilers.cc_flags
        else:
            return ""
            
    def get_mpicxx_flags(self):
        try:
            from amuse import config
            is_configured = hasattr(config, 'compilers')
        except ImportError:
            is_configured = False
        
        if is_configured:
            return config.compilers.cxx_flags
        else:
            return ""
            
    def wait_for_file(self, filename):
        for dt in [0.01, 0.01, 0.02, 0.05]:
            if os.path.exists(filename):
                return
            time.sleep(dt)
        
    def c_compile(self, objectname, string):
        root, ext = os.path.splitext(objectname)
        sourcename = root + '.c'
        
        if os.path.exists(objectname):
            os.remove(objectname)
        
        with open(sourcename, "w") as f:
            f.write(string)
            
        mpicc = self.get_mpicc_name()
        arguments = [mpicc]
        arguments.extend(self.get_mpicc_flags().split())
        arguments.extend(["-I", "lib/stopcond", "-c",  "-o", objectname, sourcename])
            
        process = subprocess.Popen(
            arguments,
            stdin = subprocess.PIPE,
            stdout = subprocess.PIPE,
            stderr = subprocess.PIPE
        )
        stdout, stderr = process.communicate()
        
        if process.returncode == 0:
            self.wait_for_file(objectname)
        
        if process.returncode != 0 or not os.path.exists(objectname):
            print "Could not compile {0}, error = {1}".format(objectname, stderr)
            raise Exception("Could not compile {0}, error = {1}".format(objectname, stderr))
    
    def cxx_compile(self, objectname, string):
        root, ext = os.path.splitext(objectname)
        sourcename = root + '.cc'
        
        if os.path.exists(objectname):
            os.remove(objectname)
            
        with open(sourcename, "w") as f:
            f.write(string)
            
        mpicxx = self.get_mpicxx_name()
        arguments = [mpicxx]
        arguments.extend(self.get_mpicxx_flags().split())
        arguments.extend(["-I", "lib/stopcond", "-c",  "-o", objectname, sourcename])
        
        process = subprocess.Popen(
            arguments,
            stdin = subprocess.PIPE,
            stdout = subprocess.PIPE,
            stderr = subprocess.PIPE
        )
        stdout, stderr = process.communicate()
        
        if process.returncode == 0:
            self.wait_for_file(objectname)
        
        if process.returncode != 0 or not os.path.exists(objectname):
            print "Could not compile {0}, error = {1}".format(objectname, stderr)
            raise Exception("Could not compile {0}, error = {1}".format(objectname, stderr))

        #~ print stdout
        #~ print stderr
            
    def c_build(self, exename, objectnames):
        
        if os.path.exists(exename):
            os.remove(exename)
            
        mpicxx = self.get_mpicxx_name()
        arguments = [mpicxx]
        arguments.extend(objectnames)
        arguments.append("-o")
        arguments.append(exename)
        
        if 'LIBS' in os.environ:
            libs = os.environ['LIBS'].split()
            arguments.extend(libs)
            
        process = subprocess.Popen(
            arguments,
            stdin = subprocess.PIPE,
            stdout = subprocess.PIPE,
            stderr = subprocess.PIPE
        )
        stdout, stderr = process.communicate()
        
        if process.returncode == 0:
            self.wait_for_file(exename)
            
        if process.returncode != 0 or not (os.path.exists(exename) or os.path.exists(exename+'.exe')):
            print "Could not compile {0}, error = {1}".format(exename, stderr)
            raise Exception("Could not build {0}, error = {1}".format(exename, stderr))

        #~ print stdout
        #~ print stderr
    
    def build_worker(self):
        
        path = os.path.abspath(self.get_path_to_results())
        codefile = os.path.join(path,"code.o")
        headerfile = os.path.join(path,"worker_code.h")
        interfacefile = os.path.join(path,"interface.o")
        self.exefile = os.path.join(path,"c_worker")
        
        self.c_compile(codefile, codestring)
        
        uc = create_c.GenerateACHeaderStringFromASpecificationClass()
        uc.specification_class = ForTestingInterface
        uc.needs_mpi = False
        header =  uc.result

        with open(headerfile, "w") as f:
            f.write(header)
        
        uc = create_c.GenerateACSourcecodeStringFromASpecificationClass()
        uc.specification_class = ForTestingInterface
        uc.needs_mpi = False
        code =  uc.result
        
        
        
        #print string
        
        self.cxx_compile(interfacefile, code)
        self.c_build(self.exefile, [interfacefile, codefile] )
    

    def setUp(self):
        super(TestCImplementationInterface, self).setUp()
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
        error = instance.set0(1)
        a_out, error = instance.get0()
        instance.stop()
        self.assertEquals(a_out, 1)
        self.assertEquals(error, 0)
        
    def test2(self):
        instance = ForTesting(self.exefile)
        print instance.grid0
        instance.grid0.a=12. | units.m
        self.assertEquals(instance.grid0.a,12.| units.m)
        instance.stop()

    def test3(self):
        instance1 = ForTesting(self.exefile)
        instance2 = ForTesting(self.exefile)
        instance1.grid0.a=12. | units.m
        instance1.grid0.new_channel_to(instance2.grid0).copy_all_attributes()
        self.assertEquals(instance2.grid0.a,12.| units.m)
        instance1.stop()
        instance2.stop()
