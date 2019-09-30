from amuse.support.interface import InCodeComponentImplementation

from amuse.test.amusetest import TestWithMPI
from amuse.test import compile_tools
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

    def build_worker(self):
        
        path = os.path.abspath(self.get_path_to_results())
        codefile = os.path.join(path,"code.o")
        headerfile = os.path.join(path,"worker_code.h")
        interfacefile = os.path.join(path,"interface.o")
        self.exefile = os.path.join(path,"c_worker")
        
        compile_tools.c_compile(codefile, codestring)
        
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

        compile_tools.cxx_compile(interfacefile, code,extra_args=['-I', path])
        compile_tools.c_build(self.exefile, [interfacefile, codefile] )
    

    def setUp(self):
        super(TestCImplementationInterface, self).setUp()
        print("building...", end=' ')
        self.check_can_compile_modules()
        try:
            self.build_worker()
        except Exception as ex:
            print(ex)
            raise
        print("done")
        
    def test1(self):
        instance = ForTestingInterface(self.exefile)
        error = instance.set0(1)
        a_out, error = instance.get0()
        instance.stop()
        self.assertEqual(a_out, 1)
        self.assertEqual(error, 0)
        
    def test2(self):
        instance = ForTesting(self.exefile)
        print(instance.grid0)
        instance.grid0.a=12. | units.m
        self.assertEqual(instance.grid0.a,12.| units.m)
        instance.stop()

    def test3(self):
        instance1 = ForTesting(self.exefile)
        instance2 = ForTesting(self.exefile)
        instance1.grid0.a=12. | units.m
        instance1.grid0.new_channel_to(instance2.grid0).copy_all_attributes()
        self.assertEqual(instance2.grid0.a,12.| units.m)
        instance1.stop()
        instance2.stop()
