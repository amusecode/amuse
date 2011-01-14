#!/usr/bin/env python

import sys
import os.path

if __name__ == "__main__":
    script_directory = sys.path[0]
    if script_directory:
        amuse_directory = os.path.normpath(os.path.join(script_directory,'../src'))
        if not os.path.exists(amuse_directory):
            print "This script must be in the bin directory of the amuse project"
        sys.path.insert(0, amuse_directory)
    else:
        print "This script must run from the command prompt"
        
    from amuse.community.support import create_c

    sys.path.append(os.getcwd())
    
    if len(sys.argv) < 3:
        print "Usage: " + sys.argv[0] + "name_of_module name_of_class"
        sys.exit(1)

    name_of_module = sys.argv[1]
    name_of_class = sys.argv[2]

    if name_of_module.endswith('.py'):
        
        if not os.path.isabs(name_of_module):
            name_of_module = os.path.join(os.getcwd(), name_of_module)
        if not os.path.exists(name_of_module):
            print "Cannot find file " + name_of_module
        if not name_of_module.startswith(amuse_directory):
            print "File "+name_of_module+" must be placed under directory"+amuse_directory
        name_of_module = name_of_module[len(amuse_directory)+1:]
        name_of_module = name_of_module[:-len('.py')]
        name_of_module = name_of_module.replace(os.sep, '.')
        
    module = __import__(name_of_module,fromlist=[name_of_class])
    class_with_legacy_functions = getattr(module, name_of_class)
    
    
    uc = create_c.MakeACHeaderStringOfAClassWithLegacyFunctions()
    uc.class_with_legacy_functions = class_with_legacy_functions
    print uc.result  

    
    
    
        
    
