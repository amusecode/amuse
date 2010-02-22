#!/usr/bin/env python

import subprocess
import os
import sys

class BuildMesa(object):
    
    def mesa_directory(self):
        return os.path.dirname(__file__)
        
        
    def get_mesa_source_from_svn(self):
        mesa_url = "https://mesa.svn.sourceforge.net/svnroot/mesa/trunk"
        revision = "1943"
        subprocess.call(['svn', 'co', '-r', revision, mesa_url, 'src'], cwd = self.mesa_directory())
        subprocess.call(['cp','-f','./mesa_reqs/makefile_header_v'+str(revision),
            './src/utils/makefile_header'], cwd = self.mesa_directory())
        
    def build_mesa(self):
        self.get_mesa_source_from_svn()
        subprocess.call(['./install'], cwd = os.path.join(self.mesa_directory(), 'src'))

    def clean_mesa(self):
        subprocess.call(['./clean'], cwd = os.path.join(self.mesa_directory(), 'src'))
        print "Finished cleaning.\n"

    def very_clean_mesa(self):
        self.clean_mesa()
        print "Also emptying the EOS, KAP, and NET caches."
        subprocess.call(['./empty_caches'], cwd = os.path.join(self.mesa_directory(), 'src'))

if __name__ == '__main__':
    if len(sys.argv) == 1:
        instance = BuildMesa()
        instance.build_mesa()
    elif sys.argv[1] == "download":
        instance = BuildMesa()
        instance.get_mesa_source_from_svn()
    elif sys.argv[1] == "build":
        instance = BuildMesa()
        instance.build_mesa()
    elif sys.argv[1] == "clean":
        instance = BuildMesa()
        instance.clean_mesa()
    elif sys.argv[1] == "veryclean":
        instance = BuildMesa()
        instance.very_clean_mesa()
    else:
        print "I'm confused: unknown argument: ", sys.argv[1]
        print "Known arguments: 'download', 'build', 'clean', 'veryclean'"
