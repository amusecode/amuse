#!/usr/bin/env python

import subprocess
import os
import sys

class BuildMesa(object):
    
    def mesa_directory(self):
        return os.path.abspath(os.path.dirname(__file__))
        
        
    def get_mesa_source_from_svn(self):
        mesa_url = "http://mesa.svn.sourceforge.net/svnroot/mesa/trunk"
        revision = "2208"
        subprocess.call(['svn', 'export', 
            '-r', revision, mesa_url, 'src'], cwd = self.mesa_directory())
        subprocess.call(['cp','-f','./mesa_reqs/makefile_header',
            './src/utils/makefile_header'], cwd = self.mesa_directory())
        subprocess.call(['cp','-f','./mesa_reqs/create_zams.f',
            './src/star/test/src/create_zams.f'], cwd = self.mesa_directory())
        subprocess.call(['cp','-f','./mesa_reqs/eos_def.f',
            './src/eos/public/eos_def.f'], cwd = self.mesa_directory())
        subprocess.call(['cp','-f','./mesa_reqs/jina_def.f',
            './src/jina/public/jina_def.f'], cwd = self.mesa_directory())
        subprocess.call(['cp','-f','./mesa_reqs/kap_def.f',
            './src/kap/public/kap_def.f'], cwd = self.mesa_directory())
        subprocess.call(['cp','-f','./mesa_reqs/net_def.f',
            './src/net/public/net_def.f'], cwd = self.mesa_directory())
        subprocess.call(['cp','-f','./mesa_reqs/star_def.f',
            './src/star/public/star_def.f'], cwd = self.mesa_directory())
        
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
    instance = BuildMesa()
    if len(sys.argv) == 1:
        instance.build_mesa()
    elif sys.argv[1] == "download":
        instance.get_mesa_source_from_svn()
    elif sys.argv[1] == "build":
        instance.build_mesa()
    elif sys.argv[1] == "rebuild":
        instance.clean_mesa()
        instance.build_mesa()
    elif sys.argv[1] == "clean":
        instance.clean_mesa()
    elif sys.argv[1] == "veryclean":
        instance.very_clean_mesa()
    else:
        print "I'm confused: unknown argument: ", sys.argv[1]
        print "Known arguments: 'download', 'build', 'rebuild', 'clean', 'veryclean'"
