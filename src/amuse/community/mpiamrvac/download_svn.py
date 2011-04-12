#!/usr/bin/env python

import subprocess
import os
import sys

class GetCodeFromSVN(object):
    revision = '149'
    username = 'studCMFA09'
    password = 'cpa9amrvac'
    url = 'https://svn.esat.kuleuven.be/amrvac/trunk'
    
    def directory(self):
        return os.path.abspath(os.path.dirname(__file__))
    
    def source_directory(self):
        return os.path.join(self.directory(), 'src')
        
    def code_directory(self):
        return os.path.join(self.source_directory(), 'mpiamrvac')
        
    def run(self):
        arguments = [
            'svn',
            'export',
            '-r',
            self.revision,
            self.url,
            '--username',
            self.username,
            '--password',
            self.password,
            self.code_directory()
        ]
        subprocess.call(arguments)
        
if __name__ == '__main__':
    instance = GetCodeFromSVN()
    instance.run()
