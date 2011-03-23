#!/usr/bin/env python

import subprocess
import os
import sys

class GetCodeFromSVN(object):
    revision = 'HEAD'
    username = 'anonymous'
    password = ''
    url = 'https://gforge.cs.vu.nl/svn/ibis/deploy/trunk'
    
    def directory(self):
        return os.path.abspath(os.path.dirname(__file__))
    
    def deploy_directory(self):
        return os.path.join(self.directory(), 'deploy')
              
    def run(self):
        arguments = [
            'svn',
            'export',
            '--trust-server-cert',
            '--non-interactive',
            '-r',
            self.revision,
            self.url,
            '--username',
            self.username,
            '--password',
            self.password,
            self.deploy_directory()
        ]
        subprocess.call(arguments)
        
if __name__ == '__main__':
    instance = GetCodeFromSVN()
    instance.run()
