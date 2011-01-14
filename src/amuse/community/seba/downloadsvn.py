#!/usr/bin/env python

import subprocess
import os
import sys
import time

class DownloadFromRemoteSvn(object):
    url = "svn://stoffel.strw.leidenuniv.nl:/data/svn/spz/SeBa/SeBaHPT"
    usermame = "amuse"
    password = "SeBa4amuse"
    
    def directory(self):
        return os.path.abspath(os.path.dirname(__file__))

    def start(self):
        arguments = ['svn', 'co']
        arguments.append('--username='+self.usermame)
        arguments.append('--password='+self.password)
        arguments.append('--no-auth-cache')
        arguments.append(self.url)
        arguments.append('src')
        subprocess.call(
            arguments, 
            cwd = self.directory()
        )
        print "Donwloaded sourcecode from svn"
        time.sleep(5) #sleeping 5 seconds, make failed because files were not added to library (assume clock failure)
    

if __name__ == '__main__':
    instance = DownloadFromRemoteSvn()
    instance.start()
