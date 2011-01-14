#!/usr/bin/env python

import subprocess
import os
import sys
import time
import urllib

class DownloadAthenaFromWebpage(object):
    url_template = "http://www.astro.princeton.edu/~jstone/downloads/athena/athena{version}.tar.gz"
    filename_template = "athena{version}.tar.gz"
    version = "4.0"
    
    def directory(self):
        return os.path.abspath(os.path.dirname(__file__))

    def src_directory(self):
        return os.path.join(self.directory(), 'src')
        
    def unpack_downloaded_file(self, filename):
        print "unpacking", filename
        arguments = ['tar', '-xf']
        arguments.append(filename)
        subprocess.call(
            arguments, 
            cwd = os.path.join(self.src_directory())
        )
        print "done"
        
        
    def start(self):
        if os.path.exists('src'):
            counter = 0
            while os.path.exists('src.{0}'.format(counter)):
                counter += 1
                if counter > 100:
                    print "too many backup directories"
                    break;
            os.rename('src', 'src.{0}'.format(counter))
        
        os.mkdir('src')
        
        url = self.url_template.format(version = self.version)
        filename = self.filename_template.format(version = self.version)
        print "downloading athena version",self.version,"from", url, "to", filename
        urllib.urlretrieve(url, filename = os.path.join(self.src_directory(),filename))
        print "downloading finished"
        self.unpack_downloaded_file(filename)
    
if __name__ == '__main__':
    instance = DownloadAthenaFromWebpage()
    if len(sys.argv) > 1:
        instance.version = sys.argv[1]
    instance.start()
