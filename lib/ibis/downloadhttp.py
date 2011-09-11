#!/usr/bin/env python

import subprocess
import os
import sys
import time
import urllib

class DownloadIbisDeployFromWebpage(object):
    url_template = "http://www.cs.vu.nl/ibis/downloads/deploy{version}.zip"
    filename_template = "deploy{version}.zip"
    unpacked_dir_template = "deploy{version}"
    version = ""
    
    def directory(self):
        return os.path.abspath(os.path.dirname(__file__))

    def unpack_downloaded_file(self, filename):
        print "unpacking", filename
        arguments = ['unzip']
        arguments.append(filename)
        subprocess.call(
            arguments, 
            cwd = self.directory()
        )
        print "done"

    def rename_unpacked_dir(self):
	dir = self.unpacked_dir_template.format(version = self.version)
	print "renaming unpacked dir", dir, "to 'deploy'"
	os.rename(dir, 'deploy')
        
    def start(self):
        if os.path.exists('deploy'):
            counter = 0
            while os.path.exists('deploy.backup.{0}'.format(counter)):
                counter += 1
                if counter > 100:
                    print "too many backup directories"
                    break;
            os.rename('deploy', 'deploy.backup.{0}'.format(counter))
        
        url = self.url_template.format(version = self.version)
        filename = self.filename_template.format(version = self.version)
        print "downloading Ibis-Deploy from", url, "to", filename
        urllib.urlretrieve(url, filename = filename)
        print "downloading finished"
        self.unpack_downloaded_file(filename)
	self.rename_unpacked_dir()
	
    
if __name__ == '__main__':
    instance = DownloadIbisDeployFromWebpage()
    if len(sys.argv) > 1:
        instance.version = sys.argv[1]
    instance.start()
