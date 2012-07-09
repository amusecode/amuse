#!/usr/bin/env python

import subprocess
import os
import sys
import time
import urllib

class DownloadIbisDeployFromWebpage(object):
    url_template = "http://www.cs.vu.nl/~niels/downloads/deploy{version}.zip"
    filename_template = "deploy{version}.zip"
    unpacked_dir_template = "deploy{version}"
    version = ""
    
    def directory(self):
        return os.path.abspath(os.path.dirname(__file__))

    def unpack_downloaded_file(self, filename):
        print "unpacking", filename
        arguments = ['unzip']
        arguments.append(filename)
        result = subprocess.call(
            arguments, 
            cwd = self.directory()
        )
	if result != 0:
	    print "error on unpacking downloaded file"
	    sys.exit(1)
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
        zipfilename = self.filename_template.format(version = self.version)
        print "downloading Ibis-Deploy from", url, "to", zipfilename
        urllib.urlretrieve(url, filename = zipfilename)
        print "downloading finished"
	if not os.path.exists(zipfilename):
		print "could not download file " + zipfilename
		sys.exit(1)

        self.unpack_downloaded_file(zipfilename)
	os.remove(zipfilename)
    
if __name__ == '__main__':
    instance = DownloadIbisDeployFromWebpage()
    if len(sys.argv) > 1:
        instance.version = sys.argv[1]
    instance.start()
