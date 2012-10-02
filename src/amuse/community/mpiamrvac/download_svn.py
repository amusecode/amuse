#!/usr/bin/env python

import subprocess
import os
import sys
import time
import urllib
from optparse import OptionParser

# mpiamrvac revisions in amuse
# 145
# 187

class GetCodeFromSVN(object):
    revision = '187'
    username = 'studCMFA09'
    password = 'cpa9amrvac'
    url = 'https://svn.esat.kuleuven.be/amrvac/trunk'
    
    def directory(self):
        return os.path.abspath(os.path.dirname(__file__))
    
    def source_directory(self):
        return os.path.join(self.directory(), 'src')
        
    def code_directory(self):
        return os.path.join(self.source_directory(), 'mpiamrvac')
        
    def start(self):
        arguments = [
            'svn',
            'co',
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
        
        

class MyFancyUrlopener(urllib.FancyURLopener):
    def retrieve(self, url, filename=None, reporthook=None, data=None):
        """retrieve(url) returns (filename, headers) for a local object
        or (tempfilename, headers) for a remote object."""
        url = urllib.unwrap(urllib.toBytes(url))
        if self.tempcache and url in self.tempcache:
            return self.tempcache[url]
        type, url1 = urllib.splittype(url)
        if filename is None and (not type or type == 'file'):
            try:
                fp = self.open_local_file(url1)
                hdrs = fp.info()
                del fp
                return urllib.url2pathname(urllib.splithost(url1)[1]), hdrs
            except IOError, msg:
                pass
        fp = self.open(url, data)
        try:
            headers = fp.info()
            code = fp.code
            if filename:
                tfp = open(filename, 'wb')
            else:
                import tempfile
                garbage, path = urllib.splittype(url)
                garbage, path = urllib.splithost(path or "")
                path, garbage = urllib.splitquery(path or "")
                path, garbage = urllib.splitattr(path or "")
                suffix = os.path.splitext(path)[1]
                (fd, filename) = tempfile.mkstemp(suffix)
                self.__tempfiles.append(filename)
                tfp = os.fdopen(fd, 'wb')
            try:
                result = filename, headers, code
                if self.tempcache is not None:
                    self.tempcache[url] = result
                bs = 1024*8
                size = -1
                read = 0
                blocknum = 0
                if reporthook:
                    if "content-length" in headers:
                        size = int(headers["Content-Length"])
                    reporthook(blocknum, bs, size)
                while 1:
                    block = fp.read(bs)
                    if block == "":
                        break
                    read += len(block)
                    tfp.write(block)
                    blocknum += 1
                    if reporthook:
                        reporthook(blocknum, bs, size)
            finally:
                tfp.close()
        finally:
            fp.close()
        del fp
        del tfp

        # raise exception if actual size does not match content-length header
        if size >= 0 and read < size:
            raise urllib.ContentTooShortError("retrieval incomplete: got only %i out "
                                       "of %i bytes" % (read, size), result)

        return result
    
class GetCodeFromHttp(object):
    url_template = "http://www.amusecode.org/codes/mpiamrvac-r{version}.tgz"
    filename_template = "mpiamrvac-r{version}.tgz"
    version = "187"
    
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
        print "downloading version",self.version,"from", url, "to", filename
        opener = MyFancyUrlopener()
        filename, httpmessage, code = opener.retrieve(url, filename = os.path.join(self.src_directory(),filename))
        if code == 404:
            os.remove(filename)
            url = self.backup_url_template.format(version = self.version)
            filename, httpmessage, code = opener.retrieve(url, filename = os.path.join(self.src_directory(),filename))
            
        print "downloading finished"
        self.unpack_downloaded_file(filename)
    
def main(must_download_from_svn = False, version = '149'):
    if must_download_from_svn:
        instance = GetCodeFromSVN()
        instance.revision = version
    else:
        instance = GetCodeFromHttp()
        instance.version = version
        
    instance.start()

        
    
def new_option_parser():
    result = OptionParser()
    result.add_option(
        "-s", "--svn", 
        default = False,
        dest="must_download_from_svn",
        help="if given will download the code from the svn repository",
        action="store_true"
    )
    result.add_option(
        "--version", 
        default = '187',
        dest="version",
        help="svn version number to download",
        type="string"
    )
    
    return result
    
if __name__ == "__main__":
    options, arguments = new_option_parser().parse_args()
    main(**options.__dict__)
    

