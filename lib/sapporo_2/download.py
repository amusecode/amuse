#!/usr/bin/env python

import subprocess
import os
import sys
import time
import urllib
import shutil
from optparse import OptionParser

        
        

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
    url = "https://github.com/treecode/sapporo2/tarball/master"
    alternative_url = "http://www.amusecode.org/codes/sapporo2-598e88c.tgz"
    filename = "master.tgz"
    
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
        
        for x in os.listdir(os.path.join(self.src_directory())):
            if x.startswith('treecode-sapporo2-'):
                subprocess.call(
                    ['mv',x, 'sapporo2-master'],
                    cwd = os.path.join(self.src_directory())
                )
                break
        print "done"
        
        
    def start(self):
        if os.path.exists('src'):
            counter = 0
            while os.path.exists('src.{0}'.format(counter)):
                counter += 1
                if counter > 5:
                    print "too many backup directories"
                    break;
            os.rename('src', 'src.{0}'.format(counter))
        
        os.mkdir('src')
        
        url = self.url
        filename = self.filename
        print "downloading master version from", url, "to", filename
        opener = MyFancyUrlopener()
        try:
            filename, httpmessage, code = opener.retrieve(url, filename = os.path.join(self.src_directory(),filename))
            if code == 404:
                os.remove(filename)
                url = self.alternative_url
                filename, httpmessage, code = opener.retrieve(url, filename = os.path.join(self.src_directory(),filename))
        except IOError as ex:
            url = self.alternative_url
            print "error in downloading, downloading master version from", url, "to", filename
            filename, httpmessage, code = opener.retrieve(url, filename = os.path.join(self.src_directory(),filename))
        print "downloading finished"
        self.unpack_downloaded_file(filename)
    
def main(must_download_from_github = False):
    if must_download_from_github:
        print "download using git is not supported yet, will download tarball instead"
        instance = GetCodeFromHttp()
    else:
        instance = GetCodeFromHttp()
        
    instance.start()

        
    
def new_option_parser():
    result = OptionParser()
    result.add_option(
        "-g", "--github", 
        default = False,
        dest="must_download_from_github",
        help="if given will download the code from the github repository using git",
        action="store_true"
    )
    return result
    
if __name__ == "__main__":
    options, arguments = new_option_parser().parse_args()
    main(**options.__dict__)
    

