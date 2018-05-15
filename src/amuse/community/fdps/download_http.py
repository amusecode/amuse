#!/usr/bin/env python

import subprocess
import os
import sys
import time
if sys.version_info <= (3,0):
    from urllib import (
            FancyURLopener, splittype, splithost, splitquery, splitattr,
            unwrap, splittype, url2pathname, ContentTooShortError,
            toBytes
            )
else:
    from urllib.request import (
            FancyURLopener, splittype, splithost, splitquery, splitattr,
            unwrap, splittype, url2pathname, ContentTooShortError,
            urlretrieve
            )
from optparse import OptionParser

        
        

class MyFancyUrlopener(FancyURLopener):
    def retrieve(self, url, filename=None, reporthook=None, data=None):
        """retrieve(url) returns (filename, headers) for a local object
        or (tempfilename, headers) for a remote object."""
        url = (
                unwrap(toBytes(url)) if sys.version_info <= (3,0)
                else unwrap(url)
                )
        if self.tempcache and url in self.tempcache:
            return self.tempcache[url]
        type, url1 = splittype(url)
        if filename is None and (not type or type == 'file'):
            try:
                fp = self.open_local_file(url1)
                hdrs = fp.info()
                del fp
                return url2pathname(splithost(url1)[1]), hdrs
            except IOError as msg:
                pass
        fp = self.open(url, data)
        try:
            headers = fp.info()
            code = fp.code
            if filename:
                tfp = open(filename, 'wb')
            else:
                import tempfile
                garbage, path = splittype(url)
                garbage, path = splithost(path or "")
                path, garbage = splitquery(path or "")
                path, garbage = splitattr(path or "")
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
            raise ContentTooShortError("retrieval incomplete: got only %i out "
                                       "of %i bytes" % (read, size), result)

        return result


class GetCodeFromHttp(object):
    url_template = "https://github.com/FDPS/FDPS/archive/{version}.zip"
    filename_template = "v{version}.zip"
    version = "4.0b"
    
    def directory(self):
        return os.path.abspath(os.path.dirname(__file__))

    def src_directory(self):
        return os.path.join(self.directory(), 'src')
        
    def unpack_downloaded_file(self, filename):
        print("unpacking", filename)
        arguments = ['unzip', '-x']
        arguments.append(filename)
        subprocess.call(
            arguments, 
            cwd = os.path.join(self.src_directory())
        )
        subprocess.call(
            ['mv', 'FDPS-{version}'.format(version = self.version), 'FDPS'],
            cwd = os.path.join(self.src_directory())
        )
        print("done")
        
        
    def start(self):
        if os.path.exists('src'):
            counter = 0
            while os.path.exists('src.{0}'.format(counter)):
                counter += 1
                if counter > 100:
                    print("too many backup directories")
                    break;
            os.rename('src', 'src.{0}'.format(counter))
        
        os.mkdir('src')
        
        url = self.url_template.format(version = self.version)
        filename = self.filename_template.format(version = self.version)
        print("downloading version",self.version,"from", url, "to", filename)
        if sys.version_info <= (3,0):
            opener = MyFancyUrlopener()
            filename, httpmessage, code = opener.retrieve(url, filename = os.path.join(self.src_directory(),filename))
            if code == 404:
                os.remove(filename)
                url = self.backup_url_template.format(version = self.version)
                filename, httpmessage, code = opener.retrieve(url, filename = os.path.join(self.src_directory(),filename))
        else:
            filename, headers = urlretrieve(url, os.path.join(self.src_directory(), filename))
                
        print("downloading finished")
        self.unpack_downloaded_file(filename)
    
def main(must_download_from_svn = False, version = '2.0e'):
    instance = GetCodeFromHttp()
    instance.version = version
    instance.start()

        
    
def new_option_parser():
    result = OptionParser()
    
    result.add_option(
        "--version", 
        default = '4.0b',
        dest="version",
        help="release to download from github",
        type="string"
    )
    
    return result
    
if __name__ == "__main__":
    options, arguments = new_option_parser().parse_args()
    main(**options.__dict__)
    

