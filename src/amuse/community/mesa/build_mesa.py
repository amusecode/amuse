#!/usr/bin/env python

import subprocess
import os
import sys
import time
import urllib

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
    url_template = "http://www.amusecode.org/codes/mesa-r{version}.tgz"
    filename_template = "mesa-r{version}.tgz"
    version = "2208"
    
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
        
class BuildMesa(object):
    
    def mesa_directory(self):
        return os.path.abspath(os.path.dirname(__file__))
        
        
    def get_mesa_source_from_svn(self):
        mesa_url = "http://mesa.svn.sourceforge.net/svnroot/mesa/trunk"
        revision = "2208"
        subprocess.call(['svn', 'export', 
            '-r', revision, mesa_url, os.path.join('src','mesa')], cwd = self.mesa_directory())
        self.patch_mesa_source()
        
    def get_mesa_source_from_http(self):
        revision = "2208"
        instance = GetCodeFromHttp()
        instance.version = "2208"
        instance.start()
        self.patch_mesa_source()
        
    def patch_mesa_source(self):
        subprocess.call(['cp','-f','./mesa_reqs/makefile_header',
            './src/mesa/utils/makefile_header'], cwd = self.mesa_directory())
        subprocess.call(['cp','-f','./mesa_reqs/create_zams.f',
            './src/mesa/star/test/src/create_zams.f'], cwd = self.mesa_directory())
        subprocess.call(['cp','-f','./mesa_reqs/eos_def.f',
            './src/mesa/eos/public/eos_def.f'], cwd = self.mesa_directory())
        subprocess.call(['cp','-f','./mesa_reqs/jina_def.f',
            './src/mesa/jina/public/jina_def.f'], cwd = self.mesa_directory())
        subprocess.call(['cp','-f','./mesa_reqs/kap_def.f',
            './src/mesa/kap/public/kap_def.f'], cwd = self.mesa_directory())
        subprocess.call(['cp','-f','./mesa_reqs/net_def.f',
            './src/mesa/net/public/net_def.f'], cwd = self.mesa_directory())
        subprocess.call(['cp','-f','./mesa_reqs/star_def.f',
            './src/mesa/star/public/star_def.f'], cwd = self.mesa_directory())
        
    def build_mesa(self):
        self.get_mesa_source_from_http()
        subprocess.call(['./install'], cwd = os.path.join(self.mesa_directory(), 'src', 'mesa'))

    def clean_mesa(self):
        subprocess.call(['./clean'], cwd = os.path.join(self.mesa_directory(), 'src', 'mesa'))
        print "Finished cleaning.\n"

    def very_clean_mesa(self):
        self.clean_mesa()
        print "Also emptying the EOS, KAP, and NET caches."
        subprocess.call(['./empty_caches'], cwd = os.path.join(self.mesa_directory(), 'src', 'mesa'))

if __name__ == '__main__':
    instance = BuildMesa()
    if len(sys.argv) == 1:
        instance.build_mesa()
    elif sys.argv[1] == "download":
        instance.get_mesa_source_from_http()
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
