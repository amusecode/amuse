#!/usr/bin/env python

import subprocess
import os
import sys
import time
import urllib.request
import urllib.parse
import urllib.error
from optparse import OptionParser
import glob


class GetCodeFromHttp(object):
    url_template = ''
    filename_template = "mesa-r{version}.zip"
    version = ""
    zip=True

    def directory(self):
        return os.path.abspath(os.path.dirname(__file__))

    def src_directory(self):
        return os.path.join(self.directory(), 'src')

    def unpack_downloaded_file(self, filename):
        print("unpacking", filename)
        if self.zip:
            arguments = ['unzip']
        else:
            arguments = ['tar','xvf']
        arguments.append(filename)
        subprocess.call(
            arguments,
            cwd=os.path.join(self.src_directory())
        )
        print("done")

    def start(self):
        try:
            os.mkdir('src')
        except FileExistsError:
            pass

        url = self.url_template.format(version=self.version)
        filename = self.filename_template.format(version=self.version)
        filepath = os.path.join(self.src_directory(), filename)
        print("downloading version", self.version, "from", url, "to", filename)
        urllib.request.urlretrieve(url, filepath)
        print("downloading finished")
        self.unpack_downloaded_file(filename)


def get_crmath():
    instance = GetCodeFromHttp()
    instance.url_template = 'https://github.com/rhdtownsend/crmath/archive/{version}.zip'
    instance.filename_template='{version}.zip'
    instance.version = 'v1.2'
    instance.start()  

    subprocess.call(
            'make',
            cwd=os.path.join(instance.src_directory(),'crmath-1.2')
        )

def get_crlibm():
    instance = GetCodeFromHttp()
    instance.url_template = 'http://www.astro.wisc.edu/~townsend/resource/download/sdk2/src/crlibm-{version}.tar.gz'
    instance.filename_template='crlibm-{version}.tar.gz'
    instance.version = '1.0beta4'
    instance.zip=False
    instance.start()  


    wd = os.path.join(instance.src_directory(),'crlibm-1.0beta4')

    subprocess.call(
            ['./configure','--enable-static=yes'],
            cwd=wd
        )

    subprocess.call(
            'make',
            cwd=wd
        )

    subprocess.call(
            ['make','install'],
            cwd=wd
        )

    subprocess.call(
        ['ar','cru', 'libcrlibm.a'] + glob.glob(wd+'/*.o') + glob.glob(wd+'/scs_lib/*.o'), 
        cwd=wd
    )

def get_mesa():
    instance = GetCodeFromHttp()

    if 'AMUSE_LOCAL_COPY_MESA' in os.environ:
        instance.url_template = os.environ['AMUSE_LOCAL_COPY_MESA']
    else:
        instance.url_template = "https://zenodo.org/record/4311514/files/mesa-r{version}.zip"

    instance.version='15140'
    instance.start()   


def main():
    get_crmath()
    get_crlibm()
    get_mesa()



if __name__ == "__main__":
    main()
