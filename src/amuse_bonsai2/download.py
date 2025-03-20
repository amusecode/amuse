#!/usr/bin/env python

import subprocess
import os
import sys
import time
import urllib.request
import urllib.parse
import urllib.error
from optparse import OptionParser


class GetCodeFromHttp(object):
    url_template = "https://github.com/treecode/Bonsai/archive/{version}.zip"
    filename_template = "{version}.zip"
    version = ""

    def directory(self):
        return os.path.abspath(os.path.dirname(__file__))

    def src_directory(self):
        return os.path.join(self.directory(), 'src')

    def unpack_downloaded_file(self, filename):
        print("unpacking", filename)
        arguments = ['tar', '-xf']
        arguments.append(filename)
        subprocess.call(
            arguments,
            cwd=os.path.join(self.src_directory())
        )
        subprocess.call(
            ['mv', 'bonsai-{version}'.format(version = self.version), 'bonsai'],
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
                    break
            os.rename('src', 'src.{0}'.format(counter))

        os.mkdir('src')

        url = self.url_template.format(version=self.version)
        filename = self.filename_template.format(version=self.version)
        filepath = os.path.join(self.src_directory(), filename)
        print("downloading version", self.version, "from", url, "to", filename)
        urllib.request.urlretrieve(url, filepath)
        print("downloading finished")
        self.unpack_downloaded_file(filename)


def main(version=''):
    instance = GetCodeFromHttp()
    instance.version = version
    instance.start()


def new_option_parser():
    result = OptionParser()
    result.add_option(
        "--version",
        default='d1c23d33263090fcf681aa727f300b31efef0f03',
        dest="version",
        help="version number to download",
        type="string"
    )
    return result


if __name__ == "__main__":
    options, arguments = new_option_parser().parse_args()
    main(**options.__dict__)
