#!/usr/bin/env python

import subprocess
import os
import urllib.request
import urllib.parse
import urllib.error
from shutil import which
from optparse import OptionParser


class GetCodeFromHttp:
    filename_template = "{version}.tar.gz"
    name = "SYCLIST"
    url_template = "https://github.com/GESEG/SYCLIST/archive/{version}.tar.gz"
    version = ""

    def directory(self):
        return os.path.abspath(os.path.dirname(__file__))

    def src_directory(self):
        return os.path.join(self.directory(), 'src')

    def unpack_downloaded_file(self, filename, name, version):
        print("unpacking", filename)
        arguments = ['tar', '-xf']
        arguments.append(filename)
        subprocess.call(
            arguments,
            cwd=os.path.join(self.src_directory())
        )
        subprocess.call(
            [
                'mv', f'{name}-{version}',
                name
            ],
            cwd=os.path.join(self.src_directory())
        )
        print("done")

    def start(self):
        if os.path.exists('src'):
            counter = 0
            while os.path.exists(f'src.{counter}'):
                counter += 1
                if counter > 100:
                    print("too many backup directories")
                    break
            os.rename('src', f'src.{counter}')

        os.mkdir('src')

        url = self.url_template.format(version=self.version)
        filename = self.filename_template.format(version=self.version)
        filepath = os.path.join(self.src_directory(), filename)
        print(
            f"downloading version {self.version}"
            f"from {url} to {filename}"
        )
        if which('wget') is not None:
            arguments = ['wget', url]
            subprocess.call(
                arguments,
                cwd=os.path.join(self.src_directory())
            )
        elif which('curl') is not None:
            arguments = ['curl', '-L', '-O', url]
            subprocess.call(
                arguments,
                cwd=os.path.join(self.src_directory())
            )
        else:
            urllib.request.urlretrieve(url, filepath)
        print("downloading finished")
        self.unpack_downloaded_file(
            filename, self.name, self.version
        )


def main(version=''):
    instance = GetCodeFromHttp()
    instance.version = version
    instance.start()


def new_option_parser():
    result = OptionParser()
    result.add_option(
        "--version",
        default='a38f080e50f157b92684df61e678f97faf24d68e',
        dest="version",
        help="SYCLIST commit hash to download",
        type="string"
    )
    return result


if __name__ == "__main__":
    options, arguments = new_option_parser().parse_args()
    main(**options.__dict__)
