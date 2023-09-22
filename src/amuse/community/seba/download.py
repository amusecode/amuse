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
    name = ["SeBa"]
    url_template = [
        "https://github.com/amusecode/SeBa/archive/{version}.tar.gz",
    ]
    version = [
        "",
    ]

    def directory(self):
        return os.path.abspath(os.path.dirname(__file__))

    def src_directory(self):
        return os.path.join(self.directory(), 'src')

    def unpack_downloaded_file(self, filename, name, version):
        print(f"unpacking {filename}")
        arguments = ['tar', '-xf']
        arguments.append(filename)
        subprocess.call(
            arguments,
            cwd=os.path.join(self.src_directory())
        )
        subprocess.call(
            [
                'mv',
                f'{name}-{version}',
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

        for i, url_template in enumerate(self.url_template):
            url = url_template.format(version=self.version[i])
            filename = self.filename_template.format(version=self.version[i])
            filepath = os.path.join(self.src_directory(), filename)
            print(f"downloading version {self.version[i]} from {url} to {filename}")
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
                filename, self.name[i], self.version[i]
            )


def main(seba_version=''):
    version = [
        seba_version,
    ]
    instance = GetCodeFromHttp()
    instance.version = version
    instance.start()


def new_option_parser():
    result = OptionParser()
    result.add_option(
        "--seba-version",
        default='2f6e7f37a53167b4b0dcd6c723dff7b5ee1aecba',
        dest="seba_version",
        help="SeBa commit hash to download",
        type="string"
    )
    return result


if __name__ == "__main__":
    options, arguments = new_option_parser().parse_args()
    main(**options.__dict__)
