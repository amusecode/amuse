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
    name = ["PeTar", "SDAR", "FDPS"]
    url_template = [
        "https://github.com/lwang-astro/PeTar/archive/{version}.tar.gz",
        "https://github.com/lwang-astro/SDAR/archive/{version}.tar.gz",
        "https://github.com/FDPS/FDPS/archive/{version}.tar.gz",
    ]
    version = [
        "",
        "",
        "",
    ]

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
                'mv', '{name}-{version}'.format(name=name, version=version),
                name
            ],
            cwd=os.path.join(self.src_directory())
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

        for i, url_template in enumerate(self.url_template):
            url = url_template.format(version=self.version[i])
            filename = self.filename_template.format(version=self.version[i])
            filepath = os.path.join(self.src_directory(), filename)
            print(
                "downloading version", self.version[i],
                "from", url, "to", filename
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
                filename, self.name[i], self.version[i]
            )


def main(petar_version='', sdar_version='', fdps_version=''):
    version = [
        petar_version,
        sdar_version,
        fdps_version,
    ]
    instance = GetCodeFromHttp()
    instance.version = version
    instance.start()


def new_option_parser():
    result = OptionParser()
    result.add_option(
        "--petar-version",
        default='1b7c279fbb94af558952c56f4fd746fcf730ce83',
        dest="petar_version",
        help="PeTar commit hash to download",
        type="string"
    )
    result.add_option(
        "--sdar-version",
        default='2d0f79f259ffd6dd972f27548698d91e7cd476d7',
        dest="sdar_version",
        help="SDAR commit hash to download",
        type="string"
    )
    result.add_option(
        "--fdps-version",
        default='3370a770cb444eab106c81a6ebbc4cb9ef5226d6',
        dest="fdps_version",
        help="FDPS commit hash to download",
        type="string"
    )
    return result


if __name__ == "__main__":
    options, arguments = new_option_parser().parse_args()
    main(**options.__dict__)
