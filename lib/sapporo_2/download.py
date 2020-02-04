#!/usr/bin/env python

import subprocess
import os
import urllib.request
import urllib.parse
import urllib.error
from optparse import OptionParser


class GetCodeFromHttp(object):
    url = "https://github.com/treecode/sapporo2/tarball/master"
    alternative_url = \
        "http://amuse.strw.leidenuniv.nl/codes/sapporo2-598e88c.tgz"
    filename = "master.tgz"

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

        for x in os.listdir(os.path.join(self.src_directory())):
            if x.startswith('treecode-sapporo2-'):
                subprocess.call(
                    ['mv', x, 'sapporo2-master'],
                    cwd=os.path.join(self.src_directory())
                )
                break
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

        url = self.url
        filename = self.filename
        filepath = os.path.join(self.src_directory(), filename)
        print("downloading sapporo2 from", url, "to", filename)
        urllib.request.urlretrieve(url, filepath)
        print("downloading finished")
        self.unpack_downloaded_file(filename)


def main(must_download_from_github=False):
    if must_download_from_github:
        print(
            "download using git is not supported yet, will download tarball "
            "instead"
        )
        instance = GetCodeFromHttp()
    else:
        instance = GetCodeFromHttp()

    instance.start()


def new_option_parser():
    result = OptionParser()
    result.add_option(
        "-g", "--github",
        default=False,
        dest="must_download_from_github",
        help="if given will download the code from the github repository "
        "using git",
        action="store_true"
    )
    return result


if __name__ == "__main__":
    options, arguments = new_option_parser().parse_args()
    main(**options.__dict__)
