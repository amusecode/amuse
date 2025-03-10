#!/usr/bin/env python

import subprocess
import os
import argparse
import urllib.request
import urllib.parse
import urllib.error
from shutil import which


class GetCodeFromHttp:
    filename_template = "{version}.tar.gz"
    name = ["rebound"]
    url_template = [
        "https://github.com/hannorein/rebound/archive/{version}.tar.gz"
    ]
    version = [
        "",
    ]

    def directory(self):
        return os.path.abspath(os.path.dirname(__file__))

    def src_directory(self):
        return os.path.join(self.directory(), "src")

    def unpack_downloaded_file(self, filename, name, version):
        print(f"unpacking {filename}")
        arguments = ["tar", "-xf"]
        arguments.append(filename)
        subprocess.call(arguments, cwd=os.path.join(self.src_directory()))
        subprocess.call(
            ["mv", f"{name}-{version}", name],
            cwd=os.path.join(self.src_directory()),
        )
        print("done")

    def start(self):
        if os.path.exists("src"):
            counter = 0
            while os.path.exists(f"src.{counter}"):
                counter += 1
                if counter > 100:
                    print("too many backup directories")
                    break
            os.rename("src", f"src.{counter}")

        os.mkdir("src")

        for i, url_template in enumerate(self.url_template):
            url = url_template.format(version=self.version[i])
            filename = self.filename_template.format(version=self.version[i])
            filepath = os.path.join(self.src_directory(), filename)
            print(f"downloading version {self.version[i]} from {url} to {filename}")
            if which("wget") is not None:
                arguments = ["wget", url]
                subprocess.call(arguments, cwd=os.path.join(self.src_directory()))
            elif which("curl") is not None:
                arguments = ["curl", "-L", "-O", url]
                subprocess.call(arguments, cwd=os.path.join(self.src_directory()))
            else:
                urllib.request.urlretrieve(url, filepath)
            print("downloading finished")
            self.unpack_downloaded_file(filename, self.name[i], self.version[i])


def new_argument_parser():
    result = argparse.ArgumentParser()
    result.add_argument(
        "--rebound-version",
        default="111bcd50de54e4cc70e1d4915918226a4adbe188",
        # default="4df2b4603058e7d767aee9c26b586b55a5f4de65",
        dest="rebound_version",
        help="Rebound commit hash to download",
        type=str,
    )
    return result


def main(version=""):
    arguments = new_argument_parser().parse_args()
    version = [
        arguments.rebound_version,
    ]
    instance = GetCodeFromHttp()
    instance.version = version
    instance.start()


if __name__ == "__main__":
    main()
