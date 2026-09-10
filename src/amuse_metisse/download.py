"""
Downloads sample MESA tracks for METISSE.

This should normally be run during installation, to a fixed directory.
"""

import os
import argparse
import urllib.request


def download(url, filename):
    if not os.path.exists(filename):
        urllib.request.urlretrieve(url, filename)


def new_argument_parser():
    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument("--directory", default="./data", help="Directory for output files")
    parser.add_argument("--verbose", "-v", action="count", default=0)
    return parser


def download_sample_metisse_tracks():
    parser = new_argument_parser()
    arguments = parser.parse_args()
    directory = arguments.directory
    os.makedirs(directory, exist_ok=True)
    download(
       "https://zenodo.org/records/17513335/files/sample_tracks_solarZ.zip?download=1",
       os.path.join(directory, "sample_tracks_solarZ.zip"),
    )
    unzip = f"unzip {os.path.join(directory, 'sample_tracks_solarZ.zip')} -d {directory}"
    os.system(unzip)
    os.remove(os.path.join(directory, "sample_tracks_solarZ.zip"))


def main():
    download_sample_metisse_tracks()


if __name__ == "__main__":
    main()
