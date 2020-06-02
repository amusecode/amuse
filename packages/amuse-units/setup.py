import sys
import os
from support.version import version, main_version
from support.classifiers import classifiers

from setuptools import setup

name = 'amuse-units'
author = 'The AMUSE team'
author_email = 'info@amusecode.org'
license_ = "Apache License 2.0"
url = 'http://www.amusecode.org/'
install_requires = [
]
description = 'The Astrophysical Multipurpose Software Environment - unit system'
with open("README.md", "r") as fh:
    long_description = fh.read()
long_description_content_type = "text/markdown"

extensions = []

all_data_files = []

packages = [
    'amuse.units',
]

package_data = {
}

setup(
    name=name,
    version=version,
    classifiers=classifiers,
    url=url,
    author_email=author_email,
    author=author,
    license=license_,
    description=description,
    long_description=long_description,
    long_description_content_type=long_description_content_type,
    install_requires=install_requires,
    python_requires=">=3.5",
    ext_modules=extensions,
    package_dir={
        'amuse.units': 'src/amuse/units'
    },
    packages=packages,
    package_data=package_data,
    data_files=all_data_files,
)
