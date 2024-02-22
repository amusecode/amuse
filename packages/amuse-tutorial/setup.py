#!/usr/bin/env python3
from support.classifiers import classifiers
from support.misc import find_data_files

from setuptools import setup

name = 'amuse-tutorial'
author = 'The AMUSE team'
author_email = 'info@amusecode.org'
license_ = "Apache License 2.0"
url = 'http://www.amusecode.org/'
install_requires = [
    'matplotlib>=2.2',
    'amuse-framework',
    'amuse-bhtree',
    'amuse-hermite',
    'amuse-seba',
    'amuse-sphray',
    'notebook',
]
description = 'The Astrophysical Multipurpose Software Environment - tutorial'
with open("README.md", "r") as fh:
    long_description = fh.read()
long_description_content_type = "text/markdown"

all_data_files = find_data_files(
    'tutorial', 'share/amuse/tutorial', '*', recursive=True
)

setup_requires = ['setuptools_scm']
use_scm_version = {
    "root": "../..",
    "relative_to": __file__,
    "version_file": "tutorial/_version.py",
}

setup(
    name=name,
    use_scm_version=use_scm_version,
    setup_requires=setup_requires,
    classifiers=classifiers,
    url=url,
    author_email=author_email,
    author=author,
    license=license_,
    description=description,
    long_description=long_description,
    long_description_content_type=long_description_content_type,
    install_requires=install_requires,
    python_requires=">=3.7",
    # cmdclass=mapping_from_command_name_to_command_class,
    data_files=all_data_files,
    scripts=["bin/amuse-tutorial"],
    packages=[],
)
