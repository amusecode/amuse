#!/usr/bin/env python3
from support.classifiers import classifiers

from setuptools import setup

import support
support.use("system")
from support.setup_codes import setup_commands

name = 'amuse-tests'
author = 'The AMUSE team'
author_email = 'info@amusecode.org'
license_ = "Apache License 2.0"
url = 'http://www.amusecode.org/'
install_requires = [
    'amuse-framework',
]
description = 'The Astrophysical Multipurpose Software Environment - tests'
with open("README.md", "r") as fh:
    long_description = fh.read()
long_description_content_type = "text/markdown"

extensions = []

all_data_files = []

packages = [
    'amuse.test.suite',
    'amuse.test.suite.ext_tests', 
    'amuse.test.suite.core_tests',
    'amuse.test.suite.compile_tests', 
    'amuse.test.suite.codes_tests',
    'amuse.test.suite.ticket_tests', 
    'amuse.test.suite.reports',
]

package_data = {
}

mapping_from_command_name_to_command_class = setup_commands()

try:
    from src.amuse.test.suite.version import version
    use_scm_version = False
    setup_requires = []
except ImportError:
    version = False
    setup_requires = ['setuptools_scm']
    use_scm_version = {
        "root": "../..",
        "relative_to": __file__,
        "write_to": "src/amuse/test/suite/version.py",
    }

setup(
    name=name,
    use_scm_version=use_scm_version,
    setup_requires=setup_requires,
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
    cmdclass=mapping_from_command_name_to_command_class,
    ext_modules=extensions,
    package_dir={
        'amuse.test.suite': 'src/amuse/test/suite',
    },
    packages=packages,
    package_data=package_data,
    data_files=all_data_files,
)
