from support.version import version, main_version
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
    'amuse-framework>=%s' % main_version,
    'amuse-bhtree>=%s' % main_version,
    'amuse-hermite>=%s' % main_version,
    'amuse-seba>=%s' % main_version,
    'amuse-sphray>=%s' % main_version,
    'notebook',
]
description = 'The Astrophysical Multipurpose Software Environment: tutorial'
with open("README.md", "r") as fh:
    long_description = fh.read()
long_description_content_type = "text/markdown"

all_data_files = find_data_files(
    'tutorial', 'share/amuse/tutorial', '*', recursive=True
)

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
    # cmdclass=mapping_from_command_name_to_command_class,
    data_files=all_data_files,
    scripts=["bin/amuse-tutorial"],
)
