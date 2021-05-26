from support.classifiers import classifiers

from setuptools import setup, find_packages
from support.setup_codes import setup_commands

name = 'amuse-framework'
author = 'The AMUSE team'
author_email = 'info@amusecode.org'
license_ = "Apache License 2.0"
url = 'http://www.amusecode.org/'
install_requires = [
    'setuptools>=41.0.0',
    'pip>=19.0.0',
    'wheel>=0.32',
    'docutils>=0.6',
    'numpy>=1.2.2',
    'pytest>=4.0',
    'h5py>=1.1.0',
]
description = 'The Astrophysical Multipurpose Software Environment'
with open("README.md", "r") as fh:
    long_description = fh.read()
long_description_content_type = "text/markdown"

extensions = []

all_data_files = []
all_data_files.append(('share/amuse', ['./config.mk', './build.py']))

packages = find_packages('src', exclude=["amuse.community.*"])
packages.append("amuse.community.interface")

package_data = {
    'amuse.rfi.tools': ['*.template'],
    'amuse': [
        '*rc'
    ]
}

mapping_from_command_name_to_command_class=setup_commands()

try:
    from src.amuse.version import version
    use_scm_version = False
    setup_requires = []
except ImportError:
    version = False
    setup_requires = ['setuptools_scm']
    use_scm_version = {
        "root": "../..",
        "relative_to": __file__,
        "write_to": "src/amuse/version.py",
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
    extras_require = {
        "MPI" : ["mpi4py>=1.1.0"]
    },
    cmdclass=mapping_from_command_name_to_command_class,
    ext_modules=extensions,
    package_dir={'': 'src'},
    packages=packages,
    package_data=package_data,
    data_files=all_data_files,
    scripts=[ "bin/amusifier" ],
)
