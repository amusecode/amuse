from setuptools import setup, find_packages
from support.setup_codes import setup_commands
from support.misc import find_data_files

name = 'amuse-devel'
author = 'The AMUSE team'
author_email = 'info@amusecode.org'
license_ = "Apache License 2.0"
url = 'http://www.amusecode.org/'
install_requires = [
    'setuptools>=41.0.0',
    'setuptools_scm',
    'pip>=19.0.0',
    'wheel>=0.32',
    'docutils>=0.6',
    'numpy>=1.2.2',
    'pytest>=4.0',
    # 'mpi4py>=1.1.0',
    'h5py>=1.1.0',
]
description = 'The Astrophysical Multipurpose Software Environment'
with open("README.md", "r") as fh:
    long_description = fh.read()
long_description_content_type = "text/markdown"
classifiers = [
    'Development Status :: 4 - Beta',
    'Environment :: Console',
    'Intended Audience :: End Users/Desktop',
    'Intended Audience :: Developers',
    'Intended Audience :: Science/Research',
    'License :: OSI Approved :: Apache Software License',
    'Operating System :: MacOS :: MacOS X',
    'Operating System :: POSIX',
    'Programming Language :: Python :: 3',
    'Programming Language :: Python :: 3.7',
    'Programming Language :: Python :: 3.8',
    'Programming Language :: Python :: 3.9',
    'Programming Language :: Python :: 3.10',
    'Programming Language :: C',
    'Programming Language :: C++',
    'Programming Language :: Fortran',
    'Topic :: Scientific/Engineering :: Astronomy',
]

extensions = []

all_data_files = []
all_data_files.append(('share/amuse', ['./config.mk']))

packages = find_packages('src')
packages.extend(
    ['amuse.examples.' + x for x in find_packages('examples')]
)

package_data = {
    'amuse.rfi.tools': ['*.template'],
    'amuse.test.suite.core_tests': [
        '*.txt', '*.dyn', '*.ini',
        '*.nemo',
        '*.dat', 'gadget_snapshot'
    ],
    'amuse.test.suite.codes_tests': [
        '*.txt', 'test_sphray_data*'
    ],
    'amuse.test.suite.ticket_tests': [
        '*.out'
    ],
    'amuse': [
        '*rc'
    ]
}

mapping_from_command_name_to_command_class = setup_commands()

setup(
    name=name,
    use_scm_version={
        "write_to": "src/amuse/version.py",
    },
    setup_requires=['setuptools_scm'],
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
    extras_require = {
        "MPI" : ["mpi4py>=1.1.0"]
    },
    cmdclass=mapping_from_command_name_to_command_class,
    ext_modules=extensions,
    package_dir={'': 'src', 'amuse.examples': 'examples'},
    packages=packages,
    package_data=package_data,
    data_files=all_data_files,
    scripts=["bin/amusifier", "bin/amuse-tutorial", ],
)
