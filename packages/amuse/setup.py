from setuptools import setup


name = 'amuse'
version = "12.0.0rc3"
author = 'The AMUSE team'
author_email = 'info@amusecode.org'
license_ = "Apache License 2.0"
url = 'http://www.amusecode.org/'
install_requires = [
    'wheel>=0.32',
    'docutils>=0.6',
    'numpy>=1.2.2',
    'nose>=0.11.1',
    'mpi4py>=1.1.0',
    'h5py>=1.1.0',
    'amuse-framework>=%s' % version,
    'amuse-athena>=%s' % version,
    'amuse-bhtree>=%s' % version,
    'amuse-brutus>=%s' % version,
    'amuse-bse>=%s' % version,
    'amuse-capreole>=%s' % version,
    'amuse-evtwin>=%s' % version,
    'amuse-fastkick>=%s' % version,
    'amuse-fi>=%s' % version,
    'amuse-fractalcluster>=%s' % version,
    'amuse-framework>=%s' % version,
    'amuse-gadget2>=%s' % version,
    'amuse-galactics>=%s' % version,
    'amuse-galaxia>=%s' % version,
    'amuse-halogen>=%s' % version,
    'amuse-hermite>=%s' % version,
    'amuse-hop>=%s' % version,
    'amuse-huayno>=%s' % version,
    'amuse-kepler>=%s' % version,
    'amuse-kepler-orbiters>=%s' % version,
    'amuse-mameclot>=%s' % version,
    'amuse-mercury>=%s' % version,
    'amuse-mmams>=%s' % version,
    'amuse-ph4>=%s' % version,
    'amuse-seba>=%s' % version,
    'amuse-secularmultiple>=%s' % version,
    'amuse-simplex>=%s' % version,
    'amuse-smalln>=%s' % version,
    'amuse-sphray>=%s' % version,
    'amuse-sse>=%s' % version,
    'amuse-twobody>=%s' % version,
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
    'Programming Language :: Python :: 2.6',
    'Programming Language :: Python :: 2.7',
    'Programming Language :: Python :: 3',
    'Programming Language :: C',
    'Programming Language :: C++',
    'Programming Language :: Fortran',
    'Topic :: Scientific/Engineering :: Astronomy',
]

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
)
