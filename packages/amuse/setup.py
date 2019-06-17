from setuptools import setup


name = 'amuse'
major_version = 12
minor_version = 0
patch_version = 0
additional_version = ".rc4"  # e.g. ".rc1"
post_version = ""  # e.g. ".post1"
version = "%i.%i.%i%s%s" % (
    major_version, minor_version, patch_version, additional_version,
    post_version
)
main_version = "%i.%i.0.rc2" % (major_version, minor_version)
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
    'amuse-framework>=%s' % main_version,
    'amuse-athena>=%s' % main_version,
    'amuse-bhtree>=%s' % main_version,
    'amuse-brutus>=%s' % main_version,
    'amuse-bse>=%s' % main_version,
    'amuse-capreole>=%s' % main_version,
    'amuse-evtwin>=%s' % main_version,
    'amuse-fastkick>=%s' % main_version,
    'amuse-fi>=%s' % main_version,
    'amuse-fractalcluster>=%s' % main_version,
    'amuse-framework>=%s' % version,
    'amuse-gadget2>=%s' % main_version,
    'amuse-galactics>=%s' % main_version,
    'amuse-galaxia>=%s' % main_version,
    'amuse-halogen>=%s' % main_version,
    'amuse-hermite>=%s' % main_version,
    'amuse-hop>=%s' % main_version,
    'amuse-huayno>=%s' % main_version,
    'amuse-kepler>=%s' % main_version,
    'amuse-kepler-orbiters>=%s' % main_version,
    'amuse-mameclot>=%s' % main_version,
    'amuse-mercury>=%s' % main_version,
    'amuse-mmams>=%s' % main_version,
    'amuse-ph4>=%s' % main_version,
    'amuse-seba>=%s' % main_version,
    'amuse-secularmultiple>=%s' % main_version,
    'amuse-simplex>=%s' % main_version,
    'amuse-smalln>=%s' % main_version,
    'amuse-sphray>=%s' % main_version,
    'amuse-sse>=%s' % main_version,
    'amuse-twobody>=%s' % main_version,
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
