from setuptools import setup
from support.version import version, main_version
from support.classifiers import classifiers

name = 'amuse'
author = 'The AMUSE team'
author_email = 'info@amusecode.org'
license_ = "Apache License 2.0"
url = 'http://www.amusecode.org/'
install_requires = [
    'matplotlib>=2.2',
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
    'amuse-phigrape>=%s' % main_version,
    'amuse-seba>=%s' % main_version,
    'amuse-secularmultiple>=%s' % main_version,
# until C++ MPI code is replaced/fixed
#    'amuse-simplex>=%s' % main_version,
    'amuse-smalln>=%s' % main_version,
    'amuse-sphray>=%s' % main_version,
    'amuse-sse>=%s' % main_version,
    'amuse-twobody>=%s' % main_version,
]
description = 'The Astrophysical Multipurpose Software Environment'
with open("README.md", "r") as fh:
    long_description = fh.read()
long_description_content_type = "text/markdown"

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
)
