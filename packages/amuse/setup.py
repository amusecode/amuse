from setuptools import setup
from support.classifiers import classifiers
from setuptools_scm import get_version

version = get_version(
    root='../..',
    relative_to=__file__,
)

name = 'amuse'
author = 'The AMUSE team'
author_email = 'info@amusecode.org'
license_ = "Apache License 2.0"
url = 'http://www.amusecode.org/'
install_requires = [
    'matplotlib>=2.2',
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
    'amuse-phigrape>=%s' % version,
    'amuse-seba>=%s' % version,
    'amuse-secularmultiple>=%s' % version,
# until C++ MPI code is replaced/fixed
#    'amuse-simplex>=%s' % version,
    'amuse-smalln>=%s' % version,
    'amuse-sphray>=%s' % version,
    'amuse-sse>=%s' % version,
    'amuse-twobody>=%s' % version,
]
description = 'The Astrophysical Multipurpose Software Environment'
with open("README.md", "r") as fh:
    long_description = fh.read()
long_description_content_type = "text/markdown"


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
)
