"""Version packaging helper

We use versioningit to automatically get a version for the packages we build from the
current git repository. When building from a tarball however, which is what we instruct
users to do, there is no git history for versioningit to inspect, and it will give an
error. There's a method that has git write the version to the pyproject.toml file when
creating the tarball using git archive, but that only works for a single pyproject.toml
file and we have many.

So, we have git write the version to a VERSION file in the root instead, and then we
have this plug-in for versioningit that checks that if we don't have a git repository.
"""
from versioningit.core import VCSDescription
from versioningit.errors import NotSdistError, NotVCSError
from versioningit.git import describe_git

from pathlib import Path


def load_version(version_file: Path):
    version = None
    if version_file.exists():
        with version_file.open('r') as f:
            for line in f:
                if line.strip().startswith('#') or line.strip() == '':
                    continue
                version = line.strip()

    return version


def get_amuse_version(project_dir, params):
    """Find the current version of AMUSE"""
    try:
        return describe_git(project_dir=project_dir, params=params)
    except (NotVCSError, NotSdistError):
        proj_dir = Path(project_dir)

        version = load_version(proj_dir.parents[0] / 'VERSION')

        if version is None:
            version = load_version(proj_dir.parents[3] / 'VERSION')

        return VCSDescription(version, 'exact', None, {})
