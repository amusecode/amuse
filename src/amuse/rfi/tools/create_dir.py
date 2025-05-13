from pathlib import Path
from shutil import copyfile, copymode, copytree
from typing import Dict

from amuse import __version__ as amuse_version


_template_dir = Path(__file__).parent / "dir_templates"


"""Map from template file to destination inside code dir, for C++."""
_templates_cxx = {
    "src_code.cc": "src/{code}.cc",
    "src_Makefile_cxx": "src/Makefile",
    "interface.cc": "interface.cc",
    "interface.py": "interface.py",
    "Makefile_cxx": "Makefile"
    }


"""Map from template file to destination inside code dir, for Fortran."""
_templates_fortran = {
    "src_code.f90": "src/{code}.f90",
    "src_Makefile_fortran": "src/Makefile",
    "interface.f90": "interface.f90",
    "interface.py": "interface.py",
    "Makefile_fortran": "Makefile"
    }


def _instantiate_template(
        tmpl_file: str, target: Path, variables: Dict[str, str]
        ) -> None:
    """Instantiates a template file at the given location.

    Args:
        tmpl_file: Name of a file in the dir_templates directory
        target: Path to the output file
        variables: Variables to substitute into the template
    """
    print(f"Instantiating {tmpl_file} to {target}")
    template = (_template_dir / tmpl_file).read_text("utf-8")
    target.write_text(template.format(**variables), "utf-8")


def _make_directory(code_dir: Path, variables: Dict[str, str]) -> None:
    """Make the basic directory tree."""
    try:
        code_dir.mkdir(parents=True)
    except FileExistsError:
        raise RuntimeError(
                f"The directory {code_dir} already exists, aborting to avoid"
                " overwriting.")

    _instantiate_template("__init__.py", code_dir / "__init__.py", variables)


_tmpls_support_cxx = {
    "configure_cxx.ac": "configure.ac",
    "config_cxx.mk.in": "config.mk.in",
}


_tmpls_support_fortran = {
    "configure_fortran.ac": "configure.ac",
    "config_fortran.mk.in": "config.mk.in",
}


_files_support = {
    "support_shared_config.guess": "config.guess",
    "support_shared_config.sub": "config.sub",
    "support_shared_install-sh": "install-sh",
    "support_shared_uninstall.sh": "uninstall.sh",
    }


def _make_support(code_dir: Path, language: str, variables: Dict[str, str]) -> None:
    """Make the support/ subdirectory with contents."""
    support_dir = code_dir / "support"
    support_dir.mkdir()

    tmpls_paths = _tmpls_support_cxx if language == "c" else _tmpls_support_fortran
    for tmpl, path in tmpls_paths.items():
        _instantiate_template(tmpl, support_dir / path.format(**variables), variables)

    # TODO: this needs to be updated when we move the codes out of community/
    shared_base = code_dir / ".." / ".." / ".." / "support" / "shared"
    shared_dir = support_dir / "shared"
    if shared_base.exists():
        # We're probably inside the AMUSE tree, so we make a symlink like for the other
        # embedded codes. This avoids having to update many copies of e.g. a broken m4
        # macro.
        shared_dir.symlink_to(shared_base)
    else:
        shared_dir.mkdir()
        for src, dst in _files_support.items():
            copyfile(_template_dir / src, shared_dir / dst)
            copymode(_template_dir / src, shared_dir / dst)

        copytree(_template_dir / "support_shared_m4", shared_dir / "m4")


def _make_wrapper(code_dir: Path, language: str, variables: Dict[str, str]) -> None:
    """Make the wrapper and the dummy code."""
    (code_dir / "src").mkdir()

    tmpls_paths = _templates_cxx if language == "c" else _templates_fortran

    for tmpl, path in tmpls_paths.items():
        _instantiate_template(tmpl, code_dir / path.format(**variables), variables)


def _make_tests(code_dir: Path, variables: Dict[str, str]) -> None:
    """Make the test."""
    (code_dir / "tests").mkdir()
    file = "test_{code}.py".format(**variables)
    _instantiate_template("tests_test_code.py", code_dir / "tests" / file, variables)


def _make_packages(code_dir: Path, code: str, variables: Dict[str, str]) -> None:
    """Make the packages/ subdirectory with contents."""
    (code_dir / "packages").mkdir()

    def make_package(pkg_type: str, suffix: str) -> None:
        pkg_dir = code_dir / "packages" / f"amuse-{code}{suffix}"
        pkg_dir.mkdir()
        (pkg_dir / code).symlink_to("../..", True)

        _instantiate_template(
                "amuse-code.amuse_deps",
                pkg_dir.parent / f"amuse-{code}{suffix}.amuse_deps", variables)

        _instantiate_template(
                f"pyproject-{pkg_type}.toml", pkg_dir / "pyproject.toml", variables)

    make_package("base", "")
    make_package("extra", "-PACKAGE")


def _variables(language: str, code: str, user_class: str) -> None:
    """Create variables to instantiate the templates with."""
    variables = {
            "framework_version": amuse_version,
            "code": code,
            "package": code.replace("_", "-"),
            "interface_class": f"{user_class}Interface",
            "user_class": user_class}

    if '+' in amuse_version:
        # We have a local version label, meaning the framework has been modified
        # locally. Versions like that can only be compared by equality, which is
        # probably good because it may be incompatible with any other versions. So we
        # required the exact version.
        variables["framework_version"] = f"=={amuse_version}"
    else:
        # We have a normal version, probably a normal release. In this case we assume
        # that future versions of the framework will continue to work.
        variables["framework_version"] = f">={amuse_version}"

    if language == "c":
        variables["includes"] = f"include_headers = [\"{code}_worker.h\"]"
    elif language == "f90":
        interface_module = variables["interface_class"]
        variables["interface_module"] = interface_module
        variables["includes"] = f"use_modules = [\"{interface_module}\"]"

    return variables


def create_code_dir(language: str, user_class: str, work_dir: str) -> None:
    """Create a directory for a new community code.

    This creates a new directory for a new community code, with instantiated templates
    to get the developer started with connecting the code to AMUSE.

    Note that the value of "code" will be used to name the directory (after having been
    lowercased) and also the class that the AMUSE user will use to interact with the
    code (with an initial capital).

    Args:
        language: Programming language, either "c" (also for C++) or "f90"
        code: Name of the new code
        work_dir: Directory inside of which the new directory will be made
    """
    code = user_class.lower()
    code_dir = Path(work_dir) / f'amuse_{code}'

    variables = _variables(language, code, user_class)

    _make_directory(code_dir, variables)
    _make_support(code_dir, language, variables)
    _make_wrapper(code_dir, language, variables)
    _make_tests(code_dir, variables)
    _make_packages(code_dir, code, variables)
