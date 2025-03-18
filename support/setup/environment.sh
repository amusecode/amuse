# Create a three-character prefix showing whether a package is installed or not.
#
is_installed() {
    installed_name=$(installed_package_name "$1")
    if is_subset "${installed_name}" "${INSTALLED_PACKAGES}" ; then
        printf 'i) '
    else
        printf '   '
    fi
}


# Detect any active Conda or venv environments
#
# This uses the following variables:
#
# CONDA_LIST - output of conda list, as set by the configure script, for improved
#              performance.
#
# This sets the following variables only if an environment is detected:
#
# ENV_TYPE - either "virtualenv" or "conda"
# ENV_NAME - name of the environment
# ENV_LIBRARY_PATH - location of /lib directory
#
detect_environments() {
    if [ "a${VIRTUAL_ENV}" != "a" ] ; then
        ENV_TYPE="virtualenv"
        ENV_NAME="${VIRTUAL_ENV}"
    fi

    if [ "a${CONDA_DEFAULT_ENV}" != "a" ] ; then
        ENV_TYPE="conda"
        ENV_NAME="${CONDA_DEFAULT_ENV}"
    fi
}


# Detect installed packages
#
# Finds installed packages, and checks that pip and wheel are installed, and if in a
# Conda environment, that they were installed using conda and not pip.
#
# This uses the following variables:
#
# ENV_TYPE - see detect_environments()
#
# This sets the following variables:
#
# INSTALLED_PACKAGES - space-separated names of installed packages
# HAVE_PYPI_PIP - if pip is installed from PyPI in a conda environment
# HAVE_PYPI_WHEEL - if wheel is installed from PyPI in a conda environment
# HAVE_PIP - if pip is installed correctly for the current environment
# HAVE_WHEEL - if wheel is installed correctly for the current environment
#
detect_installed_packages() {
    if [ "a${ENV_TYPE}" = "avirtualenv" ] ; then
        INSTALLED_PACKAGES="$(python -m pip list | tail -n +3 | cut -d ' ' -f 1)"
        HAVE_PIP=$(python -m pip list | grep '^pip')
        HAVE_WHEEL=$(python -m pip list | grep '^wheel')
    fi

    if [ "a${ENV_TYPE}" = "aconda" ] ; then
        INSTALLED_PACKAGES="$(echo "${CONDA_LIST}" | tr '^' '\n' | grep -v '^#.*' | cut -d ' ' -f 1)"

        HAVE_PYPI_PIP=$(echo "${CONDA_LIST}" | tr '^' '\n' | grep pypi | grep '^pip')
        HAVE_PYPI_WHEEL=$(echo "${CONDA_LIST}" | tr '^' '\n' | grep pypi | grep '^wheel')

        HAVE_PIP=$(echo "${CONDA_LIST}" | tr '^' '\n' | grep -v pypi | grep '^pip')
        HAVE_WHEEL=$(echo "${CONDA_LIST}" | tr '^' '\n' | grep -v pypi | grep '^wheel')
    fi

    if is_subset "sapporo_light" "${FEATURES}" ; then
        # with a dash because is_installed converts before searching
        INSTALLED_PACKAGES="${INSTALLED_PACKAGES} sapporo-light"
    fi
}


# Determine if we have the required features to build the framework
#
# This uses the following variables:
#
# FEATURES - from configure
#
# This sets the following variables:
#
# ENABLED_PACKAGES - adds amuse-framework if all requirements are met
# ENABLED_PACKAGES_TEXT - adds amuse-framework if all requirements are met
# DISABLED_PACKAGES - adds amuse-framework if features are missing
# DISABLED_PACKAGES_TEXT - adds amuse-framework if features are missing
#
check_build_framework() {
    missing_features=$(filter_out "${FEATURES}" "c c++ fortran python python-dev install mpi")

    if [ "a${missing_features}" = "a" ] ; then
        installed="$(is_installed amuse-framework)"
        ENABLED_PACKAGES="${ENABLED_PACKAGES}amuse-framework "
        ENABLED_PACKAGES_TEXT="${ENABLED_PACKAGES_TEXT}${installed}amuse-framework\n"
    else
        DISABLED_PACKAGES="${DISABLED_PACKAGES}amuse-framework "
        DISABLED_PACKAGES_TEXT="${DISABLED_PACKAGES_TEXT}amuse-framework (missing features:${COLOR_RED}${missing_features}${COLOR_END})\n"
    fi
}


# Determine if we have the required features to build Sapporo Light
#
# This uses the following variables:
#
# FEATURES - from configure
#
# This sets the following variables:
#
# ENABLED_PACKAGES - adds sapporo_light if all requirements are met
# ENABLED_PACKAGES_TEXT - adds sapporo_light if all requirements are met
# DISABLED_PACKAGES - adds sapporo_light if features are missing
# DISABLED_PACKAGES_TEXT - adds sapporo_light if features are missing
#
check_build_sapporo_light() {
    missing_features=$(filter_out "${FEATURES}" "c c++ install cuda")
    if [ "a${missing_features}" = "a" ] ; then
        installed="$(is_installed sapporo_light)"
        ENABLED_PACKAGES="${ENABLED_PACKAGES}sapporo_light "
        ENABLED_PACKAGES_TEXT="${ENABLED_PACKAGES_TEXT}${installed}sapporo_light\n"
    else
        DISABLED_PACKAGES="${DISABLED_PACKAGES}sapporo_light "
        DISABLED_PACKAGES_TEXT="${DISABLED_PACKAGES_TEXT}sapporo_light (missing features:${COLOR_RED}${missing_features}${COLOR_END})\n"
    fi
}


# Check which packages can be built with the available dependencies
#
# This uses the following variables:
#
# FEATURES - from configure
#
# This sets the following variables:
#
# ENABLED_PACKAGES - adds packages for which all requirements are met
# ENABLED_PACKAGES_TEXT - adds packages for which all requirements are met
# DISABLED_PACKAGES - adds packages for which features are missing
# DISABLED_PACKAGES_TEXT - adds packages for which features are missing
# BROKEN_PACKAGES - adds packages that are broken (issue_x dependency)
#
find_packages() {
    for code in src/amuse/community/* ; do
        for dep_file in "${code}"/packages/*.amuse_deps ; do
            # If no file matches, the loop will still run with the pattern as dep_file
            if [ ! -f "$dep_file" ] ; then
                continue
            fi
            package=$(basename "${dep_file}" .amuse_deps)
            deps=$(cat "${dep_file}")
            deps="amuse-framework gmake ${deps}"
            missing_features=$(filter_out "${FEATURES}" "${deps}")
            missing_features=$(filter_out "${ENABLED_PACKAGES}" "${missing_features}")

            if is_subset "sapporo_light" "${deps}" ; then
                NEEDS_SAPPORO_LIGHT="${NEEDS_SAPPORO_LIGHT} ${package}"
            fi

            if [ "a${missing_features}" = "a" ] ; then
                installed="$(is_installed ${package})"
                ENABLED_PACKAGES="${ENABLED_PACKAGES}${package} "
                ENABLED_PACKAGES_TEXT="${ENABLED_PACKAGES_TEXT}${installed}${package}\n"
            elif has_issue "${missing_features}" ; then
                BROKEN_PACKAGES_TEXT="${BROKEN_PACKAGES_TEXT}${package} (reference:${COLOR_RED}${missing_features}${COLOR_END})\n"
            else
                DISABLED_PACKAGES="${DISABLED_PACKAGES}${package} "
                DISABLED_PACKAGES_TEXT="${DISABLED_PACKAGES_TEXT}${package} (missing features:${COLOR_RED}${missing_features}${COLOR_END})\n"
            fi
        done
    done
}


# Analyse the environment and set variables
#
# This checks if we have an environment and if pip and wheel are available, detects
# installed packages, checks whether we can build the framework and sapporo_light and
# then finds which other packages can be installed.
#
# This uses the following variables:
#
# CONDA_LIST - output of conda list, as set by the configure script, for improved
#              performance.
#
# This sets the following variables only if an environment is detected:
#
# ENV_TYPE - either "virtualenv" or "conda"
# ENV_NAME - name of the environment
# ENV_LIBRARY_PATH - location of /lib directory
# HAVE_PYPI_PIP - if pip is installed from PyPI in a conda environment
# HAVE_PYPI_WHEEL - if wheel is installed from PyPI in a conda environment
# HAVE_PIP - if pip is installed correctly for the current environment
# HAVE_WHEEL - if wheel is installed correctly for the current environment
#
# And additionally it sets:
#
# INSTALLED_PACKAGES - space-separated names of installed packages
# ENABLED_PACKAGES - list of packages that can be installed
# ENABLED_PACKAGES_TEXT - printable description of packages that can be installed
# DISABLED_PACKAGES - list of packages that can not be installed
# DISABLED_PACKAGES_TEXT - printable description of packages that can not be installed
#
analyse_environment() {
    detect_environments
    detect_installed_packages
    check_build_framework
    check_build_sapporo_light
    find_packages
}


# Checks for any environment variables that may conflict with the build
check_shell_environment() {
    if [ "a${ENV_TYPE}" = "aconda" ] ; then
        # There's a whole lot of OMPI_ variables and we don't want to check them all one
        # by one. So we do it like this.
        ompi_vars="$(env | grep '^OMPI_' 2>/dev/null)"
        ompi_var_names="$(printf '%s' "${ompi_vars}" | cut -d '=' -f 1 | tr '\n' ' ')"
        if [ "a${ompi_vars}" != "a" ] ; then
            printf '%b\n' "${COLOR_RED}Warning:${COLOR_END} The following shell variables are set, and may"
            printf '%s\n' "interfere with the build:"
            printf '\n%s\n\n' "${ompi_vars}"
            printf '%s\n' "When installing with Conda, we use Conda-installed tools"
            printf '%s\n' "and MPI, and the OMPI_* variables are not needed and will"
            printf '%s\n' "confuse the build system. Please unset them using"
            printf '\n%s\n\n' "    unset ${ompi_var_names}"
        fi
    fi
}

