# Check whether the basic conditions for installing things are met
#
# This checks for an environment with pip and wheel, and quits with an error if they're
# not available.
#
# Args:
#   target: Target the user wants to install
#
check_install() {
    target="$1"

    if [ "a${ENV_TYPE}" = "a" ] ; then
        printf '\n%s\n\n' "Cannot install ${target}, because there is no active environment."
        print_environment_step
        exit 1
    elif [ "a${HAVE_PIP}" = "a" ] || [ "a${HAVE_WHEEL}" = "a" ] ; then
        printf '\n%s\n\n' "Cannot install ${target}, because pip and/or wheel are not available."
        print_pip_wheel_step
        exit 1
    fi
}


# Check that the AMUSE framework can be installed and error out if not
#
check_framework() {
    if ! is_subset amuse-framework "${ENABLED_PACKAGES}" ; then
        printf '%s\n' 'The AMUSE framework cannot be installed because tools or dependencies are missing.'
        printf '%s\n' 'Please run ./setup and follow the instructions to enable it.'
        exit 1
    fi
}


# Check Sapporo Light
#
check_sapporo_light() {
    if ! is_subset "sapporo_light" "${ENABLED_PACKAGES}" ; then
        printf '%s\n' 'Sapporo light cannot be installed because tools or dependencies are missing.'
        printf '%s\n' 'Please run ./setup and follow the instructions to enable it.'
        exit 1
    fi
}


# Check whether a package can be installed and error out if not
#
check_package() {
    package="$1"

    if is_subset "${package}" "${DISABLED_PACKAGES}" ; then
        printf '\n%s\n' "Package ${package} cannot be installed because tools or dependencies are missing."
        printf '%s\n\n' 'Please run ./setup and follow the instructions to enable it.'
        exit 1
    fi

    code_dir=$(code_directory "${package}")

    if [ ! -f "${code_dir}/packages/${package}.amuse_deps" ] ; then
        printf '\n%s\n' "Package ${package} was not found."
        printf '%s\n\n' 'Please run ./setup to show enabled packages.'
        exit 1
    fi
}


# Install the AMUSE framework
#
install_framework() {
    check_framework
    support/shared/uninstall.sh amuse-framework

    announce_activity install amuse-framework

    # if we're in a conda env, install the dependencies using conda first rather than
    # leaving it to pip.
    if [ "a${ENV_TYPE}" = "aconda" ] ; then
        to_install=''
        for name_ver in ${FRAMEWORK_CONDA_DEPS}  ; do
            name=$(echo "${name_ver}" | sed -e "s/'\([a-zA-Z0-9_-]*\).*/\1/")
            if ! is_subset "$name" "${INSTALLED_PACKAGES}" ; then
                to_install="${to_install} ${name_ver}"
            fi
        done
        if [ -n "${to_install}" ] ; then
            conda install -y ${to_install}
        fi
    fi

    ec_file="$(exit_code_file install amuse-framework)"
    log_file="$(log_file install amuse-framework)"

    (
        ${GMAKE} -C lib distclean && \
        ${GMAKE} -C lib install && \
        cd src && pip --no-cache-dir install .

        echo $? >"../${ec_file}"
    ) 2>&1 | tee "${log_file}"

    result=$(cat "${ec_file}")
    if [ "a${result}" = "a0" ] ; then
        INSTALLED_PACKAGES="${INSTALLED_PACKAGES} amuse-framework"
    fi

    handle_result "${result}" install amuse-framework "${log_file}"
}


# Install Sapporo Light
#
install_sapporo_light() {
    check_sapporo_light

    announce_activity install sapporo_light

    ${GMAKE} -C lib/sapporo_light distclean

    ec_file="$(exit_code_file install sapporo_light)"
    log_file="$(log_file install sapporo_light)"

    (${GMAKE} -C lib install-sapporo_light ; echo $? >"${ec_file}") 2>&1 | tee "${log_file}"

    result=$(cat "${ec_file}")
    if [ "a${result}" = "a0" ] ; then
        INSTALLED_PACKAGES="${INSTALLED_PACKAGES} sapporo-light"
    fi

    handle_result "${result}" install sapporo_light "${log_file}"
}


# Install the AMUSE framework in develop mode
#
develop_framework() {
    check_framework
    support/shared/uninstall.sh amuse-framework

    announce_activity develop amuse-framework

    ec_file="$(exit_code_file install amuse-framework)"
    log_file="$(log_file install amuse-framework)"

    (
        ${GMAKE} -C lib distclean && \
        ${GMAKE} -C lib install && \
        cd src && pip install -e .

        echo $? >"../${ec_file}"
    ) 2>&1 | tee "${log_file}"

    handle_result $(cat "$ec_file") install amuse-framework "${log_file}"
}


# Install the AMUSE framework into a packager's build environment
#
# This is set up for building Conda packages, using their recommended pip options.
#
package_framework() {
    check_framework
    ${GMAKE} -C lib distclean
    ${GMAKE} -C lib install
    (cd src && python3 -m pip install -vv --no-cache-dir --no-deps --no-build-isolation --prefix "${PREFIX}" .)
}


# Install a package
#
# This calls the package's build system to install it.
#
# If amuse-framework is not installed, then we install the framework first. Before
# installing, make distclean is run to improve reliability. If you want to take your
# chances with an incremental build and the native build system, cd into the package
# directory and call make directly.
#
# Args:
#    cmd: The type of installation, install or develop
#    package: The name of the package to install
#    brief: If set to "brief", print only a brief result, otherwise, print a full error.
#
install_package() {
    cmd="$1"
    package="$2"
    brief="$3"

    check_package "${package}"

    if ! is_subset "amuse-framework" "${INSTALLED_PACKAGES}" ; then
        save_package="${package}"
        install_framework
        package="${save_package}"
    fi

    if is_subset "${package}" "${NEEDS_SAPPORO_LIGHT}" ; then
        if ! is_subset "sapporo-light" "${INSTALLED_PACKAGES}" ; then
            save_package="${package}"
            install_sapporo_light
            package="${save_package}"
        fi
    fi

    if [ "a${package%-*-*}" != "a${package}" ] ; then
        # We are installing an amuse-code-package extension package, so we need the
        # base package as well, if it exists.
        base_package="${package%-*}"
        # If the code is e.g. CUDA-only, then there may not be a base package.
        if is_subset "${base_package}" "${EXTANT_PACKAGES}" ; then
            if ! is_subset "${base_package}" "${INSTALLED_PACKAGES}" ; then
                save_package="${package}"
                install_package "${cmd}" "${base_package}" "${brief}"
                package="${save_package}"
            fi
        fi
    fi

    save_cmd="${cmd}"
    forward_to_package "distclean" "${package}" "${brief}"
    cmd="${save_cmd}"

    forward_to_package "${cmd}" "${package}" "${brief}"

    result="$?"
    if [ "a${result}" = "a0" ] ; then
        INSTALLED_PACKAGES="${INSTALLED_PACKAGES} ${package}"
    fi

    return "${result}"
}


# Install the framework and all enabled packages
#
install_all() {
    if ! is_subset "amuse-framework" "${INSTALLED_PACKAGES}" ; then
        install_framework || exit 1
    fi

    FAILED_BUILDS=''
    if ! is_subset "sapporo-light" "${INSTALLED_PACKAGES}" ; then
        if is_subset "sapporo_light" "${ENABLED_PACKAGES}" ; then
            if ! install_sapporo_light ; then
                FAILED_BUILDS="${FAILED_BUILDS}\nsapporo_light"
            fi
        fi
    fi

    for package in ${ENABLED_PACKAGES} ; do
        installed_name=$(installed_package_name "${package}")
        if ! is_subset "${installed_name}" "${INSTALLED_PACKAGES}" ; then
            install_package install "${package}" brief
            if [ $? != '0' ] ; then
                FAILED_BUILDS="${FAILED_BUILDS}\n${package}"
            fi
        fi
    done

    printf '\n%b\n' "${COLOR_CYAN}*** Build results ***${COLOR_END}"

    if [ "${FAILED_BUILDS#*amuse-mesa-r15140}" != "${FAILED_BUILDS}" ] ; then
        printf '\n%b\n' "${COLOR_RED}amuse-mesa-r15140 failed to install${COLOR_END}"
        printf '\nMESA r15140 failed to install because it can only be installed in\n'
        printf 'develop mode. This is expected, if unfortunate. Please use\n\n'
        printf '    ./setup develop amuse-mesa-r15140\n\n'
        printf 'to install it.\n\n'

        FAILED_BUILDS="${FAILED_BUILDS#*amuse-mesa-r15140}"
    fi

    if [ "a${FAILED_BUILDS}" != "a" ] ; then
        printf '\nThe following packages failed to build:\n'
        printf '%b\n\n' "${COLOR_RED}${FAILED_BUILDS}${COLOR_END}"
        print_getting_help
        printf '\n%s\n' 'Output was saved to support/logs/.'
    else
        printf '\n%b\n\n' "${COLOR_GREEN}All enabled packages were installed successfully${COLOR_END}"
    fi
}


# Check whether the basic conditions for uninstalling things are met
#
# This checks for an environment, and quits with an error if they're not available.
#
# Args:
#   target: Target the user wants to install
#
check_uninstall() {
    target="$1"

    if [ "a${ENV_TYPE}" = "a" ] ; then
        printf '\n%s\n\n' "Cannot uninstall ${target}, because there is no active environment."
        print_environment_step
        exit 1
    fi
}


# Uninstall the AMUSE framework
#
uninstall_framework() {
    for pkg in ${INSTALLED_PACKAGES} ; do
        if [ "a${pkg#amuse-}" != "a${pkg}" ] ; then
            if [ "a${pkg}" != "aamuse-framework" ] ; then
                # We may have already uninstalled this as a dependency of another
                # package that got uninstalled, in which case we skip.
                # INSTALLED_PACKAGES is kept up-to-date as we go.
                if is_subset "${pkg}" "${INSTALLED_PACKAGES}" ; then
                    save_package="${package}"
                    uninstall_package "${pkg}" brief
                    package="${save_package}"
                fi
            fi
        fi
    done

    announce_activity uninstall amuse-framework

    ec_file="$(exit_code_file uninstall amuse-framework)"
    log_file="$(log_file uninstall amuse-framework)"

    (
        printf '%s\n\n' "Removing stray libraries, if any..." && \
        ${GMAKE} -C lib uninstall && \
        printf '\n' && \
        support/shared/uninstall.sh amuse-framework

        echo $? >"${ec_file}"
    ) 2>&1 | tee "${log_file}"

    handle_result $(cat "${ec_file}") uninstall amuse-framework "${log_file}"
}


# Uninstall Sapporo Light
#
uninstall_sapporo_light() {
    for pkg in ${INSTALLED_PACKAGES} ; do
        if [ "a${pkg#amuse-}" != "a${pkg}" ] ; then
            if is_subset "${pkg}" "${NEEDS_SAPPORO_LIGHT}" ; then
                save_package="${package}"
                uninstall_package "${pkg}" brief
                package="${save_package}"
            fi
        fi
    done

    announce_activity uninstall sapporo_light

    ec_file="$(exit_code_file uninstall sapporo_light)"
    log_file="$(log_file uninstall sapporo_light)"

    (${GMAKE} -C lib uninstall-sapporo_light ; echo $? >"${ec_file}") 2>&1 | tee "${log_file}"

    handle_result $(cat "$ec_file") uninstall sapporo_light "${log_file}"
}


# Uninstall a package
#
# This calls the shared uninstall script to uninstall the package. That script will use
# pip or conda to uninstall, as appropriate.
#
# Args:
#    package: The name of the package to uninstall
#    brief: If set to "brief", print only a brief result, otherwise, print a full error.
#
uninstall_package() {
    package="$1"
    brief="$2"

    dependents=$(extra_packages "${package}" "${INSTALLED_PACKAGES}")

    for pkg in $dependents ; do
        save_package="${package}"
        uninstall_package "${pkg}" brief
        package="${save_package}"
    done

    announce_activity uninstall "${package}"

    ec_file="$(exit_code_file uninstall ${package})"
    log_file="$(log_file uninstall ${package})"

    pkg_name=$(installed_package_name "${package}")

    (support/shared/uninstall.sh "${pkg_name}" ; echo $? >"${ec_file}") 2>&1 | tee "${log_file}"

    result=$(cat "$ec_file")

    if [ "a${result}" = "a0" ] ; then
        INSTALLED_PACKAGES=$(filter_out "${package}" "${INSTALLED_PACKAGES}")
    fi

    handle_result ${result} uninstall "${package}" "${log_file}" "${brief}"
}

