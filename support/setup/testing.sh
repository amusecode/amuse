# Check that pytest is available and error out of not
#
ensure_pytest() {
    if [ "a${PYTEST}" = "a" ] ; then
        printf '\n%b\n' "To run the tests pytest is required, but it is not installed."
        if [ "a${ENV_TYPE}" = "a" ] ; then
            printf '%s\n\n' "Please activate a conda environment or virtual environment first."
        elif [ "a${ENV_TYPE}" = "aconda" ] ; then
            printf '\n%s\n' "Please use"
            printf '\n    %b\n' "conda install pytest"
            printf '\n%s\n\n' "to install pytest, then try again."
        else
            printf '\n%s\n' "Please use"
            printf '\n    %b\n' "pip install pytest"
            printf '\n%s\n\n' "to install pytest, then try again."
        fi
        exit 1
    fi
}


# Check whether a package is installed and can therefore be tested
#
# If not, print an error message and quit.
#
check_package_installed_for_test() {
    package="$1"

    if ! is_subset "${package}" "${INSTALLED_PACKAGES}" ; then
        printf '\n%s\n' "Package ${package} is not installed, so we cannot test it."
        printf '\n%s\n' "Please install the package first using ./setup install ${package}, then try again."
        exit 1
    fi
}


# Run tests for the AMUSE framework
#
test_framework() {
    announce_activity test amuse-framework

    ec_file="$(exit_code_file test amuse-framework)"
    log_file="$(log_file test amuse-framework)"

    (
        make -C src/tests all && \
        # Tests for amuse.distributed won't be fixed as it is to be removed, disabled.
        cd src/tests && pytest --pyargs core_tests compile_tests ${PYTEST_OPTS} -k 'not TestCDistributedImplementationInterface and not TestAsyncDistributed'

        echo $? >"../${ec_file}"
    ) 2>&1 | tee "${log_file}"

    handle_result $(cat "$ec_file") test amuse-framework "${log_file}"
}


# Run tests for the AMUSE ext scripts
#
# Note that these require some of the codes to be installed, as they're mostly scenario
# tests that exercise coupled simulations.
#
test_amuse_ext() {
    announce_activity test amuse-ext

    bad_ext_tests=''
    sep=''
    for bad_test in ${BAD_EXT_TESTS} ; do
        # PyTest doesn't match Class::test_name directly
        bad_test_pt=$(echo "${bad_test}" | sed -e 's/\(.*\)::\(.*\)/(\1 and \2)/')
        bad_ext_tests="${bad_ext_tests}${sep} not ${bad_test_pt}"
        sep=' and'
    done

    ec_file="$(exit_code_file test amuse-framework)"
    log_file="$(log_file test amuse-framework)"

    (
        cd src/tests && pytest --pyargs ext_tests ticket_tests ${PYTEST_OPTS}  -k "${bad_ext_tests}"

        echo $? >"${ec_file}"
    ) 2>&1 | tee "${log_file}"

    handle_result $(cat "$ec_file") test amuse-ext "${log_file}"

    printf "\n%s\n" "The following tests were disabled because they currently fail:"
    printf "\n%b\n\n" "${COLOR_RED}${BAD_EXT_TESTS}${COLOR_END}"
    printf "%s\n\n" "This issue is tracked at https://github.com/amusecode/amuse/issues/1103"
}


# Run tests for the framework and all installed packages
#
test_all() {
    FAILED_TESTS=''
    if is_subset amuse_framework "${INSTALLED_PACKAGES}" ; then
        if ! test_framework ; then
            FAILED_TESTS="${FAILED_TESTS}\namuse-framework"
        fi

        if ! test_amuse_ext ; then
            FAILED_TESTS="${FAILED_TESTS}\namuse-ext"
        fi
    fi

    # sapporo_light does not have tests

    for package in ${INSTALLED_PACKAGES} ; do
        if ! is_subset "${package}" "amuse-framework sapporo_light" ; then
            code_dir=$(code_directory "${package}")
            if [ -f "${code_dir}/packages/${package}.amuse_deps" ] ; then
                forward_to_package test "${package}" brief
                if [ $? != '0' ] ; then
                    FAILED_TESTS="${FAILED_TESTS}\n${package}"
                fi
            fi
        fi
    done

    printf '\n%b\n' "${COLOR_CYAN}*** Test results ***${COLOR_END}"
    if [ "a${FAILED_TESTS}" != "a" ] ; then
        printf '\nThe following packages failed their tests:\n'
        printf '%b\n\n' "${COLOR_RED}${FAILED_TESTS}${COLOR_END}"
        print_getting_help
    else
        printf '\n%b\n\n' "${COLOR_GREEN}All installed packages completed their tests successfully${COLOR_END}"
    fi
}

