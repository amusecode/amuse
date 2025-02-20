# Print the base name for the log and exit code files
#
log_file_base() {
    printf '%s' "support/logs/${1}-${2}"
}


# Print the relative path to the exit code file for a command and target
#
# We show output on the terminal when installing and testing, but we also send it to a
# file so that when something goes wrong, there's something for people to send us when
# they ask for help. We also need to send the exit code to a file, otherwise we lose it
# because of the pipe to tee and we need it to see if the operation was successful. So
# this points to the file we'll use for that.
#
# Args:
#   command: The command that's being executed
#   package: The package the command is being executed for
#
exit_code_file() {
    name="$(log_file_base $1 $2)"
    printf '%s' "${name}.exit_code"
}


# Print the relative path to the log file for a command and target.
#
# See exit_code_file() above, this is the location for the log output.
#
log_file() {
    name="$(log_file_base $1 $2)"
    printf '%s' "${name}.log"
}


# Announce that we're going to do something
#
# This is a counterpart to handle_result() below
#
# Args:
#    cmd: Command that we're executing
#    target: Affected package
#
announce_activity() {
    cmd="$1"
    package="$2"

    printf '\n'
    if [ "a${cmd}" = "aclean" ] || [ "a${cmd}" = "adistclean" ] ; then
        printf '%b\n' "${COLOR_CYAN}Cleaning ${package}...${COLOR_END}"
    elif [ "a${cmd}" = "ainstall" ] ; then
        printf '%b\n' "${COLOR_CYAN}Building and installing ${package}...${COLOR_END}"
    elif [ "a${cmd}" = "adevelop" ] ; then
        printf '%b\n' "${COLOR_CYAN}Building and develop-installing ${package}...${COLOR_END}"
    elif [ "a${cmd}" = "atest" ] ; then
        printf '%b\n' "${COLOR_CYAN}Testing ${package}...${COLOR_END}"
    elif [ "a${cmd}" = "auninstall" ] ; then
        printf '%b\n' "${COLOR_CYAN}Uninstalling ${package}...${COLOR_END}"
    fi
    printf '\n'
}


# Handle the result of an installation or test step
#
# These are never supposed to fail if you install from an official distribution, but of
# course nothing is perfect. So we need to handle any errors gracefully.
#
# Args:
#   exit_code: Exit code of the process
#   command: What we did, e.g. install or test
#   package: The package we did it to
#   log_file: Location of the log file
#   brief: If set to "brief", print only a short error, otherwise the full message.
#
handle_result() {
    exit_code=$1
    cmd=$2
    package=$3
    log_file="$4"
    brief="$5"

    if [ "a${exit_code}" = "a0" ] ; then
        if [ "a${cmd}" = "ainstall" ] || [ "a${cmd}" = "adevelop" ] ; then
            printf '\n%b\n' "${COLOR_GREEN}Package ${package} was installed successfully.${COLOR_END}"
        elif [ "a${cmd}" = "atest" ] ; then
            printf '\n%b\n' "${COLOR_GREEN}Package ${package} passed its tests.${COLOR_END}"
        elif [ "a${cmd}" = "auninstall" ] ; then
            printf '\n%b\n' "${COLOR_GREEN}Package ${package} uninstall successful.${COLOR_END}"
        fi
    else
        if [ "a${cmd}" = "ainstall" ] || [ "a${cmd}" = "adevelop" ] ; then
            if [ "a${brief}" = "abrief" ] ; then
                printf '\n%b\n' "${COLOR_RED}Package ${package} failed to install correctly.${COLOR_END}"
            else
                print_install_failure "${package}" "${log_file}"
            fi
        elif [ "a${cmd}" = "atest" ] ; then
            if [ "a${brief}" = "abrief" ] ; then
                printf '\n%b\n' "${COLOR_RED}Package ${package} failed its tests.${COLOR_END}"
            else
                print_test_failure "${package}" "${log_file}"
            fi
        elif [ "a${cmd}" = "auninstall" ] ; then
            if [ "a${brief}" = "abrief" ] ; then
                printf '\n%b\n' "${COLOR_RED}Package ${package} failed to uninstall.${COLOR_END}"
            else
                print_uninstall_failure "${package}" "${log_file}"
            fi
        fi
    fi

    return "${exit_code}"
}

