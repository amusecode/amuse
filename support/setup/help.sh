. support/setup/format.sh
. support/setup/util.sh
. support/setup/dependencies.sh


print_invalid_command() {
    cmd="$1"

    printf '\n%b\n' "${COLOR_RED}\"${cmd}\" was not recognised as a valid command.${COLOR_END}"
    printf '\n%s\n' 'Try'
    printf '\n%s\n' '    ./setup help'
    printf '\n%s\n' 'to see all available commands.'
}


print_help() {
    printf '\n%b\n' "${BOLD}${COLOR_CYAN}*** AMUSE setup help ***${COLOR_END}${END_BOLD}"

    printf '%b\n' "
This setup script will help you install AMUSE and its various components from source.
It must always be called from the directory it is in, as ${BOLD}./setup${END_BOLD}.

To get started quickly, just type ${BOLD}./setup${END_BOLD}, press Enter, and follow the instructions.

Here follows an overview of all supported commands.

${BOLD}Guided setup${END_BOLD}

    ./setup

This will check whether you have an environment active, whether it is set up correctly,
and whether required dependencies are installed. It will then give suggestions for how
to accomplish these things and install AMUSE.

${BOLD}Installing AMUSE${END_BOLD}

    ./setup install all

Builds and installs the framework and then all enabled packages into the active
environment.

    ./setup install amuse-framework

Builds and installs only the framework into the active environment.

    ./setup install sapporo_light

Builds and installs the Sapporo light GPU nbody library into the active environment.

    ./setup install ${ITALIC}package1 package2...${END_ITALIC}

Builds and installs specific package(s) into the active environment. Any number of
packages can be specified.

${BOLD}Developing AMUSE${END_BOLD}

    ./setup develop amuse-framework

Builds and installs the AMUSE framework in develop mode, as an editable install, into
the active environment.

    ./setup develop ${ITALIC}package1 package2...${END_ITALIC}

Builds and installs a specific package in develop mode, as an editable install, into
the active environment.

    ./setup test all

Runs tests for the framework and all enabled packages. This requires those packages to
be installed in the active environment.

    ./setup test amuse-framework

Runs tests for the framework. This requires the framework to be installed in the active
environment.

    ./setup test amuse-ext

Runs tests for the amuse.ext scripts, and regression ('ticket') tests. This requires the
framework to be installed as well as (currently) bhtree, evtwin, fi, fractalcluster,
gadget2, galactics, halogen, hermite, hop, kepler, mesa-r14150, ph4, phigrape, seba, and
sse.

    ./setup test ${ITALIC}package1 package2...${END_ITALIC}

Runs tests for the specified package(s). This requires each package to be installed in
the active environment.

    ./setup clean

Cleans up all the compiled code for both the framework and the community codes, so that
you can (and have to!) rebuild everything from scratch.

    ./setup distclean

Like clean, but also cleans up any built-in dependencies that get built for some
community codes. Note that this does not remove downloaded files.
"
}


print_environment_step() {
    printf '%b\n' "${BOLD}${COLOR_YELLOW}* (1/4) Create and activate an environment *${COLOR_END}${END_BOLD}

To install any packages, you need to create and activate a Conda environment
or a Python virtualenv. If you have an environment into which you'd like to
install AMUSE, you should activate it now. To create a new Conda environment,
use

    conda create -n Amuse-env

Then you activate it using

    conda activate Amuse-env

You can name the environment anything you like instead of my-amuse-env.

To create a Python virtualenv, use

    python3 -m venv venv

Then you can activate it using

    source venv/bin/activate

You can put the environment anywhere you like, but do make sure that the
directory it will be in exists, or you'll get an error message.

Once you have an environment active, type

    ./setup

again to continue."
}


print_pip_wheel_step() {
    printf '%b\n' "${BOLD}${COLOR_YELLOW}* (2/4) Install pip and wheel *${COLOR_END}${END_BOLD}

Installation is disabled due to a lack of the right pip and/or wheel in the
environment. You can enable AMUSE installation by (correctly) installing pip and wheel.

To do that, use
"

    if [ "a${ENV_TYPE}" = "avirtualenv" ] ; then
        printf '    %s\n' 'python3 -m pip install pip wheel'
    elif [ "a${ENV_TYPE}" = "aconda" ] ; then
        if [ "a${HAVE_PYPI_WHEEL}" != "a" ] ; then
            printf '    %s\n' 'python3 -m pip uninstall wheel'
        fi
        if [ "a${HAVE_PYPI_PIP}" != "a" ] ; then
            printf '    %s\n' 'python3 -m pip uninstall pip'
        fi
        printf '    %s\n' 'conda install -c conda-forge pip wheel'
    fi

    printf '\n%b\n' "and then run

    ./setup

again to continue."
}


print_enable_packages_step() {
    printf '%b\n' "${BOLD}${COLOR_YELLOW}* (3/4) Enable more packages *${COLOR_END}${END_BOLD}

Some packages are disabled due to missing features. You can enable more packages by
installing additional software. Some software does require specific hardware, for
example CUDA requires an nVidia GPU to work. If you have one, see
https://amuse.readthedocs.io/en/latest/install/cuda.html for help installing CUDA to
enable those packages."

    if [ "a${ENV_TYPE}" = "aconda" ] ; then
        printf '%b\n' "
To install the dependencies using conda, use

${CONDA_CMDS}
"
    elif [ "a${ENV_TYPE}" = "avirtualenv" ] ; then
        if [ "a${MACPORTS}" != "a" ] ; then
            printf '%b\n' "
We seem to be running on a Mac with MacPorts installed. To install the dependencies
system-wide using MacPorts, use

${MACPORTS_CMDS}
"
        elif [ "a${HOMEBREW}" != "a" ] ; then
            printf '%b\n' "
We seem to be running on a Mac with Homebrew installed. To install the dependencies
system-wide using Homebrew, use

${HOMEBREW_CMDS}
"
        elif [ "a${APT}" != "a" ] ; then
            printf '%b\n' "
We seem to be running on Ubuntu or a similar Linux distribution that uses APT.
To install the dependencies system-wide using apt, use

${APT_CMDS}
"
        elif [ "a${DNF}" != "a" ] ; then
            printf '%b\n' "
We seem to be running on RedHat or a similar Linux distribution that uses DNF.
To install the dependencies system-wide using DNF, use

${DNF_CMDS}
"
        fi
    fi
}


print_install_amuse_step() {
    printf '%b\n' "${BOLD}${COLOR_YELLOW}* (4/4) Install AMUSE *${COLOR_END}${END_BOLD}

To install the AMUSE framework and all of the enabled packages into the active
${ENV_TYPE} environment ${BOLD}${ENV_NAME}${END_BOLD}, type

    ./setup install all

To install only the AMUSE framework, you can use

    ./setup install amuse-framework

To install a specific package, you can use

    ./setup install amuse-bhtree

or whichever package you want to install, as long as it's enabled.

There are some more commands available, use

    ./setup help

to show a complete overview.
"
}


print_getting_help() {
    printf '%s\n\n' 'This is not supposed to happen!'
    printf '%s\n' 'Please contact the developers on Slack or make an issue on GitHub,'
    printf '%s\n\n' 'so that we can fix this together.'
    printf '%s\n' 'https://amusecode.slack.com'
    printf '%s\n' 'https://github.com/amusecode/amuse/issues/new/choose'
}


print_install_failure() {
    package="$1"
    log_file="$2"

    printf '\n%b\n\n' "${COLOR_RED}${package} failed to build and install correctly.${COLOR_END}"
    print_getting_help
    printf '\n%s\n' "The output of the installation process was logged to $2"
}


print_test_failure() {
    package="$1"
    log_file="$2"

    printf '\n%b\n\n' "${COLOR_RED}${package} failed its tests.${COLOR_END}"
    print_getting_help
}


print_uninstall_failure() {
    package="$1"
    log_file="$2"

    printf '\n%b\n\n' "${COLOR_RED}${package} failed to uninstall correctly.${COLOR_END}"
    print_getting_help
    printf '\n%s\n' "The output of the uninstallation process was logged to $2"
}

