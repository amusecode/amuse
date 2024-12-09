define ENVIRONMENT_HELP :=

$(BOLD)* (1/4) Create and activate an environment *$(END_BOLD)

To install any packages, you need to create and activate a Conda environment
or a Python virtualenv. If you have an environment into which you'd like to
install AMUSE, you should activate it now. To create a new Conda environment,
use

    conda create -n my-amuse-env

Then you activate it using

    conda activate my-amuse-env

You can name the environment anything you like instead of my-amuse-env.

To create a Python virtualenv, use

    python -m venv ~/Envs/my-amuse-env

Then you can activate it using

    source ~/Envs/my-amuse-env/bin/activate

You can put the environment anywhere you like, but do make sure that the
directory it will be in exists, or you'll get an error message.

Once you have an environment active, type

    make configure

again to continue.

endef    # ENVIRONMENT_HELP

ENVIRONMENT_HELP := $(subst $(newline),\n,$(subst ','\'',$(ENVIRONMENT_HELP)))


define NO_PIP_WHEEL_MESSAGE :=

$(BOLD)* (2/4) Install pip and wheel *$(END_BOLD)

Installation is disabled due to a lack of the right pip and/or wheel in the
environment. You can enable AMUSE installation by correctly installing pip and wheel.

To do that, use


endef

ifeq ($(ENV_TYPE), virtualenv)
define NO_PIP_WHEEL_MESSAGE +=
    python -m pip install pip wheel

endef
endif    # ENV_TYPE is virtualenv

ifeq ($(ENV_TYPE), conda)

ifneq (,$(HAVE_PYPI_WHEEL))
define NO_PIP_WHEEL_MESSAGE +=
    python -m pip uninstall wheel

endef
endif    # HAVE_PYPI_WHEEL

ifneq (,$(HAVE_PYPI_PIP))
define NO_PIP_WHEEL_MESSAGE +=
    python -m pip uninstall pip

endef
endif    # HAVE_PYPI_PIP

define NO_PIP_WHEEL_MESSAGE +=
    conda install -c conda-forge pip wheel

endef
endif    # ENV_TYPE is conda


define NO_PIP_WHEEL_MESSAGE +=

and then run

    make configure

again to continue.

endef

NO_PIP_WHEEL_MESSAGE := $(subst $(newline),\n,$(subst ','\'',$(NO_PIP_WHEEL_MESSAGE)))


define DISABLED_PACKAGES_MESSAGE :=

$(BOLD)* (3/4) Enable more packages *$(END_BOLD)

Some packages are disabled due to missing features. You can enable more packages by
installing additional software. Some software does require specific hardware, for
example CUDA requires an nVidia GPU to work.

endef


ifeq ($(ENV_TYPE), conda)

define DISABLED_PACKAGES_MESSAGE +=

$(BOLD)Installing tools and libraries using Conda$(END_BOLD)

If you are going to install in a Conda environment, then you must use Conda to install
the required packages. You can do that using

$(CONDA_CMDS)

endef

else ifeq ($(ENV_TYPE), virtualenv)

ifneq (, $(MACPORTS))
define DISABLED_PACKAGES_MESSAGE +=

$(BOLD)Installing tools and libraries using MacPorts$(END_BOLD)

We seem to be running on a Mac with MacPorts installed. Here's how to install the
required packages for all codes system-wide using MacPorts:

$(MACPORTS_CMDS)

endef
endif


ifneq (,$(HOMEBREW))
define DISABLED_PACKAGES_MESSAGE +=

$(BOLD)Installing tools and libraries using Homebrew$(END_BOLD)

We seem to be running on a Mac with Homebrew installed. Here's how to install the
required packages for all codes system-wide using Homebrew:

$(HOMEBREW_CMDS)

endef
endif


ifneq (,$(APT))
define DISABLED_PACKAGES_MESSAGE +=

$(BOLD)Installing tools and libraries using APT$(END_BOLD)

We seem to be running on Ubuntu or a similar Linux distribution that uses APT.
Here's how to install the required packages for all codes system-wide using apt:

$(APT_CMDS)

endef
endif


ifneq (,$(DNF))
define DISABLED_PACKAGES_MESSAGE +=

$(BOLD)Installing tools and libraries using DNF$(END_BOLD)

We seem to be running on RedHat or a similar Linux distribution that uses DNF.
Here's how to install the required packages for all codes system-wide using dnf:

$(DNF_CMDS)

endef
endif

endif       # ENV_TYPE is virtualenv

ifneq (cuda, $(filter cuda, $(FEATURES)))

define DISABLED_PACKAGES_MESSAGE +=

$(BOLD)Installing CUDA$(END_BOLD)

If you have an nVidia GPU and a code that can use CUDA, then compiling with CUDA support
will greatly increase performance. Installing CUDA is a bit more complicated than the
other dependencies however. CUDA consists of two parts, the driver and the toolkit. The
driver must be installed on the system, as part of the nVidia GPU drivers. The toolkit
can be installed separately via an installer, the system package manager, or conda.

Here are a few useful links:

Installing nVidia drivers on Ubuntu: https://ubuntu.com/server/docs/nvidia-drivers-installation
NVIDIA CUDA installation guide for Linux: https://docs.nvidia.com/cuda/cuda-installation-guide-linux/index.html
CUDA on Windows Subsystem for Linux: https://docs.nvidia.com/cuda/wsl-user-guide/index.html#getting-started-with-cuda-on-wsl-2

endef

ifeq ($(ENV_TYPE), conda)

define DISABLED_PACKAGES_MESSAGE +=

To install the CUDA toolkit using conda, you should use

    conda install -c conda-forge cuda-toolkit

but note that you will still need to install the drivers separately if they are not
available already.

endef

endif    # ENV_TYPE is conda

endif    # CUDA not detected

DISABLED_PACKAGES_MESSAGE := $(subst $(newline),\n,$(subst ','\'',$(DISABLED_PACKAGES_MESSAGE)))


define INSTALL_HELP :=

$(BOLD)* (4/4) Install AMUSE *$(END_BOLD)

To install the AMUSE framework and all of the enabled packages into the active
$(ENV_TYPE) environment $(BOLD)$(ENV_NAME)$(END_BOLD), type

    make install


To install only the framework, you can use

    make install-amuse-framework

To install a specific package, you can use

    make install-amuse-bhtree

or whichever package you want to install, as long as it's enabled.

To do a development install of a specific package, you can use

    make develop-amuse-bhtree

or whichever package you want to work on.

endef

INSTALL_HELP := $(subst $(newline),\n,$(subst ','\'',$(INSTALL_HELP)))

