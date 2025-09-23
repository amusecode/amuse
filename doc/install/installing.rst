Installing AMUSE
================

To install AMUSE, we need to

- obtain a suitable Linux or macOS environment,
- install Conda, and
- install AMUSE.

(AMUSE can also be installed without Conda if needed, as described below, but we
recommend using Conda.)

If you have a Mac, then you should skip to :ref:`setting_up_macOS`, and if you're
running Linux then you can go straight to :ref:`installing_conda`. For Windows, continue
here with installing WSL.

Installing WSL
--------------

AMUSE does not run on Windows natively, so if you have a Windows computer then we need
an extra step to give it a Linux environment on top of Windows to work in. This is done
by installing Windows Subsystem for Linux, or WSL, which is done as follows:

- Open the Microsoft Store
- Install Windows Subsystem for Linux
- Restart your computer

- Open the Microsoft Store again
- Install Ubuntu

You can then open Ubuntu from the store or from the Windows menu. When you do that,
you'll see a terminal window in which you can type commands.

Setting up WSL
``````````````

The first time you do this, Ubuntu will ask you to set a username and a password.
Remember these well (or better store them in your password manager), because you'll need
them from time to time!

Once done, you will see a *prompt*, a bit of text ending with a dollar sign ``$``. This
is where you can type commands. It's a good idea to install the latest updates before
doing anything else. Type this after the prompt and press <Enter>.

.. code-block:: bash

   sudo apt update


This will ask you for your password, and then download the latest list of available
updates. We can then install them using this:

.. code-block:: bash

   sudo apt -y upgrade


Now you can continue with :ref:`installing_conda`, using the Ubuntu terminal window to
enter the instructions.

.. _setting_up_macOS:

Setting up macOS
----------------

If you have an Apple computer, then you're running macOS. AMUSE works on macOS directly,
but it does require setting up a development environment first.

Installing XCode
````````````````

The first step is to install the XCode Command Line Tools. To do that, open Terminal and
type this into the terminal window, then press return:

.. code-block:: sh

    xcode-select --install


This should make a pop-up appear, on which you can click *Install* to start the
installation. If you don't see a pop-up, click the Apple logo at the top left of your
screen, then select System Settings..., and you'll see that there's an update available
for "Command line tools for XCode". Go ahead and install the update to make the XCode
tools available.

Configuring the environment
```````````````````````````

Next, we need to make sure that AMUSE can find the files you just installed during
installation. To do that, edit the ``.zshrc`` file in your home folder and add

.. code-block:: sh

    export SDKROOT=/Library/Developer/CommandLineTools/SDKs/MacOSX.sdk


at the bottom. Then you should open a new terminal to ensure that that command has been
loaded, and then you're ready to install Conda.

If you use Bash instead of zsh, then you'll need to edit ``.bashrc`` instead. When in
doubt, you can safely edit both files to be sure.

.. _installing_conda:

Installing Conda
----------------

The next step to installing AMUSE is to install Conda, if you don't already have it
available. Conda is a package manager, a program with which you can install other
programs. It's very widely used in science and beyond, so having a working Conda setup
is very useful also outside of the world of AMUSE.

If you already have a working Conda setup, then you can continue to :ref:`installing_amuse`.

If you cannot or don't want to use Conda, see :ref:`using_a_virtualenv` below.

If you do not yet have Conda, then you can install it using the following commands in
the terminal. (Linux users can open one from the menu, Windows and macOS users will
already have one open at this point.)

To download the miniforge Conda installer, use this command:

.. code-block:: bash

   curl -L -O "https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-$(uname)-$(uname -m).sh"


You can then run the installer using

.. code-block:: bash

   bash Miniforge3-$(uname)-$(uname -m).sh


Finally, close your terminal window and open a new one to make the ``conda`` command
properly available.


.. _installing_amuse:

Installing AMUSE
----------------

To get a copy of the most recent release of AMUSE, go to AMUSE_Releases_ and look up the
most recent release. At the bottom of the description, you'll see a link **Source code
(tar.gz)**. Right-click that link and select "Copy link address", then use the ``curl``
command in your terminal to download it as above, for example:

.. code-block:: bash

   curl -L -O "https://github.com/amusecode/amuse/archive/refs/tags/v2025.9.0.tar.gz"


This ``.tar.gz`` file needs to be unpacked first (you may need to change the version if
you downloaded a newer one):

.. code-block:: bash

   tar xf v2025.9.0.tar.gz


Then we can enter the directory with the AMUSE source code:

.. code-block:: bash

   cd amuse-2025.9.0


And then you can start the installer:

.. code-block:: bash

   ./setup


From here on you can follow the instructions, using ``conda`` to create an environment
and install the dependencies.

Installing all of the AMUSE community codes will take a while. You may want to start
with just installing the framework, and install the codes as needed.

When the installer is done installing, you should have a working AMUSE setup.

If you encounter any problems, then you can ask for help in the `AMUSE
Slack <https://amusecode.slack.com>`_ or by `making an issue on
GitHub <https://github.com/amusecode/amuse/issues/new/choose>`_.


Additional packages
```````````````````

If you plan to follow the AMUSE tutorials then you'll need a few additional packages as
well. Fortunately, ``conda`` can help us here too:

.. code-block:: bash

   conda install scipy astropy jupyter pandas seaborn matplotlib


Fixing MPI on Ubuntu
````````````````````

On Ubuntu and WSL, there is a final command to run that fixes an issue with OpenMPI on
that system. This is not needed on macOS or other Linux versions.

.. code-block:: bash

   echo 'btl_tcp_if_include=lo' >>.openmpi/mca-params.conf


Using AMUSE in a new terminal
`````````````````````````````

If you close the terminal and/or want to continue working with AMUSE in a newly opened
one, then you'll first need to activate the Conda environment you made again:

.. code-block:: bash

    conda activate Amuse-env


Adding more codes
`````````````````

To access the installer, you need to enter the AMUSE source directory again

.. code-block:: bash

    cd amuse-2025.5.0


and then you can run it as before using

.. code-block:: bash

    ./setup


You should now have a working AMUSE setup. To start
using it, see :ref:`getting_started_with_amuse` or the :ref:`interactive_tutorial`


Debugging conda package installation
````````````````````````````````````

If you encounter problems with installing packages using ``conda``, or AMUSE doesn't
compile correctly, then you should check that you are using the ``conda-forge`` channel
rather than something else.

Conda can use different sources of packages, which it calls channels. Different channels
contain software packaged by different people, and packages from different channels are
often incompatible. If you type

.. code-block:: bash

    conda list


then you should see a list of packages that are installed in the active environment, and
which channel they came from. Ideally, all of them have ``conda-forge`` as the channel.

If not, then you can reinstall the package from ``conda-forge`` and see if that improves
the situation.

To reinstall a package from ``conda-forge``, use

.. code-block:: bash

    conda install -c conda-forge <package name>


If you want to combine AMUSE with another package that isn't available from conda-forge,
then you may have to install that from another channel, and hope that things work. Or
ask the maintainers of that package to add it to conda-forge and be a bit more
compatible with the rest of the world.


Alternative installation options
================================

The above instructions are the easiest way to install AMUSE, and they should work for
almost everyone wanting to use AMUSE to do astrophysics. Nevertheless, there may be
cases where you need a different setup, for example because you cannot use Conda. In
that case, you'll want one of these alternative installations.

.. _installing_from_git:

Installing from a Git repository
--------------------------------

If you plan to modify AMUSE or one of the codes in it, then you may want to install from
a local git clone instead of from a tar file. This will take more disk space and more
download time, so it shouldn't be the first option, but if you want to do it then you
can. You'll need to gave `git` installed:

.. code-block:: bash

   git clone https://github.com/amusecode/amuse.git


Then you can enter the source directory using:

.. code-block:: bash

   cd amuse


Select a version to build (use either one of these, or whichever version is relevant):

.. code-block:: bash

   git switch main                          # current development version
   git checkout checkout v2025.9.0          # tagged release

And now you can start the installer as before:

.. code-block:: bash

   ./setup


.. _using_a_virtualenv:

Using a virtualenv
------------------

In some cases, you may not want to or be able to use Conda to install AMUSE. In that
case, you can use a standard Python virtual environment (or venv for short) instead.
Unlike the `conda` command, the `pip` command that comes with virtual environments can
only install Python packages, which means that we need another package manager (such as
`apt` on Ubuntu or similar, `dnf` or Fedora or similar, or Homebrew or MacPorts on
macOS) to install the dependencies.

To install into a virtual environment, you can skip the instructions for installing
Conda (since it won't be used), and instead proceed straight away to installing AMUSE.
When the ``./setup`` command shows the instructions for making an enviroment, use the
ones for a virtual environment, and then ``./setup`` will guide you through installing
the dependencies using an appropriate external package manager and install AMUSE into
your virtual environment for you.

.. _AMUSE_Releases: https://github.com/amusecode/amuse/releases
