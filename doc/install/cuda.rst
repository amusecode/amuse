Installing CUDA
===============

If you have an nVidia GPU and are using a code with CUDA support, then compiling and
running with CUDA can give a nice performance increase. (Unfortunately, using GPUs for
computing rather than graphics is still poorly supported on other GPU brands, and the
codes in AMUSE can only use CUDA.)

Depending on your platform, installing CUDA can be a bit more complicated than the other
dependencies however. CUDA consists of two parts, the driver and the toolkit. The driver
must be installed on the system, as part of the nVidia GPU drivers. The toolkit can be
installed separately via an installer, the system package manager, or conda.

Ubuntu
------

Of the various Linux distributions, Ubuntu probably makes `installing CUDA the easiest
<https://documentation.ubuntu.com/server/how-to/graphics/install-nvidia-drivers/index.html>`_.
To see if any drivers are available for your computer, use:

.. code-block:: bash

    sudo ubuntu-drivers list --gpgpu

(The ``--gpgpu`` option tells the tool that we want to do calculations on the GPU rather
than just doing graphics, so it's important to add!)

If you have an nVidia GPU then this should list a number of driver versions. The tool
can automatically select the best one and install it using:

.. code-block:: bash

    sudo ubuntu-drivers install --gpgpu

After this, you'll need to install the CUDA toolkit, which has everything needed to
compile the codes to work with CUDA:

.. code-block:: bash

    sudo apt install nvidia-cuda-toolkit

After this, you'll need to restart the computer, and then CUDA should be available and
the installer should let you install CUDA packages..

Other Linux distributions
-------------------------

On other Linux distributions you may have to do a bit more work. The easiest way is
probably to first install the driver using either the instructions from your
distribution (if there are any) or the `ones from
nVidia <https://docs.nvidia.com/datacenter/tesla/driver-installation-guide/index.html>`_.

After that you can then install the toolkit using Conda:

.. code-block:: bash

    conda install -c conda-forge cuda-toolkit

nVidia also has a `comprehensive guide on installing CUDA on
Linux <https://docs.nvidia.com/cuda/cuda-installation-guide-linux/index.html>`_, which
makes it quite complicated but may work if you're willing to give it a go, aren't using
Ubuntu, and don't want to use Conda to install it.

Otherwise, try searching the web for some instructions, but make sure they're recent
because the best way to do it tends to change over time.


Windows
-------

There appears to be `support for using CUDA on Windows Subsystem for Linux
(WSL) <https://docs.nvidia.com/cuda/wsl-user-guide/index.html#getting-started-with-cuda-on-wsl-2>`_,
so if you're using WSL on a Windows computer with an nVidia GPU then you can give this a
try.

Like on Linux, you'll need both the nVidia driver and the CUDA toolkit. In this case,
you probably already have the Windows nVidia drivers installed, so there's no need to
install a driver. (Please don't try, the Linux drivers won't work on Windows!). Probably
the easiest way to get the toolkit is again via Conda:

.. code-block:: bash

    conda install -c conda-forge cuda-toolkit

or you can try the instructions in the link above.

