.. _cuda-setup-label:

==========================
Setting up GPGPU with CUDA
==========================

Introduction
~~~~~~~~~~~~

Here we provide help for setting up general-purpose computing on graphics processing units (GPGPU)
using CUDA. Performing (part of) the calculations on a graphics card can result 
in a significant speed-up. Several codes in AMUSE support or require GPGPU: 
phi-GRAPE (using Sapporo), ph4 (using Sapporo), Octgrav, Bonsai, HiGPUs.


Self-help script
~~~~~~~~~~~~~~~~

In the AMUSE root directory a self-help script can be found. If building or testing any of the 
codes mentioned above fails and you wonder why, it will hopefully provide you with helpful suggestions.
From a command line run the bash script `cuda_self_help`:

.. code-block:: sh

    > ./amuse-x.0/cuda_self_help


Step-by-step
~~~~~~~~~~~~

* :ref:`CUDA-capable`
* :ref:`CUDA_SDK`
* :ref:`env-vars`
* :ref:`configure-with-cuda`
* :ref:`test`


.. _CUDA-capable:

Check that your computer has a CUDA-capable Nvidia graphics card
-----------------------------------------------------------------

First determine the model of your GPU.

On Linux:

.. code-block:: sh

    > nvidia-settings -q gpus


On Mac:

   1. Click on “Apple Menu”
   2. Click on “About this Mac”
   3. Click on “More Info”
   4. Select “Graphics/Displays” under Contents list

Check whether your GPU model is listed among 
`Nvidia's CUDA-enabled GPUs <https://www.nvidia.com/object/cuda_gpus>`_.


.. _CUDA_SDK:

Check that you have installed the CUDA Toolkit (TK) and software development kit (SDK)
--------------------------------------------------------------------------------------

If not, download and install it from `CUDA Development Tools <https://developer.nvidia.com/cuda-downloads>`_.


.. _env-vars:

Set the CUDA_TK and CUDA_SDK environment variables
--------------------------------------------------

After installing the CUDA TK and SDK, make sure the environment variables CUDA_TK and CUDA_SDK are set correctly.
For shell (bash) you need to do:

.. code-block:: sh

   export CUDA_TK=/path/to/cuda_tk
   export CUDA_SDK=/path/to/cuda_sdk

'/path/to/cuda_tk' should hold directories named 'include', 'lib', and 'lib64' (where libcudart.so is located)

'/path/to/cuda_sdk' should hold a directory named 'common/inc' (where various header files are located)

We recommend you add these lines to your '.bashrc' file so that
the variables are set correctly for all sessions. If you have a
C shell you need to do a *setenv* and edit the '.cshrc file.


.. _configure-with-cuda:

Configure AMUSE with CUDA enabled
---------------------------------

AMUSE needs to be configured with the option ``--enable-cuda``. See :ref:`configuration-gpu-label`.


.. _test:

Testing
-------

Now try building for example Octgrav and run the nosetests (from AMUSE root directory),
but first re-initialize mpd (or it will remember its original environment):


.. code-block:: sh

   mpdallexit
   mpd &
   make octgrav.code
   nosetests ./test/codes_tests/test_octgrav.py

If this fails, please contact us through the `'amusecode' google group <http://groups.google.com/group/amusecode>`_, 
or on IRC at the #amuse channel on irc.freenode.net. 
