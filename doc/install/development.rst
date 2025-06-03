Installing for development
==========================

AMUSE provides a rich toolbox of astrophysical methodology, but it's in the nature of
science to try to do new things, and sometimes that requires modifying the methods too.
For computational science, that means changing the software. In case of AMUSE, if this
is what you want to do then you'll probably want to modify one of the codes included
with it, changing its C++ or Fortran code and recompiling it to try your modifications.

To do that, it's best to get a development copy of the AMUSE code using git. That will
give you a local git repository so that you can easily track your changes, and
contribute them back to AMUSE (pretty please?) if they're potentially useful to others.

To get a development copy of the AMUSE source code, you'll need ``git``, which is
available from every popular package manager. In your Conda environment,

.. code-block:: bash

    conda install git

will get you sorted. Next, you can get a local git repository with AMUSE in it using

.. code-block:: bash

    git clone https://github.com/amusecode/amuse.git


This will create a directory called ``amuse`` with the latest development version of
AMUSE in it. You can ``cd`` into that and access the installer using ``./setup`` as
before.

If you have an existing installation that's not too old, then you can probably just
switch the code you're interested in to development mode using the installer in the new
source directory:

.. code-block:: bash

    ./setup develop amuse-bhtree


This will uninstall the existing package (we'll use bhtree in this example, but you can
use any of them), and replace it with a development install. Development installs are
special because rather than copying the code into the environment, they link it. As a
result, you can try your changes right away, without having to reinstall every time.
That's really convenient.

Once you have a development installation, you can go to the code's source:

.. code-block:: bash

    cd src/amuse_bhtree

and have a look around. You'll find an ``interface.py`` that defines the code's
interface with AMUSE, the ``tests/`` directory with tests, and a directory ``src/`` with
the source code. There's also a ``bhtree_worker``, or perhaps even multiple. This file
is created from the ``bhtree`` source code together with its AMUSE wrapper using

.. code-block:: bash

   make bhtree_worker

To modify the code, open a file in ``src/`` and change it to your liking, then run the
``make`` command above in the ``src/amuse_bhtree`` directory to recompile the worker.
Then you can run your Python script (in the same environment of course!) to try it out
and see if it works.

