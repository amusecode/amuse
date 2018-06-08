===============
Obtaining AMUSE
===============

Download
--------

Go to the `Amusecode github <https://github.com/amusecode/amuse>`_ page.

Getting started
---------------

The first step in getting AMUSE to work is obtaining the AMUSE
source code. We advice you to do this even before installation of the 
prerequisite software (:ref:`prerequisite-label`). In the following 
installation instructions we assume that you will install AMUSE in a 
directory ``/amuse``. 

Releases
--------

For the official releases we provide tarballs via https://github.com/amusecode/amuse/releases.

Tarball
~~~~~~~

Obtain the tarball (e.g. release-11.2.tar.gz) from the download site and unpack it 
in the amuse directory using:

.. code-block:: sh

    > tar -xf release-11.2.tar.gz

this will make an amuse sub-directory ``amuse-release-11.2``, which we will be referring to as
the AMUSE root directory.

From here proceed by reading the  :ref:`prerequisite-label` section.

Bleeding edge
-------------

The current development version is available via git repository access 
by issuing the following command:

.. code-block:: sh

    > git clone https://github.com/amusecode/amuse.git

This will make an AMUSE root directory with the name "amuse".  
