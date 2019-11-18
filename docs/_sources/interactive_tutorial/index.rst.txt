Interactive tutorial
====================

This is the static description of the interactive tutorial for AMUSE.

The interactive tutorial is always delivered with the binary distribution
and you can start it in that distribution with:

.. code-block:: sh

    ./amuse-tutorial
    
For the source distribution you can also run the interactive tutorial 
but you will need to install IPython notebook first (`ipythonnb`_). 
Once you have installed ipython you can run the notebook with:

.. code-block:: sh

    export PYTHONPATH=$PWD/src:$PYTHONPATH
    
    cd doc/interactive_tutorial/
    
    ipython notebook --pylab inline 

.. note::
    
    You will need to run the ipython command from the 
    interactive tutorial directory.
    You also need to have amuse in your python path.

.. note::
    
    The ipython notebook code depends on the `tornado` and the `pyzmq` 
    modules. You will need to install these before 
    installing the ipython notebook.
    
.. toctree::
    :maxdepth: 1

    01-Loading_AMUSE
    02-Quantities_with_units_
    03-Generic_units
    04-Collections_of_Particles
    05-Attributes_and_functions_on_particle_collections
    06-Using_a_community_code
    07-Channels
    08-Grids


.. _ipythonnb: http://ipython.org/ipython-doc/dev/interactive/htmlnotebook.html

