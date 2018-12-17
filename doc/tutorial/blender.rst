
===========================================
Using Blender for visualising Amuse results
===========================================

.. image:: blender.jpeg
   :width: 2cm

Blender
=======

Blender is a free open source 3D content creation suite, *GPL-ed*. As it is specialised in visualisation of 3D content *and* 
scriptable in python it provides handy a tool for visualisation of Amuse data.

More information on Blender can be found on https://www.blender.org/

Starting Blender
================
We start blender from the command line, assuming all environment variables are set for amuse, and make sure mpd is running.

.. code-block:: bash

    >>blender &

After blender opened we press **CTRL-RIGHTCURSOR** three times to open a text panel next to our 3D view. In the text panel we use the Text 
menu item to open a text file:

.. image:: blenderopenscript.png
   :width: 800
   :align: center

This way we can retrieve scripts we had written before. Now, we are going to write a new script, but use our favourite editor, and once finished 
we will load the script using the method described above.

**Note:** By default, blender places a cube in the scene. You can delete it else it will obscure the sun in our example. 

The script is based on the sun-earth test in Hermite. We start with that amuse script and add a couple of lines to communicate with blender, which 
are marked by comments in the example below.

Amuse blender API
=================

To simplify communication with blender, amuse has the module amuse.ext.blender, which contains a very basic API. In our example we start with importing 
it on top of the regular stellar dynamics imports:

.. code-block:: python

    from amuse.ext.blender import blender 

To create a sphere, e.g. the sun,  with 32 segments, 32 rings and a radius of 1.0 in the current scene we type:

.. code-block:: python

    sun = blender.Primitives.sphere(segments = 32, rings = 32, radius = 1.0)

We can move and rotate our object using:

.. code-block:: python

    x,y,z = 1.0, 0.0, 0.0
    alpha = 0.2
    sun.loc = (x, y, z)
    sun.rotZ = alpha


The Hermite Sun-Earth model with blender visualisation will become like this:

.. literalinclude:: ../../examples/applications/test_blender.py

The path to this example code is {amusedir}/trunk/examples/applications/test_blender.py
      
We save this file as *myname.py*, open it in blender and run it by typing **ALT-P** (cursor in editor):

.. image:: blenderrunscript.png
   :width: 800
   :align: center 

