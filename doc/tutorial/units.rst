.. _working_with_units:

==================
Working with Units
==================

The AMUSE framework provides a unit handling library. This library
is used throughout the AMUSE framework. When interacting with a code
all data has a unit, even scaled systems have units.

==========
Quantities
==========

The basic data object is a quantity, a quantity is made up of 
a value and a unit. The value can be a single number (a scalar 
quantity) or a multi-dimensional array (a vector quantity). 

Quantities are created by the ```|(bar)``` operator. All 
quantities can be used like numbers and a lot of numpy 
functions also work on quantities.

.. code-block:: python
    
    >>> from amuse.units.si import *
    >>> from amuse.units.core import named_unit
    >>>
    >>> weigth = 80 | kg
    >>> persons = 10
    >>> print "Total weight: ", persons * weigth
    Total weight:  800 kg
    
    >>> day = named_unit("day", "d", s * 60 * 60 * 24 )
    >>> weight_loss = (0.1 | kg) / (1 | day)
    >>> print "Weight loss: ", weight_loss
    Weight loss:  0.1 1.15740740741e-05 * kg * s**-1
    >>> print "Weight loss: ", weight_loss.as_quantity_in(kg/day)
    Weight loss:  0.1 kg / d
    
===================
Working with arrays
===================

A vector quantity can be used like a python list. Take care to only 
put quantities into a vector quantity.

.. code-block:: python
    
    >>> from amuse.units.units import MSun
    >>>
    >>> masses = [] | MSun
    >>> for i in range(10):
    ...     masses.append(i**2 | MSun)
    >>> print "Masses:", masses
    Masses: [0.0, 1.0, 4.0, 9.0, 16.0, 25.0, 36.0, 49.0, 64.0, 81.0] MSun
    
.. note::

    When working with arrays, some care must be taken to ensure that
    vector quantities are created and not arrays of quantities. The
    following code will create an array of quantities
    
    >>> from amuse.units.units import MSun
    >>>
    >>> masses = []
    >>> for i in range(2):
    ...     masses.append(i**2 | MSun)
    >>> print "Masses:", masses
    Masses: [quantity<0 MSun>, quantity<1 MSun>]

    
