=================================================================
Estimating the typical mass of the most massive star in a cluster
=================================================================

In this tutorial we will estimate the typical mass of the most 
massive star for a cluster with a Salpeter initial mass function (IMF).
This tutorial also illustrates that the units won't get in the way of 
calculations.

First we import numpy and set the seed of its random number 
generator (RNG). All random numbers used in AMUSE are drawn from this 
RNG. Seeding it is of course not necessary, but it will make the 
results reproducible.

.. code-block:: python

    >>> import numpy
    >>> numpy.random.seed(123456)

We will also need units and the Salpeter IMF generator. The easiest 
way to import almost everything from the AMUSE toolbox is to do:

.. code-block:: python

    >>> from amuse.lab import *

But here we will import the required module and function manually for clarity:

.. code-block:: python

    >>> from amuse.units import units
    >>> from amuse.ext.salpeter import new_salpeter_mass_distribution

Now, we can calculate the maximum stellar mass for each of 
*number_of_cluster_realizations* realizations of an 
*number_of_stars* star cluster. The Salpeter IMF runs from 0.1 to 
125 solar mass with a slope of -2.35, by default. Suppose the 
maximum of 125 solar mass is a bit too high for our taste, so we set 
it to 100 solar mass.

.. code-block:: python

    >>> number_of_stars = 1000
    >>> number_of_cluster_realizations = 100
    >>> maxmasses = [] | units.MSun
    >>> for i in range(number_of_cluster_realizations):
    ...     maxmasses.append(max(new_salpeter_mass_distribution(
    ...         number_of_stars, 
    ...         mass_max = 100. | units.MSun
    ...     )))
    ...

Note that the way we initialize *maxmasses*, with the solar mass 
unit, forces it to be a VectorQuantity instead of a Python list of 
ScalarQuantities. This is always recommended, because it is much 
faster, and it will make sure that AMUSE always recognizes it as a 
Quantity.  

If we want to know the mean mass of the *maxmasses* VectorQuantity, 
we simply use the *mean* function of VectorQuantities:

.. code-block:: python

    >>> print "mean:  ", maxmasses.mean()
    mean:   27.4915750164 MSun

The same works for the median or the standard deviation of *maxmasses*.

.. code-block:: python

    >>> print "median:", maxmasses.median()
    >>> print "stddev:", maxmasses.std()
    median: 21.0983403429 MSun
    stddev: 19.7149800906 MSun

Slightly slower, but giving the same result, we can use the numpy functions:

.. code-block:: python

    >>> print "mean:  ", numpy.mean( maxmasses)
    >>> print "median:", numpy.median(maxmasses)
    >>> print "stddev:", numpy.std(maxmasses)
    mean:   27.4915750164 MSun
    median: 21.0983403429 1.98892e+30 * kg
    stddev: 19.7149800906 MSun

Something weird has happened to the unit of the median 
mass. The result is still correct but the unit is converted to SI 
units. This is usually caused by a multiplication of a Quantity, 
where AMUSE tries to simplify the result, cancelling out for example 
factors of kg / kg. There's no need to bother, but if it annoys you, 
it can easily be fixed by:

.. code-block:: python

    >>> print "median:", numpy.median(maxmasses).in_(units.MSun)
    median: 21.0983403429 MSun



