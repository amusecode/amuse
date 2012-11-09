==========================
Supported File Formats
==========================

.. automodule:: amuse.io

    Introduction
    ------------
    The AMUSE framework provides a generic way of reading and writing
    sets of entities (these entities can be particles, gasclouds or
    gridpoints). AMUSE provides a function to write a set to a file
    and a function to read a set from a file. These functions need
    the name of the file and the format of the data in the file. We
    will describe these functions first in this chapter. The functions
    can throw 3 kinds of exceptions, these are described next.
    
    Use
    ---
    
    To write a data set to a space separated text file do:: 
        
        >>> from amuse.io import write_set_to_file
        >>> from amuse.datamodel.core import Particles
        >>> from amuse.units import units
        >>> x = Particles(100)
        >>> x.mass = 10.0 | units.MSun
        >>> x.radius = range(1.0,101.0) | units.RSun
        >>> write_set_to_file(x, "test.csv","txt", attribute_types = [units.MSun, units.RSun])
        
    .. autofunction:: write_set_to_file
    .. autofunction:: read_set_from_file
    .. autofunction:: get_options_for_format
    
    Exceptions
    ----------
    
    .. autoclass:: UnsupportedFormatException
    .. autoclass:: CannotSaveException
    .. autoclass:: CannotLoadException
    
    
    
    Amuse HDF5 files
    ----------------
    
    AMUSE provides it's own data format based on hdf5 files. The HDF5 file storage
    provides support for storing particle set and grids in the same file. It also
    retains the particle references and can store multiple versions of the
    same set.
    
        >>> from amuse.support import io
        >>> from amuse.ic.plummer import new_plummer_model
        >>> particles = new_plummer_model(1000)
        >>> io.write_set_to_file(
                    particles,
                   'plummer.hdf5',
                   'hdf5',
            )
        >>> io.read_set_from_file('plummer.hdf5', 'hdf5')
    
    
    .. iooptions:: hdf5
    
    
    
    Text files
    ----------
    
    AMUSE has support for reading and writing text (``.txt``, ``.csv``) files.
    You can specify the txt format by entering "txt" (for space separated files) 
    or "csv" (for comma separated files)
    as the format::
    
        >>> from amuse.support import io
        >>> particles = io.read_set_from_file(
                   'plummer.txt',
                   'txt',
                   attribute_types = (units.m, units.m, units.kms),
                   attribute_names= ('x', 'y', 'vx')
            )
        >>> io.write_set_to_file(particles, 'plummer.txt', 'txt')
        
    .. iooptions:: txt
    
    By default, text files are stored with some unit information on 
    a comment line in the header. Unfortunately, amuse cannot depend 
    on this line to read the file. When reading a text file you 
    always need to specify the units of every column in the file. For
    example, to read a file with particle positions where all positions
    are stored in parsec, do::
    
        particles = io.read_set_from_file(
               'positions.txt',
               'txt',
               attribute_types = (units.parsec, units.parsec, units.parsec),
               attribute_names= ('x', 'y', 'z')
        )
    
    When writing a text file, you do not need to specify the units. If
    you don't specify the units the ``write_set_to_file`` function will
    default to the units of the quantities being saved. For reliability 
    and reproducibility we suggest to always specify the units and
    the names of the attributes to save when saving a text file::
        
        particles = io.write_set_from_file(
               particles,
               'txt',
               attribute_types = (units.parsec, units.parsec, units.parsec),
               attribute_names= ('x', 'y', 'z')
        )
        
    The amuse text file routines are very generic and since version 6.0,
    you can also specify a list of vector quantities instead of a particle
    set::
    
        particles = io.write_set_from_file(
               None,
               'txt',
               attribute_types = (units.parsec, units.parsec, units.parsec),
               attribute_names = ('x', 'y', 'z'),
               quantities = [
                    [0.1,0.2] | units.parsec, 
                    [0.3,0.4] | units.parsec, 
                    [0.5,0.6] | units.parsec, 
               ]
        )
    
    Starlab
    -------
    
    AMUSE has support for reading and writing starlab (``.dyn``) files.
    You can specify the starlab format by entering "dyn" or "starlab"
    as the format::
    
        >>> from amuse.support import io
        >>> particles = io.read_set_from_file('plummer.dyn','starlab')
        >>> io.write_set_to_file(particles, 'output.dyn', 'dyn')
    
    The starlab format support sevaral options, listed below. You can
    use these options by adding additional keyword arguments to
    the :func:`read_set_from_file` or :func:`write_set_to_file` functions. 
    For example::
    
        >>> from amuse.support import io
        >>> particles = io.read_set_from_file('plummer.dyn','starlab', must_scale = False, return_children = False)
    
    .. iooptions:: starlab

    The units of the values in the star (stellar properties)
    section of a starlab file are always in derived S.I units (solar
    mass, million years, solar luminocity etc.).
    
    The masses given in de the dynamics section of a starlab file
    are usually in *nbody* units. Some starlab tools will set the
    mass values in Solar mass units (for example the ``makemass``
    tool will return the masses in solar mass units). To read these
    files you need to set the ``dynamics_mass_units``.
    
    .. code-block:: bash
    
        > makeplummer -n 1000 > plummer1.dyn
        > cat plummer1.dyn | makemass -f 1 -x -2.0 -l 0.1 -u 20 > plummer2.dyn
        > cat plummer2.dyn | add_star -Q 0.5 -R 5 > plummer3.dyn
        > cat plummer3.dyn | scale -s > plummer4.dyn
        > cat plummer4.dyn | kira -S > plummer5.dyn
        
    The ``plummer1.dyn``, ``plummer4.dyn`` and ``plummer5.dyn`` files 
    will provide masses (and all other dynamical properties) in scaled
    *nbody* units. The ``plummer2.dyn`` and ``plummer3.dyn`` files
    will have masses in solar masses. To read each file in AMUSE, and
    return the particles with S.I. units, you need to do::
    
    
        >>> from amuse.support import io
        >>> from amuse.units import nbody_system, units
        >>> converter = nbody_system.nbody_to_si(1 | units.MSun, 1 | units.parsec)
        >>> particles1 = io.read_set_from_file('plummer1.dyn','starlab', nbody_to_si_converter = converter)
        >>> particles2 = io.read_set_from_file('plummer2.dyn','starlab', dynamics_mass_units = units.MSun, nbody_to_si_converter = converter)
        >>> particles3 = io.read_set_from_file('plummer3.dyn','starlab', dynamics_mass_units = units.MSun, nbody_to_si_converter = converter)
        >>> particles4 = io.read_set_from_file('plummer4.dyn','starlab')
        >>> particles5 = io.read_set_from_file('plummer5.dyn','starlab')
    
    .. note::
        
        No nbody converter object is needed for the last files, as the scale 
        factors given in the files will be used.
    
    The ``plummer1.dyn``, ``plummer4.dyn`` and ``plummer5.dyn`` can also
    be read in nbody units. In the following example the returned 
    particles have dynamic attributes (mass, radius, velocity, acceleration)
    in *nbody* units::
    
        >>> from amuse.support import io
        >>> particles1 = io.read_set_from_file('plummer1.dyn','starlab')
        >>> particles4 = io.read_set_from_file('plummer4.dyn','starlab', must_scale = False)
        >>> particles5 = io.read_set_from_file('plummer5.dyn','starlab', must_scale = False)
    
    
    NEMO
    ----
    
    AMUSE has support for reading and writing nemo (``.tsf``) files.
    You can specify the starlab format by entering "nemo" or "tsf"
    as the format::
    
        >>> from amuse.support import io
        >>> particles = io.read_set_from_file('plummer.tsf','nemo')
        >>> io.write_set_to_file(particles, 'output.tsf', 'tsf')
        
    The nemo format support several options, listed below. You can
    use these options by adding additional keyword arguments to
    the :func:`read_set_from_file` or :func:`write_set_to_file` functions. 
    For example::
    
        >>> from amuse.support import io
        >>> from amuse.units import nbody_system, units
        >>> converter = nbody_system.nbody_to_si(1 | units.MSun, 1 | units.parsec)
        >>> particles = io.read_set_from_file('plummer.nemo','tsf', nbody_to_si_converter = converter)
        
    .. iooptions:: nemo
        
        
    Gadget
    ------
    
    AMUSE has support for reading and writing gadget 2 files.
    You can specify the gadget format by entering "gadget"
    as the format::
    
        >>> from amuse.support import io
        >>> gas, halo, disk, bulge, stars, bndry = io.read_set_from_file('plummer.dat','gadget')
        >>> io.write_set_to_file(particles, 'output.dat', 'gadget')
    
    The gadget file format reader will return a tuple of gas, halo, disk, 
    bulge, stars and boundary particles. This is different form other
    readers which only return one set. The write_set_to_file can take
    a particles set (will be saved as halo particles) or a tuple of 
    gas, halo, disk, bulge, stars and boundary particles.
    
    The gadget format support several options, listed below. You can
    use these options by adding additional keyword arguments to
    the :func:`read_set_from_file` or :func:`write_set_to_file` functions. 
    For example (will read a gadget file with timestep information)::
    
        >>> from amuse.support import io
        >>> particles = io.read_set_from_file('plummer.nemo','gadget', has_timestep = True)
        
    .. iooptions:: gadget
        
