==================
AMUSE Style Guide
==================


This Style Guide covers the source code written in the AMUSE project. The
source of the existing codes, integrated into AMUSE, do not have to
follow this guide. Existing codes are not rewritten.

Python Code
-----------
All python code should be consistent with the 
`Python Style Guide <https://www.python.org/dev/peps/pep-0008/>`_.
We have defined some deviations from the python style guide. 
These are listed in this document.

Maximum Line Length
    Limit *most* lines to a maximum of 79 characters.
    
    Some lines may be longer than 79 characters, especially when
    this increased readability. When you need to define a class with
    a long name, having a superclass with a long name, a line with
    more than 79 characters is OK.
    
    For example::
    
        class NeedsALongNameToDescribeItsFunction (AComplicatedNameForTheSuperClass):
        
    Is better than::
    
        class NeedsALongNameToDescribeItsFunction \
            (AComplicatedNameForTheSuperClass):
    
    Or::
    
        class NeedsALongNameToDescribeItsFunction (
            AComplicatedNameForTheSuperClass):
    
    
Naming Conventions
    Naming the classes and methods in the code is very important. However,
    do not spend too much time on naming. Use long names. Use a longer
    name when a shorter name cannot be determined quickly. When you've 
    defined a very long name, it's use may be harder. When it is used 
    a lot, a refactoring to shorter name will present itself. When it is
    not used a lot, at least the class or method is described better
    in the longer name.
    
    Do not use abbreviations in names. Especially, do not use the 
    first three or four letters of a word (for example, *dir* for directory
    or *ref* for reference). Also, do not remove the vowels to shorten a
    name (for example *tmp* for temp).
    
    The following one letter variables are allowed.
    
     * Loop variables, when used as ``for x in ...``,
       names allowed are 'x' and 'y'. Example::
       
            for x in stars:
                print x.mass
        
     * Loop index variables, when used as ``for i in range(10):``, 
       only name allowed is 'i'. Thus when writing a loop in a loop, you need
       to use longer names.
       For example, to iterate over every item in a 10 by 5 table, do:s:
            
            for row in range(10):
                for column in range(5):
                    print table[row][column]
        
    These names are only allowed in short loops of maximum 6 lines. Longer
    loops need better names. Also, for looping over row indices and column
    indices use *row* and *column*.
    
    
C / C++ code
------------

All C and C++ code should be consistent with the 
`Google C++ style guid <https://google.github.io/styleguide/cppguide.html>`_

File Names
    All C++ files should have a ``.cc`` suffix.
     
Naming Conventions
    We follow the Python style for naming function names. Function names \
    are all lowercase, with underscores between words.
    

Fortran code
------------
All Fortran code written in AMUSE should be Fortran 90.  We will follow
the coding standard given in 
`Fortran Coding Standard for the Community Atmospheric Model <http://www.cesm.ucar.edu/working_groups/Software/dev_guide/dev_guide/node7.html>`_

File Names
    All fortran 90 files should have a ``.f90`` suffix.
    
Layout
    Source code should follow the "free form" fortran 90. Do not use
    the "fixed form" fortran 77 layout.
    
Indentation
    Use indentation to highlight the block structure of your
    code::
    
        FUNCTION example(arg1)
            REAL  :: arg1, example
            example = arg1 * 2.0
        END FUNCTION
        
    Use 4 character to indent the code.
    

    
    
    

    
    
