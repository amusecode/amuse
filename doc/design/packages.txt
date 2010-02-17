===============
Python Packages
===============

Like all large python projects, the AMUSE source-code is structured 
in packages, subpackages and modules. Packages and subpackages are
directories in the filesystem and modules are python files (A detailed
explanation can be found in 
`Modules <http://docs.python.org/tutorial/modules.html>`). In this
chapter the most important packages are named. 

The sourcecode of the **AMUSE Code** and **Legacy Codes** layer is 
combined in one package. The package name is **amuse**.

amuse
    Root package of the AMUSE code, does not contain test files,
    the build system or the test system.

The **amuse** and the package is further divided into three
subpackages:

amuse.support
    Contains the code of the **AMUSE Code** layer. The units system,
    data io and model and all base classes can be found in 
    this package.
    
amuse.legacy
    Contains the code the **Legacy Codes** layer. The legacy codes
    can be found in this package as well as support code for generating
    the script to legacy code messaging framework.
    
amuse.ext
    Contains **extra** and/or **extension** code. For example, making
    initial data models is not one of the main functionalities of AMUSE,
    but it is very useful to include this into the codebase.

