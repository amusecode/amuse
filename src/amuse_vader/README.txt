VADER, the Viscous Accretion Disk Evolution Resource, is a code to simulate the
evolution of viscous, thin, axisymmetric accretion disks. It has been developed
by Mark Krumholz and John Forbes [1]. For details on the code itself, please refer
to the documentation on the developers' bitbucket page [2]. This document is
concerned with the details of the AMUSE interface.

VADER has a feature where the user can define a number of functions that describe
a variety of parameters, such as the equation state parameters, source terms, and
boundary conditions. These are implemented in c, and as such VADER must be recompiled
to include new modules or edits to existing ones. The steps are as follows:

1:  Write your own functions. These must be in a file with name userFunc_{PROB}.c,
    where {PROB} is the name of your problem/set of functions, which must be placed
    in the src/prob directory of the interface. It is convenient to copy
    userFunc_none.c and edit that. Note that these functions will only be used if the
    corresponding parameters are set to true (e.g. alpha_function for userAlpha).

2:  User-defined parameters are passed to the userFunc through the params pointer,
    which is of unknown (user-definable) type. In the interface, we have fixed this
    to an array of floats of user-definable length. In an AMUSE script, the number
    of parameters must first be defined via the number_of_user_parameters parameter,
    and then each parameter is assigned by calling set_parameter(i, X), where i is
    the parameter index (starting at 0) and X is the parameter itself. AMUSE is
    agnostic as to the units of these parameters, they must be passed without unit
    and as such unit conversion must be handled by the user. 

3:  User-defined outputs (userOut in VADER) can be accessed in the interface through
    the grid_user grid. This possesses a single quantity, value, which has shape
    (nUserOut, nCells), and does not have AMUSE units (meaning that unit conversion
    is left to the user). The number of user outputs can be set through the
    number_of_user_outputs parameter. Note that currently, the grid_user must be 
    resized before it is filled.

4:  Recompiling the code is done in the root directory of AMUSE, using the following
    command:
    - make vader.code
    This automatically compiles all problems.



In order to run a simulation use the following steps:

1:  Initialize the code. The interface passes all quantities in cgs units; if 
    working in physical units, no converter is needed.	Vader has its own diagnostic
    outputs, which can be enabled via the verbosity parameter (0 through 3). Note
    that for this to show on terminal, the keyword "redirection='none'" must be
    passed. Additionally, a user-defined problem can be used by passing the
    "mode={PROB}" keyword.
    Example: viscous = Vader(redirection='none', mode='ring')

2:  Initialize a grid, either with a Keplerian, flat, or tabulated rotation curve.
    Free rotation curves are currently not avaiable through the interface.

3:  Assign parameters. Default values can be found in DOCU.txt.

4:  Assign column density and pressure distributions. Note that the pressure is
    actually vertically-integrated pressure, and has units of [pressure]*[length].
    Internal energy need only be assigned if equation_of_state_function is set to
    True.

5:  Evolve model as desired.

Details on the various parameters and (grid) properties can be found in DOCU.txt.

Note that VADER has the capability to use an algorithm called Anderson acceleration.
This is governed by a compiler parameter, AA_M (defined in src/Makefile). Whether it 
actually accelerates the code depends on the problem, as it involves additional 
overhead. This is characterized in the code paper [1], section 4. In the interface 
author's case, it was fastest to run without acceleration (AA_M=0), and this is the 
default. However, Krumholz & Forbes find that some problem only converge for AA_M>0, 
especially if backwards Euler integration is used instead of Crank-Nicholson.

[1] Krumholz, M. R. and Forbes, J. C., VADER: A Flexible, Robust, Open-Source Code
for Simulating Viscous Thin Accretion Disks, Astronomy and Computing, Vol. 11 (2015)

[2] https://bitbucket.org/krumholz/vader/src/master/
