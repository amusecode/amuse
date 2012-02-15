==============
Coupling Codes
==============

The design for coupling codes in AMUSE is based on providing the same set of functions for
every community code and using these to devise different coupling methods. As the 
coupling methods are not fixed and can change on a per problem basis the functions
to be very generic. The AMUSE library defines three sets of functions to support coupling
codes:

particle or gridpoint manipulation
    Most properties of particles (or gridpoints) can be queried and updated during the run, 
    providing a direct method of manipulating the data of a community code. Further most
    codes support removing and adding particles during the run.

stopping conditions
    Stopping conditions are designed to interrupt a code during model evolutions. Stopping
    conditions are triggered when a code encounters a predefined state (for example
    a particle escaping out of the bounding box). 
    
services
    Services are functions added to a code that use the model of that code to provide
    data for other codes. For example a smoothed particle hydrodynamics code can provide
    the state of the model at any random point (not just on the particles) which can be
    used to create a grid from an particle model.
    
    