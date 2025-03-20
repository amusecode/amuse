===============
SecularMultiple
===============

A code to compute the secular (orbit-averaged) gravitational dynamics of hierarchical multiple systems composed of nested binary orbits (simplex-type systems) with any configuration and any number of bodies. A particle can represent a binary (`is_binary = True') or a body (`is_binary = False'). The structure of the system is determined by linking to other particles with the attributes child1 and child2. Tidal interactions and relativistic corrections are included in an ad hoc fashion (tides: treating the companion as a single body, even if it is not; relativistic terms: only including binary-binary interactions).

Be aware of the following potential pitfalls.
- The code is based on an expansion of the Hamiltonian of the system, assuming the latter is hierarchical. For an inner and outer orbit, the pericenter of the outer orbit, a_out(1-e_out), should at least be larger than the apocenter of the inner orbit, a_in(1+e_out). Generally, the expansion becomes a better approximation as the ratio [a_in(1+e_out)]/[a_out(1-e_out)] decreases.
- Suborbital effects are not modeled because of orbit averaging. For example, the code does not model mean motion resonances or evection resonances.
- The code implicitly always assumes dynamical stability of the system. The user should always check for dynamical stability. Note that a system which is initially dynamically stable, can become unstable because of secular evolution. Root-finding functions can be used in this case, and are based on several popular analytic  criteria (notably, the Mardling & Aarseth 2001 criterion). 
- Internally, the code uses orbital vector to evolve the equations of motion. Externally, the user will typically input orbital elements, which are updated after running the code. NOTE: currently, the orbital conversion routines may not work properly when inputting zero orbital elements (e.g. zero eccentricity and/or zero inclination). To avoid problems, please enter tiny numbers instead (e.g. 1e-10).
- The relative accuracy of all ODE variable is specified with the code parameter `relative_tolerance’, which is by default 1e-16 (rather conservative). NOTE: the absolute tolerances are hard-coded in evolve.cpp; their values should be okay for most planetary/stellar configurations. Problems may arise if the absolute tolerances are not appropriately set; the user should modify evolve.cpp in this case.
- When the secular time-scales are short and the user-specified time-step is long, the number of secular oscillations can be extremely large such that the desired accuracy cannot be obtained using default precision. This issue can be overcome by reducing the user-specified time-step, or by adjusting the tolerances.

If you use this code for scientific publications, please refer to 2016MNRAS.459.2827H.

Adrian Hamers, January 2017