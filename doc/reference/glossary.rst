

Glossary
--------

.. glossary::
    :sorted:

    Physical Domain
      A region in space, time and objects (stars / gas clouds) that can be
      described by a set of physical models. A physical domain is a
      limited description of the physics of a stellar system. For example
      in the gravity domain the gravitational forces are taken into
      account but not the evolution of stars or the hydrodynamics of the
      interstellar medium.

    Module
      Interface to a code in AMUSE. Also, file containing python definitions.
      These definitions can be imported into the main Python program.

    Unit test
      Automatic test case for a part of the code. For example one can
      define a unit test that checks the result of one method given a
      specified argument.

    Application
      A script that combines modules into a multi-physical simulator.

    Run
      During a run the physical models are evolved to a wanted end state.
      A run can end when a certain model time has passed. For example, one
      can run a script until the age is 12 million years.

    Evolve
      Work out / develop the physical model equations. Usually, small steps
      in time are taken, the new state of the system is calculated from the
      previous state.

    Step
      An identifiable stage in a run or evolve. A code can take multiple
      steps during one evolve step. A step can be the integrator step,
      force calculation step.

    Data-model
      A single representation model for data of different physical codes.

    Set
      A collection of identifiable objects (for example all stars in a
      globular cluster). Every value in the set can be seen as belonging
      to an object.

    Grid
      A collection of identifyable locations or regions (for example a box
      of interstellar medium with a density, molecular composition). Every
      value in the collection can be seen as belonging to a specific
      location (x,y,z point).

    Code
      A computer program to solve or approximate the equations of state of
      a physical domain. Codes are usually written in Fortran or C/C++.
      For example 'sse' is a code to perform stellar evolution.

    Production Code
      A code used for modelling stellar systems and described in published
      articles.

    Toy Code
      A code with limited accuracy but that is easy to use.

    Numerical Experiment
      A numerical simulation designed to investigate a physical system.

    SPH, smoothed particle hydrodynamics
      A computational method used for simulating fluid flows.

    AMR, adaptive mesh refinement
      A numerical technique that dynamically refines the computational mesh
      in regions requiring higher spatial resolution.

    Input / output format
      The format used to read input data or write simulation results.

    Binary star
      Two stars in bound orbit.

    Star Cluster
      A group of stars. The mean distance between the stars is smaller
      than the mean distance between stars in the galaxy disks. Two types
      of star clusters exist: globular clusters, open clusters

    Dense stellar systems
      Bound conglomerate of stars with at least one 90º turnaround within a crossing time.

    Interstellar Medium
      The matter between the stars in a galaxy.

    Interstellar Cloud
      Generic name for an accumulation of gas, plasma and dust. A denser
      region of the interstellar medium.

    Molecular Cloud
      A type of interstellar cloud whose density and size permits the formation of molecules.

    Gravity Domain
      A physical domain in which gravitational interactions are taken into account.

    Hydrodynamics
      The study of the motion and behaviour of fluids, under the influence of forces such as pressure and gravity.

    Star
      Luminous spheroid of plasma held together by self-gravity.

    Stellar Evolution
      The physical processes governing the structural and chemical evolution of a star over time.

    ZAMS, zero age main sequence
      Start of the main sequence in the evolution of a star. Stable stage
      after formation of the star from a collapsing gas cloud.

    Main sequence star
      A star that derives its energy from the conversion of hydrogen into
      helium in its core.

    Metallicity
      Proportion of star matter made up of chemical elements other than
      hydrogen and helium.

    Galaxy
      A gravitationally bound system consisting of stars, stellar remnants,
      interstellar gas and dust, and dark matter.

    Star Core
      The central region of a star.

    Star Envelope
      The outer layers of a star surrounding its core.
