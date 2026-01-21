# -*- coding: ascii -*-
"""
Evolves a molecular cloud with chemical evolution

The dynamics is evolved in the sph code Fi, chemistry with Krome, assuming
an isothermal cloud (hence passive chemical evolution).

Initial condition is a smooth spherical cloud with random velocities as in
Bonnell et al. (2003)
"""
from __future__ import print_function
import numpy

from matplotlib import pyplot
from amuse.units import units, constants, nbody_system

from amuse.community.fi.interface import Fi
from amuse.community.krome.interface import Krome
from amuse.community.uclchem.interface import UCLchem

from amuse.ext.molecular_cloud import molecular_cloud
from amuse.io import write_set_to_file

gamma = 1.4
meanmwt = 1.35 | units.amu
ionization_rate=2.e-17 | units.s**-1

def get_n_T_xi(density, u):
    number_density=density/meanmwt
    temperature=((gamma - 1) * meanmwt * u / constants.kB)
    ionrate=ionization_rate*numpy.ones_like(density)
    radfield = 1.0|units.habing
    return (number_density, temperature, ionrate, radfield)

def update_chem(sph_parts, chem_parts):
    channel = sph_parts.new_channel_to(chem_parts)
    channel.transform(["number_density", "temperature", "ionrate","radfield"], get_n_T_xi, ["density", "u"])

def evolve_sph_with_chemistry(sph, chem, tend):
    sph.evolve_model(tend)
    update_chem(sph.particles, chem.particles)
    chem.evolve_model(tend)

def run_mc(N=5000, Mcloud=10000. | units.MSun, Rcloud=1. | units.parsec):
    timestep = 0.005 | units.Myr
    end_time = 0.12 | units.Myr

    conv = nbody_system.nbody_to_si(Mcloud, Rcloud)
    parts = molecular_cloud(targetN=N, convert_nbody=conv,
                            ethep_ratio=0.05, ekep_ratio=0.5).result

    rho_cloud = Mcloud / (4 / 3. * numpy.pi * Rcloud**3)

    tff = 1 / (4 * numpy.pi * constants.G * rho_cloud)**0.5
    parts.density = rho_cloud
    parts.radfield = 1.0|units.habing
    update_chem(parts, parts)

    print("Tcloud:", parts.temperature.max().in_(units.K))
    print("cloud n_H:", parts.number_density.max().in_(units.cm**-3))
    print("freefall time:", tff.in_(units.Myr))

    sph = Fi(conv)
    #chem = Krome(redirection="none")
    chem = UCLchem()
    print(parts)
    sph.parameters.self_gravity_flag = True
    sph.parameters.use_hydro_flag = True
    sph.parameters.isothermal_flag = True
    sph.parameters.integrate_entropy_flag = False
    sph.parameters.gamma = 1
    sph.parameters.verbosity = 0

    sph.parameters.timestep = timestep / 2

    sph.gas_particles.add_particles(parts)
    chem.particles.add_particles(parts)
    chem.out_species = ["OH", "OCS", "CO", "CS", "CH3OH"]
    H2_index = chem.get_index_of_species('H2')
    CO_index = chem.get_index_of_species('CO')
    tnow = sph.model_time

    sph_channel = sph.particles.new_channel_to(parts)
    chem_channel = chem.particles.new_channel_to(parts)


    f = pyplot.figure()
    pyplot.ion()
    pyplot.show()

    i = 0
    while i < (end_time / timestep + 0.5):
        evolve_sph_with_chemistry(sph, chem, i * timestep)
        tnow = sph.model_time
        sph_channel.copy()
        chem_channel.copy()
        print("done with step:", i, tnow.in_(units.Myr))
        i += 1

        n = (sph.particles.density / meanmwt).value_in(units.cm**-3)
        fh2 = chem.particles.abundances[:, H2_index]
        co = chem.particles.abundances[:, CO_index]
        
        pyplot.clf()
        pyplot.loglog(n, fh2, 'r.')
        pyplot.loglog(n, co, 'g.')
        pyplot.xlim(1.e3, 1.e6)
        pyplot.ylim(1.e-6, 1)
        pyplot.xlabel("density (cm**-3)")
        pyplot.ylabel("H_2,CO abundance")
        f.canvas.flush_events()
        particles_to_save = parts.copy()
        write_set_to_file(particles_to_save, "mol_cloud_uclchem_{}.txt".format(i), format='amuse',overwrite_file=True)
    print("done. press key to exit")
    input()


if __name__ == "__main__":
    run_mc(N=1000, Mcloud=2000. | units.MSun, Rcloud=.5 | units.parsec)
