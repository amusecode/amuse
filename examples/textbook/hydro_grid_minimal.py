import numpy
from amuse.lab import *
from amuse.units import units, constants
from amuse.ext.sph_to_grid import convert_SPH_to_grid
from amuse.units.generic_unit_converter import ConvertBetweenGenericAndSiUnits


def plot_grid(grid):
    from matplotlib import pyplot
    density = units.g / units.cm**3
    rho = grid.rho[:, :, 0].value_in(density)
    figure = pyplot.figure(figsize=(6, 6))
    plot = figure.add_subplot(1, 1, 1)
    plot.imshow(rho, origin='lower')
    figure.savefig('kelvin_helmholtz.png')
    pyplot.show()


def setup_sph_code(sph_code, N, L, rho, u):
    converter = ConvertBetweenGenericAndSiUnits(L, rho, constants.G)
    sph_code = sph_code(converter, mode='periodic')
    sph_code.parameters.periodic_box_size = 10.0 | units.parsec
    plummer = new_plummer_gas_model(N, convert_nbody=converter)    
    plummer = plummer.select(lambda r: r.length()<0.5*L,["position"])
    N = len(plummer)
    print("N=", len(plummer))
    plummer.mass = (rho * L**3) / N
    gas = Particles(N)
    gas.mass = 0.01*(rho * L**3) / N
    numpy.random.seed(12345)
    gas.x = L * numpy.random.uniform(0.0, 1.0, N)
    gas.y = L * numpy.random.uniform(0.0, 1.0, N)
    gas.z = L * numpy.random.uniform(0.0, 1.0, N)
    gas.vx = numpy.zeros(N) | units.cm / units.s
    gas.vy = numpy.zeros(N) | units.cm / units.s
    gas.vz = numpy.zeros(N) | units.cm / units.s
    gas.u = u
    if isinstance(sph_code, Fi):
        sph_code.parameters.self_gravity_flag = False
        sph_code.parameters.timestep = 0.1 | generic_unit_system.time
        gas.h_smooth = L / N**(1/3.0)
        gas.position -= 0.5 * L
        
    sph_code.gas_particles.add_particles(gas)
    sph_code.gas_particles.add_particles(plummer)
    sph_code.commit_particles()
    return sph_code


def main(N, Mtot, Rvir, t_end):

    rho = 1.14 | units.amu/units.cm**3
    u = 5.e11 | units.cm**2 / units.s**2
    sph_code = setup_sph_code(Fi, N, Rvir, rho, u)    

    grid = convert_SPH_to_grid(sph_code, (10, 10, 10), do_scale = True)
    print(grid)
    sph_code.stop()
    plot_grid(grid)
    exit()
    # the code below is not yet used, to be fixed?
    hydro = Athena(converter)
    hydro.parameters.gamma = 1.4
    hydro.parameters.courant_number=0.8
    hydro.parameters.mesh_size
    
    hydro.gas_particles.add_particles(gas)
    Etot_init = hydro.kinetic_energy \
              + hydro.potential_energy + hydro.thermal_energy
    hydro.evolve_model(t_end)
    write_set_to_file(hydro.particles, "hydro.h5", "hdf5")

    Ekin = hydro.kinetic_energy 
    Epot = hydro.potential_energy
    Eth = hydro.thermal_energy
    Etot = Ekin + Epot + Eth
    Q = (Ekin+Eth)/Epot
    dE = (Etot_init-Etot)/Etot
    com = hydro.gas_particles.center_of_mass()
    print("T=", hydro.get_time(), "M=", hydro.gas_particles.mass.sum(), end=' ') 
    print("E= ", Etot, "Q= ", Q, "dE=", dE, "CoM=", com.in_(units.RSun))

    hydro.stop()
    
def new_option_parser():
    from amuse.units.optparse import OptionParser
    result = OptionParser()
    result.add_option("-N", dest="N", type="int",default = 10000,
                      help="number of gas particles [%default]")
    result.add_option("-t", unit=units.Myr,
                      dest="t_end", type="float", default = 1|units.Myr,
                      help="end time of the simulation [%default]")
    result.add_option("-M", unit=units.MSun,
                      dest="Mtot", type="float", default = 1000|units.MSun,
                      help="Mass of the cloud [%default]")
    result.add_option("-R", unit=units.RSun,
                      dest="Rvir", type="float", default = 10|units.parsec,
                      help="Radius of the cloud [%default]")
    return result

if __name__ in ('__main__', '__plot__'):
    o, arguments  = new_option_parser().parse_args()
    main(**o.__dict__)

