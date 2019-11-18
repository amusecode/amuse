import numpy
from amuse.lab import *
#from amuse import plot
from amuse.ext import cloud
from matplotlib import pyplot
from amuse import datamodel
from amuse.ext.sph_to_grid import convert_SPH_to_grid
from amuse.community.capreole.interface import Capreole
from amuse.units.generic_unit_converter import ConvertBetweenGenericAndSiUnits

def plot_grid(grid, time= 0.0|units.day):

    pyplot.rcParams.update({'font.size': 30})
    figure = pyplot.figure(figsize=(12, 12))

    halfway = len(grid.rho[...,0,0])/2 
    rho = grid.rho[:,:,halfway].value_in(units.g/units.cm**3)
    print "Extrema density:", halfway, rho.min(), rho.max()
    max_dens = rho.max()
#    max_dens = 32

    plot = figure.add_subplot(1,1,1)
#    cax = plot.imshow(rho, interpolation='nearest', origin = 'lower', extent=[-5, 5, -5, 5], cmap="hot")
    cax = plot.imshow(rho, interpolation='bicubic', origin = 'lower', extent=[-5, 5, -5, 5], cmap="hot")
    cbar = figure.colorbar(cax, ticks=[1.e-8, 0.5*max_dens, max_dens], orientation='vertical', fraction=0.045)
    rmin = 0.0 
    rmid = "%.1f" % (0.5*max_dens)
    rmax = "%.1f" % (max_dens)
#    cbar.ax.set_yticklabels(['Low', ' ', 'High'])  # horizontal colorbar
    cbar.ax.set_yticklabels([rmin, ' ', rmax])  # horizontal colorbar
    cbar.set_label('mid-plane density [$g/cm^3$]', rotation=270)
    pyplot.xlabel("x [R$_\odot$]")
    pyplot.ylabel("y [R$_\odot$]")
    t = int(time.value_in(units.s))
    filename = "supernova_grid_T"+str(t)+".png"
    figure.savefig(filename)
#    pyplot.show()

def setup_sph_code(sph_code, N, L, rho, u):
    converter = ConvertBetweenGenericAndSiUnits(L, rho, constants.G)
    sph_code = sph_code(converter, mode = 'periodic')#, redirection = 'none')
    sph_code.parameters.periodic_box_size = 10.0 | units.parsec
    plummer = new_plummer_gas_model(N, convert_nbody=converter)    
    plummer = plummer.select(lambda r: r.length()<0.5*L,["position"])
    N = len(plummer)
    print "N=", len(plummer)
    plummer.mass = (rho * L**3) / N
    gas = Particles(N)
    gas.mass = 0.001*(rho * L**3) / N
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
    
def main(stellar_mass, stellar_radius, core_mass, core_radius, t_end, dt, resolution):
    grid_size = 10 * stellar_radius
    hydro = initialize_grid_code(resolution, grid_size)
    grid = initialize_grid(stellar_mass, stellar_radius, core_mass, core_radius, resolution, grid_size)
    cth = grid.new_channel_to(hydro.grid)
    cth.copy()
    hydro.initialize_grid()

    run_grid_code(hydro, grid, t_end, dt)
    
def initialize_grid_code(resolution, grid_size):

    converter = nbody_system.nbody_to_si(1|units.MSun, 1|units.RSun)
    instance = Athena(converter, number_of_workers=4)
    instance.initialize_code()

    instance.parameters.gamma = 5/3.0
    instance.parameters.courant_number=0.3

    instance.parameters.nx = resolution
    instance.parameters.ny = resolution
    instance.parameters.nz = resolution
        
    instance.parameters.length_x = grid_size
    instance.parameters.length_y = grid_size
    instance.parameters.length_z = grid_size
        
    instance.x_boundary_conditions = ("outflow", "outflow")
    instance.y_boundary_conditions = ("outflow", "outflow")
    instance.z_boundary_conditions = ("outflow", "outflow")
    instance.commit_parameters()

    return instance
    
def initialize_grid(stellar_mass, stellar_radius, core_mass, core_radius, resolution, grid_size):
    n = resolution
    r = grid_size.value_in(units.RSun)
    grid = datamodel.new_regular_grid((n,n,n), [r, r, r] | units.RSun)

    momentum =  units.kg / (units.s * units.m**2)
    grid_size = units.RSun
    energy_density = units.erg / grid_size**3

    supernova_energy = 1.e+51 | units.erg
    stellar_energy_density = 0.01*supernova_energy/stellar_radius**3
    supernova_energy_density = supernova_energy/core_radius**3
    stellar_density = stellar_mass/stellar_radius**3
        
    grid.rho = 1.e-10 * stellar_density
    grid.rhovx = 0.0 | momentum
    grid.rhovy = 0.0 | momentum
    grid.rhovz = 0.0 | momentum
    grid.energy = 1 | energy_density
        
    datamodel.Grid.add_global_vector_attribute("position", ["x","y","z"])

    cloud.fill_grid_with_spherical_cloud(
        grid, 
        center = [5.0, 5.0, 5.0] | units.RSun,
        radius = stellar_radius,
        rho = stellar_density,
        rhovx = 0.0 | momentum,
        rhovy = 0.0 | momentum,
        rhovz = 0.0 | momentum, 
        energy = stellar_energy_density
    )

    cloud.fill_grid_with_spherical_cloud(
        grid, 
        center = [5.0, 5.0, 5.0] | units.RSun,
        radius = core_radius,
        rho = core_mass/core_radius**3,
        rhovx = 0.0 | momentum,
        rhovy = 0.0 | momentum,
        rhovz = 0.0 | momentum, 
        energy = supernova_energy_density
        #subgridsize = 16,
    )

    return grid
         
def run_grid_code(hydro, grid, t_end, dt):
    ctg = hydro.grid.new_channel_to(grid)
    ctg.copy()
    plot_grid(grid)
        
    while hydro.model_time<t_end:
        print "Time=", hydro.model_time.in_(units.s)
        hydro.evolve_model(hydro.model_time + dt)
        ctg.copy()
        plot_grid(grid, hydro.model_time)
    hydro.stop()
        
def new_option_parser():
    from amuse.units.optparse import OptionParser
    result = OptionParser()
    result.add_option("-t", unit=units.s,
                      dest="t_end", type="float", default = 300.0|units.s,
                      help="end time of the simulation [%default]")
    result.add_option("-d", unit=units.s,
                      dest="dt", type="float", default = 50.0|units.s,
                      help="diagnostic time step [%default]")
    result.add_option("-M", unit=units.MSun,
                      dest="stellar_mass", type="float", default = 3|units.MSun,
                      help="Mass of the star [%default]")
    result.add_option("-R", unit=units.RSun,
                      dest="stellar_radius", type="float", default = 1|units.RSun,
                      help="Radius of the star [%default]")
    result.add_option("-m", unit=units.MSun,
                      dest="core_mass", type="float", default = 1.4|units.MSun,
                      help="Mass of the stellar core [%default]")
    result.add_option("-n", 
                      dest="resolution", type="int", default = 256,
                      help="Resolution of the grid [%default]")
    result.add_option("-r", unit=units.RSun,
                      dest="core_radius", type="float", default = 0.1|units.RSun,
                      help="Radius of the stellar core [%default]")
    return result

if __name__ in ('__main__', '__plot__'):
    o, arguments  = new_option_parser().parse_args()
    main(**o.__dict__)

