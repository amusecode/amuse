"""
In this script we simulate the riemann shock tube problem in 3d.

.. quote:
    The test set-up consists of two fluids of different densities 
    and pressures separated by a membrane that is then removed. The 
    resulting solution shows all three types of fluid 
    discontinuties; a shock wave moving from the high density fluid 
    to the low density one, a rarefraction (sound) wave moving in 
    the opposite direction and a contact discontinuity which marks 
    the current location of the interface between the two fluids.

We follow the directions in the paper:

:title:     A test suite for quantitative comparison of hydrodynamic 
            codes in astrophyics
:authors:   Elizabeth J. Tasker, Riccardo Brunino, Nigel L. Mitchell, 
            Dolf Michielsen, Stephen Hopton, Frazer R. Pearce, Greg L. Bryan,
            Tom Theuns
:journal:   Monthly Notices of the Royal Astronomical Society
:issue:     Volume 390, Issue 3, pages 1267-1281, November 2008
:doi:       DOI: 10.1111/j.1365-2966.2008.13836.x

See also:

:site:      http://www.astro.ufl.edu/~tasker/codecomparison2/codecomparison_riemann.html

Exact solution is based on fortran code from Frank Timmer (who based it
on code from Bruce Fryxell), see:

http://cococubed.asu.edu/code_pages/exact_riemann.shtml
"""

from amuse.support.core import late

from amuse import io
from amuse.io import text
from amuse.units import constants
from amuse.units.generic_unit_converter import ConvertBetweenGenericAndSiUnits
from amuse.units.nbody_system import nbody_to_si
from amuse.units.quantities import VectorQuantity
from amuse.units.generic_unit_system import *

from amuse.community.athena.interface import Athena
from amuse.community.capreole.interface import Capreole
from amuse.community.mpiamrvac.interface import MpiAmrVac

from amuse.community.gadget2.interface import Gadget2
from amuse.community.fi.interface import Fi
from amuse.ext.grid_to_sph import convert_grid_to_SPH

import numpy

try:
    from amuse import plot
    from matplotlib import pyplot
    IS_PLOT_AVAILABLE = True
except ImportError:
    IS_PLOT_AVAILABLE = False
    
    
from numpy import sqrt, arange, searchsorted
from optparse import OptionParser

from amuse.datamodel import Grid
from amuse.datamodel import Particles
from amuse.datamodel.grids import SamplePointsOnMultipleGrids
from amuse.datamodel.grids import SamplePointWithIntepolation
from amuse.datamodel.grids import SamplePointOnCellCenter
class CalculateExactSolutionIn1D(object):
    number_of_points = 1000
    
    rho1 = 4.0 | density
    p1 = 1.0 | mass / (length * (time**2))
    u1 = 0.0 | speed
    
    rho5 = 1.0 | density
    p5 = 0.1795 | mass / (length * (time**2))
    u5 = 0.0 | speed
    
    
    gamma = 5.0/3.0
    
    def get_post_shock_pressure_p4(self, maximum_number_of_iterations = 20, maxium_allowable_relative_error = 1e-5):
        """solve for post-shock pressure by secant method"""
        p40 = self.p1
        p41 = self.p5
        p4  = p41
        
        f0  = self.calculate_p4_from_previous_value(p40)
        for x in range(maximum_number_of_iterations):
            f1 = self.calculate_p4_from_previous_value(p41)
            if (f1 == f0):
                return p4

            p4 = p41 - (p41 - p40) * f1 / (f1 - f0)

            error = abs (p4 - p41) / p41
            if (error < maxium_allowable_relative_error):
                return p4

            p40 = p41
            p41 = p4
            f0  = f1
        
        raise Exception("solution did not converge in less than {0!r} steps.".format(maximum_number_of_iterations))
        
    def get_post_shock_density_and_velocity_and_shock_speed(self, p4):
        z  = (p4 / self.p5 - 1.0)
        c5 = (self.gamma * self.p5 / self.rho5).sqrt()
        gm1 = self.gamma - 1.0
        gp1 = self.gamma + 1.0
        gmfac1 = 0.5 * gm1 / self.gamma
        gmfac2 = 0.5 * gp1 / self.gamma

        fact = sqrt (1. + gmfac2 * z)

        u4 = c5 * z / (self.gamma * fact)
        rho4 = self.rho5 * (1.0 + gmfac2 * z) / (1.0 + gmfac1 * z)
        return u4, rho4, c5 * fact
        
        
    def calculate_p4_from_previous_value(self, p4):
        c1 = sqrt(self.gamma * self.p1 / self.rho1)
        c5 = sqrt(self.gamma * self.p5 / self.rho5)

        gm1 = self.gamma - 1.0
        gp1 = self.gamma + 1.0
        g2  = 2.0 * self.gamma
        
        z= (p4 / self.p5 - 1.0)
        fact = gm1 / g2 * (c5 / c1) * z / sqrt (1. + gp1 / g2 * z)
        fact = (1. - fact) ** (g2 / gm1)

        return self.p1 * fact - p4
        
    
    def get_solution_at_time(self, t):
        p4 = self.get_post_shock_pressure_p4()
        u4, rho4, w = self.get_post_shock_density_and_velocity_and_shock_speed(p4)
        
        #compute values at foot of rarefaction
        p3 = p4
        u3 = u4
        rho3 = self.rho1 * (p3 / self.p1)**(1. /self.gamma)
        
        c1 = sqrt (self.gamma * self.p1 / self.rho1)
        c3 = sqrt (self.gamma * p3 / rho3)
        
        xi = 0.5 | length
        xr = 1.0 | length
        xl = 0.0 | length
        
        xsh = xi + w * t
        xcd = xi + u3 * t
        xft = xi + (u3 - c3) * t
        xhd = xi - c1 * t
        
        gm1 = self.gamma - 1.0
        gp1 = self.gamma + 1.0
        
        dx = (xr - xl) / (self.number_of_points - 1)
        x = xl + dx * arange(self.number_of_points)
        
        rho = VectorQuantity.zeros(self.number_of_points, density)
        p = VectorQuantity.zeros(self.number_of_points, mass / (length * time**2))
        u = VectorQuantity.zeros(self.number_of_points, speed)
        
        for i in range(self.number_of_points):
            if x[i] < xhd:
                rho[i] = self.rho1
                p[i]   = self.p1
                u[i]   = self.u1
            elif x[i] < xft:
                u[i]   = 2. / (self.gamma + 1.0) * (c1 + (x[i] - xi) / t)
                fact   = 1. - 0.5 * gm1 * u[i] / c1
                rho[i] = self.rho1 * fact ** (2. / gm1)
                p[i]   = self.p1 * fact ** (2. * self.gamma / gm1)
            elif x[i] < xcd:
                rho[i] = rho3
                p[i]   = p3
                u[i]   = u3
            elif x[i] < xsh:
                rho[i] = rho4
                p[i]   = p4
                u[i]   = u4
            else:
                rho[i] = self.rho5
                p[i]   = self.p5
                u[i]   = self.u5
                
        return x, rho,p,u


class CalculateSolutionIn3D(object):
    number_of_workers = 1
    number_of_grid_points = 10
    gamma = 5.0/3.0
    name_of_the_code = "capreole"
    hydro_code_options = dict()#redirection="none")
    convert_generic_units = ConvertBetweenGenericAndSiUnits(
        1.0 | units.kpc, 1.0e10 | units.MSun, constants.G)
    convert_nbody_units   = nbody_to_si(
        1.0 | units.kpc, 1.0e10 | units.MSun)
    
    def __init__(self, **keyword_arguments):
        for x in keyword_arguments:
            print x, keyword_arguments[x]
            setattr(self, x, keyword_arguments[x])
            
        self.dimensions_of_mesh = (
            self.number_of_grid_points * 3, 
            self.number_of_grid_points, 
            self.number_of_grid_points
        )
        
    def new_instance_of_code(self):
        attribute = "new_instance_of_{0}_code".format(self.name_of_the_code.lower())
        return getattr(self,attribute)()
        
    def new_instance_of_athena_code(self):
        result=Athena(number_of_workers=self.number_of_workers, **self.hydro_code_options)
        result.initialize_code()
        result.parameters.gamma = self.gamma
        result.parameters.courant_number=0.3
        return result
        

    def new_instance_of_mpiamrvac_code(self):
        result=MpiAmrVac(number_of_workers=self.number_of_workers, **self.hydro_code_options)
        result.set_parameters_filename(result.default_parameters_filename)
        result.initialize_code()
        return result
        
    def new_instance_of_capreole_code(self):
        result=Capreole(number_of_workers=self.number_of_workers, **self.hydro_code_options)
        result.initialize_code()
        return result
    
    grid_hydro_codes = ("athena", "mpiamrvac", "capreole")
    
    def new_instance_of_gadget2_code(self):
        result=Gadget2(self.convert_generic_units, number_of_workers=self.number_of_workers, 
            mode = "periodic", **self.hydro_code_options)
        result.initialize_code()
        return result
    
    def new_instance_of_fi_code(self):
        result=Fi(self.convert_nbody_units, number_of_workers=self.number_of_workers, 
            mode = "periodic", **self.hydro_code_options)
        result.initialize_code()
        result.parameters.self_gravity_flag = False
        result.parameters.timestep = 0.01 | time
        return result
    
    sph_hydro_codes = ("gadget2", "fi")
    
    def set_parameters(self, instance):
        if self.name_of_the_code in self.sph_hydro_codes:
            instance.parameters.periodic_box_size = 1.0 | length
        else:
            instance.parameters.mesh_size = self.dimensions_of_mesh
            
            instance.parameters.length_x = 1 | length
            instance.parameters.length_y = 1 | length
            instance.parameters.length_z = 1 | length
            
            instance.parameters.x_boundary_conditions = ("periodic","periodic")
            instance.parameters.y_boundary_conditions = ("periodic","periodic")
            instance.parameters.z_boundary_conditions = ("periodic","periodic")
        
        result = instance.commit_parameters()
    
    def new_grid(self):
        grid = Grid.create(self.dimensions_of_mesh, [1,1,1] | length)
        self.clear_grid(grid)
        return grid
        
    
    def clear_grid(self, grid):
        density = mass / length**3
        momentum =  speed * density
        energy =  mass / (time**2 * length)

        grid.rho =  0.0 | density
        grid.rhovx = 0.0 | momentum
        grid.rhovy = 0.0 | momentum
        grid.rhovz = 0.0 | momentum
        grid.energy = 0.0 | energy
    
        return grid
    
    def initialize_grid_with_shock(self, grid):
        energy =  mass / (time**2 * length)
        
        halfway = self.dimensions_of_mesh[0]/2 - 1
        
        firsthalf = grid.x <= 0.5 | length
        secondhalf = grid.x > 0.5 | length
        grid[firsthalf].rho = 4.0  | density
        grid[firsthalf].energy = (1.0 | energy)/ (self.gamma - 1)
        grid[secondhalf].rho = 1.0  | density
        grid[secondhalf].energy = (0.1795 | energy)/ (self.gamma - 1)
    
    def set_initial_conditions(self, instance):
        if self.name_of_the_code in self.sph_hydro_codes:
            initial_grid = self.new_grid()
            self.initialize_grid_with_shock(initial_grid)
            if self.name_of_the_code == "fi":
                initial_grid.position -= 0.5 | length
            sph_particles = convert_grid_to_SPH(initial_grid, self.number_of_particles)
            instance.gas_particles.add_particles(sph_particles)
        else:
            for x in instance.itergrids():
                inmem = x.copy()
                self.clear_grid(inmem)
                self.initialize_grid_with_shock(inmem)
                from_model_to_code = inmem.new_channel_to(x)
                from_model_to_code.copy()
            
            instance.initialize_grid()
    
    def copy_results(self, instance):
        if self.name_of_the_code in self.sph_hydro_codes:
            n_samples = self.number_of_grid_points * 3
            result = Particles(n_samples)
            result.position = [(x, 0.5, 0.5) | length for x in numpy.linspace(0.0, 1.0, n_samples)]
            x, y, z = result.x, result.y, result.z
            if self.name_of_the_code == "fi":
                x -= (0.5 | length)
                y -= (0.5 | length)
                z -= (0.5 | length)
            no_speed = [0.0]*n_samples | speed
            result.rho, result.rhovx, result.rhovy, result.rhovz, result.energy = [
                self.convert_generic_units.to_generic(quantity) for quantity in 
                instance.get_hydro_state_at_point(x, y, z, no_speed, no_speed, no_speed)]
        else:
            result = []
            for x in instance.itergrids():
                result.append(x.copy())
        return result
    
    def get_solution_at_time(self, time):
        instance=self.new_instance_of_code()
        self.set_parameters(instance)
        self.set_initial_conditions(instance)
        
        print "start evolve"
        instance.evolve_model(time)
        
        print "copying results"
        result = self.copy_results(instance)
        
        print "terminating code"
        instance.stop()
        
        return result
        
    
def store_attributes(x, rho, rhovx, energy, filename):
    output = text.CsvFileText(filename = filename)
    output.quantities = (x, rho, rhovx, energy)
    output.attribute_names = ("x", "rho", "rhovx", "energy")
    output.store()


def store_attributes_of_line(grid, yindex = 0, zindex = 0, **options):
    store_attributes(
        grid.x[...,yindex,zindex],
        grid.rho[...,yindex,zindex],
        grid.rhovx[...,yindex,zindex],
        grid.energy[...,yindex,zindex],
        filename = "riemann_shock_tube_{name_of_the_code}_{number_of_grid_points}_{number_of_workers}.csv".format(**options)
    )
def new_option_parser():
    result = OptionParser()
    result.add_option(
        "-n",
        "--mesh-size", 
        dest="number_of_grid_points",
        type="int",
        default = 10,
        help="number of grid cells in the x, y and z direction"
    )
    result.add_option(
        "-w",
        "--workers", 
        dest="number_of_workers",
        type="int",
        default = 1,
        help="number of parallel workers to start"
    )
    result.add_option(
        "-c",
        "--code",
        dest="name_of_the_code",
        default="athena",
        help="name of the code to use"
    )
    result.add_option(
        "-p",
        "--number_of_particles", 
        dest="number_of_particles",
        type="int",
        default = 10000,
        help="number of particles, in case code is an SPH code"
    )
    return result

def test_riemann_shocktube_problem():
    exact = CalculateExactSolutionIn1D()
    x, rho, p, u = exact.get_solution_at_time(0.12 | time)
    
    model = CalculateSolutionIn3D()
    model.name_of_the_code = "capreole"
    model.dimensions_of_mesh = (500,1,1)
    
    grids = model.get_solution_at_time(0.12 | time)
    grid = grids[0]
    model_x = grid.x[...,0,0]
    density = grid.rho[...,0,0]
    
    index_in_model = searchsorted(model_x.value_in(length), 0.56)
    index_in_exact = searchsorted(x.value_in(length), 0.56)
    
    #store_attributes_of_line(grid, name_of_the_code = "athena-test", number_of_grid_points = 500, number_of_workers = 1)
    
    assert abs((rho[index_in_exact] - density[index_in_model])/ density[index_in_model]) < 1.e-3 |units.none
    
    
def main(**options):
    print "calculating shock using exact solution"
    exact = CalculateExactSolutionIn1D()
    xpositions, rho, p, u = exact.get_solution_at_time(0.12 | time)
    
    print "calculating shock using code"
    model = CalculateSolutionIn3D(**options)
    grids = model.get_solution_at_time(0.12 | time)
    
    print "sampling grid"
    if model.name_of_the_code in model.sph_hydro_codes:
        samples = grids
    else:
        samplepoints = [(x, 0.5, 0.5) | length for x in numpy.linspace(0.0, 1.0, 2000)]
        print len(grids)
        samples = SamplePointsOnMultipleGrids(grids, samplepoints, SamplePointOnCellCenter)
        print len(samples)
        samples.filterout_duplicate_indices()
        print len(samples)
    
    print "saving data"
    filename = "riemann_shock_tube_rho_"+model.name_of_the_code
    store_attributes(xpositions,rho,u,p,filename=filename+"_exact.csv")
    store_attributes(samples.x,samples.rho,samples.rhovx,samples.energy,filename=filename+".csv")
    
    if IS_PLOT_AVAILABLE:
        print "plotting solution"
        from prepare_figure import *
        from distinct_colours import get_distinct
        x_label = "[length]"
        y_label = "[mass/length$^3$]"
        figure = single_frame(x_label, y_label, logx=False, logy=True, xsize=14, ysize=10)
        color = get_distinct(2)
        
        plot.plot(xpositions,rho, c=color[0])
        plot.scatter(samples.x, samples.rho, c=color[1], s=100)
#        pyplot.xlim(0.3,0.7)
#        pyplot.ylim(0.5,4.5)
        pyplot.xlim(0.0,1.0)
#        pyplot.ylim(0.5,4.5)
        pyplot.xlabel("[length]")
        pyplot.ylabel("[mass/length$^3$]")
#        pyplot.savefig("riemann_shock_tube_rho_"+model.name_of_the_code+".png")
        pyplot.savefig(filename)
        pyplot.show()

if __name__ == "__main__":
    options, arguments  = new_option_parser().parse_args()
    main(**options.__dict__)
