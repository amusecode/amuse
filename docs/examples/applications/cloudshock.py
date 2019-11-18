import numpy
try:
    from matplotlib import pyplot
    IS_PLOT_AVAILABLE = True
except ImportError:
    IS_PLOT_AVAILABLE = False


from amuse.ext import cloud
from amuse.units import generic_unit_system
from amuse.community.capreole.interface import Capreole
from amuse import io

from amuse import datamodel


class CalculateCloudShock(object):
    number_of_workers = 1
    number_of_grid_points = 10
    mesh_length = 10.0 | generic_unit_system.length
    gamma = 5.0/3.0
    name_of_the_code = "athena"

    def __init__(
            self,
            number_of_grid_points=10,
            number_of_workers=1,
            name_of_the_code="capreole"):

        self.number_of_grid_points = number_of_grid_points
        self.number_of_workers = number_of_workers
        self.name_of_the_code = name_of_the_code

        self.dimensions_of_mesh = (
            self.number_of_grid_points,
            self.number_of_grid_points * 4,
            self.number_of_grid_points,
        )

    def new_instance_of_code(self):
        attribute = "new_instance_of_{0}_code".format(
            self.name_of_the_code.lower())
        return getattr(self, attribute)()

    def new_instance_of_capreole_code(self):
        result = Capreole(number_of_workers=self.number_of_workers)
        result.initialize_code()

        result.parameters.x_boundary_conditions = ("reflective", "reflective")
        result.parameters.y_boundary_conditions = ("outflow", "outflow")
        result.parameters.z_boundary_conditions = ("reflective", "reflective")

        return result

    def new_instance_of_athena_code(self):
        from amuse.community.athena.interface import Athena
        result = Athena(
            number_of_workers=self.number_of_workers, redirection="none")
        result.initialize_code()
        result.parameters.gamma = self.gamma
        result.parameters.courant_number = 0.3

        result.parameters.x_boundary_conditions = ("reflective", "reflective")
        result.parameters.y_boundary_conditions = ("outflow", "outflow")
        result.parameters.z_boundary_conditions = ("reflective", "reflective")

        return result

    def new_instance_of_mpiamrvac_code(self):
        from amuse.community.mpiamrvac.interface import MpiAmrVac
        result = MpiAmrVac(number_of_workers=self.number_of_workers)  # , redirection="none")
        result.set_parameters_filename(result.default_parameters_filename)
        result.initialize_code()
        result.parameters.maximum_number_of_grid_levels = 3
        result.parameters.spatial_discretization_method = 'tvdmu'
        result.parameters.predictor_step_discretization_method = 'tvdmu'
        result.parameters.entropy_type = 'powell'
        result.parameters.courant_number = 0.8

        result.parameters.x_boundary_conditions = ("cont", "cont")
        result.parameters.y_boundary_conditions = ("cont", "cont")
        result.parameters.z_boundary_conditions = ("cont", "cont")

        return result

    def set_parameters(self, instance):
        instance.parameters.mesh_size = self.dimensions_of_mesh

        instance.parameters.length_x = self.mesh_length
        instance.parameters.length_y = 4 * self.mesh_length
        instance.parameters.length_z = self.mesh_length

        result = instance.commit_parameters()

    def new_grid(self):
        grid = Grid.create(self.dimensions_of_mesh, [
                           1, 1, 1] | generic_unit_system.length)
        self.clear_grid(grid)
        return grid

    def initialize_grid(self, grid):
        center = self.mesh_length/2.0 * [1.0, 1.0, 1.0]
        cloud.fill_grid_with_cloud_shock(
            grid,
            center=center,
            radius=1.0 | generic_unit_system.length,
            subgridsize=4
        )

    def store_grids(self, grids, step):
        if __name__ == '__plot__':
            return

        grids_in_memory = [x.copy() for x in grids]
        io.write_set_to_file(
            grids_in_memory,
            "cloudshock_{2}_{0}_{1}.vtu".format(
                self.number_of_grid_points, step, self.name_of_the_code),
            "vtu",
            is_multiple=True
        )

    def refine_grid(self, instance):

        if hasattr(instance, 'refine_grid'):
            must_refine = True

            while must_refine:
                must_refine = instance.refine_grid()

                for x in instance.itergrids():
                    inmem = x.copy()
                    self.initialize_grid(inmem)
                    from_model_to_code = inmem.new_channel_to(x)
                    from_model_to_code.copy()

    def get_tau(self):
        rc = 1.
        xi = 10.
        cs = numpy.sqrt((self.gamma-1.0))
        cs_out = numpy.sqrt((self.gamma-1.0)*xi)
        vs = cs_out*2.7
        return (1.6*2*rc*xi**0.5/vs) | generic_unit_system.time

    def get_solution_at_time(self, time):
        instance = self.new_instance_of_code()

        self.set_parameters(instance)

        for x in instance.itergrids():
            inmem = x.copy()
            self.initialize_grid(inmem)
            from_model_to_code = inmem.new_channel_to(x)
            from_model_to_code.copy()

        self.refine_grid(instance)
        instance.initialize_grid()

        self.store_grids(instance.itergrids(), 0)

        if time > 0.0 | generic_unit_system.time:
            print "start evolve"
            dt = time / 10.0
            t = dt
            step = 1
            while t <= time:
                instance.evolve_model(t)

                print "time : ", t

                # self.store_grids(instance.itergrids(), step)

                t += dt
                step += 1

        print "sampling results"
        sample = datamodel.Grid.create(
            (1000, 4000),
            (10.0, 40) | generic_unit_system.length
        )
        sample.z = 5.0 | generic_unit_system.length

        rho, rhovx, rhovy, rhovx, rhoen = instance.get_hydro_state_at_point(
            sample.x.flatten(), sample.y.flatten(), sample.z.flatten())
        sample.rho = rho.reshape(sample.shape)
        sample.rhovx = rhovx.reshape(sample.shape)
        sample.rhovy = rhovy.reshape(sample.shape)
        sample.rhovz = rhovx.reshape(sample.shape)
        sample.energy = rhoen.reshape(sample.shape)
        print "terminating code"
        instance.stop()

        return sample


def main():
    number_of_grid_points = 40
    name_of_the_code = 'athena'
    model = CalculateCloudShock(
        number_of_grid_points=number_of_grid_points,
        number_of_workers=6,
        name_of_the_code=name_of_the_code
    )

    result = model.get_solution_at_time(0.75 * model.get_tau())

    rho = result.rho.value_in(generic_unit_system.density)

    print "done"

    if not IS_PLOT_AVAILABLE:
        return

    levels = numpy.linspace(numpy.min(rho), numpy.max(rho), 255)
    figure = pyplot.figure(figsize=(10, 10))
    plot = figure.add_subplot(1, 1, 1)
    plot.imshow(rho, origin='lower')
    figure.savefig('cloudshock_{0}_{1}.png'.format(
        name_of_the_code, number_of_grid_points))
    pyplot.show()


if __name__ == "__main__":
    main()
