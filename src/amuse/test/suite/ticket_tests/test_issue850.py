from amuse.test import amusetest

from amuse.datamodel import Particles
from amuse.units import units
from amuse.community.seba import Seba
from amuse.community.bse import Bse
from amuse.support.console import set_printing_strategy
from amuse.units import units


class TestsForIssue850(amusetest.TestCase):
    def create_stars_and_binaries(self, binary_pair=[0, 1]):
        stars = Particles(3)
        stars[0].initial_mass = 14 | units.MSun
        stars[1].initial_mass = 10 | units.MSun
        stars[2].initial_mass = 9 | units.MSun
        stars.mass = stars.initial_mass

        binary = Particles(1)
        binary.semi_major_axis = 30000 | units.RSun
        binary.eccentricity = 0.3
        binary.child1 = stars[binary_pair[0]]
        binary.child2 = stars[binary_pair[1]]

        return stars, binary

    def test_do_all_stars_evolve_no_binary(self, code=Seba):
        "Test if all stars evolve (no binaries)"
        stars, binary = self.create_stars_and_binaries(binary_pair=[0, 1])

        instance = code()
        instance.particles.add_particles(stars)

        channel_stars = instance.particles.new_channel_to(stars)
        channel_stars.copy()

        end_time = 13500. | units.Myr

        instance.evolve_model(end_time)
        channel_stars.copy()

        # All stars need to have evolved
        print(
            f"t: {end_time}, m: {stars[0].mass} {stars[1].mass} {stars[2].mass} "
        )
        for i in range(len(stars)):
            self.assertNotEqual(stars[i].mass, stars[i].initial_mass)
        instance.stop()

    def _test_do_all_stars_evolve_no_binary_bse(self):
        self.test_do_all_stars_evolve_no_binary(code=Bse)

    def test_does_seba_evolve_stars_not_in_binary(self, code=Seba):
        """
        Tests if code evolves any star that is not in a binary if a binary is
        added.
        """
        set_printing_strategy('default')
        stars, binary = self.create_stars_and_binaries(binary_pair=[0, 1])
        additional_stars = Particles(5)
        additional_stars.original_mass = [11, 20, 4, 2, 1] | units.MSun
        additional_stars.mass = additional_stars.original_mass
        stars.add_particles(additional_stars)

        instance = code()
        instance.particles.add_particles(stars)
        instance.binaries.add_particles(binary)

        channel_stars = instance.particles.new_channel_to(stars)
        channel_binary = instance.binaries.new_channel_to(binary)
        channel_stars.copy()
        channel_binary.copy()

        end_time = 13500. | units.Myr

        instance.evolve_model(end_time)
        channel_stars.copy()
        channel_binary.copy()

        # All stars need to have evolved
        print(
            f"t: {end_time}, m: {stars[0].mass} {stars[1].mass} {stars[2].mass} "
        )
        for i in range(len(stars)):
            self.assertNotEqual(stars[i].mass, stars[i].initial_mass)
        instance.stop()

    def _test_does_bse_evolve_stars_not_in_binary(self):
        self.test_does_seba_evolve_stars_not_in_binary(code=Bse)

    def test_do_all_stars_evolve_default_binary(self):
        "Test if all stars evolve (default order binary)"
        stars, binary = self.create_stars_and_binaries(binary_pair=[0, 1])

        instance = Seba()
        instance.particles.add_particles(stars)
        instance.binaries.add_particles(binary)

        channel_stars = instance.particles.new_channel_to(stars)
        channel_binary = instance.binaries.new_channel_to(binary)
        channel_stars.copy()
        channel_binary.copy()

        end_time = 13500. | units.Myr

        instance.evolve_model(end_time)
        channel_stars.copy()
        channel_binary.copy()

        # All stars need to have evolved
        print(
            f"t: {end_time}, m: {stars[0].mass} {stars[1].mass} {stars[2].mass} "
        )
        for i in range(len(stars)):
            self.assertNotEqual(stars[i].mass, stars[i].initial_mass)
        instance.stop()

    def test_do_all_stars_evolve_alternative_binary(self):
        "Test if all stars evolve (other stars in the binary)"
        stars, binary = self.create_stars_and_binaries(binary_pair=[0, 2])

        instance = Seba()
        instance.particles.add_particles(stars)
        instance.binaries.add_particles(binary)

        channel_stars = instance.particles.new_channel_to(stars)
        channel_binary = instance.binaries.new_channel_to(binary)
        channel_stars.copy()
        channel_binary.copy()

        end_time = 13500. | units.Myr

        instance.evolve_model(end_time)
        channel_stars.copy()
        channel_binary.copy()

        # All stars need to have evolved
        print(
            f"t: {end_time}, m: {stars[0].mass} {stars[1].mass} {stars[2].mass} "
        )
        for i in range(len(stars)):
            self.assertNotEqual(stars[i].mass, stars[i].initial_mass)
        instance.stop()

    def test_do_all_stars_evolve_alternative_binary_extra_stars(self):
        "Test if all stars evolve (other stars in the binary, extra stars)"
        stars, binary = self.create_stars_and_binaries(binary_pair=[0, 2])
        additional_stars = Particles(5)
        additional_stars.original_mass = [11, 20, 4, 2, 1] | units.MSun
        additional_stars.mass = additional_stars.original_mass
        stars.add_particles(additional_stars)

        instance = Seba()
        instance.particles.add_particles(stars)
        instance.binaries.add_particles(binary)

        channel_stars = instance.particles.new_channel_to(stars)
        channel_binary = instance.binaries.new_channel_to(binary)
        channel_stars.copy()
        channel_binary.copy()

        end_time = 13500. | units.Myr

        instance.evolve_model(end_time)
        channel_stars.copy()
        channel_binary.copy()

        # All stars need to have evolved
        print(
            f"t: {end_time}, m: {stars[0].mass} {stars[1].mass} {stars[2].mass} "
        )
        for i in range(len(stars)):
            self.assertNotEqual(stars[i].mass, stars[i].initial_mass)
        instance.stop()


if __name__ == "__main__":
    x = TestsForIssue850()
    x.test_do_all_stars_evolve_no_binary()
    x.test_does_seba_evolve_stars_not_in_binary()
    x.test_do_all_stars_evolve_default_binary()
    x.test_do_all_stars_evolve_alternative_binary_extra_stars()
    x.test_do_all_stars_evolve_alternative_binary()
