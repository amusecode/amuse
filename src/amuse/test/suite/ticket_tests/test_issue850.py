from amuse.test import amusetest

from amuse.datamodel import Particles
from amuse.units import units
from amuse.community.seba.interface import SeBa
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

    def test_do_all_stars_evolve_no_binary(self):
        "Test if all stars evolve (no binaries)"
        stars, binary = self.create_stars_and_binaries(binary_pair=[0, 1])

        code = SeBa()
        code.particles.add_particles(stars)

        channel_stars = code.particles.new_channel_to(stars)
        channel_stars.copy()

        end_time = 15.6001481751 | units.Myr

        code.evolve_model(end_time)
        channel_stars.copy()
        print(
            f"t: {end_time}, m: {stars[0].mass} {stars[1].mass} {stars[2].mass} "
        )
        # All stars need to have evolved
        self.assertNotEqual(stars[0].mass, stars[0].initial_mass)
        self.assertNotEqual(stars[1].mass, stars[1].initial_mass)
        self.assertNotEqual(stars[2].mass, stars[2].initial_mass)
        # self.assertEqual(a, b)

    def test_do_all_stars_evolve_default_binary(self):
        "Test if all stars evolve (default order binary)"
        stars, binary = self.create_stars_and_binaries(binary_pair=[0, 1])

        code = SeBa()
        code.particles.add_particles(stars)
        code.binaries.add_particles(binary)

        channel_stars = code.particles.new_channel_to(stars)
        channel_binary = code.binaries.new_channel_to(binary)
        channel_stars.copy()
        channel_binary.copy()

        end_time = 15.6001481751 | units.Myr

        code.evolve_model(end_time)
        channel_stars.copy()
        channel_binary.copy()
        print(
            f"t: {end_time}, m: {stars[0].mass} {stars[1].mass} {stars[2].mass} "
        )
        # All stars need to have evolved
        self.assertNotEqual(stars[0].mass, stars[0].initial_mass)
        self.assertNotEqual(stars[1].mass, stars[1].initial_mass)
        self.assertNotEqual(stars[2].mass, stars[2].initial_mass)
        # self.assertEqual(a, b)

    def test_do_all_stars_evolve_alternative_binary(self):
        "Test if all stars evolve (other stars in the binary)"
        stars, binary = self.create_stars_and_binaries(binary_pair=[0, 2])

        code = SeBa()
        code.particles.add_particles(stars)
        code.binaries.add_particles(binary)

        channel_stars = code.particles.new_channel_to(stars)
        channel_binary = code.binaries.new_channel_to(binary)
        channel_stars.copy()
        channel_binary.copy()

        end_time = 15.6001481751 | units.Myr

        code.evolve_model(end_time)
        channel_stars.copy()
        channel_binary.copy()
        print(
            f"t: {end_time}, m: {stars[0].mass} {stars[1].mass} {stars[2].mass} "
        )
        # All stars need to have evolved
        self.assertNotEqual(stars[0].mass, stars[0].initial_mass)
        self.assertNotEqual(stars[1].mass, stars[1].initial_mass)
        self.assertNotEqual(stars[2].mass, stars[2].initial_mass)
        # self.assertEqual(a, b)


if __name__ == "__main__":
    x = TestsForIssue850()
    x.test_evolve_binary()
