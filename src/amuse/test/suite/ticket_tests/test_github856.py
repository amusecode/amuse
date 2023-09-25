import os
from amuse.test import amusetest

from amuse.datamodel import new_cartesian_grid, Particles

from amuse.io import read_set_from_file, write_set_to_file


class test_github856(amusetest.TestCase):
    def test1(self):
        filename = os.path.join(self.get_path_to_results(), "github856.amuse")

        g1 = new_cartesian_grid((5, 5), 1)
        write_set_to_file(g1, filename, "amuse")
        del g1
        g2 = read_set_from_file(filename, "amuse")
        self.assertEquals(g2.get_axes_names(), "xy")

    def test2(self):
        g1 = Particles(lon=[1, 2], lat=[3, 4])
        g1.add_vector_attribute("lonlat", ["lon", "lat"])
        g2 = g1.copy()
        self.assertEquals(g2.lonlat, [[1, 3], [2, 4]])

    def test3(self):
        filename = os.path.join(self.get_path_to_results(), "github856_2.amuse")

        g1 = Particles(lon=[1, 2], lat=[3, 4])
        g1.add_vector_attribute("lonlat", ["lon", "lat"])
        write_set_to_file(g1, filename, "amuse")
        del g1
        g2 = read_set_from_file(filename, "amuse")
        self.assertEquals(g2.lonlat, [[1, 3], [2, 4]])
