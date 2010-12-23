
from amuse.support.data import base
from amuse.support.data import grids


grids.Grid.add_global_vector_attribute("position", ["x","y","z"])
grids.Grid.add_global_vector_attribute("momentum", ["rhovx","rhovy","rhovz"])