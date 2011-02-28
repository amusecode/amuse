from amuse.community import *
from amuse.test.amusetest import TestWithMPI

from interface import supportInterface
import numpy as np


if __name__ == '__main__':

    instance = supportInterface(redirection='none')
    x = np.zeros(2000)
    y = np.zeros(2000)
    z = np.zeros(2000)
    #R = instance.many_points_on_sphere(x,y,z)
    R = instance.rnd_points_on_sphere(x,y,z)
    x1 = R.x
    y1 = R.y
    z1 = R.z

    for i, v in enumerate(x):
        print x1[i], y1[i], z1[i]

    instance.stop()

