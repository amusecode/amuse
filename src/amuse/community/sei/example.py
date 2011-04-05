import os
import sys
import numpy as np
import math

from interface import SeiInterface

OMEGA = 1.0

if __name__ == '__main__':
    instance = SeiInterface()
    instance.initialization()
    impactparameter = 8.0 * 0.69336127435063 #8 Hill radii
    instance.set_state(impactparameter, 
                       np.pi * 300.0/2.0*OMEGA*impactparameter,
                       0.0,
                       0.0,
                       -3.0/2.0*OMEGA*impactparameter,
                       0.0)

    for t in np.arange(0, 6,1):
        print instance.get_state()['x'],\
              instance.get_state()['y'],\
              instance.get_state()['z'] 
        instance.evolve(t)

    instance.stop()
    
