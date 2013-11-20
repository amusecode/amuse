import numpy as np

from pylab import *

def snapshots(data,zcut):
    """ Plot 2d cut"""

    #Plot data
    figure()
    xlabel('x') #assumes dir = z
    ylabel('y')
    title('Plot of rho')
    imshow(np.transpose(data.cube[:,:,zcut]),\
               extent=[data.xmin,data.xmax,data.ymin,data.ymax])
    colorbar()
