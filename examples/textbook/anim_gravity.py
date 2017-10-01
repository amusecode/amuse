import os
###BOOKLISTSTART###
from matplotlib import pyplot
import matplotlib.animation as animation
from amuse.io import read_set_from_file
from amuse.units import units

def animate(x, y):

    def update(i):
        while i >= np: i -= np
        off = []
        for j in range(len(x[i])):
            off.append(x[i][j])
            off.append(y[i][j])
        scat.set_offsets(off)
        return scat,

    np = len(x)
    fig = pyplot.figure(figsize=(8,8))
    ax = fig.add_subplot(1,1,1)
    ax.set_xlim(-1.2, 1.2)
    ax.set_ylim(-1.2, 1.2)

    colormap = ['yellow', 'green', 'blue']
    size = [40, 20, 20]
    edgecolor = ['orange', 'green', 'blue']
    
    scat = ax.scatter(x[0], y[0], c=colormap, s=size, edgecolor=edgecolor)
    anim = animation.FuncAnimation(fig, update, interval=100)
    pyplot.show()
###BOOKLISTSTOP###

if __name__ in ('__main__'):
    
    try:
        amusedir = os.environ['AMUSE_DIR']
        dir = amusedir+'/examples/textbook/'
    except:
        print 'Environment variable AMUSE_DIR not set'
        dir = './'
    
    filename = dir+'gravity.h5'
    particles = read_set_from_file(filename, "hdf5")

    x = []
    y = []
    for si in particles.history:
        x.append(si.x.value_in(units.AU))
        y.append(si.y.value_in(units.AU))
    
    animate(x, y)
