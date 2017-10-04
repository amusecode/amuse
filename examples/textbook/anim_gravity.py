from matplotlib import pyplot
import matplotlib.animation as animation
from amuse.lab import *

fig = pyplot.figure()
ax = fig.add_subplot(1,1,1)

def animate(i):
    filename = "gravity.h5"
    particles = read_set_from_file(filename, "hdf5")
    x = []
    y = []
    for si in particles.history:
        x.append(si.x.value_in(units.AU))
        y.append(si.y.value_in(units.AU))
    ax.clear()
    ax.plot(x,y)

anim = animation.FuncAnimation(fig, animate, interval=100)
pyplot.show()

