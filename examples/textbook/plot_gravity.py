import os
###BOOKLISTSTART###
from matplotlib import pyplot
from amuse.plot import scatter, xlabel, ylabel
from amuse.io import read_set_from_file
from amuse.units import units

def plot(x, y):
    
    pyplot.figure(figsize=(8,8))

    colormap = ['yellow', 'green', 'blue']	# specific to a 3-body plot
    size = [40, 20, 20]
    edgecolor = ['orange', 'green', 'blue']
    
    for si in particles.history:
        scatter(si.x, si.y, c=colormap, s=size, edgecolor=edgecolor)
    xlabel("x")
    ylabel("y")

    save_file = 'plot_gravity.png'
    pyplot.savefig(save_file)
    print '\nSaved figure in file', save_file,'\n'
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
    
    plot(x, y)
