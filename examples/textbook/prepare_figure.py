from matplotlib import pyplot
#import seaborn 

# These are the "Tableau 20" colors as RGB.    
tableau20 = [(31, 119, 180), (174, 199, 232), (255, 127, 14), (255, 187, 120),    
             (44, 160, 44), (152, 223, 138), (214, 39, 40), (255, 152, 150),    
             (148, 103, 189), (197, 176, 213), (140, 86, 75), (196, 156, 148),    
             (227, 119, 194), (247, 182, 210), (127, 127, 127), (199, 199, 199),    
             (188, 189, 34), (219, 219, 141), (23, 190, 207), (158, 218, 229)]   


# colors
almost_black = '#262626'
blue='#3D52A1'
light_blue='#88CCEE'
cyan='#44AA99'
green='#117733'
red='#AE1C3E'
sand='#999933'
yellow='#DDCC77'
pink='#CC6677'
crimson='#882255'
violet='#AA4499'
brown='#661100'
steal='#6699CC'
rose='#AA4466'
sky_blue='#4477AA'

# -*- coding: iso-8859-1 -*-

"""
Colour-blind proof distinct colours module, based on work by Paul Tol
Pieter van der Meer, 2011
SRON - Netherlands Institute for Space Research
"""

# colour table in HTML hex format
hexcols = ['#332288', '#88CCEE', '#44AA99', '#117733', '#999933', '#DDCC77', 
           '#CC6677', '#882255', '#AA4499', '#661100', '#6699CC', '#AA4466',
           '#4477AA']

greysafecols = ['#809BC8', '#FF6666', '#FFCC66', '#64C204']

xarr = [[12], 
        [12, 6], 
        [12, 6, 5], 
        [12, 6, 5, 3], 
        [0, 1, 3, 5, 6], 
        [0, 1, 3, 5, 6, 8], 
        [0, 1, 2, 3, 5, 6, 8], 
        [0, 1, 2, 3, 4, 5, 6, 8], 
        [0, 1, 2, 3, 4, 5, 6, 7, 8], 
        [0, 1, 2, 3, 4, 5, 9, 6, 7, 8], 
        [0, 10, 1, 2, 3, 4, 5, 9, 6, 7, 8], 
        [0, 10, 1, 2, 3, 4, 5, 9, 6, 11, 7, 8]]

# get specified nr of distinct colours in HTML hex format.
# in: nr - number of colours [1..12]
# returns: list of distinct colours in HTML hex
def get_distinct(nr):

    #
    # check if nr is in correct range
    #
    
    if nr < 1 or nr > 12:
        print "wrong nr of distinct colours!"
        return

    #
    # get list of indices
    #
    
    lst = xarr[nr-1]
    
    #
    # generate colour list by stepping through indices and looking them up
    # in the colour table
    #

    i_col = 0
    col = [0] * nr
    for idx in lst:
        col[i_col] = hexcols[idx]
        i_col+=1
    return col

# displays usage information and produces example plot.
if __name__ == '__main__':
    import numpy as np
    import matplotlib.mlab as mlab
    import matplotlib.pyplot as plt

    print __doc__
    print "usage examples: "
    print "print distinct_colours.get_distinct(2)"
    print get_distinct(2)
    print "print distinct_colours.greysafecols"
    print greysafecols

    print "\ngenerating example plot: distinct_colours_example.png"
    plt.close()
    t = np.arange(0.0, 2.0, 0.01)
    n = 12
    cols = get_distinct(n)
    d = 2./n
    for i in range(n):
        s = np.sin(i*d*np.pi*t)
        plt.plot(t, s, linewidth=2.0, c=cols[i])

    plt.xlabel('time (s)')
    plt.ylabel('voltage (mV)')
    plt.title('Distinct colours example')
    plt.grid(True)
    plt.savefig("distinct_colours_example.pdf")

def figure_frame(x_label, y_label, xsize=12, ysize=10):
    figure = pyplot.figure(figsize=(xsize, ysize))
    plot = figure.add_subplot(1,1,1)
    ax = pyplot.gca()
    ax.minorticks_on() # switch on the minor ticks
    ax.locator_params(nbins=3)
    pyplot.xlabel(x_label)
    pyplot.ylabel(y_label)
    return figure, ax

from cycler import cycler
def single_frame(x_label, y_label, logx=False, logy=False, xsize=12, ysize=10,
                 ymin=-1, ymax=-1):

    pyplot.rcParams.update({'font.size': 25})
    #pyplot.rcParams['axes.color_cycle'] = [blue, green, red, sand, light_blue,
    #                                       pink, crimson, violet, brown,
    #                                       steal, rose, yellow, cyan ]
    pyplot.rcParams['axes.prop_cycle'] \
        = (cycler('color', [blue, green, red, sand, light_blue,
                            pink, crimson, violet, brown,
                            steal, rose, yellow, cyan ]))
    figure = pyplot.figure(figsize=(xsize, ysize))

    ax = pyplot.gca()
    ax.get_yaxis().get_major_formatter().set_useOffset(False)
    
    pyplot.xlabel(x_label)
    pyplot.ylabel(y_label)

    if ymax>0:
        pyplot.ylim(ymin, ymax)

    set_tickmarks(ax)

    ax.xaxis._autolabelpos = True
    ax.yaxis._autolabelpos = True

    if logx is True:
        ax.set_xscale('log')
        ax.get_xaxis().set_tick_params(pad=7)
    if logy is True:
        ax.set_yscale('log')
        ax.get_yaxis().set_tick_params(pad=7)
    return figure

def quad_frame(x_label, y_label, logx=False, logy=False, xsize=12, ysize=12):

    f, ((ax1, ax2), (ax3, ax4)) = pyplot.subplots(2, 2, sharex='col',
                                                  sharey='row', figsize=(8,8))
    set_tickmarks(ax1)
    ax1.locator_params(nbins=3)
    set_tickmarks(ax2)
    ax2.locator_params(nbins=3)
    set_tickmarks(ax3)
    ax3.locator_params(nbins=3)
    set_tickmarks(ax4)
    ax4.locator_params(nbins=3)

    ax1.xaxis.set_label_position("top")
    ax2.xaxis.set_label_position("top")
    ax1.set_xlabel(x_label)
    ax2.set_xlabel(x_label)

    ax2.yaxis.set_label_position("right")
    ax4.yaxis.set_label_position("right")
    ax2.set_ylabel(y_label)
    ax4.set_ylabel(y_label)
    return f, ax1, ax2, ax3, ax4

def _quad_frame(x_label, y_label, logx=False, logy=False, xsize=12, ysize=12):

    pyplot.rcParams.update({'font.size': 25})
    #pyplot.rcParams['axes.color_cycle'] = [blue, green, red, sand, light_blue,
    #                                       pink, crimson, violet, brown,
    #                                       steal, rose, yellow, cyan ]
    pyplot.rcParams['axes.prop_cycle'] \
        = (cycler('color', [blue, green, red, sand, light_blue,
                            pink, crimson, violet, brown,
                            steal, rose, yellow, cyan ]))
    f, ((ax1, ax2), (ax3, ax4)) = pyplot.subplots(2, 2, sharex='col',
                                                  sharey='row', figsize=(12,12))
    set_tickmarks(ax1)
    set_tickmarks(ax2)
    set_tickmarks(ax3)
    set_tickmarks(ax4)

    ax1.xaxis.set_label_position("top")
    ax2.xaxis.set_label_position("top")
    ax1.set_xlabel(x_label)
    ax2.set_xlabel(x_label)

    ax2.yaxis.set_label_position("right")
    ax4.yaxis.set_label_position("right")
    ax2.set_ylabel(y_label)
    ax4.set_ylabel(y_label)

    return f, ax1, ax2, ax3, ax4

def set_tickmarks(ax):
    ax.minorticks_on()
    return
    ax.tick_params('both', length=15, width=2, which='major')
    ax.tick_params('both', length=6, width=1, which='minor')
    ax.locator_params(nbins=3)
    ax.tick_params(axis='x', which='major', pad=20)
    ax.tick_params(axis='y', which='major', pad=20)
    ax.margins(0.25, tight=True) 
