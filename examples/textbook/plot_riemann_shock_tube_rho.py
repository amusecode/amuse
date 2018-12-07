from amuse.lab import *
from amuse.plot import *

from prepare_figure import *
from distinct_colours import get_distinct

import csv

def read_csv(filename):
    ifile  = open(filename, "r")
    reader = csv.reader(ifile)

    x = []
    rho = []
    rownum = 0
    for row in reader:
        # Save header row.
        if rownum == 0:
            header = row
        elif rownum <=2 :
            units_str = row
        else:
            colnum = 0
            x.append(float(row[0]))
            rho.append(float(row[1]))
            for col in row:
                print '%-8s: %s' % (header[colnum], col)
                colnum += 1
        rownum += 1

    ifile.close()
    return x, rho

def plot_riemann_shock_tube_rho():

    x_label = "[length]"
    y_label = "[mass/length$^3$]"
    figure = single_frame(x_label, y_label, logx=False, logy=False, xsize=14, ysize=10)
    color = get_distinct(3)
        
    x, rho = read_csv("riemann_shock_tube_problem_exact.csv")
    pyplot.plot(x,rho, c=color[0])
    x, rho = read_csv("riemann_shock_tube_rho_fiN7.csv")
    pyplot.scatter(x, rho, c=color[1], s=100, marker="o", lw=0)
    x, rho = read_csv("riemann_shock_tube_problem_athenaN2.csv")
    pyplot.scatter(x, rho, c=color[2], s=100, marker="s", lw=0)

    pyplot.xlim(0.2,0.8)

#        pyplot.savefig("riemann_shock_tube_rho_"+model.name_of_the_code+".png")
    pyplot.savefig("riemann_shock_tube_rho")
    pyplot.show()


if __name__ == "__main__":
    plot_riemann_shock_tube_rho()
    
