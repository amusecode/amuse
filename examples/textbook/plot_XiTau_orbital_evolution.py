import pickle
from amuse.lab import *
#from read_orbital_elements import read_orbital_elements
from NonConservativeMassTransfer import orbital_separation_after_mass_transfer

from prepare_figure import *
from distinct_colours import get_distinct

def read_orbital_elements(filename="orbital_elements_data.pkl"):
    with open(filename, 'rb') as infile:
        results = pickle.load(infile)
        a_in = results["semimajor_axis_binary"]
        print(list(results.keys()))
    return results["times"], results["semimajor_axis_giant"], results["eccentricity_giant"], results["giant_mass"], results["ms1_mass"], results["ms2_mass"]

def limit_range(t, t_start=0|units.day, t_end=100|units.Myr):
    istart = 0
    for ti in t:
        istart+=1
        if ti>t_start:
            break
    iend = 0
    for i in t:
        iend += 1
        if ti>t_end:
            break
    return istart, iend

def retrieve_semi_major_axis_evolution(filename, nj, 
                                       t_start=0|units.day, 
                                       t_end=100|units.Myr):
    t, a_out, e_out, M_out, m1, m2  = read_orbital_elements(filename)
    m_in = m1+m2

    istart, iend = limit_range(t, t_start, t_end)

    t_new = [] | units.day
    an_out = [] | units.AU
    a_new = [] | units.AU
    en_out = [] 
    for i in range(len(a_out)):
        if i>=istart and i<=iend:
            ai = orbital_separation_after_mass_transfer(a_out[istart], 
                                                        M_out[istart], 
                                                        M_out[i], 
                                                        m_in[istart], 
                                                        m_in[i], nj)
            a_new.append(ai)
            t_new.append(t[i])
            an_out.append(a_out[i])
            en_out.append(e_out[i])
            if t[i]>t_end:
                break
    print(i, len(t), len(an_out), len(a_new))
    return t_new, an_out, a_new, en_out


def plot_XiTau(filename, nj):

    from matplotlib import pyplot
    colors = get_distinct(4)
    figure = pyplot.figure(figsize=(16, 12))
    ax = pyplot.gca()
    ax.minorticks_on() # switch on the minor ticks
    ax.locator_params(nbins=3)
    
    t, a_out, a_new, e_out = retrieve_semi_major_axis_evolution(filename, nj)
    from matplotlib import pyplot
    #pyplot.plot(t.value_in(units.yr), e_out, c=colors[0])
    pyplot.plot(t.value_in(units.yr), a_out.value_in(units.AU), c=colors[0])
    pyplot.plot(t.value_in(units.yr), a_new.value_in(units.AU), c=colors[1])

    t, a_out, a_new, e_out = retrieve_semi_major_axis_evolution(filename, 4.0)
    pyplot.plot(t.value_in(units.yr), a_new.value_in(units.AU), c=colors[2], lw=2)
    t, a_out, a_new, e_out = retrieve_semi_major_axis_evolution(filename, 8.0)
    pyplot.plot(t.value_in(units.yr), a_new.value_in(units.AU), c=colors[3], lw=2)
    
    pyplot.xlim(0, 12)
    pyplot.xlabel("t [year]")
    pyplot.ylabel("$a_{outer}$ [AU]")
    pyplot.savefig("fig_XiTau_orbital_separation")
    pyplot.show()

def plot_HD971331(filename):
    t_end = 4000|units.day
    t_start = 8000|units.day
    from matplotlib import pyplot

    colors = get_distinct(4)
    figure = pyplot.figure(figsize=(16, 12))
    ax = pyplot.gca()
    ax.minorticks_on() # switch on the minor ticks
    ax.locator_params(nbins=3)
    
    t, a_out, e_out, M_out, m1, m2  = read_orbital_elements(filename)
    pyplot.plot(t.value_in(units.yr), a_out.value_in(units.AU), c=colors[0])

    """
    t, a_out, a_new = retrieve_semi_major_axis_evolution(filename, nj=4.0, t_end=8000|units.day)
    pyplot.plot(t.value_in(units.yr), a_new.value_in(units.AU), c='r')
    pyplot.scatter(t[-1].value_in(units.yr), a_new[-1].value_in(units.AU), s=50, c=colors[1])

    """

    t, a_out, a_new = retrieve_semi_major_axis_evolution(filename, nj=6, t_end=2000|units.day)
    pyplot.plot(t.value_in(units.yr), a_new.value_in(units.AU), c='r')
    pyplot.scatter(t[-1].value_in(units.yr), a_new[-1].value_in(units.AU), s=50, c=colors[1])

    t, a_out, a_new = retrieve_semi_major_axis_evolution(filename, nj=4, t_start=2000|units.day, t_end=5000|units.day)
    pyplot.plot(t.value_in(units.yr), a_new.value_in(units.AU), c=colors[1])
    #pyplot.scatter(t[0].value_in(units.yr), a_new[0].value_in(units.AU), s=50, c=colors[1])
    #pyplot.scatter(t[-1].value_in(units.yr), a_new[-1].value_in(units.AU), s=50, c=colors[1])

    t, a_out, a_new = retrieve_semi_major_axis_evolution(filename, nj=2, t_start=5000|units.day, t_end=8000|units.day)
    print(t.value_in(units.yr), a_new.value_in(units.AU))
    #pyplot.plot(t.value_in(units.yr), a_new.value_in(units.AU), c=colors[1])
    #pyplot.scatter(t[0].value_in(units.yr), a_new[0].value_in(units.AU), s=50, c=colors[1])
    #pyplot.scatter(t[-1].value_in(units.yr), a_new[-1].value_in(units.AU), s=50, c=colors[1])

    t, a_out, a_new = retrieve_semi_major_axis_evolution(filename, nj=-7, t_start=8000|units.day)
    pyplot.plot(t.value_in(units.yr), a_new.value_in(units.AU), c=colors[1])
    #pyplot.scatter(t[0].value_in(units.yr), a_new[0].value_in(units.AU), s=50, c=colors[1])

    pyplot.xlabel("t [year]")
    pyplot.ylabel("$a_{giant}$ [AU]")
    pyplot.ylim(0.77, 0.8201)
    pyplot.xlim(0.0, 50.0)
    pyplot.show()

def new_option_parser():
    from amuse.units.optparse import OptionParser
    result = OptionParser()
    result.add_option("-f", 
                      dest="filename",default="XiTau_orbital_evolution.pkl")
#                      dest="filename",default="HD97131_orbital_evolution.pkl")    
    result.add_option("--nj", dest="nj",type="float",default=6)
    result.add_option("-t", unit=units.day,
                      dest="t_start",type="float",default=0|units.day)
    return result

def main(filename, nj, t_start=0|units.day):
    if "XiTau" in filename:
        plot_XiTau(filename, nj)
    else:
        plot_HD971331(filename) #HD97131_orbital_parameters.pkl


if __name__ == "__main__":
    options, arguments  = new_option_parser().parse_args()
    main(**options.__dict__)


