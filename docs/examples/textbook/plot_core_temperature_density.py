import math
from amuse.lab import *
from matplotlib import pyplot
from distinct_colours import get_distinct

Second_Asymptotic_Giant_Branch = 6 | units.stellar_type
HeWhiteDwarf = 10 | units.stellar_type
    
def stellar_remnant(stellar):
    remnant = True
    if stellar.particles[0].stellar_type < HeWhiteDwarf \
       or stellar.particles[0].stellar_type > 11|units.stellar_type:
        remnant = False
    return remnant

colors = get_distinct(11)
def get_color_from_stellar_type(stype):
    st = stype.value_in(units.stellar_type)
    if st > 10: st = 10
    return colors[st]

def stellar_core_temperature_and_density(M, z=0.02, t_max=10|units.Myr):

    rho_core = []
    T_core = []
    color = []
    check = 0
    ncheck = 100

###BOOKLISTSTART1###
    stellar = MESA()
    stellar.parameters.metallicity = z
    star = stellar.particles.add_particle(Particle(mass=M))

    while not stellar_remnant(stellar) and star.age < t_max:

        stellar.evolve_model()

        nzones = star.get_number_of_zones()
        rhoc = star.get_density_profile(nzones)[0]
        Tc = star.get_temperature_profile(nzones)[0]

###BOOKLISTSTOP1###
	rho_core.append(rhoc.number)
        T_core.append(Tc.number)
        color.append(get_color_from_stellar_type(star.stellar_type))

        check += 1
        if check == ncheck:
            check = 0
            try:
                x = open('STOP')
                stop = True
            except:
                stop = False
            if stop: break

###BOOKLISTSTART2###
        print star.age.in_(units.Myr), rhoc, Tc, star.stellar_type

    stellar.stop()
###BOOKLISTSTOP2###

    return rho_core, T_core, color
    
if __name__ in ('__main__'):
    
    Mlist = [1, 10, 100] | units.MSun
    tmax = [1.2137e4, 22.6, 3.08064994] | units.Myr
    z = 0.02

    for i in range(len(Mlist)):
        M = Mlist[i]
        rhoc, Tc, color = stellar_core_temperature_and_density(M, z, tmax[i])
        size = 4*(math.log10(M.value_in(units.MSun))+1)
        print 'size =', size
        pyplot.scatter(rhoc, Tc, c=color, s=size)

    fontsize = 12
    pyplot.text(20., 1.5e7, '$1\,M_\odot$', fontsize=fontsize)
    pyplot.text(1.e2, 7.0e7, '$10\,M_\odot$', fontsize=fontsize)
    pyplot.text(4.e4, 1.6e9, '$100\,M_\odot$', fontsize=fontsize)
    
    pyplot.xlabel('core density [g/cm$^3$]')
    pyplot.ylabel('core temperature [K]')
    pyplot.xlim(1., 1.e7)
    pyplot.ylim(1.e7, 1.e10)
    pyplot.xscale('log')
    pyplot.yscale('log')

    save_file = 'plot_core_temperature_density.png'
    pyplot.savefig(save_file)
    print '\nSaved figure in file', save_file,'\n'
    #pyplot.show()
