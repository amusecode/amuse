import matplotlib.pyplot as plt
import numpy as np

from amuse.units import units
from amuse.community.mesa_r15140.interface import MESA
from amuse import datamodel

stellar_evolution = MESA()

masses=[2.0] | units.MSun
stars = datamodel.Particles(len(masses), mass=masses)

stars = stellar_evolution.native_stars.add_particles(stars)
star = stars[0]

star.evolve_one_step()

age = star.get_history('star_age') # Years

fig=plt.figure(figsize=(12*4,9*4))
ax1 = fig.add_subplot(221)
ax2 = fig.add_subplot(222)
ax3 = fig.add_subplot(223)
ax4 = fig.add_subplot(224)

logL= []
logT =[]
logTc=[]
logRhoc=[]

while age < 4.6*10**9:
    star.evolve_one_step()
    age = star.get_history('star_age') # Years

    logL.append(star.get_history('log_L'))
    logT.append(star.get_history('log_Teff'))

    logTc.append(star.get_history('log_center_T'))
    logRhoc.append(star.get_history('log_center_Rho'))

    ax1.cla()
    ax2.cla()
    ax3.cla()
    ax4.cla()

    ax1.plot(logT,logL)
    ax1.scatter(logT[-1],logL[-1],s=50)
    ax1.set_xlim(ax1.get_xlim()[::-1])
    ax1.set_xlabel(r'$\log \left(T_{eff}/K\right)$')
    ax1.set_ylabel(r'$\log \left(L/L_{\odot}\right)$')

    ax1.set_title('Model Number='+str(int(star.get_history('model_number'))))
    ax2.set_title('Star age='+str(star.get_history('star_age')))

    ax2.plot(logRhoc,logTc)
    ax2.scatter(logRhoc[-1],logTc[-1],s=50)
    ax2.set_ylabel(r'$\log \left(T_{c}/K\right)$')
    ax2.set_xlabel(r'$\log \left(\rho_{c}/\left(g\ cm^{-3}\right)\right)$')

    lgt = star.get_profile('logT')
    lgr = star.get_profile('logRho')
    ax3.set_xlabel(r'$\log \left(\rho/\left(g\ cm^{-3}\right)\right)$')
    ax3.set_ylabel(r'$\log \left(T/K\right)$')

    ax3.plot(lgr,lgt)

    m = star.get_profile('mass')
    ax4.set_yscale('log')
    for i in range(1,star.get_number_of_species()+1):
        name = star.get_name_of_species(i)
        species = star.get_profile(name)
        ax4.plot(m,species,label=name)
    ax4.legend(loc=0)   
    ax4.set_ylim(10**-5,1.0)

    ax4.set_xlabel(r'$m/M_{\odot}$')
    ax4.set_ylabel('Abundance')

    plt.pause(0.01)



