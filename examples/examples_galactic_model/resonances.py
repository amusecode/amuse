#---------------------------------------------------------------------------
# This script generates the corrotation and Lindbland resonances            | 
# as a function of radius using galactic_model.                             |
# For a complete explanation of the possible parameters and models          |
# included in galactic_model, we refer the reader to the User's manual.     |
#---------------------------------------------------------------------------

import numpy 
from amuse.units import units
from amuse.community.galactic_model.interface import BarAndSpirals3D
from matplotlib import rc 
import matplotlib.pyplot as plt

def plot_resonances(r, omega_c, ilr2, olr2, ilr4, olr4):
    rc('text', usetex=True)
    rc('font', family='serif')

    fig=plt.figure(0,figsize=(10,10))
    ax=fig.add_subplot(1,1,1)
    ax.plot(r.value_in(units.kpc), omega_c.value_in(units.kms/units.kpc), 'r-', label= 'Corotation resonance')
    ax.plot(r.value_in(units.kpc),ilr2.value_in(units.kms/units.kpc),"b--", label= 'ILRm2')
    ax.plot(r.value_in(units.kpc),olr2.value_in(units.kms/units.kpc),"g--", label= 'OLRm2')
    ax.plot(r.value_in(units.kpc),ilr4.value_in(units.kms/units.kpc),"b-", lw= 2, label= 'ILRm4')
    ax.plot(r.value_in(units.kpc),olr4.value_in(units.kms/units.kpc),"g-", lw=2, label= 'OLRm4')
    ax.legend(loc="upper right", ncol=1, shadow=False, fontsize=12)
    
   
    yticks = numpy.arange(0, 100, 10)
    xticks = numpy.arange(0, 15, 1)
    ax.set_xlim(0,15)
    ax.set_ylim(0, 100)
    ax.set_ylabel( "$\Omega$ $[kms^{-1}kpc^{-1}]$")
    ax.set_xlabel("R [Kpc]")
    plt.yticks(yticks)
    plt.xticks(xticks)
    plt.grid(True)
    plt.show()

if __name__ in('__main__', '__plot__'):
    
    MW= BarAndSpirals3D() # Axisymmetric model of the Galaxy with default values
    MW.commit_parameters()

    r= numpy.arange(0.001, 15, 0.005) |units.kpc
    phi= numpy.pi/4. #phi can take any value
    x, y, z= r*numpy.cos(phi), r*numpy.sin(phi), 0 | units.kpc
        
    circular_velocity= MW.get_velcirc(x,y,z)
    omega= circular_velocity/r # Angular velocity 
    epicyclic_frecuency= MW.get_epifreq(x,y,z)
    kappa= epicyclic_frecuency/r  # epicyclic velocity

    #Inner and Outer Lindblad resonances for 2,4 spiral arms
    ILRm2= omega - kappa/2.
    OLRm2= omega + kappa/2.
    ILRm4= omega - kappa/4.
    OLRm4= omega + kappa/4.

    
    plot_resonances(r, omega, ILRm2, OLRm2, ILRm4, OLRm4)
    
