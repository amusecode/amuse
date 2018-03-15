#from __future__ import print_function
import numpy
from amuse.units import units
from amuse.units import constants
from amuse.units import nbody_system
from amuse.ext.bridge import bridge
from amuse.community.phiGRAPE.interface import PhiGRAPE
from amuse.community.ph4.interface import ph4
from amuse.community.fi.interface import Fi
from amuse.community.bhtree.interface import BHTree
from amuse.community.gadget2.interface import Gadget2
from matplotlib import pyplot
from amuse.ic.kingmodel import new_king_model
from amuse.ext.galactics_model import new_galactics_model
from amuse.lab import new_powerlaw_mass_distribution

def make_king_model_cluster(nbodycode, N, W0, Mcluster,
                            Rcluster, parameters = []):

    converter=nbody_system.nbody_to_si(Mcluster,Rcluster)
    bodies=new_king_model(N,W0,convert_nbody=converter)

    code=nbodycode(converter)
    for name,value in parameters:
        setattr(code.parameters, name, value)
    code.particles.add_particles(bodies)
    return code

from prepare_figure import single_frame
from distinct_colours import get_distinct
from matplotlib.colors import LogNorm
def  plot_galaxy_and_stars(galaxy, stars):
    
    colors = get_distinct(3)
    single_frame('X [pc]', 'Y [pc]')
    xlim = 60
    pyplot.xlim(-xlim, xlim)
    pyplot.ylim(-xlim, xlim)
    ax = pyplot.gca()

    import numpy as np
    import pandas as pd
    from scipy import stats, integrate
    import matplotlib.pyplot as plt
    import seaborn as sns
    sns.set(color_codes=True)

    p = galaxy.select(lambda x: x<60|units.parsec,["x"])
    p = p.select(lambda x: x>-60|units.parsec,["x"])
    p = p.select(lambda y: y<60|units.parsec,["y"])
    p = p.select(lambda y: y>-60|units.parsec,["y"])
    x = p.x.value_in(units.parsec)
    y = p.y.value_in(units.parsec)
    sns.kdeplot(x, y, ax=ax)
    m = 100*numpy.sqrt(stars.mass/stars.mass.max())
    pyplot.scatter(stars.x.value_in(units.parsec), stars.y.value_in(units.parsec), c=colors[0], s=m, lw=0)
#    pyplot.show()
    pyplot.savefig("Fujii_Comparison_Figure")

from amuse.lab import new_plummer_model
def evolve_cluster_in_galaxy(N, W0, Rinit, tend, timestep, M, R):

    R_galaxy=0.1 | units.kpc
    M_galaxy=1.6e10 | units.MSun
    converter=nbody_system.nbody_to_si(M_galaxy, R_galaxy)
    galaxy=new_plummer_model(10000,convert_nbody=converter)

    print "com:", galaxy.center_of_mass().in_(units.kpc)
    print "comv:", galaxy.center_of_mass_velocity().in_(units.kms)
   
    print len(galaxy)
    galaxy_code = BHTree(converter, number_of_workers=2)
    galaxy_code.parameters.epsilon_squared = (0.01 | units.kpc)**2
    channe_to_galaxy = galaxy_code.particles.new_channel_to(galaxy)
    channe_to_galaxy.copy()
    galaxy_code.particles.add_particles(galaxy)
    inner_stars = galaxy.select(lambda r: r.length()<Rinit,["position"])
    Minner = inner_stars.mass.sum()
    print "Minner=", Minner.in_(units.MSun)
    print "Ninner=", len(inner_stars)
    vc_inner = (constants.G*Minner/Rinit).sqrt()

    converter=nbody_system.nbody_to_si(Mcluster,Rcluster)
    stars=new_king_model(N,W0,convert_nbody=converter)
    masses = new_powerlaw_mass_distribution(N, 0.1|units.MSun, 100|units.MSun, -2.35)
    stars.mass = masses
    stars.scale_to_standard(converter)
    
    stars.x += Rinit
    stars.vy += 0.8*vc_inner
    cluster_code=ph4(converter, number_of_workers=2)
    cluster_code.particles.add_particles(stars)
    channel_to_stars=cluster_code.particles.new_channel_to(stars)

    system=bridge(verbose=False)
    system.add_system(cluster_code, (galaxy_code,))
    system.add_system(galaxy_code, (cluster_code,))
    system.timestep = 0.1*timestep

    times=numpy.arange(0|units.Myr, tend, timestep)
    for i,t in enumerate(times):
        print "Time=", t.in_(units.Myr)
        channe_to_galaxy.copy()
        channel_to_stars.copy()

        inner_stars =  galaxy.select(lambda r: r.length()<Rinit,["position"])
        print "Minner=", inner_stars.mass.sum().in_(units.MSun)

        system.evolve_model(t,timestep=timestep)
    plot_galaxy_and_stars(galaxy, stars)
    galaxy_code.stop()
    cluster_code.stop()

if __name__ == "__main__":
    N=1024
    W0 = 3
    Rinit=50. | units.parsec
    timestep = 0.1 | units.Myr
#    endtime = 1.8| units.Myr
#    endtime = 1.4| units.Myr
    endtime = 2.35| units.Myr
    Mcluster = 5.e4 | units.MSun
    Rcluster = 0.8 | units.parsec
    evolve_cluster_in_galaxy(N, W0, Rinit, endtime, timestep, Mcluster, Rcluster)
        
