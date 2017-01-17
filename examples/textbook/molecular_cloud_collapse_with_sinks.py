"""
  example of molecular cloud evolution with explictly 
  split SPH and grav evolution

  Initial condition is a smooth spherical cloud with random velocities
  as in Bonnell et al. (2003)  
  
"""  

from cooling_class import SimplifiedThermalModel, SimplifiedThermalModelEvolver

import numpy
  
from matplotlib import pyplot 

from amuse.lab import *
from amuse.ext.molecular_cloud import molecular_cloud
from amuse.ext.evrard_test import body_centered_grid_unit_cube

from prepare_figure import single_frame, figure_frame, set_tickmarks
from distinct_colours import get_distinct

from hydrodynamics_class import Hydro

def make_map(sph,N=100,L=1):

    x,y=numpy.indices( ( N+1,N+1 ))

    x=L*(x.flatten()-N/2.)/N
    y=L*(y.flatten()-N/2.)/N
    z=x*0.
    vx=0.*x
    vy=0.*x
    vz=0.*x

    x=units.parsec(x)
    y=units.parsec(y)
    z=units.parsec(z)
    vx=units.kms(vx)
    vy=units.kms(vy)
    vz=units.kms(vz)

    rho,rhovx,rhovy,rhovz,rhoe=sph.get_hydro_state_at_point(x,y,z,vx,vy,vz)
    rho=rho.reshape((N+1,N+1))

    return rho

def write_output(filename, parts, conv):

    particles_nbody = ParticlesWithUnitsConverted(parts, conv.as_converter_from_nbody_to_si())

    write_set_to_file(particles_nbody, filename, "txt", attribute_names= ('rho', 'mass', 'x', 'y', 'z','vx', 'vy', 'vz'))

    return 0

def plot_hydro(time, sph, i, L=10):
    x_label = "x [pc]"
    y_label = "y [pc]"
    fig = single_frame(x_label, y_label, logx=False, logy=False, xsize=12, ysize=12)

    gas = sph.code.gas_particles
    dmp = sph.code.dm_particles
    #rho=make_map(sph,N=200,L=L)
#    pyplot.imshow(numpy.log10(1.e-5+rho.value_in(units.amu/units.cm**3)), extent=[-L/2,L/2,-L/2,L/2],vmin=1,vmax=5, origin="lower")
#    pyplot.hist2d(gas.x.value_in(units.parsec),
#                  gas.y.value_in(units.parsec),
#                  (100, 100), cmap=pyplot.cm.jet)
#    pyplot.scatter(gas.x.value_in(units.parsec),
#                   gas.y.value_in(units.parsec), s=500, lw=0, alpha=0.1, c='r')
    x = gas.x.value_in(units.parsec)
    y = gas.y.value_in(units.parsec)
    z = gas.rho.value_in(units.g/units.cm**3)
    #z = gas.u.value_in(units.g/units.cm**3),
    #cax = pyplot.tripcolor(gas.x.value_in(units.parsec),
    #                       gas.y.value_in(units.parsec),
    #                       z,
    #                       shading='flat') #, cmap=plt.cm.rainbow)
    #                       #shading='gouraud') #, cmap=plt.cm.rainbow)
    import matplotlib.tri as tri
    triang = tri.Triangulation(x, y)
    xmid = x[triang.triangles].mean(axis=1)
    ymid = y[triang.triangles].mean(axis=1)
    min_radius = 0.01
    mask = numpy.where(xmid*xmid + ymid*ymid < min_radius*min_radius, 1, 0)
    triang.set_mask(mask)
#    pyplot.tripcolor(triang, z, shading='gouraud', cmap=pyplot.cm.rainbow)
    pyplot.tripcolor(triang, z, cmap=pyplot.cm.rainbow)
    pyplot.xlim(-1, 1)
    pyplot.ylim(-1, 1)
    if len(dmp):
        m = 2.0*dmp.mass/dmp.mass.max()
        print m
        pyplot.scatter(dmp.x.value_in(units.parsec), dmp.y.value_in(units.parsec), c='k', s=m) 
    pyplot.savefig("GMCsink_"+str(i)+".png")

def run_molecular_cloud(N=100, Mcloud=100. | units.MSun, Rcloud=1. | units.parsec):

    conv = nbody_system.nbody_to_si(Mcloud,Rcloud)

    rho_cloud = 3.*Mcloud/(4.*numpy.pi*Rcloud**3)
    print rho_cloud
    tff = 0.5427/numpy.sqrt(constants.G*rho_cloud)
    print "t_ff=", tff.value_in(units.Myr), 'Myr'

    dt = 0.05 | units.Myr
    tend=1.0 | units.Myr

    parts=molecular_cloud(targetN=N,convert_nbody=conv,
            base_grid=body_centered_grid_unit_cube, seed=100).result
    
    sph = Hydro(Fi, parts)
    #sph = Hydro(Gadget2, parts)
    eps = 0.1 | units.parsec

    expected_dt = 0.2*numpy.pi*numpy.power(eps, 1.5)/numpy.sqrt(constants.G*Mcloud/N)

    print "dt_exp=", expected_dt.value_in(units.Myr)
    print "dt=", dt
    print "eps=", sph.parameters.gas_epsilon.in_(units.parsec)

    i=0
    L=6
    E0 = 0.0
    ttarget = 0.0 | units.Myr

    plot_hydro(ttarget, sph, i, L)

    while ttarget < tend:
        ttarget=float(i)*dt
        print "Evolve to time=", ttarget.in_(units.Myr)
#        print "N=", len(sph.gas_particles), len(sph.dm_particles)
#        print "Masses of dm particles:", sph.dm_particles.mass.in_(units.MSun)

        sph.evolve_model(ttarget)
        E = sph.gas_particles.kinetic_energy()+sph.gas_particles.potential_energy() + sph.gas_particles.thermal_energy()
        E_th = sph.gas_particles.thermal_energy()
        if i==0:
            E0 = E
        Eerr = (E-E0)/E0
        print 'energy=', E, 'energy_error=', Eerr, 'e_th=', E_th
        print "maximal_density:",parts.rho.max().in_(units.MSun/units.parsec**3)

        """
        filename = 'm400k_r10pc_e01_'+ str(i).zfill(2) + '.dat'
        print filename
        parts_sorted = parts.sorted_by_attribute('rho')
        write_output(filename, parts_sorted, conv)        
        """

        plot_hydro(ttarget, sph, i, L)
        i=i+1

    sph.stop()
    return parts
  
if __name__ in ("__main__","__plot__"):

    parts = run_molecular_cloud(1000, Mcloud=100. | units.MSun, Rcloud=0.5 | units.parsec)


    
