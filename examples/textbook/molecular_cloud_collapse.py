"""
  Example of molecular cloud evolution with explictly split SPH
  and gravity evolution

  Initial condition is a smooth spherical cloud with random velocities
  as in Bonnell et al. (2003)  
  
"""  

import numpy
  
from matplotlib import pyplot 

from amuse.lab import *
from amuse.ext.molecular_cloud import molecular_cloud
from amuse.ext.evrard_test import body_centered_grid_unit_cube

from prepare_figure import single_frame, figure_frame, set_tickmarks
from distinct_colours import get_distinct

def make_map(sph, N=100, L=1):

    x,y = numpy.indices((N+1,N+1))

    x = L*(x.flatten()-N/2.)/N
    y = L*(y.flatten()-N/2.)/N
    z = x*0.
    vx = 0*x
    vy = 0*x
    vz = 0*x

    x = units.parsec(x)
    y = units.parsec(y)
    z = units.parsec(z)
    vx = units.kms(vx)
    vy = units.kms(vy)
    vz = units.kms(vz)

    rho,rhovx,rhovy,rhovz,rhoe = sph.get_hydro_state_at_point(x,y,z,vx,vy,vz)
    rho = rho.reshape((N+1,N+1))

    return rho

def write_output(filename, parts, conv):

    particles_nbody \
        = ParticlesWithUnitsConverted(parts,
                                      conv.as_converter_from_nbody_to_si())

    write_set_to_file(particles_nbody, filename, "txt",
                      attribute_names=('rho', 'mass', 'x', 'y', 'z',
                                       'vx', 'vy', 'vz'))

    return 0

def plot_hydro(time, sph, i, L=10):
    x_label = "x [pc]"
    y_label = "y [pc]"
    fig = single_frame(x_label, y_label, logx=False, logy=False,
                       xsize=12, ysize=12)
    rho = make_map(sph,N=200,L=L)
    cax = pyplot.imshow(numpy.log10(1.e-5+rho.value_in(units.amu/units.cm**3)),
                        cmap="jet",
                        extent=[-L/2,L/2,-L/2,L/2],vmin=2,vmax=5)
    if i==5:
        cbar = fig.colorbar(cax, ticks=[2, 3, 4, 5], orientation='vertical', fraction=0.05)
        cbar.ax.set_yticklabels([2, " ", " ", 5])  # horizontal colorbar
        cbar.set_label('log projected density [$amu/cm^3$]', rotation=270)

    """
    rhomin = numpy.log10(rho.value_in(units.amu/units.cm**3)).min()
    rhomax = numpy.log10(rho.value_in(units.amu/units.cm**3)).max()
    cax = pyplot.imshow(numpy.log10(1.e-5+rho.value_in(units.amu/units.cm**3)),
                        extent=[-L/2,L/2,-L/2,L/2],vmin=rhomin,vmax=rhomax)
    rhomid = 0.5*(rhomin + rhomax)
    print rhomin, rhomid, rhomax
    cbar = fig.colorbar(cax, ticks=[rhomin, rhomid, rhomax], orientation='vertical', fraction=0.045)
#    cbar.ax.set_yticklabels(['Low', ' ', 'High'])  # horizontal colorbar
    low = "%.2f" % rhomin
    mid = "%.1f" % rhomid
    mx = "%.2f" % rhomax
    cbar.ax.set_yticklabels([low, mid, mx])  # horizontal colorbar
    cbar.set_label('projected density [$amu/cm^3$]', rotation=270)
    """
    
    pyplot.savefig("GMC_"+str(i)+".png")

def run_molecular_cloud(N=100, Mcloud=100. | units.MSun,
                        Rcloud=1. | units.parsec):

    conv = nbody_system.nbody_to_si(Mcloud,Rcloud)

    rho_cloud = 3.*Mcloud/(4.*numpy.pi*Rcloud**3)
    print(rho_cloud)
    tff = 0.5427/numpy.sqrt(constants.G*rho_cloud)
    print("t_ff=", tff.value_in(units.Myr), 'Myr')

    dt = 5.e-2 | units.Myr
    t_end = 1.0 | units.Myr

    parts = molecular_cloud(targetN=N,convert_nbody=conv,
                            base_grid=body_centered_grid_unit_cube,
                            seed=100).result

    sph = Fi(conv)
    sph.parameters.use_hydro_flag = True
    sph.parameters.radiation_flag = False
    sph.parameters.gamma = 1
    sph.parameters.isothermal_flag = True
    sph.parameters.integrate_entropy_flag = False
    sph.parameters.timestep = dt  
    sph.parameters.verbosity = 0
    sph.parameters.eps_is_h_flag = False    # h_smooth is constant
    eps = 0.1 | units.parsec
    sph.parameters.gas_epsilon = eps
    sph.parameters.sph_h_const = eps
    parts.h_smooth= eps

    print('eps-h flag', sph.get_eps_is_h(), sph.get_consthsm())

    expected_dt = 0.2*numpy.pi*numpy.power(eps, 1.5) \
                   / numpy.sqrt(constants.G*Mcloud/N)

    print("dt_exp=", expected_dt.value_in(units.Myr))
    print("dt=", dt)
    print("eps=", sph.parameters.gas_epsilon.in_(units.parsec))

    sph.gas_particles.add_particles(parts)

    channel_from_sph_to_parts= sph.gas_particles.new_channel_to(parts)
    channel_from_parts_to_sph= parts.new_channel_to(sph.gas_particles)
    
    i = 0
    L = 6
    E0 = 0.0
    ttarget = 0.0 | units.Myr

    plot_hydro(ttarget, sph, i, L)

    while ttarget < t_end:
        ttarget = float(i)*dt
        print(ttarget)
        sph.evolve_model(ttarget, timestep=dt)
        E = sph.gas_particles.kinetic_energy() \
             + sph.gas_particles.potential_energy() \
             + sph.gas_particles.thermal_energy()
        E_th = sph.gas_particles.thermal_energy()
        if i == 0:
            E0 = E
        Eerr = (E-E0)/E0
        print('energy=', E, 'energy_error=', Eerr, 'e_th=', E_th)
        channel_from_sph_to_parts.copy()
        plot_hydro(ttarget, sph, i, L)
        i += 1

    sph.stop()
    return parts
  
if __name__ in ("__main__","__plot__"):

    parts = run_molecular_cloud(
        1000,
        Mcloud=10000. | units.MSun,
        Rcloud=3. | units.parsec,
    )
