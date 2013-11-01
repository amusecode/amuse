"""
Fragmentation in a rotating protostar -- Boss & Bodenheimer (1979, BB79,
  http://adsabs.harvard.edu/abs/1979ApJ...234..289B)
  
  Calculates collapse of an isothermal cloud with SPH code fi. Cloud is initially 
  spherical, with an uniform density with a nonaxisymmetric perturbation of mode = 2 
  and with 50% amplitude. It rotates as a rigid body with a constant angular velocity 
  around the z-axis. The collapse leads to formation of a binary system.
  
  arguments:
    N -- number of particles
    Mcloud -- mass of the cloud
    Rcloud -- radius of the cloud
    omega -- angular velocity of the cloud
    density_perturb -- amplitude of the density perturbation
  
  Plots 10 snapshots of the simulations (bb79_rho_*.png) -- log10(rho [amu/cm**3]).
"""

import numpy
import string
import os
import sys

from matplotlib import pyplot
from matplotlib import rc

from amuse.units import nbody_system
from amuse.units import units

from amuse.ext.boss_bodenheimer import bb79_cloud
from amuse.ext.evrard_test import body_centered_grid_unit_cube

from amuse.community.fi.interface import Fi

rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
rc('text', usetex=True)

def get_grid_posvel(N=100, grid_size= 1 | units.parsec):
  """ 
  gives x,y(,z,vx,vy,vz) Cartesian coordinate grid (where z=0, vx=vy=vz=0)
  -- for plotting SPH quantities in the xy plane
  -- N points on each axis; size of the grid 'grid_size' for x and y 
  """
  
  cell_center_positions_1d = numpy.linspace(-0.5, 0.5, N)
  x,y = numpy.meshgrid(cell_center_positions_1d, cell_center_positions_1d)
  x = x.flatten()
  y = y.flatten()
    
  x *= grid_size
  y *= grid_size
  z = x * 0.0
    
  vx=0.0 * x / (1.0 | units.s)
  vy=0.0 * x / (1.0 | units.s)
  vz=0.0 * x / (1.0 | units.s)
  
  return x,y,z,vx,vy,vz

def plot_sph_rho(sph, N=100, grid_size= 1 | units.parsec, 
                 plot_name = "cloud_plot_rho.png", plot_title = "SPH density"):
  """
  plots the log10(SPH density) in the xy plane 
  input: sph -- set of sph particles 
    N -- number of grid points (in x and y)
    grid_size -- size of the grid (x and y)
    plot_name -- name for the output figure
    plot_title -- title of the plot
  """
  # getting the grid
  x,y,z,vx,vy,vz = get_grid_posvel(N, grid_size)
  
  # getting density from the SPH data
  rho,rhovx,rhovy,rhovz,rhoe=sph.get_hydro_state_at_point(x,y,z,vx,vy,vz)
  rho=rho.reshape((N,N))
  
  # plotting
  extent = (grid_size * (-0.5, 0.5, -0.5, 0.5)).value_in(units.parsec)
  pyplot.figure(figsize=(5,5.2))
  pyplot.imshow(numpy.log10(rho.value_in(units.amu/units.cm**3)),
                extent = extent)
  pyplot.title(plot_title)
  #pyplot.xlabel('x [pc]')
  #pyplot.ylabel('y [pc]')
  frame = pyplot.gca()
  frame.axes.get_xaxis().set_visible(False)
  frame.axes.get_yaxis().set_visible(False)
  cbar = pyplot.colorbar()
  cbar.set_label(r'$\log (\rho [\mathrm{amu}/\mathrm{cm}^3])$')
  pyplot.tight_layout()
  
  pyplot.savefig(plot_name)
  #pyplot.show()

def bb79_cloud_evolve(N=50000,
                      Mcloud=1. | units.MSun, 
                      Rcloud=3.2e16 | units.cm, 
                      omega=1.56e-12 | units.rad/units.s,
                      density_perturb=0.5,
                      t_total=8.3e11 | units.s):
  
  # mean density of the cloud
  rho_uni = Mcloud / (4./3.*numpy.pi*Rcloud**3)
  print " ** mean density = ", rho_uni.in_(units.g/units.cm**3)
  # free fall time of the cloud
  t_ff = numpy.sqrt(3.*numpy.pi / (32.*units.constants.G*rho_uni))
  print " ** free-fall time = ", t_ff.in_(units.yr)
  
  conv = nbody_system.nbody_to_si(Mcloud,Rcloud)
  sph = Fi(conv)
  
  # the initial conditions of BB79 
  parts_bb79 = bb79_cloud(targetN=N, omega=omega, rho_perturb=0.5, 
                          ethep_ratio=0.25, convert_nbody=conv, 
                          base_grid=body_centered_grid_unit_cube).result
  
  sph.gas_particles.add_particles(parts_bb79)
  
  sph.parameters.use_hydro_flag=True
  sph.parameters.isothermal_flag=True
  sph.parameters.integrate_entropy_flag=False
  sph.parameters.gamma=1
  sph.parameters.verbosity = 0
  sph.parameters.timestep=0.1*t_ff
  
  print " ** evolving to time: (end time = ~ {0:.3f} t_ff)".format(t_total/t_ff)
  
  # setting snapshots to be plotted
  nplot = 10
  if nplot > 1:
    plot_timestep = t_total / (nplot - 1)
  else:
    plot_timestep = t_total
  
  # evolution of the cloud
  for i in range(nplot):
    ttarget=i*plot_timestep
    t_tff = ttarget / t_ff
    
    sph.evolve_model(ttarget)
    
    plot_i = "bb79_rho_{0:03d}.png".format(i)
    tt_tff = "{0:.3f}".format(t_tff)
    title_i = "$%s\,t_{\mathrm{ff}}$" % (tt_tff)
     
    print "\t {0:.3f} t_ff -> {1}".format(t_tff,plot_i)
  
    plot_sph_rho(sph,N=300,grid_size= 0.025 | units.parsec, 
                 plot_name = plot_i, plot_title = title_i)
  
  sph.stop()
  
if __name__ in ("__main__","__plot__"):
  bb79_cloud_evolve(N=50000,
                    Mcloud=1. | units.MSun, 
                    Rcloud=3.2e16 | units.cm, 
                    omega=1.56e-12 | units.rad/units.s,
                    density_perturb=0.5,
                    t_total=8.3e11 | units.s) # ~1.5*t_ff
