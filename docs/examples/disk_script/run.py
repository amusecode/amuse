from binary import circumbinary_disc_run
from amuse.units import units

# kepler 16 
circumbinary_disc_run(    tend=40. | units.yr,         # simulation time
                          Ngas=100000,                 # number of gas particles
                          m1=0.6897  | units.MSun,     # primary mass
                          m2=0.20255 | units.MSun,     # secondary mass
                          r1=0.6489  | units.RSun,     # primary radius
                          r2=0.22623 | units.RSun,     # secondary radius
                          ecc=0.15944,                 # binary orbit eccentricity
                          Pbinary=41.08 | units.day,   # binary orbit period
                          Rmin=2.3,                    # inner edge of initial disk 
                          Rmax=36.,                    # out edge of initial disk
                          q_out=12.,                   # outer disk Toomre Q parameter
                          discfraction=0.01,           # disk mass fraction
                          Raccretion=0.1 | units.AU,   # accretion radius for sink particle
                          dt_int=0.5 | units.day,      # timestep for gas - binary grav interaction (bridge timestep)
                          Pplanet=228.776 | units.day, # period of planet (makes the r-phi map rotate with this period)
                          densitypower=1.,             # surface density powerlaw
                          eosfreq=8,                   # times between EOS updates/sink particle checks
                          mapfreq=1,                   # time between maps ( in units of dt=eosfreq*dt_int)
                          Lmap=6. | units.AU,          # size of map
                          outputfreq=100,              # output snapshot frequency (every ... dt=eosfreq*dt_int)
                          outputdir='./kepler16',      # output directory
                          label='16',                  # label for run (only for terminal output)
                          alpha=0.5,                   # viscosity alpha
                          beta=1.,                     # viscosity beta
                          balsara=False   )            # balsara viscosity switch

