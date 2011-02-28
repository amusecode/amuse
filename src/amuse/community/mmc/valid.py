import numpy as np
from amuse.community import *
from amuse.test.amusetest import TestWithMPI
from amuse.community.hermite0.interface import HermiteInterface
from amuse.community.phiGRAPE.interface import PhiGRAPEInterface
from amuse.community.bhtree.interface import BHTreeInterface

from amuse.ext.plummer import MakePlummerModel
import time as systime

import coreradius
from interface import mmcInterface
from interface import mmc

def plummer(x):
  plummer = MakePlummerModel(x)
  mass,pos,vel=plummer.new_model()

  mass=mass[0:,0]
  x=pos[0:,0]
  y=pos[0:,1]
  z=pos[0:,2]

  vx=vel[0:,0]
  vy=vel[0:,1]
  vz=vel[0:,2]
  radius=mass*0.

  return mass,radius,x,y,z,vx,vy,vz

def sort_by_radii(mass, r, vr, vt):
    indices = [r.argsort(),]
    return mass[indices], r[indices], vr[indices], vt[indices]

def crossing_time(r, vr, vt):
    R = np.sum(r)
    V = np.sqrt(np.sum(vr**2 + vt**2))
    return 2*R/V

def compare(mass, r, vr, vt, mass_, r_, vr_, vt_):
    for i, v in enumerate(mass):
        print "{0} {1} {2} {3} {4} {5} {6} {7}".format(mass[i], r[i], vr[i], vt[i], mass_[i], r_[i], vr_[i], vt_[i])

def tests_tocartesian():
  x,y,z,vx,vy,vz,ex,ey,ez = to_cartesian(np.ones(100), np.zeros(100),np.ones(100))
  for i, x in enumerate(x):
    print x,y[i],z[i],vx[i],vy[i],vz[i],ex[i],ey[i],ez[i]

def example_M67():

    #in this example we take M67 from mmc and evolve it in both
    #mmc and hermite. We compare coreradii

    mmc = mmcInterface(redirection="null")
    mmc.set_mmc_data_directory(mmc.data_directory)
    mmc.set_nt(2000)
    n_total = mmc.get_number_of_particles().n_

    mmc.set_imodel(2) #M67
    mmc.amuse_input()
    mmc.nonstandard_init()

    mmc.set_tau0(0.001)

    mass,radius,x,y,z,vx,vy,vz = plummer(n_total)
    r_, vr_, vt_ = mmc.phase_to_polar(x, y, z, vx, vy, vz)

    mmcstate  = mmc.get_state(range(1,n_total+1))

    x, y, z, vx, vy, vz,ex,ey,ez = mmc.phase_to_cartesian(mmcstate.r, mmcstate.vr, mmcstate.vt)

    mass, r_, vr_, vt_ = sort_by_radii(mass, r_, vr_, vt_)

    mmc.commit_particles()

    n_total = mmc.get_number_of_particles().n_
    mmc.set_tcrit(0)

    control = HermiteInterface(number_of_workers=3)
    control.initialize_code()                                                                    
    control.commit_parameters()                
    res = control.new_particle(mass,                                                                 
                         1/n_total * np.ones(len(mass)),
                         x, y, z,
                         vx, vy, vz)

    control_parts = res['index_of_the_particle']
    control.commit_particles()                                                                   

    mmc.set_istart(1)    

    for time_end in np.arange(0.001,1,0.001):
        n_total = mmc.get_number_of_particles().n_
        x, y, z = mmc.get_positions_3d(range(1, n_total+1))                           
        M  = mmc.get_state(range(1, n_total+1))                                        
        tcross = crossing_time(M.r, M.vr, M.vt)
        x_core,y_core,z_core,rc = coreradius.coreradius(M.mass,x,y,z)                    
        tic = systime.clock()
        mmc.evolve(time_end)          
        toc = systime.clock()
        time =  mmc.get_time().time
        tcrit_control =  mmc.get_tcrit().termination_time_units_crossing_time
        timet =  mmc.get_timet().time
        tcr = mmc.get_crossing_time().tcr

        control.evolve(time_end) 

        H = control.get_state(control_parts)
        control_xcore,control_ycore,control_zcore,control_rc = \
            coreradius.coreradius(H.mass,H.x,H.y,H.z)                    
        
        print time, time_end, rc, control_rc, toc-tic

    mmc.stop()
    control.stop()

def example_plummer():
    pass


if __name__ == '__main__':

    example_M67() # M67 from mmc

    example_plummer() # plummer from amuse
