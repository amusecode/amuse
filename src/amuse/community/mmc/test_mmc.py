from amuse.community import *
from amuse.test.amusetest import TestWithMPI
import coreradius
from .interface import mmcInterface
from .interface import mmc
from amuse.community.hermite0.interface import HermiteInterface

class mmcInterfaceTests(TestWithMPI):
    
    def test1(self):
        import numpy as np
        #instance = mmcInterface(redirection="file", redirect_file = "junk.txt")
        #instance = mmcInterface(redirection="none")
        instance = mmcInterface(redirection="null")
        instance.set_mmc_data_directory(instance.data_directory)
        instance.set_nt(10000)
        instance.set_imodel(2)#plummer
        instance.nonstandard_init()
        #instance.set_irun(10)
        instance.set_istart(2)

        n_total = instance.get_number_of_particles().n_
        instance.set_tcrit(0)
        
        R = instance.get_positions_3d(range(1,n_total))
        x = R['x']
        y = R['y']
        z = R['z']
        #V = instance.get_velocities_3d(range(1,n_total))
        #vx = V['vx']
        #vy = V['vy']
        #vz = V['vz']
        
        M  = instance.get_state(range(1,n_total))

        #        for tcrit in range(0,10000,100):                                      
        #                                                                              
        #            print instance.evolve(tcrit)                                      
        #                                                                              
        #            x, y, z = instance.get_positions_3d_(range(1,n_total))            
        #            M  = instance.get_state(range(1,n_total))                         
        #            x_core,y_core,z_core,rc = coreradius.coreradius(M.mass,x,y,z)     
        #            print tcrit, rc                                                   
        #                                                                              
        #            instance.set_tcrit(tcrit)                                         

        control = HermiteInterface()
        control.initialize_code()
        control.commit_parameters()
        control.new_particle(M.mass, np.ones(len(M.mass)), x, y, z)
        control.commit_particles()
        instance.stop()





"""

            #R = instance.get_positions_3d(range(1,n_total))
            #x = R['x']
            #y = R['y']
            #z = R['z']


        Unit_shell = instance.homogeneous_sphere_N(NT)
        R = np.diag(instance.get_state(range(NT)).r)
        print 

        f = open('data.dat','w')
        for i in R * np.matrix(Unit_shell):
            f.writelines("{0} {1} {2}\n".format(i[0,0],i[0,1],i[0,2]))

        f.close()

        NT = instance.get_number_of_particles().n_
        print NT
        instance.get_state(range(NT)).r

        Unit_shell = instance.homogeneous_sphere_N(NT)
        R = np.diag(instance.get_state(range(NT)).r)
        print 

        f = open('data2.dat','w')
        for i in R * np.matrix(Unit_shell):
            f.writelines("{0} {1} {2}\n".format(i[0,0],i[0,1],i[0,2]))

        f.close()
"""
