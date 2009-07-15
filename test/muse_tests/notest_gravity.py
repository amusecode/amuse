import unittest
import sys

from gravity.hermite0.muse_dynamics import Hermite
from gravity.BHTree.muse_dynamics import BHTree

class TestGravity(unittest.TestCase):
    def setupTwoParticles(self, gravity):
        state = gravity.dynamics_state()
        state.id = 1
        state.mass = 0.05
        state.radius = 0.25
        state.x,state.y,state.z = (0.4,0.4,0.4)
        state.vx,state.vy,state.vz = (0.2,0.2,0.2)
        gravity.add_particle(state)
        state.id = 2
        state.mass = 0.05
        state.radius = 0.25
        state.x,state.y,state.z = (1.0,1.0,1.0)
        state.vx,state.vy,state.vz = (0.2,0.4,0.4)
        gravity.add_particle(state)
        
        

class TestHermite(TestGravity):
    def test1(self):
        eps = 0.1
        gravity = Hermite(dt_dia=1.e9, eps2 = eps**2, flag_collision=0)
        gravity.setup_module()
        self.assertEquals(gravity.get_number(),0)
        self.setupTwoParticles(gravity)
        gravity.initialize_particles(0)
    
        totalEneryAtTime0 = gravity.get_kinetic_energy() + gravity.get_potential_energy()
        self.assertEquals(gravity.get_number(),2)
        self.assertAlmostEqual(gravity.get_time(), 0.0)
        gravity.evolve(0.25)
        self.assertTrue(gravity.get_time() < 0.25)
        self.assertTrue(gravity.get_time() + gravity.get_time_step() > 0.25)
        currenttime = gravity.get_time()
        gravity.evolve(gravity.get_time(), 1)
        self.assertAlmostEqual(gravity.get_time() ,currenttime)
        totalEneryAtTimeT = gravity.get_kinetic_energy() + gravity.get_potential_energy()
        self.assertAlmostEqual(totalEneryAtTimeT,totalEneryAtTime0)
        sys.stderr.write(str( gravity.get_n_steps()))

class TestBHTree(TestGravity):
    
    def test2(self):
        eps = 0.1
        gravity = BHTree(dt=0.015625,
                eps2=eps**2,
                use=1,
                theta=0.75,
                ncrit=1024,
                dt_dia=10)
        gravity.setup_module()
        self.setupTwoParticles(gravity)
        gravity.initialize_particles(0)
    
        totalEneryAtTime0 = gravity.get_kinetic_energy() + gravity.get_potential_energy()
        self.assertEquals(gravity.get_number(),2)
        self.assertAlmostEqual(gravity.get_time(), 0.0)
        gravity.evolve(0.25)
        print gravity.get_time()
        self.assertTrue(gravity.get_time() <= 0.25)
        self.assertTrue(gravity.get_time() + gravity.get_time_step() > 0.25)
        currenttime = gravity.get_time()
        gravity.evolve(gravity.get_time(), 1)
        self.assertAlmostEqual(gravity.get_time() ,currenttime)
        totalEneryAtTimeT = gravity.get_kinetic_energy() + gravity.get_potential_energy()
        self.assertAlmostEqual(totalEneryAtTimeT,totalEneryAtTime0)
        #sys.stderr.write(str( gravity.get_n_steps()))
