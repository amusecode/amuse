import unittest



class TestCollisions(unittest.TestCase):
    def createTwoParticles(self, interface):
        result = []
        state = interface.collisions_state()
        state.id = 1
        state.mass = 0.8
        state.radius = 0.955
        state.x,state.y,state.z = (3.7,0.0,0.0)
        state.vx,state.vy,state.vz = (1.0,1.0,0.0)
        state.age = 100
        result.append(state)
        state = interface.collisions_state()
        state.id = 2
        state.mass = 0.6
        state.radius = 0.65
        state.x,state.y,state.z = (5.0,0.0,1.0)
        state.vx,state.vy,state.vz = (-100.0,0.0,0.0)
        state.age = 100
        result.append(state)
        return result
    def createTwoParticlesB(self, interface):
        result = []
        state = interface.collisions_state()
        state.id = 1
        state.mass = 10.0
        state.radius = 3.7
        state.x,state.y,state.z = (0.0,0.0,1.0)
        state.vx,state.vy,state.vz = (1.0,0.0,0.0)
        state.age = 100
        result.append(state)
        state = interface.collisions_state()
        state.id = 2
        state.mass = 2
        state.radius = 1.3
        state.x,state.y,state.z = (5.0,0.0,1.0)
        state.vx,state.vy,state.vz = (-100.0,0.0,0.0)
        state.age = 100
        result.append(state)
        return result
        
        

class TestStickySpheres(TestCollisions):
    from collisions.sticky_spheres.muse_collisions import StickySpheres
    def test1(self):
        interface = self.StickySpheres()
        (star1,star2) = self.createTwoParticles(interface)
        newStar = interface.merge_stars(star1, star2)
        # note: makemeastar and stickysperes differ on radius and mass!
        self.assertAlmostEqual(newStar.radius,0.9662989,5) 
        self.assertAlmostEqual(newStar.mass,1.4,5)
        self.assertAlmostEqual(newStar.x,4.25714,5)
        self.assertAlmostEqual(newStar.y,0,5)
        self.assertAlmostEqual(newStar.z,0.428571,5)
class TestMakeMeAStar(TestCollisions):
    from collisions.mmas.muse_collisions import MakeMeAStar
    def test1(self):
        interface = self.MakeMeAStar()
        (star1,star2) = self.createTwoParticles(interface)
        newStar = interface.merge_stars(star1, star2)
        self.assertAlmostEqual(newStar.radius,2.23579,5)
        self.assertAlmostEqual(newStar.mass,1.30068,5)
        self.assertAlmostEqual(newStar.x,4.25714,5)
        self.assertAlmostEqual(newStar.y,0,5)
        self.assertAlmostEqual(newStar.z,0.428571,5)
class TestMakeMeAMassiveStar(TestCollisions):
    from collisions.mmams.muse_collisions import MMaMS
    def off2(self):
        interface = self.MMaMS()
        (star1,star2) = self.createTwoParticles(interface)
        newStar = interface.merge_stars(star1, star2)
        self.assertAlmostEqual(newStar.radius,2.23579,5)
        self.assertAlmostEqual(newStar.mass,1.30068,5)
        self.assertAlmostEqual(newStar.x,4.25714,5)
        self.assertAlmostEqual(newStar.y,0,5)
        self.assertAlmostEqual(newStar.z,0.428571,5)
        
