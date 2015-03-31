"""
    The analytic potential of the galaxy so it can be used in Bridge as an external potentail.
    Most equations are taken from B&T:
    Binney and Tremaine, Galactic Dynamics, Second Edition
"""
import numpy
from amuse.units import units, constants, quantities
from amuse.datamodel import Particle, Particles

class Abstract_Potential(object):
    def get_gravity_at_point(self, eps, x,y,z):
        """ derive the gravity from the potential """
        phi_0 = self.get_potential_at_point(eps, x,y,z)
        dpos = 0.001*(x**2+y**2+z**2).sqrt()
        phi_dx = self.get_potential_at_point(0,x+dpos,y,z) - phi_0
        phi_dy = self.get_potential_at_point(0,x,y+dpos,z) - phi_0
        phi_dz = self.get_potential_at_point(0,x,y, z+dpos) - phi_0
        return phi_dx/dpos, phi_dy/dpos, phi_dz/dpos

    def get_potential_at_point(self, eps, x, y, z):
        """ Abstract function, to be overwritten by subclass """
        pass

    def flattened_potential(self, x, y, z, a, b, mass):
        """
            Following eq. 2.69a of B&T
            a=0 gives plummer potential
            b=0 gives Kuzmin's potential for razor-thin disc
        """
        r_squared = x**2+y**2
        return constants.G * mass / (r_squared + (a + (z**2 + b**2).sqrt())**2).sqrt()

    def power_law_potential(self, r, alpha, r_0, mass_0):
        """ Following eq. 2.62 of B&T """
        rho_0 = mass_0 / (4./3. * numpy.pi * r_0**3)
        phi_0 = - constants.G * mass_0 / r_0
        v_circ_squared = 4 * numpy.pi * constants.G * rho_0 * r_0**alpha / (3 - alpha)

        if alpha == 2:
            phi_minus_phi_0 = - v_circ_squared * numpy.log(r/r_0)
        else:
            phi_minus_phi_0 = - v_circ_squared * (r_0**(2-alpha) - r**(2-alpha))/(alpha-2)

        return phi_minus_phi_0 + phi_0

    def point_mass_potential(self, r, mass):
        """ See eq. 2.34 of B&T """
        return constants.G * mass / r

class Disc_Bulge_Halo_Potential(Abstract_Potential):
    def halo_potential(self, x,y,z, Mc, Rc):
        """ TODO: Find the source for this potential -> McMillan & Portegies Zwart 2000?"""
        r=(x**2+y**2+z**2).sqrt()
        rr = (r/Rc)
        return -constants.G * (Mc/Rc)*(0.5*numpy.log(1 +rr**2) + numpy.arctan(rr)/rr)

    def get_potential_at_point(self, eps, x, y, z):
        disk = self.flattened_potential(x,y,z,
            0.0|units.kpc, 0.277|units.kpc, 1.12E+10|units.MSun)
        bulge = self.flattened_potential(x,y,z,
            3.7|units.kpc, 0.20|units.kpc, 8.07E+10|units.MSun)
        halo = self.halo_potential(x,y,z,
            Mc=5.0E+10|units.MSun, Rc=6.0|units.kpc)
        return disk + bulge + halo

class Galactic_Center_Potential_Kruijssen(Abstract_Potential):
    """
        Following Kruijssen et al 2014
        TODO: Get the real numbers for mass_0 and r_0,
        current ones are read from figure A1.
    """

    def __init__(self, beta=2.2, q=0.63, mass_0=1e8|units.MSun, r_0=42|units.parsec):
        self.q = q
        self.alpha = 3 - beta
        self.mass_0 = mass_0
        self.r_0 = r_0

    def get_potential_at_point(self, eps, x, y, z):
        r = (x**2 + y**2 + z**2/self.q**2).sqrt()
        return self.power_law_potential(r=r, alpha=self.alpha, mass_0=self.mass_0, r_0=self.r_0)

class Position_In_Potential(Abstract_Potential):
    """
        Wrapper around any other potential that has a test particle.
        Any call to get_potentail_at_point will shift the coordinates to
        put the center on the location of that test particle.
        The particle is put in a Particles set to allow channels to be used,
        however, only a single particle is allowed at a time.
    """
    def __init__(self, potential, particle=None):
        self.potential = potential

        if particle is None:
            particle = Particle()
            particle.position = [0., 0., 0.] | units.parsec
            particle.velocity = [0., 0., 0.] | units.kms

        self.particles = Particles()
        self.particles.add_particle(particle)

    @property
    def particle(self):
        return self.particles[0]

    @particle.setter
    def particle(self, particle):
        self.particles.remove_particles(self.particles)
        self.particles.add_particle(particle)

    def get_potential_at_point(self, eps, x, y, z):
        px, py, pz = self.particle.position
        return self.potential.get_potential_at_point(eps, x+px, y+py, z+pz)

