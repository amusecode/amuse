"""
    The analytic potential of the galaxy so it can be used in Bridge as an external potentail.
    Most equations are taken from B&T:
    Binney and Tremaine, Galactic Dynamics, Second Edition
"""
import numpy
from amuse.units import units, constants, quantities
from amuse.datamodel import Particle, Particles
from amuse.support.exceptions import AmuseException
from io import StringIO 

class Abstract_Potential(object):
    def get_gravity_at_point(self, eps, x,y,z):
        """ derive the gravity from the potential """
        phi_0 = self.get_potential_at_point(eps, x,y,z)
        dpos = 0.001*(x**2+y**2+z**2).sqrt()
        phi_dx = self.get_potential_at_point(0,x+dpos,y,z) - phi_0
        phi_dy = self.get_potential_at_point(0,x,y+dpos,z) - phi_0
        phi_dz = self.get_potential_at_point(0,x,y, z+dpos) - phi_0
        return -phi_dx/dpos, -phi_dy/dpos, -phi_dz/dpos

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
        return -constants.G * mass / (r_squared + (a + (z**2 + b**2).sqrt())**2).sqrt()

    def power_law_potential(self, r, alpha, r_0, mass_0):
        """ Following eq. 2.62 of B&T """
        rho_0 = mass_0 / (4./3. * numpy.pi * r_0**3)
        phi_0 = - constants.G * mass_0 / r_0
        v_circ_squared = 4 * numpy.pi * constants.G * rho_0 * r_0**alpha / (3 - alpha)

        if alpha == 2:
            phi_minus_phi_0 = - v_circ_squared * numpy.log(r/r_0)
        else:
            phi_minus_phi_0 = - v_circ_squared * (r_0**(2-alpha) - r**(2-alpha))/(alpha-2)

        return -(phi_minus_phi_0 + phi_0)

    def point_mass_potential(self, r, mass):
        """ See eq. 2.34 of B&T """
        return -constants.G * mass / r

    def point_mass_gravity(self, r, mass, unit_vector):
        """ See eq. 2.27a of B&T """
        return -constants.G * mass / r**2 * unit_vector

class Disc_Bulge_Halo_Potential(Abstract_Potential):
    def halo_potential(self, x,y,z, Mc, Rc):
        """ TODO: Find the source for this potential -> McMillan & Portegies Zwart 2000?"""
        r=(x**2+y**2+z**2).sqrt()
        rr = (r/Rc)
        return constants.G * (Mc/Rc)*(0.5*numpy.log(1 +rr**2) + numpy.arctan(rr)/rr)

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
        Following Kruijssen et al 2014, which uses the enclosed mass
        profile from Launhardt et al 2002.

        Note that the mass profile only extends to 487.9 parsec from
        the galactic center, (and only to 487.9 * 0.63 parsec in the z direction).
        Outside this range, a different potential should be used.
    """

    def __init__(self, q=0.63):
        self.q = q
        self.load_table()

    def load_table(self, rescale=True):
        table = """
            # enclosed mass profile from Launhardt et al 2002,
            # recreated and provided by Kruijssen 2014
            # radius   enclosed mass
            # (parsec) (MSun)
            0.0     0.0
            0.6261	3298000
            0.6945	3429000
            0.7751	3636000
            0.8756	3855000
            0.9712	4088000
            1.077	4335000
            1.187	4552000
            1.293	4826000
            1.408	5019000
            1.534	5374000
            1.671	5754000
            1.820	6101000
            1.994	6469000
            2.198	6995000
            2.424	7563000
            2.705	8178000
            2.964	8842000
            3.268	9561000
            3.625	1.033E7
            3.996	1.128E7
            4.406	1.232E7
            4.798	1.332E7
            5.131	1.413E7
            5.487	1.498E7
            6.013	1.558E7
            6.752	1.635E7
            7.489	1.717E7
            8.409	1.803E7
            9.327	1.894E7
            10.28	2.028E7
            11.54	2.214E7
            13.04	2.441E7
            14.91	2.718E7
            17.16	3.026E7
            19.98	3.402E7
            22.99	3.939E7
            26.13	4.516E7
            29.35	5.229E7
            32.35	5.995E7
            35.02	6.806E7
            38.61	8.036E7
            43.09	9.865E7
            47.51	1.199E8
            51.74	1.429E8
            57.75	1.789E8
            64.84	2.329E8
            71.92	2.973E8
            81.25	3.947E8
            91.79	5.188E8
            99.97	6.246E8
            109.5	7.743E8
            120.0	9.142E8
            133.9	1.068E9
            152.2	1.237E9
            179.5	1.489E9
            206.5	1.741E9
            243.4	2.056E9
            283.5	2.289E9
            332.3	2.573E9
            382.3	2.893E9
            413.8	3.098E9
            456.2	3.449E9
            487.9	3.694E9
            1e6  	3.694E9
        """
        stream = StringIO(table)
        radius, enclosed_mass = numpy.loadtxt(stream, unpack=True)
        if rescale:
            """ See footnote 18 at the bottom of page 1076 of Kruijssen """
            factor = 8.3/8.5
            radius *= factor
            enclosed_mass *= factor**2
        self.radius = radius | units.parsec
        self.enclosed_mass_profile = enclosed_mass | units.MSun

    def enclosed_mass(self, r):
        try:
            index = quantities.searchsorted(self.radius, r)
        except ValueError:
            """ This error is usually thrown when r has dimension > 1 """
            shape = r.shape
            r_flat = r.flatten()
            index = quantities.searchsorted(self.radius, r_flat)
            index.reshape(shape)

        mass_below = self.enclosed_mass_profile[index-1]
        mass_above = self.enclosed_mass_profile[index]
        radius_below = self.radius[index-1]
        radius_above = self.radius[index]

        # Linear interpolation in log space
        log_m_over_mb = numpy.log(mass_above/mass_below) * numpy.log(r/radius_below) / numpy.log(radius_above/radius_below)
        enclosed_mass = numpy.nan_to_num(numpy.exp(log_m_over_mb)) * mass_below

        return enclosed_mass

    def get_potential_at_point(self, eps, x, y, z):
        """
            Note that this potential is not entirely consistent with
            get_gravity_at_point (which should be used) because
            there a second coordinate transformation is used to flatten
            the "potential".
        """
        r = (x**2 + y**2 + z**2/self.q**2).sqrt()
        mass = self.enclosed_mass(r)
        # missing a term with the integrated contribution from outside shells
        return self.point_mass_potential(r, mass)

    def get_gravity_at_point(self, eps, x,y,z):
        """
            Overwrites the default to add a second coordinate transformation.
        """
        r = (x**2 + y**2 + z**2/self.q**2).sqrt()
        mass = self.enclosed_mass(r)
        unit_vector = []|x.unit
        for var in (x, y, z):
            unit_vector.append(var)
        unit_vector = unit_vector / (x**2 + y**2 + z**2).sqrt()
        unit_vector[2] *= 1./self.q**2
        return self.point_mass_gravity(r, mass, unit_vector)

class Position_In_Potential(Abstract_Potential):
    """
        Wrapper around any other potential that has a test particle.
        Any call to get_potential_at_point will shift the coordinates to
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

