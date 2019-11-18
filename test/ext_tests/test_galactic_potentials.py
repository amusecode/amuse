from amuse.test import amusetest
from amuse.units import units, nbody_system, constants
from amuse.ext.galactic_potentials import NFW_profile, MiyamotoNagai_profile, Plummer_profile, \
    PowerLawCutoff_profile, MWpotentialBovy2015, scipy_imported

import numpy

class TestNFWProfile(amusetest.TestCase):
    def test1(self):
        """
        See the textbook by L.Aguilar, chapter 3:
        ftp://ftp.crya.unam.mx/pub/luisfr/laguilar/GH05_Aguilar.pdf
        """
        rho0 = 1.e3|units.MSun/units.parsec**3
        rs = 6.6|units.kpc
        nfw_profile = NFW_profile(rho0, rs)
        m0 = 4.*numpy.pi*rho0*rs**3
        phi0 = -4.*numpy.pi*constants.G*rho0*rs**2
        ar0 = -phi0/rs
        
        # relative mass enclosed at rs: m_at_rs / m0 ~ 0.193147
        m_at_rs = nfw_profile.enclosed_mass(rs)
        self.assertAlmostEqual(m_at_rs/m0, 0.193147, 5)
        
        # relative mass enclosed at 5.3054*rs ~ m0
        m_eq_one = nfw_profile.enclosed_mass(5.3054*rs,)
        self.assertAlmostEqual(m_eq_one/m0, 1.0, 5)
        
        # relative mass density, rho(rs) ~ rho0/4
        rho_rs = nfw_profile.mass_density(rs)
        self.assertAlmostEqual(rho0/rho_rs, 4.0, 5)
        
        # relative gravitational potential at (r/rs)->0 is approaching 1
        phi_at_0 = nfw_profile.get_potential_at_point(0.|units.m,0.|units.m,0.|units.m,1.e-10|units.kpc)
        self.assertAlmostEqual(phi_at_0/phi0, 1.0, 4)
        
        # relative force at (r/rs)->0 is approaching -1/2
        ar_at_0 = nfw_profile.radial_force(1.e-5|units.kpc)
        self.assertAlmostEqual(ar_at_0/ar0, -0.5, 4)
        
        # relative force at (r/rs)->inf is approaching 0
        ar_at_inf = nfw_profile.radial_force(1.e10|units.kpc)
        self.assertAlmostEqual(ar_at_inf/ar0, 0.0, 5)
        
        # relative circular velocity has maximum at r/r_s=2.16258
        vc_eq_max = nfw_profile.circular_velocity(2.16258*rs)
        vc_at_r_lt_max = nfw_profile.circular_velocity(2.16248*rs)
        vc_at_r_gt_max = nfw_profile.circular_velocity(2.16268*rs)
        self.assertTrue(vc_at_r_lt_max < vc_eq_max)
        self.assertTrue(vc_at_r_gt_max < vc_eq_max)
          
class TestPlummerProfile(amusetest.TestCase):
    def test1(self):
        mass = 6.e6|units.MSun
        a = 6.|units.parsec
        plummer_profile = Plummer_profile(mass,a)
        rho0 = mass/(4./3.*numpy.pi*a**3)
        phi0 = -constants.G*mass/a
        
        # enclosed mass at R>>a is total mass
        m_tot = plummer_profile.enclosed_mass(a*1.e5)
        self.assertAlmostEqual(m_tot/mass, 1.0, 5)
        
        # mass density at the center
        rho_cen = plummer_profile.mass_density(0.|units.m)
        self.assertAlmostEqual(rho_cen/rho0, 1.0, 5)
        
        # potential at r=a is phi0/sqrt(2)
        phi_at_a = plummer_profile.get_potential_at_point(0.|units.m,0.|units.m,0.|units.m,a)
        self.assertAlmostEqual(phi_at_a/phi0*numpy.sqrt(2.), 1.0, 5)
        
class TestMiyamotoNagaiProfile(amusetest.TestCase):
    def test1(self):
        mass = 6.6e10|units.MSun
        a = 3.33|units.kpc
        b = 0.666|units.kpc
        profile = MiyamotoNagai_profile(mass,a,b)
        
        r_force = profile.force_R(0.1*a,0.2*b,1.2*a).in_(units.parsec/units.Myr**2)
        z_force = profile.force_z(0.1*a,0.2*b,1.2*a).in_(units.parsec/units.Myr**2)
        potential = profile.get_potential_at_point(0.|units.m,a*0.1,a*5.,b*0.5).in_(units.kpc**2/units.Myr**2)
        ax,ay,az = profile.get_gravity_at_point(0.|units.m,a*0.1,a*5.,b*0.5)
        density = profile.mass_density(0.1*a,-0.5*b,1.2*a).in_(units.MSun/units.kpc**3)
        vc = profile.circular_velocity_at_z0(1.0*a).in_(units.kms)
        m_r = profile.equivalent_enclosed_mass_in_plane(100.*a).in_(units.MSun)
        
        self.assertAlmostEqual(r_force, -0.263920797645|units.parsec/units.Myr**2, 12)
        self.assertAlmostEqual(z_force, -5.35763387945 |units.parsec/units.Myr**2, 11)
        self.assertAlmostEqual(potential, -0.017321166943|units.kpc**2/units.Myr**2, 12)
        self.assertAlmostEqual(ax, -0.0196231550925|units.parsec/units.Myr**2, 12)
        self.assertAlmostEqual(ay, -0.981157754625|units.parsec/units.Myr**2, 12)
        self.assertAlmostEqual(az, -0.107380572532 |units.parsec/units.Myr**2, 12)
        self.assertAlmostEqual(density, 1336672.32264|units.MSun/units.kpc**3, 5)
        self.assertAlmostEqual(vc, 149.569512197|units.kms, 9)
        self.assertAlmostEqual(m_r, 65.9857465656|1.e9*units.MSun, 9)
        
    def test2(self):
        mass = 1.984e4|units.MSun
        a = 0.0|units.parsec
        b = 6.66|units.parsec
        nm_profile = MiyamotoNagai_profile(mass,a,b)
        plummer_profile = Plummer_profile(mass,b)
        
        pot_nm = nm_profile.get_potential_at_point(0.|units.m,b*0.1,b*5.,b*0.2)
        pot_p = plummer_profile.get_potential_at_point(0.|units.m,b*0.1,b*5.,b*0.2)
        self.assertEqual(pot_nm,pot_p)
        
        ax_nm,ay_nm,az_nm = nm_profile.get_gravity_at_point(0.|units.m,b*0.1,b*5.,b*0.1)
        ax_p,ay_p,az_p = plummer_profile.get_gravity_at_point(0.|units.m,b*0.1,b*5.,b*0.1)
        print(ax_nm.in_(units.parsec/units.Myr**2), ax_p.in_(units.parsec/units.Myr**2))
        self.assertAlmostEqual(ax_nm.in_(units.parsec/units.Myr**2),ax_p.in_(units.parsec/units.Myr**2), 12)
        self.assertAlmostEqual(ay_nm.in_(units.parsec/units.Myr**2),ay_p.in_(units.parsec/units.Myr**2), 12)
        self.assertAlmostEqual(az_nm.in_(units.parsec/units.Myr**2),az_p.in_(units.parsec/units.Myr**2), 12)
        
        rho_nm = nm_profile.mass_density(b*0.,b*0.,b*6.6)
        rho_p = plummer_profile.mass_density(b*6.6)
        self.assertEqual(rho_nm,rho_p)

class TestPowerLawCutoff_profile(amusetest.TestCase):
    def test1(self):
        rho0 = 12.|units.MSun/units.parsec**3
        r0 = 1.6|units.parsec
        alpha = 1.6
        rc = 0.66|units.kpc
        power_law = PowerLawCutoff_profile(rho0,r0,alpha,rc)

        r_force = power_law.radial_force(0.1*r0).in_(units.parsec/units.Myr**2)
        potential = power_law.get_potential_at_point(0.|units.m,r0*0.1,r0*5.,r0*0.5).in_(units.kpc**2/units.Myr**2)
        ax,ay,az = power_law.get_gravity_at_point(0.|units.m,r0*0.1,r0*5.,r0*0.5)
        density = power_law.mass_density(6.6*r0).in_(units.MSun/units.parsec**3)
        vc = power_law.circular_velocity(r0).in_(units.kms)
        m_r = power_law.enclosed_mass(100.*r0).in_(units.MSun)
        
        self.assertAlmostEqual(r_force, -3.08704194743 |units.parsec/units.Myr**2, 10)
        self.assertAlmostEqual(potential, -5.91677122595e-06 |units.kpc**2/units.Myr**2, 10)
        self.assertAlmostEqual(ax, -0.00585557189412|units.parsec/units.Myr**2, 10)
        self.assertAlmostEqual(ay, -0.292778594706|units.parsec/units.Myr**2, 10)
        self.assertAlmostEqual(az, -0.0292778594706|units.parsec/units.Myr**2, 10)
        self.assertAlmostEqual(density, 0.585867989506 |units.MSun/units.parsec**3, 10)
        self.assertAlmostEqual(vc, 1.0891472277|units.kms, 10)
        self.assertAlmostEqual(m_r, 2.71756907682 |1.e5*units.MSun, 9)

    def test2(self):
        rho0 = 12.|units.MSun/units.parsec**3
        r0 = 1.6|units.parsec
        alpha = 1.6
        rc = 0.66|units.kpc
        power_law = PowerLawCutoff_profile(rho0,r0,alpha,rc)
        
        r=2*rc*numpy.arange(100001)/100000.
        
        pot=power_law.get_potential_at_point(0. | units.kpc, r, 0.*r,0*r)
        rho=power_law.mass_density( r)
        
        d=r[1]-r[0]
        rpot=r*pot
        lapl=(rpot[2:]-2*rpot[1:-1]+rpot[:-2])/d**2/r[1:-1]
        rcheck=r[1:-1]
        
        rho1=(-lapl/(4*units.pi*constants.G))
        rho2=rho[1:-1]
        d=((rho1[1:]-rho2[1:])/rho2[1:])
        # seems to be ok:        
        self.assertTrue(d.max()< 0.025)
        self.assertTrue(d.mean()< 1.e-6)
      
        
    def setUp(self):
        if not scipy_imported:
            self.skip("scipy not installed")

class TestMWpotentialBovy2015(amusetest.TestCase):
    def test1(self):
        """
        See Table 1 of Bovy 2015, http://adsabs.harvard.edu/abs/2015ApJS..216...29B
        """
        mw = MWpotentialBovy2015()
        r0 = 8.|units.kpc
        v0 = 220.|units.kms
        
        # total mass density at r=r0,z=0
        rho_r0_z0 = mw.mass_density(r0,0.|units.m,0.|units.m)
        self.assertAlmostEqual(rho_r0_z0, 0.10|units.MSun/units.parsec**3, 2)
        
        # halo mass density at r0
        rho_halo_at_r0 = mw.halo.mass_density(r0)
        self.assertAlmostEqual(rho_halo_at_r0, 0.008|units.MSun/units.parsec**3, 3)
        
        # mass enclosed in 60kpc
        mass_in_60 = mw.enclosed_mass(60.|units.kpc)
        print(mass_in_60.in_(units.MSun))
        self.assertAlmostEqual((mass_in_60/1.e11).in_(units.MSun), 4.08|units.MSun, 2)
        
        # normalization factor for bulge
        fr_r0_v0 = v0**2 / r0
        fr_r0_bulge = -mw.bulge.radial_force(r0)
        self.assertAlmostEqual(fr_r0_bulge/fr_r0_v0, 0.05, 3)
        
        # normalization factor for disk
        fr_r0_disk = -mw.disk.force_R(r0,0.|units.m,0.|units.m)
        self.assertAlmostEqual(fr_r0_disk/fr_r0_v0, 0.6, 3)
        
        # normalization factor for halo
        fr_r0_halo = -mw.halo.radial_force(r0)
        print(fr_r0_halo/fr_r0_v0)
        self.assertAlmostEqual(fr_r0_halo/fr_r0_v0, 0.35, 3)
        
    def setUp(self):
        if not scipy_imported:
            self.skip("scipy not installed")
