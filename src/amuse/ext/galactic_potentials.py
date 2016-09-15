import numpy

from amuse.units import constants, units
from amuse.support.literature import LiteratureReferencesMixIn

try:
    from scipy.special import gammainc,gamma
    scipy_imported = True
except:
    scipy_imported = False

class NFW_profile(LiteratureReferencesMixIn):
    """
    Gravitational potential of the NFW (1996) halo
    Two-power density spherical model suitable for modeling dark matter halos.
    Density--potential pair:
    * density(r) = rho0 / [r/rs * (1+r/rs)^2], where is the spherical radius
    * potential(r) = -4*pi*G*rho0*rs^2 * ln(1+r/rs)/(r/rs)
    
    .. [#] Navarro, Julio F.; Frenk, Carlos S.; White, Simon D. M., The Astrophysical Journal, Volume 490, Issue 2, pp. 493-508 (1996)
    
    :argument rho0: density parameter
    :argument rs: scale radius
    """
    def __init__(self,rho0,rs,G=constants.G):
        LiteratureReferencesMixIn.__init__(self)
        self.rho0 = rho0
        self.rs = rs
        self.G = G
        self.four_pi_rho0 = 4.*numpy.pi*self.rho0
        self.four_pi_rho0_G = self.four_pi_rho0*self.G
    
    def radial_force(self,r):
        r_rs = r/self.rs
        ar = self.four_pi_rho0_G*self.rs**3*(1./(r*self.rs+r**2)-(1./r**2)*numpy.log(1.+r_rs))
        #ar = self.four_pi_rho0_G*self.rs*((r_rs-(1.+r_rs)*numpy.log(1.+r_rs))/r_rs**2/(1.+r_rs))
        return ar
    
    def get_potential_at_point(self,eps,x,y,z):
        r = (x**2+y**2+z**2).sqrt()
        r_rs = r/self.rs
        return -1.*self.four_pi_rho0_G*self.rs**2*numpy.log(1.+r_rs)/r_rs
    
    def get_gravity_at_point(self,eps,x,y,z):
        r = (x**2+y**2+z**2).sqrt()
        fr = self.radial_force(r)
        ax = fr*x/r
        ay = fr*y/r
        az = fr*z/r
        return ax,ay,az
    
    def enclosed_mass(self,r):
        fr = self.radial_force(r)
        return -r**2/self.G*fr
    
    def circular_velocity(self,r):
        fr = self.radial_force(r)
        return (-r*fr).sqrt()
    
    def mass_density(self,r):
        r_rs = r/self.rs
        return self.rho0 / (r_rs*(1.+r_rs)**2)

class MiyamotoNagai_profile(LiteratureReferencesMixIn):
    """
    Miyamoto and Nagai (1975) axisymmetric disk.
    * potential(R,z) = -GM / sqrt(R**2 + (a+sqrt(z**2+b**2))**2)
    
    .. [#] Miyamoto, M.; Nagai, R., Astronomical Society of Japan, Publications, vol. 27, no. 4, 1975, p. 533-543 (1975)
    
    :argument mass: total mass
    :argument a: disk scale radius
    :argument b: disk scale height
    """
    def __init__(self,mass,a,b,G=constants.G):
        LiteratureReferencesMixIn.__init__(self)
        self.mass = mass
        self.a = a
        self.b = b
        self.G = G
        self.GM = self.G*self.mass
        self.a2 = self.a**2
        self.b2 = self.b**2
        
    def force_R(self,x,y,z):
        R2 = x**2+y**2
        R = R2.sqrt()
        sqrt_z2_b2 = (z**2+self.b2).sqrt()
        return -self.GM*R*(R2+(self.a+sqrt_z2_b2)**2)**(-1.5)
    
    def force_z(self,x,y,z):
        R2 = x**2+y**2
        sqrt_z2_b2 = (z**2+self.b2).sqrt()
        a_sqrt_z2_b2 = self.a+sqrt_z2_b2
        return -self.GM*z*a_sqrt_z2_b2/((R2+a_sqrt_z2_b2**2)**1.5*sqrt_z2_b2)
        
    def get_potential_at_point(self,eps,x,y,z):
        R2 = x**2+y**2
        return -self.GM/(R2+(self.a+(self.b2+z**2).sqrt())**2).sqrt()
    
    def get_gravity_at_point(self,eps,x,y,z):
        fR = self.force_R(x,y,z)
        R = (x**2+y**2).sqrt()
        ax = fR*x/R
        ay = fR*y/R
        az = self.force_z(x,y,z)
        return ax,ay,az
    
    def mass_density(self,x,y,z):
        R2 = x**2+y**2
        z2_b2 = z**2+self.b2
        sqrt_z2_b2 = z2_b2.sqrt()
        rho = self.b2*self.mass/(4.*numpy.pi) * \
            (self.a*R2+(self.a+3.*sqrt_z2_b2)*(self.a+sqrt_z2_b2)**2) / \
                ((R2+(self.a+sqrt_z2_b2)**2)**2.5*z2_b2**1.5)
        return rho
    
    def circular_velocity_at_z0(self,R):
        fR_at_z0 = self.force_R(R,0.|units.kpc,0.|units.kpc)
        return (-R*fR_at_z0).sqrt()
    
    def equivalent_enclosed_mass_in_plane(self,R):
        """
        mass, that would be enclosed in profile corresponding the disk profile in the 
        galactic plane (z=0)
        """
        fR_at_z0 = self.force_R(R,0.|units.kpc,0.|units.kpc)
        return -R**2/self.G*fR_at_z0
    
class Plummer_profile(LiteratureReferencesMixIn):
    """
    Spherically symmetric Plummer (1911) profile
    * potential(r) = -GM / sqrt(a**2 + r**2)
    * density(r) = (3M/4pi*a**3) * (1+(r/a)**2)**(-5/2)
    
    .. [#] Plummer, H. C., MNRAS, Vol. 71, p.460-470 (1911)
    
    :argument mass: total mass
    :argument a: scale radius
    """
    def __init__(self,mass,a,G=constants.G):
        LiteratureReferencesMixIn.__init__(self)
        self.mass = mass
        self.a = a
        self.G = G
        self.GM = self.G*self.mass
        self.a2 = self.a**2
    
    def radial_force(self,r):
        r2 = r**2
        return -self.GM*r*(r2+self.a2)**(-1.5)
    
    def get_gravity_at_point(self,eps,x,y,z):
        r = (x**2+y**2+z**2).sqrt()
        fr = self.radial_force(r)
        ax = fr*x/r
        ay = fr*y/r
        az = fr*z/r
        return ax, ay, az
    
    def get_potential_at_point(self,eps,x,y,z):
        r2 = x**2+y**2+z**2
        return -self.GM/(r2+self.a2).sqrt()
    
    def mass_density(self,r):
        return self.mass/(4./3.*numpy.pi*self.a**3)*(1.+(r/self.a)**2)**(-2.5)
    
    def enclosed_mass(self,r):
        fr = self.radial_force(r)
        return -r**2/self.G*fr
    
    def circular_velocity(self,r):
        fr = self.radial_force(r)
        return (-r*fr).sqrt()
    
        
class PowerLawCutoff_profile(LiteratureReferencesMixIn):
    """
    Spherically symmetric power-law density with exponential cut-off
    * density(r) = rho0*(r0/r)^alpha*exp(-(r/rc)^2)
    
    :argument rho0: density amplitude
    :argument r0: power-law scaling radius
    :argument alpha: power-law index, alpha<3
    :argument rc: cut-off radius
    """
    def __init__(self,rho0,r0,alpha,rc,G=constants.G):
        LiteratureReferencesMixIn.__init__(self)
        self.rho0 = rho0
        self.r0 = r0
        self.alpha = alpha
        self.rc = rc
        self.G = G
        self.rho0_r0_to_alpha = self.rho0*self.r0**self.alpha
        
        if 3.<=self.alpha: print "Warning: power-law index must be less than 3."
    
    def get_potential_at_point(self,eps,x,y,z):
        if scipy_imported == False:
            AmuseWarning("importing scipy failed, maybe not installed")
        r = (x**2+y**2+z**2).sqrt()
        r_rc = r/self.rc
        return -2.*numpy.pi*self.G*self.rho0_r0_to_alpha*self.rc**(3.-self.alpha)/r* \
            (gamma(1.5-0.5*self.alpha)*gammainc(1.5-0.5*self.alpha,r_rc**2)+r_rc*gamma(1.-0.5*self.alpha)*(1.-gammainc(1.-0.5*self.alpha,r_rc**2)))
    
    def radial_force(self,r):
        Mr = self.enclosed_mass(r)
        return -self.G*Mr/r**2
    
    def get_gravity_at_point(self,eps,x,y,z):
        r = (x**2+y**2+z**2).sqrt()
        fr = self.radial_force(r)
        ax = fr*x/r
        ay = fr*y/r
        az = fr*z/r
        return ax, ay, az
    
    def mass_density(self,r):
        return self.rho0_r0_to_alpha*r**(-self.alpha)*numpy.exp(-(r/self.rc)**2.)
    
    def enclosed_mass(self,r):
        """
        careful with scipy.special.gammainc :
        gammainc(a,x) = 1 / gamma(a) * integral(exp(-t) * t**(a-1), t=0..x)
        """
        if scipy_imported == False:
            AmuseWarning("importing scipy failed, maybe not installed")
        return 2.*numpy.pi*self.rho0_r0_to_alpha*self.rc**(3.-self.alpha)* \
            gamma(1.5-0.5*self.alpha)*gammainc(1.5-0.5*self.alpha,(r/self.rc)**2)
    
    def circular_velocity(self,r):
        fr = self.radial_force(r)
        return (-r*fr).sqrt()

class MWpotentialBovy2015(LiteratureReferencesMixIn):
    """
    MW-like galaxy potential consists of a bulge modeled as a power-law density
    profile that is exponentially cut-off, a Miyamoto & Nagai disk; and a
    NFW dark-matter halo. Parameters of individual components are based on
    fits to observational data. In addition to these constraints, the solar
    distance to the Galactic center is set to R0=8kpc and the circular velocity
    at the Sun to V0=220km/s.
    
    .. [#] Bovy, J; ApJSS, Volume 216, Issue 2, article id. 29, 27 pp. (2015)
    
    """
    def __init__(self):
        LiteratureReferencesMixIn.__init__(self)
        self.bulge = PowerLawCutoff_profile(2.22638e8|units.MSun/units.kpc**3, 1.|units.kpc, 1.8, 1.9|units.kpc)
        self.disk = MiyamotoNagai_profile(6.81766163214e10|units.MSun, 3.|units.kpc, 0.28|units.kpc)
        self.halo = NFW_profile(8484685.92946|units.MSun/units.kpc**3, 16.|units.kpc)
    
    def get_potential_at_point(self,eps,x,y,z):
        return self.bulge.get_potential_at_point(eps,x,y,z) + \
            self.disk.get_potential_at_point(eps,x,y,z) + \
            self.halo.get_potential_at_point(eps,x,y,z)
    
    def get_gravity_at_point(self,eps,x,y,z):
        ax_b,ay_b,az_b = self.bulge.get_gravity_at_point(eps,x,y,z)
        ax_d,ay_d,az_d = self.disk.get_gravity_at_point(eps,x,y,z)
        ax_h,ay_h,az_h = self.halo.get_gravity_at_point(eps,x,y,z)
        return ax_b+ax_d+ax_h, ay_b+ay_d+ay_h, az_b+az_d+az_h
    
    def mass_density(self,x,y,z):
        r = (x**2+y**2+z**2).sqrt()
        return self.bulge.mass_density(r)+self.disk.mass_density(x,y,z)+self.halo.mass_density(r)
    
    def circular_velocity(self,r):
        return (self.bulge.circular_velocity(r)**2+self.disk.circular_velocity_at_z0(r)**2+self.halo.circular_velocity(r)**2).sqrt()
    
    def enclosed_mass(self,r):
        return self.bulge.enclosed_mass(r)+self.disk.equivalent_enclosed_mass_in_plane(r)+self.halo.enclosed_mass(r)
    
