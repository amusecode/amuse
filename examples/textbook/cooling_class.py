import numpy

from amuse.units import units,constants
from amuse.units.quantities import zero

# Pelupessy et al. (in prep.) simple thermal model
class SimplifiedThermalModel(object):
  def __init__(self, n0=0.05 | units.cm**-3, T0=1.e4 | units.K,
                   Tmin=20 | units.K, alpha=5., reference_heating=1.e-25 | units.erg/units.s):
      self.reference_mu=(2.2 | units.amu)
      self.rho0=n0*self.reference_mu
      self.T0=T0
      self.Tmin=Tmin
      self.alpha=alpha
      self.reference_heating=reference_heating

  def equilibrium_temperature(self,rho):
    xclip=(rho/self.rho0)
    return self.Tmin+ (self.T0-self.Tmin)/(1.+numpy.log10(1.+9*xclip)**self.alpha)

  def mu(self,rho=None):
    if rho is None:
      return self.reference_mu
    else:
      return numpy.ones(numpy.shape(rho))*self.reference_mu

  def gamma(self,rho):
    return numpy.ones(numpy.shape(rho))*(self.reference_heating)

  def u_from_T(self,T):
    return constants.kB*T/self.mu()

  def T_from_u(self,u):
    return u/constants.kB*self.mu()

  def equilibrium_u(self,rho):
    return constants.kB*self.equilibrium_temperature(rho)/self.mu(rho)
    
  def tau(self,rho):
    return (constants.kB*self.equilibrium_temperature(rho)/self.gamma(rho))

  def evolve_u(self,dt,rho,u0,dudt=None):
    u_eq=self.equilibrium_u(rho)
    tau=self.tau(rho)
    if dudt is not None:
      condition1= 1.*(dudt*tau < (u0-u_eq))
      condition2= 1.-condition1
      fac=1./numpy.maximum(1-dudt/u0*tau ,1.e-5)
      u_eq=(u_eq*fac)*condition1+(u_eq+dudt*tau)*condition2
      tau=(tau*fac)*condition1+tau*condition2
    return u_eq+(u0-u_eq)*numpy.exp(-dt/tau)

  def evolve_u_radiated_energy(self,dt,rho,u0,dudt=None):
    u_eq0=self.equilibrium_u(rho)
    tau0=self.tau(rho)
    u_eq=u_eq0
    tau=tau0
    if dudt is not None:
      condition1= 1.*(dudt*tau < (u0-u_eq))
      condition2= 1.-condition1
      fac=1./numpy.maximum(1-dudt/u0*tau ,1.e-5)
      u_eq=(u_eq*fac)*condition1+(u_eq+dudt*tau)*condition2
      tau = (tau*fac)*condition1+tau*condition2
    u1=u_eq+(u0-u_eq)*numpy.exp(-dt/tau)
    rad=(u_eq-u_eq0)*dt/tau0+(u0-u_eq)*tau/tau0*(1-numpy.exp(-dt/tau))
    return u1,rad
# rad>0 -> cooling 
# rad<0 -> heating


class SimplifiedThermalModelEvolver(SimplifiedThermalModel):
        def __init__(self,particles,**kwargs):
                self.particles=particles
                SimplifiedThermalModel.__init__(self,**kwargs)
                self.radiated_energy=zero
                self.total_luminosity=zero
                self.model_time=zero
                self.umin=self.u_from_T(1. | units.K)

        def evolve_for(self, dt):
                #print " Do NOT Cool!"
                #return
                if dt>0*dt:
                    rho=self.particles.rho
                    u=self.particles.u
                    du_dt=self.particles.du_dt
                    new_u, lum= self.evolve_u_radiated_energy( dt, rho, u, du_dt )
                    self.radiated_energy+=(lum*self.particles.mass).sum()/dt
                    self.total_luminosity=(lum*self.particles.mass).sum()
                    a=numpy.where(new_u<self.umin)[0]
                    new_u[a]=self.umin
                    self.particles.u = new_u
                    # debug lines
                    nrho=numpy.isnan(rho.number).sum()
                    nu=numpy.isnan(u.number).sum()
                    ndu=numpy.isnan(du_dt.number).sum()
                    nnu=numpy.isnan(new_u.number).sum()
                    if nrho+nu+ndu+nnu>0:
                        print("nan detected in thermal evolution")
                        print(nrho,nu,ndu,nnu)
                        import pickle
                        with open("cooling_dump","w") as f:
                            pickle.dump((dt,rho,u,du_dt,new_u),f)
                        raise Exception("NaNs in thermal evolution")

        def evolve_model(self,tend):
                self.evolve_for(tend-self.model_time)
                self.model_time=tend

# COOLING
class Cooling(object):
        def __init__(self, particles): 
                self.particles = particles
                self.umin=self.u_from_T( 10. | units.K)
                self.umax=self.u_from_T( 1.e6 | units.K)

        def evolve_for(self, dt):
              #print " Do NOT Cool!"
              #return
              if dt>0*dt:
                  new_u=self.evolve_internal_energy(self.particles.u, dt,self.particles.rho/self.mu(),self.particles.du_dt)
                  a=numpy.where(new_u<self.umin)[0]
                  new_u[a]=self.umin
                  a=numpy.where(new_u>self.umax)[0]
                  new_u[a]=self.umax
                  self.particles.u = new_u
        
        def evolve_internal_energy(self, u_0, dt, n_H, du_dt_adiabatic = zero):
                function = lambda u: ((self.gerritsen_heating_function() - n_H * self.my_cooling_function(self.T_from_u(u))) / self.mu())#du_dt_adiabatic * u/u_0 + 

                u_out = self.integrate_ode(function, u_0, dt)
                return u_out
        
        def integrate_ode(self,function, x, t_end, eps = 0.01):
                """
                Integrates the given ordinary differential equation of the form:
                dx/dt = function(x)
                for a time 't_end', using the initial value 'x'.
                The routine takes small steps, such that (abs(dx) <= eps * x)
                """
                t = 0 * t_end
                while t < t_end:
                        fx = function(x)
                        dtinv=(abs(fx)/(eps*x)).amax()
                        step=t_end-t
                        if dtinv!=0*dtinv : 
                                step = min( step, 1./dtinv )
                        t += step
                        x += fx * step
                return x

        # Transforming from T to U
        def u_from_T(self, T):
                return 3.0/2.0 * constants.kB * T / self.mu()

        # Transforming from U to T
        def T_from_u(self, u):
                return 2.0/3.0 * u * self.mu() / constants.kB

        # Molecular weight
        def mu(self, X = None, Y = 0.25, Z = 0.02, x_ion = 0.1):
                """
                Compute the mean molecular weight in kg (the average weight of particles in a gas)
                X, Y, and Z are the mass fractions of Hydrogen, of Helium, and of metals, respectively.
                x_ion is the ionisation fraction (0 < x_ion < 1), 1 means fully ionised
                """
                if X is None:
                        X = 1.0 - Y - Z
                elif abs(X + Y + Z - 1.0) > 1e-6:
                        raise Exception("Error in calculating mu: mass fractions do not sum to 1.0")
                return constants.proton_mass / (X*(1.0+x_ion) + Y*(1.0+2.0*x_ion)/4.0 + Z*x_ion/2.0)
        
        # G depends on nearby sources, see 1997A&A...325..972G
        def gerritsen_heating_function(self, G_0 = 10, eps = 0.05):
                return 10.0**-24 * eps * G_0 | units.erg / units.s

        def gerritsen_cooling_function(self, T, logT = None, a = 3.24, b = 0.170): # x=1e-1
                if logT is None:
                        logT = numpy.log10(T.value_in(units.K))
                condlist = [logT <= 6.2, logT >= 6.2]
                choicelist = [10.0**-21.0 * (10**(-0.1-1.88*(5.23-logT)**4) + 10**(-a-b*(4-logT)**2)), 10.0**-22.7]
                return (units.erg*units.cm**3/units.s).new_quantity(numpy.select(condlist, choicelist))

        def my_cooling_function(self, T, logT = None, a = 3.24, b = 0.170): # x=1e-1
                if logT is None:
                        logT = numpy.log10(T.value_in(units.K))
                condlist = [logT <= 6.2, logT >= 6.2]
                choicelist = [10.0**-21.0 * (10**(-0.1-1.88*(5.23-logT)**4) + 10**(-a-b*abs(4-logT)**3)), 10.0**-22.7]
                return (units.erg*units.cm**3/units.s).new_quantity(numpy.select(condlist, choicelist))
