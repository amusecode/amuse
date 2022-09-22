# -*- coding: utf-8 -*-

import numpy
import scipy
from scipy.interpolate import BPoly, interp1d, UnivariateSpline
from numpy import exp, sqrt, pi, sin, cos, log10
from scipy.special import gamma, gammainc, hyp1f1
from scipy.integrate import ode, simps, quad
from math import factorial, sinh

#     Authors: Mark Gieles, Alice Zocchi (Surrey 2015)

class limepy:
    def __init__(self, phi0, g, **kwargs):
        r"""
        (MM, A) LIMEPY
        (Multi-Mass, Anisotropic) Lowered Isothermal Model Explorer in Python
        This code solves the models presented in Gieles & Zocchi 2015 (GZ15),
        and calculates radial profiles for some useful quantities. The models
        are defined by the distribution function (DF) of equation (1) in GZ15.
        Model parameters:
        =================
        phi0 : scalar, required
             Central dimensionless potential
        g : scalar, required
          Order of truncation (0<= g < 3.5; 0=Woolley, 1=King, 2=Wilson)
        ra : scalar, required for anisotropic models
           Anisotropy radius; default=1e8
        mj : list, required for multi-mass system
           Mean mass of each component; default=None
        Mj : list, required for multi-mass system
           Total mass of each component; default=None
        delta : scalar, optional
              Index in s_j = s x mu_j^-delta; default=0.5
              See equation (24) in GZ15
        eta : scalar, optional
            Index in ra_j = ra x mu_j^eta; default=0
            See equation (25) in GZ15
        Input for scaling:
        ==================
        G : scalar, optional
          Final scaled gravitationsl const; default=0.004302 [(km/s)^2 pc/Msun]
        M : scalar, optional
          Final scaled mass; default=10^5
        r0, rh, rv, rt : scalar, optional
                       Final scaled radius; default=rh=3
        Options:
        ========
        project : bool, optional
                Compute model properties in projection; default=False
        meanmassdef : string [global|central]
                    Definition of <m> in mu_j = m_j/<m>; default='global'
        potonly : bool, optional
                Fast solution by solving potential only; default=False
        max_step : scalar, optional
                 Maximum step size for ode output; default=1e10
        verbose : bool, optional
                Print diagnostics; default=False
        ode_atol : absolute tolerance parameter for ode solver; default=1e-7
        ode_rtol : relative tolerance parameter for ode solver; default=1e-7
        Output variables:
        =================
        All models:
        -----------
         rhat, phihat, rhohat : radius, potential and density in model units
         r, phi, rho : as above, in scaled units
         v2, v2r, v2t : total, radial and tangential mean-square velocity
         beta : anisotropy profile (equation 32, GZ15)
         mc : enclosed mass profile
         r0, rh, rv, rt : radii (King, half-mass, virial, truncation)
         K, Kr, Kt : kinetic energy: total, radial, tangential
         U, Q : potential energy, virial ratio
         A : constant in DF (equation 1, GZ15)
         volume : phase-space volume occupied by model
         nstep : number of integration steps (depends on ode_rtol & ode_atol)
         converged : bool flag to indicate whether model was solved
        Projected models:
        -----------------
         Sigma : surface (mass) density
         v2p : line-of-sight mean-square velocity
         v2R, v2T : radial and tangential component of mean-square velocity
                  on plane of the sky
        Multi-mass models:
        ------------------
         Properties of each component:
         phi0j : dimensionless central potential
         rhohatj : dimensionless density
         rhoj : density profile
         v2j : mean-square velocity profile
         v2rj, v2tj : radial and tangential component of mean-square velocity
                    profile
         r0j, raj : radii (King, anisotropy)
         Kj : kinetic energy
         Krj, Ktj : radial/tangential component of kinetic energy
        Projected multi-mass models:
        ---------------------------
         Properties of each component:
         Sigmaj : surface (mass) density
         v2pj : line-of-sight mean-square velocity profile
         v2Rj, v2Tj : radial and tangential component on the plane of the sky
                    of the mean-square velocity profile
        Examples:
        =========
        Construct a Woolley model with phi0 = 7 and print r_t/r_0 and r_v/r_h
        >>> k = limepy(7, 0)
        >>> print k.rt/k.r0, k.rv/k.rh
        >>> 19.1293426415 1.17783663028
        Construct a Michie-King model and print ratio of anisotropy radius over
        half-mass radius and the Polyachenko & Shukhman (1981) anisotropy
        parameter
        >>> a = limepy(7, 1, ra=5)
        >>> print a.ra/a.rh, 2*a.Kr/a.Kt
        >>> 1.03377960149 1.36280949941
        Create a Wilson model with phi0 = 12 in Henon/N-body units: G = M =
        r_v = 1 and print the normalisation constant A of the DF and the
        value of the DF in the centre:
        >>> w = limepy(12, 2, G=1, M=1, rv=1)
        >>> print w.A, w.df(0,0)
        >>> [ 0.00800902] [ 1303.40270676]
        Multi-mass model in physical units with r_h = 3 pc and M = 10^5 M_sun
        and print central densities of each bin over the total central density
        and the half-mass radius + half-mass radius in projection
        >>> m = limepy(7, 1, mj=[0.3,1,5], Mj=[9,3,1], rh=3, M=1e5,project=True,\
        ...: meanmassdef='central')
        >>> print m.alpha, m.rh, m.rhp
        >>> [ 0.30721416  0.14103549  0.55175035] 3.0 2.25494426759
        """

        # Set parameters and scales
        self._set_kwargs(phi0, g, **kwargs)
        self.rhoint0 = [self._rhoint(self.phi0, 0, self.ramax)]


        # In case of multi-mass model, iterate to find central densities
        # (see Section 2.2 in GZ15)

        if (self.multi):
            # Feb 2016: New method to iterate on the MF.
            # First solve isotropic model

            ratmp = self.ra*1.0

            self.ra = self.ramax*1.0
            self._init_multi(self.mj, self.Mj)

            while self.diff > self.diffcrit:
                self._poisson(True)
                if (not self.converged):
                    error = "Error: rmax reached in mf iteration"
                    raise ValueError(error)
                else:
                    self._set_alpha()
                    if self.niter > self.max_mf_iter:
                        self.converged=False
                        error = "Error: maximum number of iterations reached "
                        raise ValueError(error)

            # Now solve anisotropic multimass
            if ratmp < self.ramax:
                self.diff = 1
                self.ra = ratmp*1.0
                self._set_mass_function_variables()

                while self.diff > self.diffcrit:
                    self._poisson(True)
                    if (not self.converged):
                        error = "Error: rmax reached in mf iteration"
                        raise ValueError(error)
                    else:
                        self._set_alpha()
                        if self.niter > self.max_mf_iter:
                            self.converged=False
                            error ="Error: maximum number of iterations reached"
                            raise ValueError(error)


        self.r0 = 1.0
        if (self.multi): self.r0j = sqrt(self.s2j)*self.r0

        # ROT
        if (self.rot):
            self._poisson_rot()

        else:
            # Solve Poisson equation to get the potential
            self._poisson(self.potonly)

        if (self.multi): self.Mj = self._Mjtot

        # Optional scaling
        if (self.scale): self._scale()

        # Optional computation of model properties in projection
        if (self.project): self._project()

        # Optional output
        if (self.verbose):
            print("\n Model properties: ")
            print(" ----------------- ")
            print(" phi0 = %5.2f; g = %4.2f"%(self.phi0, self.g))
            print(" Converged = %s"%(self.converged))
            if (self.potonly):
                print(" M = %10.3f; U = %10.4f "%(self.M, self.U))
            else:
                out1 = (self.M,self.U,self.K,-self.K/self.U,2*self.Kr/self.Kt)
                frm = " M = %10.3e; U = %9.3e; K = %9.3e; Q = %6.4f; "
                frm += " 2Kr/Kt = %5.3f"
                print(frm%out1)

            out2 = (self.rv/self.rh, self.rh/self.r0)
            out2 += (self.rt/self.r0, self.ra/self.rh)
            frm = " rv/rh = %4.3f; rh/r0 = %6.3f; "
            frm += "rt/r0 = %7.3f; ra/rh = %7.3f"
            print(frm%out2)

    def _set_kwargs(self, phi0, g, **kwargs):
        """ Set parameters and scales """

        if (g<0): raise ValueError("Error: g must be larger or equal to 0")
        if (g>=3.5): raise ValueError("Error: for g>=3.5 models are infinite")

        self.model = "limepy"

        # ROT
        self.omega = 0

        self.phi0, self.g = phi0, g
        self._MS, self._RS, self._GS = None, None, None
        self.scale_radius = None
        self.scale = False
        self.project = False
        self.meanmassdef='global'
        self.maxr = 1e10
        self.max_step = self.maxr
        self.diffcrit = 1e-8
        self.max_arg_exp = 700  # Maximum argument for exponent and hyp1f1 func
        self.max_mf_iter = 100  # Maximum number of iterations to find rho0j
        self.minimum_phi = 1e-8 # Stop criterion for integrator
        self.mf_iter_index = 0.5
        self.ode_atol = 1e-7
        self.ode_rtol = 1e-7
        self.nmbin, self.delta, self.eta = 1, 0.5, 0.0

        self.G = 9.0/(4.0*pi)
        self.mu, self.alpha = numpy.array([1.0]), numpy.array([1.0])
        self.s2 = 1.0
        self.s2j = numpy.array([1.0])
        self.niter = 0

        self.potonly, self.multi, self.verbose = [False]*3
        self.ra, self.ramax = 1e8, 1e8

        self.nstep=1
        self.converged=True
        self._interpolator_set=False

        if kwargs is not None:
            for key, value in kwargs.items():
                # Check for scaling input
                if key == 'G':
                    self._GS, self.scale = value, True
                elif key == 'M':
                    self._MS, self.scale = value, True
                elif key == 'r0':
                    if self.scale_radius is None:
                        self._RS,self.scale_radius,self.scale = value,'r0', True
                    else:
                        error="Can not set scale radius to r0,already set to %s"
                        raise ValueError(error%self.scale_radius)
                elif key == 'rh':
                    if self.scale_radius is None:
                        self._RS, self.scale_radius, self.scale=value,'rh', True
                    else:
                        error="Can not set scale radius to rh,already set to %s"
                        raise ValueError(error%self.scale_radius)
                elif key == 'rv':
                    if self.scale_radius is None:
                        self._RS, self.scale_radius, self.scale=value,'rv', True
                    else:
                        error="Can not set scale radius to rv,already set to %s"
                        raise ValueError(error%self.scale_radius)
                elif key == 'rt':
                    if self.scale_radius is None:
                        self._RS, self.scale_radius, self.scale=value,'rt', True
                    else:
                        error="Can not set scale radius to rt,already set to %s"
                        raise ValueError(error%self.scale_radius)
                else:
                    # Set input parameters
                    setattr(self, key, value)

            if (self.scale):
                if self._MS is None:
                    self._MS = 1e5
                    if (self.verbose):
                        print(" No mass-scale provided, set to default M = 1e5")
                if self._RS is None:
                    self._RS, self.scale_radius = 3, 'rh'
                    if (self.verbose):
                        print(" No radius-scale provided, set to default rh = 3")
                if self._GS is None:
                    self._GS = 0.004302
                    if (self.verbose):
                        print(" No G provided, set to default: G = 0.004302")
                if (self.verbose):
                    vars=(self._GS, self._MS, self.scale_radius, self._RS)
                    print(" Model scaled to: G = %s, M = %s, %s = %s"%vars)

            if 'mj' in kwargs and 'Mj' in kwargs:
                self.multi=True
                if len(self.Mj) is not len(self.mj):
                    raise ValueError("Error: Mj and mj must have same length")
            if ('mj' not in kwargs and 'Mj' in kwargs) or \
               ('Mj' not in kwargs and 'mj' in kwargs):
                raise ValueError("Error: Supply both mj and Mj")
        self.raj = numpy.array([self.ra])

        if self.potonly and self.project:
            raise ValueError('Error: You must not use potonly and project')

        self.rot = False
        if self.omega > 0:
            print(" Warning: ROTATION PART NOT FINISHED! ")
            self.rot = True
        return

    def _logcheck(self, t, y):
        """ Logs steps and checks for final values """

        if (t>0): self.r, self._y = numpy.r_[self.r, t], numpy.c_[self._y, y]
        self.nstep+=1
        return 0 if (y[0]>self.minimum_phi) else -1

    def _set_mass_function_variables(self):
        """ Multi-mass models: Set properties for each mass bin """

        if self.meanmassdef=='global':
            Nj = self.Mj/self.mj
            self.mmean = sum(self.Mj/sum(Nj))
        elif self.meanmassdef=='central':
            self.mmean = sum(self.mj*self.alpha)    # equation (26) GZ15
        else:
            raise ValueError(" meanmass must be 'global' or 'central'")

        self.mu = self.mj/self.mmean
        self.s2j = self.mu**(-2*self.delta)         # equation (24) GZ15
        self.raj = self.ra*self.mu**self.eta        # equation (25) GZ15

        self.phi0j = self.phi0/self.s2j
        self.rhoint0 = numpy.zeros(self.nmbin)

        for j in range(self.nmbin):
            self.rhoint0[j] = self._rhoint(self.phi0j[j], 0, self.ramax)


    def _init_multi(self, mj, Mj):
        """
        Initialise parameters and arrays for multi-mass system
        (Section 2.2 in GZ15)
        """

        self.multi=True
        self.mj = numpy.array(mj)
        self.Mj = numpy.array(Mj)
        self.nmbin = len(mj)

        # Set trial value for alpha_j array, will be updated in iterations
        self.alpha = self.Mj/sum(self.Mj)
        self._set_mass_function_variables()
        self.diff = 1

    def _set_alpha(self):
        """ Set central rho_j for next iteration """

        # The power of mf_iter_index = 0.5. This is lower than the recommended
        # value of 1 as used in Da Costa & Freeman 1976 and Gunn & Griffin 1979.
        # Better convergence is reached for a smaller value, but a user may
        # experiment with larger values to gain speed.

        self.alpha *= (self.Mj/self._Mjtot)**self.mf_iter_index
        self.alpha/=sum(self.alpha)

        self._set_mass_function_variables()
        self.diff = sum((self._Mjtot/sum(self._Mjtot) -
                         self.Mj/sum(self.Mj))**2)/len(self._Mjtot)
        self.niter+=1
        self.nstep=1
        if (self.verbose):
            fracd,Mjit="", ""
            for j in range(self.nmbin):
                M1 = self._Mjtot[j]/sum(self._Mjtot)
                M2 = self.Mj[j]/sum(self.Mj)
                fracd=fracd+"%7.3f "%(( M1 - M2)/M2)
                Mjit=Mjit+"%7.3f "%(self._Mjtot[j]/sum(self._Mjtot))
            out = (self.niter, self.mf_iter_index, self.diff, self.converged*1,\
                   fracd, Mjit)
            frm = " %2i; ind=%3.1f; diff= %6.1e; conv= %s;"
            frm += " fdiff=%s; Mjtot=%s"
            print(frm%out)

    def _poisson(self, potonly):
        """ Solves Poisson equation """
        # y = [phi, u_j, U, K_j], where u = -M(<r)/G

        # Initialize
        self.r = numpy.array([0])
        self._y = numpy.r_[self.phi0, numpy.zeros(self.nmbin+1)]
        if (not potonly): self._y = numpy.r_[self._y,numpy.zeros(2*self.nmbin)]
        self._y = numpy.r_[self._y, 0]

        # Ode solving using Runge-Kutta integator of order 4(5) 'dopri5'
        # (Hairor, Norsett& Wanner 1993)
        max_step = self.maxr if (potonly) else self.max_step
        sol = ode(self._odes)
        sol.set_integrator('dopri5',nsteps=1e6,max_step=max_step,
                           atol=self.ode_atol,rtol=self.ode_rtol)
        sol.set_solout(self._logcheck)
        sol.set_f_params(potonly)
        sol.set_initial_value(self._y,0)
        sol.integrate(self.maxr)

        # Extrapolate to r_t:
        # phi(r) =~ a(r_t -r)
        # a = GM/r_t^2
        GM = -self.G*sum(sol.y[1:1+self.nmbin])
        p = 2*sol.y[0]*self.r[-1]/GM

        if (p<=0.5):
            rtfac = (1 - sqrt(1-2*p))/p
            self.rt = rtfac*self.r[-1] if (rtfac > 1) else 1.0000001*self.r[-1]
        else:
            self.rt = 1.000001*self.r[-1]

        # Set the converged flag to True if successful
        if (self.rt < self.maxr)&(sol.successful()):
            self.converged=True
        else:
            self.converged=False

        # Calculate the phase space volume occupied by the model
        dvol = (4./3*pi)**2*(self.rt**3-self.r[-1]**3)*0.5
        dvol *= (2*self._y[0,-1])**1.5
        self.volume = self._y[-1][-1]+dvol

        # Fill arrays needed if potonly=True
        self.r = numpy.r_[self.r, self.rt]
        self.rhat = self.r*1.0

        self.phihat = numpy.r_[self._y[0,:], 0]
        self.phi = self.phihat*1.0

        self._Mjtot = -sol.y[1:1+self.nmbin]/self.G

        self.M = sum(self._Mjtot)

        # Save the derivative of the potential for the potential interpolater
        dphidr = numpy.sum(self._y[1:1+self.nmbin,1:],axis=0)/self.r[1:-1]**2
        self.dphidrhat1 = numpy.r_[0, dphidr, -self.G*self.M/self.rt**2]

        self.A = self.alpha/(2*pi*self.s2j)**1.5/self.rhoint0

        if (not self.multi):
            self.mc = -numpy.r_[self._y[1,:], self._y[1,-1]]/self.G

        if (self.multi):
            self.mc = sum(-self._y[1:1+self.nmbin,:]/self.G)
            self.mc = numpy.r_[self.mc, self.mc[-1]]

        # Compute radii to be able to scale in case potonly=True
        self.U = self._y[1+self.nmbin,-1]  - 0.5*self.G*self.M**2/self.rt

        # Get half-mass radius from cubic interpolation, because the half-mass
        # radius can be used as a scale length, this needs to be done
        # accurately, a linear interpolation does not suffice. Because the
        # density array is not known yet at this point, we need a temporary
        # evaluation of density in the vicinity of r_h
        ih = numpy.searchsorted(self.mc, 0.5*self.mc[-1])-1
        rhotmp=numpy.zeros(2)
        for j in range(self.nmbin):
            phi = self.phihat[ih:ih+2]/self.s2j[j]
            rhotmp += self.alpha[j]*self._rhohat(phi, self.r[ih:ih+2], j)
        drdm = 1./(4*pi*self.r[ih:ih+2]**2*rhotmp)
        rmc_and_derivs = numpy.vstack([[self.r[ih:ih+2]],[drdm]]).T

        self.rh = BPoly.from_derivatives(self.mc[ih:ih+2], rmc_and_derivs)(0.5*self.mc[-1])
        self.rv = -0.5*self.G*self.M**2/self.U

        # Additional stuff
        if (not potonly):
            # Calculate kinetic energy (total, radial, tangential)
            self.K = numpy.sum(sol.y[2+self.nmbin:2+2*self.nmbin])
            self.Kr = numpy.sum(sol.y[2+2*self.nmbin:2+3*self.nmbin])
            self.Kt = self.K - self.Kr

            # Calculate density and velocity dispersion components
            if (not self.multi):
                self.rhohat = self._rhohat(self.phihat, self.r, 0)
                self.rho = self.rhohat*1.0
                self.v2, self.v2r, self.v2t = \
                        self._get_v2(self.phihat, self.r, self.rhohat, 0)

            # For multi-mass models, calculate quantities for each mass bin
            if (self.multi):
                for j in range(self.nmbin):
                    phi = self.phihat/self.s2j[j]
                    rhohatj = self._rhohat(phi, self.r, j)
                    v2j, v2rj, v2tj = self._get_v2(phi, self.r, rhohatj, j)
                    v2j, v2rj, v2tj = (q*self.s2j[j] for q in [v2j,v2rj,v2tj])
                    betaj = self._beta(self.r, v2rj, v2tj)

                    kj = self._y[2+self.nmbin+j,:]
                    krj = self._y[2+2*self.nmbin+j,:]
                    ktj = kj - krj

                    mcj = -numpy.r_[self._y[1+j,:], self._y[1+j,-1]]/self.G
                    rhj = numpy.interp(0.5*mcj[-1], mcj, self.r)

                    if (j==0):
                        self.rhohatj = rhohatj
                        self.rhohat = self.alpha[0] * self.rhohatj
                        self.v2j, self.v2rj, self.v2tj = v2j, v2rj, v2tj
                        self.v2 = self._Mjtot[j]*v2j/self.M
                        self.v2r = self._Mjtot[j]*v2rj/self.M
                        self.v2t = self._Mjtot[j]*v2tj/self.M

                        self.betaj = betaj
                        self.kj, self.krj, self.ktj = kj, krj, ktj
                        self.Kj, self.Krj = kj[-1], krj[-1]
                        self.ktj = self.kj - self.krj
                        self.Ktj = self.Kj - self.Krj
                        self.rhj, self.mcj = rhj, mcj
                    else:
                        self.rhohatj = numpy.vstack((self.rhohatj, rhohatj))
                        self.rhohat += self.alpha[j]*rhohatj

                        self.v2j = numpy.vstack((self.v2j, v2j))
                        self.v2rj = numpy.vstack((self.v2rj, v2rj))
                        self.v2tj = numpy.vstack((self.v2tj, v2tj))
                        self.v2 += self._Mjtot[j]*v2j/self.M
                        self.v2r += self._Mjtot[j]*v2rj/self.M
                        self.v2t += self._Mjtot[j]*v2tj/self.M

                        self.betaj = numpy.vstack((self.betaj, betaj))
                        self.kj = numpy.vstack((self.kj, kj))
                        self.krj = numpy.vstack((self.krj, krj))
                        self.ktj = numpy.vstack((self.ktj, ktj))
                        self.Kj = numpy.r_[self.Kj, kj[-1]]
                        self.Krj = numpy.r_[self.Krj, krj[-1]]
                        self.Ktj = numpy.r_[self.Ktj, ktj[-1]]
                        self.rhj = numpy.r_[self.rhj,rhj]
                        self.mcj = numpy.vstack((self.mcj, mcj))

                self.rho = self.rhohat*1.0
                self.rhoj = self.rhohatj*1.0

            # Calculate anisotropy profile (equation 32 of GZ15)
            self.beta = self._beta(self.r, self.v2r, self.v2t)

    def _rhohat(self, phi, r, j):
        """
        Wrapper for _rhoint when either: both phi or r are arrays, or both
        are scalar
        """
        if not hasattr(phi,"__len__"): phi = numpy.array([phi])
        if not hasattr(r,"__len__"): r = numpy.array([r])

        n = max([phi.size, r.size])
        rhohat = numpy.zeros(n)

        for i in range(n):
            if (phi[i]<self.max_arg_exp) or (numpy.isnan(phi[i])):
                rhohat[i] = self._rhoint(phi[i], r[i], self.raj[j])
                rhohat[i] /= self.rhoint0[j]
            else:
                # For large phi compute the ratio in one go (Section 4.1, GZ15)
                rhohat[i] = exp(phi[i]-self.phi0j[j]) if (self.multi) else 0

        return rhohat

    def _rhoint(self, phi, r, ra):
        """
        Dimensionless density integral as a function of phi and r (scalar)
        """

        # Isotropic case first (equation 8, GZ15)
        rhoint = exp(phi)*gammainc(self.g + 1.5, phi)

        # Anisotropic case, add r-dependent part explicitly (equation 11, GZ15)
        if (self.ra < self.ramax) and (phi>0) and (r>0):
            p, g = r/ra, self.g
            p2 = p**2
            g3, g5, fp2 = g+1.5, g+2.5, phi*p2

            func = hyp1f1(1, g5, -fp2)  if fp2 < self.max_arg_exp else g3/fp2
            rhoint += p2*phi**(g+1.5)*func/gamma(g5)
            rhoint /= (1+p2)
        return rhoint

    def _get_v2(self, phi, r, rho, j):
        v2 = numpy.zeros(r.size)
        v2r, v2t = numpy.zeros(r.size), numpy.zeros(r.size)
        for i in range(r.size-1):
            v2[i], v2r[i], v2t[i] = self._rhov2int(phi[i], r[i], self.raj[j])
            v2[i] /= rho[i]*self.rhoint0[j]
            v2r[i] /= rho[i]*self.rhoint0[j]
            v2t[i] /= rho[i]*self.rhoint0[j]

        return v2, v2r, v2t

    def _rhov2int(self, phi, r, ra):
        """Compute dimensionless pressure integral for phi, r """

        # Isotropic case first (equation 9, GZ15)
        rhov2r = exp(phi)*gammainc(self.g + 2.5, phi)
        rhov2  = 3*rhov2r
        rhov2t = 2*rhov2r

        # Add anisotropy, add parts depending explicitly on r
        # (see equations 12, 13, and 14 of GZ15)
        if (ra < self.ramax) and (r>0) and (phi>0):
            p, g = r/ra, self.g
            p2 = p**2
            p12 = 1+p2
            g3, g5, g7, fp2 = g+1.5, g+2.5, g+3.5, phi*p2

            P1 = p2*phi**g5/gamma(g7)
            H1 = hyp1f1(1, g7, -fp2) if fp2 <self.max_arg_exp else g5/fp2
            H2 = hyp1f1(2, g7, -fp2) if fp2 <self.max_arg_exp else g5*g3/fp2**2

            rhov2r += P1*H1
            rhov2r /= p12

            rhov2t /= p12
            rhov2t += 2*P1*(H1/p12 + H2)
            rhov2t /= p12

            rhov2 = rhov2r + rhov2t
        return rhov2, rhov2r, rhov2t

    def _beta(self, r, v2r, v2t):
        """ Calculate the anisotropy profile """

        beta = numpy.zeros(r.size)
        if (self.ra < self.ramax):
            c = (v2r>0.)
            # Equation  (32), GZ15
            beta[c] = 1.0 - 0.5*v2t[c]/v2r[c]
        return beta

    def _odes(self, x, y, potonly):
        """ Solve ODEs """
        # y = [phi, u_j, U, K_j], where u = -M(<r)/G
        if (self.multi):
            derivs = [numpy.sum(y[1:1+self.nmbin])/x**2] if (x>0) else [0]
            for j in range(self.nmbin):
                phi = y[0]/self.s2j[j]
                derivs.append(-9.0*x**2*self.alpha[j]*self._rhohat(phi, x, j))

            dUdx  = 2.0*pi*numpy.sum(derivs[1:1+self.nmbin])*y[0]/9.
        else:
            derivs = [y[1]/x**2] if (x>0) else [0]
            derivs.append(-9.0*x**2*self._rhohat(y[0], x, 0))

            dUdx  = 2.0*pi*derivs[1]*y[0]/9.

        derivs.append(dUdx)

        # Calculate pressure tensor components for all the mass bins
        if (not potonly): #dK_j/dx
            rhov2j, rhov2rj = [], []
            for j in range(self.nmbin):
                rv2, rv2r, rv2t = self._rhov2int(y[0]/self.s2j[j],
                                                 x, self.raj[j])
                P = self.alpha[j]*self.s2j[j]*2*pi*x**2*rv2/self.rhoint0[j]
                rhov2j.append(P)

                Pr = self.alpha[j]*self.s2j[j]*2*pi*x**2*rv2r/self.rhoint0[j]
                rhov2rj.append(Pr)

            for j in range(self.nmbin):
                derivs.append(rhov2j[j])

            for j in range(self.nmbin):
                derivs.append(rhov2rj[j])

        dVdvdr = (4*pi)**2*x**2 * (2*y[0])**1.5/3 if (x>0) and (y[0]>0) else 0
        derivs.append(dVdvdr)
        return derivs


    # ROT: UNDER CONSTRUCTIONS  ##############
    def _poisson_rot(self):
        self._init_rot()
        diff_rot = 1e99
        while (diff_rot > self.diff_rot_crit):
            # Compute density
            rho_l = self._rhohat_rot()
            self.rho_l = rho_l
            # Compute potential

            diff_rot = 0

    # def _phil_rot(self):
    #     phi_l = numpy.zeros(len(self.r)*nl).reshape(len(self.r), nl)

    #     phi_l[:,0] = self.phi0*sqrt(2)

    def legendre_Ul(self, l, theta):
        il = l/2.0
        dl = 1.0/2**l
        m = numpy.arange(il+1)
        cm = gamma(2*l-2*m+1)*(-1)**m/(gamma(m+1)*gamma(l-m+1)*gamma(l-2*m+1))
        gm = cos(theta)**(l-2*m)
        Pl = dl*numpy.sum(cm*gm)
        return sqrt((2*l+1)/2)*Pl

    def _rhohat_rot(self):
        for ip in range(len(self.r)):
            for it in range(self.ntheta):
                self.rhon_g[ip][it] = self._rhoint_rot(self.phin_g[ip][it], self.r[ip], self.theta[it])

        nl = int(self.lmax/2+1)
        rho_l = numpy.zeros(len(self.r)*nl).reshape(len(self.r), nl)

        for il in range(nl):
            l = 2*il
            for ir in range(len(self.r)):
                Ul = numpy.zeros(len(self.theta))
                for it in range(self.ntheta):
                    Ul[it] = self.legendre_Ul(l, self.theta[it])
                rho_l[ir][il] = -2.0 * simps(self.rhon_g[ir]*Ul, x=cos(self.theta))

        return rho_l

    def _rhoint_rot(self, phi, r, theta):
        # Compute eq 3 for scalar phi, r and theta
        g, om = self.g, self.omega
        Q = 3*sqrt(2)*r*sin(theta)
        # Eq 3
        integ = quad(lambda x: exp(phi-x)*gammainc(g, phi-x)*numpy.sinh(om*Q*x**0.5), 0, phi)[0]

        return integ/(om*Q*gamma(1.5))

    def _init_rot(self):
        self.lmax = 6
        self.diff_rot_crit = 1e-3

        # Grid params
        self.ntheta = (2*self.lmax + 1)

        # Compute spherical model first
        self._poisson(True)
        phin = self.phi
        self.theta = numpy.linspace(1e-4, 0.5*pi, self.ntheta)

        # Set up grid
        self.r_g, self.theta_g = numpy.meshgrid(self.r, self.theta, indexing ='ij')
        self.phin_g, tmp = numpy.meshgrid(self.phi, self.theta, indexing ='ij')
        self.rhon_g = self.phin_g*0.0

        # Fudge TBD
        self.r[0] = 0.5*self.r[1]

        # Make copy for iterations
        self.phin_prev_g = self.phin_g*1.0
        return
    ################

    def _setup_phi_interpolator(self):
        """ Setup interpolater for phi, works on scalar and arrays """

        # Generate piecewise 3th order polynomials to connect the discrete
        # values of phi obtained from from Poisson, using dphi/dr
        self._interpolator_set = True

        if (self.scale):
            phi_and_derivs = numpy.vstack([[self.phi],[self.dphidr1]]).T
        else:
            phi_and_derivs = numpy.vstack([[self.phihat],[self.dphidrhat1]]).T

        self._phi_poly = BPoly.from_derivatives(self.r,phi_and_derivs)

    def _scale(self):
        """ Scales the model to the units set in the input: GS, MS, RS """

        Mstar = self._MS/self.M
        Gstar = self._GS/self.G
        if (self.scale_radius=='r0'): Rstar = self._RS/self.r0
        if (self.scale_radius=='rh'): Rstar = self._RS/self.rh
        if (self.scale_radius=='rv'): Rstar = self._RS/self.rv
        if (self.scale_radius=='rt'): Rstar = self._RS/self.rt
        v2star =  Gstar*Mstar/Rstar

        # Update the scales that define the system (see Section 2.1.2 of GZ15)
        self.G *= Gstar
        self.rs = Rstar
        self.s2 *= v2star
        self.s2j *= v2star

        # Anisotropy radii
        self.ra, self.raj, self.ramax = (q*Rstar for q in [self.ra,
                                                           self.raj,
                                                           self.ramax])

        # Scale all variable needed when run with potonly=True
        self.r, self.r0, self.rt = (q*Rstar for q in [self.rhat,
                                                      self.r0, self.rt])
        self.rh, self.rv = (q*Rstar for q in [self.rh,self.rv])

        self.M *= Mstar
        self.phi = self.phihat * v2star
        self.dphidr1 = self.dphidrhat1 * v2star/Rstar
        self.mc *= Mstar
        self.U *= Mstar*v2star
        self.A *= Mstar/(v2star**1.5*Rstar**3)
        self.volume *= v2star**1.5*Rstar**3


        # Scale density, velocity dispersion components, kinetic energy
        self.rho = self.rhohat*Mstar/Rstar**3

        if (not self.potonly):
            self.rho = self.rhohat*Mstar/Rstar**3
            self.v2, self.v2r, self.v2t = (q*v2star for q in [self.v2,
                                                          self.v2r,self.v2t])
            self.K,self.Kr,self.Kt=(q*Mstar*v2star for q in [self.K, self.Kr,
                                                             self.Kt])

            if (self.multi):
                self.rhoj = self.rhohatj * Mstar/Rstar**3
                for j in range(self.nmbin):
                    self.rhoj[j] *= self.alpha[j]
                self.Mj *= Mstar # 28/1/19 Thanks to William

                self.mcj *= Mstar
                self.rhj *= Rstar
                self.v2j,self.v2rj,self.v2tj=(q*v2star for q in
                                              [self.v2j, self.v2rj,self.v2tj])

                self.kj *= Mstar*v2star
                self.Krj *= Mstar*v2star
                self.Ktj *= Mstar*v2star
                self.Kj *= Mstar*v2star

    def _tonp(self, q):
        q = numpy.array([q]) if not hasattr(q,"__len__") else numpy.array(q)
        return q

    def _project(self):
        """
        Compute projected mass density (Sigma) and projected <v2> profiles
        """

        # Initialise the projected quantities:
        # R is the projected (2d) distance from the center, Sigma is the
        # projected density, v2p is the line-of-sight velocity dispersion,
        # v2R and v2T are the radial and tangential velocity dispersion
        # components projected on the plane of the sky

        # Initialise some arrays
        R = self.r
        Sigma = numpy.zeros(self.nstep)
        v2p = numpy.zeros(self.nstep)
        v2R = numpy.zeros(self.nstep)
        v2T = numpy.zeros(self.nstep)
        mcp = numpy.zeros(self.nstep)

        if (self.multi):
            Sigmaj = numpy.zeros((self.nmbin, self.nstep))
            v2pj = numpy.zeros((self.nmbin, self.nstep))
            v2Rj = numpy.zeros((self.nmbin, self.nstep))
            v2Tj = numpy.zeros((self.nmbin, self.nstep))
            mcpj = numpy.zeros((self.nmbin, self.nstep)) # TBD

        # Project model properties for each R
        for i in range(self.nstep-1):
            c = (self.r >= R[i])
            r = self.r[c]
            z = sqrt(abs(r**2 - R[i]**2)) # avoid small neg. values

            Sigma[i] = 2.0*abs(simps(self.rho[c], x=z))

            betaterm1 = 1 if i==0 else 1-self.beta[c]*R[i]**2/self.r[c]**2

            # eq 41 in paper has a small mistake:  (1-R^2)/r^2 should be (1-R^2/r^2) as below
            betaterm2 = 1 - self.beta[c] if i==0 else 1-self.beta[c]*(1-R[i]**2/self.r[c]**2)

            v2p[i] = abs(2.0*simps(betaterm1*self.rho[c]*self.v2r[c], x=z))
            v2p[i] /= Sigma[i]

            v2R[i] = abs(2.0*simps(betaterm2*self.rho[c]*self.v2r[c], x=z))
            v2R[i] /= Sigma[i]

            v2T[i] = abs(simps(self.rho[c]*self.v2t[c], x=z))
            v2T[i] /= Sigma[i]

            # Cumulative mass in projection
            if (i>0):
                x = self.r[i-1:i+1]
                mcp[i] = mcp[i-1] + 2*pi*simps(x*Sigma[i-1:i+1], x=x)
            mcp[-1] = mcp[-2]

            # Radius containing half the mass in projection
            self.rhp = numpy.interp(0.5*mcp[-1], mcp, self.r)

            if (self.multi):
                for j in range(self.nmbin):
                    Sigmaj[j,i] = 2.0*simps(self.rhoj[j,c], x=z)
                    if (i==0):
                        betaterm1 = 1
                        betaterm2 = 1 -self.betaj[j,c]
                    else:
                        betaterm1 = 1-self.betaj[j,c]*R[i]**2/self.r[c]**2
                        # eq 41 in paper has a small mistake:  (1-R^2)/r^2 should be (1-R^2/r^2) as below
                        betaterm2 = 1-self.betaj[j,c]*(1-R[i]**2/self.r[c]**2)

                    v2int = simps(betaterm1*self.rhoj[j,c]*self.v2rj[j,c], x=z)
                    v2pj[j,i] = abs(2.0*v2int)
                    v2pj[j,i] /= Sigmaj[j,i]

                    v2int = simps(betaterm2*self.rhoj[j,c]*self.v2rj[j,c], x=z)
                    v2Rj[j,i] = abs(2.0*v2int)
                    v2Rj[j,i] /= Sigmaj[j,i]

                    v2Tj[j,i] = abs(simps(self.rhoj[j,c]*self.v2tj[j,c], x=z))
                    v2Tj[j,i] /= Sigmaj[j,i]

        self.R, self.Sigma = R, Sigma
        self.v2p, self.v2R, self.v2T = v2p, v2R, v2T
        self.mcp = mcp

        # TBD: Compute half-mass radii in projection

        if (self.multi):
            self.Sigmaj = Sigmaj
            self.v2pj, self.v2Rj, self.v2Tj =  v2pj, v2Rj, v2Tj


        return

    def get_Paz(self, az_data, R_data, jns):
        """
        Computes probability of line of sight acceleration at projected R : P(az|R)
        """

        # Under construction !!!

        # Return P(az|R)

        az_data = abs(az_data) # Consider only positive values

        # Assumes units of az [m/s^2] if self.G ==0.004302, else models units
        # Conversion factor from [pc (km/s)^2/Msun] -> [m/s^2]
        az_fac = 1./3.0857e10 if (self.G==0.004302) else 1

        if (R_data < self.rt):
            nz = self.nstep                   # Number of z value equal to number of r values
            zt = sqrt(self.rt**2 - R_data**2) # maximum z value at R

            z = numpy.logspace(log10(self.r[1]), log10(zt), nz)

            spl_Mr = UnivariateSpline(self.r, self.mc, s=0, ext=1)  # Spline for enclosed mass

            r = sqrt(R_data**2 + z**2)                        # Local r array
            az = self.G*spl_Mr(r)*z/r**3                     # Acceleration along los
            az[-1] = self.G*spl_Mr(self.rt)*zt/self.rt**3    # Ensure non-zero final data point

            az *= az_fac # convert to [m/s^2]
            az_spl = UnivariateSpline(z, az, k=4, s=0, ext=1) # 4th order needed to find max (can be done easier?)

            zmax = az_spl.derivative().roots()  # z where az = max(az), can be done without 4th order spline?
            azt = az[-1]                        # acceleration at the max(z) = sqrt(r_t**2 - R**2)

            # Setup spline for rho(z)
            if jns == 0 and self.nmbin == 1:
                rho = self.rho
            else:
                rho = self.rhoj[jns]

            rho_spl = UnivariateSpline(self.r, rho, ext=1, s=0)
            rhoz = rho_spl(sqrt(z**2 + R_data**2))
            rhoz_spl = UnivariateSpline(z, rhoz, ext=1, s=0)

            # Now compute P(a_z|R)
            # There are 2 possibilities depending on R:
            #  (1) the maximum acceleration occurs within the cluster boundary, or
            #  (2) max(a_z) = a_z,t (this happens when R ~ r_t)

            nr, k = nz, 3 # bit of experimenting

            # Option (1): zmax < max(z)
            if len(zmax)>0:
                zmax = zmax[0] # Take first entry for the rare cases with multiple peaks
                # Set up 2 splines for the inverse z(a_z) for z < zmax and z > zmax
                z1 = numpy.linspace(z[0], zmax, nr)
                z2 = (numpy.linspace(zmax, z[-1], nr))[::-1] # Reverse z for ascending az

                z1_spl = UnivariateSpline(az_spl(z1), z1, k=k, s=0, ext=1)
                z2_spl = UnivariateSpline(az_spl(z2), z2, k=k, s=0, ext=1)

            # Option 2: zmax = max(z)
            else:
                zmax = z[-1]
                z1 = numpy.linspace(z[0], zmax, nr)
                z1_spl = UnivariateSpline(az_spl(z1), z1, k=k, s=0, ext=1)

            # Maximum acceleration along this los
            azmax = az_spl(zmax)

            # Now determine P(az_data|R)
            if (az_data < azmax):
                z1 = max([z1_spl(az_data), z[0]]) # first radius where az = az_data
                Paz = rhoz_spl(z1)/abs(az_spl.derivatives(z1)[1])

                if (az_data> azt):
                    # Find z where a_z = a_z,t
                    z2 = z2_spl(az_data)
                    Paz += rhoz_spl(z2)/abs(az_spl.derivatives(z2)[1])

                # Normalize to 1
                Paz /= rhoz_spl.integral(0, zt)
                self.z = z
                self.az = az
                self.Paz = Paz
                self.azmax = azmax
                self.zmax = zmax
            else:
                self.Paz = 0
        else:
            self.Paz = 0

        return


    def interp_phi(self, r):
        """ Returns interpolated potential at r, works on scalar and arrays """

        if not hasattr(r,"__len__"): r = numpy.array([r])
        if (not self._interpolator_set): self._setup_phi_interpolator()

        phi = numpy.zeros([r.size])
        inrt = (r<self.rt)
        # Use 3th order polynomials to interp, using phi'
        if (sum(inrt)>0): phi[inrt] = self._phi_poly(r[inrt])
        return phi

    def df(self, *arg):
        """
        Returns the value of the normalised DF at a given position in phase
        space, can only be called after solving Poisson's equation
        Arguments can be:
          - r, v                   (isotropic single-mass models)
          - r, v, j                (isotropic multi-mass models)
          - r, v, theta, j         (anisotropic models)
          - x, y, z, vx, vy, vz, j (all models)
        Here j specifies the mass bin, j=0 for single mass
        Works with scalar and array input
        """

        if (len(arg)<2) or (len(arg)==5) or (len(arg)==6) or (len(arg)>7):
            raise ValueError("Error: df needs 2, 3, 4 or 7 arguments")

        if (len(arg)<=3) and (self.ra<self.ramax):
            raise ValueError("Error: model is anisotropy, more input needed")

        if len(arg) == 2:
            r, v = (self._tonp(q) for q in arg)
            j = 0

        if len(arg) == 3:
            r, v = (self._tonp(q) for q in arg[:-1])
            j = arg[-1]

        if len(arg) == 4:
            r, v, theta = (self._tonp(q) for q in arg[:-1])
            j = arg[-1]

        if len(arg) < 7: r2, v2 = r**2, v**2

        if len(arg) == 7:
            x, y, z, vx, vy, vz = (self._tonp(q) for q in arg[:-1])
            j = arg[-1]
            r2 = x**2 + y**2 + z**2
            v2 = vx**2 + vy**2 + vz**2
            r, v = sqrt(r2), sqrt(v2)

        # Interpolate potential and calculate escape velocity as a
        # function of the distance from the center
        phi = self.interp_phi(r)
        vesc2 = 2.0*phi                  # Note: phi > 0

        DF = numpy.zeros([max(r.size, v.size)])

        if (len(r) == len(v2)):
            # Bitwise & needed to allow for possibility that r and v2 are arrays
            c = (r<self.rt) & (v2<vesc2)

        if (len(r) == 1) and (len(v2) > 1):
            c = (v2<vesc2)
            if (r>self.rt): c=False

        if (len(r) > 1) and (len(v2) == 1):
            # Bitwise & needed to allow for possibility that r and v2 are arrays
            c = (r<self.rt) & (numpy.zeros(len(r))+v2<vesc2)

        if (sum(c)>0):
            # Compute the DF: equation (1), GZ15

            E = (phi-0.5*v2)/self.s2j[j]  # Dimensionless positive energy
            DF[c] = exp(E[c])

            if (self.g>0): DF[c] *= gammainc(self.g, E[c])

            if (len(arg)==7): J2 = v2*r2 - (x*vx + y*vy + z*vz)**2
            if (len(arg)==4): J2 = sin(theta)**2*v2*r2
            if (len(arg)<=3): J2 = numpy.zeros(len(c))

            DF[c] *= exp(-J2[c]/(2*self.raj[j]**2*self.s2j[j]))

            DF[c] *= self.A[j]

        else:
	        DF = numpy.zeros(max(len(r),len(v)))

        return DF
