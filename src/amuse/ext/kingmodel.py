import math
import numpy

from amuse.support import exceptions
from amuse.units import nbody_system

from amuse.support import data
class MakeKingModel(object):
    def __init__(self, number_of_particles, W0, convert_nbody = None, do_scale = False, 
            beta = 0.0, seed = None, verbose = False):
        self.number_of_particles = number_of_particles
        self.convert_nbody = convert_nbody
        self.do_scale = do_scale
        self.verbose = verbose
        self.beta = beta
        self.W0 = W0
        self.beta_w0 = beta*W0
        self.scale_fac = math.exp(self.beta_w0)

        self.YMAX = 4.0 # Note: make sure YMAX is a float.
        self.NG = 1000
        self.g_integral = self.compute_g_integral()
        self.v33 = self.compute_v33()
        
        # // profile
        self.NM = 2500
        self.rr=[]; self.d=[]; self.v2=[]; self.psi=[]; self.zm=[]
        self.NINDX = 100
        self.index=[]
        
        numpy.random.seed(seed)
        
    def compute_v33(self):
        v33 = []
    	for i in range(self.NG+1):
            v33.append(self.scale_fac * math.pow(((self.YMAX/self.NG) * i), 3) / 3.0)
        return v33
    
    def compute_g_integral(self):
        dy = self.YMAX/self.NG
        g = [0.0]
        for i in range(self.NG):
            g.append(g[i]+self.trapzd_gaus2(i*dy,(i+1)*dy))
        return g
            
    # Trapezoid Rule for integration, after NR.
    def trapzd_gaus2(self, a, b):
        eps = 1.e-10
        max_steps = 50
        previous_sum = -1.0e30
        for i in range(max_steps):
            sum = self.trapint_gaus2(a, b, i, previous_sum);
            if (abs(sum - previous_sum) < eps * abs(previous_sum)): return sum
            previous_sum = sum;
        return 0.0
        
    # Integrate func from a to b in n steps, using the trapezoid rule.
    def trapint_gaus2(self, a, b, n, previous_sum):
        if (n == 0):
            return 0.5 * (b - a) * (self.gaus2(a) + self.gaus2(b))
        else:
            base = 2**(n-1)
            dx = (b - a) / base
            x = a + 0.5 * dx
            sum = 0.0
            for i in range(base):
                sum = sum + self.gaus2(x)
                x = x + dx
            return 0.5 * (previous_sum + (b - a)/base * sum);
            
    def gaus2(self, x):
        return x*x*math.exp(-x*x)
        
    def gg_integral(self, y, g4_flag):
        # Array g contains the indefinite integral (from 0) of y^2 exp(-y^2) dy,
        # up to a maximum value of YMAX, i.e.
        #
        #	g[i]  =  integral{0 to Y} y^2 exp(-y^2) dy,
        #
        # where Y  =  i*YMAX/NG (i = 0,...,NG).
        #
        # Function computes  g2  =  integral{0 to y} y^2 exp(-y^2) dy
        #               and  g4  =  integral{0 to y} y^4 exp(-y^2) dy.
        if (y >= self.YMAX):
            g2 = self.g_integral[self.NG]
        else:
            yindx = self.NG * y / self.YMAX
            intyindx = int(yindx)
            g2 = self.g_integral[intyindx] + (yindx-intyindx)*(self.g_integral[intyindx+1]-self.g_integral[intyindx])
        if (g4_flag):
            return (g2, 1.5*g2 - 0.5*y**3*math.exp(-y*y))
        else:
            return (g2, 0.0)
    
    # Return density d and local velocity dispersion v2,
    # given scaled potential psi (< 0).
    def get_dens_and_vel(self, psi, v2_flag):
        if (psi >= -self.beta_w0): 
            if (v2_flag): return (0.0, 1.0)
            else: return 0.0
        #  // Distribution function is
        #  //
        #  //	f(r,v) = A exp(-psi) (exp(-v^2 / 2 sig^2) - exp(psi - beta*psi0)),
        #  //
        #  // where psi = phi/sig^2 and psi0 = -w0.
        #  //
        #  // Density is the integral of 4 pi v^2 f, velocity dispersion v2
        #  // is (integral of 4 pi v^4 f) / density.
        p_max = -psi - self.beta_w0	# // v_max^2 / (2 sig^2)
        y_max = math.sqrt(p_max)
        ep = math.exp(-psi)
        (g2, g4) = self.gg_integral(y_max, v2_flag)
        dens = ep * g2 - self.scale_fac * y_max * p_max / 3
        #  // Note that this expression for dens excludes an overall factor
        #  // of 8 pi sqrt(2) sig^3 A.  Scaling to the central density
        #  // is handled elsewhere (in rhs).
        if (v2_flag):
            v2 = 2 * (ep * g4 - 0.2 * self.scale_fac * y_max * p_max * p_max) / dens
            return (dens, v2)
        else:
            return dens
        #  // The unit of v2 is sig^2.
    
    def rhs(self, y, x):
        # //  Define RHS of ODE, for use by rk4.
        ypr = []
        ypr.append(y[1])
        if (x <= 0):
            d = 1
        else:
            d = (self.get_dens_and_vel(y[0]/x, False)) * self.dc_inverse
            # // d is rho/rho_0
            # // False suppresses vel
        ypr.append(9 * x * d)
        return ypr
    
    # //  Runge-Kutta-4 method 
    def step(self, y, x, dx, N):
        # Doesn't look very efficient...
        dy1=[]; dy2=[]; dy3=[]; dy4=[]
        y1=[]; y2=[]; y3=[]; y4=[]
        dydx = self.rhs(y, x)
        for i in range(N):
            dy1.append(dx*dydx[i])
            y1.append(y[i] + 0.5*dy1[i])
        dydx = self.rhs(y1, x+0.5*dx)
        for i in range(N):
            dy2.append(dx*dydx[i])
            y2.append(y[i] + 0.5*dy2[i])
        dydx = self.rhs(y2, x+0.5*dx)
        for i in range(N):
            dy3.append(dx*dydx[i])
            y3.append(y[i] + dy3[i])
        dydx = self.rhs(y3, x+dx)
        for i in range(N):
            dy4.append(dx*dydx[i])
            y[i] = y[i] + (dy1[i] + 2*dy2[i] + 2*dy3[i] + dy4[i])/6.0
        return y
    
    def rk4(self, x, x_next, y, N, dx):
        # // Integrate the y(x) from x to x_next, using an integration step of dx.
        # // On entry y is y at the initial x.  On exit, x is replaced by x_next
        # // and y is y(x).
        while (x < x_next - .5*dx): # // NOTE!!  Beware of rounding error!
            y = self.step(y, x, dx, N)
            x = x + dx
        return (x, y)
    
    def poisson(self):
        # //       Self-contained 1-D (spherical) Poisson's equation solver.
        # //       Currently knows about normal and lowered King models.
        # //
        # //        Input:  nmax is the maximum number of points allowed
        # //                w0 is the dimensionless central potential
        # //                iout allows messages if nonzero
        #
        # //        Output: x   is scaled radius (r/rc)
        # //                d   is scaled density (1 at center)
        # //                v2  is scaled velocity dispersion (1 at center)
        # //                psi is scaled potential (-W0 at center)
        # //                zm  is cumulative mass (scaling from x, d scalings)
        # //                nprof is the actual number of points generated
        # //                v20 is the central 3-D velocity dispersion (unit = sig^2)
        RLIN = 0.25
        NLIN = 105
        RMAX = 1e4
        
        nmax = self.NM
        psi0 = - abs(self.W0)
        #  // Initialize at center of cluster.
        y = [0.0, psi0]
        xn = 0.0
        self.rr.append(0.0)
        self.psi.append(psi0)
        self.zm.append(0.0)
        #  // Establish density scaling factor.
        (density, velocity_dispersion) = self.get_dens_and_vel(psi0, True)
        self.d.append(density)
        self.v2.append(velocity_dispersion)
        self.dc_inverse = 1./density
        fac = math.pow(10, (math.log10(RMAX/RLIN) / (nmax-NLIN)))
        #// 	Poisson's equation is:
        #//
        #// 		(1/r^2) d/dr (r^2 dphi/dr)  =  4 pi G rho,
        #//
        #// 	where r is radius, phi is potential, and rho is density, given
        #// 	(for equal-mass stars) by
        #//
        #// 		rho	=  {integral (v < ve)} 4 pi v^2 f(v) dv,
        #//
        #// 	where ve is the escape velocity,
        #//
        #// 		ve^2	=  -2 phi.
        #//
        #// 	The (3-D) velocity distribution is
        #//
        #// 		f(v)	=  A (exp(-v^2 / 2 sig^2)
        #// 					 - exp(-ve^2 / 2 sig^2)),
        #//
        #// 	where sig^2 is a 1-D velocity dispersion (not quite the
        #// 	central velocity dispersion, except in the limit of infinite
        #// 	central potential).  In King's (1966) paper, he uses
        #// 	j^2 = 1 / (2 sig^2).
        #//
        #// 	Following King, we define the core radius rc by
        #//
        #// 		rc^2	=  9 sig^2 / (4 pi G rho0)
        #//
        #// 	and the dimensionless depth as
        #//
        #// 		W0	=  -phi0 / sig^2,
        #//
        #// 	where rho0 and phi0 are the central density and potential,
        #// 	respectively.
        #//
        #// 	We then scale as follows:
        #//
        #// 		x	=  r / rc
        #//
        #// 		d	=  rho / rho0
        #//
        #// 		psi	=  phi / sig^2,
        #//
        #// 	to obtain
        #//
        #// 		(x psi)''  =  9 x d,
        #//
        #// 	where ' = d/dx.
        #//
        #// 	We integrate this ODE from the cluster center (x = 0, d = 1,
        #// 	psi = -W0) to the tidal radius (d = 0) by defining
        #//
        #//		y(0)	=  (x psi)
        #//		y(1)	=  y(0)'
        #//
        #//	We cover the first RLIN core radii linearly with NLIN points;
        #//	the remaining coverage is logarithmic, out to RMAX core radii,
        #//	if necessary.  We stop when d <= 0.
        stopped = False
        for i in range(1,nmax+1):
            xo = xn
            if (i <= NLIN):
                xn = (RLIN * i) / NLIN
            else:
                xn = fac * xo
            dx = 0.051*(xn-xo)
            (xo, y) = self.rk4(xo, xn, y, 2, dx)
            # //  N.B. Remember that y(1) is x*psi and xo is updated by step.
            xn = xo
            self.rr.append(xn)
            self.psi.append(y[0] / xn)
            (density, velocity_dispersion) = self.get_dens_and_vel(self.psi[i], True)
            self.d.append(density)
            self.v2.append(velocity_dispersion)
            if (density < 0):
                # 	// Density is negative, calculation is over.
                # 	// Interpolate to the tidal radius.
                self.rr[i] = self.rr[i-1] + (self.rr[i] - self.rr[i-1]) / (0.1 - self.d[i]/self.d[i-1])
                self.d[i] = 0
                self.v2[i] = 0
            self.zm.append(self.rr[i] * y[1] - y[0])
            if (density <= 0):
                stopped = True
                break
        if not stopped: i = nmax
        nprof = i
        #  // Scale d and v2 to their central values.  Save v2_0 (unit = sig^2).
        v2_0 = self.v2[0]
        d_0 = self.d[0]
        c_zm = 4.0 * math.pi / 9.0
        self.d[:] = [x/d_0 for x in self.d]
        self.v2[:] = [x/v2_0 for x in self.v2]
        self.zm[:] = [x*c_zm for x in self.zm]
        return nprof, v2_0
    
    def coordinates_from_spherical(self, radius, theta, phi):
        x = radius * numpy.sin( theta ) * numpy.cos( phi )
        y = radius * numpy.sin( theta ) * numpy.sin( phi )
        z = radius * numpy.cos( theta )
        return [x,y,z]
    
    def setpos(self):
        # // Obtain a random position for body b from the King profile
        # // and return the scaled potential at that location.
        
        #  //  Choose radius randomly from the mass distribution.
        rno = numpy.random.uniform()
        i = int(self.NINDX * rno)
        found_index = False
        for i1 in range(self.index[i],self.index[i+1]+2): #(i1 = self.indx[i]; i1 <= indx[i+1]+1; i1++)
            if (self.zm[i1] > rno):
                found_index = True
                break
        if (not found_index):
            raise exceptions.AmuseException("makeking: error in getpos")
        rfac = (rno - self.zm[i1-1]) / (self.zm[i1] - self.zm[i1-1])
        radius = self.rr[i1-1] + rfac * (self.rr[i1] - self.rr[i1-1])
        potential = self.psi[i1-1] + rfac * (self.psi[i1] - self.psi[i1-1])
        
        #  //  Angular position random.
        theta = numpy.arccos(numpy.random.uniform(-1.0, 1.0))
        phi = numpy.random.uniform(0.0, 2.0*math.pi)
        return self.coordinates_from_spherical(radius, theta, phi), potential
    
    def setvel(self, potential):
        #// Obtain a random velocity for body b from the King profile
        #// given its scaled potential.
        
        #    // Array v33[] contains the second term in the integral for the density,
        #    // namely exp(beta*W0) * v_esc^3 / 3 (scaling v^2 by 2 sig^2, as usual).
        #    // As with g[], the array indices run from 0 to NG, spanning a range 0 to
        #    // YMAX, i.e. the cumulative distribution function for v is
        #    //
        #    //		exp(-p) * g[i] - v33[i],
        #    //
        #    // where y = i*YMAX/NG (i = 0,...,NG) and v = sqrt(2)*sig*y (sig = 1 here).

        #    //  Choose speed randomly from the distribution at this radius.
        v = 0
        if (potential < -self.beta_w0):
            pfac = math.exp(-potential)
            ymax = math.sqrt(-potential-self.beta_w0)
            #	// Will obtain v by bisection.  Determine maximum possible
            #	// range in the index i.
            il = 0
            iu = int(((self.NG/self.YMAX) * math.sqrt(-potential))) #	// Binning OK for W0 < 16,
            #   				        		// *only* if beta >= 0.
            if (iu > self.NG): 
                iu = self.NG
            rl = 0
            ru = pfac * self.g_integral[iu] - self.v33[iu]
            rno = numpy.random.uniform(0.0, ru)
            while (iu - il > 1):
                im = (il + iu) / 2
                rm = pfac * self.g_integral[im] - self.v33[im]
                if (rm > rno):
                    iu = im
                    ru = rm
                else:
                    il = im
                    rl = rm
            #	// Maximum possible range of il here (for beta = 0) is
            #	//	0 to NG*sqrt(-p)/YMAX.
            #	// Maximum possible value of v (for beta = 0) is the local
            #	//      escape speed, sqrt(-2*p).
            v = (self.YMAX/self.NG) * math.sqrt(2.0) * (il + (rno - rl)/(ru - rl))
        #    //  Direction is random.
        theta = numpy.arccos(numpy.random.uniform(-1.0, 1.0))
        phi = numpy.random.uniform(0.0, 2.0*math.pi)
        return self.coordinates_from_spherical(v, theta, phi)
        
    def makeking(self):
        #// Create a King model, and optionally initialize an N-body system
        #// with total mass = n, core radius = 1.
        if (self.W0 > 16): 
            raise exceptions.AmuseException("makeking: must specify w0 < 16")
        #    // Compute the cluster density/velocity/potential profile
        (nprof, v20) = self.poisson()
        zm = self.zm
        d = self.d
        rr = self.rr
        v2 = self.v2
        psi = self.psi
        if not (len(zm) == nprof+1 and len(d) == nprof+1 and 
            len(rr) == nprof+1 and len(v2) == nprof+1 and 
            len(psi) == nprof+1):
            print len(zm), len(d), len(rr), len(v2), len(psi), nprof+1
            raise exceptions.AmuseException("Error in result of Poisson")
        #    // Determine statistics and characteristic scales of the King model.
        rho0 = 1 / zm[nprof]#	 // Central density for total mass = 1
        #    // Unit of velocity = sig, where rc^2 = 9 sig^2 / (4 pi G rho0)
        sig = math.sqrt(4.0 * math.pi * rho0 / 9.0)# // This 3 was v20 in the f77 version...
        #					 // v20 is central vel. disp. / sig^2
        #    // Scale the zm array to unit total mass.
        inv_total_mass = 1.0 / zm[nprof]
        zm[:] = [x*inv_total_mass for x in zm]
        #    // Index the mass distribution, and determine the core mass and
        #    // the half-mass radius.
        
        #    // By construction, rr[indx[j]] and rr[indx[j+1]] bracket the
        #    // radius containing a fraction j / NINDX of the total mass.
        self.index.append(0)
        dz = 1.0/self.NINDX
        z = dz
        for j in range(1,nprof):
            if (rr[j] < 1): jcore = j
            if (zm[j] < 0.5): jhalf = j
            if (zm[j] > z):
                self.index.append(j - 1)
                z = z + dz
        self.index.append(nprof)
        if not (len(self.index)==self.NINDX+1):
            raise exceptions.AmuseException("Error in length of indx")
        zmcore = zm[jcore] + (zm[jcore+1] - zm[jcore]) * (1 - 
            rr[jcore]) / (rr[jcore+1] - rr[jcore])
        rhalf = rr[jhalf] + (rr[jhalf+1] - rr[jhalf]) * (0.5 -  
            zm[jhalf]) / (zm[jhalf+1] - zm[jhalf])
        #    // Compute the kinetic and potential energies, and determine the
        #    // virial radius and ratio.
        kin = 0; pot =0
        for i in range(nprof):
            kin = kin + (zm[i+1] - zm[i]) * (v2[i+1] + v2[i])
            pot = pot - (zm[i+1] - zm[i]) * (zm[i+1] + zm[i]) / (rr[i+1] + rr[i])
        kin *= 0.25*sig*sig*v20
        rvirial = -0.5/pot
        #    // Initialize the N-body system.
        if self.verbose:
            print " King model, w0 = ",self.W0,", Rt/Rc = ",rr[nprof],", Rh/Rc = ",rhalf,", Mc/M = ", zmcore
            #    // Write essential model information
            print "initial_mass", 1.0
            print "initial_rtidal_over_rvirial",  rr[nprof] / (0.25/kin)
        #    // Assign positions and velocities. Note that it may actually
        #    // be preferable to do this in layers instead.
        masses = numpy.zeros((self.number_of_particles,1)) + (1.0 / self.number_of_particles)
        positions = []
        velocities = []
        #    // Convenient to have the "unscaled" system 
        #    // be as close to standard units as possible, so rescale position
        #    // and velocity with 'xfac' and 'vfac' to force
        #    // the virial radius to 1.  (Steve, 9/04)
        xfac = 1.0/rvirial
        vfac = 1.0/math.sqrt(xfac)
        for i in range(self.number_of_particles):
            (position, potential) = self.setpos()
            velocity = self.setvel(potential)
            #	// Unit of length = rc.
            #	// Unit of velocity = sig.
            position[:] = [xfac*comp for comp in position]
            velocity[:] = [vfac*comp*sig for comp in velocity]
            positions.append(position)
            velocities.append(velocity)
        #    // System is in virial equilibrium in a consistent set of units
        #    // with G, core radius, and total mass = 1.
        
        return (masses, positions, velocities)
   
    @property
    def result(self):
        masses, positions, velocities = self.makeking()
        result = data.Particles(self.number_of_particles)
        result.mass = nbody_system.mass.new_quantity(masses)
        result.position = nbody_system.length.new_quantity(positions)
        result.velocity = nbody_system.speed.new_quantity(velocities)
        result.move_to_center()
        if self.do_scale:
            result.scale_to_standard()
        
        if not self.convert_nbody is None:
            result = data.ParticlesWithUnitsConverted(result, self.convert_nbody.as_converter_from_si_to_nbody())
            result = result.copy_to_memory()
        
        return result
    

def new_king_model(number_of_particles, W0, *list_arguments, **keyword_arguments):
    """
    Create a King model with the given number of particles and King dimensionless 
    depth W0. Returns a set of particles with equal mass and positions and velocities 
    distributed to fit a King distribution model. The model is centered around the 
    origin. Positions and velocities are optionally scaled such that the kinetic and 
    potential energies are 0.25 and -0.5 in nbody-units, respectively.

    :argument number_of_particles: Number of particles to include in the King model
    :argument W0: King dimensionless depth, allowed range: < 0, 16 ]
    :argument convert_nbody:  When given will convert the resulting set to SI units
    :argument do_scale: scale the result to exact nbody units (M=1, K=0.25, U=-0.5)
    :argument beta:  Steve's rescaling parameter (< 1) [0]. Models with b > 0 are just 
    :argument seed:  Seed for the random number generator
        rescaled King models; models with b < 0 approach isothermal spheres as 
        b --> -infinity.
    :argument verbose: Be verbose (output is suppressed by default) [False]
    """
    uc = MakeKingModel(number_of_particles, W0, *list_arguments, **keyword_arguments)
    return uc.result
