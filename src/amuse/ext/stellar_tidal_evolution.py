# Equations for tidal evolution of binary/planetary orbits
# From Hansen 2010, based on Eggleton et al. 1998

import threading
import numpy
import math
from scipy import integrate
from amuse.lab import *
from amuse.support import literature

#2010ApJ...723..285H
def sigma_planet():
    sig_p = 3.4e-7 
    return 5.9e-54 * sig_p | units.g**-1 *units.cm**-2 * units.s**-1

#2010ApJ...723..285H
def sigma_star():
    sig_s = 7.8e-8 
    #sig_s = 1e+9 # Earth
    return 6.4e-59 * sig_s  | units.g**-1 *units.cm**-2 * units.s**-1

def angular_frequency(Ms, Mp, a):
    return (constants.G*(Ms+Mp)/(a**3))**(0.5)

def J_orb(Ms, Mp, a, e):
    return Ms*Mp*(constants.G*a*(1-e**2)/(Ms + Mp))**0.5
   
def interp(x, *args, **kwds):
     if type(x) in (float, int):
         return numpy.interp([x], *args, **kwds).item()
     else :
         return numpy.interp(x, *args, **kwds)

class TidalEvolution(literature.LiteratureReferencesMixIn):
    """
        Tidal evolution between a planet(esimal) and a star.
        Based on
        .. [#] ** 2010ApJ...723..285H
        .. [#] ** 
    """
    def __init__(self, central_particle=Particles()):
        literature.LiteratureReferencesMixIn.__init__(self)
        self.current_time = 0 | units.s
        self.pericenter_interaction_factor = 4
        self.central_particle = central_particle
        if not hasattr(self.central_particle, "gyration_radius_sq"):
            self.central_particle.gyration_radius_sq = 0.2
        if hasattr(self.central_particle, "radius"):
            self.central_particle.old_radius = self.central_particle.radius
        self.orbiters = Particles(0)
        if not hasattr(self.central_particle, "Omega"):
            self.central_particle.Omega = 2.6e-6|units.s**-1
        self.all_merged_orbiters = Particles()
        
    @property
    def particles(self):
        return ParticlesSuperset([self.central_particle, self.orbiters])

    def set_current_time(self, time):
        self.current_time = time
    def get_semimajor_axis(self):
        return self.semimajor_axis
    def get_eccentricity(self):
        return self.eccentricity
    def get_pericenter_interaction_factor(self):
        return self.pericenter_interaction_factor

    def orbital_evolution_time_scale(self):
        e = self.orbiters[0].eccentricity
        a = self.orbiters[0].semimajor_axis
        m = self.orbiters[0].mass
        r = self.orbiters[0].radius
        M = self.central_particle[0].mass
        R = self.central_particle[0].radius
        O = self.central_particle[0].Omega
        Op = angular_frequency(M, m, a)
        dt =  abs(e/self.edot_star(a, M, R, m, r, O, Op, e))
        return dt

    def add_particles(self, p):
        self.orbiters.add_particles(p)

    def delete_particles(self, p):
        self.orbiters.removel_particles(p)

    def evolve_model(self, time):
        M = self.central_particle.mass
        R = self.central_particle.radius
        Os = self.central_particle.Omega

        interacting_bodies = self.orbiters.select(lambda a, e: a*(1-e)<self.pericenter_interaction_factor*R,["semimajor_axis", "eccentricity"])
        if len(interacting_bodies):
            print "N tidal:", len(self.orbiters), "of which N=", len(interacting_bodies), "tidally interacting" 
        self.orbiters_with_error = Particles()

        """
        #Interestingly enough, integrate.quadpack is not thread-safe.
        # AvanE and SPZ, 7 Jan 2014
        #
        threads = []
        nproc = 3
        Nper_proc = max(1, len(interacting_bodies)/nproc)
        print "npp=", Nper_proc
        for offset in range(0, len(interacting_bodies), Nper_proc):
            subset = interacting_bodies[offset:offset+Nper_proc].copy()
            print "lss=", len(subset)
            thread = threading.Thread(target=self.evolve_multiple_orbiters, args=[subset, time, M, R, Os])
            threads.append(thread)
        for ti in threads:
            ti.start()
        for ti in threads:
            ti.join()

        """
        for pi in interacting_bodies:
            self.evolve_individual_orbiter(pi, time, M, R, Os)


        print "Central_particle:", self.central_particle.Omega
        self.current_time = time
        if len(interacting_bodies):
            print "Post tidal interaction", len(interacting_bodies)
        if len(self.orbiters_with_error)>0:
            print "Error in N=", len(self.orbiters_with_error), "orbiters."
            print self.orbiters_with_error

        merged_orbiters = interacting_bodies.select(lambda a, e, r: a*(1-e)<R+r,["semimajor_axis", "eccentricity", "radius"])
        self.all_merged_orbiters = Particles()
        if len(merged_orbiters)>0:
            print "Merged orbiters N= ", len(merged_orbiters)
            print merged_orbiters
            self.all_merged_orbiters.add_particles(merged_orbiters-self.orbiters_with_error)

    def contains_nan(self, dO):
        has_nan = False
        for xi in dO:
            if math.isnan(xi):
                has_nan = True
        return has_nan

    def evolve_multiple_orbiters(self, interacting_bodies, time, M, R, Os):
        for pi in interacting_bodies:
            print "ev=", pi.key
            self.evolve_individual_orbiter(pi, time, M, R, Os)

    def evolve_individual_orbiter(self, oi, time, M, R, Os):
        dOmega = zero
        m = oi.mass
        r = oi.radius
        current_time = self.current_time
        while current_time<time:
            #print "Time=", oi.name, current_time, oi.semi_major_axis, oi.eccentricity, self.central_particle.Omega
            a = oi.semimajor_axis
            e = oi.eccentricity
            Op = angular_frequency(M, m, a)

            dt = time-current_time
            dt = min(dt, abs(e/self.edot_star(a, M, R, m, r, Os, Op, e)) )
            t_end = current_time + dt
            
            da = integrate.quad(lambda x: self.adot_star(a, M, R, m, r, Os, Op, e).value_in(units.RSun/units.s), current_time.value_in(units.s), t_end.value_in(units.s))
            de = integrate.quad(lambda x: self.edot_star(a, M, R, m, r, Os, Op, e).value_in(units.s**-1), current_time.value_in(units.s), t_end.value_in(units.s))
            dO = integrate.quad(lambda x: self.Omegadot_star(a, M, R, m, r, Os, Op, e).value_in(units.s**-2), current_time.value_in(units.s), t_end.value_in(units.s))
            if self.contains_nan(da):
                print "NAN's detected in da", da
                self.orbiters_with_error.add_particles(oi.as_set())
                break
            if self.contains_nan(de):
                print "NAN's detected, de", de
                self.orbiters_with_error.add_particles(oi.as_set())
                break
            if self.contains_nan(dO):
                print "NAN's detectedm dO", dO
                self.orbiters_with_error.add_particles(oi.as_set())
                break
            oi.semimajor_axis += da[0] | units.RSun
            oi.eccentricity += de[0]
            if oi.eccentricity<0:
                oi.eccentricity = 0
            dOmega += dO[0] | units.s**-1
            ######oi.age += dt
            current_time += dt
            #print "Time=", oi.name, current_time, oi.semi_major_axis, oi.eccentricity, self.central_particle.Omega
        self.central_particle.Omega += dOmega

#        self.central_particle.Omega = self.central_particle.Omega * (self.central_particle.old_radius/self.central_particle.radius)**2

    def J_star(self, Ms, Rs, Omega_s):
        k2s = self.central_particle[0].gyration_radius_sq
        Is = k2s*Ms*Rs*Rs
        return Is*Omega_s

    def J_planet(self, Mp, Rp, k2p, Omega_p):
        k2p = 0.2 # depends on the planet
        Ip = k2p*Mp*Rp*Rp
        return Ip*Omega_p

    # timescale for particle with mass Mb
    def tidal_timescale(self, Ma, Mb, Rb, a, sigma):
        denominator = (9*Ma*(Ma + Mb)*(Rb**10)*sigma)
        t_tidal = float("infinity") | units.s
        if not denominator==zero:
            t_tidal = Mb*(a**8) / denominator
        return t_tidal

    def Tp(self, Ms, Mp, Rp, a):
        return self.tidal_timescale(Ms, Mp, Rp, a, sigma_planet())

    def Ts(self, Ms, Mp, Rs, a):
        return self.tidal_timescale(Mp, Ms, Rs, a, sigma_star())

    def adot_planet(self, a, Ms, Rs, Mp, Rp, Omega_s, Omega_p, e):
        T_p = self.Tp(Ms, Mp, Rp, a)
        omega = angular_frequency(Ms, Mp, a)
        return -(a/T_p) * (self.f1(e) - Omega_p/omega * self.f2(e))

    def adot_star(self, a, Ms, Rs, Mp, Rp, Omega_s, Omega_p, e):
        T_s = self.Ts(Ms, Mp, Rs, a)
        omega = angular_frequency(Ms, Mp, a)
        adot = -(a/T_s) * (self.f1(e) - Omega_s/omega * self.f2(e))
        return adot

    def edot_planet(self, a, Ms, Rs, Mp, Rp, Omega_s, Omega_p, e):
        T_p =  self.Tp(Ms, Mp, Rp, a)
        omega = angular_frequency(Ms, Mp, a)
        return -9./2. * e/T_p * (self.f3(e) - 11./18. * Omega_p/omega * self.f4(e))

    def edot_star(self, a, Ms, Rs, Mp, Rp, Omega_s, Omega_p, e):
        T_s =  self.Ts(Ms, Mp, Rs, a)
        omega = angular_frequency(Ms, Mp, a)
        edot = -9./2. * e/T_s * (self.f3(e) - 11./18. * Omega_s/omega * self.f4(e))
        return edot

    def Omegadot_planet(self, a, Ms, Rs, Mp, Rp, Omega_s, Omega_p, e):
        T_p =  self.Tp(Ms, Mp, Rp, a)
        omega = angular_frequency(Ms, Mp, a)
        gamma = J_orb(Ms, Mp, a, e)/self.J_planet(Mp, Rp, Omega_p)
        return gamma/2. * Omega_p/T_p * (self.f5(e) - Omega_p/omega * self.f6(e))
    
    def Omegadot_star(self, a, Ms, Rs, Mp, Rp, Omega_s, Omega_p, e):
        T_s =  self.Ts(Ms, Mp, Rs, a)
        omega = angular_frequency(Ms, Mp, a)
        gamma = J_orb(Ms, Mp, a, e)/self.J_star(Ms, Rs, Omega_s)
        return gamma/2. * Omega_s/T_s * (self.f5(e) - Omega_s/omega * self.f6(e))

    def f1(self, e):
        res = 1
        if (e > 0):
            res = (1 + 31./2. *e*e + 255./8. *(e**4) + 185./16. *(e**6) + 25./64. *(e**8)) \
                /math.pow((1 - e*e), 7.5) 
        return res

    def f2(self, e):
        res = 1
        if (e > 0):
            res = (1 + 15./2. *e*e + 45./8. *math.pow(e,4) + 5./16. *math.pow(e, 6)) \
                /math.pow((1 - e*e), 6) 
        return res

    def f3(self, e):
        res = 1
        if (e > 0):
            res = (1 + 15./4. * e*e + 15./8. *math.pow(e,4) + 5./64. *math.pow(e, 6)) \
                /math.pow((1 - e*e), 6.5) 
        return res

    def f4(self, e):
        res = 1
        if (e > 0):
            res = (1 + 3./2. *e*e + 1./8. *math.pow(e,4)) /math.pow((1 - e*e), 5) 
        return res

    def f5(self, e):
        res = 1
        if (e > 0):
            res = (1 + 15./2. * e*e + 45./8. *math.pow(e,4) + 5./16. *math.pow(e, 6)) \
                /math.pow((1 - e*e), 6.5) 
        return res

    def f6(self, e):
        res = 1
        if (e > 0):
            res = (1 + 3 * e*e + 3./8. *math.pow(e,4)) /math.pow((1 - e*e), 5) 
        return res

#import unittest
from amuse.test.amusetest import TestCase
class TestTidalInteraction(TestCase):
    def test_remove_orbiters(self):
        M = 1 | units.MSun
        m = 1 | units.MJupiter
        a = 1 | units.AU
        e = 0.99
        Omega_s = 2.6e-6|units.s**-1

        star = Particles(1)
        star.mass = M
        star.Omega = Omega_s
        star.stellar_type = 1|units.stellar_type
        stellar = SeBa()
        stellar.particles.add_particles(star)

        channel_from_se_to_framework = stellar.particles.new_channel_to(star)
        channel_from_se_to_framework.copy_attributes(["age", "mass", "radius", "luminosity", "temperature", "stellar_type"])

        tidal = TidalEvolution(star)
        planet = Particles(2)
        planet.mass = [1, 100] * m
        planet.radius = [0.001, 0.001] |units.RSun
        planet.semimajor_axis = [1, 5.2] * a
        planet.eccentricity = [e, 0.1]
        tidal.add_particles(planet)
        channel_from_tc_to_framework = tidal.central_particle.new_channel_to(star)
        channel_from_to_to_framework = tidal.orbiters.new_channel_to(planet)
        channel_from_framework_to_tc = star.new_channel_to(tidal.central_particle)
        channel_from_framework_to_to = planet.new_channel_to(tidal.orbiters)

        dt = 1|units.Myr
        time = 0*dt
        He_WD = 10 | units.stellar_type
        while star.stellar_type<He_WD:
            print "T=", stellar.model_time
            stellar.particles.evolve_one_step()
            time = stellar.particles.age
            adiabatic_expansion_factor = star[0].mass/stellar.particles[0].mass
            tidal.central_particle.mass = stellar.particles[0].mass
            star[0].Omega = star[0].Omega * (star[0].radius/stellar.particles[0].radius)**2
            planet.semimajor_axis *= adiabatic_expansion_factor

            channel_from_se_to_framework.copy_attributes(["age", "mass", "radius", "luminosity", "temperature", "stellar_type"])

            channel_from_framework_to_tc.copy_attributes(["mass", "radius", "Omega"])
            channel_from_framework_to_to.copy_attributes(["semimajor_axis"])


            tidal.central_particle.Omega = star[0].Omega
            tidal.central_particle.radius = star[0].radius
            tidal.evolve_model(time)

            channel_from_to_to_framework.copy_attributes(["semimajor_axis", "eccentricity"])

            if len(tidal.orbiters_with_error)>0:
                print "Remove orbiter with error:", len(tidal.orbiters_with_error)
                tidal.orbiters.remove_particles(tidal.orbiters_with_error)
                planet.remove_particle(tidal.orbiters_with_error)

            if len(tidal.all_merged_orbiters)>0:
                print "Merged planets/asteroids: N=", len(tidal.all_merged_orbiters)
                print "removed:", tidal.all_merged_orbiters
                tidal.orbiters.remove_particles(tidal.all_merged_orbiters)
                planet.remove_particle(tidal.all_merged_orbiters)
            print "Remaining orbiters:", tidal.orbiters
            print "N=", len(tidal.orbiters)
        self.assertEquals(len(tidal.orbiters), 1)

def tidal_interaction(M, m, a, e, Omega_s, tend):

    star = Particles(1)
    star.mass = M
    star.Omega = Omega_s
    star.stellar_type = 1|units.stellar_type
    stellar = SeBa()
    stellar.particles.add_particles(star)

    channel_from_se_to_framework = stellar.particles.new_channel_to(star)
    channel_from_se_to_framework.copy_attributes(["age", "mass", "radius", "luminosity", "temperature", "stellar_type"])

    tidal = TidalEvolution(star)
    planet = Particles(1)
    planet.mass = 1*m
    planet.radius = 0.001 |units.RSun
    planet.semimajor_axis = 1.*a
    planet.eccentricity = e
    tidal.add_particles(planet)
    channel_from_tc_to_framework = tidal.central_particle.new_channel_to(star)
    channel_from_to_to_framework = tidal.orbiters.new_channel_to(planet)
    channel_from_framework_to_tc = star.new_channel_to(tidal.central_particle)
    channel_from_framework_to_to = planet.new_channel_to(tidal.orbiters)

#    bodies = ParticlesSuperset([star, planet])

#    dt = 1|units.Myr
    dt = tidal.orbital_evolution_time_scale()
    time = zero
    while time<tend:
        dt_se = stellar.particles[0].time_step
        dt = min(dt, dt_se)
        dt = max(1|units.Myr, dt)
        print "dt_tidal=", dt
        time += dt
        stellar.evolve_model(time)

        adiabatic_expansion_factor = star[0].mass/stellar.particles[0].mass
        tidal.central_particle.mass = stellar.particles[0].mass
        star[0].Omega = star[0].Omega * (star[0].radius/stellar.particles[0].radius)**2
        #expand planetary orbit due to stellar mass loss
        planet.semimajor_axis *= adiabatic_expansion_factor

        channel_from_se_to_framework.copy_attributes(["age", "mass", "radius", "luminosity", "temperature", "stellar_type"])

        channel_from_framework_to_tc.copy_attributes(["mass", "radius", "Omega"])
        channel_from_framework_to_to.copy_attributes(["semimajor_axis"])


        tidal.central_particle.Omega = star[0].Omega
        tidal.central_particle.radius = star[0].radius
        tidal.evolve_model(time)

        channel_from_to_to_framework.copy_attributes(["semimajor_axis", "eccentricity"])#, "merged_with_central_star"])

        if len(tidal.all_merged_orbiters)>0:
            print "Merged planets/asteroids: N=", len(tidal.all_merged_orbiters)
            print "removed:", tidal.all_merged_orbiters
            tidal.orbiters.remove_particles(tidal.all_merged_orbiters)
            planet.remove_particle(tidal.all_merged_orbiters)
        print "Remaining orbiters:", tidal.orbiters
        print "N=", len(tidal.orbiters)

        if len(tidal.orbiters)==0:
            return

    print "current time=", tidal.current_time, tidal.orbiters.semimajor_axis, tidal.orbiters.eccentricity

def new_option_parser():
#    from optparse import OptionParser
    from amuse.units.optparse import OptionParser
    result = OptionParser()
    result.add_option("-M", unit=units.MSun,
                      dest="M",  type="float",  default = 1|units.MSun,
                      help="stellar mass [%default]")
    result.add_option("-m", unit=units.MJupiter,
                      dest="m",  type="float",  default = 0.001|units.MJupiter,
                      help="planet mass [%default]")
    result.add_option("-a", unit=units.RSun, 
                      dest="a",  type="float",  default = 1|units.RSun,
                      help="planet semi major axis [%default]")
    result.add_option("-e", 
                      dest="e",  type="float", default = 0.6,
                      help="planet eccentricity [%default]")
    result.add_option("-t",  unit=units.Myr,
                      dest="tend",  type="float", default = 1|units.Myr,
                      help="end time of integration [%default]")
    result.add_option("-O", unit=units.s**-1,
                      dest="Omega_sun", type="float", default = 2.6e-6|units.s**-1,
                      help="Stellar angular something [%default]")
    return result

if __name__ in ('__main__', '__plot__'):
    set_printing_strategy("custom", #nbody_converter = converter, 
                          preferred_units = [units.MSun, units.RSun, units.Myr], 
                          precision = 11, prefix = "", 
                          separator = " [", suffix = "]")

    o, arguments  = new_option_parser().parse_args()
    tidal_interaction(o.M, o.m, o.a, o.e, o.Omega_sun, o.tend)

