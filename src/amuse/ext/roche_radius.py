"""
    Various equations to calculate the size of the Roche lobe of a system.

    There are three ways to use this module:
    1. From the command line:
        python roche_radius.py -a 10 -e 0.5 -m 10

    2. Using the direct functions:
        q = 10
        print eggleton_formula(q)

        e = 0.5
        A = sepinsky_A_parameter(e)
        print sepinsky_formula(q, A)

    3. Using the Roche_Orbit class:
        orbit = Roche_Orbit(mass_2=10|units.MSun, eccentricity=0.5)

        print orbit.eggleton_roche_radius()
        print orbit.sepinsky_roche_radius()

    The functions also accept arrays instead of single numbers:
        orbit.mass_1 = [1.0, 2.0, 3.0] | units.MSun

        print orbit.sepinsky_roche_radius()
"""
import numpy
import math
from amuse.units.optparse import OptionParser
from amuse.units import units, constants, quantities
from amuse.support.console import set_printing_strategy

""" Equation 47-52 of Sepinsky """

def low_q_low_A(q, A):
    log_q = numpy.log10(q)

    return 1 + 0.11 * (1-A) - 0.05 * (1-A) * numpy.exp(-(0.5 * (1 + A) + log_q)**2)

def high_q_low_A(q, A):
    log_q = numpy.log10(q)

    return 1.226 - 0.21 * A - 0.15 * (1-A) * numpy.exp((0.25 * A - 0.3) * (log_q)**1.55)

def low_q_medium_A(q, A):
    log_q = numpy.log10(q)
    log_A = numpy.log10(A)

    g_0 = 0.9978 - 0.1229 * log_A - 0.1273 * (log_A)**2
    g_1 = 0.001 + 0.02556 * log_A
    g_2 = 0.0004 + 0.0021 * log_A
    return g_0 + g_1 * log_q * g_2 * (log_q)**2

def high_q_medium_A(q, A):
    log_q = numpy.log10(q)
    log_A = numpy.log10(A)

    h_0 = 1.0071 -  0.0907 * log_A - 0.0495 * (log_A)**2
    h_1 = -0.004 -  0.163 * log_A - 0.214 * (log_A)**2
    h_2 = 0.00022 -  0.0108 * log_A - 0.02718 * (log_A)**2
    return h_0 + h_1 * log_q * h_2 * (log_q)**2

def low_q_high_A(q, A):
    log_q = numpy.log10(q)
    log_A = numpy.log10(A)

    num_0 = 6.3014 * (log_A)**1.3643
    den_0 = numpy.exp(2.3644 * (log_A)**0.70748) - 1.4413 * numpy.exp(-0.0000184 * (log_A)**-4.5693)
    i_0 = num_0 / den_0

    den_1 = 0.0015 * numpy.exp(8.84 * (log_A)**0.282) + 15.78
    i_1 = log_A / den_1

    num_2 = 1 + 0.036 * numpy.exp(8.01 * (log_A)**0.879)
    den_2 = 0.105 * numpy.exp(7.91 * (log_A)**0.879)
    i_2 = num_2 / den_2

    den_3 = 1.38 * numpy.exp(-0.035 * (log_A)**0.76) + 23.0 * numpy.exp(-2.89 * (log_A)**0.76)
    i_3 = 0.991 / den_3

    return i_0 + i_1 * numpy.exp(-i_2 * (log_q + i_3)**2)

def high_q_high_A(q, A):
    log_q = numpy.log10(q)
    log_A = numpy.log10(A)

    num_0 = 1.895 * (log_A)**0.837
    den_0 = numpy.exp(1.636 * (log_A)**0.789) - 1
    j_0 = num_0 / den_0

    num_1 = 4.3 * (log_A)**0.98
    den_1 = numpy.exp(2.5 * (log_A)**0.66) + 4.7
    j_1 = num_1 / den_1

    den_2 = 8.8 * numpy.exp(-2.95 * (log_A)**0.76) + 1.64 * numpy.exp(-0.03 * (log_A)**0.76)
    j_2 = 1.0 / den_2

    j_3 = 0.256 * numpy.exp(-1.33 * (log_A)**2.9) * (5.5 * numpy.exp(1.33 * (log_A)**2.9) + 1)

    return j_0 + j_1 * numpy.exp(-j_2 * (log_q)**j_3)

functions = {'low': {'low':low_q_low_A, 'medium':low_q_medium_A, 'high':low_q_high_A },
             'high': {'low':high_q_low_A, 'medium':high_q_medium_A, 'high':high_q_high_A }}

def pick_formula(q, A):
    log_q = numpy.log10(q)
    log_A = numpy.log10(A)

    q_regime = 'low' if log_q <= 0 else 'high'
    A_regime = 'low' if log_A <= -0.1 else 'medium' if -0.1 < log_A <= 0.2 else 'high'

    function = functions[q_regime][A_regime]

    return function(q, A)

vec_pick_formula = numpy.vectorize(pick_formula)

def sepinsky_formula(q=1, A=1):
    """ The correction to the Eggleton Roche radius following Sepinsky, Willems and Kalogera 2007 """
    return vec_pick_formula(q, A)

def sepinsky_A_parameter(eccentricity=0.0, angular_velocity_ratio=1.0, true_anomaly=numpy.pi):
    numerator = angular_velocity_ratio**2 * (1.0 + eccentricity)**4
    denominator = (1.0 + eccentricity * numpy.cos(true_anomaly))**3
    return numerator / denominator

def eggleton_formula(mass_ratio):
    """ Use the Eggleton formula for the given mass ratio.

        The Eggleton formula assumes a circular corotating system, which means:
        eccentricity = 0
        true_anomaly = pi
        angular_velocity_ratio = 1
    """

    two_third = mass_ratio**(2.0/3.0)
    one_third = mass_ratio**(1.0/3.0)
    return 0.49 * two_third / ( 0.6 * two_third + numpy.log(1.0 + one_third))

def separation(semimajor_axis, eccentricity, true_anomaly):
    """ Return the orbital separation in the same units as the semimajor axis"""
    numerator = semimajor_axis * (1.0-eccentricity**2)
    denominator = 1.0 + eccentricity * numpy.cos(true_anomaly)
    return numerator / denominator

class Roche_Orbit(object):
    """
        A set of orbital parameters that allows the calculation of the Roche radius.

        See:
        Eggleton 1983
        Sepinsky, Willems and Kalogera 2007
    """

    def __init__(self, mass_1=1.0|units.MSun, mass_2=1.0|units.MSun, eccentricity=0.0, true_anomaly=0.0, angular_velocity_ratio=1.0, semimajor_axis=1.0|units.RSun):
        self.mass_1 =  mass_1
        self.mass_2 = mass_2
        self.eccentricity = eccentricity
        self.true_anomaly = true_anomaly
        self.angular_velocity_ratio = angular_velocity_ratio
        self.semimajor_axis = semimajor_axis

    @property
    def mass_ratio(self):
        return self.mass_1 / self.mass_2

    @mass_ratio.setter
    def mass_ratio(self, value):
        self.mass_1 = value * self.mass_2

    @property
    def total_mass(self):
        return self.mass_1 + self.mass_2

    @property
    def A(self):
        return sepinsky_A_parameter(self.eccentricity, self.angular_velocity_ratio, self.true_anomaly)

    @property
    def period(self):
        return (4.0 * numpy.pi**2 * self.semimajor_axis**3 / (constants.G * self.total_mass)).sqrt()

    @period.setter
    def period(self, period):
        self.semimajor_axis = (period/(2.0 * numpy.pi) * (self.total_mass * constants.G).sqrt())**(2.0/3.0)

    def sepinsky_over_eggleton(self):
        return sepinsky_formula(self.mass_ratio, self.A)

    def eggleton_roche_over_separation(self):
        return eggleton_formula(self.mass_ratio)

    def eggleton_roche_radius(self):
        """ The Roche radius assumes a curcular orbit with the current separation.

            Note that this is not really correct for non-circular orbits.
        """
        return self.eggleton_roche_over_separation() * self.separation()

    def sepinsky_roche_radius(self):
        return self.sepinsky_over_eggleton() * self.eggleton_roche_radius()

    def separation(self):
        return separation(self.semimajor_axis, self.eccentricity, self.true_anomaly)

def new_option_parser():
    parser = OptionParser(description="Calculate the Roche radius for a given orbit.")
    parser.add_option("-a", dest="semimajor_axis", type="float", default = 1.0, unit=units.AU, help="The orbit semimajor axis [%default %unit]")
    parser.add_option("-p", dest="period", type="float", default = 0.0, unit=units.day, help="The orbital period, which sets the semimajor axis [%unit]")
    parser.add_option("-e", dest="eccentricity", type="float", default = 0.0, help="The orbit eccentricity [%default]")
    parser.add_option("-M", dest="mass_1", type="float", default = 1.0, unit=units.MSun, help="The mass of the primary (the object that has the Roche lobe) [%default %unit]")
    parser.add_option("-m", dest="mass_2", type="float", default = 1.0, unit=units.MSun, help="The mass of the secondary (the object which causes the Roche lobe) [%default %unit]")
    parser.add_option("-n", dest="true_anomaly", type="float", default = 0.0, help="The true anomaly, the angle between the objects location and it's periastron location [%default]")
    parser.add_option("-f", dest="angular_velocity_ratio", type="float", default = 1.0, help="The rotational angular velocity of object 1, in units of the orbital angular velocity at periastron [%default]")
    return parser

def create_orbit_from_options():
    options, args  = new_option_parser().parse_args()
    options = options.__dict__

    period = options.pop("period")

    orbit = Roche_Orbit(**options)

    if period > quantities.zero:
        orbit.period = period

    return orbit

def print_results(orbit):
    set_printing_strategy("custom", preferred_units = [units.MSun, units.RSun, units.Myr], precision = 7)

    if orbit.A == 1.0:
        print "This is a circular, corotating orbit, so the eggleton formula is correct."
    else:
        print "Warning: This is not a circular, corotating orbit, so the eggleton formula is not correct."

    print "Roche radius for: M =", orbit.mass_1, "m =", orbit.mass_2, "a =", orbit.semimajor_axis, "e =", orbit.eccentricity
    print

    print "Eggleton Roche radius =", orbit.eggleton_roche_radius()
    print "Sepinsky Roche radius =", orbit.sepinsky_roche_radius()

if __name__ == '__main__':

    orbit = create_orbit_from_options()
    print_results(orbit)

