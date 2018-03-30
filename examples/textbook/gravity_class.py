import os.path
import math
import numpy
from amuse.lab import *
from amuse.ext.protodisk import ProtoPlanetaryDisk
from amuse.ext.evrard_test import uniform_unit_sphere
import time

from initialize_sstars import *
from calculate_orbital_parameters import calculate_orbital_parameters

from amuse.community.mercury.interface import Mercury
from amuse.community.mikkola.interface import Mikkola
from amuse.units.optparse import OptionParser

from amuse.community.adaptb.interface import Adaptb

class Gravity:
    def __init__(self, gravity_code, particles):
        self.converter=nbody_system.nbody_to_si(1|units.MSun,1|units.AU)
        self.code = gravity_code(self.converter)
        self.code.initialize_code()
        self.get_gravity_at_point = self.code.get_gravity_at_point
        self.get_potential_at_point = self.code.get_potential_at_point
        if hasattr(self.code.parameters, "bs_tolerance"):
            self.code.parameters.bs_tolerance  = 1.0e-10
            self.code.parameters.word_length = 512
        if hasattr(self.code.parameters, "lightspeed"):
            self.code.parameters.lightspeed = 0 | units.kms

        for pi in particles:
            self.code.particles.add_particles(pi)
        if isinstance(self.code, Mercury):
            self.code.commit_particles() 
    
        self.channel_to_framework = []
        self.channel_from_framework = []
        for pi in particles:
            self.channel_to_framework.append(self.code.particles.new_channel_to(pi))
            self.channel_from_framework.append(pi.new_channel_to(self.code.particles))

    @property
    def model_time(self):
        return self.code.model_time
    @property
    def stop(self):
        return self.code.stop
    @property
    def particles(self):
        return self.code.particles

    def evolve_model(self, model_time):
        for gi in self.channel_from_framework:
            gi.copy()
            #gi.copy_attributes("mass")
        self.code.evolve_model(model_time)
        for gi in self.channel_to_framework:
            gi.copy()

def print_diagnostics(model_time, particles, converter):
    bh = particles[0]
    stars = particles[1]
    for si in stars:
        a, e, m0, m1 = calculate_orbital_parameters(bh[0], si, converter)
        print "Orbital elements for star ", si.name, model_time.in_(units.yr), a.in_(units.AU), e, bh[0].mass.in_(units.MSun), si.mass.in_(units.MSun)

def main(t_end=1687|units.yr, n_steps=1, filename=None):
    black_hole, stars = initialize_sstars(2012|units.yr, S_name, S_a_arcsec, S_ecc, S_inc, S_omra, S_Omega, S_tperi, S_Period)

    gravity = Gravity(Mercury, [black_hole, stars])
#    gravity = Gravity(Mikkola, [black_hole, stars])
    
    print_diagnostics(gravity.model_time, [black_hole, stars], gravity.converter) 
    if filename:
        write_set_to_file(gravity.particles, filename, "hdf5")

    dt = t_end/float(n_steps)
    while gravity.model_time < t_end:
        gravity.evolve_model(gravity.model_time+dt)
        if filename:
            write_set_to_file(gravity.particles, filename, 'hdf5')
        print_diagnostics(gravity.model_time, [black_hole, stars], gravity.converter) 
    gravity.stop()

def new_option_parser():
    result = OptionParser()
    result.add_option("-n", dest="n_steps", type="int", default = 1,
                      help="number of diagnostics time steps [10]")
    result.add_option("-f", dest="filename", default = None,
                      help="write output filename")
    result.add_option("-t", unit=units.Myr, 
                      dest="t_end", type="float", default = 0.000001|units.Myr,
                      help="end time of the simulation [0.0000001] %unit")
    return result

if __name__ in ('__main__', '__plot__'):
    o, arguments  = new_option_parser().parse_args()
    main(**o.__dict__)
