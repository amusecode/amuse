"""
   Nbody integration of N particles with a Salpeter initial mass
   function between Mmin and Mmax and with stellar evolution with
   metallicity z.
"""
import numpy
from matplotlib import pyplot
from amuse.lab import *
from amuse.ext.LagrangianRadii import LagrangianRadii
from prepare_figure import single_frame
from distinct_colours import get_distinct

class Gravity:
        def __init__(self, gravity_code, particles):

            self.particles = particles
            self.converter = nbody_system.nbody_to_si(1|units.MSun,
                                                      1|units.parsec)
            self.code = gravity_code(self.converter)
            self.code.particles.add_particles(self.particles)
            self.channel_to_framework \
                    = self.code.particles.new_channel_to(self.particles)
            self.channel_from_framework \
                    = self.particles.new_channel_to(self.code.particles)

        def evolve_model(self, time):
            self.channel_from_framework
            self.code.evolve_model(time)
            self.channel_to_framework

        @property
        def model_time(self):            
            return self.code.model_time
        @property
        def particles(self):
            return self.code.particles
        @property
        def kinetic_energy(self):
            return self.code.kinetic_energy
        @property
        def potential_energy(self):
            return self.code.potential_energy
        @property
        def stop(self):
            return self.code.stop

def generate_initial_conditions(N, W0, Rvir, Mmin, Mmax):

    numpy.random.seed(123)

    masses = new_salpeter_mass_distribution(N, Mmin, Mmax)
    Mtot_init = masses.sum()
    converter=nbody_system.nbody_to_si(Mtot_init,Rvir)
    bodies = new_king_model(N, W0,convert_nbody=converter)
    bodies.mass = masses.copy()
    bodies.scale_to_standard(convert_nbody=converter)

    return bodies

def run_only_gravity(bodies, t_end):

    Mtot_init = bodies.mass.sum()
    time = [] | units.Myr
    Lr25 = [] |units.parsec
    Lr50 = [] |units.parsec
    Lr75 = [] |units.parsec
    converter = nbody_system.nbody_to_si(bodies.mass.sum(), 1|units.parsec)
    gravity = ph4(converter)
    gravity.particles.add_particles(bodies)
    channel_from_gravity = gravity.particles.new_channel_to(bodies)

    dt = 0.1|units.Myr
    while True:
        time.append(gravity.model_time)
        Lr25.append(LagrangianRadii(gravity.particles)[5])
        Lr50.append(LagrangianRadii(gravity.particles)[6])
        Lr75.append(LagrangianRadii(gravity.particles)[7])

        gravity.evolve_model(time[-1]+dt)
        channel_from_gravity.copy_attributes(["x", "y", "z", "vx", "vy", "vz"])

        print("G: T=", time[-1], "M=", bodies.mass.sum(), \
              "(dM[SE]=", bodies.mass.sum()/Mtot_init, ")")
        if time[-1] >= t_end:
            break
    gravity.stop()
    return time, Lr25, Lr50, Lr75

def run_sequential_gravity_and_stellar(bodies, t_end):
    Mtot_init = bodies.mass.sum()
    time = [] | units.Myr
    Lr25 = [] |units.parsec
    Lr50 = [] |units.parsec
    Lr75 = [] |units.parsec
    converter = nbody_system.nbody_to_si(bodies.mass.sum(), 1|units.parsec)
    gravity = ph4(converter)
    gravity.particles.add_particles(bodies)
    channel_from_gravity = gravity.particles.new_channel_to(bodies)
    channel_to_gravity = bodies.new_channel_to(gravity.particles)

    stellar = SSE()
    stellar.parameters.metallicity = 0.02
    stellar.particles.add_particle(bodies)
    channel_from_stellar = stellar.particles.new_channel_to(bodies)
    
    dt = 0.1|units.Myr
    while  True:
            
        stellar.evolve_model(gravity.model_time + dt/2)
        channel_from_stellar.copy_attributes(["mass"])
        channel_to_gravity.copy_attributes(["mass"])

        gravity.evolve_model(gravity.model_time+dt)
        channel_from_gravity.copy_attributes(["x", "y", "z", "vx", "vy", "vz"])

        stellar.evolve_model(gravity.model_time + dt)
        channel_from_stellar.copy_attributes(["mass"])
        channel_to_gravity.copy_attributes(["mass"])

        time.append(gravity.model_time)
        Lr25.append(LagrangianRadii(gravity.particles)[5])
        Lr50.append(LagrangianRadii(gravity.particles)[6])
        Lr75.append(LagrangianRadii(gravity.particles)[7])
        
        print("GS: T=", time[-1], "M=", bodies.mass.sum(), \
              "(dM[SE]=", bodies.mass.sum()/Mtot_init, ")")
        if time[-1] >= t_end:
            break

    gravity.stop()
    stellar.stop()
    return time, Lr25, Lr50, Lr75

def run_event_driven_gravity_and_stellar(bodies, t_end):

    Mtot_init = bodies.mass.sum()
    time = [] | units.Myr
    Lr25 = [] |units.parsec
    Lr50 = [] |units.parsec
    Lr75 = [] |units.parsec
    converter = nbody_system.nbody_to_si(bodies.mass.sum(), 1|units.parsec)
    gravity = ph4(converter, number_of_workers=2)
    gravity.particles.add_particles(bodies)

    stellar = SSE()
    stellar.parameters.metallicity = 0.02
    stellar.particles.add_particle(bodies)

    channel_from_gravity = gravity.particles.new_channel_to(
            bodies, attributes=["mass", "x", "y", "z", "vx", "vy", "vz"])
    channel_from_stellar = stellar.particles.new_channel_to(
            bodies, attributes=["mass"])
    channel_from_stellar_to_gravity = stellar.particles.new_channel_to(
            gravity.particles, attributes=["mass"])
    
    while True:

        dt = 0.5*stellar.particles.time_step.min()
        stellar.evolve_model(gravity.model_time + dt/2)
        channel_from_stellar_to_gravity.copy()

        dt = 0.5*stellar.particles.time_step.min()
        gravity.evolve_model(stellar.model_time + dt)
        channel_from_gravity.copy()

        stellar.evolve_model(gravity.model_time)
        channel_from_stellar.copy()
        
        time.append(gravity.model_time)
        Lr25.append(LagrangianRadii(gravity.particles)[5])
        Lr50.append(LagrangianRadii(gravity.particles)[6])
        Lr75.append(LagrangianRadii(gravity.particles)[7])

        stellar.evolve_model()
        
        print("GSE: T=", time[-1], "M=", bodies.mass.sum(), \
              "(dM[SE]=", bodies.mass.sum()/Mtot_init, ")")

        if time[-1] >= t_end:
            break

    gravity.stop()
    stellar.stop()
    return time, Lr25, Lr50, Lr75

def main(N, W0, t_end, Rvir, Mmin, Mmax):
    bodies = generate_initial_conditions(N, W0, Rvir, Mmin, Mmax)
    print(numpy.sort(bodies.mass.value_in(units.MSun)))

    x_label = "t [Myr]"
    y_label = "R [pc]"
    fig = single_frame(x_label, y_label, logx=False, logy=False,
                       xsize=14, ysize=12)
    color = get_distinct(4)

    time, Lr25, Lr50, Lr75 = run_only_gravity(bodies.copy(), t_end)
    pyplot.plot(time.value_in(units.Myr), Lr25.value_in(units.parsec),
                c=color[0], label= 'without mass loss')
    pyplot.plot(time.value_in(units.Myr), Lr50.value_in(units.parsec),
                c=color[0])
    pyplot.plot(time.value_in(units.Myr), Lr75.value_in(units.parsec),
                c=color[0])

    time, Lr25, Lr50, Lr75 \
            = run_sequential_gravity_and_stellar(bodies.copy(), t_end)
    pyplot.plot(time.value_in(units.Myr), Lr25.value_in(units.parsec),
                c=color[1], label= 'with mass loss')
    pyplot.plot(time.value_in(units.Myr), Lr50.value_in(units.parsec),
                c=color[1])
    pyplot.plot(time.value_in(units.Myr), Lr75.value_in(units.parsec),
                c=color[1])
    
    time, Lr25, Lr50, Lr75 \
            = run_event_driven_gravity_and_stellar(bodies.copy(), t_end)
    pyplot.plot(time.value_in(units.Myr), Lr25.value_in(units.parsec),
                c=color[2], label= 'event driven')
    pyplot.plot(time.value_in(units.Myr), Lr50.value_in(units.parsec),
                c=color[2])
    pyplot.plot(time.value_in(units.Myr), Lr75.value_in(units.parsec),
                c=color[2])

    pyplot.legend(loc="upper left", ncol=1, shadow=False, fontsize=24)

    save_file = 'gravity_stellar_comparison.png'
    pyplot.savefig(save_file)
    print('\nSaved figure in file', save_file,'\n')
    pyplot.show()
    
def new_option_parser():
    from amuse.units.optparse import OptionParser
    result = OptionParser()
    result.add_option("-N", dest="N", type="int",default = 1000,
                      help="number of stars [1000]")
    result.add_option("-M", unit=units.MSun, dest="Mmax", type="float",
                      default = 100|units.MSun,
                      help="maximun stellar mass [100] MSun")
    result.add_option("-m", unit=units.MSun, dest="Mmin", type="float",
                      default = 1|units.MSun,
                      help="minimum stellar mass [1] MSun")
    result.add_option("-R", unit=units.parsec, dest="Rvir", type="float",
                      default = 3.0|units.parsec,
                      help="cluster virial radius [3] in parsec")
    result.add_option("-t", unit=units.Myr,
                      dest="t_end", type="float", default = 10.0|units.Myr,
                      help="end time of the simulation [10] Myr")
    result.add_option("-W", dest="W0", type="float", default = 3.0,
                      help="dimensionless depth of the King potential (W0) [%default]")
    return result

if __name__ in ('__main__', '__plot__'):
    o, arguments  = new_option_parser().parse_args()
    main(**o.__dict__)


