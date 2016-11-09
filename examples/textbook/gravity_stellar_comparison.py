"""
   Nbody integration of N particles with a Salpeter initial mass
   function between Mmin and Mmax and with stellar evolution with
   metalicity z.
"""
from amuse.lab import *
from amuse.ext.LagrangianRadii import LagrangianRadii

class Gravity:
        def __init__(self, gravity_code, particles):

            self.particles = particles
            self.converter = nbody_system.nbody_to_si(1|units.MSun, 1|units.parsec)
            self.code = gravity_code(self.converter)
            self.code.particles.add_particles(self.particles)
            self.channel_to_framework = self.code.particles.new_channel_to(self.particles)
            self.channel_from_framework = self.particles.new_channel_to(self.code.particles)

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

def main(N=10, W0=7.0, t_end=10, Rvir=1, Mmin=0.1, Mmax= 100, z=0.02):
    t_end = t_end | units.Myr
    Rvir = Rvir | units.parsec
    Mmin = Mmin | units.MSun
    Mmax = Mmax | units.MSun

    import numpy
#    numpy.random.seed(1)

    masses = new_salpeter_mass_distribution(N, Mmin, Mmax)
    Mtot_init = masses.sum()
    converter=nbody_system.nbody_to_si(Mtot_init,Rvir)
    bodies = new_king_model(N, W0,convert_nbody=converter)
    bodies.mass = masses.copy()
    bodies.scale_to_standard(convert_nbody=converter)

#    gravity_only = Gravity(ph4, bodies)

#    stars = bodies.copy()

#    numpy.random.seed(1)
#    masses = new_salpeter_mass_distribution(N, Mmin, Mmax)
    stars = new_king_model(N, W0,convert_nbody=converter)
    stars.mass = masses.copy()
    stars.position = bodies.position.copy()
    stars.velocity = bodies.velocity.copy()
#    stars.scale_to_standard(convert_nbody=converter)
    
    stellar = SSE()
    stellar.parameters.metallicity = z
    stellar.particles.add_particle(stars)
    channel_from_stellar = stellar.particles.new_channel_to(stars)

#    gravity_stellar = Gravity(ph4, stars)

    converter = nbody_system.nbody_to_si(1|units.MSun, 1|units.parsec)
    gravity_only = ph4(converter)
    gravity_only.particles.add_particles(bodies)
    gravity_stars = ph4(converter)
    gravity_stars.particles.add_particles(stars)
    channel_to_gravity = stars.new_channel_to(gravity_stars.particles)

    
    Etot_init = gravity_only.kinetic_energy + gravity_only.potential_energy
    dE_gr = 0 | Etot_init.unit
    dE_se = 0 | Etot_init.unit
    time = 0.0 | t_end.unit
    dt = 0.1|units.Myr
    Lr_stars25 = [] |units.parsec
    Lr_bodies25 = [] |units.parsec
    Lr_stars50 = [] |units.parsec
    Lr_bodies50 = [] |units.parsec
    Lr_stars75 = [] |units.parsec
    Lr_bodies75 = [] |units.parsec
    t = [] | units.Myr
    t.append(time)
    Lr_stars25.append(LagrangianRadii(gravity_stars.particles)[5])
    Lr_bodies25.append(LagrangianRadii(gravity_only.particles)[5])
    Lr_stars50.append(LagrangianRadii(gravity_stars.particles)[6])
    Lr_bodies50.append(LagrangianRadii(gravity_only.particles)[6])
    Lr_stars75.append(LagrangianRadii(gravity_stars.particles)[7])
    Lr_bodies75.append(LagrangianRadii(gravity_only.particles)[7])
    while time < t_end:
        time = time+dt

        Etot_gr = gravity_only.kinetic_energy + gravity_only.potential_energy
        gravity_only.evolve_model(time)
        dE_gr += (gravity_only.kinetic_energy + gravity_only.potential_energy-Etot_gr)

        gravity_stars.evolve_model(time)

        stellar.evolve_model(time)
        Etot_se = gravity_only.kinetic_energy + gravity_only.potential_energy
        channel_from_stellar.copy_attributes(["mass"])
        channel_to_gravity.copy_attributes(["mass"])
        print stellar.particles.mass.sum(), gravity_only.particles.mass.sum(), gravity_stars.particles.mass.sum()


        t.append(time)
        Lr_stars25.append(LagrangianRadii(gravity_stars.particles)[5])
        Lr_bodies25.append(LagrangianRadii(gravity_only.particles)[5])
        Lr_stars50.append(LagrangianRadii(gravity_stars.particles)[6])
        Lr_bodies50.append(LagrangianRadii(gravity_only.particles)[6])
        Lr_stars75.append(LagrangianRadii(gravity_stars.particles)[7])
        Lr_bodies75.append(LagrangianRadii(gravity_only.particles)[7])
        #print "Time=", time, Lr_stars, Lr_bodies

    Ekin = gravity_only.kinetic_energy 
    Epot = gravity_only.potential_energy
    Etot = Ekin + Epot
    dE = Etot_init-Etot
    Mtot = bodies.mass.sum()
    print "T=", time, 
    print "M=", Mtot, "(dM[SE]=", Mtot/Mtot_init, ")",
    print "E= ", Etot, "Q= ", Ekin/Epot,
    print "dE/E=", (Etot_init-Etot)/Etot,
    print "(dE[gr]/E=", dE_gr/Etot, ",", 
    print "dE[se]/E=", (Etot_init-Etot-dE_gr)/Etot, ")"
    Etot_init -= dE

    gravity_only.stop()
    gravity_stars.stop()
    stellar.stop()
    return t, Lr_stars25, Lr_bodies25, Lr_stars50, Lr_bodies50, Lr_stars75, Lr_bodies75
    
def new_option_parser():
    from amuse.units.optparse import OptionParser
    result = OptionParser()
    result.add_option("-N", dest="N", type="int",default = 1000,
                      help="number of stars [100]")
    result.add_option("-M", dest="Mmax", type="float",default = 100,
                      help="maximal stellar mass [100] MSun")
    result.add_option("-m", dest="Mmin", type="float",default = 1,
                      help="minimal stellar mass [0.1] MSun")
    result.add_option("-R", dest="Rvir", type="float",default = 3.0,
                      help="cluser virial radius [1] in parsec")
    result.add_option("-t", dest="t_end", type="float", default = 10.0,
                      help="end time of the simulation [1] Myr")
    result.add_option("-W", dest="W0", type="float", default = 3.0,
                      help="Dimension-less depth of the King potential (W0) [7.0]")
    result.add_option("-z", dest="z", type="float", default = 0.02,
                      help="metalicity [0.02]")
    return result

if __name__ in ('__main__', '__plot__'):
    o, arguments  = new_option_parser().parse_args()
    time, Lr_stars25, Lr_bodies25, Lr_stars50, Lr_bodies50, Lr_stars75, Lr_bodies75 = main(**o.__dict__)
    from matplotlib import pyplot
    from prepare_figure import single_frame
    from distinct_colours import get_distinct

    x_label = "t [Myr]"
    y_label = "R [pc]"
    fig = single_frame(x_label, y_label, logx=False, logy=False, xsize=14, ysize=12)
    color = get_distinct(4)

#    import numpy
#    Lr_stars = numpy.matrix.transpose(numpy.matrix(Lr_stars))
#    Lr_bodies = numpy.matrix.transpose(numpy.matrix(Lr_bodies))
    
#    print Lr_stars[6], Lr_bodies[6]
    pyplot.plot(time.value_in(units.Myr), Lr_stars25.value_in(units.parsec), c=color[0], label= 'with mass loss')
    pyplot.plot(time.value_in(units.Myr), Lr_bodies25.value_in(units.parsec), c=color[1], label= 'without mass loss')

    pyplot.plot(time.value_in(units.Myr), Lr_stars50.value_in(units.parsec), c=color[0])
    pyplot.plot(time.value_in(units.Myr), Lr_bodies50.value_in(units.parsec), c=color[1])

    pyplot.plot(time.value_in(units.Myr), Lr_stars75.value_in(units.parsec), c=color[0])
    pyplot.plot(time.value_in(units.Myr), Lr_bodies75.value_in(units.parsec), c=color[1])
    
    pyplot.legend(loc="upper left", ncol=1, shadow=False, fontsize=24)

    #pyplot.show()
    pyplot.savefig("gravity_stellar_comparison")

