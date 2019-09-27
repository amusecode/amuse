from amuse.lab import *

def new_option_parser():
    from amuse.units.optparse import OptionParser
    result = OptionParser()
    result.add_option("-N", dest="Nsph", type="int", default = 20000,
                      help="number of gas particles [%default]")
    result.add_option("-M", unit=units.MSun,
                      dest="mass", type="float", default = 5.5|units.MSun,
                      help="mass of the star [%default]")
    result.add_option("-R", unit=units.AU,
                      dest="radius", type="float", default = 0.42|units.AU,
                      help="Size of star [%default]")
    return result

if __name__ in ('__main__', '__plot__'):
    o, arguments = new_option_parser().parse_args()

###BOOKLISTSTART1###
    stellar = MESA()
    stellar.particles.add_particle(Particle(mass = o.mass))
    XiTau = stellar.particles[0]
    while XiTau.radius<=o.radius:
        print "Time=", stellar.model_time.in_(units.Myr), stellar.particles.radius.in_(units.AU), stellar.particles.stellar_type
        stellar.evolve_model()

    target_core_mass = XiTau.core_mass
    print XiTau
    print "Core mass:", XiTau.core_mass.in_(units.MSun)
    sph_model = convert_stellar_model_to_SPH(
        XiTau, 
        o.Nsph,
        with_core_particle = True,
        target_core_mass  = target_core_mass,
        do_store_composition = False,
        base_grid_options=dict(type="fcc")
    )
    core_particle = sph_model.core_particle.as_set()
    gas_particles = sph_model.gas_particles

    print "Ngas=", len(gas_particles), "Ncore=", core_particle
    write_set_to_file(core_particle, "Hydro_PrimaryStar_XiTau.amuse", "amuse", append_to_file=False)
    write_set_to_file(gas_particles, "Hydro_PrimaryStar_XiTau.amuse", "amuse")
    stellar.stop()
###BOOKLISTSTOP1###


