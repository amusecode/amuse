from amuse.lab import *

def hydro_sink_particles(sinks, bodies):

    all_lost_particles = Particles()
    for s in sinks:
        xs,ys,zs=s.x,s.y,s.z
        radius_squared = s.radius**2
        insink=bodies.select_array(lambda x,y,z: (x-xs)**2+(y-ys)**2+(z-zs)**2 < radius_squared,['x','y','z'])  
        print "Ninsink=", len(insink)
        if len(insink)==0:
            return insink

        cm=s.position*s.mass
        print "sink=", s.position, s.mass, cm
        p=s.velocity*s.mass
        s.mass+=insink.total_mass()
        s.position=(cm+insink.center_of_mass()*insink.total_mass())/s.mass
        s.velocity=(p+insink.total_momentum())/s.mass
        print "insinkstar=", s.mass, s.position, s.velocity
        all_lost_particles.add_particles(insink)
    return all_lost_particles

def main(N, Mtot, Rvir, rsink):
    converter=nbody_system.nbody_to_si(Mtot, Rvir)
    bodies = new_plummer_gas_model(N, convert_nbody=converter)

    sink = Particles(1)
    sink.mass = 0 | units.MSun
    sink.radius = rsink
    sink.position = (0, 0, 0) | units.AU
    sink.velocity = (0, 0, 0) | units.kms

    accreted = hydro_sink_particles(sink, bodies)
    print "Particles in sink: N=", len(accreted), " M=", sink.mass
    print "sink position=", sink.position.as_quantity_in(units.AU)
    print "sink velocity=", sink.velocity.as_quantity_in(units.kms)
    
def new_option_parser():
    from amuse.units.optparse import OptionParser
    result = OptionParser()
    result.add_option("-N", dest="N", type="int",default = 100,
                      help="number of sph particles [100]")
    result.add_option("-M", unit=units.MSun,
                      dest="Mtot", type="float", default = 1|units.MSun,
                      help="Mass of molcular cloud [%default]")
    result.add_option("-R", unit=units.AU,
                      dest="Rvir", type="float", default = 100|units.AU,
                      help="Radius of cloud [%default]")
    result.add_option("-r", unit=units.AU,
                      dest="rsink", type="float", default = 100|units.AU,
                      help="Radius of the sink [%default]")
    return result

if __name__ in ('__main__', '__plot__'):
    o, arguments  = new_option_parser().parse_args()
    main(**o.__dict__)

