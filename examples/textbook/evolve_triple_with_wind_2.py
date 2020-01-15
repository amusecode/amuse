"""
   N-body integration of N partcles with a Salpeter initial mass
   function between Mmin and Mmax and with stellar evolution with
   metallicity z.
"""
import sys
import math, numpy
from optparse import OptionParser

import matplotlib
#matplotlib.use('Agg')
from matplotlib import pyplot as plt
from prepare_figure import single_frame
from distinct_colours import get_distinct

from amuse.units import units, constants, nbody_system
from amuse.units.quantities import zero
from amuse.datamodel import Particle, Particles
from amuse.support.console import set_printing_strategy
from amuse.io import store

from amuse.ext.orbital_elements import new_binary_from_orbital_elements
from amuse.ext.orbital_elements import orbital_elements_from_binary

from amuse.community.symple.interface import symple	# symplectic
from amuse.community.huayno.interface import Huayno	# symplectic
from amuse.community.smalln.interface import SmallN	# time reversible
from amuse.community.hermite.interface import Hermite	# not symplectic
from amuse.community.seba.interface import SeBa
from amuse.community.sse.interface import SSE

def orbital_period(a, Mtot):
    return 2*numpy.pi*(a**3/(constants.G*Mtot)).sqrt()

def semimajor_axis(P, Mtot):
    return (constants.G*Mtot*P**2/(4*numpy.pi**2))**(1./3)
    
def get_orbital_elements_of_triple(stars):
    inner_binary = stars[0]+stars[1]
    outer_binary = Particles(1)
    outer_binary[0].mass = inner_binary.mass.sum()
    outer_binary[0].position = inner_binary.center_of_mass()
    outer_binary[0].velocity = inner_binary.center_of_mass_velocity()
    outer_binary.add_particle(stars[2])
    M1, M2, ain, ein, ta_in, inc_in, lan_in, aop_in \
        = orbital_elements_from_binary(inner_binary, G=constants.G)
    M12, M3, aout, eout, ta_out, outc_out, lan_out, aop_out \
        = orbital_elements_from_binary(outer_binary, G=constants.G)
    return ain, ein, aout, eout
    
def evolve_triple_with_wind(M1, M2, M3, Pora, Pin_0, ain_0, aout_0,
                            ein_0, eout_0, t_end, nsteps, scheme, integrator,
                            t_stellar, dt_se, dtse_fac, interp):

    import random
    from amuse.ext.solarsystem import get_position

    numpy.random.seed(42)

    print("Initial masses:", M1, M2, M3)
    triple = Particles(3)
    triple[0].mass = M1
    triple[1].mass = M2
    triple[2].mass = M3
    stellar = SeBa()
    stellar.particles.add_particles(triple)
    channel_from_stellar = stellar.particles.new_channel_to(triple)

    # Evolve to t_stellar.
    
    stellar.evolve_model(t_stellar)
    channel_from_stellar.copy_attributes(["mass"])
    M1 = triple[0].mass
    M2 = triple[1].mass
    M3 = triple[2].mass
    print("t=", stellar.model_time.in_(units.Myr))
    print("M=", stellar.particles.mass.in_(units.MSun))
    print("R=", stellar.particles.radius.in_(units.RSun))
    print("L=", stellar.particles.luminosity.in_(units.LSun))
    print("T=", stellar.particles.temperature.in_(units.K))
    print("Mdot=", \
        -stellar.particles.wind_mass_loss_rate.in_(units.MSun/units.yr))

    # Start the dynamics.
    # Inner binary:
    
    tmp_stars = Particles(2)
    tmp_stars[0].mass = M1
    tmp_stars[1].mass = M2

    if Pora == 1:
        ain_0 = semimajor_axis(Pin_0, M1+M2)
    else:
        Pin_0 = orbital_period(ain_0, M1+M2)
    print('Pin =', Pin_0)
        
    print('ain_0 =', ain_0)
    print('M1+M2 =', M1+M2)
    print('Pin_0 =', Pin_0.value_in(units.day), '[day]')
    #print 'semi:', semimajor_axis(Pin_0, M1+M2).value_in(units.AU), 'AU'
    #print 'period:', orbital_period(ain_0, M1+M2).value_in(units.day), '[day]'
    
    dt_init = 0.01*Pin_0
    ma = 180
    inc = 60
    aop = 180
    lon = 0
    r,v = get_position(M1, M2, ein_0, ain_0, ma, inc, aop, lon, dt_init)
    tmp_stars[1].position = r
    tmp_stars[1].velocity = v
    tmp_stars.move_to_center()

    # Outer binary:
    
    r,v = get_position(M1+M2, M3, eout_0, aout_0, 0, 0, 0, 0, dt_init)
    tertiary = Particle()
    tertiary.mass = M3
    tertiary.position = r
    tertiary.velocity = v
    tmp_stars.add_particle(tertiary)
    tmp_stars.move_to_center()

    triple.position = tmp_stars.position
    triple.velocity = tmp_stars.velocity

    Mtriple = triple.mass.sum()
    Pout = orbital_period(aout_0, Mtriple)

    print("T=", stellar.model_time.in_(units.Myr))
    print("M=", stellar.particles.mass.in_(units.MSun))
    print("Pout=", Pout.in_(units.Myr))
    print('tK =', ((M1+M2)/M3)*Pout**2*(1-eout_0**2)**1.5/Pin_0)

    converter = nbody_system.nbody_to_si(triple.mass.sum(), aout_0)

    if integrator == 0:
        gravity = Hermite(converter)
        gravity.parameters.timestep_parameter = 0.01
    elif integrator == 1:
        gravity = SmallN(converter)
        gravity.parameters.timestep_parameter = 0.01
        gravity.parameters.full_unperturbed = 0
    elif integrator == 2:
        gravity = Huayno(converter)
        gravity.parameters.inttype_parameter = 20
        gravity.parameters.timestep = (1./256)*Pin_0
    else:
        gravity = symple(converter)
        gravity.parameters.integrator = 10
        #gravity.parameters.timestep_parameter = 0.
        gravity.parameters.timestep = (1./128)*Pin_0

    print(gravity.parameters)

    gravity.particles.add_particles(triple)
    channel_from_framework_to_gd = triple.new_channel_to(gravity.particles)
    channel_from_gd_to_framework = gravity.particles.new_channel_to(triple)
    
    Etot_init = gravity.kinetic_energy + gravity.potential_energy
    Etot_prev = Etot_init

    gravity.particles.move_to_center()

    # Note: time = t_diag = 0 at the start of the dynamical integration.
    
    dt_diag = t_end/float(nsteps)
    t_diag = dt_diag
    time = 0.0 | t_end.unit
    t_se = t_stellar + time

    print('t_end =', t_end)
    print('dt_diag =', dt_diag)

    ain, ein, aout, eout = get_orbital_elements_of_triple(triple)
    print("Triple elements t=",  time,  \
        "inner:", triple[0].mass, triple[1].mass, ain, ein, \
        "outer:", triple[2].mass, aout, eout)

    t = [time.value_in(units.Myr)]
    Mtot = triple.mass.sum()
    mtot = [Mtot.value_in(units.MSun)]
    smai = [ain/ain_0] 
    ecci = [ein/ein_0]
    smao = [aout/aout_0] 
    ecco = [eout/eout_0]

    if interp:
        
        # Create arrays of stellar times and masses for interpolation.

        times = [time]
        masses = [triple.mass.copy()]
        while time < t_end:
            time += dt_se
            stellar.evolve_model(t_stellar+time)
            channel_from_stellar.copy_attributes(["mass"])
            times.append(time)
            masses.append(triple.mass.copy())

        time = 0.0 | t_end.unit
        print('\ntimes:', times, '\n')

    # Evolve the system.
    
    def advance_stellar(t_se, dt):
        E0 = gravity.kinetic_energy + gravity.potential_energy
        t_se += dt

        if interp:
            t = t_se-t_stellar
            i = int(t/dt_se)
            mass = masses[i] + (t-times[i])*(masses[i+1]-masses[i])/dt_se
            triple.mass = mass
            #print 't_se =', t_se, 'masses =', mass
        else:
            stellar.evolve_model(t_se)
            channel_from_stellar.copy_attributes(["mass"])

        channel_from_framework_to_gd.copy_attributes(["mass"])
        return t_se, gravity.kinetic_energy + gravity.potential_energy - E0

    def advance_gravity(tg, dt):
        tg += dt
        gravity.evolve_model(tg)
        channel_from_gd_to_framework.copy()
        return tg

    while time < t_end:

        if scheme == 1:

            # Advance to the next diagnostic time.
            
            dE_se = zero
            dt = t_diag - time

            if dt > 0|dt.unit:
                time = advance_gravity(time, dt)

        elif scheme == 2:
            
            # Derive dt from Pin using dtse_fac.
            
            dt = dtse_fac*Pin_0
            if time + dt > t_diag: dt = t_diag - time

            if dt > 0|dt.unit:
                t_se, dE_se = advance_stellar(t_se, dt)
                time = advance_gravity(time, dt)
            
        elif scheme == 3:
            
            # Derive dt from Pin using dtse_fac.
            
            dt = dtse_fac*Pin_0
            if time + dt > t_diag: dt = t_diag - time

            if dt > 0|dt.unit:
                time = advance_gravity(time, dt)
                t_se, dE_se = advance_stellar(t_se, dt)
            
        elif scheme == 4:
            
            # Derive dt from Pin using dtse_fac.
            
            dt = dtse_fac*Pin_0
            if time + dt > t_diag: dt = t_diag - time

            if dt > 0|dt.unit:
                t_se, dE_se = advance_stellar(t_se, 0.5*dt)
                time = advance_gravity(time, dt)
                t_se, dE_se2 = advance_stellar(t_se, 0.5*dt)
                dE_se += dE_se2
            
        elif scheme == 5:

            # Use the specified dt_se.
            
            dE_se = zero
            dt = dt_se
            if time + dt > t_diag: dt = t_diag - time

            if dt > 0|dt.unit:

                # For use with symple only: set up average mass loss.
    
                channel_from_stellar.copy_attributes(["mass"])
                m0 = triple.mass.copy()
                stellar.evolve_model(t_se+dt)
                channel_from_stellar.copy_attributes(["mass"])
                t_se = stellar.model_time
                m1 = triple.mass
                dmdt = (m1-m0)/dt
                for i in range(len(dmdt)):
                    gravity.set_dmdt(i, dmdt[i])

                time = advance_gravity(time, dt)

        else:

            print('unknown option')
            sys.exit(0)

        if time >= t_diag:
            
            t_diag = time + dt_diag

            Ekin = gravity.kinetic_energy 
            Epot = gravity.potential_energy
            Etot = Ekin + Epot
            dE = Etot_prev - Etot
            Mtot = triple.mass.sum()
            print("T=", time, end=' ') 
            print("M=", Mtot, "(dM[SE]=", Mtot/Mtriple, ")", end=' ')
            print("E= ", Etot, "Q= ", Ekin/Epot, end=' ')
            print("dE=", (Etot_init-Etot)/Etot, "ddE=", (Etot_prev-Etot)/Etot, end=' ') 
            print("(dE[SE]=", dE_se/Etot, ")")
            Etot_init -= dE
            Etot_prev = Etot
            ain, ein, aout, eout = get_orbital_elements_of_triple(triple)
            print("Triple elements t=",  t_stellar + time,  \
                "inner:", triple[0].mass, triple[1].mass, ain, ein, \
                "outer:", triple[2].mass, aout, eout)

            t.append(time.value_in(units.yr))
            mtot.append(Mtot.value_in(units.MSun))
            smai.append(ain/ain_0)
            ecci.append(ein/ein_0)
            smao.append(aout/aout_0)
            ecco.append(eout/eout_0)

            if eout > 1 or aout <= zero:
                print("Binary ionized or merged")
                break

    gravity.stop()
    stellar.stop()

    return t, mtot, smai, ecci, smao, ecco

def main(M1, M2, M3, Pora, Pin, ain, aout, ein, eout,
         t_end, nsteps, scheme, integrator,
         t_stellar, dt_se, dtse_fac, interp, show):

    two_frames = False
    plot_ae = True

    color = get_distinct(4)
    if two_frames:
        plt.figure(figsize=(10, 8))
    else:
        plt.figure(figsize=(12, 8))

    if scheme == 5:
        if integrator != 3:
            print('Warning: scheme = 5 forces integrator = 3')
        integrator = 3
    srange = [1,scheme]		# 1 = no mass loss; other = mass loss
    				# assume scheme > 1
    
    i = 0
    lw = [1,2]
    
    srange = [1, scheme]

    for s in srange:
        time, mtot, ai, ei, ao, eo \
            = evolve_triple_with_wind(M1, M2, M3,
                                      Pora, Pin, ain, aout,
                                      ein, eout,
                                      t_end, nsteps,
                                      s, integrator,
                                      t_stellar, dt_se,
                                      dtse_fac, interp)
        if i == 0:
            if two_frames: plt.subplot(1,2,1)
            plt.plot(time, ai, c=color[0], linewidth=lw[i],
                     label='inner, no mass loss')
            plt.plot(time, ao, c=color[3], linewidth=lw[i],
                     label='outer, no mass loss')
            plt.xlabel('time (yr)')
            plt.ylabel('$a/a_0$')
            if two_frames:
                plt.subplot(1,2,2)
                if plot_ae:
                    plt.plot(ai, ei, c=color[0], linewidth=lw[i])
                    plt.plot(ao, eo, c=color[3], linewidth=lw[i])
                    plt.xlabel('$a/a_0$')
                    plt.ylabel('$e/e_0$')
                else:
                    plt.plot(time, mtot, c=color[0], linewidth=lw[i])
                    plt.xlabel('time (yr)')
                    plt.ylabel('M')
            i = 1
        else:
            if two_frames: plt.subplot(1,2,1)
            plt.plot(time, ai, c=color[1], linewidth=lw[i],
                     label='inner, mass loss')
            plt.plot(time, ao, c=color[2], linewidth=lw[i],
                     label='outer, mass loss')
            if two_frames:
                plt.subplot(1,2,2)
                if plot_ae:
                    plt.plot(ai, ei, c=color[1], linewidth=lw[i])
                    plt.plot(ao, eo, c=color[2], linewidth=lw[i])
                else:
                    plt.plot(time, mtot, c=color[1], linewidth=lw[i])

    if two_frames: plt.subplot(1,2,1)
    plt.legend(loc='best')

    integrators = ['hermite', 'smalln', 'huayno', 'symple']
    label = integrators[integrator]
    label += ' integrator, stellev scheme= {:d}'.format(scheme)
    save_file \
        = 'evolve_triple_with_wind_t={:.3f}'.format(t_end.value_in(units.Myr)) \
        	+'_i={:d}_s={:d}'.format(integrator, scheme)
    if scheme < 5:
        label += ', dtse_fac = {:.3f}'.format(dtse_fac)
        save_file += '_dtsefac={:.3f}'.format(dtse_fac)
    else:
        label += ', dt_se = {:.1f}'.format(dt_se.value_in(units.yr))
        save_file += '_dtse={:.1f}'.format(dt_se.value_in(units.yr))
    save_file += '.png'

    if two_frames:
        plt.tight_layout()
        plt.subplots_adjust(top=0.88)
    #plt.suptitle(label, y=0.97, fontsize=15)

    ax = plt.gca()
    ax.minorticks_on()		# switch on the minor ticks
    ax.tick_params(axis='both', which='both', direction='in')
    ax.locator_params(nbins=3)
    ax.get_yaxis().get_major_formatter().set_useOffset(False)

    plt.savefig(save_file, dpi=300)
    print('\nSaved figure in file', save_file,'\n')
    if show: plt.show()
    
def new_option_parser():
    
    from amuse.units.optparse import OptionParser
    result = OptionParser()
    result.add_option("--ain", unit=units.AU,
                      dest="ain", type="float", default = 0.63|units.AU,
                      help="orbital separation [%default]")
    result.add_option("--aout", unit=units.AU,
                      dest="aout", type="float", default = 100|units.AU,
                      help="orbital separation [%default]")
    result.add_option("--dtse", unit=units.Myr,
                      dest="dt_se", type="float", default = 1.e-3|units.Myr,
                      help="stellar mass-loss time step [%default]")
    result.add_option("--dtse_fac",
                      dest="dtse_fac", type="float", default = 0.1,
                      help="stellar mass-loss time step fraction [%default]")
    result.add_option("--ein",
                      dest="ein", type="float", default = 0.2,
                      help="orbital eccentricity [%default]")
    result.add_option("--eout",
                      dest="eout", type="float", default = 0.6,
                      help="orbital eccentricity [%default]")
    result.add_option("-i",
                      dest="integrator", type="int", default = 2,
                      help="integration scheme [%default]")
    result.add_option("-I",
                      dest="interp", action="store_false", default = True,
                      help="interpolate stellar evolution [%default]")
    result.add_option("--M1", unit=units.MSun,
                      dest="M1", type="float", default = 60|units.MSun,
                      help="Primary mass [%default]")
    result.add_option("--M2", unit=units.MSun,
                      dest="M2", type="float", default = 30|units.MSun,
                      help="secondary mass [%default]")
    result.add_option("--M3", unit=units.MSun,
                      dest="M3", type="float", default = 20|units.MSun,
                      help="secondary mass [%default]")
    result.add_option("-n",
                      dest="nsteps", type="int", default = 1000,
                      help="number of data points [%default]")
    result.add_option("--Pin", unit=units.day,
                      dest="Pin", type="float", default = 19|units.day,
                      help="orbital period [%default]")
    result.add_option("--Pora",
                      dest="Pora", type="int", default = 1,
                      help="period (1) or semimajor axis (2) [%default]")
    result.add_option("-s",
                      dest="scheme", type="int", default = 3,
                      help="stellar integration method [%default]")
    result.add_option("-S",
                      dest="show", action="store_false", default = True,
                      help="show plot on display [%default]")
    result.add_option("-t", unit=units.Myr,
                      dest="t_end", type="float", default = 1.e-3|units.Myr,
                      help="end time of the dynamical simulation [%default]")
    result.add_option("--ts", unit=units.Myr,
                      dest="t_stellar", type="float", default = 4.|units.Myr,
                      help="stellar evolution time [%default]")
    return result

if __name__ in ('__main__', '__plot__'):
    
    #set_printing_strategy("custom", 
    #                      preferred_units = [units.MSun, units.AU, units.Myr], 
    #                      precision = 12, prefix = "", 
    #                      separator = " [", suffix = "]")

    o, arguments  = new_option_parser().parse_args()
    print(o.__dict__)
    main(**o.__dict__)
