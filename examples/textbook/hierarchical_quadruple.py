"""
Simple script to run a `2+2' quadruple system with SecularMultiple.
The two `inner' binaries are denoted with `A' and `B'; the wider outer binary (`superorbit') is denoted with `C'. 
Orbital parameters can be provided with command line arguments.
Note: setting N_output to a large value will slow down the script due to Python overhead, but will make nicer-looking plots.

Adrian Hamers, December 2017
"""

import numpy

from amuse.community.secularmultiple.interface import SecularMultiple
from amuse.units import quantities,units,constants
from amuse.datamodel import Particles
from matplotlib import pyplot

###BOOKLISTSTART1###
def initialize_multiple_system(N_bodies, masses, semimajor_axis, eccentricity, inclination, argument_of_pericenter, longitude_of_ascending_node):

    N_binaries = N_bodies-1
    particles = Particles(N_bodies+N_binaries)
    for index in range(N_bodies):
        particle = particles[index]
        particle.mass = masses[index]
        particle.is_binary = False
        particle.radius = 1.0 | units.RSun
        particle.child1 = None
        particle.child2 = None

    for index in range(N_binaries):
        particle = particles[index+N_bodies]
        particle.is_binary = True
        particle.semimajor_axis = semimajor_axis[index]
        particle.eccentricity = eccentricity[index]
        particle.inclination = inclination[index]
        particle.argument_of_pericenter = argument_of_pericenter[index]
        particle.longitude_of_ascending_node = longitude_of_ascending_node[index]
        
        # Specify the `2+2' hierarchy:
        
        if index==0:
            particle.child1 = particles[0]
            particle.child2 = particles[1]
        elif index==1:
            particle.child1 = particles[2]
            particle.child2 = particles[3]
        elif index==2:
            particle.child1 = particles[4]
            particle.child2 = particles[5]        
    binaries = particles[particles.is_binary]

    return particles, binaries
###BOOKLISTSTOP1###

def evolve_quadruple(N_output, end_time, m1, m2, m3, m4, aA, aB, aC, eA, eB, eC, iA, iB, iC, ApA, ApB, ApC, LANA, LANB, LANC):

    masses = [m1, m2, m3, m4]
    semimajor_axis = [aA, aB, aC]
    eccentricity = [eA, eB, eC]
    inclination = numpy.deg2rad([iA, iB, iC])
    argument_of_percienter = numpy.deg2rad([ApA, ApB, ApC])
    longitude_of_ascending_node = numpy.deg2rad([LANA, LANB, LANC])
    print(longitude_of_ascending_node)
    
    N_bodies = 4
    N_binaries = N_bodies-1
    particles, binaries = initialize_multiple_system(N_bodies, masses, semimajor_axis, eccentricity, inclination, argument_of_percienter, longitude_of_ascending_node)

    code = SecularMultiple()
    code.particles.add_particles(particles)

    channel_from_particles_to_code = particles.new_channel_to(code.particles)
    channel_from_code_to_particles = code.particles.new_channel_to(particles)
    channel_from_particles_to_code.copy()

    ### set up some arrays for plotting ###
    print_smas_AU = [[] for x in range(N_binaries)]
    print_rps_AU = [[] for x in range(N_binaries)]
    print_parent_is_deg = [[] for x in range(N_binaries)]
    print_times_Myr = []

    time = 0.0|units.yr
    output_time_step = end_time/float(N_output)
    while time <= end_time:
        time += output_time_step
        code.evolve_model(time)

        channel_from_code_to_particles.copy()
        print('='*50)
        print('t/Myr',time.value_in(units.Myr))
        print('e',binaries.eccentricity)
        print('i/deg', numpy.rad2deg(binaries.inclination))
        print('AP/deg', \
            numpy.rad2deg(binaries.argument_of_pericenter))  
        print('LAN/deg', \
            numpy.rad2deg(binaries.longitude_of_ascending_node))
            
        ### write to output arrays ###
        print_times_Myr.append(time.value_in(units.Myr))
        for index_binary in range(N_binaries):
            print_smas_AU[index_binary].append( binaries[index_binary].semimajor_axis.value_in(units.AU) )
            print_rps_AU[index_binary].append( binaries[index_binary].semimajor_axis.value_in(units.AU)*(1.0 - binaries[index_binary].eccentricity) )
            print_parent_is_deg[index_binary].append( numpy.rad2deg(binaries[index_binary].inclination_relative_to_parent) )

    ### compute the `canonical' maximum eccentricity/periapsis distance that applies in the quadrupole-order test-particle limit if the `outer' binary is replaced by a point mass ###
    print(inclination[0],inclination[2],longitude_of_ascending_node[0],longitude_of_ascending_node[2])
    i_AC_init = compute_mutual_inclination(inclination[0],inclination[2],longitude_of_ascending_node[0],longitude_of_ascending_node[2])
    i_BC_init = compute_mutual_inclination(inclination[1],inclination[2],longitude_of_ascending_node[1],longitude_of_ascending_node[2])
    
    canonical_rp_min_A_AU = (semimajor_axis[0]*(1.0 - numpy.sqrt( 1.0 - (5.0/3.0)*numpy.cos(i_AC_init)**2 ) )).value_in(units.AU)
    canonical_rp_min_B_AU = (semimajor_axis[1]*(1.0 - numpy.sqrt( 1.0 - (5.0/3.0)*numpy.cos(i_BC_init)**2 ) )).value_in(units.AU)

    data = print_times_Myr,print_smas_AU,print_rps_AU,print_parent_is_deg,canonical_rp_min_A_AU,canonical_rp_min_B_AU
    return data

def compute_mutual_inclination(INCL_k,INCL_l,LAN_k,LAN_l):
    cos_INCL_rel = numpy.cos(INCL_k)*numpy.cos(INCL_l) + numpy.sin(INCL_k)*numpy.sin(INCL_l)*numpy.cos(LAN_k-LAN_l)
    return numpy.arccos(cos_INCL_rel)

def plot_function(data):
    print_times_Myr,print_smas_AU,print_rps_AU,print_parent_is_deg,canonical_rp_min_A_AU,canonical_rp_min_B_AU = data

    N_binaries = len(print_smas_AU)

    pyplot.rc('text',usetex=True)
    pyplot.rc('legend',fancybox=True)
            
    linewidth=4
    dlinewidth=2
    fig=pyplot.figure(figsize=(10,9))
    plot1=fig.add_subplot(2,1,1,yscale="log")
    plot2=fig.add_subplot(2,1,2)

    from distinct_colours import get_distinct

    colors = get_distinct(4)
    labels = ["$A$","$B$","$C$"]
    labels_i = ["$i_{AC}$","$i_{BC}$",None]
    for index_binary in range(N_binaries):
        label = labels[index_binary]
        label_i = labels_i[index_binary]        
        color = colors[index_binary]
        
        plot1.plot(print_times_Myr,print_smas_AU[index_binary],color=color,linestyle='dashed',linewidth=dlinewidth)
        plot1.plot(print_times_Myr,print_rps_AU[index_binary],color=color,linewidth=linewidth,label=label)
        
        plot2.plot(print_times_Myr,print_parent_is_deg[index_binary],color=color,linewidth=linewidth,label=label_i)

    plot1.axhline(y = canonical_rp_min_A_AU, color= colors[0],linestyle='dotted',linewidth=dlinewidth)
    plot1.axhline(y = canonical_rp_min_B_AU, color= colors[1],linestyle='dotted',linewidth=dlinewidth)

    handles,labels = plot1.get_legend_handles_labels()
    plot1.legend(handles,labels,loc="upper right",fontsize=12)

    handles,labels = plot2.get_legend_handles_labels()
    plot2.legend(handles,labels,loc="lower right",fontsize=12)

    plot1.set_xlabel("t [Myr]",fontsize=18)
    plot2.set_xlabel("t [Myr]",fontsize=18)

    plot1.set_ylabel("$a_i [\mathrm{AU}]$",fontsize=18)
    plot2.set_ylabel("$i_{kl} [\mathrm{deg}]$",fontsize=18)

    plot1.set_xlim(0.0,print_times_Myr[-1])
    plot2.set_xlim(0.0,print_times_Myr[-1])

    plot1.tick_params(axis='both', which ='major', labelsize = 18)
    plot2.tick_params(axis='both', which ='major', labelsize = 18)

    fig.savefig("figure.eps")

    pyplot.show()

def new_option_parser():
    from amuse.units.optparse import OptionParser
    result = OptionParser()
    result.add_option("--end_time", unit=units.Myr,
                      dest="end_time", type="float", default = 5.0|units.Myr,
                      help="integration time [%default]")
    result.add_option("--N_output", 
                      dest="N_output", type="int", default = 400,
                      help="number of output steps [%default]")
    result.add_option("--m1", unit=units.MSun,
                      dest="m1", type="float", default = 1.0|units.MSun,
                      help="mass of object 1 [%default]")
    result.add_option("--m2", unit=units.MSun,
                      dest="m2", type="float", default = 0.8|units.MSun,
                      help="mass of object 2 [%default]")
    result.add_option("--m3", unit=units.MSun,
                      dest="m3", type="float", default = 1.1|units.MSun,
                      help="mass of object 3 [%default]")
    result.add_option("--m4", unit=units.MSun,
                      dest="m4", type="float", default = 0.9|units.MSun,
                      help="mass of object 4 [%default]")
    result.add_option("--aA", unit=units.AU,
                      dest="aA", type="float", default = 1.0|units.AU,
                      help="semimajor axis of orbit A [%default]")
    result.add_option("--aB", unit=units.AU,
                      dest="aB", type="float", default = 1.2|units.AU,
                      help="semimajor axis of orbit B [%default]")
    result.add_option("--aC", unit=units.AU,
                      dest="aC", type="float", default = 100.0|units.AU,
                      help="semimajor axis of orbit C (the `superorbit') [%default]")
    result.add_option("--eA",
                      dest="eA", type="float", default = 0.1,
                      help="eccentricity of orbit A [%default]")
    result.add_option("--eB",
                      dest="eB", type="float", default = 0.1,
                      help="eccentricity of orbit B [%default]")
    result.add_option("--eC",
                      dest="eC", type="float", default = 0.3,
                      help="eccentricity of orbit C (the `superorbit') [%default]")
    result.add_option("--iA",
                      dest="iA", type="float", default = 75.0,
                      help="inclination of orbit A in degrees [%default]")
    result.add_option("--iB",
                      dest="iB", type="float", default = 80.0,
                      help="inclination of orbit B in degrees [%default]")
    result.add_option("--iC",
                      dest="iC", type="float", default = 0.001,
                      help="inclination of orbit C (the `superorbit') in degrees [%default]")
    result.add_option("--ApA",
                      dest="ApA", type="float", default = 10.0,
                      help="argument of periapsis of orbit A in degrees [%default]")
    result.add_option("--ApB",
                      dest="ApB", type="float", default = 30.0,
                      help="argument of periapsis of orbit B in degrees [%default]")
    result.add_option("--ApC",
                      dest="ApC", type="float", default = 60.0,
                      help="argument of periapsis of orbit C (the `superorbit') in degrees [%default]")
    result.add_option("--LANA",
                      dest="LANA", type="float", default = 0.001,
                      help="longitude of the ascending node of orbit A in degrees [%default]")
    result.add_option("--LANB",
                      dest="LANB", type="float", default = 0.001,
                      help="longitude of the ascending node of orbit B in degrees [%default]")
    result.add_option("--LANC",
                      dest="LANC", type="float", default = 0.001,
                      help="longitude of the ascending node of orbit C (the `superorbit') in degrees [%default]")
                      
    return result
    

if __name__ in ('__main__', '__plot__'):
    o, arguments  = new_option_parser().parse_args()
    data = evolve_quadruple(**o.__dict__)
    
    plot_function(data)
