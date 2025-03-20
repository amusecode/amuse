import numpy
import amuse.community.twobody.twobody as twobody


from amuse.units.units import *
if __name__=='__main__':
    pos_earth=numpy.array([ 8.418982185410142E-01, 
                            5.355823303978186E-01, 
                            2.327960005926782E-05]) | AU 
    vel_earth=numpy.array([-9.488931818313919E-03, 
                            1.447515189957170E-02, 
                            3.617712172296458E-07]) | AUd
    pos_moon=numpy.array([8.426656721530955E-01, 
                          5.331110650437484E-01, 
                         -6.837900390288286E-05]) | AU
    vel_moon=numpy.array([-8.943352544499154E-03, 
                           1.467795416516487E-02, 
                           4.840393580601162E-05]) | AUd
    mass_earth=5.9742e24 | kg
    mass_moon=7.3477e22 | kg 
    radius_earth=6371 | km
    radius_moon=1737.1 | km

    cmpos=(pos_earth*mass_earth+pos_moon*mass_moon)/(mass_earth+mass_moon)
    cmvel=(vel_earth*mass_earth+vel_moon*mass_moon)/(mass_earth+mass_moon)
    
    pos_earth=pos_earth-cmpos
    vel_earth=vel_earth-cmvel
    pos_moon=pos_moon-cmpos
    vel_moon=vel_moon-cmvel

    nb=twobody.TwoBody()
    
    earth,err=nb.new_particle(mass_earth.value_in(kg),radius_earth.value_in(m),  
      pos_earth[0].value_in(m),pos_earth[1].value_in(m),pos_earth[2].value_in(m), 
      vel_earth[0].value_in(m/s),vel_earth[1].value_in(m/s),vel_earth[2].value_in(m/s))
    moon,err=nb.new_particle(mass_moon.value_in(kg),radius_moon.value_in(m),  
      pos_moon[0].value_in(m),pos_moon[1].value_in(m),pos_moon[2].value_in(m), 
      vel_moon[0].value_in(m/s),vel_moon[1].value_in(m/s),vel_moon[2].value_in(m/s))

    nb.evolve_model( (27.321582| day ).value_in(s))
    moonstate,err=nb.get_state(moon)
    print(pos_moon[0].value_in(km), (moonstate['x'] | m).value_in(km))
    print(pos_moon[1].value_in(km), (moonstate['y'] | m).value_in(km))
    print(pos_moon[2].value_in(km), (moonstate['z'] | m).value_in(km))
    earthstate,err=nb.get_state(earth)
    print(pos_earth[0].value_in(km), (earthstate['x'] | m).value_in(km))
    print(pos_earth[1].value_in(km), (earthstate['y'] | m).value_in(km))
    print(pos_earth[2].value_in(km), (earthstate['z'] | m).value_in(km))
