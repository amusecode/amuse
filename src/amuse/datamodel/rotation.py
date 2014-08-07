import numpy

def new_rotation_matrix(phi, theta, psi):
    """
    Return the rotation matrix, to rotate positions, around the x-axis (phi), y-axis (theta) and z-axis (psi).
    See wikipedia for reference
    """
    cos = numpy.cos
    sin = numpy.sin
    
    return numpy.array( (
        (cos(theta)*cos(psi), -cos(phi)*sin(psi) + sin(phi)*sin(theta)*cos(psi),  sin(phi)*sin(psi) + cos(phi)*sin(theta)*cos(psi)) ,
        (cos(theta)*sin(psi),  cos(phi)*cos(psi) + sin(phi)*sin(theta)*sin(psi), -sin(phi)*cos(psi) + cos(phi)*sin(theta)*sin(psi)) ,
        (-sin(theta)        ,  sin(phi)*cos(theta)                             ,  cos(phi)*cos(theta))
    ) )

def rotated(positions, phi, theta, psi):
    """
    Return the positions, rotated by phi, theta and psi around the x, y and z axes
    """
#    print phi, theta, psi
    rotation_matrix = new_rotation_matrix(phi, theta, psi)
#    print "RT=", rotation_matrix
    return positions.dot(rotation_matrix.transpose())

def rotate(particles, phi, theta, psi):
    """
    Rotate the positions and the velocities around 0,0,0.
    """
    particles.position = rotated(particles.position,  phi, theta, psi)
    particles.velocity = rotated(particles.velocity,  phi, theta, psi)

def add_spin(particles, omega):
    """
    Add solid-body rotation to the velocity of the particles, relative to the 
    center-of-mass position.
    """
    if not omega.is_vector():
        omega = omega * [0.0, 0.0, 1.0]
    particles.velocity += omega.cross(particles.position - particles.center_of_mass())

#___________________________________________________________

def from_inertial_to_rotating_frame(particles, initial_phase, omega, time):
    '''
    This function transforms from an inertial system to a rotating one.
    particles must be defined in an inertial frame
   '''

    angle= initial_phase + omega*time
    C1= particles.vx + omega*particles.y
    C2= particles.vy - omega*particles.x
    rotating_x = particles.x*numpy.cos(angle) + particles.y*numpy.sin(angle)
    rotating_y = -particles.x*numpy.sin(angle) + particles.y*numpy.cos(angle) 
    rotating_z = particles.z
    rotating_vx = C1*numpy.cos(angle) + C2*numpy.sin(angle) 
    rotating_vy = C2*numpy.cos(angle) - C1*numpy.sin(angle)
    rotating_vz = particles.vz
    return rotating_x, rotating_y, rotating_z, rotating_vx, rotating_vy, rotating_vz


def creation_particle_set_in_rotating_frame(particles, initial_phase, omega, time, stellar_evolution=False):
    '''
    This function creates a particle set
    in a rotating frame. Input:
    particles: set located in an inertial frame
    initial_phase= initial angle of the rotating frame
    omega= pattern speed of the rotating frame
    stellar evolution: flag. If it is true other attributes of 
    particles are copied to the new set.
    '''
    rot_x, rot_y, rot_z, rot_vx, rot_vy, rot_vz= from_inertial_to_rotating_frame(particles, initial_phase, omega, time) 
    no_inertial_system= particles.copy()
    no_inertial_system.x = rot_x
    no_inertial_system.y = rot_y
    no_inertial_system.z = rot_z
    no_inertial_system.vx = rot_vx
    no_inertial_system.vy = rot_vy
    no_inertial_system.vz = rot_vz
    if stellar_evolution== True:
        no_inertial_system.age= particles.age
        no_inertial_system.mass= particles.mass
        no_inertial_system.radius= particles.radius
        no_inertial_system.luminosity= particles.luminosity
        no_inertial_system.temperature= particles.temperature
        no_inertial_system.stellar_type= particles.stellar_type
   
    return no_inertial_system 


def in_rotating_frame(particles, initial_phase, omega, time):
    '''
    This function places a particle set into a rotating system that moves
    with a pattern speed equals to omega. Input:
    particles: set located in an inertial frame
    initial_phase= initial angle of the rotating frame
    omega= pattern speed of the rotating frame
    stellar evolution: flag. If it is true other attributes of 
    particles are updated to the new set.
    
    '''    
    rot_x, rot_y, rot_z, rot_vx, rot_vy, rot_vz= from_inertial_to_rotating_frame(particles, initial_phase, omega, time)     
    particles.x= rot_x
    particles.y= rot_y
    particles.z= rot_z
    particles.vx= rot_vx
    particles.vy= rot_vy
    particles.vz= rot_vz

def from_rotating_to_inertial_frame(particles, initial_phase, omega, time):
    '''
    This function transforms from a rotating system to an inertial one
    particles are defined in a rotating frame
    '''
    angle= initial_phase + omega*time
    C1= part_noin.vx - part_noin.y*omega
    C2= part_noin.vy + part_noin.x*omega
  
    inertial.x= particles.x*numpy.cos(angle)-particles.y*numpy.sin(angle)
    inertial.y= particles.x*numpy.sin(angle)+particles.y*numpy.cos(angle)
    inertial.z= particles.z
    inertial.vx= C1*numpy.cos(angle) - C2*numpy.sin(angle)
    ienrtial.vy= C1*numpy.sin(angle) + C2*numpy.cos(angle)
    inertial.vz= particles.vz
    return inertial.x, inertial.y, inertial.z, inertial.vx, inertial.vy, inertial.vz

#_______________________________________

def creation_particle_set_in_inertial_frame(particles, initial_phase, omega, time, stellar_evolution=False):
    '''
    This function creates a particle set
    in an inertial frame. Input:
    particles: set located in a rotating frame
    initial_phase= initial angle of the rotating frame
    omega= pattern speed of the rotating frame
    stellar evolution: flag. If it is true other attributes of 
    particles are copied to the new set.
    '''
    inertial_x, inertial_y, inertial_z, inertial_vx, inertial_vy, inertial_vz= from_rotating_to_inertial_frame(particles, initial_phase, omega, time) 

    inertial_system= particles.copy()
    inertial_system.x = rot_x
    inertial_system.y = rot_y
    inertial_system.z = rot_z
    inertial_system.vx = rot_vx
    inertial_system.vy = rot_vy
    inertial_system.vz = rot_vz
    if stellar_evolution== True:
        inertial_system.age= particles.age
        inertial_system.mass= particles.mass
        inertial_system.radius= particles.radius
        inertial_system.luminosity= particles.luminosity
        inertial_system.temperature= particles.temperature
        inertial_system.stellar_type= particles.stellar_type
   
    return inertial_system 

def in_inertial_frame(particles, initial_phase, omega, time):
    '''
    This function places a particle set into a rotating system that moves
    with a pattern speed equals to omega. Input:
    particles: set located in a rotating frame
    initial_phase= initial angle of the rotating frame
    omega= pattern speed of the rotating frame
    stellar evolution: flag. If it is true other attributes of 
    particles are updated to the new set.
    
    '''    
    inertial_x, inertial_y, inertial_z, inertial_vx, inertial_vy, inertial_vz= from_rotating_to_inertial_frame(particles, initial_phase, omega, time)      
    particles.x= inertial_x
    particles.y= inertial_y
    particles.z= inertial_z
    particles.vx= inertial_vx
    particles.vy= inertial_vy
    particles.vz= inertial_vz
