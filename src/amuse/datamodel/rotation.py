import numpy
from amuse.units.trigo import cos, sin

def new_rotation_matrix(phi, theta, psi):
    """
    Return the rotation matrix, to rotate positions, around the x-axis (phi), y-axis (theta) and z-axis (psi).
    See wikipedia for reference
    """
    
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

