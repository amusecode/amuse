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
    rotation_matrix = new_rotation_matrix(phi, theta, psi)
    units = positions.unit
    value = positions.value_in(units)
    
    result = []
    for x in value:
        result.append(numpy.dot(rotation_matrix, x))
    return units.new_quantity(result)

def rotate(particles, phi, theta, psi):
    """
    Rotate the positions and the velocities around 0,0,0.
    """
    particles.position = rotated(particles.position,  phi, theta, psi)
    particles.velocity = rotated(particles.velocity,  phi, theta, psi)

