import numpy

from amuse.units import units

from amuse.datamodel import Particle
from amuse.datamodel import Particles
_solsysdat= \
 [['MERCURY',1.66013679527193009E-07,20.,5.43, \
 -3.83966017419175965E-01, -1.76865300855700736E-01, 2.07959213998758705E-02, \
  5.96286238644834141E-03, -2.43281292146216750E-02,-2.53463209848734695E-03, \
  0., 0., 0.],
 ['VENUS',2.44783833966454430E-06,20.,5.24, \
  6.33469157915745540E-01, 3.49855234102151691E-01,-3.17853172088953667E-02, \
 -9.84258038001823571E-03,  1.76183746921837227E-02, 8.08822351013463794E-04, \
  0., 0., 0.],
  ['EARTHMOO',3.04043264264672381E-06,20.,5.52, \
  2.42093942183383037E-01, -9.87467766698604366E-01, -4.54276292555233496E-06, \
  1.64294055023289365E-02,  4.03200725816140870E-03,  1.13609607260006795E-08, \
  0., 0., 0.],
 ['MARS',3.22715144505386530E-07,20.,3.94, \
  2.51831018120174499E-01,  1.52598983115984788E+00,  2.57781137811807781E-02, \
 -1.32744166042475433E-02,  3.46582959610421387E-03,  3.98930013246952611E-04, \
  0., 0., 0.],
 ['JUPITER',9.54791938424326609E-04,3.,1.33, \
  4.84143144246472090E+00, -1.16032004402742839E+00, -1.03622044471123109E-01, \
  1.66007664274403694E-03,  7.69901118419740425E-03, -6.90460016972063023E-05, \
  0., 0., 0.],
 ['SATURN',2.85885980666130812E-04,3.,0.70, \
  8.34336671824457987E+00,  4.12479856412430479E+00, -4.03523417114321381E-01, \
 -2.76742510726862411E-03,  4.99852801234917238E-03,  2.30417297573763929E-05, \
  0., 0., 0.],
 ['URANUS',4.36624404335156298E-05,3.,1.30, \
  1.28943695621391310E+01, -1.51111514016986312E+01, -2.23307578892655734E-01, \
  2.96460137564761618E-03,  2.37847173959480950E-03, -2.96589568540237556E-05, \
  0., 0., 0.],
 ['NEPTUNE',5.15138902046611451E-05,3.,1.76, \
  1.53796971148509165E+01, -2.59193146099879641E+01,  1.79258772950371181E-01, \
  2.68067772490389322E-03,  1.62824170038242295E-03, -9.51592254519715870E-05, \
  0., 0., 0.],
 ['PLUTO',7.39644970414201173E-09,3.,1.1, \
 -1.15095623952731607E+01, -2.70779438829451422E+01,  6.22871533567077229E+00, \
  2.97220056963797431E-03, -1.69820233395912967E-03, -6.76798264809371094E-04, \
  0., 0., 0.]]

def _planets_only(define_mercury_attributes = False):
    data = numpy.array([tuple(entry) for entry in _solsysdat], dtype=[('name','S10'), 
        ('mass','<f8'), ('celimit','<f8'), ('density','<f8'), 
        ('x','<f8'), ('y','<f8'), ('z','<f8'), 
        ('vx','<f8'), ('vy','<f8'), ('vz','<f8'), 
        ('Lx','<f8'), ('Ly','<f8'), ('Lz','<f8')])
    
    planets = Particles(len(_solsysdat))
    planets.name = data['name']
    planets.mass = units.MSun.new_quantity(data['mass'])
    density = (units.g/units.cm**3).new_quantity(data['density'])
    planets.radius = ((planets.mass/density) ** (1/3.0)).as_quantity_in(units.km)
    for attribute in ['x', 'y', 'z']:
        setattr(planets, attribute, units.AU.new_quantity(data[attribute]))
    for attribute in ['vx', 'vy', 'vz']:
        setattr(planets, attribute, units.AUd.new_quantity(data[attribute]).as_quantity_in(units.km / units.s))
    
    if define_mercury_attributes:
        planets.density = density
        angular_momentum_unit = units.MSun * units.AU**2/units.day
        for attribute in ['Lx', 'Ly', 'Lz']:
            setattr(planets, attribute, angular_momentum_unit.new_quantity(data[attribute]).as_quantity_in(units.J * units.s))
        planets.celimit = units.none.new_quantity(data['celimit'])
    
    return planets

def new_solar_system_for_mercury():
    """
    Create initial conditions for the symplectic integrator Mercury, describing 
    the solar system. Returns a tuple consisting of two particle sets. The first 
    set contains the central particle (sun) and the second contains the planets 
    and Pluto (the 'orbiters'). The positions and velocities are in heliocentric 
    coordinates.
    
    Defined attributes sun: 
    name, mass, radius, j2, j4, j6, Lx, Ly, Lz
    
    Defined attributes orbiters: 
    name, mass, radius, density, x, y, z, vx, vy, vz, Lx, Ly, Lz, celimit
    """
    planets = _planets_only(define_mercury_attributes = True)
    centre = Particles(1)
    centre.name = 'SUN'
    centre.mass = 1.0 | units.MSun
    centre.radius = 0.0000001 | units.AU
    centre.j2 = .0001|units.AU**2
    centre.j4 = .0|units.AU**4
    centre.j6 = .0|units.AU**6
    centre.angularmomentum = [0.0, 0.0, 0.0] | units.MSun * units.AU**2/units.day
    return centre, planets

def new_solar_system():
    """
    Create initial conditions describing the solar system. Returns a single 
    particle set containing the sun, planets and Pluto. The model is centered at 
    the origin (center-of-mass(-velocity) coordinates).
    
    Defined attributes: 
    name, mass, radius, x, y, z, vx, vy, vz
    """
    sun = Particle()
    sun.name = 'SUN'
    sun.mass = 1.0 | units.MSun
    sun.radius = 1.0 | units.RSun
    planets = _planets_only()
    
    particles = Particles()
    particles.add_particle(sun)
    particles.add_particles(planets)
    particles.move_to_center()
    return particles
