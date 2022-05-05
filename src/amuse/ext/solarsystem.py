import numpy

from amuse.units import units, nbody_system, constants
from amuse.datamodel import Particles, Particle


pi_over_180 = numpy.pi/180.

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
    planets.name = list(data['name'])
    print(planets.name.dtype)
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
    centre.angular_momentum = [0.0, 0.0, 0.0] | units.MSun * units.AU**2/units.day
    return centre, planets


def new_kepler():
  from amuse.community.kepler.interface import Kepler
  converter = nbody_system.nbody_to_si(1|units.MSun,1|units.AU)
  kepler = Kepler(converter)
  kepler.initialize_code()
  return kepler

def get_position(mass_sun, mass_planet, ecc, semi, mean_anomaly, incl, argument, longitude, delta_t=0.|units.day):
  """
  cartesian position and velocity from orbital elements,
  where the orbit is evolved from given mean_anomaly 
  by time delta_t
  argument -- argument of perihelion
  longitude -- longitude of ascending node
  """
  kepler = new_kepler()
  kepler.initialize_from_elements(mass=(mass_sun+mass_planet),
                                  semi=semi,
                                  ecc=ecc,
                                  mean_anomaly=mean_anomaly)
  kepler.transform_to_time(time=delta_t)
  r = kepler.get_separation_vector()
  v = kepler.get_velocity_vector()
  
  kepler.stop()
  
  a1 = ([numpy.cos(longitude), -numpy.sin(longitude), 0.0], [numpy.sin(longitude), numpy.cos(longitude), 0.0], [0.0, 0.0, 1.0])
  a2 = ([1.0, 0.0, 0.0], [0.0, numpy.cos(incl), -numpy.sin(incl)], [0.0, numpy.sin(incl), numpy.cos(incl)])
  a3 = ([numpy.cos(argument), -numpy.sin(argument), 0.0], [numpy.sin(argument), numpy.cos(argument), 0.0], [0.0, 0.0, 1.0])
  A = numpy.dot(numpy.dot(a1,a2),a3)
  r_vec = numpy.dot(A,numpy.reshape(r,3,'F'))
  v_vec = numpy.dot(A,numpy.reshape(v,3,'F'))
  
  # for relative vectors
  r[0] = r_vec[0]
  r[1] = r_vec[1]
  r[2] = r_vec[2]
  v[0] = v_vec[0]
  v[1] = v_vec[1]
  v[2] = v_vec[2]
  
  return r,v

def get_sun_and_planets(delta_JD=0.|units.day):
  """
  eight planets of the Solar System
  as for JD = 2457099.500000000 = A.D. 2015-Mar-18 00:00:00.0000 (CT)
  http://ssd.jpl.nasa.gov/horizons.cgi
  """
  planets = Particles(8)
  
  # mass
  planets.mass = [3.302e23,
                  48.685e23,
                  5.97219e24,
                  6.4185e23,
                  1898.13e24,
                  5.68319e26,
                  86.8103e24,
                  102.41e24] | units.kg
  #radius
  planets.radius = [2439.7,
                    6051.8,
                    6378.1,
                    3396.2,
                    71492,
                    60268,
                    25559,
                    24764] | units.km
  
  # eccentricity
  planets_ecc = [2.056263501026885E-01,
                 6.756759719005901E-03,
                 1.715483324953308E-02,
                 9.347121362500883E-02,
                 4.877287772914470E-02,
                 5.429934603664216E-02,
                 4.911406962716518E-02,
                 8.494660388602767E-03]
  
  # semi-major axis
  planets_semi = [3.870989725156447E-01,
                  7.233252880006816E-01,
                  1.000816989613834E+00,
                  1.523624142457679E+00,
                  5.203543088590996E+00,
                  9.547316304899041E+00,
                  1.915982879739036E+01,
                  2.997013749028780E+01] | units.AU
  
  # mean anomaly [degrees]
  planets_mean_anomaly = [2.256667460183225E+02,
                          3.096834722926926E+02,
                          6.970055236286768E+01,
                          5.013506750245609E+01,
                          1.213203242081277E+02,
                          1.423311616732398E+02,
                          2.079860620353052E+02,
                          2.712246916734600E+02]
  planets_mean_anomaly = numpy.array(planets_mean_anomaly) * pi_over_180
  
  # inclination [IN degrees]
  planets_inclination = [7.004026765179669E+00,
                         3.394480103844425E+00,
                         3.563477431351056E-03,
                         1.848403408106458E+00,
                         1.303457729562742E+00,
                         2.488017444885577E+00,
                         7.728000142736371E-01,
                         1.767720502209091E+00]
  planets_inclination = numpy.array(planets_inclination) * pi_over_180
  
  # Longitude of Ascending Node [OM degrees]
  planets_longitude = [4.831163083479358E+01,
                       7.663982595051040E+01,
                       1.775515437672556E+02,
                       4.951282677064384E+01,
                       1.005036717671826E+02,
                       1.135683875842263E+02,
                       7.388411509910506E+01,
                       1.317497218434830E+02]
  planets_longitude = numpy.array(planets_longitude) * pi_over_180
  
  # Argument of Perihelion [W degrees]
  planets_argument = [2.916964171964058E+01,
                      5.469102797401222E+01,
                      2.877495001117996E+02,
                      2.865420083537150E+02,
                      2.740725976811202E+02,
                      3.398666856578898E+02,
                      9.666856264946740E+01,
                      2.951871807292030E+02]
  planets_argument = numpy.array(planets_argument) * pi_over_180
  
  planets.name = ['Mercury',
                  'Venus',
                  'Earth',
                  'Mars',
                  'Jupiter',
                  'Saturn',
                  'Uranus',
                  'Neptune']
  
  ### to compare with JPL, mass of the Sun needs to be rescaled
  #mg_nasa = 1.32712440018e20 | (units.m**3 / units.s**2)
  #g_nasa = 6.67259e-11 | (units.m**3 / units.kg / units.s**2)
  #ms = mg_nasa / g_nasa
  
  sun = Particle()
  sun.name = 'Sun'
  #sun.mass = ms
  sun.mass = 1.0 | units.MSun
  sun.position = [0.,0.,0.] | units.AU
  sun.velocity = [0.,0.,0.] | units.kms
  
  # get the position and velocity vectors relative to sun 
  # by evolving in Kepler
  for i,ecc_i in enumerate(planets_ecc):
    r, v = get_position(sun.mass,
                        planets[i].mass,
                        planets_ecc[i],
                        planets_semi[i],
                        planets_mean_anomaly[i],
                        planets_inclination[i],
                        planets_longitude[i],
                        planets_argument[i],
                        delta_t=delta_JD)
    planets[i].position = r
    planets[i].velocity = v
    
  return sun, planets

def solar_system_in_time(time_JD=2457099.5|units.day):
  """
  Initial conditions of Solar system --
  particle set with the sun + eight planets,
  at the center-of-mass reference frame.

  Defined attributes: 
  name, mass, radius, x, y, z, vx, vy, vz
  """
  time_0 = 2457099.5 | units.day
  delta_JD = time_JD-time_0
  sun, planets = get_sun_and_planets(delta_JD=delta_JD)
  
  solar_system = Particles()
  solar_system.add_particle(sun)
  solar_system.add_particles(planets)
  solar_system.move_to_center()
  
  ### to compare with JPL, relative positions and velocities need to be corrected for the
  # Sun's vectors with respect to the barycenter
  #r_s = (3.123390770608490E-03, -4.370830943817017E-04, -1.443425433116342E-04) | units.AU
  #v_s = (3.421633816761503E-06,  5.767414405893875E-06, -8.878039607570240E-08) | (units.AU / units.day)
  #print sun
  #print planets.position.in_(units.AU) + r_s
  #print planets.velocity.in_(units.AU/units.day) + v_s
  
  return solar_system

def old_new_solar_system():
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

def new_solar_system(Julian_date=-1|units.day):
    if Julian_date<0|units.day:
        return old_new_solar_system()
    else:
        return solar_system_in_time(Julian_date)

def new_option_parser():
    from amuse.units.optparse import OptionParser
    result = OptionParser()
    result.add_option("-d", dest="Julian_date", unit=units.day,
                      type=float, default = 2438871.5|units.day,
                      help="julian date [%default]")
    return result

if __name__ in ('__main__', '__plot__'):
  o, arguments  = new_option_parser().parse_args()
  solar_system = new_solar_system(o.Julian_date)
  print(solar_system)


