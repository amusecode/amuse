from amuse.support.data import core
from amuse.support.units import nbody_system
from amuse.support.units import units
import numpy as np

#mercury orbiter:
#mass, radius, x, y, z, vx, vy, vz, Lx, Ly, Lz, celimit

solsysdat= \
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

class solarsystem(object):
    def __init__(self):
        pass

    def new_solarsystem(self):
        planets = core.Particles(9)
        for i, d in enumerate(solsysdat):
            planets[i].mass = d[1] | units.MSun
            planets[i].radius = 1e-12 | units.AU #dummy
            planets[i].density = d[3] | units.g/units.cm**3
            planets[i].x = d[4] | units.AU
            planets[i].y = d[5] | units.AU
            planets[i].z = d[6] | units.AU
            planets[i].vx = d[7] | units.AUd
            planets[i].vy = d[8] | units.AUd
            planets[i].vz = d[9] | units.AUd
            planets[i].Lx = d[10] | units.MSun * units.AU**2/units.day
            planets[i].Ly = d[11] | units.MSun * units.AU**2/units.day
            planets[i].Lz = d[12] | units.MSun * units.AU**2/units.day
            planets[i].celimit = d[2] | units.none

        centre = core.Particles(1)
        centre.mass = 1.0 | units.MSun
        centre.radius = 0.01 | units.AU
        centre.j2 = .0001|units.AU**2
        centre.j4 = .0|units.AU**4
        centre.j6 = .0|units.AU**6
        centre.angularmomentum = [0.0, 0.0, 0.0] | units.MSun * units.AU**2/units.day

        return centre, planets
