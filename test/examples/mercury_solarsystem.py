import numpy
from amuse.legacy.mercury.interface import Mercury

try:
    from matplotlib import pyplot
    HAS_MATPLOTLIB = True
except ImportError:
    HAS_MATPLOTLIB = False

solarsystem= \
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

def planetplot():
  instance=Mercury()
  instance.initialize_code()
  ids=dict()
  xpos=dict()
  ypos=dict()
  for x in solarsystem:
    pid,err=instance.new_orbiter(x[1],x[3], \
     x[4],x[5],x[6],x[7],x[8],x[9],x[10],x[11],x[12],x[2])
    ids[x[0]]=pid
    xpos[x[0]]=[]
    ypos[x[0]]=[]
    xpos[x[0]].append( x[4])
    ypos[x[0]].append( x[5])

  instance.commit_particles()

  t_end=365.25*20
  time=0

  while time<t_end:
    time=time+8
    err=instance.evolve(time)
    for p in solarsystem:
      mass,dens,x,y,z,vx,vy,vz,sx,sy,sz,celimit,err= \
        instance.get_orbiter_state(ids[p[0]])
      xpos[p[0]].append(x)
      ypos[p[0]].append(y)

  instance.stop()

  pyplot.plot(xpos['MERCURY'],ypos['MERCURY'])
  pyplot.plot(xpos['VENUS'],ypos['VENUS'])
  pyplot.plot(xpos['EARTHMOO'],ypos['EARTHMOO'])
  pyplot.plot(xpos['MARS'],ypos['MARS'])
  pyplot.plot(xpos['JUPITER'],ypos['JUPITER'])
  pyplot.xlim(-6.0, 6.0)
  pyplot.ylim(-6.0, 6.0)
  pyplot.savefig('solarsytem.png')

def energyplot():
  instance=Mercury()
  instance.initialize_code()
  for x in solarsystem:
    pid,err=instance.new_orbiter(x[1],x[3], \
     x[4],x[5],x[6],x[7],x[8],x[9],x[10],x[11],x[12],x[2])
  instance.commit_particles()

  t_end=365.25*1000000
  time=0

  e0,err=instance.get_total_energy()
  energy=[0.]
  times=[0.]
  while time<t_end:
    time=time+8*10000.
    err=instance.evolve(time)
    if(err!=0):
      print "err:", err
    e,err=instance.get_total_energy()
    t,err=instance.get_time()
    energy.append((e0-e)/e0)
    times.append(t/365.25)
  
  instance.stop()

  pyplot.plot(times,energy)
  pyplot.xlim(0,t_end/365.25)
  pyplot.ylim(-1.e-6,1.e-6)
  pyplot.savefig('solarsytem_energy.png')
  
if __name__=='__main__':
  print "this may take a while.."
  energyplot()
