
import numpy
from amuse.community.mercury.interface import Mercury


try:
    from matplotlib import pyplot
    HAS_MATPLOTLIB = True
except ImportError:
    HAS_MATPLOTLIB = False


instance=Mercury()
instance.initialize_code()
mass=3.04043264264672381E-06
dens=5.52
x=2.42093942183383037E-01
y=-9.87467766698604366E-01
z=-4.54276292555233496E-06
vx=1.64294055023289365E-02
vy=4.03200725816140870E-03
vz=1.13609607260006795E-08
sx=sy=sz=0.
celimit=20.
pid,err=instance.new_orbiter(mass,dens,x,y,z,vx,vy,vz,sx,sy,sz,celimit)  
instance.commit_particles()

t_end=365.25*1000
time=0
xx=[x]
yy=[y]

while time<t_end:
    time=time+800
    err=instance.evolve(time)
    mass,dens,x,y,z,vx,vy,vz,sx,sy,sz,celimit,err=instance.get_orbiter_state(pid)
    xx.append(x)
    yy.append(y)

instance.stop()
