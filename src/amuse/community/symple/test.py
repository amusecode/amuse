from amuse.lab import *
from interface import symple
from matplotlib import pyplot as plt
import sys, math

n = 2

g = symple(redirection='none')
g.initialize_code()
g.parameters.set_defaults()
print g.parameters

# Make a binary.

p = new_plummer_model(n)
r = 0.5|nbody_system.length
vcirc = (0.25*nbody_system.G*p[0].mass/r).sqrt()
v = 0.5*vcirc
p[0].position = [r, zero, zero]
p[1].position = -p[0].position
p[0].velocity = [zero, v, zero]
p[1].velocity = -p[0].velocity

g.particles.add_particles(p)
g.commit_particles()

E0 = g.kinetic_energy + g.potential_energy
M = p[0].mass+p[1].mass
GM = (nbody_system.G)*M
a = -0.5*nbody_system.G*p[0].mass*p[1].mass/E0
P = 2*math.pi*(a**3/GM).sqrt()
L = 2*p[0].mass*(p[0].x*p[0].vy)
print 'L =', L
mu = p[0].mass*p[1].mass/M
print 'mu =', mu
ee = E0/mu
print 'ee =', ee
h = L/mu
print 'h =', h
e = math.sqrt(1+2*ee*h**2/GM**2)
print 'E =', E0
print 'a =', a
print 'P =', P
print 'e =', e

# Include mass loss.

mdotfac = 0.
for i in range(n):
    g.set_dmdt(i, -mdotfac*p[i].mass/(1.e3|nbody_system.time))

g.parameters.integrator = 10
g.parameters.timestep = P/128
g.parameters.eta = 0.25			# timestep and eta are mutually exclusive
print g.parameters

time = 0.|nbody_system.time
dt = P/64
t_end = 100*P

t = [time.number]
x = [(g.particles[0].x-g.particles[1].x).number]
y = [(g.particles[0].y-g.particles[1].y).number]
e = [0.]
while time < t_end:
    time += dt
    g.evolve_model(time)
    t.append(time.number)
    x.append((g.particles[0].x-g.particles[1].x).number)
    y.append((g.particles[0].y-g.particles[1].y).number)
    E1 = g.kinetic_energy + g.potential_energy
    dE = E1/E0 - 1
    e.append(dE)
    
print 'dE/E =', dE
print ''

if 0:
    plt.figure(figsize=(6,6))
    plt.plot(x, y)
    plt.xlim(-1.2, 1.2)
    plt.ylim(-1.2, 1.2)
    plt.axes().set_aspect('equal')
else:
    plt.plot(t, e)
    
plt.show()

g.stop()
