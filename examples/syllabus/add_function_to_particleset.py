from amuse.lab import *

def function(body):
    return 0.2*body[0].mass

bodies = Particles(1)
bodies.mass = 1|units.MSun
bodies.add_global_function_attribute("new_function", function)
print(bodies.new_function().in_(units.MSun))
