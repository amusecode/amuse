from amuse.community.arepo import Arepo
from amuse.units import nbody_system
# import unit seconds as s

# Check code runs without errors
x = Arepo(redirection="none")
x.initialize_code()
print(x.get_position(21057))  # This is the ID of P[0]
print(x.get_position(21060))

#print("Evolving")
#x.evolve_model(0.00001)
#print("Evolving another step")
#x.evolve_model(0.00002)

#x.run_sim()
#END_TIME = 1.0 | s
#x.evolve_model(END_TIME)

# x.(evolve for a single timestep)

# amuse tests to check for diverging behaviour
