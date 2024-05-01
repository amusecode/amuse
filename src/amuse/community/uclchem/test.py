from interface import Uclchem
species='H H2'
dict = "{'outspecies': 2}"
chem = Uclchem()
test = chem.sim_cloud(outSpeciesIn=species, dictionary=dict)
print(test)