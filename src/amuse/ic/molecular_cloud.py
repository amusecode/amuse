"""
Generate a molecular cloud of SPH particles
"""

from amuse.ext.molecular_cloud import molecular_cloud

def new_molecular_cloud(
    nf=32, power=-3., targetN=10000, ethep_ratio=0.01,
    convert_nbody=None, ekep_ratio=1., seed=None, base_grid=None,
):
    return molecular_cloud(
        nf=nf, power=power, targetN=targetN, ethep_ratio=ethep_ratio,
        convert_nbody=convert_nbody, ekep_ratio=ekep_ratio, seed=seed,
        base_grid=base_grid,
    ).result
