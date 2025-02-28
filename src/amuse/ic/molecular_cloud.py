"""
Generate a molecular cloud of SPH particles
"""

import warnings
from amuse.ext.molecular_cloud import molecular_cloud


def new_molecular_cloud(
    nf=32,
    power=-3.0,
    target_number_of_particles=None,
    ethep_ratio=0.01,
    convert_nbody=None,
    ekep_ratio=1.0,
    seed=None,
    base_grid=None,
    targetN=None,
):
    """
    Creates a molecular cloud of SPH particles
    """
    if targetN is not None:
        warnings.warn(
            "targetN is deprecated, use target_number_of_particles instead",
            category=FutureWarning,
        )
        if target_number_of_particles is not None and target_number_of_particles != targetN:
            raise ValueError(
                "targetN and target_number_of_particles have different values, "
                "this is only allowed if one of them is None"
            )
        target_number_of_particles = targetN
    if target_number_of_particles is None:
        raise ValueError("target_number_of_particles must be set")
    return molecular_cloud(
        nf=nf,
        power=power,
        targetN=target_number_of_particles,
        ethep_ratio=ethep_ratio,
        convert_nbody=convert_nbody,
        ekep_ratio=ekep_ratio,
        seed=seed,
        base_grid=base_grid,
    ).result
