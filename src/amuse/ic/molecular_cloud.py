"""
Generate a molecular cloud of SPH particles
"""

import warnings
from amuse.ext.molecular_cloud import molecular_cloud


def new_molecular_cloud(
    nf=32,
    power=-3.0,
    target_number_of_stars=None,
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
            "targetN is deprecated, use target_number_of_stars instead",
            category=FutureWarning,
        )
        if target_number_of_stars != targetN:
            raise ValueError(
                "targetN and target_number_of_stars have different values, "
                "this is only allowed if targetN is None and target_number_of_stars "
                "is not None"
            )
        targetN = target_number_of_stars
    if target_number_of_stars is None:
        raise ValueError("target_number_of_stars must be set")
    return molecular_cloud(
        nf=nf,
        power=power,
        targetN=targetN,
        ethep_ratio=ethep_ratio,
        convert_nbody=convert_nbody,
        ekep_ratio=ekep_ratio,
        seed=seed,
        base_grid=base_grid,
    ).result
