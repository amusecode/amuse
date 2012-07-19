#ifndef    STDFUNC
#   define STDFUNC

#include "stdinc.h"
#include "starlab_constants.h"
#include "star_support.h"

// Standard functions for various implementations.

//		Tidal capture routines.
       real tf2_energy_diss(const real eta, const stellar_type);
       real tf3_energy_diss(const real eta, const stellar_type);
       real lineair_interpolation(const real, const real, const real,
				  const real, const real);

//		Super nova utilities.
       real post_sn_cm_velocity(const real, const real, const real,
                                const real, const real, const real,
                                const real,
                                const real, const real, const real);

       real post_sn_semi_major_axis(const real, const real, const real,
                                    const real, const real, const real,
                                    const real,
                                    const real, const real, const real);
       real post_sn_eccentricity(const real, const real, const real,
                                 const real, const real, const real, 
                                 const real,
                                 const real, const real, const real);
       real random_angle(const real, const real);
       real eccentric_anomaly(const real, const real);
       real random_eccentric_anomaly(const real);
       real random_separation(const real, const real);
       real paczynski_distribution(const real, const real);
       real random_paczynski_velocity(const real);
       real maxwellian(const real, const real);
       real random_maxwellian_velocity(const real);
       real gravitational_focussed_velocity(const real, const real,
                                            const real, const real,
                                            const real);
       real gf_velocity_distribution_extremum(const real, const real,
                                            const real, const real);
//       real random_focussed_maxwellian_velocity(const real, const real);
       real random_focussed_maxwellian_velocity(const real, const real,
					        const real, const real,
						const real);
       real gauss(const real, const real);
       real gauss();
       real cross_section(const real, const real, const real,  const real);

       real eddington_limit(const real radius,
			    const real dt,
			    const real mu=1);

       real kinetic_energy(const real, const real);
       real potential_energy(const real, const real, const real);
       real velocity_at_infinity(const real,
                                 const real,
                                 const real,
                                 const real);
       void fool(char *);

//		Timescale utilities.
       real turn_off_mass(const real);
       real main_sequence_time(const real);
       real zero_age_main_sequnece_radius(const real);
       real kelvin_helmholds_timescale(const real, const real, const real);
       real nucleair_evolution_timescale(const real, const real);
       real dynamic_timescale(const real, const real);

       real roche_radius(const real, const real, const real);
//		Integration utilities

#endif
