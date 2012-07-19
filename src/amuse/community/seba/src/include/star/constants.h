
       //=======================================================//    _\|/_
      //  __  _____           ___                    ___       //      /|\
     //  /      |      ^     |   \  |         ^     |   \     //          _\|/_
    //   \__    |     / \    |___/  |        / \    |___/    //            /|\
   //       \   |    /___\   |  \   |       /___\   |   \   // _\|/_
  //     ___/   |   /     \  |   \  |____  /     \  |___/  //   /|\     
 //                                                       //            _\|/_
//=======================================================//              /|\

/*
 *  constants.h: header file for constants.
 *           
 *.............................................................................
 *    version 1:  May 1993   Simon F. Portegies Zwart
 *    version 2:
 *.............................................................................
 *     This file includes:
 *  1) definition of physical constants, assumed numbers and pre set
 *     numbers.
 *
 *.............................................................................
 */
#ifndef     _CONSTANTS
#   define  _CONSTANTS

#include "stdinc.h"       

#define static_cast(type, expr) 	(type)(expr)
#define const_cast(type, expr) 		(type)(expr)
#define reintepret_cast(type, expr) 	(type)(expr)
#define dynamic_cast(type, expr) 	(type)(expr)

#define SEBA_VERSION 2.0

#define TRUE 1
#define FALSE 0

// Stellar remnant constants.
// Remnant lifetimes.
//#define  SPI_TIME	0.0005		     // Spiral in duration time.
#define  POST_AGB_TIME  0.005                // time to lose giant envelope.

//Assymetric supernovae kicks
#define MEAN_KICK 450.0			// km/s
#define VAR_KICK  90.0			// km/s

//#define BINARY_UPDATE_TIMESTEP 0.9      // Timestep increase for binary   
#define BINARY_UPDATE_TIMESTEP 0.1      // Timestep increase for binary   
#define TEN_PERCENT 0.1
//#define STABLE_ZETA 15
//#define J_SPI_FR 0.5			//Pols 93 advises 1/3

// Model choices.
//#define exp_HE_lifetime FALSE		// exponential he star lifetime fits.

enum binary_history  {ms_ms = 0, bs_ms, bs_bs, he_ms, heN_ms, 
                      wr_ms, he_he, rscvn, 
		      wuma, wdxb, 
                      lmxb, mmxb, hmxb, spi, spi_2, 
                      wd_ms, wd_he, wd_wd, wd_ns, wdXns, 
                      pr_ms, pr_he, pr_st, prXst, pr_pr, prXpr,
                      ns_ms, ns_he, ns_st, nsXst, ns_ns, nsXns, st_st, 
                      no_of_type};

enum Be_binaries    {be_ms=0, be_he, be_wd, be_ns, no_of_be_binary};

enum mathematical_constant {one_third,
			    two_third,
			    pi,
			    two_pi
};

enum physics_constants{G, gravitational_constant,
		       C, speed_of_light,
		       Myear, million_years,
		       days, seconds_per_day,
		       km_per_s, kilometer_per_second,
		       kilometer, kilometer_in_centimeters,
		       nucleair_efficiency
};

enum astronomical_scale_parameter {PC, parsec,
				   AU, astronomical_unit
				};

enum solar_parameter {Msun, solar_mass,
		      Rsun, solar_radius,
		      Lsun, solar_luminosity,
		      Tsun, solar_temperature,
		      Zsun, solar_metalicity,
		      energy_to_mass_in_internal_units
};

enum pulsar_initial_conditions {pulsar_magnetic_field,
				pulsar_pulse_period,
				kanonical_neutron_star_radius,
				kanonical_neutron_star_mass,
				maximum_neutron_star_mass,
				minimum_neutron_star_mass
};
    
enum stellar_mass_limits {low_mass_star_mass_limit,
			  medium_mass_star_mass_limit,
			  massive_star_mass_limit,
			  upper_ZAMS_mass_for_degenerate_core,
			  minimum_main_sequence,
			  helium_dwarf_mass_limit,
			  carbon_dwarf_mass_limit,
			  Chandrasekar_mass,
//			  neutron_star,
//			  minimum_neutron_star,
//			  maximum_neutron_star,
			  helium2neutron_star,
			  COcore2black_hole,
			  super_giant2neutron_star,
			  super_giant2black_hole,
			  maximum_main_sequence,
			  minimum_helium_star,
			  helium_star_lifetime_fraction,
			  helium_star_final_core_fraction,
};

enum boolean_parameter {hyper_critical,
			super_giant_disintegration,
			proto_star_to_binary,
			impulse_kick_for_black_holes,
			use_angular_momentum_tidal,
			use_common_envelope_gamma_gamma,
			use_common_envelope_alpha_alpha
		     };

enum super_nova_kick_distribution {internally_decided_velocity_kick,
                                   no_velocity_kick,
				   Maxwellian_velocity_kick,
				   Paczynski_velocity_kick,
				   delta_function_velocity_kick
                                  };

enum accretion_parameter {black_hole_accretion_limit,
			  neutron_star_accretion_limit,
			  white_dwarf_accretion_limit,
			  thermo_nuclear_flash
		      };

enum model_parameter {star_formation_efficiency,
		      star_formation_timescale,
		      magnetic_mass_limit,
		      magnetic_braking_exponent,
		      corotation_eccentricity,
		      tidal_circularization_radius,
		      core_overshoot,
		      hydrogen_fraction,
		      common_envelope_efficiency,
		      envelope_binding_energy,
		      specific_angular_momentum_loss,
		      dynamic_mass_transfer_gamma,
		      non_massive_star_envelope_fraction_lost,
		      massive_star_envelope_fraction_lost,
		      relaxation_driven_mass_loss_constant,
		      time_independent_mass_loss_law,
		      massive_star_mass_loss_law,
		      Darwin_Riemann_instability_factor,
		      homogeneous_sphere_gyration_radius_sq,
		      radiative_star_gyration_radius_sq,
		      convective_star_gyration_radius_sq,
		      rejuvenation_exponent,
		      spiral_in_time,
                     };

enum observational_parameter {B_emission_star_mass_limit,
			      Barium_star_mass_limit,
			      Blue_straggler_mass_limit
                     };


enum safety_parameter {timestep_factor,
		       maximum_binary_update_time_fraction,
		       minimum_timestep,
		       minimum_mass_step,
		       maximum_timestep,
		       maximum_recursive_calls,
               tiny,
                      };

enum dynamics_update_parameter {stellar_mass_update_limit,
				semi_major_axis_update_limit,
   			        binary_update_time_fraction
				};

static
class stellar_evolution_constants {  // Easy to have a name for compiling.
  public:

  real mathematics(mathematical_constant);
  real physics(physics_constants);

  real super_nova_kick(super_nova_kick_distribution pk =
		       internally_decided_velocity_kick,
		       const real v_disp = 600);
  
  real parameters(solar_parameter);
  real parameters(astronomical_scale_parameter);
  real parameters(pulsar_initial_conditions);
  real parameters(stellar_mass_limits);
  
  int use_common_envelope_method();

  bool parameters(boolean_parameter);
  real parameters(accretion_parameter);
  real parameters(model_parameter);

  real parameters(observational_parameter);

  real safety(safety_parameter);

  real star_to_dyn(dynamics_update_parameter);

} cnsts;

static
class stellar_model_constants {  // Easy to have a name for compiling.

  public:
  real a(int, real);
  real b(int, real);
  real c(int, real);
} smc;


#endif		// _CONSTANTS




