#include "constants.h"
#include "stdfunc.h"

// constructor for constant class.
// Currently there are not functionalities.

//stellar_evolution_constants::stellar_evolution_constants() {
//  //Dummy function in order to allow compilation.
//}


real stellar_evolution_constants::mathematics(mathematical_constant pm) {

    switch(pm) {
	case one_third:                        return 0.33333333333333333333;
             break;
	case two_third:                        return 0.66666666666666666666;
             break;
	case pi:                               return 3.14159265358979323846;
	  break;                               // of the orde of 1
	case two_pi:                           return 6.28318530717958647692;
             break;
        default:
	     cerr << "\nNo recognized option in "
		     "stellar_evolution_constants::"
		     "parameters(mathematical_constant "
	          << pm << ")"
		  << endl;
             exit(1);
    }
}

real stellar_evolution_constants::physics(physics_constants pp) {

    // physics parameters 
    switch(pp) {
	case gravitational_constant:
	case G:                             return 6.67e-8;    // [cgs]
             break;
        case speed_of_light:					    
        case C:                             return 2.9979e+10; // [cm/s]
             break;
        case million_years:
        case Myear:                         return 3.15e+13;   // [s]
             break;
        case seconds_per_day:
        case days:                          return 8.6400e+4;   // [s]
             break;
        case kilometer_per_second:
        case km_per_s:                      return 1.0e+5;      // [cm/s]
             break;
        case kilometer_in_centimeters:
        case kilometer:                     return 1.0e+5;   // [cm]
             break;
	case nucleair_efficiency:                      return 0.007;
	      break;                          // nuc. energy production eff
	                                      // Delta E = 0.007 Mc^2
        default:
	     cerr << "\nNo recognized option in "
		     "stellar_evolution_constants::"
		     "parameters(physics_constants "
	          << pp << ")"
		  << endl;
             exit(1);
    }
}

real stellar_evolution_constants::super_nova_kick(
				  super_nova_kick_distribution pk,
				  const real v_disp){
  
  // Super Nova Kick functions.
  // kick velicity imparted to the remnant during the supernova event.
  // changing the kick velocity requires recompilation.

  //	Paczynsky's velocity distribution with a dispersion
  //	of 600 km/s is prefered by Hartman.

  //	Gauss kick velocity distribution.
  //	Lyne, A.G., \& Lorimer, D.R., 1994, Nature 369, 127
  //	propose v_k = 450 +- 90 km/s
  //	this is comparable to a one dimensional gaussion
  //	of 260 km/s.

  //    Hansen & Phinney (1997, MNRAS 291, 569) kick velocity distribution.
  //    They prever a Maxwellian with a velocity dispersion of
  //    270 km/s.
  
  // selected kick distribution imparted to a newly formed neutron star
  // in a supernova. 
    switch(pk) {
        case internally_decided_velocity_kick:
        case no_velocity_kick:              return 0;
             break;                                
	case Maxwellian_velocity_kick:      return
					    random_maxwellian_velocity(v_disp);
             break;
	case Paczynski_velocity_kick:       return
					    random_paczynski_velocity(v_disp);
             break;
	case delta_function_velocity_kick:  return v_disp;
             break;
    default:
	     cerr << "\nNo recognized option in "
		     "stellar_evolution_constants::"
		     "super_nova_kick(super_nova_kick_distribution "
	          << pk << ")"
		  << endl;
             exit(1);
    }
}  

real stellar_evolution_constants::parameters(astronomical_scale_parameter pa) {

    // Astronomical distance scale parameters in centimeters [cm]
    switch(pa) {
	case PC:
	case parsec:                             return 3.0857e+18;
             break;                                
	case AU:                      
	case astronomical_unit:                  return 1.496e+13; 
             break;                                
        default:
	     cerr << "\nNo recognized option in "
		     "stellar_evolution_constants::"
		     "parameters(astronomical_scale_parameter "
	          << pa << ")"
		  << endl;
             exit(1);
    }
}

real stellar_evolution_constants::parameters(solar_parameter ps) {

    // Solar parameters in cgs.
    switch(ps) {
	case solar_mass:
        case Msun:                      return 1.989e+33; //[gram]
             break;                                
	case solar_radius:
	case Rsun:                      return 6.96e+10;//[centimeter]
             break;                                
	case solar_luminosity:
	case Lsun:                      return 3.862e+33; // [erg/s]
             break;                                
	case solar_temperature:
	case Tsun:                      return 5770;     // [Kelvin]
             break;
        case Zsun:
        case solar_metalicity:          return 0.02;
	     break;
        case energy_to_mass_in_internal_units:         return 
				   cnsts.physics(nucleair_efficiency)
				 * pow(cnsts.physics(C), 2)
                                 * cnsts.parameters(solar_mass)
                                 / (cnsts.parameters(solar_luminosity) *
				    cnsts.physics(Myear));
	     break;
        default:
	     cerr << "\nNo recognized option in "
		     "stellar_evolution_constants::"
		     "parameters(solar_parameter "
	          << ps << ")"
		  << endl;
             exit(1);
    }
}

real stellar_evolution_constants::parameters(pulsar_initial_conditions pp) {

    // Netron star basic parameters.
    switch(pp) {
	case pulsar_magnetic_field:              return 12;     // [log Gs]
             break;                                
	case pulsar_pulse_period:                return 0.1;    // [s]
             break;                                
        case kanonical_neutron_star_radius:      return 1.5e-5; // [cm]
             break;                               
        case kanonical_neutron_star_mass:        return 1.34;   // [msun]
             break;                                
	case maximum_neutron_star_mass:          return 1.5;    // [Msun]
             break;                                
	case minimum_neutron_star_mass:          return 0.0925; // [Msun]
             break;                              // Newtonian polytrope 
        default:
	     cerr << "\nNo recognized option in "
		     "stellar_evolution_constants::"
		     "parameters(pulsar_initial_conditions "
	          << pp << ")"
		  << endl;
             exit(1);
    }
}

real stellar_evolution_constants::parameters(stellar_mass_limits pm) {

    // Switch statement for certain limits. All values are in
    // solar units [Msun].
    switch(pm) {
	case low_mass_star_mass_limit:           return 1.5;
             break;                                
	case medium_mass_star_mass_limit:        return 15;
             break;                               
	case massive_star_mass_limit:            return 25;
             break;                               
	case upper_ZAMS_mass_for_degenerate_core: return 2.3;
             break;                                
	case minimum_main_sequence:              return 0.075;
             break;                              // for pop III: 0.095 [Msun] 
	case helium_dwarf_mass_limit:            return 0.45;
             break;                                
	case carbon_dwarf_mass_limit:            return 1.2;
             break;
	case Chandrasekar_mass:                  return 1.44;
             break;                                
	case helium2neutron_star:                return 2.2;
             break;                                
        case COcore2black_hole:                  return 5; // was 10
                                                      // (SPZ+GN: 27 Jul 2000)
             break;                                
	case super_giant2neutron_star:           return 8;
             break;                                
        case super_giant2black_hole:             return 25; // was 40
             break;                                
	case maximum_main_sequence:              return 100;
             break;                               
	case minimum_helium_star:                return 0.33;
             break;
	case helium_star_lifetime_fraction:      return 0.9;
             break;
	case helium_star_final_core_fraction:    return 0.80;
             break;
        default:
	     cerr << "\nNo recognized option in "
		     "stellar_evolution_constants::"
		     "parameters(stellar_mass_limits "
	          << pm << ")"
		  << endl;
             exit(1);
    }
}

int stellar_evolution_constants::use_common_envelope_method() {
  int cc_parameter = 1; // default (use alpha-gamma)
  //  cc_parameter = 2;     // default (use gamma-gamma)
  //  cc_parameter = 3;     // default (use alpha-alpha)
  return cc_parameter;
}

bool stellar_evolution_constants::parameters(boolean_parameter pb) {

    switch(pb) {                                
	case hyper_critical:                        return false;
             break;                                // 10^8 x Eddington
                                                   // see: Chevalier 1993
        case super_giant_disintegration:               return false;
             break;
        case proto_star_to_binary:                     return false;
             break;						    
        case impulse_kick_for_black_holes:             return true;
             break;				
        case use_angular_momentum_tidal:               return false;
	     break;
        case use_common_envelope_gamma_gamma:
             if (cnsts.use_common_envelope_method()==2) 
	       return true;
	     else
	       return false;
	     break;
        case use_common_envelope_alpha_alpha:
             if(cnsts.use_common_envelope_method()==3) 
	       return true;
	     else
	       return false;
	     break;
        default:
	     cerr << "\nNo recognized option in "
		     "stellar_evolution_constants::"
		     "parameters(boolean_parameter "
	          << pb << ")"
		  << endl;
             exit(1);
    }
}
    
real stellar_evolution_constants::parameters(accretion_parameter pa) {

    switch(pa) {                        // Mass of accretion to core.
	case black_hole_accretion_limit:        return 0.1;
             break;                             // 
	case neutron_star_accretion_limit:      return 0.05;
             break;                             //
	case white_dwarf_accretion_limit:       return 0.01;
             break;                             //
        case thermo_nuclear_flash:              return 1;  // used to be 0.05;
             break;                             // mass fractioncauses flash
        default:
	     cerr << "\nNo recognized option in "
		     "stellar_evolution_constants::"
		     "parameters(accretion_parameter "
	          << pa << ")"
		  << endl;
             exit(1);

      }
}


real stellar_evolution_constants::parameters(model_parameter pm) {

    switch(pm) {
	case star_formation_efficiency:                return 1.0;
             break;                             // [%] 
	case star_formation_timescale:                 return 1.0; // [Myr]
             break;                                 // if(<0) abs times KH timescale 
	case magnetic_mass_limit:                      return 0.7;
             break;                             // [msun] magnetic if(M<this)
	case magnetic_braking_exponent:                return 2.5;
	      break;                        
	case corotation_eccentricity:                  return 0.001;
	      break;                          // corotation if(ecc<this)
	case tidal_circularization_radius:             return 5.0;
	      break;                          
	case core_overshoot:                           return 0.125;
	      break;                          // overshoot fraction.
        case hydrogen_fraction:                        return 0.7;
              break;                          // X
 	case common_envelope_efficiency:               return 4;
	      break;                          // (alpha_ce)
        case envelope_binding_energy:                  return 0.5;
	      break;                          // (lambda)
	case specific_angular_momentum_loss:           return 3.;
	      break;                          // (beta)
        case dynamic_mass_transfer_gamma:               return 1.75;
              break;			      // (gamma)
        case non_massive_star_envelope_fraction_lost:  return 0.03;
	      break;                          
        case massive_star_envelope_fraction_lost:      return 0.9;
	      break;               
        case relaxation_driven_mass_loss_constant:     return 1;           
             break;						       
        case massive_star_mass_loss_law:               return 6.8;
	      break;                          
        case time_independent_mass_loss_law:           return 1;
	      break;
        case Darwin_Riemann_instability_factor:        return
	                                       cnsts.mathematics(one_third);
	      break;
        case homogeneous_sphere_gyration_radius_sq:    return 0.4;
	      break;
        case radiative_star_gyration_radius_sq:        return 0.03;
	      break;
        case convective_star_gyration_radius_sq:       return 0.2;
	      break;
        case rejuvenation_exponent:                    return 1;
	      break;
        case spiral_in_time:                           return 0.0005; // Myr
	      break;
        default:
	     cerr << "\nNo recognized option in "
		     "stellar_evolution_constants::"
		     "parameters(model_parameter "
	          << pm << ")"
		  << endl;
             exit(1);
    }
}


real stellar_evolution_constants::parameters(observational_parameter pm) {

    switch(pm) {
        case B_emission_star_mass_limit:               return 0.1;
             break; 
        case Barium_star_mass_limit:                   return 0.01;
	      break;                        
        case Blue_straggler_mass_limit:                return 0;
	      break;                        
        default:
	     cerr << "\nNo recognized option in "
		     "stellar_evolution_constants::"
		     "parameters(observational_parameter "
	          << pm << ")"
		  << endl;
             exit(1);
    }
}

// Safetly parameters means that changing these may cause
// unexpected errors and changes the bahavious of the program.
real stellar_evolution_constants::safety(safety_parameter ps) {

  switch(ps) {
    case timestep_factor:                      return 0.01;       
          break;                             // 0.01 ok, but big for AM CVns
    case maximum_binary_update_time_fraction: return 1.0;// 0.9;
          break;                            // see also star_to_dyn
    case minimum_timestep:                   return 1.e-11; // Do not change!
          break;                                   // == 10*HubbleT*precision
    case minimum_mass_step:                  return 1.e-5;
          break;                                   // Tricky used to be 1.e-7
    case maximum_timestep:                   return 1000; // was 1000  
          break;                                   // 2.5 works fine but slow
    case maximum_recursive_calls:            return 1000;
          break;
      case tiny:                             return 5e-13;
          break;
    case number_of_steps:                   return 10;
        break;
        default:
	     cerr << "\nNo recognized option in "
		     "stellar_evolution_constants::"
		     "parameters(safety_parameter "
	          << ps << ")"
		  << endl;
             exit(1);
  } 
}

real stellar_evolution_constants::star_to_dyn(dynamics_update_parameter dup) {

    switch(dup) {
        case stellar_mass_update_limit:                return 0.001;
             break; 
        case semi_major_axis_update_limit:             return 0.001;
	      break;                        
        case binary_update_time_fraction:              return 0.001; //was 0.05
             break;
        default:
	     cerr << "\nNo recognized option in "
		     "stellar_evolution_constants::"
		     "star_to_dyn(dynamica_update_parameter "
	          << dup << ")"
		  << endl;
             exit(1);
    }
}

local real ap(real zeta, real a, real b=0, real c=0, real d=0, real e=0) {

  real ai = a + zeta*(b + zeta*(c + zeta*(d + zeta*e)));

  return ai;
}

real stellar_model_constants::a(int index, real z) {

  real zeta = log10(z/cnsts.parameters(solar_metalicity));

  real a = 0;
  switch(index) {
     case 1: a = ap(zeta, 1593.890, 2053.038, 1231.226, 232.7785);
             break;
     case 2: a = ap(zeta, 2.706708E+3, 1.483131E+3, 5.772723E+2, 7.411230E+1);  
             break;
     case 3: a = ap(zeta, 1.466143E+2, -1.048442E+2, -6.795374E+1, -1.391127E+1);
             break;
     case 4: a = ap(zeta, 4.141960E-2, 4.564888E-2, 2.958542E-2, 5.571483E-3);
             break;
     case 5: a = ap(zeta, 3.426349E-1);
             break;
     case 6: a = ap(zeta, 1.949814E+1, 1.758178E+0, -6.008212E+0, -4.470533E+0); 
             break;
     case 7: a = ap(zeta, 4.903830E+0);
             break;
     case 8: a = ap(zeta, 5.212154E-2, 3.166411E-2, -2.750074E-3, -2.271549E-3); 
             break;
     case 9: a = ap(zeta, 1.312179E+0, -3.294936E-1, 9.231860E-2, 2.610989E-2);
             break;
     case 10: a = ap(zeta, 8.073972E-1);
              
             break;
     case 11: {real a11 = ap(zeta, 1.031538E+0, -2.434480E-1, 7.732821E+0, 
			     6.460705E+0, 1.374484E+0);
               a = a11 * smc.a(14, z);
     }
             break;
     case 12: {real a12 = ap(zeta, 1.043715E+0, -1.577474E+0, -5.168234E+0, 
			     -5.596506E+0, -1.299394E+0);
               a = a12 * smc.a(14, z);
     }
             break;
     case 13: a = ap(zeta, 7.859573E+2, -8.542048E+0, -2.642511E+1, 
		     -9.585707E+0);
             break;
     case 14: a = ap(zeta, 3.858911E+3, 2.459681E+3, -7.630093E+1, 
		     -3.486057E+2, -4.861703E+1);
             break;
     case 15: a = ap(zeta, 2.888720E+2, 2.952979E+2, 1.850341E+2, 3.797254E+1); 
             break;
     case 16: a = ap(zeta, 7.196580E+0, 5.613746E-1, 3.805871E-1, 8.398728E-2); 
             break;
     case 17: {real s = log10(z);
               real log_a17 = max(0.097 - 0.1072*(s+3),
				 max(0.097, 
				     min(0.1461, 0.1461 + 0.1237*(s+2))));
              a = pow(10., log_a17);
     }
             break; 
     case 18: {real a18 = ap(zeta, 2.187715E-1, -2.154437E+0, -3.768678E+0, 
		            -1.975518E+0, -3.021475E-1);
              a = a18 * smc.a(20, z);
     }
             break;
     case 19: {real a19 = ap(zeta, 1.466440E+0, 1.839725E+0, 6.442199E+0, 
			        4.023635E+0, 6.957529E-1);
              a = a19 * smc.a(20, z);
     }
             break;
     case 20: a = ap(zeta, 2.652091E+1, 8.178458E+1, 1.156058E+2, 7.633811E+1,
		     1.950698E+1);
             break;
     case 21: a = ap(zeta, 1.472103E+0, -2.947609E+0, -3.312828E+0, 
		     -9.945065E-1);
             break;
     case 22: a = ap(zeta, 3.071048E+0, -5.679941E+0, -9.745523E+0, 
		     -3.594543E+0);
             break;
     case 23: a = ap(zeta, 2.617890E+0, 1.019135E+0, -3.292551E-2, 
		     -7.445123E-2);
             break;
     case 24: a = ap(zeta, 1.075567E-2, 1.773287E-2, 9.610479E-3, 1.732469E-3);
             break;
     case 25: a = ap(zeta, 1.476246E+0, 1.899331E+0, 1.195010E+0, 3.035051E-1);
             break;
     case 26: a = ap(zeta, 5.502535E+0, -6.601663E-2, 9.968707E-2, 3.599801E-2);
             break;
     case 27: a = ap(zeta, 9.511033E+1, 6.819618E+1, -1.045625E+1, -1.474939E+1);
             break;
     case 28: a = ap(zeta, 3.113458E+1, 1.012033E+1, -4.650511E+0, -2.463185E+0);
             break;
     case 29: {real a29 = ap(zeta, 1.413057E+0, 4.578814E-1, -6.850581E-2, 
		                  -5.588658E-2);
              a = pow(a29, smc.a(32, z)); 
     }
             break;
     case 30: a = ap(zeta, 3.910862E+1, 5.196646E+1, 2.264970E+1, 2.873680E+0);
             break;
     case 31: a = ap(zeta, 4.597479E+0, -2.855179E-1, 2.709724E-1);
             break;
     case 32: a = ap(zeta, 6.682518E+0, 2.827718E-1, -7.294429E-2);
             break;
     case 33: a = max(0.6355 - 0.4192*zeta, 
		      max(1.25, 
			  min(1.4, 1.5135 + 0.3769*zeta)));
             break;
     case 34: a = ap(zeta, 1.910302E-1, 1.158624E-1, 3.348990E-2, 2.599706E-3);
             break;
     case 35: a = ap(zeta, 3.931056E-1, 7.277637E-2, -1.366593E-1, -4.508946E-2);
             break;
     case 36: a = ap(zeta, 3.267776E-1, 1.204424E-1, 9.988332E-2, 2.455361E-2);
             break;
     case 37: a = ap(zeta, 5.990212E-1, 5.570264E-2, 6.207626E-2, 1.777283E-2);
             break;

     case 38: a = ap(zeta, 7.330122E-1, 5.192827E-1, 2.316416E-1, 8.346941E-3);
             break;
     case 39: a = ap(zeta, 1.172768E+0, -1.209262E-1, -1.193023E-1, -2.859837E-2);
             break;
     case 40: a = ap(zeta, 3.982622E-1, -2.296279E-1, -2.262539E-1, -5.219837E-2);
             break;
     case 41: a = ap(zeta, 3.571038E+0, -2.223625E-2, -2.611794E-2, -6.359648E-3);
             break;
     case 42: {real a42 = ap(zeta, 1.9848E+0, 1.1386E+0, 3.5640E-1);
              a = min(1.25, max(1.1, a42));
     }
             break;
     case 43: a = ap(zeta, 6.300E-2, 4.810E-2, 9.840E-3);
             break;
     case 44: {real a44 = ap(zeta, 1.200E+0, 2.450E+0);
              a = min(1.3, max(0.45, a44));
     }
             break;
     case 45: a = ap(zeta, 2.321400E-1, 1.828075E-3, -2.232007E-2, -3.378734E-3);
             break;
     case 46: a = ap(zeta, 1.163659E-2, 3.427682E-3, 1.421393E-3, -3.710666E-3);
             break;
     case 47: a = ap(zeta, 1.048020E-2, -1.231921E-2, -1.686860E-2, -4.234354E-3);
             break;
     case 48: a = ap(zeta, 1.555590E+0, -3.223927E-1, -5.197429E-1, -1.066441E-1);
             break;
     case 49: {real a49 = ap(zeta, 9.7700E-2, -2.3100E-1, -7.5300E-2);
              a = max(a49, 0.145);
     }
             break;
     case 50: {real a50 = ap(zeta, 2.4000E-1, 1.8000E-1, 5.9500E-1);
              a = min(a50, 0.306 + 0.053*zeta);
     }
             break;
     case 51: {real a51 = ap(zeta, 3.3000E-1, 1.3200E-1, 2.1800E-1);
              a = min(a51, 0.3625 + 0.062*zeta);
     }
             break;
     case 52: {real a52 = ap(zeta, 1.1064E+0, 4.1500E-1, 1.8000E-1);
       			a = max(a52, 0.9);
                if(z>0.01)
					a = min(a, 1.0);
	 }
             break;
     case 53: {real a53 = ap(zeta, 1.1900E+0, 3.7700E-1, 1.7600E-1);
				a = max(a53, 1.0);
	            if(z>0.01)
					a = min(a, 1.1);
     }
             break;
     case 54: a = ap(zeta, 3.855707E-1, -6.104166E-1, 5.676742E+0, 1.060894E+1,
		    5.284014E+0);
             break;
     case 55: a = ap(zeta, 3.579064E-1, -6.442936E-1, 5.494644E+0, 1.054952E+1,
		    5.280991E+0);
             break;
     case 56: a = ap(zeta, 9.587587E-1, 8.777464E-1, 2.017321E-1);
             break;
     case 57: a = max(0.6355 - 0.4192*zeta, max(1.25, 
			 min(1.4, 1.5135 + 0.3769*zeta)));
             break;
     case 58: a = ap(zeta, 4.907546E-1, -1.683928E-1, -3.108742E-1, -7.202918E-2);
             break;
     case 59: a = ap(zeta, 4.537070E+0, -4.465455E+0, -1.612690E+0, -1.623246E+0);
             break;
     case 60: a = ap(zeta, 1.796220E+0, 2.814020E-1, 1.423325E+0, 3.421036E-1);
             break;
     case 61: a = ap(zeta, 2.256216E+0, 3.773400E-1, 1.537867E+0, 4.396373E-1);
             break;
     case 62: {real a62 = ap(zeta, 8.4300E-2, -4.7500E-2, -3.5200E-2);
              a = max(0.065, a62);
     }
             break;
     case 63: {a = ap(zeta, 7.3600E-2, 7.4900E-2, 4.4260E-2);
              if(z<0.004)
                  a = min(0.055, a);
    }
             break;
     case 64: {real a64 = ap(zeta, 1.3600E-1, 3.5200E-2);
              a = max(0.091, min(0.121, a64));
              if (smc.a(68, z) >= smc.a(66, z)) {
                  a =  smc.a(58, z)*pow(smc.a(66, z), smc.a(60, z))
                  / (smc.a(59, z) + pow(smc.a(66, z), smc.a(61, z)));
                                   
              }
     }
             break;
     case 65: {a = ap(zeta, 1.564231E-3, 1.653042E-3, -4.439786E-3, 
		     -4.951011E-3, -1.216530E-03);
     }
             break;
     case 66: {real a66 = ap(zeta, 1.4770E+0, 2.9600E-1);
             a = max(0.8, 
		     min(0.8 - 2.0*zeta, 
			 max(a66, 
			     min(1.6, -0.308 - 1.046*zeta))));
     }
             break;
     case 67: a = ap(zeta, 5.210157E+0, -4.143695E+0, -2.120870E+0);
             break;
     case 68: {real a68 = ap(zeta, 1.1160E+0, 1.6600E-1);
             a68 = max(0.9, min(a68, 1.0));
             a = min(a68, smc.a(66, z));
     }
             break;

     case 69: a = ap(zeta, 1.071489E+0, -1.164852E-1, -8.623831E-2, 
		     -1.582349E-2);
             break;
     case 70: a = ap(zeta, 7.108492E-1, 7.935927E-1, 3.926983E-1, 3.622146E-2);
             break;
     case 71: a = ap(zeta, 3.478514E+0, -2.585474E-2, -1.512955E-2, 
		     -2.833691E-3);
             break;
     case 72: {real a72 = ap(zeta, 9.132108E-1, -1.653695E-1, 0.0,  3.636784E-2);
             if(z>0.01)
	       a = max(a72, 0.95);
	     else
	       a = a72;
     }
             break;
     case 73: a = ap(zeta, 3.969331E-3, 4.539076E-3, 1.720906E-3, 
		     1.897857E-4);
             break;
     case 74: {real a74 = ap(zeta, 1.600E+0, 7.640E-1, 3.322E-1);
             a = max(1.4, min(a74, 1.6));
     }
             break;
     case 75: {real a75 = ap(zeta, 8.109E-1, -6.282E-1);
             a = max(max(1.0, 
			 min(a75, 1.27)), 0.6355 - 0.4192*zeta);
     }
             break;
     case 76: {real a76 = ap(zeta, 1.192334E-2, 1.083057E-2, 
			     1.230969E+0, 1.551656E+0);
             a = max(a76, 
		     -0.1015564 - 0.2161264*zeta - 0.05182516*pow(zeta, 2));
     }
             break;
     case 77: {real a77 = ap(zeta, -1.668868E-1, 5.818123E-1, -1.105027E+1, -1.668070E+1);
             a = max(-0.3868776 - 0.5457078*zeta - 0.1463472*pow(zeta, 2),
		     min(0.0, a77));
     }
             break;
     case 78: {real a78 = ap(zeta, 7.615495E-1, 1.068243E-1, 
			          -2.011333E-1, -9.371415E-2);
             a = max(0.0, min(a78, 7.454 + 9.046*zeta));
     }
             break;
     case 79: {real a79 = ap(zeta, 9.409838E+0, 1.522928E+0);
             a = min(a79, max(2.0, - 13.3 - 18.6*zeta));
     }
             break;
     case 80: {real a80 = ap(zeta, -2.7110E-1, -5.7560E-1, -8.3800E-2);
             a = max(0.0585542, a80);
     }
             break;
     case 81: {real a81 = ap(zeta, 2.4930E+0, 1.1475E+0);
             a = min(1.5, max(0.4, a81));
     }
             break;

     default:  cerr << "Not a valid index for "
		    << "constant.stellar_model_constant::a(int i= " 
		    << index << ")" << endl;
  }

  return a;
}


local real bp(real zeta, real a, real b=0, real c=0, real d=0, real e=0) {

  // According to HPT2000 ap == bp (see Appendix A)
  return ap(zeta, a, b, c, d, e);
}

real stellar_model_constants::b(int index, real z) {

  real zeta = log10(z/cnsts.parameters(solar_metalicity));

  real b = 0;
  switch(index) {
     case 1: {real b1 = bp(zeta, 3.9700E-1, 2.8826E-1, 5.2930E-1);
              b = min(b1, 0.54);}
             break;
     case 2: {real s = log10(z);
             real b2 = pow(10, -4.6739 - 0.9394*s);
	     b = min(max(b2, -0.04167+55.67*z), 0.4771 -9329.21*pow(z, 2.94));
     }
	     break;
     case 3: {real s = log10(z);
             real b3 = max(-0.1451, -2.2794 - s*(1.5175 + s*0.254));
             b = pow(10, b3);
	     if(z>0.004) 
	       b = max(b, 0.7307 + 14265.1*pow(z, 3.395));
             }
	     break;
     case 4: {real b4 = bp(zeta, 9.960283E-1, 8.164393E-1, 
			  2.383830E+0, 2.223436E+0, 8.638115E-1);
             b = b4 + 0.1231572*pow(zeta, 5);
     }
             break;
     case 5: b = bp(zeta, 2.561062E-1, 7.072646E-2, -5.444596E-2, -5.798167E-2,
		    -1.349129E-2);
             break;
     case 6: {real b6 = bp(zeta, 1.157338E+0, 1.467883E+0, 4.299661E+0, 
			  3.130500E+0, 6.992080E-1);
             b = b6 + 0.01640687*pow(zeta, 5);
     }
             break;
     case 7: b = bp(zeta, 4.022765E-1, 3.050010E-1, 9.962137E-1, 7.914079E-1,
		    1.728098E-1);
             break;

     case 8: cerr << "Unknown index for b8"<<endl;
          break;
          
     case 9: b = bp(zeta, 2.751631E+3, 3.557098E+2);
             break;
     case 10: b = bp(zeta, -3.820831E-2, 5.872664E-2);
              break;
     case 11: {real b11 = bp(zeta, 1.071738E+2, -8.970339E+1, -3.949739E+1);
              b = pow(b11, 2);
     }
             break;
     case 12: b = bp(zeta, 7.348793E+2, -1.531020E+2, -3.793700E+1);
             break;
     case 13: {real b13 = bp(zeta, 9.219293E+0, -2.005865E+0, -5.561309E-1);
              b = pow(b13, 2);
     }
             break;
     case 14: {real b14 = bp(zeta, 2.917412E+0, 1.575290E+0, 5.751814E-1);
              b = pow(b14, smc.b(15,z));
     }
             break;
     case 15: b = bp(zeta, 3.629118E+0, -9.112722E-1, 1.042291E+0);
             break;
     case 16: {real b16 = bp(zeta, 4.916389E+0, 2.862149E+0, 7.844850E-1);
              b = pow(b16, smc.b(15,z));
     }
	      break;
     case 17: b = 1;
              if (zeta>-1) 
		b = 1 - 0.3880523*pow(zeta+1, 2.862149);
	     break;
     case 18: b = bp(zeta, 5.496045E+1, -1.289968E+1, 6.385758E+0);
             break;
     case 19: b = bp(zeta, 1.832694E+0, -5.766608E-2, 5.696128E-2);
             break;
     case 20: b = bp(zeta, 1.211104E+2);
             break;
     case 21: b = bp(zeta, 2.214088E+2, 2.187113E+2, 1.170177E+1, -2.635340E+1);
             break;
     case 22: b = bp(zeta, 2.063983E+0, 7.363827E-1, 2.654323E-1, -6.140719E-2);
             break;
     case 23: b = bp(zeta, 2.003160E+0, 9.388871E-1, 9.656450E-1, 2.362266E-1);
             break;
     case 24: {real b24 = bp(zeta, 1.609901E+1, 7.391573E+0, 
			    2.277010E+1, 8.334227E+0);
              b = pow(b24, smc.b(28, z));
     }
             break;
     case 25: b = bp(zeta, 1.747500E-1, 6.271202E-2, -2.324229E-2, -1.844559E-2);
             break;
     case 26: b = 5 - 0.09138012*pow(z, -0.3671407);
             break;
     case 27: {real b27 = bp(zeta, 2.752869E+0, 2.729201E-2, 
			    4.996927E-1, 2.496551E-1);
              b = pow(b27, 2*smc.b(28, z));
     }
             break;
     case 28: b = bp(zeta, 3.518506E+0, 1.112440E+0, -4.556216E-1, -2.179426E-1);
             break;


     case 29: b = bp(zeta, 1.626062E+2, -1.168838E+1, -5.498343E+0);
             break;
     case 30: b = bp(zeta, 3.336833E-1, -1.458043E-1, -2.011751E-2);
             break;
     case 31: {real b31 = bp(zeta, 7.425137E+1, 1.790236E+1, 
			    3.033910E+1, 1.018259E+1);
              b = pow(b31, smc.b(33, z));
     }
             break;
     case 32: b = bp(zeta, 9.268325E+2, -9.739859E+1, -7.702152E+1, 
		     -3.158268E+1);
             break;
     case 33: b = bp(zeta, 2.474401E+0, 3.892972E-1);
             break;
     case 34: {real b34 = bp(zeta, 1.127018E+1, 1.622158E+0, 
			    -1.443664E+0, -9.474699E-1);
              b = pow(b34, smc.b(33, z));
     }
             break;
     case 35: cerr << "Unknown index for b35"<<endl;
	     break;
     case 36: {real b36 = bp(zeta, 1.445216E-1, -6.180219E-2, 
			    3.093878E-2, 1.567090E-2);
              b = pow(b36, 4);
     }
             break;
     case 37: {real b37 = bp(zeta, 1.304129E+0, 1.395919E-1, 
			    4.142455E-3, -9.732503E-3);
              b = 4*b37;
     }
             break;
     case 38: {real b38 = bp(zeta, 5.114149E-1, -1.160850E-2);
              b = pow(b38, 4);
     }
             break;
     case 39: b = bp(zeta, 1.314955E+2, 2.009258E+1, -5.143082E-1, -1.379140E+0);
             break;
     case 40: {real b40 = bp(zeta, 1.823973E+1, -3.074559E+0, -4.307878E+0);
             b = max(b40, 1.);
     }
             break;
     case 41: {real b41 = bp(zeta, 2.327037E+0, 2.403445E+0, 
			    1.208407E+0, 2.087263E-1);
             b = pow(b41, smc.b(42, z));
     }
             break;
     case 42: b = bp(zeta, 1.997378E+0, -8.126205E-1);
             break;
     case 43: b = bp(zeta, 1.079113E-1, 1.762409E-2, 1.096601E-2, 3.058818E-3);
             break;
     case 44: {real b44 = bp(zeta, 2.327409E+0, 6.901582E-1, 
			    -2.158431E-1, -1.084117E-1);
              b = pow(b44, 5);
     }
             break;
     case 45: {real p = zeta+1;
              if(p<=0)
		b = 1;
	      else
		b = 1 - (2.47162*p - 5.401682*pow(p, 2) + 3.247361*pow(p, 3));
     }
	      break;
     case 46: {b = bp(zeta, 2.214315E+0, -1.975747E+0);
     }
             break;
     case 47: {real p = zeta+1;
              b = 1.127733*p + 0.2344416*pow(p, 2) - 0.3793726*pow(p, 3);
     }
	      break;
     case 48: b = bp(zeta, 5.072525E+0, 1.146189E+1, 6.961724E+0, 1.316965E+0);
             break;
     case 49: b = bp(zeta, 5.139740E+0);
             break;
     case 51: {real b51 = bp(zeta, 1.125124E+0, 1.306486E+0, 
			    3.622359E+0, 2.601976E+0, 3.031270E-1);
             b = b51 - 0.1343798*pow(zeta, 5);
     }
             break;
     case 52: b = bp(zeta, 3.349489E-1, 4.531269E-3, 1.131793E-1, 2.300156E-1,
		     7.632745E-2);
             break;
     case 53: {real b53 = bp(zeta, 1.467794E+0, 2.798142E+0, 9.455580E+0, 
			    8.963904E+0, 3.339719E+0);
             b = b53 + 0.4426929*pow(zeta, 5);
     }
             break;
     case 54: b = bp(zeta, 4.658512E-1, 2.597451E-1, 9.048179E-1, 7.394505E-1,
		     1.607092E-1);
             break;
     case 55: {real b55 = bp(zeta, 1.0422E+0, 1.3156E-1, 4.5000E-2);
             b = min(b55, 0.99164 - 743.123*pow(z, 2.83));
     }
             break;
     case 56: {real b56 = bp(zeta, 1.110866E+0, 9.623856E-1, 2.735487E+0, 
			    2.445602E+0, 8.826352E-1);
              b = b56 + 0.1140142*pow(zeta, 5);
     }
             break;
     case 57: {real b57 = bp(zeta, -1.584333E-1, -1.728865E-1, -4.461431E-1, 
		     -3.925259E-1, -1.276203E-1);
              b = b57 - 0.01308728*pow(zeta, 5);
     }
             break;
     default:  cerr << "Not a valid index for "
		    << "constant.stellar_model_constant::b(int i= " 
		    << index << ")" << endl;
     }

return b;
}




local real cp(real zeta, real a, real b=0, real c=0, real d=0, real e=0) {
  // According to Tout96 cp == ap (== bp) of HPT 2000 
  return ap(zeta, a, b, c, d, e);
}

real stellar_model_constants::c(int index, real z) {

  real zeta = log10(z/cnsts.parameters(solar_metalicity));

  real c = 0;
  switch(index) {
     case 1: c = cp(zeta, 3.970417E-01, -3.2913574E-01, 3.4776688E-01, 
			3.7470851E-01, 9.011915E-02);
             break;
     case 2: c = cp(zeta, 8.527626E+00,-2.441225973E+01, 5.643597107E+01, 
			3.706152575E+01, 5.4562406E+00);  
             break;
     case 3: c = cp(zeta,2.5546E-04, -1.23461E-03, -2.3246E-04, 
			4.5519E-04, 1.6176E-04 );
             break;
     case 4: c = cp(zeta, 5.432889E+00, -8.62157806E+00, 1.344202049E+01, 
			1.451584135E+01, 3.39793084E+00);
             break;
     case 5: c = cp(zeta, 5.563579E+00,-1.032345224E+01, 1.944322980E+01, 
			1.897361347E+01, 4.16903097E+00);
             break;
     case 6: c = cp(zeta, 7.8866060E-01, -2.90870942E+00,  6.54713531E+00,
			4.05606657E+00, 5.3287322E-01); 
             break;
     case 7: c = cp(zeta, 5.86685E-03, -1.704237E-02, 3.872348E-02, 
			2.570041E-02, 3.83376E-03);
             break;
     case 8: c = cp(zeta, 1.715359E+00, 6.2246212E-01, -9.2557761E-01, 
			-1.16996966E+00, -3.0631491E-01); 
             break;
     case 9: c = cp(zeta, 6.597788E+00, -4.2450044E-01,-1.213339427E+01,
			-1.073509484E+01, -2.51487077E+00);
             break;
     case 10: c = cp(zeta,1.008855000E+01, -7.11727086E+00,-3.167119479E+01, 
			-2.424848322E+01,-5.33608972E+00 );
             break;
     case 11: c = cp(zeta, 1.012495E+00, 3.2699690E-01, -9.23418E-03, 
			-3.876858E-02, -4.12750E-03);
             break;
     case 12: c = cp(zeta,7.490166E-02, 2.410413E-02, 7.233664E-02, 
			3.040467E-02, 1.97741E-03 );
 	     break;
     case 13: c = cp(zeta,1.077422E-02 );
             break;
     case 14: c = cp(zeta, 3.082234E+00, 9.447205E-01, -2.15200882E+00, 
			-2.49219496E+00, -6.3848738E-01);
             break;
     case 15: c = cp(zeta,1.784778E+01, -7.4534569E+00,-4.896066856E+01,
			-4.005386135E+01, -9.09331816E+00 ); 
             break;
     case 16: c = cp(zeta,2.2582E-04, -1.86899E-03, 3.88783E-03, 
			1.42402E-03,-7.671E-05 ); 
             break;
  }
  return c;
}     


