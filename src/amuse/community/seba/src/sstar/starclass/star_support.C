//
// star_support: Helper functions to aid in setting up and manipulating
//

#include "star_support.h"

char* type_string(stellar_type tpe) {
   
      switch(tpe) {
         case Static_Star:		return "static_star";
         case SPZDCH_Star:		return "SPZDCH_star";
         case NAS:			return "not_a_star";
         case Proto_Star:		return "proto_star";
         case Planet:		        return "planet";
         case Brown_Dwarf:		return "brown_dwarf";
         case Main_Sequence: 		return "main_sequence";
         case Hyper_Giant:		return "hyper_giant";
         case Hertzsprung_Gap:		return "hertzsprung_gap";
         case Sub_Giant:		return "sub_giant";
         case Horizontal_Branch:	return "horizontal_branch";
         case Super_Giant:		return "super_giant";
         case Carbon_Star:		return "carbon_star";
         case Helium_Star:		return "helium_star";
         case Helium_Giant:		return "helium_giant";
         case Helium_Dwarf:		return "helium_dwarf";
         case Carbon_Dwarf:		return "carbon_dwarf";
         case Oxygen_Dwarf:		return "oxygen_dwarf";
         case Thorn_Zytkow:		return "thorne_zytkow";
         case Xray_Pulsar:		return "xray_pulsar";
         case Radio_Pulsar:		return "radio_pulsar";
         case Neutron_Star:		return "neutron_star";
         case Black_Hole:		return "black_hole";
         case Disintegrated:		return "Disintegrated";
         case Double:			return "binary"; 
         default:			return "unknown_stellar_type";
      }
   }

char* type_short_string(stellar_type tpe) {

      switch(tpe) {
         case Static_Star:              return "SS";
         case SPZDCH_Star:              return "PZHs";
         case NAS:                      return "nas";
         case Proto_Star:		return "ps";
         case Planet:                   return "pl";
         case Brown_Dwarf:              return "bd";
         case Main_Sequence:            return "ms";
         case Hyper_Giant:              return "hy";
         case Hertzsprung_Gap:          return "hg";
         case Sub_Giant:                return "gs";
         case Horizontal_Branch:        return "hb";
         case Super_Giant:              return "sg";
         case Carbon_Star:              return "co";
         case Helium_Star:              return "he";
         case Helium_Giant:		return "gh";
         case Helium_Dwarf:             return "hd";
         case Carbon_Dwarf:             return "cd";
         case Oxygen_Dwarf:             return "od";
         case Thorn_Zytkow:             return "TZO";
         case Xray_Pulsar:              return "xp";
         case Radio_Pulsar:             return "rp";
         case Neutron_Star:             return "ns";
         case Black_Hole:               return "bh";
         case Disintegrated:		return "di";
         case Double:			return "bin"; 
         default:                       return "?";
      }
   }

char* type_string(stellar_type_summary tpe) {
   
      switch(tpe) {
         case ZAMS:			return "ZAMS";
         case Early_Giant: 		return "early_giant";
         case Late_Giant:		return "late_giant";
         case Helium_Remnant:		return "helium_remnant";
         case White_Dwarf:		return "white_dwarf";
         case Neutron_Remnant:		return "neutron_remnant";
         case Inert_Remnant:		return "inert_remnant";
         case Unspecified:		return "unspecified";
         case Undefined:		return "undefined";
         default:			return "unknown_stellar_type_summary";
      }
   }

char* type_short_string(stellar_type_summary tpe) {

      switch(tpe) {
	case ZAMS:                     return "ms";
       case Early_Giant:              return "ge";
       case Late_Giant:               return "gl";
       case Helium_Remnant:           return "dh";
       case White_Dwarf:              return "dw";
       case Neutron_Remnant:          return "rn";
       case Inert_Remnant:            return "ri";
       case Unspecified:              return "unspecified";
       case Undefined:                return "undefined";
       default:                       return "unknown_stellar_type_summary";
     }
}

stellar_type_summary summarize_stellar_type(stellar_type tpe) {
   
      switch(tpe) {
         case Static_Star:	        
         case SPZDCH_Star:	        
         case Proto_Star:	        
         case Planet:		        
         case Brown_Dwarf:		return Unspecified;
	   
         case Main_Sequence: 		return ZAMS;
	 case Hyper_Giant:
         case Hertzsprung_Gap:		
         case Sub_Giant:		return Early_Giant;

         case Horizontal_Branch:
         case Super_Giant:
         case Thorn_Zytkow:		return Late_Giant;

         case Carbon_Star:	
         case Helium_Star:	
         case Helium_Giant:		return Helium_Remnant;
	   
         case Helium_Dwarf:
         case Carbon_Dwarf:
         case Oxygen_Dwarf:		return White_Dwarf;

         case Xray_Pulsar:
         case Radio_Pulsar:
         case Neutron_Star:             return Neutron_Remnant;


         case Black_Hole:		return Inert_Remnant;


         case NAS:		
         case Disintegrated:	
         case Double:		
         default:			return Undefined;
      }
   }

char* type_string(spectral_class class_tpe) {

      switch(class_tpe) {
         case O5:			return "O5";
         case O6:			return "O6";
         case O7:			return "O7";
         case O8:			return "O8";
         case O9:			return "O9";
         case O95:			return "O9.5";
         case B0:			return "B0";
         case B05:			return "B0.5";
         case B1:			return "B1";
         case B2:			return "B2";
         case B3:			return "B3";
         case B5:			return "B5";
         case B6:			return "B6";
         case B7:			return "B7";
         case B8:			return "B8";
         case B9:			return "B9";
         case B95:			return "B9.5";
         case A0:			return "A0";
         case A1:			return "A1";
         case A2:			return "A2";
         case A3:			return "A3";
         case A4:			return "A4";
         case A5:			return "A5";
         case A7:			return "A7";
         case F0:			return "F0";
         case F2:			return "F2";
         case F3:			return "F3";
         case F5:			return "F5";
         case F6:			return "F6";
         case F7:			return "F7";
         case F8:			return "F8";
         case G0:			return "G0";
         case G1:			return "G1";
         case G2:			return "G2";
         case G5:			return "G5";
         case K0:			return "K0";
         case K5:			return "K5";
         case M0:			return "M0";
         case M5:			return "M5";
         case M8:			return "M8";
         case he:			return "he";
         case wd:			return "wd";
         case ns:			return "ns";
         case bh:			return "bh";
         case bd:			return "bd";
         case di:			return "di";
         case bin:			return "bin";
         default:			return " ";
      }
   }


char* type_short_string(spectral_class class_tpe) {

      switch(class_tpe) {
         case O5:
         case O6:
         case O7:
         case O8:
         case O9:
         case O95:			return "O";
         case B0:
         case B05:
         case B1:
         case B2:
         case B3:
         case B5:
         case B6:
         case B7:
         case B8:
         case B9:
         case B95:			return "B";
         case A0:
         case A1:
         case A2:
         case A3:
         case A4:
         case A5:
         case A7:			return "A";
         case F0:
         case F2:
         case F3:
         case F5:
         case F6:
         case F7:
         case F8:			return "F";
         case G0:
         case G1:
         case G2:
         case G5:			return "G";
         case K0:
         case K5:			return "K";
         case M0:
         case M5:
         case M8:			return "M";
         case he:			return "O";
         case wd:			return "O";
         case ns:
         case bh:
         case bd:
         case di:
         case bin:
         default:			return "N";
      }
   }

char* type_string(luminosity_class lum_c) {

    switch(lum_c) {
	case I:            return "I";
	case II:           return "II";
	case III:          return "III";
	case IV:           return "IV";
	case V:            return "V";
	default:           return "?";
    }
}

char* type_string(star_type_spec spec) {

      switch(spec) {
         case NAC:			return "";
         case Emission:			return "emission";
         case Barium:			return "Barium";
         case Blue_Straggler:		return "blue_straggler";
         case Rl_filling:               return "Roche_lobe_filling";
         case Runaway:			return "runaway";
         case Merger:			return "merger";
         case Accreting:		return "accreting";
         case Dsntgr:		        return "disintegrated";
         default:			return "?";
      }
   }

char* type_short_string(star_type_spec spec) {

      switch(spec) {
         case NAC:                      return "";
         case Emission:                 return "e";
         case Barium:                   return "Ba";
         case Blue_Straggler:           return "Bs";
         case Rl_filling:               return "Rlof";
         case Runaway:                  return "rnway";
         case Merger:                   return "mrgr";
         case Accreting:                return "accr";
         case Dsntgr:                   return "dsntgr";
         default:                       return "?";
      }
   }

//      On stellar temperature determined.
//spectral_class get_spectral_class(const real mass) {
spectral_class get_spectral_class(const real temperature) {
//	Morton, D.C., 1967, in Stellar Astronomy, Gordon and Breach
//	Schience Publishers, NY, (eds. H-Y. Chiu, R.L., Warasila and J.T., Remo),
//	volume 2, p. 157.
//	Values indicated with MB are from Mihallas, D. \& Binney, J., 1981 
//	Galactic Astronomy (second edition), (W.H. Freeman and Company NY),
//	Table 3.5, p. 111.


     real log_temp = log10(temperature);
     spectral_class spectral_type;
        if(log_temp>=5.574) spectral_type = O5;
        else if(log_temp>=4.562) spectral_type = O6;
        else if(log_temp>=4.553) spectral_type = O7;
        else if(log_temp>=4.544) spectral_type = O8;
        else if(log_temp>=4.535) spectral_type = O9;
        else if(log_temp>=4.506) spectral_type = O95;
        else if(log_temp>=4.490) spectral_type = B0;
        else if(log_temp>=4.418) spectral_type = B05;
        else if(log_temp>=4.354) spectral_type = B1;
        else if(log_temp>=4.312) spectral_type = B2;
        else if(log_temp>=4.253) spectral_type = B3;
        else if(log_temp>=4.193) spectral_type = B5;
        else if(log_temp>=4.164) spectral_type = B6;
        else if(log_temp>=4.134) spectral_type = B7;
        else if(log_temp>=4.079) spectral_type = B8;
        else if(log_temp>=4.029) spectral_type = B9;
        else if(log_temp>=4.000) spectral_type = B95;
        else if(log_temp>=4.000) spectral_type = A0;
        else if(log_temp>=3.982) spectral_type = A1;
        else if(log_temp>=3.969) spectral_type = A2;
        else if(log_temp>=3.958) spectral_type = A3;
        else if(log_temp>=3.936) spectral_type = A4;
        else if(log_temp>=3.929) spectral_type = A5;
        else if(log_temp>=3.914) spectral_type = A7;
        else if(log_temp>=3.876) spectral_type = F0;
        else if(log_temp>=3.860) spectral_type = F2;
        else if(log_temp>=3.845) spectral_type = F3;
        else if(log_temp>=3.833) spectral_type = F5;
        else if(log_temp>=3.818) spectral_type = F6;
        else if(log_temp>=3.804) spectral_type = F7;
        else if(log_temp>=3.793) spectral_type = F8;
        else if(log_temp>=3.777) spectral_type = G0;
        else if(log_temp>=3.770) spectral_type = G1;
        else if(log_temp>=3.763) spectral_type = G2;
        else if(log_temp>=3.748) spectral_type = G5; //MB
        else if(log_temp>=3.708) spectral_type = K0; //MB
        else if(log_temp>=3.623) spectral_type = K5; //MB
        else if(log_temp>=3.568) spectral_type = M0; //MB
        else if(log_temp>=3.477) spectral_type = M5; //MB
        else if(log_temp>=3.398) spectral_type = M8; //MB
        else spectral_type = M8;

    return spectral_type;
}

luminosity_class get_luminosity_class(const real temperature,
				      const real lum) {

  luminosity_class lum_c = V;
  real log_temp = log10(temperature);
  real log_lum  = log10(lum);

  if (log_temp<=4.5 && log_lum>=lum_class_limit(log_temp, I))
    lum_c = I;
  else if (log_temp<=4.4 && log_lum>=lum_class_limit(log_temp,II))
    lum_c = II;
  else if (log_temp<=4.2 && log_lum>=lum_class_limit(log_temp,III))
    lum_c = III;
  else if (log_temp<=4.1 && log_lum>=lum_class_limit(log_temp,IV))
    lum_c = IV;
  
  return lum_c;
}

real lum_class_limit(const real log_temp,
		     luminosity_class lmc) {

    real cfit1, cfit2, cfit3, cfit4;
    
    switch(lmc) {
	case I:         cfit1 = -30.49;
	                cfit2 = 24.21;
	                cfit3 = -6.01;
		        cfit4 = 0.532;
	break;
	case II:        cfit1 = 173.3;
	                cfit2 = -112.2;
	                cfit3 = 23.78;
	                cfit4 = -1.591;
	break;
	case III:       cfit1 = -251.3;
	                cfit2 = 239.8;
	                cfit3 = -73.16;
	                cfit4 = 7.266;
	break;
	case IV:        cfit1 = -89.86;
	                cfit2 = 83.60;
	                cfit3 = -25.62;
	                cfit4 = 2.61;
	break;
	default:        cfit1 = cfit2 = cfit3 = cfit4 = 0;
    }
    
    return cfit1 + log_temp*(cfit2 + log_temp*(cfit3 + log_temp*(cfit4)));
}


stellar_type extract_stellar_type_string(char* star_type_string) {

     stellar_type type = NAS;

     if (!strcmp(star_type_string, "planet")) 
        type = Planet;
     else if (!strcmp(star_type_string, "proto_star")) 
        type = Proto_Star;
     else if (!strcmp(star_type_string, "brown_dwarf")) 
        type = Brown_Dwarf;
     else if (!strcmp(star_type_string, "main_sequence"))
        type = Main_Sequence;
     else if (!strcmp(star_type_string, "hyper_giant")) 
        type = Hyper_Giant;
     else if (!strcmp(star_type_string, "hertzsprung_gap"))
        type = Hertzsprung_Gap;
     else if (!strcmp(star_type_string, "sub_giant"))
        type = Sub_Giant;
     else if (!strcmp(star_type_string, "horizontal_branch"))
        type = Horizontal_Branch;
     else if (!strcmp(star_type_string, "super_giant"))
        type = Super_Giant;
     else if (!strcmp(star_type_string, "thorne_zytkow"))
        type = Thorn_Zytkow;
     else if (!strcmp(star_type_string, "carbon_star"))
        type = Carbon_Star;
     else if (!strcmp(star_type_string, "helium_star"))
        type = Helium_Star;
     else if (!strcmp(star_type_string, "helium_giant"))
        type = Helium_Giant;
     else if (!strcmp(star_type_string, "helium_dwarf"))
        type = Helium_Dwarf;
     else if (!strcmp(star_type_string, "carbon_dwarf"))
        type = Carbon_Dwarf;
     else if (!strcmp(star_type_string, "oxygen_dwarf"))
        type = Oxygen_Dwarf;
     else if (!strcmp(star_type_string, "xray_pulsar"))
        type = Xray_Pulsar;
     else if (!strcmp(star_type_string, "radio_pulsar"))
        type = Radio_Pulsar;
     else if (!strcmp(star_type_string, "neutron_star"))
        type = Neutron_Star;
     else if (!strcmp(star_type_string, "black_hole"))
        type = Black_Hole;
     else if (!strcmp(star_type_string, "Disintegrated"))
        type = Disintegrated;
     else if (!strcmp(star_type_string, "SPZDCH_star")) 
        type = SPZDCH_Star;
     else if (!strcmp(star_type_string, "static_star")) 
        type = Static_Star;
     else
        type = NAS;

     return type;
   }


stellar_type_summary extract_stellar_type_summary_string(char* star_type_string) {

     stellar_type_summary type = ZAMS;

     if (!strcmp(star_type_string, "zams")) 
        type = ZAMS;
     else if (!strcmp(star_type_string, "early_giant")) 
        type = Early_Giant;
     else if (!strcmp(star_type_string, "late_giant"))
        type = Late_Giant;
     else if (!strcmp(star_type_string, "helium_remnant")) 
        type = Helium_Remnant;
     else if (!strcmp(star_type_string, "white_dwarf"))
        type = White_Dwarf;
     else if (!strcmp(star_type_string, "inert_remnant"))
        type = Inert_Remnant;
     else if (!strcmp(star_type_string, "unspecified"))
        type = Unspecified;
     else if (!strcmp(star_type_string, "undefined"))
        type = Undefined;
     else if (!strcmp(star_type_string, "SPZDCH_star"))
        type = Unspecified;
     else if (!strcmp(star_type_string, "proto_star")) 
        type = Unspecified;
     else if (!strcmp(star_type_string, "static_star"))
        type = Unspecified;
     else
        type = no_of_star_type_summ;

     return type;
   }

star_type_spec extract_stellar_spec_summary_string(char* star_spec_string) {

     star_type_spec type = NAC;

     if (!strcmp(star_spec_string, "emission")) 
        type = Emission;
     else if (!strcmp(star_spec_string, "Barium")) 
        type = Barium;
     else if (!strcmp(star_spec_string, "blue_straggler"))
        type = Blue_Straggler;
     else if (!strcmp(star_spec_string, "Roche_lobe_filling")) 
        type = Rl_filling;
     else if (!strcmp(star_spec_string, "runaway"))
        type = Runaway;
     else if (!strcmp(star_spec_string, "merger"))
        type = Merger;
     else if (!strcmp(star_spec_string, "accreting"))
        type = Accreting;
     else if (!strcmp(star_spec_string, "disintegrated"))
        type = Dsntgr;
     else
        type = NAC;

     return type;
   }

 
char* type_string(mass_transfer_type type) {
   
  switch(type) {
    case Unknown:		return "Unknown";
    case Nuclear: 		return "Nuclear";
    case AML_driven:		return "AML_driven";
    case Thermal:		return "Thermal";
    case Dynamic:		return "Dynamic";
    default:			return "unknown_mass_transfer_type";
  }
}

char* type_short_string(mass_transfer_type type) {
   
  switch(type) {
    case Unknown:		return "?";
    case Nuclear: 		return "nuc";
    case AML_driven:		return "aml";
    case Thermal:		return "th";
    case Dynamic:		return "dyn";
    default:			return "???";
  }
}


bool post_supernova_star(stellar_type tpe) {

  if (tpe == Xray_Pulsar  ||
      tpe == Radio_Pulsar ||
      tpe == Neutron_Star ||
      tpe == Black_Hole)
    return true;

  return false;
}

supernova_type type_of_supernova(stellar_type progenitor) {

  supernova_type ns_type = NAT;
    switch(progenitor) {
      case Super_Giant: ns_type = SN_II;
			break;
      case Hyper_Giant: ns_type = SN_Ic;
			break;
      case Helium_Giant: ns_type = SN_Ib;
			break;
      case Thorn_Zytkow: ns_type = SN_II;
			break;
      case Carbon_Dwarf: 
      case Helium_Dwarf:
      case Oxygen_Dwarf:
			ns_type = SN_Ia;
			break;
// Does not work properly because the argument is the progenitor which
// in the case of a neutron star is also a post_supernova_star.
//      case Xray_pulsar: 
//      case Radio_Pulsar: 
//      case Neutron_Star: ns_type = SN_IV;
//			break;
      default:
			ns_type = NAT;
    };
	
  return ns_type;
}

char *type_string(supernova_type sn_type) {

  switch(sn_type) {
    case NAT:    return "nat";
    break;
    case SN_Ia:    return "Ia";
    break;
    case SN_Ib:    return "Ib";
    break;
    case SN_Ic:    return "Ic";
    break;
    case SN_II:    return "II";
    break;
    case SN_IIL:    return "IIL";
    break;
    case SN_IIP:    return "IIP";
    break;
    case SN_IV:    return "IV";
    break;
    default:       return "Unknown supernova_type";
  }
}

bool remmant(stellar_type tpe) {

  if (tpe == Helium_Dwarf ||
      tpe == Carbon_Dwarf ||
      tpe == Oxygen_Dwarf ||
      tpe == Xray_Pulsar  ||
      tpe == Radio_Pulsar ||
      tpe == Neutron_Star ||
      tpe == Black_Hole)
    return true;
  
  return false;
}

#if 0
real main_sequence_mass_loss() {

  static time[] = {   ,   , };
  static mass[] = {   ,   , };
  

}

//              Standard lineair interpolation routine.
real lineair_interpolation(const real x,
			   const real x1, const real x2,
			   const real y1, const real y2) {

#endif






