//
// triple_star.C
//

#include "double_star.h"
#include "single_star.h"
#include "main_sequence.h"


//              Mass transfer utilities.
/*
real get_relative_mass() {
      return get_primary()->get_relative_mass() + get_secondary()->get_relative_mass();
   }
*/

star* double_star::reduce_mass(const real mdot) {

cerr <<"void double_star::reduce_mass(mdot="<<mdot<<")"<<endl;

        real envelope_mass = get_primary()->get_envelope_mass()
                          + get_secondary()->get_envelope_mass();

        if (envelope_mass<=0) 
cerr <<"Really in trouble! envelope_mass = "<<envelope_mass<<endl;

/*
        if (envelope_mass<=mdot) {
           get_primary()->reduce_mass(primary->get_envelope_mass());
           secondary->reduce_mass(secondary->get_envelope_mass());
           return;
        }
*/
cerr<<"reduce_envelope of both components with the fraction: "
    <<mdot/envelope_mass<<endl;


      return this;
     }

star* double_star::subtrac_mass_from_donor(const real dt, real& mdot) {

cerr <<"real double_star::subtrac_mass_from_donor(dt= "<<dt<<")"<<endl;
cerr<<"double_star not allowed to fill Roche-lobe"<<endl;
exit(1);

/*
        real mdot = relative_mass*dt/root->get_donor_timescale();
//        mdot = mass_ratio_mdot_limit(mdot);

        real evelope_mass = primary->get_envelope_mass()
                          + secondary->get_envelope_mass();

*/
        //return 0;
        mdot = 0;
 
      return this;
     }

real double_star::wind_velocity() {
    cerr<<"double_star::wind_velocity()"<<endl;
    real v_wind = 2.5*sqrt(cnsts.physics(G)*cnsts.parameters(solar_mass)*get_total_mass()
			   / (get_effective_radius()*cnsts.parameters(solar_radius)))
                  / cnsts.physics(km_per_s);
    return v_wind;

}
void double_star::stellar_wind(const real dt) {
cerr<<"void double_star::stellar_wind(const real "<<dt<<")"<<endl;
      perform_wind_loss(dt);
}


real double_star::temperature() {
//cerr<<"real double_star::temperature()"<<endl;

      real temp = get_primary()->temperature()
                + get_secondary()->temperature();

      return temp;
   }

real double_star::bolometric_correction() {
cerr<<"real double_star::bolometric_correction()"<<endl;

      real bc = get_primary()->bolometric_correction()
              + get_secondary()->bolometric_correction();

      return bc;
   } 

real double_star::mass_ratio_mdot_limit(real mdot) {
cerr<<"real double_star::mass_ratio_mdot_limit(mdot="<<mdot<<")"<<endl;

        real accretor_mass = get_companion()->get_total_mass();
           if (accretor_mass<get_total_mass()) {
              real mdot_max = get_total_mass() - accretor_mass;
              if (mdot>mdot_max) mdot = mdot_max;
           }

        return mdot;
     }

real double_star::kelvin_helmholds_timescale() {
cerr<<"real double_star::kelvin_helmholds_timescale()"<<endl;

        return 1;
     }

real double_star::nucleair_evolution_timescale() {
cerr<<"real double_star::nucleair_evolution_timescale()"<<endl;


        return cnsts.physics(nucleair_efficiency);
     }

real double_star::dynamic_timescale() {
cerr<<"real double_star::dynamic_timescale()"<<endl;

        return 1;
     }


real double_star::accretion_limit(const real mdot, const real dt) {
cerr<<"real double_star::accretion_limit(mdot="<<mdot<<", dt="<<dt<<")"<<endl;

        real r_rl = get_binary()->roche_radius(this);
cerr <<"rl "<<r_rl ;
        real prim_mdot = get_primary()->accretion_limit(mdot, dt);
        real sec_mdot  = get_secondary()->accretion_limit(mdot, dt);
cerr<<"mdot_p="<<prim_mdot<<" mdot_s="<<sec_mdot<<endl;
        real accretion = prim_mdot + sec_mdot;

        return accretion;

}
/*
//		Should be inplemented in base_element.h
void double_star::adjust_accretor_radius(const real mdot, const real dt) {
cerr<<"d double_star::adjust_accretor_radius(mdot="<<mdot
    <<", dt="<<dt<<")"<<endl;

      primary->adjust_accretor_radius(mdot, dt);
      secondary->adjust_accretor_radius(mdot, dt);
   }
*/

real double_star::mass_transfer_timescale(mass_transfer_type &type) {
cerr<<"real double_star::mass_transfer_timescale()"<<endl;

      real mtt_s = get_secondary()->mass_transfer_timescale(type);
      real mtt_p = get_primary()->mass_transfer_timescale(type);

      real mtt = 1./(1./mtt_p + 1./mtt_s);
      
      return mtt;
}

real double_star::zeta_adiabatic() {
    cerr<<"real double_star::zeta_adiabatic()"<<endl;
    return 0;
}

real double_star::zeta_thermal() {
    cerr<<"real double_star::zeta_thermal()"<<endl;
    return 0;
}

real double_star::add_mass_to_accretor(real mdot, bool hydrogen, const real dt) {
cerr<<"real double_star::add_mass_to_accretor(mdot="<<mdot
    <<", dt="<<dt<<")"<<endl;

      if (bin_type==Merged || bin_type==Disrupted) 
         return 0;

      if (mdot<=0) {
         cerr << "double_star::add_mass_to_accretor(mdot=" << mdot << ")"<<endl;
         cerr << "mdot (mdot=" << mdot << ", dt="<<dt
              <<"): mdot smaller than zero!" << endl;
         cerr << "Action: No action!" << endl;
         return 0;
      }

      mass_transfer_type type = Unknown;
      real mtt_p = get_primary()->mass_transfer_timescale(type);
      real mtt_s = get_secondary()->mass_transfer_timescale(type);
      real mtt_b = mass_transfer_timescale(type);

      real mdot_p = mdot*mtt_b/mtt_p;
      real mdot_s = mdot*mtt_b/mtt_s;

      mdot_p = get_primary()->add_mass_to_accretor(mdot_p, dt, hydrogen);
      mdot_s = get_secondary()->add_mass_to_accretor(mdot_s, dt, hydrogen);

      real mdot_left = mdot - (mdot_p + mdot_s);

//              Spiral in anyway (temporary treatment)
       real r_lobe = Starlab::min(get_primary()->get_effective_radius()/semi,
                         get_secondary()->get_effective_radius()/semi);
       real a_spi = semi*(get_primary()->get_total_mass()
                          /(get_primary()->get_total_mass()+mdot))
                          /(1. + (2.*mdot_left
                   /      (cnsts.parameters(common_envelope_efficiency)
		   *       cnsts.parameters(envelope_binding_energy)
		   *       r_lobe*get_secondary()->get_total_mass())));

       if (a_spi <= get_primary()->get_effective_radius()
                 + get_secondary()->get_effective_radius())
          merge_elements(get_primary(), get_secondary());
       else
          semi = a_spi;

      return mdot_p + mdot_s;

}
real double_star::accrete_from_stellar_wind(const real mdot, const real dt) {

    // cerr<<"real double_star::accrete_from_stellar_wind(mdot="<<mdot
    //     <<", dt="<<dt<<")"<<endl;

    //		Just try to prevent error in star::accrete_from_stellar_wind()

      if(bin_type!=Merged && bin_type!=Disrupted) {
      real mdot_p = get_primary()->accrete_from_stellar_wind(mdot, dt);
      real mdot_s = mdot - mdot_p;
      if (mdot_s>0) 
         get_secondary()->accrete_from_stellar_wind(mdot_s, dt);
      }

      return 0;		// Added by Steve 8/4/97 to satisfy DEC C++
  }

void double_star::set_rotational_velocity(const real v) {
    cerr<<"void double_star::set_rotational_velocity(v="<<v<<")"<<endl;
}

void double_star::set_previous_radius(const real r) {
    //cerr<<"void double_star::set_previous_radius(r="<<r<<")"<<endl;
}

void double_star::adjust_triple_after_wind_loss(star * donor,
                                                const real md_dot,
                                                const real dt) {

cerr <<"void double_star::adjust_triple_after_wind_loss(d="
     <<donor<<", m.="<<md_dot<<" dt="<<dt<<"): "<<semi;
     if (bin_type!=Merged && bin_type!=Disrupted) {
        star* accretor = donor->get_companion();

     real M_old = get_total_mass() + md_dot;
     real old_donor_mass = donor->get_total_mass() + md_dot;
     real old_accretor_mass = accretor->get_total_mass();
cerr<<"m's: "<<M_old<<" "<<old_donor_mass<<" "<<old_accretor_mass
             <<" "<<md_dot<<endl;

     real ma_dot = accretor->accrete_from_stellar_wind(md_dot, dt);

     real M_new = get_total_mass();
     real new_donor_mass = donor->get_total_mass();
     real new_accretor_mass = accretor->get_total_mass();
cerr<<"M's: "<<M_new<<" "<<new_donor_mass<<" "<<new_accretor_mass
             <<" "<<md_dot<<endl;

cerr<<"mdot: "<<md_dot<<" "<<ma_dot<<endl;
     if(md_dot>0 && ma_dot>=0) {
        real alpha = 1 - ma_dot/md_dot;
        semi *= pow(pow(new_donor_mass/old_donor_mass, 1-alpha)
             *  new_accretor_mass/old_accretor_mass, -2)
             *  M_old/M_new;
     }
     }
     else  {
//        base_element * old_d = donor;
        donor->reduce_mass(md_dot);
//        if (old_d==get_primary()) donor=primary;
//        else donor=secondary;
     }
cerr<<" a-->"<<semi<<endl;
}


// What is the gyration radius of a binary?
real double_star::gyration_radius_sq() {
  cerr << "double_star::gyration_radius_sq()"<<endl;
  cerr << "Unknown: return 0"<<endl; 

  return 0; // Dummy
}
