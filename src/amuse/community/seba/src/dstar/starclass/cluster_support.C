
#include "cluster_support.h"
//#include "general_support.h"
//#include "star_support.h"

/*-----------------------------------------------------------------------------
 *  make_profile --  make profile for whole n-binary cluster
 *-----------------------------------------------------------------------------
 */
void  make_profile(initial_cluster& c_init,
                   cluster_profile& c_prof, profile prof,
                   star_type_spec s_prof = Emission) {

      double_profile binary;
      double_init b_init;

      if (!c_init.field) {
      	cerr<<"assuming solar metallicity in make_profile in cluster_support.C"<<endl;
         real m_to = turn_off_mass(c_init.end_time, 0.02);
         cout << "Turn off mass (" << type_string(get_spectral_class(m_to))  
              << ") = " << turn_off_mass(c_init.end_time, 0.02) 
              << " Solar mass."<< endl;
      }

      for (int i=0; i<c_init.no_of_binaries; i++) {
          b_init = c_init.random_initial_conditions();
          make_profile(i+1, c_init.start_time, binary, b_init);
          c_prof.enhance_cluster_profile(binary, prof, s_prof);
      }

}    

/*-----------------------------------------------------------------------------
 *  cluster_profile -- Support routines for cluster simulation.
 *-----------------------------------------------------------------------------
 */

cluster_profile::cluster_profile() {

        no_of_init_bins=no_of_fin_bins=0;
        no_of_runners=no_of_mergers=no_of_singles=0;

        for (int s=MIN_TYPE_NUMBER; s<MAX_TYPE_NUMBER; s++) {
            for (int p=MIN_TYPE_NUMBER; p<MAX_TYPE_NUMBER; p++) {
                init_bins[p][s] = 0;
                fin_bins[p][s] = 0;
                bins_stt[p][s] = 0;
                bins_spt[p][s] = 0;
                for (int t = NAC; t<no_of_spec_type; t++)
                    bins_spa[p][s][t] = 0;
            }
            singles[s] = 0;
            mergers[s] = 0;
            single_stt[s] = 0;
            single_spt[s] = 0;
            merge_stt[s] = 0;
            merge_spt[s] = 0;
            runner_stt[s] = 0;
            runner_spt[s] = 0;
            for (int t = NAC; t<no_of_spec_type; t++) {
                 single_spa[s][t] = 0;
                 merge_spa[s][t] = 0;
                 runner_spa[s][t] = 0;
            }
        }

    }

/*-----------------------------------------------------------------------------
 *  enhance_cluster_profile -- Support routines for cluster profiler.
 *-----------------------------------------------------------------------------
 */
void cluster_profile::enhance_cluster_profile(double_profile& binary,
                                              profile prof,
                                              star_type_spec s_prof = Emission) {

        if (prof==star_type) 
           star_type_cluster_profile(binary);
        else if (prof==spectral_type) 
           spectral_type_cluster_profile(binary);
        else if (prof==spectral_addition) 
           spectral_addition_cluster_profile(binary, s_prof);
        else {
           cerr << "cluster_profile::enhance_cluster_profile()" <<endl;
           cerr << "profile type unknown" << endl;
        }
    }

void cluster_profile::enhance_cluster_profile(double_profile& binary) {

      star_type_cluster_profile(binary);
      spectral_type_cluster_profile(binary);
      for (int prof = Emission; prof<no_of_spec_type; prof++)
          spectral_addition_cluster_profile(binary, (star_type_spec)prof);
    }

void cluster_profile::enhance_cluster_profile(star_state& single) {

      star_type_cluster_profile(single);
      spectral_type_cluster_profile(single);
      for (int prof = Emission; prof<no_of_spec_type; prof++)
          spectral_addition_cluster_profile(single, (star_type_spec)prof);
    }

void cluster_profile::star_type_cluster_profile(double_profile& binary) {

        init_bins[binary.init.primary.type][binary.init.secondary.type]++;
        no_of_init_bins++;
        
        if (binary.final.type==Merged) { 
           merge_stt[binary.final.primary.type]++; 
        }
        else if (binary.final.type==Disrupted) {
           runner_stt[binary.final.primary.type]++; 
           runner_stt[binary.final.secondary.type]++; 
        }
        else {
           bins_stt[binary.final.primary.type]
                   [binary.final.secondary.type]++;
        } 
    } 

void cluster_profile::star_type_cluster_profile(star_state& single) {

      if (single.class_spec[Merger])
         merge_stt[single.type]++;
      else if (single.class_spec[Runaway]) 
         runner_stt[single.type]++;
      else
         single_stt[single.type]++;
   }

void cluster_profile::spectral_type_cluster_profile(double_profile& binary) {

//put_state(binary.init);
        init_bins[binary.init.primary.class_tpe]
                 [binary.init.secondary.class_tpe]++;
       
        if (binary.final.type==Merged) {
           merge_spt[binary.final.primary.class_tpe]++;
        }
        else if (binary.final.type==Disrupted) {
           runner_spt[binary.final.primary.class_tpe]++;
           runner_spt[binary.final.secondary.class_tpe]++;
        }
        else {
           bins_spt[binary.final.primary.class_tpe]
                   [binary.final.secondary.class_tpe]++;
        }
    }

void cluster_profile::spectral_type_cluster_profile(star_state& single) {

      if (single.class_spec[Merger])
         merge_spt[single.class_tpe]++;
      else if (single.class_spec[Runaway])
         runner_spt[single.class_tpe]++;
      else
         single_spt[single.class_tpe]++;
      no_of_singles++;
    }

void cluster_profile::spectral_addition_cluster_profile(double_profile& binary,
                                                  star_type_spec s_prof) {

//        if (s_prof==Runaway || s_prof==Merger) {
//           spectral_addition_single_profile(binary, s_prof);
//           return ;
//        }

//		For spectral_addition_cluster_profile()
// 		init_bins contains te final systems which containd
//		a spectral_addition star as primary.
        if (binary.final.primary.class_spec[s_prof]) {
           init_bins[binary.final.primary.class_tpe]
                 [binary.final.secondary.class_tpe]++;
           no_of_init_bins++;
        }
       
        if (binary.final.type==Disrupted) {
           if (binary.final.primary.class_spec[s_prof]) {
              runner_spa[binary.final.primary.class_tpe][s_prof]++;
           }
           if (binary.final.secondary.class_spec[s_prof]) {
              runner_spa[binary.final.secondary.class_tpe][s_prof]++;
           }
        }
        else if (binary.final.type==Merged) {
           if (binary.final.primary.class_spec[s_prof]) {
              merge_spa[binary.final.primary.class_tpe][s_prof]++;
           }
        }
        else {
//              fin_bins contains te final systems which containd
//              a spectral_addition star as secondary.
           if (binary.final.primary.class_spec[s_prof]) {
              bins_spa[binary.final.primary.class_tpe]
                      [binary.final.secondary.class_tpe][s_prof]++;
           }
           else if (binary.final.secondary.class_spec[s_prof]) {
              bins_spa[binary.final.secondary.class_tpe]
                      [binary.final.primary.class_tpe][s_prof]++;
           }
        }
    }

void cluster_profile::spectral_addition_cluster_profile(star_state& single,
                                                  star_type_spec prof) {
      if (single.class_spec[prof]) {
         if (single.class_spec[Merger])
            merge_spa[single.class_tpe][prof]++;
         else if (single.class_spec[Runaway])
            runner_spa[single.class_tpe][prof]++;
         else
            single_spa[single.class_tpe][prof]++;
      }
   }


void cluster_profile::spectral_addition_single_profile(double_profile& binary,
                                                  star_type_spec s_prof) {

//              For spectral_addition_single_profile()
//              fin_bins primary (row) contain stellar type information.
//	  	secondary contains spectral information.	
        if (binary.final.primary.class_spec[s_prof]) 
           bins_spa[binary.final.primary.type]
                   [binary.final.secondary.class_tpe][s_prof]++;
        else if (binary.final.secondary.class_spec[s_prof]) 
           bins_spa[binary.final.secondary.type]
                   [binary.final.primary.class_tpe][s_prof]++;
   }

void cluster_profile::print_profile() {


     print_profile(Main_Sequence);
     print_profile(O5);
//     print_profile(O0_O9);
     for (int type=Emission; type<no_of_spec_type; type++) 
         print_profile((star_type_spec)type,star_type);

   }

void cluster_profile::print_profile(profile prof) {

        if (prof==star_type) {
           cout << "Stellar type data:" << endl;
           print_profile(Main_Sequence);
        }
        else if (prof==spectral_type) {
           cout << "Spectral type data:" << endl;
           print_profile(O5);
//           print_profile(O0_O9);
        }
        else if (prof==spectral_addition) {
           print_profile(Emission, prof);
        }
        else {
           cerr << "cluster_profile::print_profile()" <<endl;
           cerr << "profile type unknown" << endl;
        }
    }

void cluster_profile::print_profile(stellar_type type) {

//		Initialize output counters.
//        stellar_type s, p, v_min, v_max;
        int s, p;
        int v_min = Main_Sequence;
        int v_max = no_of_stellar_type;

        cout << endl << endl;
//		Print final conditions.
        cout <<"\t";
        for (s=v_min; s<v_max; s++) {
//          if (s_tot[s]) 
            cout << type_short_string((stellar_type)s) <<"\t";
        }

        cout << endl;
        for (s=v_min; s<v_max; s++) {
//          if (s_tot[s]) {
            cout << type_short_string((stellar_type)s) << "\t";
            for (p=v_min; p<v_max; p++) {
//              if (p_tot[p]) 
                cout << bins_stt[p][s] << "\t";
            }
            cout << endl;
        }

        cout << endl;
//              Print final conditions.
        cout <<"\t";
        for (s=v_min; s<v_max; s++)
            cout << type_short_string((stellar_type)s) <<"\t";
        cout << endl;
        cout << "runaway\t";
        for (s=v_min; s<v_max; s++)
            cout << runner_stt[s] << "\t";
        cout << endl;
        cout << "mergers\t";
        for (s=v_min; s<v_max; s++) 
            cout << merge_stt[s] << "\t";
        cout << endl;
        cout << "singles\t";
        for (s=v_min; s<v_max; s++) 
            cout << single_stt[s] << "\t";
        cout << endl;
    }



void cluster_profile::print_profile(spectral_class type) {

//              Initialize output counters.
        // spectral_class s, p, v_min, v_max;
        int s, p;
        int v_min = O5; //O0_O9;
        int v_max = no_spectral_class;

//              Print final conditions.
        cout << endl;
        cout <<"\t";
        for (s=v_min; s<v_max; s++)
            cout << type_string((spectral_class)s) <<"\t";
        cout << endl;
        for (s=v_min; s<v_max; s++) {
            cout << type_string((spectral_class)s) << "\t";
            for (p=v_min; p<v_max; p++)
                cout << bins_spt[p][s] << "\t";
            cout << endl;
        }


        cout << endl;
//              Print final conditions for single stars.
        cout <<"\t";
        for (s=v_min; s<v_max; s++)
            cout << type_string((spectral_class)s) <<"\t";
        cout << endl;
        cout << "runaway\t";
        for (s=v_min; s<v_max; s++)
            cout << runner_spt[s] << "\t";
        cout << endl;
        cout << "mergers\t";
        for (s=v_min; s<v_max; s++)
            cout << merge_spt[s] << "\t";
        cout << endl;
        cout << "singles\t";
        for (s=v_min; s<v_max; s++)
            cout << single_spt[s] << "\t";
        cout << endl;
    }

void cluster_profile::print_profile(star_type_spec type, profile prof) {

        cout << endl;
        cout << type_string(type) << ": " << endl;
        if (type==Runaway || type==Merger) {
//           print_single_profile(type);
           return ;
        }

//              Initialize output counters.
        // spectral_class s, p, v_min, v_max;
        int s, p;
        int v_min = O5; //O0_O9;
        int v_max = no_spectral_class;

//              Print final conditions.
        cout <<"\t";
        for (s=v_min; s<v_max; s++)
            cout << type_string((spectral_class)s) <<"\t";
        cout << endl;
        for (s=v_min; s<v_max; s++) {
            cout << type_string((spectral_class)s) << "\t";
            for (p=v_min; p<v_max; p++)
                cout << bins_spa[p][s][type] << "\t";
            cout << endl;
        }

        if (type!=Rl_filling && type!=Accreting) {
        cout << endl;
//              Print final conditions.
        cout <<"\t";
        for (s=v_min; s<v_max; s++)
            cout << type_string((spectral_class)s) <<"\t";
        cout << endl;
        cout << "runaway\t";
        for (s=v_min; s<v_max; s++)
            cout << runner_spa[s][type] << "\t";
        cout << endl;
        cout << "mergers\t";
        for (s=v_min; s<v_max; s++)
            cout << merge_spa[s][type] << "\t";
        cout << endl;
        cout << "singles\t";
        for (s=v_min; s<v_max; s++)
            cout << single_spa[s][type] << "\t";
        cout << endl;
        }
    }

void cluster_profile::print_single_profile(star_type_spec type) {

//              Initialize output counters.
        //spectral_class p, p_min, p_max;
        //stellar_type s, s_min, s_max;
        int s, p;
        int p_min = O5; //O0_O9;
        int p_max = no_spectral_class;
        int s_min = Main_Sequence;
        int s_max = no_of_stellar_type;


//              Print initial conditions.
        cout <<"\t";
        for (s=s_min; s<s_max; s++)
            cout << type_short_string((stellar_type)s) <<"\t";
        cout << endl;
        for (p=p_min; p<p_max; p++) {
            cout << type_string((spectral_class)p) << "\t";
            for (s=s_min; s<s_max; s++)
                cout << bins_spa[s][p][type] << "\t";
            cout << endl;
        }
}

/*-----------------------------------------------------------------------------
 *  initial_cluster -- Initialization for cluster and field systems. 
 *-----------------------------------------------------------------------------
 */
initial_cluster::initial_cluster() {

    start_time = 0;
    end_time = 12000;
    field = FALSE;

    start_id = 1;
    n_steps  = 10;
    no_of_singles  = 1000;
    no_of_binaries = 100;

    r_core = 0.5;
    r_halfm = 5;
    r_tidal = 50;
    rho_core = 1.e+5;
    v_disp = 30;

    m_min = 0.08;
    m_max = 100.0;
//		Scalo (1986)
    m_alpha = 2.7;
//	Used for model D9.
//    m_alpha = 2.27;
//    m_alpha = 2.35;

    q_min = 0.3;
    q_max = 0.95;
    q_alpha = 0.7;

    a_min = 0.6;
    a_max = 1.e+4;
    a_alpha = 0.0; //1.3;

    e_min = 0.035;
    e_max = .7;
    e_alpha = 1;

}

double_init initial_cluster::random_initial_conditions() {

//      Should be Sake's initial conditions for spectroskopic
//      binary stars.
//      Literature: Hogeveen S. PhD Thesis Amsterdam 1991.
//
//              Currently standard initialization of
//              Pols 1993.
//
    double_init init;

//              Start time.
    init.start_time = start_time;
    if (!field)
       init.end_time = end_time;
    else
       init.end_time = randinter(start_time, end_time);
    init.n_steps = n_steps;

//              Primary mass.
    real constant = pow(m_min, 1-m_alpha) - pow(m_max, 1-m_alpha);
    real random = randinter(0., 1.);
    init.mass_prim = pow(random*constant
                   + pow(m_max, 1-m_alpha), 1/(1-m_alpha));

//              Eccentricity.
    random = randinter(0., 1.);  // Heggie 1975
    init.eccentricity = sqrt(random); // F(e) = 2e

//              Mass ratio
    random = randinter(0.125, 1.);
    real q = init.q = 1./pow(random, 1/3.) - 1;

//              Semi major axis
    a_min = 0.49*pow(q, 2./3.)/(0.6*pow(q, 2./3.) + log(1 + pow(q, 1/3.)));
    a_min = 2*pow(init.mass_prim, 0.7)/a_min;
    a_max = 5000*a_min;
    constant = log(a_max) - log(a_min);
    random = randinter(0., 1.);
    init.semi = a_min*exp(random*constant);


    return init;
}

double_init random_initial_conditions(initial_cluster& c) {

    extern real randinter(real, real);

    real random, constant;

    double_init init;

//              Start time.
    if (!c.field)
       init.end_time = c.end_time;
    else
       init.end_time = randinter(c.start_time, c.end_time);
    init.n_steps = c.n_steps;

//              Primary mass.
    constant = pow(c.m_min, 1-c.m_alpha) - pow(c.m_max, 1-c.m_alpha);
    random = randinter(0., 1.);
    init.mass_prim = pow(random*constant
                   + pow(c.m_max, 1-c.m_alpha), 1/(1-c.m_alpha));

//              Eccentricity.
    random = randinter(0., 1.);  // Heggie 1975
    init.eccentricity = sqrt(random); // F(e) = 2e

//              Mass ratio
    random = randinter(0.125, 1.);
    real q = init.q = 1./pow(random, 1/3.) - 1;

//              Semi major axis
    c.a_min = 0.49*pow(q, 2./3.)/(0.6*pow(q, 2./3.) + log(1 + pow(q, 1/3.)));
    c.a_min = 2*pow(init.mass_prim, 0.7)/c.a_min;
    c.a_max = 5000*c.a_min;
    constant = log(c.a_max) - log(c.a_min);
    random = randinter(0., 1.);
    init.semi = c.a_min*exp(random*constant);

    return init;

}

real next_output_time(int n_out, int n_tot, real t1, real t2) {
//cerr<<"real next_output_time(nj="<<n_out<<", nt="<<n_tot<<", t1="<<t1
//                                        <<", t2="<<t2<<")"<<endl;
cerr<<"assuming solar metallicity (3x) in next_output_time in cluster_support.C"<<endl;
      real m1 = turn_off_mass(t1, 0.02);
      real m2 = turn_off_mass(t2, 0.02);
      real dm = (m1-m2)/n_tot;
      real dt = (t2-t1)/n_tot;

      real m_out = m1 - dm*(n_out+1);
//      real t_out = t1+dt*(n_out+1);
   
      return main_sequence_time(m_out, 0.02);

   }

     
