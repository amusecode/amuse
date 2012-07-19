
//// rdc_movie: read input snapshot and reduce it to simple readable
////            ASCII for movie program.
////          
//// Options:     none
//-----------------------------------------------------------------------------
//   version 1:  Sept 1998   Simon Portegies Zwart   spz@grape.c.u-tokyo.ac.jp
//                                                   University of Tokyo
//.............................................................................
//   non-local functions: 
//-----------------------------------------------------------------------------

#include "single_star.h"
#include "sstar_to_dyn.h"
#include "dyn.h"

#ifndef TOOLBOX

#define NNEBULAE 10

typedef struct nebulae {
  int id;
  real time;
  vec pos;
};
nebulae nebula[NNEBULAE];

#define NSN_REMNANT 4
typedef struct sn_remnants {
  int id;
  bool bh;
  real time;
  vec pos;
};
sn_remnants sn_remnant[NSN_REMNANT];

#define NCOLLISIONS 10
typedef struct collisions {
  int starid;
  char name[25];
  real time;
};
collisions collision[NCOLLISIONS];

enum filter_type {Visual=0, Radio, X_rays}; 

local filter_type get_filter_type(char fltr) {

  switch(fltr) {
    case 'V': return Visual;
              break;
    case 'R': return Radio;
              break;
    case 'X': return X_rays;
              break;
    default: cerr << "No such filter"<<endl;
             exit(-1);
  };

}

void new_camera_position(vec &v_new,
			 vec first_pos,
			 vec last_pos,
			 int nsteps,
			 int nsnap) {

  real fstep = (1.0*nsnap)/(1.0*nsteps);
  v_new[0] = first_pos[0] + fstep * (last_pos[0] - first_pos[0]);
  v_new[1] = first_pos[1] + fstep * (last_pos[1] - first_pos[1]);
  v_new[2] = first_pos[2] + fstep * (last_pos[2] - first_pos[2]);
  
}

#if 0
void new_camera_position(vec &v_new,
			 vec v_old,
			 real theta_rotation,
			 real phi_rotation,
			 int nsteps,
			 int nsnap) {

  
  real theta = nsnap*theta_rotation/nsteps;
  real phi   = nsnap*phi_rotation/nsteps;

  PRC(nsnap);PRC(theta);PRL(phi);
  real r = abs(v_old);
  v_new[0] = r * sin(theta);
  v_new[1] = r * sin(theta) * cos(phi);
  v_new[2] = r * cos(theta);

  cerr << "vector = " << v_new << endl;

}
#endif

local void print_camera_on_star(dyn* b, vec com_pos,
				int camera_on_star_id,
                                real r_star,
				real aperture, int blur_samples) {

    vec cam_view = b->get_pos() + b->get_vel();
    vec cam_pos = b->get_pos() + r_star*b->get_vel();

    // Look at the clusters CoM for now.
//    cam_view[0] = cam_pos;
    cam_view[0] = 0;
    cam_view[1] = 0;
    cam_view[2] = 0;
    
    vec normal;
    normal[0] = sqrt(pow(abs(cam_pos), 2)
              -      pow(cam_pos[1], 2)
              -      pow(cam_pos[2], 2));
    normal[1] = 1;
    normal[2] = (-cam_pos[0] * normal[0] - cam_pos[1])/cam_pos[2];

    cout << "// Normal to camera " << endl;
    cout << "   #declare normal_to_camera = < " << normal[0] << ", "
	                                        << normal[1] << ", "
	                                        << normal[2] << " >" << endl;
    
    cout << "// camera located on star #" << camera_on_star_id << endl;
    cout << "camera {" << endl
	 << "   location < " << cam_pos[0] << ", "
                             << cam_pos[1] << ", "
                             << cam_pos[2] << " >" << endl
	 << "   look_at  < " << cam_view[0] << ", "
	                     << cam_view[1] << ", "
	                     << cam_view[2] << " >" << endl
	 << "   blur_samples " << blur_samples << endl; 
    if (aperture>0)
	cout << "   focal_point < " << cam_view[0] << ", "
	                            << cam_view[1] << ", "
	                            << cam_view[2] << " >" << endl
	     << "   aperture " << aperture << endl;
    cout << "}" << endl << endl;

}


local bool print_camera_position_recursive(dyn* b, vec cam_pos,
					   int camera_on_star_id,
                                           real r_star,
					   real aperture, int blur_samples) {

  if (b->get_oldest_daughter()) {
    
    vec com_pos  = b->get_pos() - cam_pos;


    for_all_daughters(dyn, b, bb)
      if (bb->n_leaves() >= 2) 
	return print_camera_position_recursive(bb, com_pos,
                                               camera_on_star_id, r_star,
                                               aperture, blur_samples);

      else if (bb->get_index()==camera_on_star_id) {
	print_camera_on_star(bb, com_pos, 
                             camera_on_star_id, r_star,
                             aperture, blur_samples);

	return true;
      }

  }
  return false;
}

local void print_filename_counter(int counter, ostream& s) {

    if (counter>=10000) { 
      cout << "\nToo namy filenames in print_povray_header()" << endl;
      exit(1);
    }
    else if (counter<10)
      s << "000" << counter;
    else if (counter<100)
      s << "00" << counter;
    else if (counter<1000)
      s << "0" << counter;
    else
      s << counter;
}

local void print_pl_nebula(vec pos, real scale) {

    cout << "object { Pl_Nebula scale " << scale
	 << " translate < " << pos[0] << ", "
                            << pos[1] << ", "
		            << pos[2] << " > }"
	 << endl;
}

local void print_sn_nebula(vec pos, real scale) {

  cout << "object { SN_Remnant scale " << scale
       << " translate < " << pos[0] << ", "
                          << pos[1] << ", "
		          << pos[2] << " > }"
       << endl;
}

local void print_some_data(dyn *b) {

  if (b->get_oldest_daughter()) {

    real time = b->get_starbase()->conv_t_dyn_to_star(b->get_real_system_time());
    
    int nts=0, nbs=0, nss=0;
    for_all_daughters(dyn, b, bb)
      if (bb->n_leaves() > 2)
	nts++;
      else if (bb->n_leaves() >= 2)
	nbs++;
      else 
	nss++;
  }
}

local void print_hertzsprung_Russell_diagram(dyn* b, vec cam_pos) {

  if (b->get_oldest_daughter()) {

    print_some_data(b);

    int stype_s[no_of_star_type_summ];
    for (int i=0; i<no_of_star_type_summ; i++) 
      stype_s[i]=0;
    
    cout << "#declare Stellar_HRD = union {" << endl;
    cout << "   object { Y_Axis } " << endl;
    cout << "   object { X_Axis }\n " << endl;
    
    for_all_leaves(dyn, b, bi) {
      
      star_type_spec tpe_class = NAC;
      spectral_class star_class;
      stellar_type stype = NAS;
      stellar_type_summary sstype = ZAMS;
      real t_cur, t_rel, m_rel, m_env, m_core, co_core; 
      real T_eff, L_eff, p_rot, b_fld;
	//real T_eff=0, L_eff=0;
      if (bi->get_use_sstar()) {
	stype = bi->get_starbase()->get_element_type();
	sstype = summarize_stellar_type(stype);
	stype_s[dynamic_cast(int, sstype)]++;
	
	T_eff = bi->get_starbase()->temperature();
	L_eff = bi->get_starbase()->get_luminosity();
      }
      else if (bi->get_star_story()) {
	extract_story_chapter(stype, t_cur, t_rel, m_rel, m_env, 
			      m_core, co_core,
			      T_eff, L_eff, p_rot, b_fld,
			      *bi->get_star_story());
	sstype = summarize_stellar_type(stype);
	stype_s[dynamic_cast(int, sstype)]++;
	star_class = get_spectral_class(T_eff);
#if 0
	if (find_qmatch(bi->get_star_story(), "Type")) {
	cerr <<"Reading"<<endl;
	stype = extract_stellar_type_string(
                getsq(bi->get_star_story(), "Type"));
	sstype = summarize_stellar_type(stype);
	stype_s[dynamic_cast(int, sstype)]++;
	T_eff = getrq(bi->get_star_story(), "T_eff");
	star_class = get_spectral_class(T_eff);
	L_eff = getrq(bi->get_star_story(), "L_eff");
#endif	
      }
      else {
	cout << "    No stellar information found for: ";
	bi->pretty_print_node(cout);
	return;
      }
      
      real xt_pos = (4.78 - log10(T_eff))/4.78;
      real yl_pos = (log10(L_eff) + 1)/7.;

      if (xt_pos>0&&xt_pos<1 && yl_pos>0&&yl_pos<1)
	cout << "   object { Red_Sphere translate < "
	     << xt_pos << ", "
	     << yl_pos << ", "
	     << "0 > }" << endl;
    }

    cout << "\n} // End Stellar_HRD" << endl;
    cout << "object { Stellar_HRD translate < "
	 << cam_pos[0]+0.5 << ", "
	 << cam_pos[1]-1 << ", "
	 << cam_pos[2] + 2
	 << " > }\n" << endl;

  }
}

local void add_collision_effect(dyn *bi, vec pos, real time, 
				real scale) { 

  for (int i=0; i<NCOLLISIONS; i++) {
    if(bi->get_name()) {
      if(!strcmp(bi->get_name(), collision[i].name) && 
	 time >= collision[i].time && time < collision[i].time+0.125) {
	 cerr << "Adding collision at time = " << time << endl;
	 cerr << "object { OStar scale "
	      << 4 * scale * pow(5 * (0.125 - (time-collision[i].time)), 3)
	      << " translate < "
	      << pos[0] << ", " << pos[1] << ", " << pos[2] << " > }"
	      << endl;
	 if(collision[i].starid == 0) { // RLOF
	   cout << "object { KStar scale ";
	 } 
	 else { // Collision
	   cout << "object { OStar scale ";
	 }
	   cout << 4 * scale * pow(5 * (0.125 - (time-collision[i].time)), 3)
	     << " translate < "
	       << pos[0] << ", " << pos[1] << ", " << pos[2] << " > }"
		 << endl;
       }
     }
  }
}

local void print_star(dyn *bi, vec pos,
		      real scale_L, filter_type filter) {

  // To solar radii
  //  vec pos = bi->get_pos() - dc_pos;
//  pos[0] = bi->get_starbase()->conv_r_dyn_to_star(pos[0]);
//  pos[1] = bi->get_starbase()->conv_r_dyn_to_star(pos[1]);
//  pos[2] = bi->get_starbase()->conv_r_dyn_to_star(pos[2]);

  real time = bi->get_starbase()->conv_t_dyn_to_star(bi->get_real_system_time());

  // And now to parsec
  real Rsun_per_parsec = cnsts.parameters(solar_radius)
                       / cnsts.parameters(parsec);
//  pos[0] *= Rsun_per_parsec;
//  pos[1] *= Rsun_per_parsec;
//  pos[2] *= Rsun_per_parsec;

     star_type_spec tpe_class = NAC;
     spectral_class star_class;
     stellar_type stype = NAS;
     stellar_type_summary sstype = ZAMS;
     real t_cur, m_rel, m_env, m_core, co_core, T_eff, L_eff, p_rot, b_fld;
     real t_rel=0, R_eff=0;
     real M_tot, U, B, V, R, I;	

     if (bi->get_star_story()) {
//       if (find_qmatch(bi->get_star_story(), "Type")) {
	  // cerr << "Reading"<< endl;
	  // put_dyn(bi, cerr);
	  stype = extract_stellar_type_string(
		      getsq(bi->get_star_story(), "Type"));
	  M_tot = getrq(bi->get_star_story(), "M_rel");
	  T_eff = getrq(bi->get_star_story(), "T_eff");
	  star_class = get_spectral_class(T_eff);
	  L_eff = getrq(bi->get_star_story(), "L_eff");
//	}
//       extract_story_chapter(stype, t_cur, t_rel, m_rel, m_env, 
//			     m_core, co_core,
//			     T_eff, L_eff, p_rot, b_fld,
//			     *bi->get_star_story());

//       PRL(T_eff);
//       T_eff *= 1000;

//       M_tot = m_env + m_core;
       sstype = summarize_stellar_type(stype);
       star_class = get_spectral_class(T_eff);
       
//       PRC(L_eff);PRC(T_eff);PRL(M_tot);

       ltm_to_ubvri(log10(L_eff), log10(T_eff), M_tot,
		     U, B, V, R, I);
//       PRC(U);PRC(B);PRC(V);PRC(R);PRC(I);
       
       if (find_qmatch(bi->get_star_story(), "Class"))
	 tpe_class = extract_stellar_spec_summary_string(
             getsq(bi->get_star_story(), "Class"));
       if (L_eff>0)
          R_eff = pow(T_eff/cnsts.parameters(solar_temperature), 2)
	       / sqrt(L_eff);
     }
     else if (bi->get_use_sstar()) {

        put_dyn(bi, cerr);
       	stype = bi->get_starbase()->get_element_type();
	M_tot  = bi->get_starbase()->conv_m_dyn_to_star(bi->get_mass());
        t_cur = bi->get_starbase()->get_current_time();
        t_rel = bi->get_starbase()->get_relative_age();
        T_eff = bi->get_starbase()->temperature();
        L_eff = bi->get_starbase()->get_luminosity();
        star_class = get_spectral_class(T_eff);
	R_eff = bi->get_starbase()->get_effective_radius();
	ltm_to_ubvri(log10(L_eff), log10(T_eff), M_tot,
		     U, B, V, R, I);

//       PRC(L_eff);PRC(T_eff);PRL(M_tot);
//       PRC(U);PRC(B);PRC(V);PRC(R);PRC(I);

     }
     else {
       cout << "    No stellar information found for: ";
       bi->pretty_print_node(cout);
       return;
     }

     int known_at;
     bool is_known = false;
     bool should_be_known = false;
     
     if (remnant(stype) && t_rel <= 1) 
       should_be_known = true;
     
     if (remnant(stype)) {
       for (int i=0; i<NNEBULAE; i++)
	 if (bi->get_index() == nebula[i].id) {
	   is_known = true;
	   known_at = i;
	 }

     }

     if (should_be_known) {
       if (!is_known)
	 for (int i=0; i<NNEBULAE; i++) 
	   if (nebula[i].id < 0) {
	     nebula[i].id  = bi->get_index();
	     nebula[i].pos = pos;
	     nebula[i].time = t_rel;
	     known_at = i;
	   }
       
       if (stype== Helium_Dwarf || Carbon_Dwarf || Oxygen_Dwarf) 
	 print_pl_nebula(nebula[known_at].pos, scale_L * (0.5 + t_rel));
       else {
#if 0
	 print_sn_nebula(nebula[known_at].pos, scale_L * (0.5 + t_rel));
	 if (t_rel<0.01) 
	   cout << "object { single_star \n"
		<< "         finish { ambient <"
		<< 0.4 << ", " << 0.6 << ", " << 0.7 << ">}\n" 
		<< "         scale " << scale_L * pow(100 * (0.01-t_rel), 3)
		<< " translate < "
		<< pos[0] << ", " << pos[1] << ", " << pos[2]-0.5 << " > }"
		<< endl;
#endif
	 
       }
     }
     else if(is_known)
       nebula[known_at].id = -1;

     //  Add SN remnant
     for (int i=0; i<NSN_REMNANT; i++)
       if(bi->get_index() == sn_remnant[i].id &&
	  t_cur >= sn_remnant[i].time && t_cur < sn_remnant[i].time+0.025) {
	 if (sn_remnant[i].bh)
	   print_sn_nebula(pos, scale_L * (0.5 + t_rel));
	 else
	   print_pl_nebula(pos, scale_L * (0.5 + t_rel));
	 if(t_cur >= sn_remnant[i].time && t_cur < sn_remnant[i].time+0.01) {
	   cout << "object { single_star \n"
		<< "         finish { ambient <";
	   if (sn_remnant[i].bh)
	     cout << 0.4 << ", " << 0.6 << ", " << 0.7 << ">}\n";
	   else
	     cout << 0.56 << ", " << 0.14 << ", " << 0.14 << ">}\n";
	   cout << "         scale " << scale_L * pow(100 * (0.01-t_rel), 3)
		<< " translate < "
		<< pos[0] << ", " << pos[1] << ", " << pos[2]-0.5 << " > }"
		<< endl;
	 }
       }

  add_collision_effect(bi, pos, time, scale_L); 
#if 0
     //  Add collisions
cerr << "Add collisions"<<endl;
     for (int i=0; i<NCOLLISIONS; i++) {
       PRC(i);PRC(bi->get_index());PRL(collision[i].starid);
       if(bi->get_index() == collision[i].starid && 
	  t_cur >= collision[i].time && t_cur < collision[i].time+0.125) {
	 cout << "object { OStar scale "
	      << scale_L * pow(5 * (0.125 - (t_cur-collision[i].time)), 3)
	      << " translate < "
	      << pos[0] << ", " << pos[1] << ", " << pos[2]-0.5 << " > }"
	      << endl;
       }
     }
#endif
     
     // The eye works as an exponential detector.
     real logL_eff = scale_L * sqrt(log10(1 + L_eff));
#if 0       
     cout << "object { "
	  << type_short_string(star_class) << "Star scale "
	  << logL_eff << " translate < "
	  << pos[0] << ", " << pos[1] << ", " << pos[2] << " > }"
	  << endl;

     if (tpe_class == Accreting)
       cout << "object { Accretion_Disc scale "
	    << logL_eff << " translate < "
	    << pos[0] << ", " << pos[1] << ", " << pos[2] << " > }"
	    << endl;
#endif

     real G;
     real apparent_size = -1;
     if (filter == Visual) {

	// Transform UBVRI into RGB on a rather clumsy way.
	real UM100 = -7.38974,  UM1 = 5.17096,  UM01 = 17.9608;
	real BM100 = -6.23550,  BM1 = 5.06279,  BM01 = 16.7472;
	real VM100 = -5.91088,  VM1 = 4.46967,  VM01 = 15.0682;
	real RM100 = -15.9099,  RM1 = 4.13746,  RM01 = 13.7278;
	real IM100 = -25.9089,  IM1 = 3.81839,  IM01 = 12.0133;

//	RM100 = -16; RM01 = 14;
	RM100 = -6; RM01 = 11.4;
	R = Starlab::max(0., Starlab::min(1., (R-RM100)/(RM01-RM100)));
//	R = (R-RM100)/(RM01-RM100);
//	PRC(R);

	VM100 = -5.8; VM01 = 13;
	G = 1 - Starlab::min(1., Starlab::max(0., sqrt(abs((V-VM01)/(VM01-VM100)))));
//	G = 1 - Starlab::min(1., Starlab::max(0., sqrt(abs((V-VM1)/(VM01-VM100)))));
//	G = 1 - sqrt(abs((V-VM1)/(VM01-VM100)));
//	PRC(G);

//	BM100 = -7; BM01 = 17;
	BM100 = -7; BM01 = 15;
	B = pow(Starlab::max(0., Starlab::min(1., (BM01-B)/(BM01-BM100))), 2);
//	B = (BM01-B)/(BM01-BM100);
//	PRL(B);

//	apparent_size = 0.1 * scale_L * (VM01-V)/(VM01-VM100);
	apparent_size = 0.1 * scale_L * logL_eff;
     }
     else {

       real VM100 = -5.91088,  VM1 = 4.46967,  VM01 = 15.0682;
       apparent_size = 0.1 *scale_L * M_tot;

       switch(stype) {
          case Carbon_Star:  
          case Helium_Star:  
          case Helium_Giant: R = 0; B=0; // green
	    G = 1 - Starlab::min(1., Starlab::max(0., sqrt(abs((V-VM1)/(VM01-VM100)))));
	    apparent_size = 0.1 * scale_L * (VM01-V)/(VM01-VM100);
	                      break;
          case Carbon_Dwarf: 
          case Helium_Dwarf: 
          case Oxygen_Dwarf: G=0; B=0; // Red
	    R = 1 - Starlab::min(1., Starlab::max(0., sqrt(abs((V-VM1)/(VM01-VM100)))));
	    apparent_size = 0.1 * scale_L * (VM01-V)/(VM01-VM100);
	                      break;
          case Thorn_Zytkow: R = 1; G=0; B=0;
	    apparent_size = 0.1 * scale_L * (VM01-V)/(VM01-VM100);
	                      break;
          case Xray_Pulsar:  
          case Radio_Pulsar: 
          case Neutron_Star: R = 1; G=0; B=0; // Blue
	                       B = -0.3 * log10(Starlab::min(1., p_rot));
			       apparent_size = 0.1 * scale_L * B;
 	                       break;
          case Black_Hole:   R = 1; G=1; B=1;
	                      break;
          default:
	       R = 0; G=0; B=0;
	    apparent_size = -1;
	                      break;
       
       };
     }
     
  if(apparent_size>0) {
     cout << "object { single_star \n"
          << "         finish { ambient <"
          << R << ", " << G << ", " << B << ">}\n" 
	  << "         scale " << apparent_size << " translate < "
	  << pos[0] << ", " << pos[1] << ", " << pos[2] << " > }"
	  << endl;
   }

     if (tpe_class == Accreting)
       cout << "object { Accretion_Disc scale "
	    << 0.05*logL_eff << " translate < "
	    << pos[0] << ", " << pos[1] << ", " << pos[2] << " > }"
	    << endl;

}

local void print_node(dyn *bi, vec pos, real mass_scale,
		      real mmax)  {

  real time = getrq(bi->get_root()->get_dyn_story(), "real_system_time");

  real apparent_size = mass_scale * sqrt(bi->get_mass());
  real mmin = 0;
  real m = bi->get_mass();
  real R = sqrt((mmax-m)/(mmax-mmin));
  real G = 0.5;
  real B = sqrt(m/(mmax-mmin));
  if(m>0.0125) { // Red
    R = 1;
    G = 0;
    B = 0;
  }
  else if(m<0.003125) {
    R = 0.8;
    G = 0.8;
    B = 0.8;
  }
  else { // Yellow
    R = 0.6;
    G = 0.8;
    B = 0.196078;
//    R = 2.50;
//    G = 2.00;
//    B = 0.00;
  }
#if 0
  vec cam_pos;
    cam_pos[0] = 1;
    cam_pos[1] = 0;
    cam_pos[2] = 0;
  real V=100;
  real v = (V-pos[3])*tan(abs(cam_pos-pos)/pos[3]);
  R = 1;G=0;B=0;
  cout << "object { simple_star \n"
       << "         finish { ambient <"
       << R << ", " << G << ", " << B << ">}\n";
  cout << "         scale " << apparent_size; 
  cout << " translate < "
       << pos[0]-v << ", " << pos[1] << ", " << pos[2] << " > }"
       << endl;
  R=0;G=0;B=1;
  cout << "object { simple_star \n"
       << "         finish { ambient <"
       << R << ", " << G << ", " << B << ">}\n";
  cout << "         scale " << apparent_size; 
  cout << " translate < "
       << pos[0]+v << ", " << pos[1] << ", " << pos[2] << " > }"
       << endl;
#endif


  R=G=B=1;
  cout << "object { simple_star \n"
       << "         finish { ambient <"
       << R << ", " << G << ", " << B << ">}\n";
  if(time>=25) 
    cout << "         scale " << 1./(1+0.2*(25-20))  *apparent_size; 
  else if(time>=20) 
    cout << "         scale " << 1./(1+0.2*(time-20))  *apparent_size; 
  else
    cout << "         scale " << apparent_size; 
  cout << " translate < "
       << pos[0] << ", " << pos[1] << ", " << pos[2] << " > }"
       << endl;

//  pos[2] -= apparent_size;
  add_collision_effect(bi, pos, time, mass_scale); 
//  pos[2] += apparent_size;

#if 0
     //  Add collisions
     for (int i=0; i<NCOLLISIONS; i++) {
//       if(bi->get_index() == collision[i].starid && 
       if(bi->get_name()) {
//       PRC(bi->get_index());PRC(bi->get_name());PRL(collision[i].name);
       if(!strcmp(bi->get_name(), collision[i].name) && 
	  time >= collision[i].time && time < collision[i].time+0.125) {
	 cerr << "Adding collision at time = " << time << endl;
	 cerr << "object { OStar scale "
	      << 2 *mass_scale * pow(5 * (0.125 - (time-collision[i].time)), 3)
	      << " translate < "
	      << pos[0] << ", " << pos[1] << ", " << pos[2]-apparent_size << " > }"
	      << endl;
	 cout << "object { OStar scale "
	      << 2 *mass_scale * pow(5 * (0.125 - (time-collision[i].time)), 3)
	      << " translate < "
	      << pos[0] << ", " << pos[1] << ", " << pos[2]-apparent_size << " > }"
	      << endl;
       }
     }
     }
#endif

//  cout << "object { single_star \n"
//       << "         finish { ambient <"
//       << R << ", " << G << ", " << B << ">}\n" 
//       << "         scale " << apparent_size << " translate < "
//       << pos[0] << ", " << pos[1] << ", " << pos[2] << " > }"
//       << endl;
}

local int print_povray_binary_recursive(dyn *b,
					vec dc_pos, 
					real mass_limit, real number_limit,
					bool povray, real scale_L,
					filter_type filter) {

  int nb = 0;
  if (b->get_oldest_daughter()) {
    
    vec r_com  = b->get_pos() - dc_pos;

    real m_tot = b->get_starbase()->conv_m_dyn_to_star(b->get_mass());

    for_all_daughters(dyn, b, bb)
      if (bb->n_leaves() >= 2) 
	nb += print_povray_binary_recursive(bb, dc_pos,
					    mass_limit, number_limit,
					    povray, scale_L, filter);
      else {
	print_star(bb, bb->get_pos()-r_com, scale_L, filter);
      }

  }
  return nb;
}

local void print_povtime(real time,
			 vec pos,
			 real scale=1,
			 real depth=0.25) {

    int p = cout.precision(LOW_PRECISION);
    cout << "text { ttf \"timrom.ttf\" ";
    if(time==0)
      cout << "\"Time = " << time << " \" ";
    else if(time>0)
      cout << "\"Time = " << time/73.387 << " Myr \" ";
    else
      cout << "\"Time = " << -time << " N-body \" ";

    cout << depth << ", 0\n"
	 << "       pigment { Red }\n"
         << "       translate < " << pos[0] << ", "
	                  << pos[1] << ", "
	                  << pos[2] << " >\n"
         << "       scale " << scale
	 << "}\n" << endl;
    cout.precision(p);
}

local void print_start_text(vec cam_pos) {
}

local void print_povray_stars(dyn *b, real mass_limit,
			      real number_limit,
			      bool povray,
			      real scale_L,
			      filter_type filter) {

cerr << "Print STARS"<<endl;
  for (int i=0; i<NNEBULAE; i++) 
    nebula[i].id  = -1;
    
  bool cod = false;

  vec dc_pos = 0;
  bool try_com = false;
  if(abs(dc_pos) == 0) {
    if (find_qmatch(b->get_dyn_story(), "density_center_pos")) {
      
      if (getrq(b->get_dyn_story(), "density_center_time")
	  != b->get_real_system_time()) {
	warning("mkpovfile: neglecting out-of-date density center");
	try_com = true;
      } else
	cod = true;
      
      dc_pos = getvq(b->get_dyn_story(), "density_center_pos");
    }

    if (try_com && find_qmatch(b->get_dyn_story(), "com_pos")) {

      if (getrq(b->get_dyn_story(), "com_time")
	  != b->get_real_system_time()) {
	warning("lagrad: neglecting out-of-date center of mass");
      } else
	dc_pos = getvq(b->get_dyn_story(), "com_pos");
    }
  }

  // For now: put denxity center in geometric origin
  dc_pos = 0;
  
  int ns=0, nb=0;
  for_all_daughters(dyn, b, bi) 
    if (bi->is_leaf()                    &&
	((bi->get_index()>0 && bi->get_index() <= number_limit) ||
	 bi->get_starbase()->get_total_mass() >= mass_limit)) {
	
      print_star(bi, bi->get_pos() - dc_pos, scale_L, filter);
      ns++;
    }

  for_all_daughters(dyn, b, bi) {
    if (!bi->is_leaf())                   
      nb += print_povray_binary_recursive(bi, dc_pos,
					  mass_limit, number_limit,
					  povray, scale_L, filter);
  }
  
if (!povray && ns+nb==0) 
    cout << "            (none)\n";
}

local void print_povray_bodies(dyn *b, real mass_limit,
			      real number_limit, real mmax, 
			      real scale_M) {


/*
  real mmin = VERY_LARGE_NUMBER, mmax = -1;
  for_all_leaves(dyn, b, bi) {
    real m = bi->get_mass();
    mmin = min(mmin, m);
    mmax = max(mmax, m);
  } 
*/
    
  int ns=0;
  for_all_leaves(dyn, b, bi) 
    if(bi->get_mass() >= mass_limit) {
	
      print_node(bi, bi->get_pos(), scale_M, mmax);
      ns++;
    }
  
  if (ns==0) 
    cout << "//            (none)\n";

}

void print_povray_header(dyn* b, vec cam_pos,
                         int camera_on_star_id, 
			 real aperture, real gamma, 
			 int blur_samples,
                         int counter, real scale_L,
			 int horizontal, int vertical,
			 bool print_hrd) {

  cout << "\n\n// End of POVRay scenary." << endl;
  cout << "\n\n// Start new POVRay scenary." << endl;

  cout << "//       filename = starlab_";
  print_filename_counter(counter, cout);
  cout << ".pov\n";
    
  cout << "#include \"colors.inc\"" << endl
       << "#include \"shapes.inc\"" << endl
       << "#include \"textures.inc\"" << endl
       << "#include \"glass.inc\"" << endl;
  cout << "#include \"astro.inc\"" << endl;

  cerr<< endl;

  cout << "#version 3.0" << endl << endl;

  cout << "global_settings { " << endl
       << "  assumed_gamma " << gamma << endl 
       << "  max_trace_level 15" << endl
       << "  ambient_light White" << endl
       << "}"
       << endl << endl;

  if (camera_on_star_id<=0) {

    vec normal;
    normal[1] = 1;

//    cout << "// Normal to camera " << endl;
//    cout << "   #declare normal_to_camera = < " << normal[0] << ", "
//	                                        << normal[1] << ", "
//	                                        << normal[2] << " >" << endl;
    cout  << "#declare camera_location  = <0, 0, -1>" << endl;
    cout << "#declare normal_to_camera = -z" << endl;
    cout << "#declare camera_up        = y" << endl;
    
    cout << "camera {" << endl
	 << "   location < " << cam_pos[0] << ", "
                             << cam_pos[1] << ", "
                             << cam_pos[2] << " >" << endl
	 << "   look_at  <0.0, 0.0, " << cam_pos[2]+5 << ">" << endl
	 << "   blur_samples " << blur_samples << endl; 
    if (aperture>0)
      cout << "   focal_point <0, 0, " << cam_pos[2]+5 << ">" << endl
	   << "   aperture " << aperture << endl;
    cout << "}" << endl << endl;

    //
    cout << "light_source { < " << cam_pos[0] << ", "
                                << cam_pos[1] << ", "
                                << cam_pos[2] << " > White }" << endl;

//    cout << "#object { Initial_text }\n" << endl;
    //
    
//    real time = b->get_starbase()->conv_t_dyn_to_star(b->get_real_system_time());
//    real time = b->get_real_system_time();
    real time = -1*getrq(b->get_dyn_story(), "real_system_time");
//    real time = -(20 + 0.015625*(counter-160));
    vec time_pos;
    if(vertical<=400) {         //320x240
      time_pos[0] = -10;
      time_pos[1] = 6;
    }
    else {        //640x480
      time_pos[0] = -16;
      time_pos[1] = -11;
    }
//    time_pos[0] = -13;
//    time_pos[1] = -9.5;
//    time_pos[2] = 1;
    time_pos[0] = -7;
    time_pos[1] = -5;
    time_pos[2] = 1;
//
//    real scale = 0.2;
//    real depth = 0.25;
    real scale = 0.1;
    real depth = 0.125;
    print_povtime(time, time_pos, scale, depth);    
      
    if (print_hrd)
      print_hertzsprung_Russell_diagram(b, cam_pos);
    
  }
  else {
    cam_pos = 0;
    if(!print_camera_position_recursive(b, cam_pos, camera_on_star_id, scale_L,
					aperture, blur_samples)) {
      cout << "No star id = " << camera_on_star_id << " found." << endl;
      cout << " in print_camera_position_recursive() " << endl;
      exit(-1);
    }
  }

  cout << "// Here follow the definitions of the stars..." << endl
       << endl;
}

void print_mpegplayer_param_file(ostream& s,
				 int first_frame /* = 1   */,
				 int last_frame  /* = 100 */,
				 int horizontal  /* = 120 */,
				 int vertical    /* = 90  */,
				 int GOP_size    /* = 10  */
				) {
				 

     s <<"# mpeg_encode parameter file\n"
	  << "PATTERN         IBBBPBBBBP\n"
	  << "OUTPUT          starlab.mpg\n"
	  << endl;

     s << "YUV_SIZE        "<<vertical<<"x"<<horizontal<<"\n"
	  << "BASE_FILE_FORMAT        PPM\n"
	  << "INPUT_CONVERT   *\n"
	  << "GOP_SIZE        "<<GOP_size<<"\n" 
	  << "SLICES_PER_FRAME  1\n"
	  << endl;

     s << "INPUT_DIR       .\n"
	  << "INPUT\n"
	  << "starlab_*.ppm   [";
          print_filename_counter(first_frame, s);
      s   << "-";
          print_filename_counter(last_frame-1, s);
      s   << "]\n"
	  << "END_INPUT\n" << endl;

    s << "FORCE_ENCODE_LAST_FRAME\n"
//	 << "PIXEL           HALF\n"
	 << "PIXEL           WHOLE\n"
	 << "RANGE           10\n"
	 << endl;

    s << "PSEARCH_ALG     LOGARITHMIC\n"
	 << "BSEARCH_ALG     CROSS2\n"
	 << endl;

    s << "IQSCALE         8\n"
	 << "PQSCALE         10\n"
	 << "BQSCALE         25\n"
	 << endl;

    s << "REFERENCE_FRAME ORIGINAL" << endl;
}

void rdc_and_wrt_movie(dyn *b, bool povray, real scale_L, real mmax, 
		       char fltr) {

    strcpy(collision[0].name, "255a+255b"); 
    collision[0].starid = 0; 
    collision[0].time = 2.37956;
    strcpy(collision[1].name, "576a+576b"); 
    collision[1].starid = 0;
    collision[1].time = 3.79467;
    strcpy(collision[2].name, "463a+463b"); 
    collision[2].starid = 0;
    collision[2].time = 8.01511;
    strcpy(collision[3].name, "1747a+1747b"); 
  collision[3].starid = 0;
  collision[3].time = 8.97172;
  strcpy(collision[4].name, "1a+218b"); 
  collision[4].starid = 1;
  collision[4].time = 9.21207;
  strcpy(collision[5].name, "1b<+2>"); 
  collision[5].starid = 1;
  collision[5].time = 14.0345;
  strcpy(collision[6].name, "1a<+2>"); 
  collision[6].starid = 1;
  collision[6].time = 14.6999;
  strcpy(collision[7].name, "6a+6b"); 
  collision[7].starid = 0;
  collision[7].time = 15.9929;
  strcpy(collision[8].name, "1b<+3>"); 
  collision[8].starid = 1;
  collision[8].time = 22.9681;
  strcpy(collision[9].name, "14a+14b"); 
  collision[9].starid = 26;
  collision[9].time = 34.916;

#if 0
  // Fill in the collision database.
    strcpy(collision[0].name, "1+167"); 
    collision[0].starid = 167; 
    collision[0].time = 20.3479;
    strcpy(collision[1].name, "1<+2>"); 
    collision[1].starid = 69;
    collision[1].time = 21.1829;
    strcpy(collision[2].name, "30+64"); 
    collision[2].starid = 30;
    collision[2].time = 21.7115;
    strcpy(collision[3].name, "1<+3>"); 
  collision[3].starid = 4;
  collision[3].time = 23.6172;
  strcpy(collision[4].name, "1<+4>"); 
  collision[4].starid = 1435;
  collision[4].time = 25.2399;
  strcpy(collision[5].name, "155+938"); 
  collision[5].starid = 155;
  collision[5].time = 25.8208;
  strcpy(collision[6].name, "23+34"); 
  collision[6].starid = 23;
  collision[6].time = 26.3169;
  strcpy(collision[7].name, "6+403"); 
  collision[7].starid = 6;
  collision[7].time = 26.4171;
  strcpy(collision[8].name, "102+483"); 
  collision[8].starid = 102;
  collision[8].time = 26.901;
  strcpy(collision[9].name, "26+1007"); 
  collision[9].starid = 26;
  collision[9].time = 27.817;
  strcpy(collision[10].name, "29+173"); 
  collision[10].starid = 29;
  collision[10].time = 28.2324;
  strcpy(collision[11].name, "1<+6>"); 
  collision[11].starid = 2;
  collision[11].time = 28.234;
  strcpy(collision[12].name, "1<+7>"); 
  collision[12].starid = 1917;
  collision[12].time = 28.7764;
  strcpy(collision[13].name, "1<+8>");  
  collision[13].starid = 46;
  collision[13].time = 29.3036;
  strcpy(collision[14].name, "1<+9>");  
  collision[14].starid = 1562;
  collision[14].time = 29.6476;
  strcpy(collision[15].name, "1<+10>");  
  collision[15].starid = 0;
  collision[15].time = 30.3179;
  strcpy(collision[16].name, "86+123");  
  collision[16].starid = 0;
  collision[16].time = 31.8339;
  strcpy(collision[17].name, "1<+11>");  
  collision[17].starid = 0;
  collision[17].time = 31.9796;
  strcpy(collision[18].name, "2+149");
  collision[18].starid = 0;
  collision[18].time = 32.3984;
  strcpy(collision[19].name, "1<+12>");
  collision[19].starid = 0;
  collision[19].time = 32.5874;
  strcpy(collision[20].name, "1<+13>");
  collision[20].starid = 0;
  collision[20].time = 32.7505;
  strcpy(collision[21].name, "15+19");
  collision[21].starid = 0;
  collision[21].time = 32.9756;
  strcpy(collision[22].name, "1<+14>");
  collision[22].starid = 0;
  collision[22].time = 33.4833;
  strcpy(collision[23].name, "1<+15>");
  collision[23].starid = 0;
  collision[23].time = 33.5243;
  strcpy(collision[24].name, "1<+16>");
  collision[24].starid = 0;
  collision[24].time = 33.7629;
  strcpy(collision[25].name, "1<+17>");
  collision[25].starid = 0;
  collision[25].time = 34.4345;
  strcpy(collision[26].name, "1<+18>");
  collision[26].starid = 0;
  collision[26].time = 34.4508;
  strcpy(collision[27].name, "1<+19>");
  collision[27].starid = 0;
  collision[27].time = 34.6628;
  strcpy(collision[28].name, "1<+20>");
  collision[28].starid = 0;
  collision[28].time = 34.7393;
#endif








#if 0
  sn_remnant[0].bh = true;
  sn_remnant[0].id = 4441;
  sn_remnant[0].time = 4.42;
  sn_remnant[1].bh = true;
  sn_remnant[1].id = 1;
  sn_remnant[1].time = 5.35;
  sn_remnant[1].bh = false;
  sn_remnant[2].id = 2056;
  sn_remnant[2].time = 6.14;
  sn_remnant[1].bh = false;
  sn_remnant[3].id = 3193;
  sn_remnant[3].time = 6.26;
#endif

    filter_type filter = get_filter_type(fltr);
    
    real mass_limit = 0;
    real number_limit = VERY_LARGE_NUMBER;
    if(b->get_use_sstar())
      print_povray_stars(b, mass_limit, number_limit, povray, scale_L, filter);
    else
      print_povray_bodies(b, mass_limit, number_limit, mmax, scale_L);
    
}

#else

//-----------------------------------------------------------------------------
//  main  --  driver to use  compute_mass_radii() as a tool
//-----------------------------------------------------------------------------
enum SF_base_type {No_base=0, StTr_station, SW_ISD}; 

main(int argc, char ** argv)
{
    bool  c_flag = false;      // if TRUE: a comment given on command line
    bool  v_flag = false;      // if TRUE: a comment given on command line
    bool  B_flag = false;  
    bool  print_hrd = false;
    bool  Stellar = true;

    char filter = 'V';
    real aperture    = -1;
    vec first_campos, last_campos, cam_pos;
    first_campos[2] = last_campos[2] = -10;
    real theta_rotate = 0;
    real phi_rotate   = 0;
    int blur_samples = 10;
    int camera_on_star_id = -1;
    int nstart = 0;

    char * mpeg_paramfile = "paramfile.mpeg";
    int horizontal =  90;
    int vertical   = 120;
    int first_frame = 1;
    int last_frame = 0;
    real scale_L   = 1;

    real gamma = 2.2;

    SF_base_type base_type = No_base;
    int nsteps = 32767;
    int nskip = 0; 
    int nint_skip = 0;
    
    char  *comment;
    check_help();

    extern char *poptarg;
    int c;
    char* param_string = "a:Bb:F:f:g:H:hI:L:N:n:W:X:x:Y:y:Z:z:P:T:S:s:c:";

    while ((c = pgetopt(argc, argv, param_string)) != -1)
	switch(c)
	    {
	    case 'a': aperture = atof(poptarg);
	              break;
//	    case 'B': base_type = (SF_base_type)atoi(poptarg);
//	              break;
	    case 'B': Stellar = !Stellar;
	              break;
	    case 'b': blur_samples = atoi(poptarg);
	              break;
	    case 'F': filter = *poptarg;
	              break;
	    case 'f': mpeg_paramfile = poptarg;
	              break;
	    case 'g': gamma = atof(poptarg);
	              break;
	    case 'H': horizontal = atoi(poptarg);
	              break;
	    case 'h': print_hrd = true;
	              break;
	    case 'I': camera_on_star_id = atoi(poptarg);
	              break;
	    case 'W': vertical = atoi(poptarg);
	              break;
	    case 'L': scale_L *= atof(poptarg);
	              break;
	    case 'P': phi_rotate = atof(poptarg);
	              phi_rotate *= cnsts.mathematics(pi)/180;
	              break;
	    case 'T': theta_rotate = atof(poptarg);
	              theta_rotate *= cnsts.mathematics(pi)/180;
	              break;
	    case 'N': nsteps = atoi(poptarg);
	              break;
	    case 'n': nstart = atoi(poptarg);
	              break;
	    case 'X': first_campos[0] = atof(poptarg);
	              break;
	    case 'Y': first_campos[1] = atof(poptarg);
	              break;
	    case 'Z': first_campos[2] = atof(poptarg);
	              break;
	    case 'x': last_campos[0] = atof(poptarg);
	              break;
	    case 'y': last_campos[1] = atof(poptarg);
	              break;
	    case 'z': last_campos[2] = atof(poptarg);
	              break;
	    case 'S': nskip = atoi(poptarg);
	              break;
	    case 's': nint_skip = atoi(poptarg);
	              break;
	    case 'c': c_flag = true;
		      comment = poptarg;
		      break;
	    case 'v': v_flag = true;
		      break;
            case '?': params_to_usage(cout, argv[0], param_string);
		      get_help();
		      exit(1);
	    }            

    dyn *b;

    int GOP_size = Starlab::min(1000, nsteps);
    remove(mpeg_paramfile);
    ofstream os(mpeg_paramfile, ios::app|ios::out);
    if (!os) cerr << "\nerror: couldn't create file "
		  << mpeg_paramfile << endl;

    print_mpegplayer_param_file(os, first_frame, GOP_size,
				horizontal, vertical, GOP_size);
    os.close();

  int id;
  real time, mmax = -1;
  vec pos;
    
    int nread = nstart;
    int j=0;
    for (int i = 0; i < nskip; i++) {
	cerr << " skipping snapshot " << i << endl;
	if (!forget_node()) exit(1);
	if (j == nint_skip) {
	    nread++;
	    j=0;
	}
	else
	    j++;
    }
    
    int nsnap = 0;	
    while (b = get_dyn()) {

      if(nsnap==0)
	for_all_leaves(dyn, b, bi) {
	  mmax = Starlab::max(mmax, bi->get_mass());
	}

    if(Stellar) {
       b->set_use_sstar(true);
       real T_start = 0;
       if(b->get_starbase()->get_stellar_evolution_scaling()) {
	 //addstar(b,                             // Note that T_start and
//		 T_start,                       // Main_Sequence are
//		 Main_Sequence,                 // defaults. They are
//		 true);                         // ignored if a star
       }
       else {
	 cerr << "No stellar evolution scaling present"<<endl;
	 cerr << "Continue with node prescription"<<endl;
       }
//       b->get_starbase()->print_stellar_evolution_scaling(cerr);
     }
      //b->set_use_sstar(true); 
      // Do not make new stars but make sure the starbase is searched.
      //b->set_use_sstar(false);

	nsnap ++;
	nread ++;
	cerr << " reading snapshot " << nread << ", and "
	     << nsnap << " done." << endl;
	
	if (camera_on_star_id<=0) 
	  new_camera_position(cam_pos, first_campos, last_campos,
			      nsteps, nsnap);
	
	print_povray_header(b, cam_pos, camera_on_star_id,
			    aperture, gamma, blur_samples, nread, scale_L,
			    horizontal, vertical, print_hrd);

       if (base_type == StTr_station) {
	  // Plot the starbase at the origin.
	  cout << "#include \"starbase.inc\"" << endl;
	  cout << "object { Base_Station }\n" << endl;
	}
       else if (base_type == SW_ISD) {
	 cout << "object { StarWars_Scene }\n" << endl;
       }
	
	rdc_and_wrt_movie(b, true, scale_L, mmax, filter);
	last_frame++;

	if (v_flag)
	  put_dyn(b);
	
	rmtree(b);

	for (int i = 0; i < nint_skip; i++){
	  nsnap ++;
	  cerr << " skipping snapshot " << nsnap << endl;
	  if (!forget_node()) exit(1);
	}

	if (nsnap==nsteps)
	  break;
	
    }

}

#endif

// endof: mkpovfile.C

