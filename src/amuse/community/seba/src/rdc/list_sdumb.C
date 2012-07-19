
//// list_snap:  Print times of all snapshots in the input stream.
////
//// Options:    none


#include "single_star.h"
#include "sstar_to_dyn.h"
#include "hdyn.h"

#ifndef TOOLBOX

#if 0
local void combine_ubvri(real Up, real Bp, real Vp, real Rp, real Ip,
                         real Us, real Bs, real Vs, real Rs, real Is,
                         real &U, real &B, real &V, real &R, real &I) {

    real ln10g25 = 2.3025851/2.5;

    U = -2.5 * log10(exp(-ln10g25*Up) + exp(-ln10g25*Us));
    B = -2.5 * log10(exp(-ln10g25*Bp) + exp(-ln10g25*Bs));
    V = -2.5 * log10(exp(-ln10g25*Vp) + exp(-ln10g25*Vs));
    R = -2.5 * log10(exp(-ln10g25*Rp) + exp(-ln10g25*Rs));
    I = -2.5 * log10(exp(-ln10g25*Ip) + exp(-ln10g25*Is));

}

#endif

local bool star_is_bound(hdyn* b, vec& cod_vel) {

  bool bound = false;
  if (0.5*square(b->get_vel()-cod_vel) + b->get_pot() < 0) 
    bound = true;

    return bound;

}

local void get_UBVRI_of_star(hdyn *bi, vec pos,
	real &Us, real &Bs, real &Vs, real &Rs, real &Is) {


  // To solar radii
  //  vec pos = bi->get_pos() - dc_pos;
  pos[0] = bi->get_starbase()->conv_r_dyn_to_star(pos[0]);
  pos[1] = bi->get_starbase()->conv_r_dyn_to_star(pos[1]);
  pos[2] = bi->get_starbase()->conv_r_dyn_to_star(pos[2]);

  real time = bi->get_starbase()->conv_t_dyn_to_star(bi->get_system_time());

  // And now to parsec
  real Rsun_per_parsec = cnsts.parameters(solar_radius)
                       / cnsts.parameters(parsec);
  pos[0] *= Rsun_per_parsec;
  pos[1] *= Rsun_per_parsec;
  pos[2] *= Rsun_per_parsec;
  

     star_type_spec tpe_class = NAC;
     spectral_class star_class;
     stellar_type stype = NAS;
     stellar_type_summary sstype = ZAMS;
     real t_cur, m_rel, m_env, m_core, mco_core, T_eff, L_eff, p_rot, b_fld;
     real t_rel=0, R_eff=0;
     real M_tot;
     if (bi->get_use_sstar()) {
       	stype = bi->get_starbase()->get_element_type();
	M_tot  = bi->get_starbase()->conv_m_dyn_to_star(bi->get_mass());
        t_cur = bi->get_starbase()->get_current_time();
        t_rel = bi->get_starbase()->get_relative_age();
        T_eff = bi->get_starbase()->temperature();
        L_eff = bi->get_starbase()->get_luminosity();
        star_class = get_spectral_class(T_eff);
	R_eff = bi->get_starbase()->get_effective_radius();
	ltm_to_ubvri(log10(L_eff), log10(T_eff), M_tot,
		     Us, Bs, Vs, Rs, Is);

     }
     else if (bi->get_star_story()) {

       extract_story_chapter(stype, t_cur, t_rel, 
			     m_rel, m_env, m_core, mco_core,
			     T_eff, L_eff, p_rot, b_fld,
			     *bi->get_star_story());

       M_tot = m_env + m_core;
       sstype = summarize_stellar_type(stype);
       star_class = get_spectral_class(T_eff);
       
       ltm_to_ubvri(log10(L_eff), log10(T_eff), M_tot,
		     Us, Bs, Vs, Rs, Is);
       
       if (find_qmatch(bi->get_star_story(), "Class"))
	 tpe_class = extract_stellar_spec_summary_string(
             getsq(bi->get_star_story(), "Class"));
       if (L_eff>0)
          R_eff = 
	    * sqrt(L_eff)/pow(T_eff/cnsts.parameters(solar_temperature), 2);
     }
     else {
       cout << "    No stellar information found: " << endl;
       return;
     }
}

local void print_star(hdyn *bi, bool bound, vec pos, vec vel,
		      real &Up, real &Bp, real &Vp, real &Rp, real &Ip, 
		      bool verbose) { 


  int id = bi->get_index();

  //  vec cod_vel = 0;
  //  bool bound = star_is_bound(bi, cod_vel);

  // To solar radii
  //  vec pos = bi->get_pos() - dc_pos;
  pos[0] = bi->get_starbase()->conv_r_dyn_to_star(pos[0]);
  pos[1] = bi->get_starbase()->conv_r_dyn_to_star(pos[1]);
  pos[2] = bi->get_starbase()->conv_r_dyn_to_star(pos[2]);

  real time = bi->get_starbase()->conv_t_dyn_to_star(bi->get_system_time());

  // And now to parsec
  real Rsun_per_parsec = cnsts.parameters(solar_radius)
                       / cnsts.parameters(parsec);
  pos[0] *= Rsun_per_parsec;
  pos[1] *= Rsun_per_parsec;
  pos[2] *= Rsun_per_parsec;

//  real to_Rsun_Myr = cnsts.physics(km_per_s) * cnsts.physics(Myear)
//               / cnsts.parameters(solar_radius);
                       
//      real to_Rsun_Myr = cnsts.physics(km_per_s) * cnsts.physics(Myear)
//	               / cnsts.parameters(solar_radius);
//                       
//      real to_dyn      = bi->get_starbase()->conv_r_star_to_dyn(1)
//                       / bi->get_starbase()->conv_t_star_to_dyn(1);
//      vel = vel/(to_Rsun_Myr * to_dyn);
                       
  real to_km_per_second = cnsts.parameters(solar_radius)
                        / (cnsts.physics(km_per_s) * cnsts.physics(Myear));
  real to_star      = bi->get_starbase()->conv_r_dyn_to_star(1)/
	              bi->get_starbase()->conv_t_dyn_to_star(1);
  to_km_per_second = to_star*to_km_per_second;

//	PRC(to_star);PRL(to_km_per_second);
  vel[0] *= to_km_per_second;
  vel[1] *= to_km_per_second;
  vel[2] *= to_km_per_second;

     star_type_spec tpe_class = NAC;
     spectral_class star_class;
     stellar_type stype = NAS;
     stellar_type_summary sstype = ZAMS;
     real t_cur, m_rel, m_env, m_core, mco_core, T_eff, L_eff, p_rot, b_fld;
     real t_rel=0, R_eff=0;
     real M_tot, Us, Bs, Vs, Rs, Is;	
     if (bi->get_use_sstar()) {
       	stype = bi->get_starbase()->get_element_type();
	M_tot  = bi->get_starbase()->conv_m_dyn_to_star(bi->get_mass());
        t_cur = bi->get_starbase()->get_current_time();
        t_rel = bi->get_starbase()->get_relative_age();
        T_eff = bi->get_starbase()->temperature();
        L_eff = bi->get_starbase()->get_luminosity();
        star_class = get_spectral_class(T_eff);
	R_eff = bi->get_starbase()->get_effective_radius();
	ltm_to_ubvri(log10(L_eff), log10(T_eff), M_tot,
		     Us, Bs, Vs, Rs, Is);

     }
     else if (bi->get_star_story()) {

       extract_story_chapter(stype, t_cur, t_rel, 
			     m_rel, m_env, m_core, mco_core,
			     T_eff, L_eff, p_rot, b_fld,
			     *bi->get_star_story());

       M_tot = m_env + m_core;
       sstype = summarize_stellar_type(stype);
       star_class = get_spectral_class(T_eff);
       
       ltm_to_ubvri(log10(L_eff), log10(T_eff), M_tot,
		     Us, Bs, Vs, Rs, Is);
       
       if (find_qmatch(bi->get_star_story(), "Class"))
	 tpe_class = extract_stellar_spec_summary_string(
             getsq(bi->get_star_story(), "Class"));
       if (L_eff>0)
          R_eff = pow(T_eff/cnsts.parameters(solar_temperature), 2)
	        * sqrt(L_eff);
     }
     else {
       cout << "    No stellar information found for: ";
       bi->pretty_print_node(cout);
       return;
     }

     real U, B, V, R, I;
     combine_ubvri(Up, Bp, Vp, Rp, Ip,
 	           Us, Bs, Vs, Rs, Is,
                   U, B, V, R, I);

     if(verbose)
       cout << " Time= " << time << " id= " << id << " b= " << bound 
	    << " type= " << stype << " m= " << M_tot << " R= " << R_eff
	    << " L= " << L_eff 
	    << " T_eff= " << T_eff 
	    << " r= "  << pos[0] << " " << pos[1] << " " << pos[2] 
	    << " v= "  << vel[0] << " " << vel[1] << " " << vel[2] 
	    << " ubvri= "<<  U << " " << B << " " << V << " " << R << " "
	    << I	<< " :: ";
     else {
       cout << time <<" "<< id << bound <<" "
	    <<" "<< stype <<" "
	    << M_tot <<" "<< R_eff <<" "<< L_eff  <<" "<< T_eff 
	    <<" "<< pos[0] << " " << pos[1] << " " << pos[2] 
	    <<" "<< vel[0] << " " << vel[1] << " " << vel[2] 
	    <<" "<< U << " " << B << " " << V << " " << R << " "
	    << I	<< " ";
       //       PRC(id);PRC(pos[0]);PRC(pos[1]);PRL(pos[2]);
     }

     Up=U;
     Bp=B;
     Vp=V;
     Ip=I;
}

local int print_binary_recursive(hdyn *b, vec dc_pos, vec dc_vel,
	real &U, real &B, real &V, real &R, real &I, bool verbose) { 

  int nb = 0;
  if (b->get_oldest_daughter()) {
    
    // Not relative to COM SPZ@27 Aug 2004
    //    vec r_com  = b->get_pos();
    //    vec v_com  = b->get_vel();
    vec r_com  = b->get_pos() - dc_pos;
    vec v_com  = b->get_vel() - dc_vel;
//    vec r_com  = dyn_something_relative_to_root(b, &dyn::get_pos) -dc_pos;
//    vec v_com  = dyn_something_relative_to_root(b, &dyn::get_vel) -dc_vel;

    real m_tot = b->get_starbase()->conv_m_dyn_to_star(b->get_mass());
    //    PRC(m_tot);PRC(r_com);PRL(v_com);

    for_all_daughters(hdyn, b, bb)
      if (bb->n_leaves() >= 2) {
	if(verbose) cout << "\nBinary: ";
	else cout << "\n2 ";
	U=B=V=R=I=VERY_LARGE_NUMBER;
	nb += print_binary_recursive(bb, dc_pos, dc_vel, U, B, V, R, I,
				     verbose);
      }
      else {
	if (bb->get_parent()==bb->get_root()) {
	  U=B=V=R=I=VERY_LARGE_NUMBER;
	  if(verbose) cout << "\nStar:: ";
	  else cout << "\n1 ";
	}

	// Specific for star clusters which are not at the COM.
	// Such as Arches star cluster models.
	//        print_star(bb, bb->get_pos()-r_com, bb->get_vel()-v_com,
	//           U, B, V, R, I, verbose);
	vec nul = 0;
	bool bound = star_is_bound(bb, nul);
	// includ correction for density center position.
        print_star(bb, bound, bb->get_pos()-dc_pos, bb->get_vel()-dc_vel,
	           U, B, V, R, I, verbose);
//	if (bb->get_parent()==bb->get_root())
//	   cerr << endl;
      }

  }
  return nb;
}

local void print_integrated_cluster(hdyn *b, vec dc_pos,
	                            real r_min, real r_max,
	                            real v_min, real v_max,
	real &U, real &B, real &V, real &R, real &I,
	int project) { 


  int nb = 0;
  real rcom;
  vec r;
  real m_tot = 0;
  if (b->get_oldest_daughter()) {

    vec r_com  = b->get_pos() - dc_pos;

    real Uc, Bc, Vc, Rc, Ic;
    Uc=Bc=Vc=Rc=Ic=VERY_LARGE_NUMBER;
    real Us, Bs, Vs, Rs, Is;
    real U, B, V, R, I;
    U=B=V=R=I=VERY_LARGE_NUMBER;
    for_all_leaves(hdyn, b, bb) {
	r = bb->get_pos()-dc_pos;
	if(project<=2)
          rcom = r[project];
	else
	  rcom = abs(r);
      if (rcom>r_min && rcom<r_max) {

       m_tot += bb->get_mass();
       Us=Bs=Vs=Rs=Is=VERY_LARGE_NUMBER;
       get_UBVRI_of_star(bb, dc_pos, Us, Bs, Vs, Rs, Is);
       if (Vs<v_min && Vs>v_max) {
       combine_ubvri(Uc, Bc, Vc, Rc, Ic,
 	             Us, Bs, Vs, Rs, Is,
                     U, B, V, R, I);
       Uc=U;
       Bc=B;
       Vc=V;
       Rc=R;
       Ic=I;
}
}
    }

    real L_tot=0, T_eff=0;
    cout << " Cluster: m= " << m_tot << " R= " << r_max << " L= " << L_tot
	 << " T_eff= " << T_eff 
	 << " r= "  << dc_pos[0] << " " << dc_pos[1] << " " << dc_pos[2] 
         << " ubvri= "<<  Uc << " " << Bc << " " << Vc << " " << Rc << " "
         << Ic	<< " :: ";
	
  }
}

main(int argc, char ** argv)
{
    check_help();

    bool C_flag = false;

    real r_max = VERY_LARGE_NUMBER;
    real r_min = 0;
    real v_max = -VERY_LARGE_NUMBER;
    real v_min = VERY_LARGE_NUMBER;

    bool verbose = false;
    int project = 0; // project along x, y, z or no projection [3]

    extern char *poptarg;
    int c;
    char* param_string = "CV:vR:r:p:";

    while ((c = pgetopt(argc, argv, param_string)) != -1)
	switch(c)
	    {
	    case 'C': C_flag = true;
	              break;
	    case 'V': v_max = atof(poptarg);
	              break;
		      //	    case 'v': v_min = atof(poptarg);
		      //	              break;
	    case 'v': verbose = !verbose;
	              break;
	    case 'R': r_max = atof(poptarg);
	              break;
	    case 'r': r_min = atof(poptarg);
	              break;
	    case 'p': project = atoi(poptarg);
	              break;
            case '?': params_to_usage(cerr, argv[0], param_string);
	              get_help();
                      exit(1);
	    }

    hdyn *b = NULL;
    int count = 0;

    bool cod, try_com = false;
    vec dc_pos = 0;
    vec dc_vel = 0;
    while (b = get_hdyn()) {

	real rd_min = b->get_starbase()->conv_r_star_to_dyn(r_min);
	real rd_max = b->get_starbase()->conv_r_star_to_dyn(r_max);

	if(verbose)
	  cout << "\n\nTime = "
	       << b->get_starbase()->conv_t_dyn_to_star(b->get_system_time())
	       << endl;

	compute_max_cod(b, dc_pos, dc_vel);

	real U,B,V,R,I;
	U=B=V=R=I=VERY_LARGE_NUMBER;

	PRC(dc_pos);PRL(dc_vel);
	dc_pos = 0;
	dc_vel = 0;
	if (!C_flag)
	   print_binary_recursive(b, dc_pos, dc_vel, U, B, V, R, I,
				  verbose);
	else
  	   print_integrated_cluster(b, dc_pos, rd_min, rd_max, v_min,
v_max,
				    U, B,
                                    V, R, I, project);

       rmtree(b);
    }
}

#endif
