//
// gravity.C
//
// Version 1999/2/9 Jun Makino
//
// Change in the calling format for apply_vf, from function name to
// add explicit address operator. This was necessary to pass g++ 2.8.1
// 
// Version 1999/1/1 Jun Makino  --- cleaned up and some comments  added.
//
// The gravitational force calculation package designed for TREE
// implementation of SPH/NBODY programs
//
//
// Structure:
//
//  calculate_gravity (top level routine)
//      calculate_uncorrcted_gravity (calculate gravity with const softening)
//          setup_tree
//	    set_cm_quantities_for_default_tree
//	    calculate_gravity_using_tree,eps2_for_gravity(NON-GRAPE)
//          evaluate_gravity_using_default_tree_and_list(GRAPE)
//      corrcted_gravity (apply SPH form-factor correction)
//

#ifndef NOGRAPHICS
#define GRAPHICS
#endif

#include "stdinc.h"
#include "BHtree.h"

#ifdef TREE
void set_cm_quantities_for_default_tree();
#endif    
extern "C" double cpusec();
void clear_tree_counters();
void print_tree_counters();

#ifndef REAL_GRAVITY
void real_system::calculate_gravity()
{
    // This is peudo gravity: just harmonic force
    apply_vf(real_particle::set_harmonic_gravity,1.0);
    
}
#endif
#ifdef REAL_GRAVITY
void real_system::calculate_gravity()
{
    if (use_self_gravity){
	calculate_uncorrected_gravity();
	correct_gravity();
    }
}
#endif

void evaluate_gravity_using_default_tree_and_list(real theta2,
					  real eps2,
					  int ncrit);

void real_system::calculate_uncorrected_gravity()
{
#ifndef TREE    
  calculate_uncorrected_gravity_direct();
#else
//    cerr << "Call setup tree, cpu = " <<cpusec() << endl;
    setup_tree();
//    cerr << "Call set cm, cpu = " <<cpusec() << endl;
    set_cm_quantities_for_default_tree();
//    cerr << "Call evaluate_gravity, cpu = " <<cpusec() << endl;
    clear_tree_counters();
    apply_vf(&real_particle::clear_acc_phi_gravity);
    //    real (befrac::*func_ptr3)(const real, const real, const real)
    //       func_ptr = &befrac::n_br_wd;
    //    void calculate_vf2R = &real_particle::calculate_gravity_using_tree;
    //    apply_vf2R(calculate_vf2R,
    //	     eps2_for_gravity, theta_for_tree*theta_for_tree);
    
    apply_vf(&real_particle::calculate_gravity_using_tree,
	       eps2_for_gravity, theta_for_tree*theta_for_tree);
#endif
    
//    print_tree_counters();
//    cerr << "Exit evaluate_gravity, cpu = " <<cpusec() << endl;
    
}

void real_system::calculate_uncorrected_gravity_direct()
{
    apply_vf(&real_particle::clear_acc_phi_gravity);
    int i, j;
    real_particle * pi;
    real_particle * pj;
    for(i = 0,  pi = &(pb[0]); i<n-1; i++,pi++){
	for(j = i+1,  pj = pi+1; j<n; j++,pj++){
	    accumulate_mutual_gravity(*pi, *pj, eps2_for_gravity);
	}
    }
}

#if 0
template <class T>
void accumulate_mutual_gravity(T & p1,
			       T & p2,
			       real eps2)
#endif    
void accumulate_mutual_gravity(real_particle & p1,
			       real_particle & p2,
			       real eps2)
{
    vec dx = p1.pos-p2.pos;
    double r2inv = 1/(dx*dx+eps2);
    double rinv  = sqrt(r2inv);
    double r3inv = r2inv*rinv;
    p1.phi_gravity -= p2.mass*rinv;
    p2.phi_gravity -= p1.mass*rinv;
    p1.acc_gravity -= p2.mass*r3inv*dx;
    p2.acc_gravity += p1.mass*r3inv*dx;
}

#ifdef SPH
#if 1
// The next function perform the gravity correction with the simplest possible
// method, namely assuming that one particle is uniform density sphere and the
// other a point mass.
//
void apply_sph_correction_for_mutual_gravity(real_particle & p1,
					     real_particle & p2,
					     real eps2,
					     int correct_both)
{
    vec dx = p1.pos-p2.pos;
    real r2 = dx*dx;
    real scaledist = (p1.h + p2.h)*1;
    real scaledist2 = scaledist * scaledist ;

    if (r2 > scaledist2) return;
    real hsinv2 = 1.0/(scaledist2+eps2);
    real hsinv = sqrt(hsinv2);
    real phi_h = hsinv;
    real fcoef = hsinv2*hsinv;
    
    double r2inv = 1/(r2+eps2);
    double r = sqrt(r2);
    double rinv  = sqrt(r2inv);
    double r3inv = r2inv*rinv;
    real phifact = rinv-hsinv+0.5*fcoef*(r2-scaledist2);
    real ffact = r3inv - fcoef;

    p1.phi_gravity += p2.mass*phifact;
    p1.acc_gravity += p2.mass*ffact*dx;
    if(correct_both){
	p2.phi_gravity += p1.mass*phifact;
	p2.acc_gravity -= p1.mass*ffact*dx;
    }
}
#else


// Correction using spline kernel (assuming eps2 = 0)
//
void apply_real_correction_for_mutual_gravity(real_particle & p1,
					     real_particle & p2,
					     real eps2,
					     int correct_both)
{
    vec dx = p1.pos-p2.pos;
    real r2 = dx*dx;
    real scaledist = (p1.h + p2.h);
    real scaledist2 = scaledist * scaledist ;
    //    PRC(r2); PRL(scaledist2);
    if (r2 > scaledist2) return;
    real u2 = r2/scaledist2;
    real u = sqrt(u);
    real f, g;
    real scale3 = scaledist2*scaledist;
    double r2inv = 1/(r2+eps2);
    double r = sqrt(r2);
    double rinv  = sqrt(r2inv);
    double r3inv = r2inv*rinv;
    if  (u < 1){
	f = (1.4 - ((0.05*u-0.015)*u2+1.0/3.0)*u2*2)/scaledist;
	g = ((0.5*u-1.2)*u2+4.0/3.0)/scale3;
    }else{
	f = -rinv/15.0+(1.6-(((-u/30+0.3)*u-1)+4.0/3.0)*u2)/scaledist;
	g = r3inv*((((-u/6+1.2)*u-3)+8.0/3.0)*u2*u-1.0/15.0);
    }
    real phifact = rinv-f;
    real ffact = r3inv - g;
    //    PRC(phifact/rinv); PRC(ffact/r3inv); PRL(dx);
    p1.phi_gravity += p2.mass*phifact;
    p1.acc_gravity += p2.mass*ffact*dx;
    if(correct_both){
	p2.phi_gravity += p1.mass*phifact;
	p2.acc_gravity -= p1.mass*ffact*dx;
    }
}
#endif
#endif

#ifdef SPH
void real_particle::correct_gravity(real eps2)
{
    //    cerr << "correct_gravity for "; PRC(index); PRL(eps2);
    //    PRC(acc_gravity); PRL(phi_gravity);
    for(int j = 0; j<nnb; j++){
	real_particle * pj = pnb[j];
	if (pj->index > index){
	    apply_real_correction_for_mutual_gravity(*this, *pj, eps2, 1);
	}
    }
    //    PRC(acc_gravity); PRL(phi_gravity);
}
#else
void real_particle::correct_gravity(real eps2)
{
}
#endif
void real_system::correct_gravity()
{
    apply_vf(&real_particle::correct_gravity, eps2_for_gravity);
}

