#include "mmas.h"
#include <gsl/gsl_errno.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_odeiv.h>
#include <gsl/gsl_interp.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_roots.h>

struct merge_stars_consistently_params {
    mmas *mmas_ptr;
    double e_tot;
    int n_shells;
    double m_ejecta;
};

double merge_stars_consistently_eq(double f_heat, void *params) {
    merge_stars_consistently_params *p = (struct merge_stars_consistently_params*)params;
    int status = p->mmas_ptr->shock_heating(f_heat);
    if (status < 0) return -1e99;
    p->mmas_ptr->merge_stars(1.0, p->n_shells);
    real energy_p = p->mmas_ptr->compute_stellar_energy(p->mmas_ptr->get_product());
    real vesc = sqrt(p->mmas_ptr->get_product().star_mass / p->mmas_ptr->get_product().star_radius);
    real energy_ej = p->m_ejecta * vesc*vesc;
    PRL(energy_p); PRL(p->e_tot); PRL(energy_ej);
    cerr << " -~~*~~- \n";
    return (energy_p + energy_ej - p->e_tot); 
}

void mmas::merge_stars_consistently(int n_shells, int flag_do_shock_heating) {
    compute_extra();
    
    merge_stars_consistently_params p;
    p.mmas_ptr = this;
    
    real f_lost = mass_loss()/100.0;
    p.m_ejecta = (model_a->star_mass + model_b->star_mass)*f_lost;
    
    
    real energy_a = compute_stellar_energy(get_model_a());
    real energy_b = compute_stellar_energy(get_model_b());
    p.e_tot = energy_a + energy_b;
    p.n_shells = 1000;
    
    // specifies minimal and maximal heating
    real f_heat_min = 0.1;
    real f_heat_max = 10;
    
    gsl_function F;
    F.function = &merge_stars_consistently_eq;
    F.params   = &p;
    
    const gsl_root_fsolver_type *T = gsl_root_fsolver_brent;
    //   const gsl_root_fsolver_type *T = gsl_root_fsolver_falsepos;
    //   const gsl_root_fsolver_type *T = gsl_root_fsolver_bisection;
    gsl_root_fsolver *s            = gsl_root_fsolver_alloc(T);
    
    /* ----------- SOLVER ------------ */
    
    fprintf(stderr, "\n----------------------------------------------------\n");
    fprintf(stderr, "Computing f_heat using %s method:\n", 
    gsl_root_fsolver_name (s));
    
    // Calculating lower bound test interval
    double f1;// = merge_stars_consistently_eq(f_heat_min, F.params);
    do {
        f1 = merge_stars_consistently_eq(f_heat_min, F.params);
        if (f1 == -1e99) f_heat_min *= 2.0;
    } while (f1 == -1e99 && f_heat_min < 1.0);
    
    // Calculating upper bound test interval
    double f2;
    do {
        f2 = merge_stars_consistently_eq(f_heat_max, F.params);
        if (f2 == -1e99) f_heat_max /= 2.0;
    } while (f2 == -1e99 && f_heat_max > 1.0);
    if (f1 == -1e99 or f2 == -1e99) f2 = f1; // No shock heating
    
    if (f1*f2 < 0.0 and flag_do_shock_heating)    // if heating rate is appropriate, solve for f_heat, otherwise assume no heating
    {
        fprintf(stderr,"\n\n\n  Solving for f_heat ... \n\n\n");
        gsl_root_fsolver_set (s, &F, f_heat_min, f_heat_max);
        fprintf (stderr, "%5s [%9s, %9s] %9s %9s\n",
            "iter", "lower", "upper", "root", "err(est)");
        
        int status;
        int iter = 0, iter_max = 1000;
        double rel_err = 1.0e-3;
        do {
            status      = gsl_root_fsolver_iterate (s);
            double r    = gsl_root_fsolver_root    (s);
            double x_lo = gsl_root_fsolver_x_lower (s);
            double x_hi = gsl_root_fsolver_x_upper (s);
            
            status = gsl_root_test_interval (x_lo, x_hi, 0, rel_err);
            
            if (status == GSL_SUCCESS)
                fprintf (stderr, "Converged in %u iterations:\n", iter+1);
            
            fprintf (stderr, "***************************************\n");
            fprintf (stderr, "%5d [%.7f, %.7f] %.7f %.7f\n",
                iter, x_lo, x_hi, r, x_hi - x_lo);
            fprintf (stderr, "***************************************\n");
            iter++;
        } while (status == GSL_CONTINUE && iter < iter_max);
        fprintf(stderr, "----------------------------------------------------\n");
        
        double f_heat = gsl_root_fsolver_root(s);
        
        /* ----------- SOLVER ------------ */
        
        cerr << endl; 
        cerr << " **** Solved energy consistensy **** \n";
        PRL(f_heat);
        cerr << " **** Now I build hi-resolution model **** \n";
        cerr << endl;
        
        shock_heating(f_heat);
    }
    else
    {
        if (flag_do_shock_heating)
            fprintf(stderr,"\n\n\n  Failed to solve for f_heat ... Assuming no shock-heating.  \n\n\n");
        else
            fprintf(stderr,"\n\n\n  Assuming no shock-heating.  \n\n\n");
    }
    gsl_root_fsolver_free (s);
    merge_stars(1.0, n_shells);
}
