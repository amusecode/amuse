#include "types.h"

void Particle::set_ODE_quantities(double delta_time)
{
    for (int i=0; i<3; i++)
    {        
        dspin_vec_dt[i] = 0.0;
        de_vec_dt[i] = 0.0;
        dh_vec_dt[i] = 0.0;
    }
    
    e = norm3(e_vec);
    h = norm3(h_vec);
    spin_vec_norm = norm3(spin_vec);

    for (int i=0; i<3; i++)
    {        
        e_vec_unit[i] = e_vec[i]/e;
        h_vec_unit[i] = h_vec[i]/h;
    }
    
    e_p2 = e*e;
    j_p2 = 1.0 - e_p2;
    j = sqrt(j_p2);
    j_p3 = j*j_p2;
    j_p4 = j*j_p3;
    j_p5 = j*j_p4;

   
    a = h*h*child1_mass_plus_child2_mass/( CONST_G*child1_mass_times_child2_mass*child1_mass_times_child2_mass*j_p2 );

    /* assuming constant semimajor axis with mass loss */
    //double factor_h_vec = (child1_mass_dot_external/(2.0*child1_mass))*(1.0 + child2_mass/child1_mass_plus_child2_mass) + (child2_mass_dot_external/(2.0*child2_mass))*(1.0 + child1_mass/child1_mass_plus_child2_mass); 
    /* assuming constant SPECIFIC orbital angular momentum with mass loss, i.e. a(m1+m2)=const. */
    double factor_h_vec = child1_mass_dot_external/child1_mass + child2_mass_dot_external/child2_mass - (child1_mass_dot_external + child2_mass_dot_external)/child1_mass_plus_child2_mass;

    /* assuming constant spin angular momentum of the body */
    double factor_spin_vec = - (mass_dot_external/mass + 2.0*radius_dot_external/radius);

    /* set external time derivatives if appropriate */
    double h_vec_dot_external[3] = {h_vec_x_dot_external,h_vec_y_dot_external,h_vec_z_dot_external};
    double e_vec_dot_external[3] = {e_vec_x_dot_external,e_vec_y_dot_external,e_vec_z_dot_external};
    double spin_vec_dot_external[3] = {spin_vec_x_dot_external,spin_vec_y_dot_external,spin_vec_z_dot_external};

    for (int i=0; i<3; i++)
    {        
        dspin_vec_dt[i] = spin_vec[i]*factor_spin_vec + spin_vec_dot_external[i];
        de_vec_dt[i] = e_vec_dot_external[i];
        dh_vec_dt[i] = h_vec[i]*factor_h_vec + h_vec_dot_external[i];
    }
    dmass_dt = mass_dot_external;
    dradius_dt = radius_dot_external + radius_ddot_external*delta_time;

}

void Particle::reset_ODE_quantities()
{
    for (int i=0; i<3; i++)
    {        
        dspin_vec_dt[i] = 0.0;
        de_vec_dt[i] = 0.0;
        dh_vec_dt[i] = 0.0;
    }
}
