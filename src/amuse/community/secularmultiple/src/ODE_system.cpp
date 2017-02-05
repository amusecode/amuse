#include "types.h"
#include "evolve.h"
#include "ODE_system.h"

int compute_y_dot(realtype time, N_Vector y, N_Vector y_dot, void *data_)
{
	UserData data;
	data = (UserData) data_;
    ParticlesMap *particlesMap = data->particlesMap;

    extract_ODE_variables(particlesMap, y, true);

    /****************************
     * compute right-hand sides *
     * **************************/

    double hamiltonian = 0.0;
    ParticlesMapIterator it_p;
    std::vector<int>::iterator it_parent_p,it_parent_q;
    
    for (it_p = particlesMap->begin(); it_p != particlesMap->end(); it_p++)
    {
        Particle *P_p = (*it_p).second;

        if (P_p->is_binary == 1)
        {
            
            /* Newtonian gravitational point mass dynamics */
            
            /* binary pairs */
            for (it_parent_p = P_p->parents.begin(); it_parent_p != P_p->parents.end(); it_parent_p++)
            {
                int i = std::distance(P_p->parents.begin(), it_parent_p);
                Particle *P_q = (*particlesMap)[(*it_parent_p)];
                int connecting_child_in_parent_q = P_p->connecting_child_in_parents[i];
                hamiltonian += compute_EOM_binary_pairs(particlesMap,P_p->index,P_q->index,connecting_child_in_parent_q,false);
                
                /* binary triplets */
                for (it_parent_q = P_q->parents.begin(); it_parent_q != P_q->parents.end(); it_parent_q++)
                {
                    int j = std::distance(P_q->parents.begin(), it_parent_q);
                    Particle *P_u = (*particlesMap)[(*it_parent_q)];
                    int connecting_child_in_parent_u = P_q->connecting_child_in_parents[j];
                    hamiltonian += compute_EOM_binary_triplets(particlesMap,P_p->index,P_q->index,P_u->index,connecting_child_in_parent_q,connecting_child_in_parent_u,false);
                    //printf("cross applied %d %d %d %d %d\n",P_p->index,P_q->index,P_u->index,connecting_child_in_parent_q,connecting_child_in_parent_u);

                }
            }
            
            /* Pairwise PN corrections */
            if (P_p->include_pairwise_1PN_terms == 1)
            {
                hamiltonian += compute_EOM_pairwise_1PN(particlesMap,P_p->index,false);
            }
            if (P_p->include_pairwise_25PN_terms == 1)
            {
                hamiltonian += compute_EOM_pairwise_25PN(particlesMap,P_p->index,false);
            }
            
            /* tidal friction (ad hoc) */
            Particle *P_child1 = (*particlesMap)[P_p->child1];
            Particle *P_child2 = (*particlesMap)[P_p->child2];

            if (P_child1->include_tidal_friction_terms == 1 || P_child1->include_tidal_bulges_precession_terms == 1 || P_child1->include_rotation_precession_terms == 1)
            {
                if (P_child1->tides_method == 0)
                {
                    compute_EOM_equilibrium_tide_BO(particlesMap,P_p->index,P_child1->index,P_child2->index,P_child1->include_tidal_friction_terms,P_child1->include_tidal_bulges_precession_terms,P_child1->include_rotation_precession_terms,P_child1->minimum_eccentricity_for_tidal_precession);
                }
                else if (P_child1->tides_method == 1)
                {
                    compute_EOM_equilibrium_tide_BO_full(particlesMap,P_p->index,P_child1->index,P_child2->index,P_child1->include_tidal_friction_terms,P_child1->include_tidal_bulges_precession_terms,P_child1->include_rotation_precession_terms,P_child1->minimum_eccentricity_for_tidal_precession);
                }
                else if (P_child1->tides_method == 2)
                {
                    compute_EOM_equilibrium_tide(particlesMap,P_p->index,P_child1->index,P_child2->index,P_child1->include_tidal_friction_terms,P_child1->include_tidal_bulges_precession_terms,P_child1->include_rotation_precession_terms,P_child1->minimum_eccentricity_for_tidal_precession);
                }
            }
            if (P_child2->include_tidal_friction_terms == 1 || P_child2->include_tidal_bulges_precession_terms == 1 || P_child2->include_rotation_precession_terms == 1)
            {
                if (P_child2->tides_method == 0)
                {
                    compute_EOM_equilibrium_tide_BO(particlesMap,P_p->index,P_child2->index,P_child1->index,P_child2->include_tidal_friction_terms,P_child2->include_tidal_bulges_precession_terms,P_child2->include_rotation_precession_terms,P_child2->minimum_eccentricity_for_tidal_precession);
                }
                else if (P_child2->tides_method == 1)
                {
                    compute_EOM_equilibrium_tide_BO_full(particlesMap,P_p->index,P_child2->index,P_child1->index,P_child2->include_tidal_friction_terms,P_child2->include_tidal_bulges_precession_terms,P_child2->include_rotation_precession_terms,P_child2->minimum_eccentricity_for_tidal_precession);
                }
                else if (P_child2->tides_method == 2)
                {
                    compute_EOM_equilibrium_tide(particlesMap,P_p->index,P_child2->index,P_child1->index,P_child2->include_tidal_friction_terms,P_child2->include_tidal_bulges_precession_terms,P_child2->include_rotation_precession_terms,P_child2->minimum_eccentricity_for_tidal_precession);
                }

            }
        }
            
    }
    
    write_ODE_variables_dots(particlesMap,y_dot);

    data->hamiltonian = hamiltonian;

    return 0;
}

double compute_EOM_binary_pairs(ParticlesMap *particlesMap, int inner_binary_index, int outer_binary_index, int connecting_child_in_outer_binary, bool compute_hamiltonian_only)
{

    /* stop if no triple terms are to be computed for the binary pair */
    if ((include_quadrupole_order_terms == false) && (include_octupole_order_binary_pair_terms == false) && (include_hexadecupole_order_binary_pair_terms == false) && (include_dotriacontupole_order_binary_pair_terms == false) )
    {
        return 0.0;
    }
    

    /*********************
     * preamble          *
     ********************/
    Particle *inner_binary = (*particlesMap)[inner_binary_index];
    Particle *outer_binary = (*particlesMap)[outer_binary_index];
    
    Particle *P_child1 = (*particlesMap)[inner_binary->child1];
    Particle *P_child2 = (*particlesMap)[inner_binary->child2];
    Particle *P_sibling;
    if (connecting_child_in_outer_binary==1)
    {
        P_sibling = (*particlesMap)[outer_binary->child2];
    }
    else if (connecting_child_in_outer_binary==2)
    {
        P_sibling = (*particlesMap)[outer_binary->child1];
    }
//    printf("compute_EOM_binary_pairs inner_binary_index %d outer_binary_index %d connecting_child_in_outer_binary %d P_sibling %d sibling_mass %g\n",inner_binary_index,outer_binary_index,connecting_child_in_outer_binary,P_sibling->index,P_sibling->mass);

    double e_in = inner_binary->e;
    double e_in_p2 = inner_binary->e_p2;
    double e_in_p4 = e_in_p2*e_in_p2;
    double e_out = outer_binary->e;
    double e_out_p2 = outer_binary->e_p2;
    
    double *e_in_vec = inner_binary->e_vec;
    double *e_out_vec = outer_binary->e_vec;
    double *h_in_vec = inner_binary->h_vec;
    double *h_out_vec = outer_binary->h_vec;
    
    double *e_in_vec_unit = inner_binary->e_vec_unit;
    double *e_out_vec_unit = outer_binary->e_vec_unit;
    double *h_in_vec_unit = inner_binary->h_vec_unit;
    double *h_out_vec_unit = outer_binary->h_vec_unit;
    
    double h_in = inner_binary->h;
    double h_out = outer_binary->h;
    
    double j_in = inner_binary->j;
    double j_in_p2 = inner_binary->j_p2;
    double j_in_p3 = inner_binary->j_p3;
    double j_out = outer_binary->j;
    double j_out_p2 = outer_binary->j_p2;
    double j_out_p3 = outer_binary->j_p3;
    double j_out_p4 = outer_binary->j_p4;
    double j_out_p5 = outer_binary->j_p5;
    double j_out_p6 = j_out*j_out_p5;
    double j_out_p7 = j_out*j_out_p6;
    double j_out_p8 = j_out*j_out_p7;
    double j_out_p9 = j_out*j_out_p8;
    double j_out_p10 = j_out*j_out_p9;
    double j_out_p11 = j_out*j_out_p10;
    double j_out_p13 = j_out_p2*j_out_p11;
    
    double j_out_p2_inv = 1.0/j_out_p2;
    double j_out_p5_inv = 1.0/j_out_p5;
    double j_out_p7_inv = 1.0/j_out_p7;
    double j_out_p9_inv = 1.0/j_out_p9;
    double j_out_p11_inv = 1.0/j_out_p11;
    double j_out_p13_inv = 1.0/j_out_p13;
    
    double j_in_vec[3],j_out_vec[3];
    for (int i=0; i<3; i++)
    {
        j_in_vec[i] = j_in*h_in_vec_unit[i];
        j_out_vec[i] = j_out*h_out_vec_unit[i];
    }
    
    double a_in = inner_binary->a;
    double a_out = outer_binary->a;
    
    double m1 = P_child1->mass;
    double m2 = P_child2->mass;
    double m3 = P_sibling->mass;
        
    double m1_plus_m2 = inner_binary->child1_mass_plus_child2_mass;
    double m1_minus_m2 = inner_binary->child1_mass_minus_child2_mass;
    double m1_times_m2 = inner_binary->child1_mass_times_child2_mass;

    double A_quad = c_1div8*CONST_G*(a_in*a_in/(a_out*a_out*a_out))*m1_times_m2*m3/m1_plus_m2;
    double A_oct = A_quad*c_15div8*(a_in/a_out)*fabs(m1_minus_m2)/m1_plus_m2;
    double A_hd = 0.0;
    double A_tc = 0.0;

    if (include_quadrupole_order_terms == false)
    {
        A_quad = 0.0;
    }
    if (include_octupole_order_binary_pair_terms == false)
    {
        A_oct = 0.0;
    }
    if (include_hexadecupole_order_binary_pair_terms == true)
    {
        A_hd = c_3div1024*CONST_G*(a_in*a_in*a_in*a_in/(a_out*a_out*a_out*a_out*a_out))*(m1_times_m2*m3*(m1*m1 - m1_times_m2 + m2*m2)/(m1_plus_m2*m1_plus_m2*m1_plus_m2));
    }
    if (include_dotriacontupole_order_binary_pair_terms == true)
    {
        A_tc = -c_105div4096*CONST_G*(a_in*a_in*a_in*a_in*a_in/(a_out*a_out*a_out*a_out*a_out*a_out))*(m1_times_m2*m3*fabs(m1_minus_m2)*(m1*m1 + m2*m2)/(m1_plus_m2*m1_plus_m2*m1_plus_m2*m1_plus_m2));
    }

    double Lambda_in = h_in/j_in;
    double Lambda_out = h_out/j_out;

    double e_in_vec_dot_e_out_vec = dot3(e_in_vec,e_out_vec);
    double j_in_vec_dot_j_out_vec = dot3(j_in_vec,j_out_vec);
    double e_in_vec_dot_j_out_vec = dot3(e_in_vec,j_out_vec);
    double j_in_vec_dot_e_out_vec = dot3(j_in_vec,e_out_vec);
    
    double e_in_vec_dot_e_out_vec_p2 = e_in_vec_dot_e_out_vec*e_in_vec_dot_e_out_vec;
    double j_in_vec_dot_j_out_vec_p2 = j_in_vec_dot_j_out_vec*j_in_vec_dot_j_out_vec;
    double e_in_vec_dot_j_out_vec_p2 = e_in_vec_dot_j_out_vec*e_in_vec_dot_j_out_vec;
    double j_in_vec_dot_e_out_vec_p2 = j_in_vec_dot_e_out_vec*j_in_vec_dot_e_out_vec;

    double j_in_vec_dot_j_out_vec_p4 = j_in_vec_dot_j_out_vec_p2*j_in_vec_dot_j_out_vec_p2;
    double e_in_vec_dot_j_out_vec_p4 = e_in_vec_dot_j_out_vec_p2*e_in_vec_dot_j_out_vec_p2;

    /* dotriacontupole */
    double e_in_vec_dot_e_out_vec_p3 = e_in_vec_dot_e_out_vec*e_in_vec_dot_e_out_vec_p2;
    double j_in_vec_dot_j_out_vec_p3 = j_in_vec_dot_j_out_vec*j_in_vec_dot_j_out_vec_p2;
    double e_in_vec_dot_j_out_vec_p3 = e_in_vec_dot_j_out_vec*e_in_vec_dot_j_out_vec_p2;
    

    /***************************
     * compute the Hamiltonian *
     **************************/

    double f1 = (1.0-6.0*e_in_p2)*j_out_p2 + 15.0*e_in_vec_dot_j_out_vec_p2 - 3.0*j_in_vec_dot_j_out_vec_p2;
	double f2 = (1.0-8.0*e_in_p2)*j_out_p2 + 35.0*e_in_vec_dot_j_out_vec_p2 - 5.0*j_in_vec_dot_j_out_vec_p2;
	double f3 = -10.0*e_in_vec_dot_j_out_vec*j_in_vec_dot_e_out_vec*j_in_vec_dot_j_out_vec;
    double f4,f5,f6,f7,f8,f9,f10,f11,f12; /* hexadecupole */
    double g,g1,g2,g3,h1,h2,h3,h4,h5; /* dotriacontupole */
    if (include_hexadecupole_order_binary_pair_terms == true)
    {
        f4 = -6.0 + e_out_p2 + 40.0*e_in_p2*(1.0 + 8.0*e_out_p2) - 20.0*e_in_p4*(8.0+15.0*e_out_p2);
        f5 = -2.0*j_out_p2 - e_in_p2*j_out_p2 + 21.0*e_in_vec_dot_j_out_vec_p2;
        f6 = (1.0 - 10.0*e_in_p2)*(4.0 + 3.0*e_out_p2);
        f7 = 8.0 + 6.0*e_out_p2 + e_in_p2*(6.0 + 29.0*e_out_p2);
        f8 = j_out_p2 + 13.0*e_in_p2*j_out_p2 - 7.0*j_in_vec_dot_j_out_vec_p2;
        f9 = -2.0 - 3.0*e_out_p2 + 4.0*e_in_p2*(5.0 + 3.0*e_out_p2);
        f10 = j_out_p2 - e_in_p2*j_out_p2 + 7.0*e_in_vec_dot_j_out_vec_p2;
        f11 = 2.0 + e_out_p2;
        f12 = 3.0*f4*j_out_p4 + 420.0*e_in_vec_dot_e_out_vec_p2*j_out_p2*f5 \
            - 5880.0*j_out_p2*e_in_vec_dot_e_out_vec*j_in_vec_dot_e_out_vec*e_in_vec_dot_j_out_vec*j_in_vec_dot_j_out_vec \
            + 5.0*( 28.0*j_out_p2*f6*e_in_vec_dot_j_out_vec_p2 - 6.0*j_out_p2*f7*j_in_vec_dot_j_out_vec_p2 \
                - 12.0*j_out_p2*f8*j_in_vec_dot_e_out_vec_p2 + 98.0*j_out_p2*f9*e_in_vec_dot_j_out_vec_p2 \
                - 441.0*f11*e_in_vec_dot_j_out_vec_p4 + 42.0*f11*f10*j_in_vec_dot_j_out_vec_p2 \
                - 21.0*f11*j_in_vec_dot_j_out_vec_p4);
    }
    if (include_dotriacontupole_order_binary_pair_terms == true)
    {
        h1 = (1.0 - 4.0*e_in_p2)*(8.0 + e_out_p2);
        h2 = 8.0 + 3.0*e_out_p2;
        h3 = -8.0 + e_out_p2 - 4.0*e_in_p4*(80.0 + 179.0*e_out_p2) + e_in_p2*(64.0 + 748.0*e_out_p2);
        h4 = -8.0 - 19.0*e_out_p2 + 6.0*e_in_p2*(16.0 + 5.0*e_out_p2);
        h5 = 8.0 + e_out_p2 - 2.0*e_in_p2*(16.0 + 29.0*e_out_p2); 
        
        g1 = (-26.0 + 15.0*e_in_p2)*j_out_p2 + 18.0*j_in_vec_dot_j_out_vec_p2 + 99.0*e_in_vec_dot_j_out_vec_p2;
        g2 = h1*j_out_p2 + 9.0*h2*e_in_vec_dot_j_out_vec_p2 + 6.0*j_out_p2*j_in_vec_dot_e_out_vec_p2 - 3.0*h2*j_in_vec_dot_j_out_vec_p2;
        g3 = h3*j_out_p4 - 693.0*h2*e_in_vec_dot_j_out_vec_p4 + 42.0*e_in_vec_dot_j_out_vec_p2*(h4*j_out_p2 + 9.0*h2*j_in_vec_dot_j_out_vec_p2) \
            + 14.0*h5*j_in_vec_dot_j_out_vec_p2*j_out_p2 - 21.0*h2*j_in_vec_dot_j_out_vec_p4 \
            - 28.0*j_in_vec_dot_e_out_vec_p2*j_out_p2*( (1.0 + 23.0*e_in_p2)*j_out_p2 - 9.0*j_in_vec_dot_j_out_vec_p2 );
        
        g = -3024.0*e_in_vec_dot_e_out_vec_p2*e_in_vec_dot_j_out_vec*j_in_vec_dot_e_out_vec*j_in_vec_dot_j_out_vec*j_out_p2 \
            + 28.0*j_out_p2*e_in_vec_dot_e_out_vec_p3*g1 + 28.0*e_in_vec_dot_j_out_vec*j_in_vec_dot_e_out_vec*j_in_vec_dot_j_out_vec*g2 \
            + e_in_vec_dot_e_out_vec*g3;

//        f12 = 3.0*f4*j_out_p4 + 420.0*e_in_vec_dot_e_out_vec_p2*j_out_p2*f5 \
            - 5880.0*j_out_p2*e_in_vec_dot_e_out_vec*j_in_vec_dot_e_out_vec*e_in_vec_dot_j_out_vec*j_in_vec_dot_j_out_vec \
            + 5.0*( 28.0*j_out_p2*f6*e_in_vec_dot_j_out_vec_p2 \
                - 6.0*j_out_p2*f7*j_in_vec_dot_j_out_vec_p2 );
    }
    
    double binary_pair_hamiltonian = A_quad*j_out_p5_inv*f1 - A_oct*j_out_p7_inv*( e_in_vec_dot_e_out_vec*f2 + f3 ) \
        + A_hd*j_out_p11_inv*f12 + A_tc*j_out_p13_inv*g;


    if (compute_hamiltonian_only == true)
    {
        return binary_pair_hamiltonian;
    }

    
    /****************************************
     * compute gradients of the Hamiltonian *
     ***************************************/    
    double grad_e_in_vec_phi[3],    grad_j_in_vec_phi[3];
    double grad_e_out_vec_phi[3],   grad_j_out_vec_phi[3];
    
    double grad_j_in_vec_f1[3],     grad_j_in_vec_f2[3],        grad_j_in_vec_f3[3];
    double grad_j_out_vec_f1[3],    grad_j_out_vec_f2[3],       grad_j_out_vec_f3[3];    
    double grad_e_in_vec_f1[3],     grad_e_in_vec_f2[3],        grad_e_in_vec_f3[3];
    double grad_e_out_vec_f3[3];

    /* triacontadipole */
    double grad_j_in_vec_g1[3],     grad_j_in_vec_g2[3],        grad_j_in_vec_g3[3];
    double grad_j_out_vec_g1[3],    grad_j_out_vec_g2[3],       grad_j_out_vec_g3[3];    
    double grad_e_in_vec_g1[3],     grad_e_in_vec_g2[3],        grad_e_in_vec_g3[3];
    double grad_e_out_vec_g2[3],    grad_e_out_vec_g3[3];


    for (int i=0; i<3; i++)
    {
        /* separate terms occurring in the gradients */
        if (include_quadrupole_order_terms == true)
        {
            grad_j_in_vec_f1[i] = -6.0*j_in_vec_dot_j_out_vec*j_out_vec[i];
            grad_j_out_vec_f1[i] = -6.0*j_in_vec_dot_j_out_vec*j_in_vec[i] + 2.0*(1.0-6.0*e_in_p2)*j_out_vec[i] \
                + 30.0*e_in_vec_dot_j_out_vec*e_in_vec[i];
            grad_e_in_vec_f1[i] = -12.0*j_out_p2*e_in_vec[i] + 30.0*e_in_vec_dot_j_out_vec*j_out_vec[i];
        }
        if (include_octupole_order_binary_pair_terms == true)
        {
            grad_j_in_vec_f2[i] = -10.0*j_in_vec_dot_j_out_vec*j_out_vec[i];
            grad_j_in_vec_f3[i] = -10.0*e_in_vec_dot_j_out_vec*( j_in_vec_dot_j_out_vec*e_out_vec[i] + j_in_vec_dot_e_out_vec*j_out_vec[i] );
            grad_j_out_vec_f2[i] = -10.0*j_in_vec_dot_j_out_vec*j_in_vec[i] + 2.0*(1.0-8.0*e_in_p2)*j_out_vec[i] \
                + 70.0*e_in_vec_dot_j_out_vec*e_in_vec[i];
            grad_j_out_vec_f3[i] = -10.0*j_in_vec_dot_e_out_vec*( j_in_vec_dot_j_out_vec*e_in_vec[i] + e_in_vec_dot_j_out_vec*j_in_vec[i] );
            grad_e_in_vec_f2[i] = -16.0*j_out_p2*e_in_vec[i] + 70.0*e_in_vec_dot_j_out_vec*j_out_vec[i];
            grad_e_in_vec_f3[i] = -10.0*j_in_vec_dot_e_out_vec*j_in_vec_dot_j_out_vec*j_out_vec[i];
            grad_e_out_vec_f3[i] = -10.0*e_in_vec_dot_j_out_vec*j_in_vec_dot_j_out_vec*j_in_vec[i];
        }
        if (include_dotriacontupole_order_binary_pair_terms == true)
        {
            grad_j_in_vec_g1[i] = 36.0*j_in_vec_dot_j_out_vec*j_out_vec[i];
            grad_j_in_vec_g2[i] = 12.0*j_out_p2*j_in_vec_dot_e_out_vec*e_out_vec[i] - 6.0*h2*j_in_vec_dot_j_out_vec*j_out_vec[i];
            grad_j_in_vec_g3[i] = 756.0*e_in_vec_dot_j_out_vec_p2*h2*j_in_vec_dot_j_out_vec*j_out_vec[i] + 28.0*h5*j_in_vec_dot_j_out_vec*j_out_p2*j_out_vec[i] \
                - 84.0*h2*j_in_vec_dot_j_out_vec_p3*j_out_vec[i] - 56.0*j_in_vec_dot_e_out_vec*j_out_p2*((1.0 + 23.0*e_in_p2)*j_out_p2 - 9.0*j_in_vec_dot_j_out_vec_p2)*e_out_vec[i] \
                + 504.0*j_in_vec_dot_e_out_vec_p2*j_out_p2*j_in_vec_dot_j_out_vec*j_out_vec[i];

            grad_j_out_vec_g1[i] = 2.0*(-26.0 + 15.0*e_in_p2)*j_out_vec[i] + 36.0*j_in_vec_dot_j_out_vec*j_in_vec[i] + 198.0*e_in_vec_dot_j_out_vec*e_in_vec[i];
            grad_j_out_vec_g2[i] = 2.0*h1*j_out_vec[i] + 18.0*h2*e_in_vec_dot_j_out_vec*e_in_vec[i] + 12.0*j_in_vec_dot_e_out_vec_p2*j_out_vec[i] \
                - 6.0*h2*j_in_vec_dot_j_out_vec*j_in_vec[i];
            grad_j_out_vec_g3[i] = 4.0*h3*j_out_p2*j_out_vec[i] - 2772.0*h2*e_in_vec_dot_j_out_vec_p3*e_in_vec[i] \
                + 84.0*e_in_vec_dot_j_out_vec*( h4*j_out_p2 + 9.0*h2*j_in_vec_dot_j_out_vec_p2 )*e_in_vec[i] \
                + 42.0*e_in_vec_dot_j_out_vec_p2*( 2.0*h4*j_out_vec[i] + 18.0*h2*j_in_vec_dot_j_out_vec*j_in_vec[i] ) \
                + 28.0*h5*( j_out_p2*j_in_vec_dot_j_out_vec*j_in_vec[i] + j_in_vec_dot_j_out_vec_p2*j_out_vec[i] ) \
                - 84.0*h2*j_in_vec_dot_j_out_vec_p3*j_in_vec[i] \
                - 28.0*j_in_vec_dot_e_out_vec_p2*( 2.0*j_out_vec[i]*( (1.0 + 23.0*e_in_p2)*j_out_p2 - 9.0*j_in_vec_dot_j_out_vec_p2 ) \
                    + j_out_p2*( 2.0*(1.0 + 23.0*e_in_p2)*j_out_vec[i] - 18.0*j_in_vec_dot_j_out_vec*j_in_vec[i]) );

            grad_e_in_vec_g1[i] = 30.0*j_out_p2*e_in_vec[i] + 198.0*e_in_vec_dot_j_out_vec*j_out_vec[i];
            grad_e_in_vec_g2[i] = -8.0*(8.0 + e_out_p2)*j_out_p2*e_in_vec[i] + 18.0*e_in_vec_dot_j_out_vec*h2*j_out_vec[i];
            grad_e_in_vec_g3[i] = j_out_p4*( -16.0*e_in_p2*(80.0 + 179.0*e_out_p2) + 2.0*(64.0 + 748.0*e_out_p2) )*e_in_vec[i] \
                - 2772.0*h2*e_in_vec_dot_j_out_vec_p3*j_out_vec[i] + 84.0*e_in_vec_dot_j_out_vec*( h4*j_out_p2 + 9.0*h2*j_in_vec_dot_j_out_vec_p2 )*j_out_vec[i] \
                + 504.0*e_in_vec_dot_j_out_vec_p2*j_out_p2*(16.0 + 5.0*e_out_p2)*e_in_vec[i] - 56.0*j_in_vec_dot_j_out_vec_p2*j_out_p2*(16.0 + 29.0*e_out_p2)*e_in_vec[i] \
                - 1288.0*j_in_vec_dot_e_out_vec_p2*j_out_p4*e_in_vec[i];
        
            grad_e_out_vec_g2[i] = 2.0*(1.0 - 4.0*e_in_p2)*j_out_p2*e_out_vec[i] + 54.0*e_in_vec_dot_j_out_vec_p2*e_out_vec[i] \
                + 12.0*j_out_p2*j_in_vec_dot_e_out_vec*j_in_vec[i] - 18.0*j_in_vec_dot_j_out_vec_p2*e_out_vec[i];
            grad_e_out_vec_g3[i] = j_out_p4*( 2.0 + 1496.0*e_in_p2 -1432.0*e_in_p4 )*e_out_vec[i] - 4158.0*e_in_vec_dot_j_out_vec_p4*e_out_vec[i] \
                + 42.0*e_in_vec_dot_j_out_vec_p2*( (-38.0 + 60.0*e_in_p2)*j_out_p2 + 54.0*j_in_vec_dot_j_out_vec_p2 )*e_out_vec[i] \
                + 14.0*j_in_vec_dot_j_out_vec_p2*j_out_p2*(2.0 - 116.0*e_in_p2)*e_out_vec[i] - 126.0*j_in_vec_dot_j_out_vec_p4*e_out_vec[i] \
                - 56.0*j_in_vec_dot_e_out_vec*j_out_p2*( (1.0 + 23.0*e_in_p2)*j_out_p2 - 9.0*j_in_vec_dot_j_out_vec_p2 )*j_in_vec[i];
        }
            
        /* complete gradients */
        grad_j_in_vec_phi[i] = 0.0;
        grad_j_out_vec_phi[i] = 0.0;
        grad_e_in_vec_phi[i] = 0.0;
        grad_e_out_vec_phi[i] = 0.0;
        
        if (include_quadrupole_order_terms == true)
        {
            grad_j_in_vec_phi[i] += A_quad*j_out_p5_inv*grad_j_in_vec_f1[i];
            grad_j_out_vec_phi[i] += -5.0*A_quad*j_out_p7_inv*j_out_vec[i]*f1 + A_quad*j_out_p5_inv*grad_j_out_vec_f1[i];
            grad_e_in_vec_phi[i] += A_quad*j_out_p5_inv*grad_e_in_vec_f1[i];
        }    
        if (include_octupole_order_binary_pair_terms == true)
        {
            grad_j_in_vec_phi[i] += -A_oct*j_out_p7_inv*( e_in_vec_dot_e_out_vec*grad_j_in_vec_f2[i] + grad_j_in_vec_f3[i] );
            grad_j_out_vec_phi[i] += 7.0*A_oct*j_out_p9_inv*j_out_vec[i]*( e_in_vec_dot_e_out_vec*f2 + f3 ) \
                - A_oct*j_out_p7_inv*( e_in_vec_dot_e_out_vec*grad_j_out_vec_f2[i] + grad_j_out_vec_f3[i] );
            grad_e_in_vec_phi[i] += -A_oct*j_out_p7_inv*( e_out_vec[i]*f2 + e_in_vec_dot_e_out_vec*grad_e_in_vec_f2[i] \
                + grad_e_in_vec_f3[i] );
            grad_e_out_vec_phi[i] += -A_oct*j_out_p7_inv*( e_in_vec[i]*f2 + grad_e_out_vec_f3[i] );
        }
        if (include_hexadecupole_order_binary_pair_terms == true)
        {
            grad_j_in_vec_phi[i] += A_hd*j_out_p11_inv*( \
                - 5880.0*j_out_p2*e_in_vec_dot_e_out_vec*e_in_vec_dot_j_out_vec*(j_in_vec_dot_j_out_vec*e_out_vec[i] + j_in_vec_dot_e_out_vec*j_out_vec[i]) \
                + 5.0*( -12.0*j_out_p2*f7*j_in_vec_dot_j_out_vec*j_out_vec[i] - 12.0*j_out_p2*(2.0*f8*j_in_vec_dot_e_out_vec*e_out_vec[i] \
                        - 14.0*j_in_vec_dot_e_out_vec_p2*j_in_vec_dot_j_out_vec*j_out_vec[i]) \
                    + 84.0*f11*f10*j_in_vec_dot_j_out_vec*j_out_vec[i] \
                    - 84.0*f11*j_in_vec_dot_j_out_vec_p2*j_in_vec_dot_j_out_vec*j_out_vec[i] ) \
                );
            grad_j_out_vec_phi[i] += -11.0*A_hd*j_out_p11_inv*j_out_p2_inv*f12*j_out_vec[i] \
                + A_hd*j_out_p11_inv*(12.0*f4*j_out_p2*j_out_vec[i] \
                    + 420.0*e_in_vec_dot_e_out_vec_p2*(2.0*f5*j_out_vec[i] + j_out_p2*(-4.0*j_out_vec[i] - 2.0*e_in_p2*j_out_vec[i] + 42.0*e_in_vec_dot_j_out_vec*e_in_vec[i])) \
                    - 5880.0*e_in_vec_dot_e_out_vec*j_in_vec_dot_e_out_vec*(2.0*e_in_vec_dot_j_out_vec*j_in_vec_dot_j_out_vec*j_out_vec[i] \
                        + j_out_p2*j_in_vec_dot_j_out_vec*e_in_vec[i] + j_out_p2*e_in_vec_dot_j_out_vec*j_in_vec[i]) \
                    + 5.0*( \
                        + 56.0*f6*(e_in_vec_dot_j_out_vec_p2*j_out_vec[i] + j_out_p2*e_in_vec_dot_j_out_vec*e_in_vec[i]) \
                        - 12.0*f7*(j_in_vec_dot_j_out_vec_p2*j_out_vec[i] + j_out_p2*j_in_vec_dot_j_out_vec*j_in_vec[i]) \
                        - 12.0*j_in_vec_dot_e_out_vec_p2*(2.0*f8*j_out_vec[i] + j_out_p2*(2.0*j_out_vec[i] + 26.0*e_in_p2*j_out_vec[i] - 14.0*j_in_vec_dot_j_out_vec*j_in_vec[i]) ) \
                        + 196.0*f9*(e_in_vec_dot_j_out_vec_p2*j_out_vec[i] + j_out_p2*e_in_vec_dot_j_out_vec*e_in_vec[i]) \
                        - 1764.0*f11*e_in_vec_dot_j_out_vec_p2*e_in_vec_dot_j_out_vec*e_in_vec[i] \
                        + 42.0*f11*( j_in_vec_dot_j_out_vec_p2*(2.0*j_out_vec[i] - 2.0*e_in_p2*j_out_vec[i] + 14.0*e_in_vec_dot_j_out_vec*e_in_vec[i]) \
                            + 2.0*f10*j_in_vec_dot_j_out_vec*j_in_vec[i] ) - 84.0*f11*j_in_vec_dot_j_out_vec_p2*j_in_vec_dot_j_out_vec*j_in_vec[i] ) \
                );
            grad_e_in_vec_phi[i] += A_hd*j_out_p11_inv*( \
                + 240.0*j_out_p4*(1.0 + 8.0*e_out_p2 - e_in_p2*(8.0 + 15.0*e_out_p2))*e_in_vec[i] + 840.0*e_in_vec_dot_e_out_vec*j_out_p2*f5*e_out_vec[i] \
                + 420.0*e_in_vec_dot_e_out_vec_p2*j_out_p2*(-2.0*j_out_p2*e_in_vec[i] + 42.0*e_in_vec_dot_j_out_vec*j_out_vec[i]) \
                - 5880.0*j_out_p2*j_in_vec_dot_e_out_vec*j_in_vec_dot_j_out_vec*( e_in_vec_dot_j_out_vec*e_out_vec[i] + e_in_vec_dot_e_out_vec*j_out_vec[i] ) \
                + 5.0*( \
                    + 28.0*j_out_p2*(4.0 + 3.0*e_out_p2)*( -20.0*e_in_vec_dot_j_out_vec_p2*e_in_vec[i] + 2.0*(1.0 - 10.0*e_in_p2)*e_in_vec_dot_j_out_vec*j_out_vec[i] ) \
                    - 12.0*j_out_p2*(6.0 + 29.0*e_out_p2)*j_in_vec_dot_j_out_vec_p2*e_in_vec[i] - 312.0*j_out_p4*j_in_vec_dot_e_out_vec_p2*e_in_vec[i] \
                    + 98.0*j_out_p2*(8.0*e_in_vec_dot_j_out_vec_p2*(5.0 + 3.0*e_out_p2)*e_in_vec[i] + 2.0*e_in_vec_dot_j_out_vec*f9*j_out_vec[i]) \
                    - 1764.0*f11*e_in_vec_dot_j_out_vec_p2*e_in_vec_dot_j_out_vec*j_out_vec[i] \
                    + 42.0*f11*j_in_vec_dot_j_out_vec_p2*(-2.0*j_out_p2*e_in_vec[i] + 14.0*e_in_vec_dot_j_out_vec*j_out_vec[i]) ) \
                );
            grad_e_out_vec_phi[i] += A_hd*j_out_p11_inv*( \
                + 6.0*j_out_p4*(1.0 + 320.0*e_in_p2 - 300.0*e_in_p4)*e_out_vec[i] + 840.0*e_in_vec_dot_e_out_vec*j_out_p2*f5*e_in_vec[i] \
                - 5880.0*j_out_p2*e_in_vec_dot_j_out_vec*j_in_vec_dot_j_out_vec*(j_in_vec_dot_e_out_vec*e_in_vec[i] + e_in_vec_dot_e_out_vec*j_in_vec[i]) \
                + 5.0*( \
                    + 168.0*j_out_p2*(1.0 - 10.0*e_in_p2)*e_in_vec_dot_j_out_vec_p2*e_out_vec[i] - 6.0*j_out_p2*j_in_vec_dot_j_out_vec_p2*(12.0 + 58.0*e_in_p2)*e_out_vec[i] \
                    - 24.0*j_out_p2*f8*j_in_vec_dot_e_out_vec*j_in_vec[i] + 98.0*j_out_p2*e_in_vec_dot_j_out_vec_p2*(-6.0 + 24.0*e_in_p2)*e_out_vec[i] \
                    - 882.0*e_in_vec_dot_j_out_vec_p4*e_out_vec[i] + 84.0*f10*j_in_vec_dot_j_out_vec_p2*e_out_vec[i] - 42.0*j_in_vec_dot_j_out_vec_p4*e_out_vec[i]) \
                );
        }
        if (include_dotriacontupole_order_binary_pair_terms == true)
        {
            grad_j_in_vec_phi[i] += A_tc*j_out_p13_inv*( \
                - 3024.0*e_in_vec_dot_e_out_vec_p2*e_in_vec_dot_j_out_vec*j_out_p2*( j_in_vec_dot_j_out_vec*e_out_vec[i] + j_in_vec_dot_e_out_vec*j_out_vec[i] ) \
                + 28.0*j_out_p2*e_in_vec_dot_e_out_vec_p3*grad_j_in_vec_g1[i] + 28.0*e_in_vec_dot_j_out_vec*( j_in_vec_dot_j_out_vec*g2*e_out_vec[i] \
                    + j_in_vec_dot_e_out_vec*g2*j_out_vec[i] + j_in_vec_dot_e_out_vec*j_in_vec_dot_j_out_vec*grad_j_in_vec_g2[i] ) \
                    + e_in_vec_dot_e_out_vec*grad_j_in_vec_g3[i] );
            grad_j_out_vec_phi[i] += -13.0*A_tc*j_out_p13_inv*j_out_p2_inv*g*j_out_vec[i] + A_tc*j_out_p13_inv*( \
                - 3024.0*e_in_vec_dot_e_out_vec_p2*j_in_vec_dot_e_out_vec*( j_in_vec_dot_j_out_vec*j_out_p2*e_in_vec[i] + e_in_vec_dot_j_out_vec*j_out_p2*j_in_vec[i] \
                    + 2.0*e_in_vec_dot_j_out_vec*j_in_vec_dot_j_out_vec*j_out_vec[i] ) + 28.0*e_in_vec_dot_e_out_vec_p3*( 2.0*g1*j_out_vec[i] + j_out_p2*grad_j_out_vec_g1[i] ) \
                    + 28.0*j_in_vec_dot_e_out_vec*( j_in_vec_dot_j_out_vec*g2*e_in_vec[i] + e_in_vec_dot_j_out_vec*g2*j_in_vec[i] \
                        + e_in_vec_dot_j_out_vec*j_in_vec_dot_j_out_vec*grad_j_out_vec_g2[i] ) + e_in_vec_dot_e_out_vec*grad_j_out_vec_g3[i] );
            grad_e_in_vec_phi[i] += A_tc*j_out_p13_inv*( \
                - 3024.0*j_in_vec_dot_e_out_vec*j_in_vec_dot_j_out_vec*j_out_p2*( 2.0*e_in_vec_dot_e_out_vec*e_in_vec_dot_j_out_vec*e_out_vec[i] \
                    + e_in_vec_dot_e_out_vec_p2*j_out_vec[i] ) + 84.0*j_out_p2*e_in_vec_dot_e_out_vec_p2*g1*e_out_vec[i] \
                + 28.0*j_out_p2*e_in_vec_dot_e_out_vec_p3*grad_e_in_vec_g1[i] + 28.0*j_in_vec_dot_e_out_vec*j_in_vec_dot_j_out_vec*( g2*j_out_vec[i] \
                    + e_in_vec_dot_j_out_vec*grad_e_in_vec_g2[i] ) + g3*e_out_vec[i] + e_in_vec_dot_e_out_vec*grad_e_in_vec_g3[i] );
            grad_e_out_vec_phi[i] += A_tc*j_out_p13_inv*( \
                - 3024.0*e_in_vec_dot_j_out_vec*j_in_vec_dot_j_out_vec*j_out_p2*( 2.0*j_in_vec_dot_e_out_vec*e_in_vec_dot_e_out_vec*e_in_vec[i] \
                    + e_in_vec_dot_e_out_vec_p2*j_in_vec[i] ) + 84.0*j_out_p2*e_in_vec_dot_e_out_vec_p2*g1*e_in_vec[i] \
                + 28.0*e_in_vec_dot_j_out_vec*j_in_vec_dot_j_out_vec*( g2*j_in_vec[i] + j_in_vec_dot_e_out_vec*grad_e_out_vec_g2[i] ) \
                + g3*e_in_vec[i] + e_in_vec_dot_e_out_vec*grad_e_out_vec_g3[i] );
        }
    }
    
    double j_in_vec_cross_grad_j_in_vec_phi[3],             j_in_vec_cross_grad_e_in_vec_phi[3];
    double j_out_vec_cross_grad_j_out_vec_phi[3],           j_out_vec_cross_grad_e_out_vec_phi[3];
    
    double e_in_vec_cross_grad_e_in_vec_phi[3],             e_in_vec_cross_grad_j_in_vec_phi[3];
    double e_out_vec_cross_grad_e_out_vec_phi[3],           e_out_vec_cross_grad_j_out_vec_phi[3];
    

    cross3(j_in_vec,        grad_j_in_vec_phi,              j_in_vec_cross_grad_j_in_vec_phi);
    cross3(j_in_vec,        grad_e_in_vec_phi,              j_in_vec_cross_grad_e_in_vec_phi);
    cross3(e_in_vec,        grad_e_in_vec_phi,              e_in_vec_cross_grad_e_in_vec_phi);
    cross3(e_in_vec,        grad_j_in_vec_phi,              e_in_vec_cross_grad_j_in_vec_phi);
    
    cross3(j_out_vec,       grad_j_out_vec_phi,             j_out_vec_cross_grad_j_out_vec_phi);
    cross3(j_out_vec,       grad_e_out_vec_phi,             j_out_vec_cross_grad_e_out_vec_phi);
    cross3(e_out_vec,       grad_e_out_vec_phi,             e_out_vec_cross_grad_e_out_vec_phi);
    cross3(e_out_vec,       grad_j_out_vec_phi,             e_out_vec_cross_grad_j_out_vec_phi);
    

    for (int i=0; i<3; i++)
    {
        inner_binary->de_vec_dt[i] += (-1.0/(Lambda_in))*( e_in_vec_cross_grad_j_in_vec_phi[i] \
            + j_in_vec_cross_grad_e_in_vec_phi[i] );
        inner_binary->dh_vec_dt[i] += -1.0*( j_in_vec_cross_grad_j_in_vec_phi[i] \
            + e_in_vec_cross_grad_e_in_vec_phi[i] );

        outer_binary->de_vec_dt[i] += (-1.0/(Lambda_out))*( e_out_vec_cross_grad_j_out_vec_phi[i] \
            + j_out_vec_cross_grad_e_out_vec_phi[i] );
        outer_binary->dh_vec_dt[i] += -1.0*( j_out_vec_cross_grad_j_out_vec_phi[i] \
            + e_out_vec_cross_grad_e_out_vec_phi[i] );  
    }

    return binary_pair_hamiltonian;
    
    if (1==0)
    {
        printf("e_in %g %g %g\n",e_in_vec[0],e_in_vec[1],e_in_vec[2]);
        printf("e_out %g %g %g\n",e_out_vec[0],e_out_vec[1],e_out_vec[2]);    
        printf("h_in %g %g %g\n",h_in_vec[0],h_in_vec[1],h_in_vec[2]);
        printf("h_out %g %g %g\n",h_out_vec[0],h_out_vec[1],h_out_vec[2]);    

        printf("grad1 %g %g %g\n",grad_e_in_vec_f1[0],grad_e_in_vec_f1[1],grad_e_in_vec_f1[2]);
        printf("grad2 %g %g %g\n",grad_e_in_vec_f2[0],grad_e_in_vec_f2[1],grad_e_in_vec_f2[2]);    
        printf("grad3 %g %g %g\n",grad_e_in_vec_f3[0],grad_e_in_vec_f3[1],grad_e_in_vec_f3[2]);    

        printf("de_in_dt %g %g %g\n",inner_binary->de_vec_dt[0],inner_binary->de_vec_dt[1],inner_binary->de_vec_dt[2]);
        printf("de_out_dt %g %g %g\n",outer_binary->de_vec_dt[0],outer_binary->de_vec_dt[1],outer_binary->de_vec_dt[2]);    
        printf("dh_in_dt %g %g %g\n",inner_binary->dh_vec_dt[0],inner_binary->dh_vec_dt[1],inner_binary->dh_vec_dt[2]);
        printf("dh_out_dt %g %g %g\n",outer_binary->dh_vec_dt[0],outer_binary->dh_vec_dt[1],outer_binary->dh_vec_dt[2]); 
    }
}



double compute_EOM_binary_triplets(ParticlesMap *particlesMap, int binary_A_index, int binary_B_index, int binary_C_index, int connecting_child_in_binary_B_to_binary_A, int connecting_child_in_binary_C_to_binary_B, bool compute_hamiltonian_only)
{
    
    if (include_octupole_order_binary_triplet_terms == false)
    {
        return 0.0;
    }

    /*********************
     * preamble          *
     ********************/

    Particle *binary_A = (*particlesMap)[binary_A_index];
    Particle *binary_B = (*particlesMap)[binary_B_index];
    Particle *binary_C = (*particlesMap)[binary_C_index];

    Particle *binary_A_child1 = (*particlesMap)[binary_A->child1];
    Particle *binary_A_child2 = (*particlesMap)[binary_A->child2];

    Particle *binary_B_child1 = (*particlesMap)[binary_B->child1];
    Particle *binary_B_child2 = (*particlesMap)[binary_B->child2];

    Particle *binary_C_child1 = (*particlesMap)[binary_C->child1];
    Particle *binary_C_child2 = (*particlesMap)[binary_C->child2];
   
    double B_ijB = 0.0;
    
    if (connecting_child_in_binary_B_to_binary_A==1)
    {
        B_ijB = binary_B_child2->mass/binary_B->mass;
    }
    else if (connecting_child_in_binary_B_to_binary_A==2)
    {
        B_ijB = -binary_B_child1->mass/binary_B->mass;
    }

    double M_C_CS_B = 0.0;
    
    if (connecting_child_in_binary_C_to_binary_B==1)
    {
        M_C_CS_B = binary_C_child2->mass;
    }
    else if (connecting_child_in_binary_C_to_binary_B==2)
    {
        M_C_CS_B = binary_C_child1->mass;
    }

    double e_A = binary_A->e;
    double e_B = binary_B->e;
    double e_C = binary_C->e;
    double e_A_p2 = binary_A->e_p2;
    double e_B_p2 = binary_B->e_p2;
    double e_C_p2 = binary_C->e_p2;
    
    double *e_A_vec = binary_A->e_vec;
    double *e_B_vec = binary_B->e_vec;
    double *e_C_vec = binary_C->e_vec;
    
    double *h_A_vec = binary_A->h_vec;
    double *h_B_vec = binary_B->h_vec;
    double *h_C_vec = binary_C->h_vec;

    double *e_A_vec_unit = binary_A->e_vec_unit;
    double *e_B_vec_unit = binary_B->e_vec_unit;
    double *e_C_vec_unit = binary_C->e_vec_unit;
    
    double *h_A_vec_unit = binary_A->h_vec_unit;
    double *h_B_vec_unit = binary_B->h_vec_unit;
    double *h_C_vec_unit = binary_C->h_vec_unit;
    
    double *j_A_vec_unit = h_A_vec_unit;
    double *j_B_vec_unit = h_B_vec_unit;    
    double *j_C_vec_unit = h_C_vec_unit;    

    double h_A = binary_A->h;
    double h_B = binary_B->h;
    double h_C = binary_C->h;
    
    double j_A = binary_A->j;
    double j_A_p2 = binary_A->j_p2;
    double j_B = binary_B->j;
    double j_C = binary_C->j;
    double j_C_p2 = binary_C->j_p2;
    double j_C_p4 = binary_C->j_p4;
    double j_C_p7 = j_C*j_C_p2*j_C_p4;
    double j_C_p9 = j_C_p7*j_C_p2;
        
    double j_C_p7_inv = 1.0/j_C_p7;
    double j_C_p9_inv = 1.0/j_C_p9;
    
    double j_A_vec[3],j_B_vec[3],j_C_vec[3];
    for (int i=0; i<3; i++)
    {
        j_A_vec[i] = j_A*h_A_vec_unit[i];
        j_B_vec[i] = j_B*h_B_vec_unit[i];
        j_C_vec[i] = j_C*h_C_vec_unit[i];
    }
    
    double a_A = binary_A->a;
    double a_B = binary_B->a;
    double a_C = binary_C->a;    
    
    double M_A1 = binary_A_child1->mass;
    double M_A2 = binary_A_child2->mass;

    double A_cross = -(c_9div32*M_A1*M_A2*B_ijB*M_C_CS_B/(M_A1 + M_A2))*(a_A*a_A*a_B/(a_C*a_C*a_C*a_C));
    double Lambda_A = h_A/j_A;
    double Lambda_B = h_B/j_B;
    double Lambda_C = h_C/j_C;

    double e_A_vec_dot_e_B_vec = dot3(e_A_vec,e_B_vec);
    double e_B_vec_dot_e_C_vec = dot3(e_B_vec,e_C_vec);
    double e_A_vec_dot_e_C_vec = dot3(e_A_vec,e_C_vec);

    double e_A_vec_dot_j_C_vec = dot3(e_A_vec,j_C_vec);
    double e_B_vec_dot_j_C_vec = dot3(e_B_vec,j_C_vec);
    double e_B_vec_dot_j_A_vec = dot3(e_B_vec,j_A_vec);
    double e_C_vec_dot_j_A_vec = dot3(e_C_vec,j_A_vec);
    double j_A_vec_dot_j_C_vec = dot3(j_A_vec,j_C_vec);
    
    double e_A_vec_dot_j_C_vec_p2 = e_A_vec_dot_j_C_vec*e_A_vec_dot_j_C_vec;
    double j_A_vec_dot_j_C_vec_p2 = j_A_vec_dot_j_C_vec*j_A_vec_dot_j_C_vec;
    
    /***************************
     * compute the Hamiltonian *
     **************************/
    
    double f1 = j_C_p2*(1.0 - 6.0*e_A_p2) + 25.0*e_A_vec_dot_j_C_vec_p2 - 5.0*j_A_vec_dot_j_C_vec_p2;
    double f0 = -10.0*e_A_vec_dot_e_B_vec*e_A_vec_dot_e_C_vec*j_C_p2 + 50.0*e_A_vec_dot_e_C_vec*e_A_vec_dot_j_C_vec*e_B_vec_dot_j_C_vec \
        + 2.0*e_C_vec_dot_j_A_vec*e_B_vec_dot_j_A_vec*j_C_p2 - 10.0*e_B_vec_dot_j_C_vec*e_C_vec_dot_j_A_vec*j_A_vec_dot_j_C_vec \
        + e_B_vec_dot_e_C_vec*f1;
            
    double binary_triplet_hamiltonian = A_cross*j_C_p7_inv*f0;
    
    if (compute_hamiltonian_only == true)
    {
        return binary_triplet_hamiltonian;
    }

    /****************************************
     * compute gradients of the Hamiltonian *
     ***************************************/
    double grad_e_A_vec_H[3],     grad_j_A_vec_H[3];
    double grad_e_B_vec_H[3],     grad_j_B_vec_H[3];
    double grad_e_C_vec_H[3],     grad_j_C_vec_H[3];    
    
    for (int i=0; i<3; i++)
    {
        
        /* gradient w.r.t. e_A */
        grad_e_A_vec_H[i] = A_cross*j_C_p7_inv*( \
            - 10.0*j_C_p2*(e_A_vec_dot_e_C_vec*e_B_vec[i] + e_A_vec_dot_e_B_vec*e_C_vec[i]) \
            + 50.0*e_B_vec_dot_j_C_vec*(e_A_vec_dot_j_C_vec*e_C_vec[i] + e_A_vec_dot_e_C_vec*j_C_vec[i]) \
            + e_B_vec_dot_e_C_vec*(50.0*e_A_vec_dot_j_C_vec*j_C_vec[i] - 12.0*j_C_p2*e_A_vec[i]) );
        
        /* gradient w.r.t. j_A */
        grad_j_A_vec_H[i] = A_cross*j_C_p7_inv*( \
            + 2.0*e_B_vec_dot_j_A_vec*j_C_p2*e_C_vec[i] + 2.0*e_C_vec_dot_j_A_vec*j_C_p2*e_B_vec[i] \
            - 10.0*e_B_vec_dot_j_C_vec*(j_A_vec_dot_j_C_vec*e_C_vec[i] + e_C_vec_dot_j_A_vec*j_C_vec[i]) \
            - 10.0*e_B_vec_dot_e_C_vec*j_A_vec_dot_j_C_vec*j_C_vec[i] );
            
        /* gradient w.r.t. e_B */
        grad_e_B_vec_H[i] = A_cross*j_C_p7_inv*( \
            - 10.0*e_A_vec_dot_e_C_vec*j_C_p2*e_A_vec[i] + 50.0*e_A_vec_dot_e_C_vec*e_A_vec_dot_j_C_vec*j_C_vec[i] \
            + 2.0*e_C_vec_dot_j_A_vec*j_C_p2*j_A_vec[i] - 10.0*e_C_vec_dot_j_A_vec*j_A_vec_dot_j_C_vec*j_C_vec[i] \
            + f1*e_C_vec[i] );
            
        /* gradient w.r.t. j_B */
        grad_j_B_vec_H[i] = 0.0;

        /* gradient w.r.t. e_C */
        grad_e_C_vec_H[i] = A_cross*j_C_p7_inv*( \
            - 10.0*e_A_vec_dot_e_B_vec*j_C_p2*e_A_vec[i] + 50.0*e_A_vec_dot_j_C_vec*e_B_vec_dot_j_C_vec*e_A_vec[i] \
            + 2.0*e_B_vec_dot_j_A_vec*j_C_p2*j_A_vec[i] - 10.0*e_B_vec_dot_j_C_vec*j_A_vec_dot_j_C_vec*j_A_vec[i] \
            + f1*e_B_vec[i] );

        /* gradient w.r.t. j_C */
        grad_j_C_vec_H[i] = -7.0*A_cross*j_C_p9_inv*f0*j_C_vec[i] + A_cross*j_C_p7_inv*( \
            - 20.0*e_A_vec_dot_e_B_vec*e_A_vec_dot_e_C_vec*j_C_vec[i] \
            + 50.0*e_A_vec_dot_e_C_vec*(e_B_vec_dot_j_C_vec*e_A_vec[i] + e_A_vec_dot_j_C_vec*e_B_vec[i]) \
            + 4.0*e_C_vec_dot_j_A_vec*e_B_vec_dot_j_A_vec*j_C_vec[i] \
            - 10.0*e_C_vec_dot_j_A_vec*(j_A_vec_dot_j_C_vec*e_B_vec[i] + e_B_vec_dot_j_C_vec*j_A_vec[i]) \
            + e_B_vec_dot_e_C_vec*(2.0*(1.0 - 6.0*e_A_p2)*j_C_vec[i] + 50.0*e_A_vec_dot_j_C_vec*e_A_vec[i] \
                - 10.0*j_A_vec_dot_j_C_vec*j_A_vec[i]) );
            
    }
    
    double j_A_vec_cross_grad_j_A_vec_H[3],                   j_A_vec_cross_grad_e_A_vec_H[3];
    double j_B_vec_cross_grad_j_B_vec_H[3],                   j_B_vec_cross_grad_e_B_vec_H[3];    
    double j_C_vec_cross_grad_j_C_vec_H[3],                   j_C_vec_cross_grad_e_C_vec_H[3];        
    
    double e_A_vec_cross_grad_e_A_vec_H[3],                   e_A_vec_cross_grad_j_A_vec_H[3];
    double e_B_vec_cross_grad_e_B_vec_H[3],                   e_B_vec_cross_grad_j_B_vec_H[3];
    double e_C_vec_cross_grad_e_C_vec_H[3],                   e_C_vec_cross_grad_j_C_vec_H[3];
    
    cross3(j_A_vec,             grad_j_A_vec_H,               j_A_vec_cross_grad_j_A_vec_H);
    cross3(j_A_vec,             grad_e_A_vec_H,               j_A_vec_cross_grad_e_A_vec_H);
    cross3(j_B_vec,             grad_j_B_vec_H,               j_B_vec_cross_grad_j_B_vec_H);
    cross3(j_B_vec,             grad_e_B_vec_H,               j_B_vec_cross_grad_e_B_vec_H);
    cross3(j_C_vec,             grad_j_C_vec_H,               j_C_vec_cross_grad_j_C_vec_H);
    cross3(j_C_vec,             grad_e_C_vec_H,               j_C_vec_cross_grad_e_C_vec_H);
    
    cross3(e_A_vec,             grad_e_A_vec_H,               e_A_vec_cross_grad_e_A_vec_H);
    cross3(e_A_vec,             grad_j_A_vec_H,               e_A_vec_cross_grad_j_A_vec_H);    
    cross3(e_B_vec,             grad_e_B_vec_H,               e_B_vec_cross_grad_e_B_vec_H);
    cross3(e_B_vec,             grad_j_B_vec_H,               e_B_vec_cross_grad_j_B_vec_H);    
    cross3(e_C_vec,             grad_e_C_vec_H,               e_C_vec_cross_grad_e_C_vec_H);
    cross3(e_C_vec,             grad_j_C_vec_H,               e_C_vec_cross_grad_j_C_vec_H);
    
    for (int i=0; i<3; i++)
    {
        binary_A->de_vec_dt[i] += (-1.0/(Lambda_A))*( e_A_vec_cross_grad_j_A_vec_H[i] \
            + j_A_vec_cross_grad_e_A_vec_H[i] );
        binary_A->dh_vec_dt[i] += -1.0*( j_A_vec_cross_grad_j_A_vec_H[i] \
            + e_A_vec_cross_grad_e_A_vec_H[i] );

        binary_B->de_vec_dt[i] += (-1.0/(Lambda_B))*( e_B_vec_cross_grad_j_B_vec_H[i] \
            + j_B_vec_cross_grad_e_B_vec_H[i] );
        binary_B->dh_vec_dt[i] += -1.0*( j_B_vec_cross_grad_j_B_vec_H[i] \
            + e_B_vec_cross_grad_e_B_vec_H[i] );

        binary_C->de_vec_dt[i] += (-1.0/(Lambda_C))*( e_C_vec_cross_grad_j_C_vec_H[i] \
            + j_C_vec_cross_grad_e_C_vec_H[i] );
        binary_C->dh_vec_dt[i] += -1.0*( j_C_vec_cross_grad_j_C_vec_H[i] \
            + e_C_vec_cross_grad_e_C_vec_H[i] );
    }

    return binary_triplet_hamiltonian;
}



double compute_EOM_pairwise_1PN(ParticlesMap *particlesMap, int binary_index, bool compute_hamiltonian_only)
{
//    printf("compute_EOM_pairwise_1PN\n");
    Particle *binary = (*particlesMap)[binary_index];
    double e = binary->e;
    double a = binary->a;
    double *e_vec_unit = binary->e_vec_unit;
    double *h_vec_unit = binary->h_vec_unit;
    double j = binary->j;
    double j_p2 = binary->j_p2;
    double m1 = binary->child1_mass;
    double m2 = binary->child2_mass;
    double mt = m1+m2;

    double hamiltonian_1PN = -3.0*CONST_G_P2*m1*m2*mt/(a*a*CONST_C_LIGHT_P2*j);
    if (compute_hamiltonian_only == true)
    {
        return hamiltonian_1PN;
    }
    
    double q_vec_unit[3];
    cross3(h_vec_unit,e_vec_unit,q_vec_unit);
    
    double GMdiva = CONST_G*mt/a;
    double Z_1PN = 3.0*sqrt(GMdiva)*GMdiva/(a*CONST_C_LIGHT_P2*j_p2);
    for (int i=0; i<3; i++)
    {
        binary->de_vec_dt[i] += e*Z_1PN*q_vec_unit[i];
    }
    
    return hamiltonian_1PN;
}

double compute_EOM_pairwise_25PN(ParticlesMap *particlesMap, int binary_index, bool compute_hamiltonian_only)
{
    Particle *binary = (*particlesMap)[binary_index];
    double e = binary->e;
    double e_p2 = e*e;
    double a = binary->a;
    double *e_vec_unit = binary->e_vec_unit;
    double *h_vec_unit = binary->h_vec_unit;
    double j = binary->j;
    double j_p2 = binary->j_p2;
    double j_p4 = binary->j_p4;
    double m1 = binary->child1_mass;
    double m2 = binary->child2_mass;
    double mt = m1+m2;

    double a_p3 = a*a*a;
    double GMdiva = CONST_G*mt/a;
    double c_common = CONST_G_P3*m1*m2/(CONST_C_LIGHT_P5*a_p3*j_p4);
    double f_e = 1.0 + c_121div304*e_p2;
    double f_h = 1.0 + c_7div8*e_p2;

    double de_dt = -c_304div15*c_common*mt*e*f_e/(a*j);
    double dh_dt = -c_32div5*c_common*m1*m2*sqrt(GMdiva)*f_h;

    for (int i=0; i<3; i++)
    {
        binary->de_vec_dt[i] += de_dt*e_vec_unit[i];
        binary->dh_vec_dt[i] += dh_dt*h_vec_unit[i];
    }
    
    return 0.0; // N/A
}

double compute_EOM_equilibrium_tide_BO(ParticlesMap *particlesMap, int binary_index, int star_index, int companion_index, int include_tidal_friction_terms, int include_tidal_bulges_precession_terms, int include_rotation_precession_terms, double minimum_eccentricity_for_tidal_precession)
/* Barker & Ogilvie (2009; http://adsabs.harvard.edu/abs/2009MNRAS.395.2268B) */
/* PLUS terms associated with tidal bulges and rotation in only the Z direction */

/* NOTE: in SecularMultiple, the h-vector is defined as the orbital angular momentum vector,
 * NOT the SPECIFIC orbital angular momentum vector. Compared to the notation used by Eggleton,
 * h_vec_SecularMultiple = mu*h_vec_Eggleton where mu = m*M/(m+M) is the reduced mass.
 * In particular, note the line `star->dspin_vec_dt[i] += -dh_vec_dt_star[i]/I;' */
{
//    printf("tides BO\n");
//    printf("TIDES %d %d %d\n",binary_index,star_index,companion_index);
//    printf("minimum_eccentricity_for_tidal_precession %g\n",minimum_eccentricity_for_tidal_precession);
    Particle *binary = (*particlesMap)[binary_index];
    Particle *star = (*particlesMap)[star_index];
    Particle *companion = (*particlesMap)[companion_index];
    
    /* orbit quantities */
    double e = binary->e;
    double e_p2 = binary->e_p2;
    double a = binary->a;
    double h = binary->h;
    double *e_vec = binary->e_vec;
    double *h_vec = binary->h_vec;
    double *e_vec_unit = binary->e_vec_unit;
    double *h_vec_unit = binary->h_vec_unit;
    double j = binary->j;
    double j_p2 = binary->j_p2;
    double j_p3 = binary->j_p3;
    double j_p4 = binary->j_p4;
    double j_p8 = j_p4*j_p4;
    double j_p10 = j_p2*j_p8;  
    double j_p13 = j_p3*j_p10; 
    double j_p4_inv = 1.0/j_p4; 
    double j_p10_inv = 1.0/j_p10; 
    double j_p13_inv = 1.0/j_p13;
    double P_orb = compute_orbital_period(binary);
    double n = 2.0*M_PI/P_orb; /* mean motion */
    

    /* stellar properties */
    double *spin_vec = star->spin_vec;
    double M = star->mass;
    double m = companion->mass;
    double R = star->radius;
    double k_AM = star->tides_apsidal_motion_constant;
    double t_V = star->tides_viscous_time_scale;
    double tau = 3.0*(1.0 + 1.0/(2.0*k_AM))*R*R*R/(CONST_G*M*t_V);
    double rg = star->tides_gyration_radius;
    double I = rg*M*R*R; // moment of intertia
    
    double R_div_a = R/a;
	double R_div_a_p5 = pow(R_div_a,5.0);
    double t_f_inv = 3.0*k_AM*tau*n*n*(m/M)*R_div_a_p5;

    double f_tides1 = f_tides1_function_BO(e_p2,j_p10_inv,j_p13_inv);
    double f_tides2 = f_tides2_function_BO(e_p2,j_p10_inv,j_p13_inv);
    double f_tides3 = f_tides3_function_BO(e_p2,j_p10_inv,j_p13_inv);
    double f_tides4 = f_tides4_function_BO(e_p2,j_p10_inv,j_p13_inv);
    double f_tides5 = f_tides5_function_BO(e_p2,j_p10_inv,j_p13_inv);

    double spin_vec_dot_e_vec = dot3(spin_vec,e_vec);
    double spin_vec_dot_h_vec = dot3(spin_vec,h_vec);

    double q_vec_unit[3];
    cross3(h_vec_unit,e_vec_unit,q_vec_unit);

    double spin_vec_dot_e_vec_unit = dot3(spin_vec,e_vec_unit);
    double spin_vec_dot_h_vec_unit = dot3(spin_vec,h_vec_unit);
    double spin_vec_dot_q_vec_unit = dot3(spin_vec,q_vec_unit);
    double mu = m*M/(m+M);
    double C = m*k_AM*R_div_a_p5/(mu*n);
    
    double Z_rot = 0.0;
    double Z_TB = 0.0;
    if (include_rotation_precession_terms == 1)
    {
        Z_rot = C*c_1div2*j_p4_inv*(2.0*spin_vec_dot_h_vec_unit*spin_vec_dot_h_vec_unit - spin_vec_dot_q_vec_unit*spin_vec_dot_q_vec_unit - spin_vec_dot_e_vec_unit*spin_vec_dot_e_vec_unit);
//        printf("BO Z_rot %g\n",Z_rot);
    }    
    if (include_tidal_bulges_precession_terms == 1)
    {
        Z_TB = C*15.0*n*n*(mu/M)*f_tides2;
//        printf("BO Z_TB %g\n",Z_TB);
    }
    double Z = Z_rot + Z_TB;
        
    double dh_vec_dt_star_i;
   
    for (int i=0; i<3; i++)
    {            
        if (include_tidal_friction_terms == 1)
        {
            dh_vec_dt_star_i = -t_f_inv*( (h/(2.0*n))*(spin_vec_dot_e_vec*f_tides5*e_vec[i] - spin_vec[i]*f_tides3) \
                + h_vec[i]*(f_tides4 - spin_vec_dot_h_vec*f_tides2/(2.0*n*h)) );
            binary->dh_vec_dt[i] += dh_vec_dt_star_i;

            binary->de_vec_dt[i] += -(t_f_inv/h)*( spin_vec_dot_e_vec*f_tides2*h_vec[i]/(2.0*n) \
                + 9.0*e_vec[i]*(f_tides1*h - c_11div18*spin_vec_dot_h_vec*f_tides2/n) );
            star->dspin_vec_dt[i] += -dh_vec_dt_star_i/I;
        }
        if (include_rotation_precession_terms == 1 || include_tidal_bulges_precession_terms == 1)
        {
            if (e >= minimum_eccentricity_for_tidal_precession)
            {
                binary->de_vec_dt[i] += Z*e*q_vec_unit[i];
            }
        }
    }

    return 0;
}

double compute_EOM_equilibrium_tide_BO_full(ParticlesMap *particlesMap, int binary_index, int star_index, int companion_index, int include_tidal_friction_terms, int include_tidal_bulges_precession_terms, int include_rotation_precession_terms, double minimum_eccentricity_for_tidal_precession)
/* Barker & Ogilvie (2009; http://adsabs.harvard.edu/abs/2009MNRAS.395.2268B) */
/* PLUS terms associated with tidal bulges and rotation in the X, Y & Z directions */

/* NOTE: in SecularMultiple, the h-vector is defined as the orbital angular momentum vector,
 * NOT the SPECIFIC orbital angular momentum vector. Compared to the notation used by Eggleton,
 * h_vec_SecularMultiple = mu*h_vec_Eggleton where mu = m*M/(m+M) is the reduced mass.
 * In particular, note the line `star->dspin_vec_dt[i] += -dh_vec_dt_star[i]/I;' */
{
//    printf("tides BO full\n");
//    printf("TIDES %d %d %d\n",binary_index,star_index,companion_index);
    Particle *binary = (*particlesMap)[binary_index];
    Particle *star = (*particlesMap)[star_index];
    Particle *companion = (*particlesMap)[companion_index];
    
    /* orbit quantities */
    double e = binary->e;
    double e_p2 = binary->e_p2;
    double a = binary->a;
    double h = binary->h;
    double *e_vec = binary->e_vec;
    double *h_vec = binary->h_vec;
    double *e_vec_unit = binary->e_vec_unit;
    double *h_vec_unit = binary->h_vec_unit;
    double j = binary->j;
    double j_p2 = binary->j_p2;
    double j_p3 = binary->j_p3;
    double j_p4 = binary->j_p4;
    double j_p8 = j_p4*j_p4;
    double j_p10 = j_p2*j_p8;  
    double j_p13 = j_p3*j_p10; 
    double j_p4_inv = 1.0/j_p4; 
    double j_p10_inv = 1.0/j_p10; 
    double j_p13_inv = 1.0/j_p13;
    double P_orb = compute_orbital_period(binary);
    double n = 2.0*M_PI/P_orb; /* mean motion */
    

    /* stellar properties */
    double *spin_vec = star->spin_vec;
    double M = star->mass;
    double m = companion->mass;
    double R = star->radius;

    double k_AM = star->tides_apsidal_motion_constant;
    double t_V = star->tides_viscous_time_scale;
    double tau = 3.0*(1.0 + 1.0/(2.0*k_AM))*R*R*R/(CONST_G*M*t_V);
    double rg = star->tides_gyration_radius;
    double I = rg*M*R*R; // moment of intertia
    
    double R_div_a = R/a;
	double R_div_a_p5 = pow(R_div_a,5.0);
    double t_f_inv = 3.0*k_AM*tau*n*n*(m/M)*R_div_a_p5;

    double f_tides1 = f_tides1_function_BO(e_p2,j_p10_inv,j_p13_inv);
    double f_tides2 = f_tides2_function_BO(e_p2,j_p10_inv,j_p13_inv);
    double f_tides3 = f_tides3_function_BO(e_p2,j_p10_inv,j_p13_inv);
    double f_tides4 = f_tides4_function_BO(e_p2,j_p10_inv,j_p13_inv);
    double f_tides5 = f_tides5_function_BO(e_p2,j_p10_inv,j_p13_inv);

    double spin_vec_dot_e_vec = dot3(spin_vec,e_vec);
    double spin_vec_dot_h_vec = dot3(spin_vec,h_vec);

    double q_vec_unit[3];
    cross3(h_vec_unit,e_vec_unit,q_vec_unit);

    double spin_vec_dot_e_vec_unit = dot3(spin_vec,e_vec_unit);
    double spin_vec_dot_h_vec_unit = dot3(spin_vec,h_vec_unit);
    double spin_vec_dot_q_vec_unit = dot3(spin_vec,q_vec_unit);
    double mu = m*M/(m+M);
    double C = m*k_AM*R_div_a_p5/(mu*n);
    
    double C_rot;
    double X_rot = 0.0;
    double Y_rot = 0.0;
    double Z_rot = 0.0;
    double Z_TB = 0.0;
    if (include_rotation_precession_terms == 1)
    {
        C_rot = C*j_p4_inv*spin_vec_dot_h_vec_unit;
        X_rot = -C_rot*spin_vec_dot_e_vec_unit;
        Y_rot = -C_rot*spin_vec_dot_q_vec_unit;
        Z_rot = C*c_1div2*j_p4_inv*(2.0*spin_vec_dot_h_vec_unit*spin_vec_dot_h_vec_unit - spin_vec_dot_q_vec_unit*spin_vec_dot_q_vec_unit - spin_vec_dot_e_vec_unit*spin_vec_dot_e_vec_unit);
    }    
    if (include_tidal_bulges_precession_terms == 1)
    {
        Z_TB = C*15.0*n*n*(mu/M)*f_tides2;
    }
    
    double X = X_rot;
    double Y = Y_rot;
    double Z = Z_rot + Z_TB;
    

    double dh_vec_dt_star_i;
   
    for (int i=0; i<3; i++)
    {
            
        if (include_tidal_friction_terms == 1)
        {
            dh_vec_dt_star_i = -t_f_inv*( (h/(2.0*n))*(spin_vec_dot_e_vec*f_tides5*e_vec[i] - spin_vec[i]*f_tides3) \
                + h_vec[i]*(f_tides4 - spin_vec_dot_h_vec*f_tides2/(2.0*n*h)) );
            binary->dh_vec_dt[i] += dh_vec_dt_star_i;

            binary->de_vec_dt[i] += -(t_f_inv/h)*( spin_vec_dot_e_vec*f_tides2*h_vec[i]/(2.0*n) \
                + 9.0*e_vec[i]*(f_tides1*h - c_11div18*spin_vec_dot_h_vec*f_tides2/n) );
            star->dspin_vec_dt[i] += -dh_vec_dt_star_i/I;
        }
        if (include_rotation_precession_terms == 1 || include_tidal_bulges_precession_terms == 1)
        {
            if (e >= minimum_eccentricity_for_tidal_precession)
            {
                binary->de_vec_dt[i] += e*(Z*q_vec_unit[i] - Y*h_vec_unit[i]);
                
                dh_vec_dt_star_i = h*(-X*q_vec_unit[i] + Y*e_vec_unit[i]);
                binary->dh_vec_dt[i] += dh_vec_dt_star_i;
                star->dspin_vec_dt[i] += -dh_vec_dt_star_i/I;
                
            }
        }
    }

    return 0;
}


double f_tides1_function_BO(double e_p2, double j_p10_inv, double j_p13_inv)
{
    return j_p13_inv*(1.0 + e_p2*(c_15div4 + e_p2*(c_15div8 + e_p2*c_5div64)));
}
double f_tides2_function_BO(double e_p2, double j_p10_inv, double j_p13_inv)
{
    return j_p10_inv*(1.0 + e_p2*(c_3div2 + e_p2*c_1div8));
}
double f_tides3_function_BO(double e_p2, double j_p10_inv, double j_p13_inv)
{
    return j_p10_inv*(1.0 + e_p2*(c_9div2 + e_p2*c_5div8));
}
double f_tides4_function_BO(double e_p2, double j_p10_inv, double j_p13_inv)
{
    return j_p13_inv*(1.0 + e_p2*(c_15div2 + e_p2*(c_45div8 + e_p2*c_5div16)));
}
double f_tides5_function_BO(double e_p2, double j_p10_inv, double j_p13_inv)
{
    return j_p10_inv*(3.0 + c_1div2*e_p2);
}  


double compute_EOM_equilibrium_tide(ParticlesMap *particlesMap, int binary_index, int star_index, int companion_index, int include_tidal_friction_terms, int include_tidal_bulges_precession_terms, int include_rotation_precession_terms, double minimum_eccentricity_for_tidal_precession)

/* Equilibrium tide in vector form adopted from Eggleton & Kisseleva 1998 */
 
/* NOTE: in SecularMultiple, the h-vector is defined as the orbital angular momentum vector,
 * NOT the SPECIFIC orbital angular momentum vector. Compared to the notation used by Eggleton,
 * h_vec_SecularMultiple = mu*h_vec_Eggleton where mu = m*M/(m+M) is the reduced mass.
 * In particular, note the line `star->dspin_vec_dt[i] += -dh_vec_dt_star[i]/I;' */

{
//    printf("tides EK \n");
//    printf("TIDES %d %d %d\n",binary_index,star_index,companion_index);
    Particle *binary = (*particlesMap)[binary_index];
    Particle *star = (*particlesMap)[star_index];
    Particle *companion = (*particlesMap)[companion_index];
    
    /* orbit quantities */
    double e = binary->e;
    double e_p2 = binary->e_p2;
    double a = binary->a;
    double h = binary->h;
    double *e_vec = binary->e_vec;
    double *h_vec = binary->h_vec;
    double *e_vec_unit = binary->e_vec_unit;
    double *h_vec_unit = binary->h_vec_unit;
    double j = binary->j;
    double j_p2 = binary->j_p2;
    double j_p3 = binary->j_p3;
    double j_p4 = binary->j_p4;
    double j_p8 = j_p4*j_p4;
    double j_p10 = j_p2*j_p8;  
    double j_p13 = j_p3*j_p10; 
    double j_p4_inv = 1.0/j_p4;
    double j_p10_inv = 1.0/j_p10; 
    double j_p13_inv = 1.0/j_p13;
    double P_orb = compute_orbital_period(binary);
    double n = 2.0*M_PI/P_orb; /* mean motion */
    
    /* stellar properties */
    double *spin_vec = star->spin_vec;
    double M = star->mass;
    double m = companion->mass;
    double mu = m*M/(m+M);
    double R = star->radius;
    double k_AM = star->tides_apsidal_motion_constant;
    double t_V = star->tides_viscous_time_scale;
    double rg = star->tides_gyration_radius;
    double I = rg*M*R*R; // moment of intertia
    
    double R_div_a = R/a;
	double R_div_a_p5 = pow(R_div_a,5.0);
    double R_div_a_p8 = pow(R_div_a,8.0);
    double t_f_inv = (9.0/t_V)*R_div_a_p8*((M+m)*m/(M*M))*(1.0 + 2.0*k_AM)*(1.0 + 2.0*k_AM);
    
    double q_vec_unit[3];
    cross3(h_vec_unit,e_vec_unit,q_vec_unit);

    double spin_vec_dot_e_vec_unit = dot3(spin_vec,e_vec_unit);
    double spin_vec_dot_h_vec_unit = dot3(spin_vec,h_vec_unit);
    double spin_vec_dot_q_vec_unit = dot3(spin_vec,q_vec_unit);
        
    double V,W,X,Y,Z;
    VWXYZ_tides_function(
        include_tidal_friction_terms, include_tidal_bulges_precession_terms, include_rotation_precession_terms, minimum_eccentricity_for_tidal_precession, \
        t_f_inv,k_AM, \
        m, M, mu, n, R_div_a_p5, \
        e, e_p2, j_p4_inv, j_p10_inv, j_p13_inv, \
        spin_vec_dot_h_vec_unit, spin_vec_dot_e_vec_unit, spin_vec_dot_q_vec_unit, \
        &V, &W, &X, &Y, &Z);

    double dh_vec_dt_star[3];
    
    for (int i=0; i<3; i++)
    {
        dh_vec_dt_star[i] = h*( Y*e_vec_unit[i] - W*h_vec_unit[i] - X*q_vec_unit[i] );
        binary->dh_vec_dt[i] += dh_vec_dt_star[i];

        binary->de_vec_dt[i] += e*( -V*e_vec_unit[i] - Y*h_vec_unit[i] + Z*q_vec_unit[i] );
        
        star->dspin_vec_dt[i] += -dh_vec_dt_star[i]/I; /* conservation of total angular momentum (orbit+spin) */
        
    }

    return 0;
}

double f_tides1_function(double e_p2, double j_p10_inv, double j_p13_inv)
{
    return j_p13_inv*(1.0 + e_p2*(c_15div4 + e_p2*(c_15div8 + e_p2*c_5div64)));
}
double f_tides2_function(double e_p2, double j_p10_inv, double j_p13_inv)
{
    return j_p10_inv*(1.0 + e_p2*(c_3div2 + e_p2*c_1div8));
}
double f_tides3_function(double e_p2, double j_p10_inv, double j_p13_inv)
{
    return j_p13_inv*(1.0 + e_p2*(c_15div2 + e_p2*(c_45div8 + e_p2*c_5div16)));
}
double f_tides4_function(double e_p2, double j_p10_inv, double j_p13_inv)
{
    return j_p10_inv*(1.0 + e_p2*(3.0 + e_p2*c_3div8));
}
double f_tides5_function(double e_p2, double j_p10_inv, double j_p13_inv)
{
    return j_p10_inv*(1.0 + e_p2*(c_9div2 + e_p2*c_5div8));
}


int VWXYZ_tides_function
(
    int include_tidal_friction_terms, int include_tidal_bulges_precession_terms, int include_rotation_precession_terms, double minimum_eccentricity_for_tidal_precession, \
    double t_f_inv, double k_AM, \
    double m, double M, double mu, double n, double R_div_a_p5, \
    double e, double e_p2, double j_p4_inv, double j_p10_inv,double j_p13_inv, \
    double spin_vec_dot_h_vec_unit, double spin_vec_dot_e_vec_unit,double spin_vec_dot_q_vec_unit, \
    double* V, double* W, double* X, double* Y, double* Z
)
{
   
    *V = 0.0;
    *W = 0.0;
    *X = 0.0;
    *Y = 0.0;
    *Z = 0.0;

    double f2 = f_tides2_function(e_p2,j_p10_inv,j_p13_inv); /* needed for both pure tidal dissipation and pure tidal bulges terms */
    
    if (include_tidal_friction_terms == 1)
    {
        double f1 = f_tides1_function(e_p2,j_p10_inv,j_p13_inv);
        double f3 = f_tides3_function(e_p2,j_p10_inv,j_p13_inv);
        double f4 = f_tides4_function(e_p2,j_p10_inv,j_p13_inv);
        double f5 = f_tides5_function(e_p2,j_p10_inv,j_p13_inv);
    
        *V += 9.0*t_f_inv*(f1    - c_11div18*(spin_vec_dot_h_vec_unit/n)*f2);
        *W += t_f_inv*(f3        - (spin_vec_dot_h_vec_unit/n)*f4);
        *X +=  -t_f_inv*spin_vec_dot_q_vec_unit*f5/(2.0*n);
        *Y += t_f_inv*spin_vec_dot_e_vec_unit*f2/(2.0*n);
    }

    if (e < minimum_eccentricity_for_tidal_precession)
    {
        return 0;
    }


    if ((include_tidal_bulges_precession_terms == 1) || (include_rotation_precession_terms) == 1)
    {
        double C = m*k_AM*R_div_a_p5/(mu*n);
    
        if (include_tidal_bulges_precession_terms == 1)
        {
            *Z += C*15.0*n*n*(mu/M)*f2;
//            printf("include_tidal_bulges_precession_terms Z %g\n",*Z);
        }    
        if (include_rotation_precession_terms == 1)
        {
            double C_XY = -C*spin_vec_dot_h_vec_unit*j_p4_inv;
            *X += C_XY*spin_vec_dot_e_vec_unit;
            *Y += C_XY*spin_vec_dot_q_vec_unit;            
            *Z += C*c_1div2*j_p4_inv*(2.0*spin_vec_dot_h_vec_unit*spin_vec_dot_h_vec_unit - spin_vec_dot_q_vec_unit*spin_vec_dot_q_vec_unit - spin_vec_dot_e_vec_unit*spin_vec_dot_e_vec_unit);
//            printf("include_rotation_precession_terms Z %g\n",*Z);
        }
    }
    return 0;
}

double compute_orbital_period(Particle *particle)
{
	double a = particle->a;
	double total_mass = particle->child1_mass_plus_child2_mass;
	return 2.0*M_PI*sqrt(a*a*a/(CONST_G*total_mass));
}

void extract_ODE_variables(ParticlesMap *particlesMap, N_Vector &y, bool reset_ODE_quantities)
{
    ParticlesMapIterator it_p;
    int k=1;
    int k_component;

    double spin_vec[3],e_vec[3],h_vec[3];
    for (it_p = particlesMap->begin(); it_p != particlesMap->end(); it_p++)
    {
        Particle *P_p = (*it_p).second;
        if (P_p->is_binary == 0) // particle is a body
        {
            for (k_component=0; k_component<3; k_component++)
            {
                P_p->spin_vec[k_component] = Ith(y,k + k_component);
            }
            
            k=k+3;
        }
        if (P_p->is_binary == 1) // particle is a binary
        {
            for (k_component=0; k_component<3; k_component++)
            {
                P_p->e_vec[k_component] = Ith(y,k + k_component);
                P_p->h_vec[k_component] = Ith(y,k + k_component + 3);
            }
            
            k=k+6;
        }
        P_p->set_ODE_quantities();
        
        if (reset_ODE_quantities == true)
        {
            P_p->reset_ODE_quantities();
        }
        
    }
}

void write_ODE_variables_dots(ParticlesMap *particlesMap, N_Vector &y_dot)
{
    ParticlesMapIterator it_p;
    int k=1;
    int k_component;

    double spin_vec[3],e_vec[3],h_vec[3];
    for (it_p = particlesMap->begin(); it_p != particlesMap->end(); it_p++)
    {
        Particle *P_p = (*it_p).second;
        if (P_p->is_binary == 0) // particle is a body
        {
            for (k_component=0; k_component<3; k_component++)
            {
                Ith(y_dot,k + k_component) = P_p->dspin_vec_dt[k_component];
            }
            
            k=k+3;
        }
        if (P_p->is_binary == 1) // particle is a binary
        {
            for (k_component=0; k_component<3; k_component++)
            {
                Ith(y_dot,k + k_component)      = P_p->de_vec_dt[k_component];
                Ith(y_dot,k + k_component + 3)  = P_p->dh_vec_dt[k_component];
            }
            
            k=k+6;
        }
    }
}












int root_finding_functions(realtype t, N_Vector y, realtype *root_functions, void *data_)
{

	UserData data;
	data = (UserData) data_;
    ParticlesMap *particlesMap = data->particlesMap;
    int N_root_finding = data->N_root_finding;
    
    extract_ODE_variables(particlesMap, y, false);

    double large_quantity = 1.0e10;
    for (int i=0; i<N_root_finding; i++)
    {
        root_functions[i] = large_quantity;
    }

    ParticlesMapIterator it_p;
    std::vector<int>::iterator it_parent;
    
    int i_root = 0;
    
    for (it_p = particlesMap->begin(); it_p != particlesMap->end(); it_p++)
    {
        Particle *P_p = (*it_p).second;
        if (P_p->is_binary == 1)
        {
            if (P_p->check_for_secular_breakdown == 1)
            {
                double hamiltonian=0.0;
                for (it_parent = P_p->parents.begin(); it_parent != P_p->parents.end(); it_parent++)
                {
                    int i = std::distance(P_p->parents.begin(), it_parent);
                    Particle *P_q = (*particlesMap)[(*it_parent)];
                    int connecting_child_in_parent = P_p->connecting_child_in_parents[i];
                    hamiltonian += compute_EOM_binary_pairs(particlesMap,P_p->index,P_q->index,connecting_child_in_parent,false);
                }

                double AM_time_scale = compute_AM_time_scale(P_p);
                double orbital_period = compute_orbital_period(P_p);
                root_functions[i_root] = 1.0 - AM_time_scale/orbital_period;

                i_root++;
            }
            
            if (P_p->check_for_dynamical_instability == 1)
            {
                if (P_p->parent != -1)
                {
                    Particle *P_parent = (*particlesMap)[P_p->parent];
                    Particle *P_child1 = (*particlesMap)[P_p->child1];
                    Particle *P_child2 = (*particlesMap)[P_p->child2];
                
                    double a_out = P_parent->a;
                    double e_out = P_parent->e;
                    double a_in = P_p->a;
                    double e_in = P_p->e;
                    double M_p = P_p->mass;
                    double ra_in = a_in*(1.0+e_in);
                    double rp_out = a_out*(1.0-e_out);
                    

                    Particle *P_sibling = (*particlesMap)[P_p->sibling];
        
                    if (P_p->dynamical_instability_criterion == 0) /* for mass ratios on the order of unity */
                    {
                        /* Mardling & Aarseth 2001 */
                        double q_out = P_sibling->mass/M_p;
                        double rel_INCL = 0.0;
                        get_inclination_relative_to_parent(P_p->index,&rel_INCL);
                        root_functions[i_root] = rp_out - a_in*2.8*pow( (1.0+q_out)*(1.0+e_out)/sqrt(1.0-e_out),2.0/5.0)*(1.0 - 0.3*rel_INCL/M_PI);
                    }
                    else if (P_p->dynamical_instability_criterion > 0)
                    /* in case of a central dominant particle
                     * m1 is the `inner' mass; m2 is the `outer' mass */
                    {
                        int central_particle_index = P_p->dynamical_instability_central_particle;
                        Particle *P_central_particle = (*particlesMap)[central_particle_index];
                        double central_particle_mass = P_central_particle->mass;
                        double m2 = P_sibling->mass;
                        double m1;
                        if (P_p->child1 == central_particle_index)
                        {
                            m1 = P_child2->mass;
                        }
                        else if (P_p->child2 == central_particle_index)
                        {
                            m1 = P_child1->mass;
                        }
                        
                        int central_particle_parent;
                        std::vector<int>::iterator it_C_parent;
                        for (it_C_parent = P_central_particle->parents.begin(); it_C_parent != P_central_particle->parents.end(); it_C_parent++)
                        {
                            central_particle_parent = *it_C_parent;
                            if (P_p->child1 == central_particle_parent)
                            {
                                m1 = P_child2->mass;
                            }
                            else if (P_p->child2 == central_particle_parent)
                            {
                                m1 = P_child1->mass;
                            }
                        }

                        double mu1 = m1/central_particle_mass;
                        double mu2 = m2/central_particle_mass;
                        double R_H = c_1div2*(a_in+a_out)*pow( c_1div3*(mu1+mu2), c_1div3 );
                        
                        if (P_p->dynamical_instability_criterion == 1)
                        {
                            double K = P_p->dynamical_instability_K_parameter;
                            root_functions[i_root] = (rp_out - ra_in) - K*R_H;
                        }

                        else if (P_p->dynamical_instability_criterion == 2)
                        /* Petrovich 2015 */
                        {
                            root_functions[i_root] = rp_out/ra_in - ( 2.4*pow( max(mu1,mu2), c_1div3)*sqrt(a_out/a_in) + 1.15 );
                        }
                        
                    }

                    if (root_functions[i_root] <= 0.0)
                    {
                        P_p->dynamical_instability_has_occurred = 1;
                    }


                }
                i_root++;                
            }
            if (P_p->check_for_physical_collision_or_orbit_crossing == 1)
            {
                Particle *P_child1 = (*particlesMap)[P_p->child1];
                Particle *P_child2 = (*particlesMap)[P_p->child2];
                
                double cross_section = 0.0;
                cross_section_function(P_child1,&cross_section);
                cross_section_function(P_child2,&cross_section);

                double periapse_distance = P_p->a*(1.0 - P_p->e);
                root_functions[i_root] = 1.0 - periapse_distance/cross_section;
                
//                printf("root finding check_for_physical_collision_or_orbit_crossing %d %g %g %g\n",P_p->index,P_p->a,cross_section, root_functions[i_root]);
                if (root_functions[i_root] >= 0.0)
                {
                    P_p->physical_collision_or_orbit_crossing_has_occurred = 1;
//                    printf("root finding check_for_physical_collision_or_orbit_crossing %d %g %g %g\n",P_p->index,P_p->a,cross_section, root_functions[i_root]);
                }                
                
                i_root++;
            }
            if (P_p->check_for_minimum_periapse_distance == 1)
            {
                Particle *P_child1 = (*particlesMap)[P_p->child1];
                Particle *P_child2 = (*particlesMap)[P_p->child2];
                
                double cross_section = P_p->check_for_minimum_periapse_distance_value;
                double periapse_distance = P_p->a*(1.0 - P_p->e);
                root_functions[i_root] = 1.0 - periapse_distance/cross_section;
                
//                printf("root finding check_for_minimum_periapse_distance %d %g %g %g\n",P_p->index,P_p->a,cross_section, root_functions[i_root]);
                if (root_functions[i_root] >= 0.0)
                {
                    P_p->minimum_periapse_distance_has_occurred = 1;
//                    printf("root finding check_for_minimum_periapse_distance %d %g %g %g\n",P_p->index,P_p->a,cross_section, root_functions[i_root]);
                }                

                i_root++;
            }
        }
        else /* P_p not a binary */
        {

            if (P_p->check_for_RLOF_at_pericentre == 1)
            {
                if (P_p->parent != -1)
                {
                    Particle *P_parent = (*particlesMap)[P_p->parent];
                    Particle *P_sibling = (*particlesMap)[P_p->sibling];
                    
                    double a = P_parent->a;
                    double e = P_parent->e;
                    double rp = a*(1.0 - e);
                    double subject_mass = P_p->mass;
                    double companion_mass = P_sibling->mass;
    
                    double spin_angular_frequency = P_p->spin_vec_norm;
                    double orbital_angular_frequency_periapse = sqrt( CONST_G*(subject_mass + companion_mass)*(1.0 + e)/(rp*rp*rp) );
                    double f = spin_angular_frequency/orbital_angular_frequency_periapse;
            
                    double roche_radius_pericenter;
                    if (P_p->check_for_RLOF_at_pericentre_use_sepinsky_fit == 0)
                    {
                        roche_radius_pericenter = roche_radius_pericenter_eggleton(rp, subject_mass/companion_mass);
                    }
                    else
                    {
                        roche_radius_pericenter = roche_radius_pericenter_sepinsky(rp, subject_mass/companion_mass, e, f);
                    }
                    
                    root_functions[i_root] = 1.0 - P_p->radius/roche_radius_pericenter;
//                    printf("check_for_RLOF_at_pericentre rp %g roche_radius_pericenter %g R %g\n", rp, roche_radius_pericenter, P_p->radius);
                    
                    if (root_functions[i_root] <= 0.0)
                    {
                        P_p->RLOF_at_pericentre_has_occurred = 1;
                    }
                }
                
                i_root++;
            }
        }
    }
    return 0;
}

void cross_section_function(Particle *p, double *cross_section)
{ 
    if (p->is_binary==0)
    {
        *cross_section += p->radius;
    }
    else
    {
        *cross_section += p->a*(1.0 + p->e);
    }
}
double compute_AM_time_scale(Particle *P_p)
{
    double e = P_p->e;
    double e_p2 = P_p->e_p2;
    double de_dt = dot3(P_p->e_vec_unit,P_p->de_vec_dt);    
    double AM_time_scale = ((1.0-e_p2)/(e*fabs(de_dt)));
    
    return AM_time_scale;
}
         

double roche_radius_pericenter_eggleton(double rp, double q)
{
    /* 2007ApJ...660.1624S Eqs. (45) */    
    /* q is defined as m_primary/m_secondary */
    double q_pow_one_third = pow(q,c_1div3);
    double q_pow_two_third = q_pow_one_third*q_pow_one_third;
    return rp*0.49*q_pow_two_third/(0.6*q_pow_two_third + log(1.0 + q_pow_one_third));
}
double roche_radius_pericenter_sepinsky(double rp, double q, double e, double f)
{
    /* 2007ApJ...660.1624S Eqs. (47)-(52) */
    double log_q = log10(q);
    double A = f*f*(1.0 + e); // assumes pericenter
    double log_A = log10(A);

    double R_L_pericenter_eggleton = roche_radius_pericenter_eggleton(rp,q);
    double ratio = 0.0; // this is R_L divided by R_L_pericenter_eggleton

    if (log_q < 0.0)
    {
        if (log_A <= -0.1)
        {
            double c = 0.5*(1.0+A) + log_q;
            ratio = 1.0 + 0.11*(1.0-A) - 0.05*(1.0-A)*exp(-c*c);
        }
        if ((log_A > -0.1) && (log_A < 0.2))
        {
            double g_0 = 0.9978 - 0.1229*log_A - 0.1273*log_A*log_A;
            double g_1 = 0.001 + 0.02556*log_A;
            double g_2 = 0.0004 + 0.0021*log_A;
            ratio = g_0 + g_1*log_q * g_2*log_q*log_q;
        }
        if (log_A >= 0.2)
        {
            double num_0 = 6.3014*pow(log_A,1.3643);
            double den_0 = exp(2.3644*pow(log_A,0.70748)) - 1.4413*exp(-0.0000184*pow(log_A,-4.5693));
            double i_0 = num_0/den_0;

            double den_1 = 0.0015*exp(8.84*pow(log_A,0.282)) + 15.78;
            double i_1 = log_A/den_1;

            double num_2 = 1.0 + 0.036*exp(8.01*pow(log_A,0.879));
            double den_2 = 0.105*exp(7.91*pow(log_A,0.879));
            double i_2 = num_2/den_2;

            double den_3 = 1.38*exp(-0.035*pow(log_A,0.76)) + 23.0*exp(-2.89*pow(log_A,0.76));
            double i_3 = 0.991/den_3;

            double c = log_q + i_3;
            ratio = i_0 + i_1*exp(-i_2*c*c);
        }
    }
    if (log_q >= 0.0)
    {
        if (log_A <= -0.1)
        {
            ratio = 1.226 - 0.21*A - 0.15*(1.0-A)*exp( (0.25*A - 0.3)*pow(log_q,1.55) );
        }
        if ((log_A > -0.1) && (log_A < 0.2))
        {
            double log_A_p2 = log_A*log_A;
            double h_0 = 1.0071 - 0.0907*log_A - 0.0495*log_A_p2;
            double h_1 = -0.004 - 0.163*log_A - 0.214*log_A_p2;
            double h_2 = 0.00022 - 0.0108*log_A - 0.02718*log_A_p2;
            ratio = h_0 + h_1*log_q + h_2*log_q*log_q;
        }
        if (log_A >= 0.2)
        {
            double num_0 = 1.895*pow(log_A,0.837);
            double den_0 = exp(1.636*pow(log_A,0.789)) - 1.0;
            double j_0 = num_0/den_0;

            double num_1 = 4.3*pow(log_A,0.98);
            double den_1 = exp(2.5*pow(log_A,0.66)) + 4.7;
            double j_1 = num_1/den_1;

            double den_2 = 8.8*exp(-2.95*pow(log_A,0.76)) + 1.64*exp(-0.03*pow(log_A,0.76));
            double j_2 = 1.0/den_2;

//            double j_3 = 0.256*exp(-1.33*pow(log_A,2.9))*( 5.5*exp(1.33*pow(log_A,2.9)) + 1.0 );
            double j_3 = 0.256*(5.5 + exp(-1.33*pow(log_A,2.9)));

            ratio = j_0 + j_1*exp(-j_2*pow(log_q,j_3));
            
//            printf("log_A %g\n",log_A);
//            printf("1 %g %g %g \n",num_0,den_0,j_0);
//            printf("2 %g %g %g \n",num_1,den_1,j_1);            
//            printf("2 %g %g %g \n",den_2,j_2,j_3);            
//            printf("ratio %g %g %g \n",ratio);            
        }
    }

    if (ratio == 0.0)
    {
        printf("unrecoverable error occurred in function roche_radius_pericenter_sepinsky in ODE_system.cpp\n");
        printf("log_q %g log_A %g ratio %g\n",log_q,log_A,ratio);
        printf("rp %g q %g e %g f %g\n",rp,q,e,f);
        exit(-1);
    }
    
    return ratio*R_L_pericenter_eggleton;
}

