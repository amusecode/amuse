/*
*/

#include "types.h"
#include "newtonian.h"
#include "../interface.h" /* for parameters */
#include <stdio.h>

double compute_orbital_period(Particle *particle)
{
	double a = particle->a;
	double total_mass = particle->child1_mass_plus_child2_mass;
	return 2.0*M_PI*sqrt(a*a*a/(CONST_G*total_mass));
}

double compute_EOM_binary_pairs(ParticlesMap *particlesMap, int inner_binary_index, int outer_binary_index, int connecting_child_in_outer_binary, bool compute_hamiltonian_only)
{
    /* last checked 23-06-15 */


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
    
    /* set alpha = +1 */
    double m1 = P_child1->mass;
    double m2 = P_child2->mass;
    double m3 = P_sibling->mass;
//    printf("m1 %g m2 %g m3 %g\n",m1,m2,m3);
//    printf("a_in %g a_out %g\n",a_in,a_out);
        
    double m1_plus_m2 = inner_binary->child1_mass_plus_child2_mass;
    double m1_minus_m2 = inner_binary->child1_mass_minus_child2_mass;
    double m1_times_m2 = inner_binary->child1_mass_times_child2_mass;

    double A_quad = c_1div8*CONST_G*(a_in*a_in/(a_out*a_out*a_out))*m1_times_m2*m3/m1_plus_m2;
    //double A_oct = A_quad*c_15div8*(a_in/a_out)*m1_minus_m2/m1_plus_m2;
    double A_oct = A_quad*c_15div8*(a_in/a_out)*fabs(m1_minus_m2)/m1_plus_m2; /* changed 06-08-15 */
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
//        A_hd = c_3div1024*CONST_G*(a_in*a_in*a_in*a_in/(a_out*a_out*a_out*a_out*a_out))*(m1_times_m2*m3*(m1*m1*m1 + m2*m2*m2)/(m1_plus_m2*m1_plus_m2*m1_plus_m2*m1_plus_m2));
        /* the above two expressions are mathematically identical */
//        printf("hd true %g %g %g\n",A_quad,A_oct,A_hd);
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


#ifdef IGNORE
    if (include_octupole_order_binary_pair_terms == true)
    {
        printf("include_octupole_order_binary_pair_terms\n");
    }
    if (include_hexadecupole_order_binary_pair_terms == true)
    {
        printf("include_hexadecupole_order_binary_pair_terms\n");
    }
    if (include_dotriacontupole_order_binary_pair_terms == true)
    {
        printf("include_dotriacontupole_order_binary_pair_terms\n");
    }
#endif

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
            
        //printf("testtttt %g %g \n",inner_binary->de_vec_dt[i],inner_binary->dh_vec_dt[i]);
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
    /* last checked 23-06-15 */
    
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
   
   /* set alpha = +1 */
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
//    double j_B_p2 = binary_B->j_p2;
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
//    cross_term_hamiltonian *= -1.0;
    
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

