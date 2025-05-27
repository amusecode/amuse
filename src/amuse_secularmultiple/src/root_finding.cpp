/*
*/

#include "types.h"
#include "evolve.h"
#include "root_finding.h"
#include <stdio.h>


int root_finding_functions(realtype time, N_Vector y, realtype *root_functions, void *data_)
{

	UserData data;
	data = (UserData) data_;
    ParticlesMap *particlesMap = data->particlesMap;
    int N_root_finding = data->N_root_finding;
    double start_time = data->start_time;
    double delta_time = time - start_time;
    
    extract_ODE_variables(particlesMap, y, delta_time, false); // do not reset ODE quantities

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
                //printf("sb %g\n",root_functions[i_root]);

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
                    
//                    printf("check_for_dynamical_instability p %d par %d c1 %d c2 %d \n",P_p->index,P_parent->index,P_child1->index,P_child2->index);
//                    printf("a_in %g a_out %g\n",a_in,a_out);
    
                    Particle *P_sibling = (*particlesMap)[P_p->sibling];
        
                    if (P_p->dynamical_instability_criterion == 0) /* for mass ratios on the order of unity */
                    {
                        /* Mardling & Aarseth 2001 */
                        double q_out = P_sibling->mass/M_p;
                        double rel_INCL = 0.0;
                        get_inclination_relative_to_parent(P_p->index,&rel_INCL);
                        root_functions[i_root] = rp_out - a_in*2.8*pow( (1.0+q_out)*(1.0+e_out)/sqrt(1.0-e_out),2.0/5.0)*(1.0 - 0.3*rel_INCL/M_PI);
                        //printf("di MA %g\n",root_functions[i_root]);
                        //printf("di MA %g %g\n",P_sibling->mass,M_p);
                    }
                    if (P_p->dynamical_instability_criterion == 1) /* for S-type test particles in binaries */
                    {
                        /* Wiegert & Holman 1999 */
                        double m1 = M_p;
                        double m2 = P_sibling->mass;
                        double mu = m2/(m1+m2);
                        double e = e_out;
                        root_functions[i_root] = (0.464 - 0.38*mu - 0.631*e + 0.586*mu*e + 0.15*e*e - 0.198*mu*e*e) - a_in/a_out;
                        //printf("di WH %g %g %g %g %g %g\n",root_functions[i_root],m1,m2,e,a_in/a_out,(0.464 - 0.38*mu - 0.631*e + 0.586*mu*e + 0.15*e*e - 0.198*mu*e*e));
                        
                           
                        //printf("di MA %g %g\n",P_sibling->mass,M_p);
                    }
                    else if (P_p->dynamical_instability_criterion > 1)
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
                            
//                        else /* the central particle is not a direct child of P_p */
//                        {
//                            if (P_child1->is_binary == true)
//                            {
//                                m1 = P_child2->mass;
//                            }
//                            else if (P_child2->is_binary == true)
//                            {
//                                m1 = P_child1->mass;
//                            }
//                            else
//                            {
//                                printf("error in root finding function dynamical_stability_criterion > 0: the system should be `fully nested'; exiting\n");
//                                exit(-1);
//                            }

                        double mu1 = m1/central_particle_mass;
                        double mu2 = m2/central_particle_mass;
                        double R_H = c_1div2*(a_in+a_out)*pow( c_1div3*(mu1+mu2), c_1div3 );
                        
                        if (P_p->dynamical_instability_criterion == 2)
                        {
                            double K = P_p->dynamical_instability_K_parameter;
                            root_functions[i_root] = (rp_out - ra_in) - K*R_H;
                            //printf("rf %g %g %g %g %g\n",central_particle_mass,m1,m2,K,root_functions[i_root]);
                        }

                        else if (P_p->dynamical_instability_criterion == 3)
                        /* Petrovich 2015 */
                        {
                            root_functions[i_root] = rp_out/ra_in - ( 2.4*pow( max(mu1,mu2), c_1div3)*sqrt(a_out/a_in) + 1.15 );
                        }
                        
                    }

                    if (root_functions[i_root] <= 0.0)
                    {
                        P_p->dynamical_instability_has_occurred = 1;
//                        printf("stop\n");
                    }

                        //printf("di P15a %d %d %g %g %g\n",P_p->index,central_particle_index,m1,m2,central_particle_mass);
                        //printf("di P15b %g %g %g \n",a_in,a_out,root_functions[i_root]);

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

int read_root_finding_data(ParticlesMap *particlesMap, int *roots_found)
{
    ParticlesMapIterator it_p;
    
    int i_root = 0;
    
    for (it_p = particlesMap->begin(); it_p != particlesMap->end(); it_p++)
    {
        Particle *P_p = (*it_p).second;
        if (P_p->is_binary == 1)
        {
            if (P_p->check_for_secular_breakdown == 1)
            {
                if FOUND_ROOT
                {
                    P_p->secular_breakdown_has_occurred = 1;
                }
                i_root++;

            }
            if (P_p->check_for_dynamical_instability == 1)
            {
                if FOUND_ROOT
                {
                    P_p->dynamical_instability_has_occurred = 1;
                }
                i_root++;                
            }
            if (P_p->check_for_physical_collision_or_orbit_crossing == 1)
            {
                if FOUND_ROOT
                {
                    P_p->physical_collision_or_orbit_crossing_has_occurred = 1;
                }
                i_root++;                
            }
            if (P_p->check_for_minimum_periapse_distance == 1)
            {
                if FOUND_ROOT
                {
                    P_p->minimum_periapse_distance_has_occurred = 1;
                }
                i_root++;                
            }
        }
        else /* P_p not a binary */
        {
            if (P_p->check_for_RLOF_at_pericentre == 1)
            {
                if FOUND_ROOT
                {
                    P_p->RLOF_at_pericentre_has_occurred = 1;
                }
                i_root++;
            }
        }
    }
    return 0;
}

int check_for_initial_roots(ParticlesMap *particlesMap)
{
    ParticlesMapIterator it_p;
    
    int N_root_found = 0;
    
    for (it_p = particlesMap->begin(); it_p != particlesMap->end(); it_p++)
    {
        Particle *P_p = (*it_p).second;
        if (P_p->is_binary == 1)
        {
            if (P_p->check_for_secular_breakdown == 1)
            {
                if (P_p->secular_breakdown_has_occurred == 1)
                {
                    N_root_found++;
                }

            }
            if (P_p->check_for_dynamical_instability == 1)
            {
                if (P_p->dynamical_instability_has_occurred == 1)
                {
                    N_root_found++;
                }
            }
            if (P_p->check_for_physical_collision_or_orbit_crossing == 1)
            {
                if (P_p->physical_collision_or_orbit_crossing_has_occurred == 1)
                {
                    N_root_found++;
                }
            }
            if (P_p->check_for_minimum_periapse_distance == 1)
            {
                if (P_p->minimum_periapse_distance_has_occurred == 1)
                {
                    N_root_found++;
                }
            }            
        }
        else /* P_p not a binary */
        {
            if (P_p->check_for_RLOF_at_pericentre == 1)
            {
                if (P_p->RLOF_at_pericentre_has_occurred == 1)
                {
                    N_root_found++;
                }
            }
        }
    }
    return N_root_found;
}
