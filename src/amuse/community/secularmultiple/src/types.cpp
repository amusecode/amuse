#include "types.h"

void Particle::set_ODE_quantities()
{

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
    
    child1_mass_plus_child2_mass = child1_mass + child2_mass;
    child1_mass_minus_child2_mass = child1_mass - child2_mass;
    child1_mass_times_child2_mass = child1_mass*child2_mass;
    
    a = h*h*child1_mass_plus_child2_mass/( CONST_G*child1_mass_times_child2_mass*child1_mass_times_child2_mass*j_p2 );
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

int determine_binary_parents_levels_and_masses(ParticlesMap *particlesMap, int *N_bodies, int *N_binaries, int *N_root_finding)
{
    *N_bodies = 0;
    *N_binaries = 0;
    *N_root_finding = 0;
    
    /* determine parent for each particle */
    ParticlesMapIterator it_p,it_q;
    for (it_p = particlesMap->begin(); it_p != particlesMap->end(); it_p++)
    {
        Particle *P_p = (*it_p).second;
        
        P_p->parent = -1;
        
        if (P_p->is_binary == 1)
        {
            (*N_binaries)++;
            
            /* root finding */
            if (P_p->check_for_secular_breakdown == 1)
            {
                (*N_root_finding)++;
            }
            if (P_p->check_for_dynamical_instability == 1)
            {
                (*N_root_finding)++;
            }
            if (P_p->check_for_physical_collision_or_orbit_crossing == 1)
            {
                (*N_root_finding)++;
            }
            if (P_p->check_for_minimum_periapse_distance == 1)
            {
                (*N_root_finding)++;
            }

            /* parents and siblings */
            for (it_q = particlesMap->begin(); it_q != particlesMap->end(); it_q++)
            {
                Particle *P_q = (*it_q).second;
                
                if (P_q->index == P_p->child1)
                {
                    P_q->parent = P_p->index;
                    P_q->sibling = P_p->child2;
                }
                if (P_q->index == P_p->child2)
                {
                    P_q->parent = P_p->index;
                    P_q->sibling = P_p->child1;
                }
            }
        }
        else
        {
            (*N_bodies)++;
            
            /* root finding */
            if (P_p->check_for_RLOF_at_pericentre == 1)
            {
                (*N_root_finding)++;
            }
            
        }
    }

    /* determine levels and set of parents for each particle */
    int highest_level = 0;
    for (it_p = particlesMap->begin(); it_p != particlesMap->end(); it_p++)
    {
        Particle *P_p = (*it_p).second;

        P_p->connecting_child_in_parents.clear();
        P_p->parents.clear();
        P_p->level=0;
        
        int child = P_p->index;
        int parent = P_p->parent;

        if (parent != -1) /* if parent == -1, P_p is the `top' binary, for which level=0 */
        {
            while (parent != -1) /* search parents until reaching the top binary */
            {
                for (it_q = particlesMap->begin(); it_q != particlesMap->end(); it_q++)
                {
                    Particle *P_q = (*it_q).second;
                    if (P_q->index == parent)
                    {
                        if (child==P_q->child1)
                        {
                            P_p->connecting_child_in_parents.push_back(1);
                        }
                        else if (child==P_q->child2)
                        {
                            P_p->connecting_child_in_parents.push_back(2);
                        }
                        P_p->parents.push_back(parent);
                        P_p->level++;
                        
                        child = P_q->index;
                        parent = P_q->parent;
                        break;
                    }
                }
            }
        }
        highest_level = max(P_p->level,highest_level);
    }
    
    /* set binary masses -- to ensure this happens correctly, do this from highest level to lowest level */
    int level=highest_level;
    while (level > -1)
    {
        for (it_p = particlesMap->begin(); it_p != particlesMap->end(); it_p++)
        {
            Particle *P_p = (*it_p).second;
            if ((P_p->is_binary == 1) && (P_p->level == level))
            {
                Particle *P_p_child1 = (*particlesMap)[P_p->child1];
                Particle *P_p_child2 = (*particlesMap)[P_p->child2];
                
                /* these quantities are used in ODE_system.cpp */
                P_p->child1_mass = P_p_child1->mass;
                P_p->child2_mass = P_p_child2->mass;
                P_p->mass = P_p->child1_mass + P_p->child2_mass;

            }
        }
        level--;
    }
}
