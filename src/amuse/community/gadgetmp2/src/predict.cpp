/* ################################################################################## */
/* ###                                                                            ### */
/* ###                                 Gadgetmp2                                  ### */
/* ###                                                                            ### */
/* ###   Original: Gadget2 in the version used in Amuse                           ### */
/* ###   Author: Gadget2 and Amuse contributors                                   ### */
/* ###                                                                            ### */
/* ###   Modified: July 2020                                                      ### */
/* ###   Author: Thomas Schano                                                    ### */
/* ###                                                                            ### */
/* ###   Changes are intended to enable precise calculations in                   ### */
/* ###   non periodic small domain simulations in which comoving parts            ### */
/* ###   are simulated in std precision                                           ### */
/* ###                                                                            ### */
/* ################################################################################## */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#ifndef NOMPI
#include <mpi.h>
#endif
#include <gsl/gsl_math.h>

//#include "allvars.hpp"
#include "proto.hpp"


/*! \file predict.c
 *  \brief drift particles by a small time interval
 *
 *  This function contains code to implement a drift operation on all the
 *  particles, which represents one part of the leapfrog integration scheme.
 */


/*! This function drifts all particles from the current time to the future:
 *  time0 - > time1
 *
 *  If there is no explicit tree construction in the following timestep, the
 *  tree nodes are also drifted and updated accordingly. Note: For periodic
 *  boundary conditions, the mapping of coordinates onto the interval
 *  [0,All.BoxSize] is only done before the domain decomposition, or for
 *  outputs to snapshot files.  This simplifies dynamic tree updates, and
 *  allows the domain decomposition to be carried out only every once in a
 *  while.
 */
void gadgetmp2::move_particles(int time0, int time1)
{
    int i, j;
    my_float dt_drift, dt_gravkick, dt_hydrokick, dt_entr;
    double t0, t1;


    t0 = second();

    if(All.ComovingIntegrationOn)
    {
        dt_drift = get_drift_factor(time0, time1);
        dt_gravkick = get_gravkick_factor(time0, time1);
        dt_hydrokick = get_hydrokick_factor(time0, time1);
    }
    else
    {
        dt_drift = dt_gravkick = dt_hydrokick = (time1 - time0) * All.Timebase_interval;
    }

    for(i = 0; i < NumPart; i++)
    {
        for(j = 0; j < 3; j++)
            P[i].Pos[j] += P[i].Vel[j] * dt_drift;

        if(P[i].Type == 0)
        {
            for(j = 0; j < 3; j++)
                SphP[i].VelPred[j] += P[i].GravAccel[j] * dt_gravkick + SphP[i].HydroAccel[j] * dt_hydrokick;
            SphP[i].Density *= exp(-SphP[i].DivVel * dt_drift);
            SphP[i].Hsml *= exp(const_0_333333333333 * SphP[i].DivVel * dt_drift);

            if(SphP[i].Hsml < All.MinGasHsml)
                SphP[i].Hsml = All.MinGasHsml;

            dt_entr = (time1 - (P[i].Ti_begstep + P[i].Ti_endstep) / const_2) * All.Timebase_interval;

            SphP[i].Pressure = (SphP[i].Entropy + SphP[i].DtEntropy * dt_entr) * pow(SphP[i].Density, const_GAMMA);

#ifdef MORRIS97VISC
            SphP[i].Alpha += SphP[i].DAlphaDt * dt_drift;
#endif
        }
    }

    /* if domain-decomp and tree are not going to be reconstructed, update dynamically.  */
    if(All.NumForcesSinceLastDomainDecomp < All.TotNumPart * All.TreeDomainUpdateFrequency)
    {
        for(i = 0; i < Numnodestree; i++)
            for(j = 0; j < 3; j++)
                Nodes[All.MaxPart + i].u_d_s[j] += Extnodes[All.MaxPart + i].vs[j] * dt_drift;

        force_update_len();

        force_update_pseudoparticles();
    }

    t1 = second();

    All.CPU_Predict += timediff(t0, t1);
}
