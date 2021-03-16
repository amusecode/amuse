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
#include <stdlib.h>
#include <string.h>
#include <math.h>
#ifndef NOMPI
#include <mpi.h>
#endif

//#include "allvars.hpp"
#include "proto.hpp"


/*! \file init.c
 *  \brief Code for initialisation of a simulation from initial conditions
 */



/*! This function is used to find an initial smoothing length for each SPH
 *  particle. It guarantees that the number of neighbours will be between
 *  desired_ngb-MAXDEV and desired_ngb+MAXDEV. For simplicity, a first guess
 *  of the smoothing length is provided to the function density(), which will
 *  then iterate if needed to find the right smoothing length.
 */
void gadgetmp2::setup_smoothinglengths(void)
{
    int i, no, p;

    if(RestartFlag == 0)
    {

        for(i = 0; i < N_gas; i++)
        {
            no = Father[i];

            while(10 * All.DesNumNgb * P[i].Mass > Nodes[no].u_d_mass)
            {
                p = Nodes[no].u.d.father;

                if(p < 0)
                    break;

                no = p;
            }
#ifndef TWODIMS
            SphP[i].Hsml =
                    pow(const_3 / (const_4 * const_PI) * All.DesNumNgb * P[i].Mass / Nodes[no].u_d_mass,  const_1 / const_3) * Nodes[no].len;
#else
            SphP[i].Hsml =
                    pow(const_1 / (const_PI) * All.DesNumNgb * P[i].Mass / Nodes[no].u.d.mass, const_1 / const_2) * Nodes[no].len;
#endif
        }
    }

    density();
}

