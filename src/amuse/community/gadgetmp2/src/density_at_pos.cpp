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

//#include "allvars.hpp"
#include "proto.hpp"

const int debug = 0;
/*void hydro_state_evaluate(my_float h, my_float pos[3], my_float vel[3], my_float *numngb,
    my_float *dhsml_out, my_float *rho_out, my_float *rhov_out, my_float *rhov2_out, my_float
    *rhoe_out);
*/
void gadgetmp2::hydro_state_at_point(my_float pos[3], my_float vel[3], my_float *h_out, my_float *ngb_out,
my_float *dhsml_out, my_float *rho_out, my_float *rhov_out, my_float *rhov2_out, my_float *rhoe_out)
{
    my_float low, up, h, dhsml, low_ngb, up_ngb, ngb;
    my_float rho, rhov[3], rhov2, rhoe;
    int i, iter;

    up = low = SphP[0].Hsml;
    for(i = 1; i < N_gas; i++){
        if (low > SphP[i].Hsml)
            low = SphP[i].Hsml;
        else if (up < SphP[i].Hsml)
            up = SphP[i].Hsml;
    }
#ifndef NOMPI
    //MPI_Allreduce(MPI_IN_PLACE, &low, 1, MPI_DOUBLE, MPI_MIN, GADGET_WORLD);
    low=mpi_all_min(low);
    //MPI_Allreduce(MPI_IN_PLACE, &up,  1, MPI_DOUBLE, MPI_MAX, GADGET_WORLD);
    up=mpi_all_max(up);
#endif
    low *= const_1_26;
    up /= const_1_26;

    iter = 0;
    do {
        low /= const_1_26;
        hydro_state_evaluate(low, pos, vel, &low_ngb, &dhsml, &rho, rhov, &rhov2, &rhoe);
#ifndef NOMPI
        //MPI_Allreduce(MPI_IN_PLACE, &low_ngb, 1, MPI_DOUBLE, MPI_SUM, GADGET_WORLD);
        low_ngb=mpi_all_sum(low_ngb);
#endif
        if (iter > MAXITER)
            endrun(3210);
        if (debug) printf("%d - Searching for lower h boundary: %f (ngb: %f)\n",iter, low.toDouble(), low_ngb.toDouble());
        iter++;
    } while (low_ngb > All.DesNumNgb);
    iter = 0;
    do {
        up *= const_1_26;
        hydro_state_evaluate(up, pos, vel, &up_ngb, &dhsml, &rho, rhov, &rhov2, &rhoe);
#ifndef NOMPI
        //MPI_Allreduce(MPI_IN_PLACE, &up_ngb, 1, MPI_DOUBLE, MPI_SUM, GADGET_WORLD);
        up_ngb=mpi_all_sum(up_ngb);
#endif
        if (iter > MAXITER)
            endrun(3211);
        if (debug) printf("%d - Searching for upper h boundary: %f (ngb: %f)\n",iter, up.toDouble(), up_ngb.toDouble());
        iter++;
    } while (up_ngb < All.DesNumNgb);

    iter = 0;
    ngb = All.DesNumNgb + 2*All.MaxNumNgbDeviation; // Makes sure first evaluation of condition is true:
    while (fabs(All.DesNumNgb - ngb) > All.MaxNumNgbDeviation) {
        h = pow(const_0_5 * (pow(low, const_3) + pow(up, const_3)), const_1 / const_3);
        hydro_state_evaluate(h, pos, vel, &ngb, &dhsml, &rho, rhov, &rhov2, &rhoe);
#ifndef NOMPI
        //MPI_Allreduce(MPI_IN_PLACE, &ngb, 1, MPI_DOUBLE, MPI_SUM, GADGET_WORLD);
        ngb=mpi_all_sum(ngb);
#endif

        if (ngb > All.DesNumNgb){
            up = h;
            if (up <= All.MinGasHsml)
                break;
        } else
            low = h;

        if(iter > MAXITER)
            endrun(3212);
        if (debug) printf("%d - Searching for h: %f (ngb: %f)\n",iter, h.toDouble(), ngb.toDouble());
        iter++;
    }

    if (h <= All.MinGasHsml)
        h = All.MinGasHsml;

    //    ngb = dhsml = rho = rhov2 = rhoe = rhov[0] = rhov[1] = rhov[2] = 0;
    hydro_state_evaluate(h, pos, vel, &ngb, &dhsml, &rho, rhov, &rhov2, &rhoe);
#ifndef NOMPI
    //MPI_Allreduce(MPI_IN_PLACE, &ngb,   1, MPI_DOUBLE, MPI_SUM, GADGET_WORLD);
    ngb=mpi_all_sum(ngb);
    //MPI_Allreduce(MPI_IN_PLACE, &dhsml, 1, MPI_DOUBLE, MPI_SUM, GADGET_WORLD);
    dhsml=mpi_all_sum(dhsml);
    //MPI_Allreduce(MPI_IN_PLACE, &rho,   1, MPI_DOUBLE, MPI_SUM, GADGET_WORLD);
    rho=mpi_all_sum(rho);
    //MPI_Allreduce(MPI_IN_PLACE, &rhov,  3, MPI_DOUBLE, MPI_SUM, GADGET_WORLD);
    rhov[0]=mpi_all_sum(rhov[0]);
    rhov[1]=mpi_all_sum(rhov[1]);
    rhov[2]=mpi_all_sum(rhov[2]);
    //MPI_Allreduce(MPI_IN_PLACE, &rhov2, 1, MPI_DOUBLE, MPI_SUM, GADGET_WORLD);
    rhov2=mpi_all_sum(rhov2);
    //MPI_Allreduce(MPI_IN_PLACE, &rhoe,  1, MPI_DOUBLE, MPI_SUM, GADGET_WORLD);
    rhoe=mpi_all_sum(rhoe);
#endif

    *h_out      = h;
    *ngb_out    = ngb;
    *dhsml_out  = dhsml;
    *rho_out    = rho;
    *rhov2_out  = rhov2;
    *rhoe_out   = rhoe;
    rhov_out[0] = rhov[0];
    rhov_out[1] = rhov[1];
    rhov_out[2] = rhov[2];
}

void gadgetmp2::hydro_state_evaluate(my_float h, my_float pos[3], my_float vel[3],
my_float *numngb_out, my_float *dhsml_out, my_float *rho_out, my_float *rhov_out,
my_float *rhov2_out, my_float *rhoe_out)
{
    int j, n, startnode, numngb, numngb_inbox;
    my_float h2, fac, hinv, hinv3, hinv4;
    my_float rho, rhov[3], rhov2, rhoe, wk, dwk;
    my_float dx, dy, dz, r, r2, u, mass_j;
    my_float dvx, dvy, dvz;
    my_float weighted_numngb, dhsmlrho;

    h2 = h * h;
    hinv = const_1 / h;
#ifndef  TWODIMS
    hinv3 = hinv * hinv * hinv;
#else
    hinv3 = hinv * hinv / boxSize_Z;
#endif
    hinv4 = hinv3 * hinv;

    rho.setZero(); rhov2.setZero(); rhoe.setZero(); rhov[0].setZero(); rhov[1].setZero(); rhov[2].setZero();
    weighted_numngb.setZero();
    dhsmlrho.setZero();

    startnode = All.MaxPart;
    numngb = 0;
    do
    {
        numngb_inbox = ngb_treefind_variable(&pos[0], h, &startnode);

        for(n = 0; n < numngb_inbox; n++)
        {
            j = Ngblist[n];

            dx = pos[0] - P[j].Pos[0];
            dy = pos[1] - P[j].Pos[1];
            dz = pos[2] - P[j].Pos[2];

            r2 = dx * dx + dy * dy + dz * dz;

            if(r2 < h2)
            {
                numngb++;

                r = sqrt(r2);

                u = r * hinv;

                if(u < const_0_5)
                {
                    wk = hinv3 * (KERNEL_COEFF_1 + KERNEL_COEFF_2 * (u - const_1) * u * u);
                    dwk = hinv4 * u * (KERNEL_COEFF_3 * u - KERNEL_COEFF_4);
                }
                else
                {
                    wk = hinv3 * KERNEL_COEFF_5 * (const_1 - u) * (const_1 - u) * (const_1 - u);
                    dwk = hinv4 * KERNEL_COEFF_6 * (const_1- u) * (const_1 - u);
                }

                mass_j = P[j].Mass;

                fac = mass_j * wk;
                rho += fac;

                weighted_numngb += NORM_COEFF * wk / hinv3;

                dhsmlrho += -mass_j * (NUMDIMS * hinv * wk + u * dwk);

                dvx = vel[0] - SphP[j].VelPred[0];
                dvy = vel[1] - SphP[j].VelPred[1];
                dvz = vel[2] - SphP[j].VelPred[2];

                rhov[0] -= fac * dvx;
                rhov[1] -= fac * dvy;
                rhov[2] -= fac * dvz;

                rhoe  += fac * SphP[j].Entropy;

                rhov2 += fac * (dvx*dvx + dvy*dvy + dvz*dvz);
            }
        }
    }
    while(startnode >= 0);

    *rho_out    = rho;
    *rhoe_out   = rhoe;
    *rhov2_out  = rhov2;
    *dhsml_out  = dhsmlrho;
    *numngb_out = weighted_numngb;
    rhov_out[0] = rhov[0];
    rhov_out[1] = rhov[1];
    rhov_out[2] = rhov[2];
}
