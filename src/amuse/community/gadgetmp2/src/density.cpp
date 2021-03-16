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


/*! \file density.c
 *  \brief SPH density computation and smoothing length determination
 *
 *  This file contains the "first SPH loop", where the SPH densities and
 *  some auxiliary quantities are computed.  If the number of neighbours
 *  obtained falls outside the target range, the correct smoothing
 *  length is determined iteratively, if needed.
 */



/*! This function computes the local density for each active SPH particle,
 *  the number of neighbours in the current smoothing radius, and the
 *  divergence and curl of the velocity field.  The pressure is updated as
 *  well.  If a particle with its smoothing region is fully inside the
 *  local domain, it is not exported to the other processors. The function
 *  also detects particles that have a number of neighbours outside the
 *  allowed tolerance range. For these particles, the smoothing length is
 *  adjusted accordingly, and the density computation is executed again.
 *  Note that the smoothing length is not allowed to fall below the lower
 *  bound set by MinGasHsml.
 */
void gadgetmp2::density(void)
{
    long long ntot, ntotleft;
    int *noffset, *nbuffer, *nsend, *nsend_local, *numlist, *ndonelist;
    int i, j, n, ndone, npleft, maxfill, source, iter = 0;
    int level, ngrp, sendTask, recvTask, place, nexport;
    double  tstart, tend, tstart_ngb = 0, tend_ngb = 0;
    my_float dt_entr;
    double sumt, sumcomm, timengb, sumtimengb;
    double timecomp = 0, timeimbalance = 0, timecommsumm = 0, sumimbalance;
#ifndef NOMPI
    MPI_Status status;
#endif

    noffset = new int[NTask];	/* offsets of bunches in common list */
    nbuffer = new int[NTask];
    nsend_local = new int[NTask];
    nsend = (int*)malloc(sizeof(int) * NTask * NTask);
    ndonelist = new int[NTask];

    for(n = 0, NumSphUpdate = 0; n < N_gas; n++)
    {
        SphP[n].Left.setZero(); SphP[n].Right.setZero();

        if(P[n].Ti_endstep == All.Ti_Current)
            NumSphUpdate++;
    }

    numlist = new int[NTask * NTask];
#ifndef NOMPI
    MPI_Allgather(&NumSphUpdate, 1, MPI_INT, numlist, 1, MPI_INT, GADGET_WORLD);
#else
    numlist[0] = NumSphUpdate;
#endif
    for(i = 0, ntot = 0; i < NTask; i++)
        ntot += numlist[i];
    free(numlist);



    /* we will repeat the whole thing for those particles where we didn't
     * find enough neighbours
     */
    do
    {
        i = 0;			/* beginn with this index */
        ntotleft = ntot;		/* particles left for all tasks together */

        while(ntotleft > 0)
        {
            for(j = 0; j < NTask; j++)
                nsend_local[j] = 0;

            /* do local particles and prepare export list */
            tstart = second();
            for(nexport = 0, ndone = 0; i < N_gas && nexport < All.BunchSizeDensity - NTask; i++)
                if(P[i].Ti_endstep == All.Ti_Current)
                {
                    ndone++;

                    for(j = 0; j < NTask; j++)
                        Exportflag[j] = 0;

                    density_evaluate(i, 0);

                    for(j = 0; j < NTask; j++)
                    {
                        if(Exportflag[j])
                        {
                            //DensDataIn[nexport].Pos[0] = P[i].Pos[0];
                            //DensDataIn[nexport].Pos[1] = P[i].Pos[1];
                            //DensDataIn[nexport].Pos[2] = P[i].Pos[2];
                            DensDataIn->set_init_Pos0(P[i].Pos[0],nexport);
                            DensDataIn->set_init_Pos1(P[i].Pos[1],nexport);
                            DensDataIn->set_init_Pos2(P[i].Pos[2],nexport);
                            //DensDataIn[nexport].Vel[0] = SphP[i].VelPred[0];
                            //DensDataIn[nexport].Vel[1] = SphP[i].VelPred[1];
                            //DensDataIn[nexport].Vel[2] = SphP[i].VelPred[2];
                            DensDataIn->set_init_Vel0(SphP[i].VelPred[0],nexport);
                            DensDataIn->set_init_Vel1(SphP[i].VelPred[1],nexport);
                            DensDataIn->set_init_Vel2(SphP[i].VelPred[2],nexport);
                            //DensDataIn[nexport].Hsml = SphP[i].Hsml;
                            DensDataIn->set_init_Hsml(SphP[i].Hsml,nexport);
                            //DensDataIn[nexport].Index = i;
                            DensDataIn->set_Index(i,nexport);
                            //DensDataIn[nexport].Task = j;
                            DensDataIn->set_Task(j,nexport);
                            nexport++;
                            nsend_local[j]++;
                        }
                    }
                }
            tend = second();
            timecomp += timediff(tstart, tend);

            qsort(DensDataIn, nexport, densdata_in::get_size(), dens_compare_key);

            for(j = 1, noffset[0] = 0; j < NTask; j++)
                noffset[j] = noffset[j - 1] + nsend_local[j - 1];

            tstart = second();

#ifndef NOMPI
            MPI_Allgather(nsend_local, NTask, MPI_INT, nsend, NTask, MPI_INT, GADGET_WORLD);
#else
            nsend[0] = nsend_local[0];
#endif
            tend = second();
            timeimbalance += timediff(tstart, tend);


            /* now do the particles that need to be exported */

            for(level = 1; level < (1 << PTask); level++)
            {
                tstart = second();
                for(j = 0; j < NTask; j++)
                    nbuffer[j] = 0;
                for(ngrp = level; ngrp < (1 << PTask); ngrp++)
                {
                    maxfill = 0;
                    for(j = 0; j < NTask; j++)
                    {
                        if((j ^ ngrp) < NTask)
                            if(maxfill < nbuffer[j] + nsend[(j ^ ngrp) * NTask + j])
                                maxfill = nbuffer[j] + nsend[(j ^ ngrp) * NTask + j];
                    }
                    if(maxfill >= All.BunchSizeDensity)
                        break;

                    sendTask = ThisTask;
                    recvTask = ThisTask ^ ngrp;

                    if(recvTask < NTask)
                    {
                        if(nsend[ThisTask * NTask + recvTask] > 0 || nsend[recvTask * NTask + ThisTask] > 0)
                        {
                            /* get the particles */
#ifndef NOMPI
                            MPI_Sendrecv(DensDataIn->get_buff_start(noffset[recvTask]),
                                         nsend_local[recvTask] * densdata_in::get_size(), MPI_BYTE,
                                         recvTask, TAG_DENS_A,
                                         DensDataGet->get_buff_start(nbuffer[ThisTask]),
                                         nsend[recvTask * NTask + ThisTask] * densdata_in::get_size(),
                                    MPI_BYTE, recvTask, TAG_DENS_A, GADGET_WORLD, &status);
#else
                            fprintf(stderr, "NO MPI, SO NO SENDING");
                            exit(1);
#endif
                        }
                    }

                    for(j = 0; j < NTask; j++)
                        if((j ^ ngrp) < NTask)
                            nbuffer[j] += nsend[(j ^ ngrp) * NTask + j];
                }
                tend = second();
                timecommsumm += timediff(tstart, tend);


                tstart = second();
                for(j = 0; j < nbuffer[ThisTask]; j++)
                    density_evaluate(j, 1);
                tend = second();
                timecomp += timediff(tstart, tend);

                /* do a block to explicitly measure imbalance */
                tstart = second();
#ifndef NOMPI
                MPI_Barrier(GADGET_WORLD);
#endif
                tend = second();
                timeimbalance += timediff(tstart, tend);

                /* get the result */
                tstart = second();
                for(j = 0; j < NTask; j++)
                    nbuffer[j] = 0;
                for(ngrp = level; ngrp < (1 << PTask); ngrp++)
                {
                    maxfill = 0;
                    for(j = 0; j < NTask; j++)
                    {
                        if((j ^ ngrp) < NTask)
                            if(maxfill < nbuffer[j] + nsend[(j ^ ngrp) * NTask + j])
                                maxfill = nbuffer[j] + nsend[(j ^ ngrp) * NTask + j];
                    }
                    if(maxfill >= All.BunchSizeDensity)
                        break;

                    sendTask = ThisTask;
                    recvTask = ThisTask ^ ngrp;

                    if(recvTask < NTask)
                    {
                        if(nsend[ThisTask * NTask + recvTask] > 0 || nsend[recvTask * NTask + ThisTask] > 0)
                        {
                            /* send the results */
#ifndef NOMPI
                            MPI_Sendrecv(DensDataResult->get_buff_start(nbuffer[ThisTask]),
                                         nsend[recvTask * NTask + ThisTask] * densdata_out::get_size(),
                                    MPI_BYTE, recvTask, TAG_DENS_B,
                                    DensDataPartialResult->get_buff_start(noffset[recvTask]),
                                    nsend_local[recvTask] * densdata_out::get_size(),
                                    MPI_BYTE, recvTask, TAG_DENS_B, GADGET_WORLD, &status);
#else
                            fprintf(stderr, "NO MPI, SO NO SENDING");
                            exit(1);
#endif

                            /* add the result to the particles */
                            for(j = 0; j < nsend_local[recvTask]; j++)
                            {
                                source = j + noffset[recvTask];
                                //place = DensDataIn[source].Index;
                                place = DensDataIn->read_Index(source);

                                //SphP[place].NumNgb += DensDataPartialResult[source].Ngb;
                                SphP[place].NumNgb += DensDataPartialResult->read_re_init_Ngb(source);
                                //SphP[place].Density += DensDataPartialResult[source].Rho;
                                SphP[place].Density += DensDataPartialResult->read_re_init_Rho(source);
                                //SphP[place].DivVel += DensDataPartialResult[source].Div;
                                SphP[place].DivVel += DensDataPartialResult->read_re_init_Div(source);

                                //SphP[place].DhsmlDensityFactor += DensDataPartialResult[source].DhsmlDensity;
                                SphP[place].DhsmlDensityFactor += DensDataPartialResult->read_re_init_DhsmlDensity(source);

                                //SphP[place].Rot[0] += DensDataPartialResult[source].Rot[0];
                                //SphP[place].Rot[1] += DensDataPartialResult[source].Rot[1];
                                //SphP[place].Rot[2] += DensDataPartialResult[source].Rot[2];
                                SphP[place].Rot[0] +=  DensDataPartialResult->read_re_init_Rot0(source);
                                SphP[place].Rot[1] +=  DensDataPartialResult->read_re_init_Rot1(source);
                                SphP[place].Rot[2] +=  DensDataPartialResult->read_re_init_Rot2(source);
                            }
                        }
                    }

                    for(j = 0; j < NTask; j++)
                        if((j ^ ngrp) < NTask)
                            nbuffer[j] += nsend[(j ^ ngrp) * NTask + j];
                }
                tend = second();
                timecommsumm += timediff(tstart, tend);

                level = ngrp - 1;
            }

#ifndef NOMPI
            MPI_Allgather(&ndone, 1, MPI_INT, ndonelist, 1, MPI_INT, GADGET_WORLD);
#else
            ndonelist[0] = ndone;
#endif
            for(j = 0; j < NTask; j++)
                ntotleft -= ndonelist[j];
        }



        /* do final operations on results */
        tstart = second();
        for(i = 0, npleft = 0; i < N_gas; i++)
        {
            if(P[i].Ti_endstep == All.Ti_Current)
            {
                {
                    SphP[i].DhsmlDensityFactor =
                            const_1 / (const_1 + SphP[i].Hsml * SphP[i].DhsmlDensityFactor / (NUMDIMS * SphP[i].Density));

                    SphP[i].CurlVel = sqrt(SphP[i].Rot[0] * SphP[i].Rot[0] +
                            SphP[i].Rot[1] * SphP[i].Rot[1] +
                            SphP[i].Rot[2] * SphP[i].Rot[2]) / SphP[i].Density;

                    SphP[i].DivVel /= SphP[i].Density;

                    dt_entr = (All.Ti_Current - (P[i].Ti_begstep + P[i].Ti_endstep) / const_2) * All.Timebase_interval;

                    SphP[i].Pressure =
                            (SphP[i].Entropy + SphP[i].DtEntropy * dt_entr) * pow(SphP[i].Density, const_GAMMA);
                }


                /* now check whether we had enough neighbours */

                if(SphP[i].NumNgb < (All.DesNumNgb - All.MaxNumNgbDeviation) ||
                        (SphP[i].NumNgb > (All.DesNumNgb + All.MaxNumNgbDeviation)
                         && SphP[i].Hsml > (1.01 * All.MinGasHsml)))
                {
                    /* need to redo this particle */
                    npleft++;

                    if(SphP[i].Left > 0 && SphP[i].Right > 0)
                        if((SphP[i].Right - SphP[i].Left) < 1.0e-3 * SphP[i].Left)
                        {
                            /* this one should be ok */
                            npleft--;
                            P[i].Ti_endstep = -P[i].Ti_endstep - 1;	/* Mark as inactive */
                            continue;
                        }

                    if(SphP[i].NumNgb < (All.DesNumNgb - All.MaxNumNgbDeviation))
                        SphP[i].Left = dmax(SphP[i].Hsml, SphP[i].Left);
                    else
                    {
                        if(SphP[i].Right != "0")
                        {
                            if(SphP[i].Hsml < SphP[i].Right)
                                SphP[i].Right = SphP[i].Hsml;
                        }
                        else
                            SphP[i].Right = SphP[i].Hsml;
                    }

                    if(iter >= MAXITER - 10)
                    {
                        printf
                                ("i=%d task=%d ID=%d Hsml=%g Left=%g Right=%g Ngbs=%g Right-Left=%g\n   pos=(%g|%g|%g)\n",
                                 i, ThisTask, (int) P[i].ID, SphP[i].Hsml.toDouble(), SphP[i].Left.toDouble(), SphP[i].Right.toDouble(),
                                 ((my_float) SphP[i].NumNgb).toDouble(), (SphP[i].Right - SphP[i].Left).toDouble(), P[i].Pos[0].toDouble(), P[i].Pos[1].toDouble(),
                                P[i].Pos[2].toDouble());
                        fflush(stdout);
                    }

                    if(SphP[i].Right > const_0 && SphP[i].Left > const_0)
                        SphP[i].Hsml = pow(const_0_5 * (pow(SphP[i].Left,const_3) + pow(SphP[i].Right,const_3)), const_1 / const_3);
                    else
                    {
                        if(SphP[i].Right == const_0 && SphP[i].Left == const_0)
                            endrun(8188);	/* can't occur */

                        if(SphP[i].Right == const_0 && SphP[i].Left > const_0)
                        {
                            if(P[i].Type == 0 && fabs(SphP[i].NumNgb - All.DesNumNgb) < const_0_5 * All.DesNumNgb)
                            {
                                SphP[i].Hsml *=
                                        const_1- (SphP[i].NumNgb -
                                                  All.DesNumNgb) / (NUMDIMS * SphP[i].NumNgb) * SphP[i].DhsmlDensityFactor;
                            }
                            else
                                SphP[i].Hsml *= const_1_26;
                        }

                        if(SphP[i].Right > const_0 && SphP[i].Left == const_0)
                        {
                            if(P[i].Type == 0 && fabs(SphP[i].NumNgb - All.DesNumNgb) < const_0_5 * All.DesNumNgb)
                            {
                                my_float d=(SphP[i].NumNgb -
                                            All.DesNumNgb) / (NUMDIMS * SphP[i].NumNgb) * SphP[i].DhsmlDensityFactor;
                                if(d<.99)
                                    SphP[i].Hsml *= const_1 - d;
                                else
                                    SphP[i].Hsml /= const_1_26;
                            }
                            else
                                SphP[i].Hsml /= const_1_26;
                        }
                    }

                    if(SphP[i].Hsml <= const_0)
                    {
                        printf
                                ("** i=%d task=%d ID=%d Hsml=%g Left=%g Right=%g Ngbs=%g dhsml=%g\n ",
                                 i, ThisTask, (int) P[i].ID, SphP[i].Hsml.toDouble(), SphP[i].Left.toDouble(), SphP[i].Right.toDouble(),
                                 ((my_float) SphP[i].NumNgb).toDouble(), SphP[i].DhsmlDensityFactor.toDouble());
                        fflush(stdout);
                    }


                    if(SphP[i].Hsml < All.MinGasHsml)
                        SphP[i].Hsml = All.MinGasHsml;
                }
                else
                    P[i].Ti_endstep = -P[i].Ti_endstep - 1;	/* Mark as inactive */
            }
        }
        tend = second();
        timecomp += timediff(tstart, tend);


        numlist = new int[NTask * NTask];

#ifndef NOMPI
        MPI_Allgather(&npleft, 1, MPI_INT, numlist, 1, MPI_INT, GADGET_WORLD);
#else
        numlist[0] = npleft;
#endif
        for(i = 0, ntot = 0; i < NTask; i++)
            ntot += numlist[i];
        free(numlist);

        if(ntot > 0)
        {
            if(iter == 0)
                tstart_ngb = second();

            iter++;

            if(iter > 0 && ThisTask == 0)
            {
                printf("ngb iteration %d: need to repeat for %d%09d particles.\n", iter,
                       (int) (ntot / 1000000000), (int) (ntot % 1000000000));
                fflush(stdout);
            }

            if(iter > MAXITER)
            {
                printf("failed to converge in neighbour iteration in density()\n");
                fflush(stdout);
                endrun(1155);
            }
        }
        else
            tend_ngb = second();
    }
    while(ntot > 0);


    /* mark as active again */
    for(i = 0; i < NumPart; i++)
        if(P[i].Ti_endstep < 0)
            P[i].Ti_endstep = -P[i].Ti_endstep - 1;

    free(ndonelist);
    free(nsend);
    free(nsend_local);
    free(nbuffer);
    free(noffset);


    /* collect some timing information */
    if(iter > 0)
        timengb = timediff(tstart_ngb, tend_ngb);
    else
        timengb = 0;

#ifndef NOMPI
    MPI_Reduce(&timengb, &sumtimengb, 1, MPI_DOUBLE, MPI_SUM, 0, GADGET_WORLD);
    MPI_Reduce(&timecomp, &sumt, 1, MPI_DOUBLE, MPI_SUM, 0, GADGET_WORLD);
    MPI_Reduce(&timecommsumm, &sumcomm, 1, MPI_DOUBLE, MPI_SUM, 0, GADGET_WORLD);
    MPI_Reduce(&timeimbalance, &sumimbalance, 1, MPI_DOUBLE, MPI_SUM, 0, GADGET_WORLD);
#else
    sumtimengb = timengb;
    sumt = timecomp;
    sumcomm = timecommsumm;
    sumimbalance = timeimbalance;
#endif

    if(ThisTask == 0)
    {
        All.CPU_HydCompWalk += sumt / NTask;
        All.CPU_HydCommSumm += sumcomm / NTask;
        All.CPU_HydImbalance += sumimbalance / NTask;
        All.CPU_EnsureNgb += sumtimengb / NTask;
    }
}



/*! This function represents the core of the SPH density computation. The
 *  target particle may either be local, or reside in the communication
 *  buffer. it computes the cumulative and rated density of neighbors
 */
void gadgetmp2::density_evaluate(int target, int mode)
{
    int j, n, startnode, numngb, numngb_inbox;
    my_float h, h2, fac, hinv, hinv3, hinv4;
    my_float rho, divv, wk, dwk;
    my_float dx, dy, dz, r, r2, u, mass_j;
    my_float dvx, dvy, dvz, rotv[3];
    my_float weighted_numngb, dhsmlrho;
    my_float pos[3], vel[3];

    if(mode == 0)
    {
        pos[0] = P[target].Pos[0];
        pos[1] = P[target].Pos[1];
        pos[2] = P[target].Pos[2];
        vel[0] = SphP[target].VelPred[0];
        vel[1] = SphP[target].VelPred[1];
        vel[2] = SphP[target].VelPred[2];
        h = SphP[target].Hsml;
    }
    else
    {
        //pos = DensDataGet[target].Pos;
        pos[0] =  DensDataGet->read_re_init_Pos0(target);
        pos[1] =  DensDataGet->read_re_init_Pos1(target);
        pos[2] =  DensDataGet->read_re_init_Pos2(target);
        //vel = DensDataGet[target].Vel;
        vel[0] =  DensDataGet->read_re_init_Vel0(target);
        vel[1] =  DensDataGet->read_re_init_Vel1(target);
        vel[2] =  DensDataGet->read_re_init_Vel2(target);
        //h = DensDataGet[target].Hsml;
        h =  DensDataGet->read_re_init_Hsml(target);
    }

    h2 = h * h;
    hinv = const_1 / h;
#ifndef  TWODIMS
    hinv3 = hinv * hinv * hinv;
#else
    hinv3 = hinv * hinv / boxSize_Z;
#endif
    hinv4 = hinv3 * hinv;

    rho.setZero(); divv.setZero(); rotv[0].setZero(); rotv[1].setZero(); rotv[2].setZero();
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
                    dwk = hinv4 * KERNEL_COEFF_6 * (const_1 - u) * (const_1 - u);
                }

                mass_j = P[j].Mass;

                rho += mass_j * wk;

                weighted_numngb += NORM_COEFF * wk / hinv3;

                dhsmlrho += -mass_j * (NUMDIMS * hinv * wk + u * dwk);

                if(r > const_0)
                {
                    fac = mass_j * dwk / r;

                    dvx = vel[0] - SphP[j].VelPred[0];
                    dvy = vel[1] - SphP[j].VelPred[1];
                    dvz = vel[2] - SphP[j].VelPred[2];

                    divv -= fac * (dx * dvx + dy * dvy + dz * dvz);

                    rotv[0] += fac * (dz * dvy - dy * dvz);
                    rotv[1] += fac * (dx * dvz - dz * dvx);
                    rotv[2] += fac * (dy * dvx - dx * dvy);
                }
            }
        }
    }
    while(startnode >= 0);

    if(mode == 0)
    {
        SphP[target].NumNgb = weighted_numngb;
        SphP[target].Density = rho;
        SphP[target].DivVel = divv;
        SphP[target].DhsmlDensityFactor = dhsmlrho;
        SphP[target].Rot[0] = rotv[0];
        SphP[target].Rot[1] = rotv[1];
        SphP[target].Rot[2] = rotv[2];
    }
    else
    {
        //DensDataResult[target].Rho = rho;
        DensDataResult->set_init_Rho(rho,target);
        //DensDataResult[target].Div = divv;
        DensDataResult->set_init_Div(divv,target);
        //DensDataResult[target].Ngb = weighted_numngb;
        DensDataResult->set_init_Ngb(weighted_numngb,target);
        //DensDataResult[target].DhsmlDensity = dhsmlrho;
        DensDataResult->set_init_DhsmlDensity(dhsmlrho,target);
        //DensDataResult[target].Rot[0] = rotv[0];
        //DensDataResult[target].Rot[1] = rotv[1];
        //DensDataResult[target].Rot[2] = rotv[2];
        DensDataResult->set_init_Rot0(rotv[0],target);
        DensDataResult->set_init_Rot1(rotv[1],target);
        DensDataResult->set_init_Rot2(rotv[2],target);
    }
}




/*! This routine is a comparison kernel used in a sort routine to group
 *  particles that are exported to the same processor.
 */
int gadgetmp2::dens_compare_key(const void *a, const void *b)
{
    if(((struct densdata_in *) a)->read_Task() < (((struct densdata_in *) b)->read_Task()))
        return -1;

    if(((struct densdata_in *) a)->read_Task() > (((struct densdata_in *) b)->read_Task()))
        return +1;

    return 0;
}
