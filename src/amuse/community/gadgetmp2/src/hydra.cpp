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

/*! \file hydra.c
 *  \brief Computation of SPH forces and rate of entropy generation
 *
 *  This file contains the "second SPH loop", where the SPH forces are
 *  computed, and where the rate of change of entropy due to the shock heating
 *  (via artificial viscosity) is computed.
 */


static my_float hubble_a, atime, hubble_a2, fac_mu, fac_vsic_fix, a3inv, fac_egy;


/*! This function is the driver routine for the calculation of hydrodynamical
 *  force and rate of change of entropy due to shock heating for all active
 *  particles .
 */
void gadgetmp2::hydro_force(void)
{
    long long ntot, ntotleft;
    int i, j, k, n, ngrp, maxfill, source, ndone;
    int *nbuffer, *noffset, *nsend_local, *nsend, *numlist, *ndonelist;
    int level, sendTask, recvTask, nexport, place;
    my_float soundspeed_i;
    double tstart, tend, sumt, sumcomm;
    double timecomp = 0, timecommsumm = 0, timeimbalance = 0, sumimbalance;
#ifndef NOMPI
    MPI_Status status;
#endif

    if(All.ComovingIntegrationOn)
    {
        /* Factors for comoving integration of hydro */
        hubble_a = All.Omega0 / (All.Time * All.Time * All.Time)
                + (1 - All.Omega0 - All.OmegaLambda) / (All.Time * All.Time) + All.OmegaLambda;

        hubble_a = All.Hubble * sqrt(hubble_a);
        hubble_a2 = All.Time * All.Time * hubble_a;

        fac_mu = pow(All.Time, 3 * (const_GAMMA - 1) / 2) / All.Time;

        fac_egy = pow(All.Time, 3 * (const_GAMMA - 1));

        fac_vsic_fix = hubble_a * pow(All.Time, 3 * const_GAMMA_MINUS1);

        a3inv = 1 / (All.Time * All.Time * All.Time);
        atime = All.Time;
    }
    else
        hubble_a = hubble_a2 = atime = fac_mu = fac_vsic_fix = a3inv = fac_egy = const_1;


    /* `NumSphUpdate' gives the number of particles on this processor that want a force update */
    for(n = 0, NumSphUpdate = 0; n < N_gas; n++)
    {
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


    noffset = new int[NTask];	/* offsets of bunches in common list */
    nbuffer = new int[NTask];
    nsend_local = new int[NTask];
    nsend = new int[NTask * NTask];
    ndonelist = new int[NTask];


    i = 0;			/* first particle for this task */
    ntotleft = ntot;		/* particles left for all tasks together */

    while(ntotleft > 0)
    {
        for(j = 0; j < NTask; j++)
            nsend_local[j] = 0;

        /* do local particles and prepare export list */
        tstart = second();
        for(nexport = 0, ndone = 0; i < N_gas && nexport < All.BunchSizeHydro - NTask; i++)
            if(P[i].Ti_endstep == All.Ti_Current)
            {
                ndone++;

                for(j = 0; j < NTask; j++)
                    Exportflag[j] = 0;

                hydro_evaluate(i, 0);

                for(j = 0; j < NTask; j++)
                {
                    if(Exportflag[j])
                    {
                        HydroDataIn->set_init_Pos0(P[i].Pos[0],nexport);
                        HydroDataIn->set_init_Pos1(P[i].Pos[1],nexport);
                        HydroDataIn->set_init_Pos2(P[i].Pos[2],nexport);
                        HydroDataIn->set_init_Vel0(SphP[i].VelPred[0],nexport);
                        HydroDataIn->set_init_Vel1(SphP[i].VelPred[1],nexport);
                        HydroDataIn->set_init_Vel2(SphP[i].VelPred[2],nexport);
                        HydroDataIn->set_init_Hsml(SphP[i].Hsml,nexport);
                        HydroDataIn->set_init_Mass(P[i].Mass,nexport);
                        HydroDataIn->set_init_DhsmlDensityFactor(SphP[i].DhsmlDensityFactor,nexport);
                        HydroDataIn->set_init_Density(SphP[i].Density,nexport);
                        HydroDataIn->set_init_Pressure(SphP[i].Pressure,nexport);
                        HydroDataIn->set_Timestep(P[i].Ti_endstep - P[i].Ti_begstep,nexport);

                        /* calculation of F1 */
                        soundspeed_i = sqrt(const_GAMMA * SphP[i].Pressure / SphP[i].Density);
                        HydroDataIn->set_init_F1(
                                                    fabs(SphP[i].DivVel) / (fabs(SphP[i].DivVel) +
                                                    SphP[i].CurlVel + const_0_0001 * soundspeed_i /
                                                    SphP[i].Hsml / fac_mu)
                                                    ,nexport);

                        HydroDataIn->set_Index(i,nexport);
                        HydroDataIn->set_Task(j,nexport);
#ifdef MORRIS97VISC
                        HydroDataIn->set_init_Alpha(SphP[i].Alpha,nexport);
#endif
                        nexport++;
                        nsend_local[j]++;
                    }
                }
            }
        tend = second();
        timecomp += timediff(tstart, tend);

        qsort(HydroDataIn, nexport, hydrodata_in::get_size(), hydro_compare_key);

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
                if(maxfill >= All.BunchSizeHydro)
                    break;

                sendTask = ThisTask;
                recvTask = ThisTask ^ ngrp;

                if(recvTask < NTask)
                {
                    if(nsend[ThisTask * NTask + recvTask] > 0 || nsend[recvTask * NTask + ThisTask] > 0)
                    {
#ifndef NOMPI
                        /* get the particles */
                        MPI_Sendrecv(HydroDataIn->get_buff_start(noffset[recvTask]),
                                     nsend_local[recvTask] * hydrodata_in::get_size(), MPI_BYTE,
                                     recvTask, TAG_HYDRO_A,
                                     HydroDataGet->get_buff_start(nbuffer[ThisTask]),
                                     nsend[recvTask * NTask + ThisTask] * hydrodata_in::get_size(), MPI_BYTE,
                                recvTask, TAG_HYDRO_A, GADGET_WORLD, &status);
#else
                        fprintf(stderr, "NOT SUPPORTED");
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

            /* now do the imported particles */
            tstart = second();
            for(j = 0; j < nbuffer[ThisTask]; j++)
                hydro_evaluate(j, 1);
            tend = second();
            timecomp += timediff(tstart, tend);

            /* do a block to measure imbalance */
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
                if(maxfill >= All.BunchSizeHydro)
                    break;

                sendTask = ThisTask;
                recvTask = ThisTask ^ ngrp;

                if(recvTask < NTask)
                {
                    if(nsend[ThisTask * NTask + recvTask] > 0 || nsend[recvTask * NTask + ThisTask] > 0)
                    {
#ifndef NOMPI
                        /* send the results */
                        MPI_Sendrecv(HydroDataResult->get_buff_start(nbuffer[ThisTask]),
                                     nsend[recvTask * NTask + ThisTask] * hydrodata_out::get_size(),
                                MPI_BYTE, recvTask, TAG_HYDRO_B,
                                HydroDataPartialResult->get_buff_start(noffset[recvTask]),
                                nsend_local[recvTask] * hydrodata_out::get_size(),
                                MPI_BYTE, recvTask, TAG_HYDRO_B, GADGET_WORLD, &status);
#else
                        fprintf(stderr, "NOT SUPPORTED");
                        exit(1);
#endif
                        /* add the result to the particles */
                        for(j = 0; j < nsend_local[recvTask]; j++)
                        {
                            source = j + noffset[recvTask];
                            place = HydroDataIn->read_Index(source);

                            SphP[place].HydroAccel[0] +=  HydroDataPartialResult->read_re_init_Acc0(source);
                            SphP[place].HydroAccel[1] +=  HydroDataPartialResult->read_re_init_Acc1(source);
                            SphP[place].HydroAccel[2] +=  HydroDataPartialResult->read_re_init_Acc2(source);

                            SphP[place].DtEntropy += HydroDataPartialResult->read_re_init_DtEntropy(source);

                            if(SphP[place].MaxSignalVel < HydroDataPartialResult->read_re_init_MaxSignalVel(source))
                                SphP[place].MaxSignalVel = HydroDataPartialResult->read_re_init_MaxSignalVel(source);
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

    free(ndonelist);
    free(nsend);
    free(nsend_local);
    free(nbuffer);
    free(noffset);



    /* do final operations on results */
    tstart = second();

    for(i = 0; i < N_gas; i++)
        if(P[i].Ti_endstep == All.Ti_Current)
        {
            SphP[i].DtEntropy *= const_GAMMA_MINUS1 / (hubble_a2 * pow(SphP[i].Density, const_GAMMA_MINUS1));
#ifdef SPH_BND_PARTICLES
            if(P[i].ID == 0)
            {
                SphP[i].DtEntropy = const_0;
                for(k = 0; k < 3; k++)
                    SphP[i].HydroAccel[k] = const_0;
            }
#endif
        }

    tend = second();
    timecomp += timediff(tstart, tend);

    /* collect some timing information */

#ifndef NOMPI
    MPI_Reduce(&timecomp, &sumt, 1, MPI_DOUBLE, MPI_SUM, 0, GADGET_WORLD);
    MPI_Reduce(&timecommsumm, &sumcomm, 1, MPI_DOUBLE, MPI_SUM, 0, GADGET_WORLD);
    MPI_Reduce(&timeimbalance, &sumimbalance, 1, MPI_DOUBLE, MPI_SUM, 0, GADGET_WORLD);
#else
    sumt = timecomp;
    sumcomm = timecommsumm;
    sumimbalance = timeimbalance;
#endif
    if(ThisTask == 0)
    {
        All.CPU_HydCompWalk += sumt / NTask;
        All.CPU_HydCommSumm += sumcomm / NTask;
        All.CPU_HydImbalance += sumimbalance / NTask;
    }
}


/*! This function is the 'core' of the SPH force computation. A target
 *  particle is specified which may either be local, or reside in the
 *  communication buffer.
 */
void gadgetmp2::hydro_evaluate(int target, int mode)
{
    int j, k, n, timestep, startnode, numngb;
    my_float pos[3], vel[3];
    my_float mass, h_i, dhsmlDensityFactor, rho, pressure, f1, f2;
#ifdef MORRIS97VISC
    my_float alpha_visc, alpha_visc_j;
#endif
    my_float acc[3], dtEntropy, maxSignalVel;
    my_float dx, dy, dz, dvx, dvy, dvz;
    my_float h_i2, hinv, hinv4;
    my_float p_over_rho2_i, p_over_rho2_j, soundspeed_i, soundspeed_j;
#ifdef MONAGHAN83VISC
    my_float soundspeed_ij, h_ij;
#endif
    my_float hfc, dwk_i, vdotr, vdotr2, visc, mu_ij, rho_ij, vsig;
    my_float h_j, dwk_j, r, r2, u, hfc_visc;

#ifndef NOVISCOSITYLIMITER
    my_float dt;
#endif

    if(mode == 0)
    {
        pos[0] = P[target].Pos[0];
        pos[1] = P[target].Pos[1];
        pos[2] = P[target].Pos[2];
        vel[0] = SphP[target].VelPred[0];
        vel[1] = SphP[target].VelPred[1];
        vel[2] = SphP[target].VelPred[2];
        h_i = SphP[target].Hsml;
        mass = P[target].Mass;
        dhsmlDensityFactor = SphP[target].DhsmlDensityFactor;
        rho = SphP[target].Density;
        pressure = SphP[target].Pressure;
        timestep = P[target].Ti_endstep - P[target].Ti_begstep;
        soundspeed_i = sqrt(const_GAMMA * pressure / rho);
#ifdef MORRIS97VISC
        alpha_visc = SphP[target].Alpha;
#else
        f1 = fabs(SphP[target].DivVel) /
                (fabs(SphP[target].DivVel) + SphP[target].CurlVel +
                 const_0_0001 * soundspeed_i / SphP[target].Hsml / fac_mu);
#endif
    }
    else
    {
        //pos = HydroDataGet[target].Pos;
        pos[0]       =  HydroDataGet->read_re_init_Pos0(target);
        pos[1]       =  HydroDataGet->read_re_init_Pos1(target);
        pos[2]       =  HydroDataGet->read_re_init_Pos2(target);
        //vel = HydroDataGet[target].Vel;
        vel[0]       =  HydroDataGet->read_re_init_Vel0(target);
        vel[1]       =  HydroDataGet->read_re_init_Vel1(target);
        vel[2]       =  HydroDataGet->read_re_init_Vel2(target);
        //h_i = HydroDataGet[target].Hsml;
        h_i =  HydroDataGet->read_re_init_Hsml(target);
        //mass = HydroDataGet[target].Mass;
        mass =  HydroDataGet->read_re_init_Mass(target);
        //dhsmlDensityFactor = HydroDataGet[target].DhsmlDensityFactor;
        dhsmlDensityFactor =  HydroDataGet->read_re_init_DhsmlDensityFactor(target);
        //rho = HydroDataGet[target].Density;
        rho =  HydroDataGet->read_re_init_Density(target);
        //pressure = HydroDataGet[target].Pressure;
        pressure =  HydroDataGet->read_re_init_Pressure(target);
        //timestep = HydroDataGet[target].Timestep;
        timestep =  HydroDataGet->read_Timestep(target);
        soundspeed_i = sqrt(const_GAMMA * pressure / rho);
        //f1 = HydroDataGet[target].F1;
        f1 =  HydroDataGet->read_re_init_F1(target);
#ifdef MORRIS97VISC
        //alpha_visc = HydroDataGet[target].Alpha;
        alpha_visc =  HydroDataGet->read_re_init_Alpha(target);
#endif
    }


    /* initialize variables before SPH loop is started */
    acc[0] = acc[1] = acc[2] = dtEntropy = const_0;
    maxSignalVel = const_0;

    p_over_rho2_i = pressure / (rho * rho) * dhsmlDensityFactor;
    h_i2 = h_i * h_i;

    /* Now start the actual SPH computation for this particle */
    startnode = All.MaxPart;
    do
    {
        numngb = ngb_treefind_pairs(&pos[0], h_i, &startnode);

        for(n = 0; n < numngb; n++)
        {
            j = Ngblist[n];

            dx = pos[0] - P[j].Pos[0];
            dy = pos[1] - P[j].Pos[1];
            dz = pos[2] - P[j].Pos[2];
#ifdef MORRIS97VISC
            alpha_visc_j = SphP[j].Alpha;
#endif
            r2 = dx * dx + dy * dy + dz * dz;
            h_j = SphP[j].Hsml;
            if(r2 < h_i2 || r2 < h_j * h_j)
            {
                r = sqrt(r2);
                if(r > const_0)
                {
                    p_over_rho2_j = SphP[j].Pressure / (SphP[j].Density * SphP[j].Density);
                    soundspeed_j = sqrt(const_GAMMA * p_over_rho2_j * SphP[j].Density);
                    dvx = vel[0] - SphP[j].VelPred[0];
                    dvy = vel[1] - SphP[j].VelPred[1];
                    dvz = vel[2] - SphP[j].VelPred[2];
                    vdotr = dx * dvx + dy * dvy + dz * dvz;

                    if(All.ComovingIntegrationOn)
                        vdotr2 = vdotr + hubble_a2 * r2;
                    else
                        vdotr2 = vdotr;

                    if(r2 < h_i2)
                    {
                        hinv = const_1 / h_i;
#ifndef  TWODIMS
                        hinv4 = hinv * hinv * hinv * hinv;
#else
                        hinv4 = hinv * hinv * hinv / boxSize_Z;
#endif
                        u = r * hinv;
                        if(u < const_0_5)
                            dwk_i = hinv4 * u * (KERNEL_COEFF_3 * u - KERNEL_COEFF_4);
                        else
                            dwk_i = hinv4 * KERNEL_COEFF_6 * (const_1 - u) * (const_1 - u);
                    }
                    else
                    {
                        dwk_i = const_0;
                    }

                    if(r2 < h_j * h_j)
                    {
                        hinv = const_1 / h_j;
#ifndef  TWODIMS
                        hinv4 = hinv * hinv * hinv * hinv;
#else
                        hinv4 = hinv * hinv * hinv / boxSize_Z;
#endif
                        u = r * hinv;
                        if(u < const_0_5)
                            dwk_j = hinv4 * u * (KERNEL_COEFF_3 * u - KERNEL_COEFF_4);
                        else
                            dwk_j = hinv4 * KERNEL_COEFF_6 * (const_1 - u) * (const_1 - u);
                    }
                    else
                    {
                        dwk_j = const_0;
                    }

                    if(soundspeed_i + soundspeed_j > maxSignalVel)
                        maxSignalVel = soundspeed_i + soundspeed_j;

                    if(vdotr2 < const_0)	/* ... artificial viscosity */
                    {
#ifndef MONAGHAN83VISC
                        mu_ij = fac_mu * vdotr2 / r;	/* note: this is negative! */
#else
                        h_ij = const_0_5 * (h_i + h_j);
                        mu_ij = fac_mu * h_ij * vdotr2 / (r2 + const_0_0001 * h_ij * h_ij);
#endif
                        vsig = soundspeed_i + soundspeed_j - const_3 * mu_ij;



                        if(vsig > maxSignalVel)
                            maxSignalVel = vsig;

                        rho_ij = const_0_5 * (rho + SphP[j].Density);

#ifdef MORRIS97VISC
                        visc = 0.25 * (alpha_visc + alpha_visc_j) * vsig * (-mu_ij) / rho_ij;
#else

                        f2 = fabs(SphP[j].DivVel) / (fabs(SphP[j].DivVel) + SphP[j].CurlVel + const_0_0001 * soundspeed_j / fac_mu / SphP[j].Hsml);

#ifndef MONAGHAN83VISC
                        visc = const_0_25 * All.ArtBulkViscConst * vsig * (-mu_ij) / rho_ij * (f1 + f2);
#else
                        soundspeed_ij = (soundspeed_i+soundspeed_j) * 0.5;

                        visc = ((-All.ArtBulkViscConst) * soundspeed_ij * mu_ij + All.ArtBulkViscBeta * mu_ij * mu_ij) / rho_ij;
#endif //MONAGHAN83VISC

#endif //MORRIS97VISC

                        /* .... end artificial viscosity evaluation */
#ifndef NOVISCOSITYLIMITER
                        /* make sure that viscous acceleration is not too large */
                        dt = imax(timestep, (P[j].Ti_endstep - P[j].Ti_begstep)) * All.Timebase_interval;
                        if(dt > const_0 && (dwk_i + dwk_j) < const_0)
                        {
                            visc = dmin(visc, const_0_5 * fac_vsic_fix * vdotr2 /
                                        (const_0_5 * (mass + P[j].Mass) * (dwk_i + dwk_j) * r * dt));
                        }
#endif
                    }
                    else
                        visc = const_0;

                    p_over_rho2_j *= SphP[j].DhsmlDensityFactor;

                    hfc_visc = const_0_5 * P[j].Mass * visc * (dwk_i + dwk_j) / r;

                    hfc = hfc_visc + P[j].Mass * (p_over_rho2_i * dwk_i + p_over_rho2_j * dwk_j) / r;

                    acc[0] -= hfc * dx;
                    acc[1] -= hfc * dy;
                    acc[2] -= hfc * dz;
                    dtEntropy += const_0_5 * hfc_visc * vdotr2;
                }
            }
        }
    }
    while(startnode >= 0);

    /* Now collect the result at the right place */
    if(mode == 0)
    {
        for(k = 0; k < 3; k++)
            SphP[target].HydroAccel[k] = acc[k];
        SphP[target].DtEntropy = dtEntropy;
        SphP[target].MaxSignalVel = maxSignalVel;
    }
    else
    {
        HydroDataResult->set_init_Acc0(acc[0],target);
        HydroDataResult->set_init_Acc1(acc[1],target);
        HydroDataResult->set_init_Acc2(acc[2],target);
        HydroDataResult->set_init_DtEntropy(dtEntropy,target);
        HydroDataResult->set_init_MaxSignalVel(maxSignalVel,target);
    }
}




/*! This is a comparison kernel for a sort routine, which is used to group
 *  particles that are going to be exported to the same CPU.
 */
int gadgetmp2::hydro_compare_key(const void *a, const void *b)
{
    if(((hydrodata_in *) a)->read_Task() < (((hydrodata_in *) b)->read_Task()))
        return -1;
    if(((hydrodata_in *) a)->read_Task() > (((hydrodata_in *) b)->read_Task()))
        return +1;
    return 0;
}
