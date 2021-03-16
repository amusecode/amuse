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
#include <float.h>
#ifndef NOMPI
#include <mpi.h>
#endif

//#include "allvars.hpp"
#include "proto.hpp"

/*! \file gravtree.c
 *  \brief main driver routines for gravitational (short-range) force computation
 *
 *  This file contains the code for the gravitational force computation by
 *  means of the tree algorithm. To this end, a tree force is computed for
 *  all active local particles, and particles are exported to other
 *  processors if needed, where they can receive additional force
 *  contributions. If the TreePM algorithm is enabled, the force computed
 *  will only be the short-range part.
 */

/*! This function computes the gravitational forces for all active
 *  particles.  If needed, a new tree is constructed, otherwise the
 *  dynamically updated tree is used.  Particles are only exported to other
 *  processors when really needed, thereby allowing a good use of the
 *  communication buffer.
 */
void gadgetmp2::gravity_tree(void)
{
    long long ntot;
    int numnodes, nexportsum = 0;
    int i, j, iter = 0;
    int *numnodeslist, maxnumnodes, nexport, *numlist, *nrecv, *ndonelist;
    double tstart, tend, timetree = 0, timecommsumm = 0, timeimbalance = 0, sumimbalance;
    double ewaldcount;
    double costtotal, ewaldtot, *costtreelist, *ewaldlist;
    double maxt, sumt, *timetreelist, *timecommlist;
    double sumcomm;
    my_float fac, plb, plb_max;

#ifndef NOGRAVITY
    int *noffset, *nbuffer, *nsend, *nsend_local;
    long long ntotleft;
    int ndone, maxfill, ngrp;
    int k, place;
    int level, sendTask, recvTask;
    my_float ax, ay, az;
#ifndef NOMPI
    MPI_Status status;
#endif // NOMPI
#endif

    /* set new softening lengths */
    if(All.ComovingIntegrationOn)
        set_softenings();


    /* contruct tree if needed */
    tstart = second();
    if(TreeReconstructFlag)
    {
        if(ThisTask == 0)
            printf("Tree construction.\n");

        force_treebuild(NumPart);

        TreeReconstructFlag = 0;

        if(ThisTask == 0)
            printf("Tree construction done.\n");
    }
    tend = second();
    All.CPU_TreeConstruction += timediff(tstart, tend);

    costtotal =0;
    ewaldcount = 0;

    /* Note: 'NumForceUpdate' has already been determined in find_next_sync_point_and_drift() */
    numlist =  new int[NTask * NTask]; // just NTask is needed!
#ifndef NOMPI
    MPI_Allgather(&NumForceUpdate, 1, MPI_INT, numlist, 1, MPI_INT, GADGET_WORLD);
#else
    numlist[0] = NumForceUpdate;
#endif
    for(i = 0, ntot = 0; i < NTask; i++)
        ntot += numlist[i]; // sum all particles to be transfered
    free(numlist);


#ifndef NOGRAVITY
    if(ThisTask == 0)
        printf("Begin tree force.\n");

    noffset =  new int[NTask];	/* offsets of bunches in common list */
    nbuffer =  new int[NTask];
    nsend_local =  new int[NTask]; //particles to send to specific Task
    nsend =  new int[NTask * NTask];
    ndonelist =  new int[NTask];

    i = 0;			/* beginn with this index */
    ntotleft = ntot;		/* particles left for all tasks together */

    while(ntotleft > 0)
    {
        iter++;

        for(j = 0; j < NTask; j++)
            nsend_local[j] = 0;

        /* do local particles and prepare export list */
        tstart = second();
        for(nexport = 0, ndone = 0; i < NumPart && nexport < All.BunchSizeForce - NTask; i++)
            if(P[i].Ti_endstep == All.Ti_Current)
            {
                ndone++;

                for(j = 0; j < NTask; j++)
                    Exportflag[j] = 0;
#ifndef PMGRID
                costtotal += force_treeevaluate(i, 0, &ewaldcount); // decide to export to which task
#else
                costtotal += force_treeevaluate_shortrange(i, 0);
#endif
                for(j = 0; j < NTask; j++)
                {
                    if(Exportflag[j] != 0) //collect data to export into intermediate buffer
                    {
                        /*		    for(k = 0; k < 3; k++)
                         *		      GravDataGet[nexport].u.Pos[k] = P[i].Pos[k];
                         */
                        //GravDataGet[nexport].u0 = P[i].Pos[0];
                        GravDataGet->set_init_u0(P[i].Pos[0],nexport);
                        //GravDataGet[nexport].u1 = P[i].Pos[1];
                        GravDataGet->set_init_u1(P[i].Pos[1],nexport);
                        //GravDataGet[nexport].u2 = P[i].Pos[2];
                        GravDataGet->set_init_u2(P[i].Pos[2],nexport);
#ifdef UNEQUALSOFTENINGS
                        //GravDataGet[nexport].Type = P[i].Type;
                        GravDataGet->set_Type(P[i].Type,nexport);
#ifdef ADAPTIVE_GRAVSOFT_FORGAS
                        if(P[i].Type == 0)
                            //GravDataGet[nexport].Soft = SphP[i].Hsml;
                            GravDataGet->set_init_Soft(SphP[i].Hsml,nexport);
#endif
#endif
                        //GravDataGet[nexport].OldAcc = P[i].OldAcc;
                        GravDataGet->set_init_OldAcc(P[i].OldAcc,nexport);
                        GravDataIndexTable[nexport].Task = j;
                        GravDataIndexTable[nexport].Index = i;
                        GravDataIndexTable[nexport].SortIndex = nexport;
                        nexport++;
                        nexportsum++;
                        nsend_local[j]++;
                    }
                }
            }
        tend = second();
        timetree += timediff(tstart, tend);

        qsort(GravDataIndexTable, nexport, sizeof(struct gravdata_index), grav_tree_compare_key); //generate pattern to setup real buffer

        for(j = 0; j < nexport; j++)
            gravdata_in::lcopy(GravDataIn, j, GravDataGet, GravDataIndexTable[j].SortIndex);
        //               GravDataIn[j] = GravDataGet[GravDataIndexTable[j].SortIndex];   //this builts the real transfer buffer

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

#ifndef NOMPI
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
                if(maxfill >= All.BunchSizeForce)
                    break;

                sendTask = ThisTask;
                recvTask = ThisTask ^ ngrp;

                if(recvTask < NTask)
                {
                    if(nsend[ThisTask * NTask + recvTask] > 0 || nsend[recvTask * NTask + ThisTask] > 0)
                    {
                        /* get the particles */
                        MPI_Sendrecv(GravDataIn->get_buff_start(noffset[recvTask]),
                                     nsend_local[recvTask] * gravdata_in::get_size(), MPI_BYTE,
                                     recvTask, TAG_GRAV_A,
                                     GravDataGet->get_buff_start(nbuffer[ThisTask]),
                                     nsend[recvTask * NTask + ThisTask] * gravdata_in::get_size(), MPI_BYTE,
                                recvTask, TAG_GRAV_A, GADGET_WORLD, &status);
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
            {
#ifndef PMGRID
                costtotal += force_treeevaluate(j, 1, &ewaldcount);
#else
                costtotal += force_treeevaluate_shortrange(j, 1);
#endif
            }
            tend = second();
            timetree += timediff(tstart, tend);

            tstart = second();
            MPI_Barrier(GADGET_WORLD);
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
                if(maxfill >= All.BunchSizeForce)
                    break;

                sendTask = ThisTask;
                recvTask = ThisTask ^ ngrp;
                if(recvTask < NTask)
                {
                    if(nsend[ThisTask * NTask + recvTask] > 0 || nsend[recvTask * NTask + ThisTask] > 0)
                    {
                        /* send the results */
                        MPI_Sendrecv(GravDataResult->get_buff_start(nbuffer[ThisTask]),
                                     nsend[recvTask * NTask + ThisTask] * gravdata_in::get_size(),
                                MPI_BYTE, recvTask, TAG_GRAV_B,
                                GravDataOut->get_buff_start(noffset[recvTask]),
                                nsend_local[recvTask] * gravdata_in::get_size(),
                                MPI_BYTE, recvTask, TAG_GRAV_B, GADGET_WORLD, &status);

                        /* add the result to the particles */
                        for(j = 0; j < nsend_local[recvTask]; j++)
                        {
                            place = GravDataIndexTable[noffset[recvTask] + j].Index;

                            /*			  for(k = 0; k < 3; k++)
                                 *			    P[place].GravAccel[k] += GravDataOut[j + noffset[recvTask]].u.Acc[k];
                                 */
                            //P[place].GravAccel[0] += GravDataOut[j + noffset[recvTask]].u0;
                            P[place].GravAccel[0] += GravDataOut->read_re_init_u0(j + noffset[recvTask]);
                            //P[place].GravAccel[1] += GravDataOut[j + noffset[recvTask]].u1;
                            P[place].GravAccel[1] += GravDataOut->read_re_init_u1(j + noffset[recvTask]);
                            //P[place].GravAccel[2] += GravDataOut[j + noffset[recvTask]].u2;
                            P[place].GravAccel[2] += GravDataOut->read_re_init_u2(j + noffset[recvTask]);

                            //P[place].GravCost += GravDataOut[j + noffset[recvTask]].Ninteractions;
                            P[place].GravCost += GravDataOut->read_Ninteractions(j + noffset[recvTask]);
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
#endif
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

    /* now add things for comoving integration */

#ifndef PERIODIC
#ifndef PMGRID
    if(All.ComovingIntegrationOn)
    {
        fac = const_0_5 * All.Hubble * All.Hubble * All.Omega0 / All.G;

        for(i = 0; i < NumPart; i++)
            if(P[i].Ti_endstep == All.Ti_Current)
                for(j = 0; j < 3; j++)
                    P[i].GravAccel[j] += fac * P[i].Pos[j];
    }
#endif
#endif

    for(i = 0; i < NumPart; i++)
        if(P[i].Ti_endstep == All.Ti_Current)
        {

            ax = P[i].GravAccel[0];
            ay = P[i].GravAccel[1];
            az = P[i].GravAccel[2];
            P[i].OldAcc = sqrt(ax * ax + ay * ay + az * az);
        }


    if(All.TypeOfOpeningCriterion == 1)
        All.ErrTolTheta = 0;	/* This will switch to the relative opening criterion for the following force computations */

    /*  muliply by G */
    for(i = 0; i < NumPart; i++)
        if(P[i].Ti_endstep == All.Ti_Current)
            for(j = 0; j < 3; j++)
                P[i].GravAccel[j] *= All.G;


    /* Finally, the following factor allows a computation of a cosmological simulation
                     *     with vacuum energy in physical coordinates */
#ifndef PERIODIC
#ifndef PMGRID
    if(All.ComovingIntegrationOn == 0)
    {
        fac = All.OmegaLambda * All.Hubble * All.Hubble;

        for(i = 0; i < NumPart; i++)
            if(P[i].Ti_endstep == All.Ti_Current)
                for(j = 0; j < 3; j++)
                    P[i].GravAccel[j] += fac * P[i].Pos[j];
    }
#endif
#endif

    if(ThisTask == 0)
        printf("tree is done.\n");

#else /* gravity is switched off */

    for(i = 0; i < NumPart; i++)
        if(P[i].Ti_endstep == All.Ti_Current)
            for(j = 0; j < 3; j++)
                P[i].GravAccel[j] = 0;

#endif




    /* Now the force computation is finished */

    /*  gather some diagnostic information */

    timetreelist =  new double[NTask];
    timecommlist = new double[NTask];
    costtreelist = new double[NTask];
    numnodeslist = new int[NTask];
    ewaldlist = new double[NTask];
    nrecv =  new int[NTask];

    numnodes = Numnodestree;
#ifndef NOMPI
    MPI_Gather(&costtotal, 1, MPI_DOUBLE, costtreelist, 1, MPI_DOUBLE, 0, GADGET_WORLD);
    MPI_Gather(&numnodes, 1, MPI_INT, numnodeslist, 1, MPI_INT, 0, GADGET_WORLD);
    MPI_Gather(&timetree, 1, MPI_DOUBLE, timetreelist, 1, MPI_DOUBLE, 0, GADGET_WORLD);
    MPI_Gather(&timecommsumm, 1, MPI_DOUBLE, timecommlist, 1, MPI_DOUBLE, 0, GADGET_WORLD);
    MPI_Gather(&NumPart, 1, MPI_INT, nrecv, 1, MPI_INT, 0, GADGET_WORLD);
    MPI_Gather(&ewaldcount, 1, MPI_DOUBLE, ewaldlist, 1, MPI_DOUBLE, 0, GADGET_WORLD);
    MPI_Reduce(&nexportsum, &nexport, 1, MPI_INT, MPI_SUM, 0, GADGET_WORLD);
    MPI_Reduce(&timeimbalance, &sumimbalance, 1, MPI_DOUBLE, MPI_SUM, 0, GADGET_WORLD);
#else
    costtreelist[0] = costtotal;
    numnodeslist[0] = numnodes;
    timetreelist[0] = timetree;
    timecommlist[0] = timecommsumm;
    nrecv[0] = NumPart;
    ewaldlist[0] = ewaldcount;
    nexport  = nexportsum;
    sumimbalance = timeimbalance;
#endif
    if(ThisTask == 0)
    {
        All.TotNumOfForces += ntot;

        fprintf(FdTimings, "Step= %d  t= %g  dt= %g \n", All.NumCurrentTiStep, All.Time, All.TimeStep);
        fprintf(FdTimings, "Nf= %d%09d  total-Nf= %d%09d  ex-frac= %g  iter= %d\n",
                (int) (ntot / 1000000000), (int) (ntot % 1000000000),
                (int) (All.TotNumOfForces / 1000000000), (int) (All.TotNumOfForces % 1000000000),
                (nexport / ((my_float) ntot)).toDouble(), iter);
        /* note: on Linux, the 8-byte integer could be printed with the format identifier "%qd", but doesn't work on AIX */

        fac = NTask / ((my_float) All.TotNumPart);

        for(i = 0, maxt = timetreelist[0], sumt = 0, plb_max = const_0,
            maxnumnodes = 0, costtotal = 0, sumcomm = 0, ewaldtot = 0; i < NTask; i++)
        {
            costtotal += costtreelist[i];

            sumcomm += timecommlist[i];

            if(maxt < timetreelist[i])
                maxt = timetreelist[i];
            sumt += timetreelist[i];

            plb = nrecv[i] * fac;

            if(plb > plb_max)
                plb_max = plb;

            if(numnodeslist[i] > maxnumnodes)
                maxnumnodes = numnodeslist[i];

            ewaldtot += ewaldlist[i];
        }
        fprintf(FdTimings, "work-load balance: %g  max=%g avg=%g PE0=%g\n",
                (maxt / (sumt / NTask)), maxt, (sumt / NTask), timetreelist[0]);
        fprintf(FdTimings, "particle-load balance: %g\n", plb_max.toDouble());
        fprintf(FdTimings, "max. nodes: %d, filled: %g\n", maxnumnodes,
                (maxnumnodes / (All.TreeAllocFactor * All.MaxPart)));
        fprintf(FdTimings, "part/sec=%g | %g  ia/part=%g (%g)\n", (ntot / (sumt + 1.0e-20)),
                (ntot / (maxt * NTask)), ( (costtotal)) / ntot), ((( ewaldtot) / ntot));
        fprintf(FdTimings, "\n");

        fflush(FdTimings);

        All.CPU_TreeWalk += sumt / NTask;
        All.CPU_Imbalance += sumimbalance / NTask;
        All.CPU_CommSum += sumcomm / NTask;
    }

    free(nrecv);
    free(ewaldlist);
    free(numnodeslist);
    free(costtreelist);
    free(timecommlist);
    free(timetreelist);
}



/*! This function sets the (comoving) softening length of all particle
 *  types in the table All.SofteningTable[...].  We check that the physical
 *  softening length is bounded by the Softening-MaxPhys values.
 */
void gadgetmp2::set_softenings(void)
{
    int i;

    if(All.ComovingIntegrationOn)
    {
        if(All.SofteningGas * All.Time > All.SofteningGasMaxPhys)
            All.SofteningTable[0] = All.SofteningGasMaxPhys / All.Time;
        else
            All.SofteningTable[0] = All.SofteningGas;

        if(All.SofteningHalo * All.Time > All.SofteningHaloMaxPhys)
            All.SofteningTable[1] = All.SofteningHaloMaxPhys / All.Time;
        else
            All.SofteningTable[1] = All.SofteningHalo;

        if(All.SofteningDisk * All.Time > All.SofteningDiskMaxPhys)
            All.SofteningTable[2] = All.SofteningDiskMaxPhys / All.Time;
        else
            All.SofteningTable[2] = All.SofteningDisk;

        if(All.SofteningBulge * All.Time > All.SofteningBulgeMaxPhys)
            All.SofteningTable[3] = All.SofteningBulgeMaxPhys / All.Time;
        else
            All.SofteningTable[3] = All.SofteningBulge;

        if(All.SofteningStars * All.Time > All.SofteningStarsMaxPhys)
            All.SofteningTable[4] = All.SofteningStarsMaxPhys / All.Time;
        else
            All.SofteningTable[4] = All.SofteningStars;

        if(All.SofteningBndry * All.Time > All.SofteningBndryMaxPhys)
            All.SofteningTable[5] = All.SofteningBndryMaxPhys / All.Time;
        else
            All.SofteningTable[5] = All.SofteningBndry;
    }
    else
    {
        All.SofteningTable[0] = All.SofteningGas;
        All.SofteningTable[1] = All.SofteningHalo;
        All.SofteningTable[2] = All.SofteningDisk;
        All.SofteningTable[3] = All.SofteningBulge;
        All.SofteningTable[4] = All.SofteningStars;
        All.SofteningTable[5] = All.SofteningBndry;
    }

    for(i = 0; i < 6; i++)
        All.ForceSoftening[i] = 2.8 * All.SofteningTable[i];

    All.MinGasHsml = All.MinGasHsmlFractional * All.ForceSoftening[0];
}


/*! This function is used as a comparison kernel in a sort routine. It is
 *  used to group particles in the communication buffer that are going to
 *  be sent to the same CPU.
 */
int gadgetmp2::grav_tree_compare_key(const void *a, const void *b)
{
    if(((struct gravdata_index *) a)->Task < (((struct gravdata_index *) b)->Task))
        return -1;

    if(((struct gravdata_index *) a)->Task > (((struct gravdata_index *) b)->Task))
        return +1;

    return 0;
}
