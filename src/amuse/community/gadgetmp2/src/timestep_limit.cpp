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
/* ################################################################################## */
/* ###                                                                            ### */
/* ###                              timestep_limit.c                              ### */
/* ###                                                                            ### */
/* ###   Original: timestep.c (public version of Gadget 2)                        ### */
/* ###   Author: Volker Springel                                                  ### */
/* ###                                                                            ### */
/* ###   Modified: February 2011                                                  ### */
/* ###   Author: Fabrice Durier                                                   ### */
/* ###                                                                            ### */
/* ################################################################################## */


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#ifndef NOMPI
#include <mpi.h>
#endif
#include <time.h>
#include <assert.h>

//#include "allvars.hpp"
#include "proto.hpp"

/*! \file timestep_limit.c
 *
 *  This file is a modified version of timestep.c, which include the routines
 *  for 'kicking' particles in momentum space and assigning new timesteps.
 *
 *  This file now includes the Timestep Limiter routine for communication of
 *  the new timestep between neighbouring SPH particles
 */


static my_float fac1, fac2, fac3, hubble_a, atime, a3inv;
static my_float dt_displacement = 0;

static int NDone;
static long long NTotDone;


/*! This function advances the system in momentum space, i.e. it does apply
 *  the 'kick' operation after the forces have been computed. Additionally, it
 *  assigns new timesteps to particles following the Limiter criterion.
 *  At start-up, a half-timestep is carried out, as well as at the end of the simulation.
 *  In between, the half-step kick that ends the previous timestep and
 *  the half-step kick for the new timestep are combined into one operation.
 */

void gadgetmp2::advance_and_find_timesteps(void)
{
    long long ntot, ntotleft;
    int *noffset, *nbuffer, *nsend, *nsend_local, *numlist, *ndonelist;
    int n, ndone, maxfill, source;
    int level, ngrp, sendTask, recvTask, place, nexport;
    double t0, t1;
    my_float sumt = 0, sumimbalance = 0, sumcomm = 0;
    double timecomp = 0, timeimbalance = 0, timecommsumm = 0;
    #ifndef NOMPI
    MPI_Status status;
    #endif

    int CptLimit = 0;
    int shrinkcount = 0, shrinktot = 0;

    int i, j, no, ti_step, ti_min, tend, tstart;
    my_float dt_entr, dt_entr2, dt_gravkick, dt_hydrokick, dt_gravkick2, dt_hydrokick2;
    my_float minentropy, aphys;
    my_float dv[3];

    #ifdef FLEXSTEPS
    int ti_grp;
    #endif
    #if defined(PSEUDOSYMMETRIC) && !defined(FLEXSTEPS)
    my_float apred, prob;
    int ti_step2;
    #endif
    #ifdef MAKEGLASS
    my_float disp, dispmax, globmax, dmean, fac, disp2sum, globdisp2sum;
    #endif
    #ifdef MORRIS97VISC
    my_float soundspeed, tau, f_fac;
    #endif

    noffset = new int[NTask];  /* offsets of bunches in common list */
    nbuffer = new int[NTask];
    nsend_local = new int[NTask];
    nsend = new int[NTask * NTask];
    ndonelist = new int[NTask];

    if(ThisTask == 0)
    {
        fprintf(stdout,"Start time-stepping evaluation ...\n");
        fflush(stdout);
    }

    if(All.ComovingIntegrationOn)
    {
        fac1 = 1 / (All.Time * All.Time);
        fac2 = 1 / std::pow(All.Time, 3 * GAMMA - 2);
        fac3 = std::pow(All.Time, 3 * (1 - GAMMA) / 2.0);
        hubble_a = All.Omega0 / (All.Time * All.Time * All.Time)
        + (1 - All.Omega0 - All.OmegaLambda) / (All.Time * All.Time) + All.OmegaLambda;

        hubble_a = All.Hubble * sqrt(hubble_a);
        a3inv = 1 / (All.Time * All.Time * All.Time);
        atime = All.Time;
    }
    else
        fac1 = fac2 = fac3 = hubble_a = a3inv = atime = 1;

    if(Flag_FullStep || dt_displacement == 0)
        find_dt_displacement_constraint(hubble_a * atime * atime);

    #ifdef MAKEGLASS
    for(i = 0, dispmax = 0, disp2sum = 0; i < NumPart; i++)
    {
        for(j = 0; j < 3; j++)
        {
            P[i].GravPM[j] *= -1;
            P[i].GravAccel[j] *= -1;
            P[i].GravAccel[j] += P[i].GravPM[j];
            P[i].GravPM[j] = 0;
        }

        disp = sqrt(P[i].GravAccel[0] * P[i].GravAccel[0] +
        P[i].GravAccel[1] * P[i].GravAccel[1] + P[i].GravAccel[2] * P[i].GravAccel[2]);

        disp *= 2.0 / (3 * All.Hubble * All.Hubble);

        disp2sum += disp * disp;

        if(disp > dispmax)
            dispmax = disp;
    }

    #ifndef NOMPI
    //MPI_Allreduce(&dispmax, &globmax, 1, MPI_DOUBLE, MPI_MAX, GADGET_WORLD);
    globmax=mpi_all_max(dispmax);
    //MPI_Allreduce(&disp2sum, &globdisp2sum, 1, MPI_DOUBLE, MPI_SUM, GADGET_WORLD);
    globdisp2sum=mpi_all_sum(disp2sum);
    #else
    globmax = dispmax;
    globdisp2sum = disp2sum;
    #endif
    dmean = pow(P[0].Mass / (All.Omega0 * 3 * All.Hubble * All.Hubble / (8 * M_PI * All.G)), 1.0 / 3);

    if(globmax > dmean)
        fac = dmean / globmax;
    else
        fac = 1.0;

    if(ThisTask == 0)
    {
        printf("\nglass-making:  dmean= %g  global disp-maximum= %g  rms= %g\n\n",
               dmean.toDouble(), globmax.toDouble(), sqrt(globdisp2sum / All.TotNumPart).toDouble());
        fflush(stdout);
    }

    for(i = 0, dispmax = 0; i < NumPart; i++)
    {
        for(j = 0; j < 3; j++)
        {
            P[i].Vel[j] = 0;
            P[i].Pos[j] += fac * P[i].GravAccel[j] * 2.0 / (3 * All.Hubble * All.Hubble);
            P[i].GravAccel[j] = 0;
        }
    }
    #endif



    /* Now assign new timesteps and kick */

    #ifdef FLEXSTEPS
    if((All.Ti_Current % (4 * All.PresentMinStep)) == 0)
        if(All.PresentMinStep < TIMEBASE)
            All.PresentMinStep *= 2;

        for(i = 0; i < NumPart; i++)
        {
            if(P[i].Ti_endstep == All.Ti_Current)
            {
                ti_step = get_timestep(i, &aphys, 0);

                /* make it a power 2 subdivision */
                ti_min = TIMEBASE;
                while(ti_min > ti_step)
                    ti_min >>= 1;
                ti_step = ti_min;

                if(ti_step < All.PresentMinStep)
                    All.PresentMinStep = ti_step;
            }
        }

        ti_step = All.PresentMinStep;

        #ifndef NOMPI
        MPI_Allreduce(&ti_step, &All.PresentMinStep, 1, MPI_INT, MPI_MIN, GADGET_WORLD);
        #else
        All.PresentMinStep = ti_step;
        #endif
        #endif


        for(i = 0; i < NumPart; i++)
        {
            if(P[i].Ti_endstep == All.Ti_Current)
            {
                ti_step = get_timestep(i, &aphys, 0);

                /* make it a power 2 subdivision */
                ti_min = TIMEBASE;
                while(ti_min > ti_step)
                    ti_min >>= 1;
                ti_step = ti_min;

                #ifdef FLEXSTEPS
                ti_grp = P[i].FlexStepGrp % All.PresentMaxStep;
                ti_grp = (ti_grp / All.PresentMinStep) * All.PresentMinStep;
                ti_step = ((P[i].Ti_endstep + ti_grp + ti_step) / ti_step) * ti_step - (P[i].Ti_endstep + ti_grp);
                #else

                #ifdef PSEUDOSYMMETRIC
                if(P[i].Type != 0)
                {
                    if(P[i].Ti_endstep > P[i].Ti_begstep)
                    {
                        apred = aphys + ((aphys - P[i].AphysOld) / (P[i].Ti_endstep - P[i].Ti_begstep)) * ti_step;
                        if(fabs(apred - aphys) < 0.5 * aphys)
                        {
                            ti_step2 = get_timestep(i, &apred, -1);
                            ti_min = TIMEBASE;
                            while(ti_min > ti_step2)
                                ti_min >>= 1;
                            ti_step2 = ti_min;

                            if(ti_step2 < ti_step)
                            {
                                get_timestep(i, &apred, ti_step);
                                prob =
                                ((apred - aphys) / (aphys - P[i].AphysOld) * (P[i].Ti_endstep -
                                P[i].Ti_begstep)) / ti_step;
                                if(prob < get_random_number(P[i].ID))
                                    ti_step /= 2;
                            }
                            else if(ti_step2 > ti_step)
                            {
                                get_timestep(i, &apred, 2 * ti_step);
                                prob =
                                ((apred - aphys) / (aphys - P[i].AphysOld) * (P[i].Ti_endstep -
                                P[i].Ti_begstep)) / ti_step;
                                if(prob < get_random_number(P[i].ID + 1))
                                    ti_step *= 2;
                            }
                        }
                    }
                    P[i].AphysOld = aphys;
                }
                #endif

                #ifdef SYNCHRONIZATION
                if(ti_step > (P[i].Ti_endstep - P[i].Ti_begstep))	/* timestep wants to increase */
                {
                    if(((TIMEBASE - P[i].Ti_endstep) % ti_step) > 0)
                        ti_step = P[i].Ti_endstep - P[i].Ti_begstep;	/* leave at old step */
                }
                #endif
                #endif /* end of FLEXSTEPS */

                if(All.Ti_Current == TIMEBASE)	/* we here finish the last timestep. */
                    ti_step = 0;

                if((TIMEBASE - All.Ti_Current) < ti_step)	/* check that we don't run beyond the end */
                    ti_step = TIMEBASE - All.Ti_Current;

                #ifdef TIMESTEP_LIMITER
                /*  ------------------------------------------------------------------------------------------------------------------  */
                /* 	  Main change in the time intregration module:  */
                /* 	  loop over all active gas particles to communicate  */
                /* 	  their NEXT timestep to the neighbours */
                /*  ------------------------------------------------------------------------------------------------------------------  */

                P[i].Ti_sizestep = ti_step;
            }
        }
        t0 = second();
        NumSphUpdate = 0;

        if(All.NumCurrentTiStep > 0)
        {
            for(n = 0; n < N_gas; n++)
                if(P[n].Type == 0)
                {
                    if(P[n].Ti_endstep == All.Ti_Current)
                        NumSphUpdate++;
                }
        }
        t1 = second();
        timecomp += timediff(t0, t1);

        t0 = second();
        numlist = new int[NTask * NTask];

        #ifndef NOMPI
        MPI_Allgather(&NumSphUpdate, 1, MPI_INT, numlist, 1, MPI_INT, GADGET_WORLD);
        #else
        numlist[0] = NumSphUpdate;
        #endif
        for(i = 0, ntot = 0; i < NTask; i++)
            ntot += numlist[i];
        free(numlist);
        t1 = second();
        timeimbalance += timediff(t0, t1);

        if(ThisTask == 0)
        {
            fprintf(stdout,"\n      ---> number of timestep active gas particles:    %6lld \n", ntot);
            fflush(stdout);
        }

        if(ntot)
            NTotDone = 1;
        else
            NTotDone = 0;

        /* loop for time-step limiter */
        while(NTotDone > 0)
        {
            NDone = 0;

            i = 0;                        /* first particle for this task */
            ntotleft = ntot;              /* particles left for all tasks together */

            while(ntotleft > 0)
            {
                for(j = 0; j < NTask; j++)
                    nsend_local[j] = 0;

                /* do local particles and prepare export list */
                t0 = second();
                for(nexport = 0, ndone = 0; i < N_gas && nexport < All.BunchSizeTime - NTask; i++)
                    if(P[i].Type == 0)
                    {
                        if(P[i].Ti_endstep == All.Ti_Current)
                        {
                            ndone++;

                            for(j = 0; j < NTask; j++)
                                Exportflag[j] = 0;

                            /* -------------------------------- */
                            NDone += time_limiter_evaluate(i, 0);
                            /* -------------------------------- */

                            for(j = 0; j < NTask; j++)
                            {
                                if(Exportflag[j])
                                {
                                    //TimeDataIn[nexport].Pos[0] =    P[i].Pos[0];
                                    TimeDataIn->set_init_Pos0(P[i].Pos[0],nexport);
                                    //TimeDataIn[nexport].Pos[1] =    P[i].Pos[1];
                                    TimeDataIn->set_init_Pos1(P[i].Pos[1],nexport);
                                    //TimeDataIn[nexport].Pos[2] =    P[i].Pos[2];
                                    TimeDataIn->set_init_Pos2(P[i].Pos[2],nexport);
                                    //TimeDataIn[nexport].Hsml   = SphP[i].Hsml;
                                    TimeDataIn->set_init_Hsml(SphP[i].Hsml,nexport);
                                    //TimeDataIn[nexport].Size   =    P[i].Ti_sizestep;
                                    TimeDataIn->set_Size(P[i].Ti_sizestep,nexport);
                                    //TimeDataIn[nexport].Begin  =    P[i].Ti_endstep;
                                    TimeDataIn->set_Begin(P[i].Ti_endstep,nexport);
                                    //TimeDataIn[nexport].Index  = i;
                                    TimeDataIn->set_Index(i,nexport);
                                    //TimeDataIn[nexport].Task   = j;
                                    TimeDataIn->set_Task(j,nexport);
                                    nexport++;
                                    nsend_local[j]++;
                                }
                            }
                        }
                    }
                    t1 = second();
                    timecomp += timediff(t0, t1);

                    qsort(TimeDataIn, nexport, timedata_in::get_size(), time_compare_key);

                    for(j = 1, noffset[0] = 0; j < NTask; j++)
                        noffset[j] = noffset[j - 1] + nsend_local[j - 1];

                    t0 = second();

                    #ifndef NOMPI
                    MPI_Allgather(nsend_local, NTask, MPI_INT, nsend, NTask, MPI_INT, GADGET_WORLD);
                    #else
                    nsend[0] = nsend_local[0];
                    #endif
                    t1 = second();
                    timeimbalance += timediff(t0, t1);

                    /* now do the particles that need to be exported */

                    for(level = 1; level < (1 << PTask); level++)
                    {
                        t0 = second();
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
                            if(maxfill >= All.BunchSizeTime)
                                break;

                            sendTask = ThisTask;
                            recvTask = ThisTask ^ ngrp;

                            #ifndef NOMPI
                            if(recvTask < NTask)
                            {
                                if(nsend[ThisTask * NTask + recvTask] > 0 || nsend[recvTask * NTask + ThisTask] > 0)
                                {
                                    /* get the particles */
                                    MPI_Sendrecv(TimeDataIn->get_buff_start(noffset[recvTask]),
                                                 nsend_local[recvTask] * timedata_in::get_size(), MPI_BYTE,
                                                 recvTask, 0,
                                                 TimeDataGet->get_buff_start(nbuffer[ThisTask]),
                                                 nsend[recvTask * NTask + ThisTask] * timedata_in::get_size(),
                                                 MPI_BYTE, recvTask, 0, GADGET_WORLD, &status);
                                }
                            }
                            #endif

                            for(j = 0; j < NTask; j++)
                                if((j ^ ngrp) < NTask)
                                    nbuffer[j] += nsend[(j ^ ngrp) * NTask + j];
                        }
                        t1 = second();
                        timecommsumm += timediff(t0, t1);

                        t0 = second();
                        /* --------------------------------- */
                        for(j = 0; j < nbuffer[ThisTask]; j++)
                            NDone += time_limiter_evaluate(j, 1);
                        /* --------------------------------- */
                        t1 = second();
                        timecomp += timediff(t0, t1);

                        /* do a block to explicitly measure imbalance */
                        t0 = second();
                        #ifndef NOMPI
                        MPI_Barrier(GADGET_WORLD);
                        #endif
                        t1 = second();
                        timeimbalance += timediff(t0, t1);

                        level = ngrp - 1;
                    }

                    t0 = second();

                    #ifndef NOMPI
                    MPI_Allgather(&ndone, 1, MPI_INT, ndonelist, 1, MPI_INT, GADGET_WORLD);
                    #else
                    ndonelist[0] = ndone;
                    #endif
                    for(j = 0; j < NTask; j++)
                        ntotleft -= ndonelist[j];
                    t1 = second();
                    timeimbalance += timediff(t0, t1);
            }

            t0 = second();
            numlist = new int[NTask * NTask];

            #ifndef NOMPI
            MPI_Allgather(&NDone, 1, MPI_INT, numlist, 1, MPI_INT, GADGET_WORLD);
            #else
            numlist[0] = NDone;
            #endif
            for(i = 0, NTotDone = 0; i < NTask; i++)
                NTotDone += numlist[i];
            free(numlist);
            t1 = second();
            timeimbalance += timediff(t0, t1);

            if(ThisTask == 0)
            {
                fprintf(stdout," %3d) ---> number of timestep shrinked gas neighbors:  %6lld \n", CptLimit++, NTotDone);
                fflush(stdout);
            }
        }

        /*   do final operations on results */
        t0 = second();
        for(i = 0; i < NumPart; i++)
        {
            if(P[i].Ti_endstep == All.Ti_Current)
            {
                ti_step =  P[i].Ti_sizestep;

                /*  ------------------------------------------------------------------------------------------------------------------  */
                #endif

                tstart  = (P[i].Ti_begstep + P[i].Ti_endstep) / 2;	/* midpoint of old step */
                tend    =  P[i].Ti_endstep + ti_step / 2;	        /* midpoint of new step */

                if(All.ComovingIntegrationOn)
                {
                    dt_entr = (tend - tstart) * All.Timebase_interval;
                    dt_entr2 = (tend - P[i].Ti_endstep) * All.Timebase_interval;
                    dt_gravkick = get_gravkick_factor(tstart, tend);
                    dt_hydrokick = get_hydrokick_factor(tstart, tend);
                    dt_gravkick2 = get_gravkick_factor(P[i].Ti_endstep, tend);
                    dt_hydrokick2 = get_hydrokick_factor(P[i].Ti_endstep, tend);
                }
                else
                {
                    dt_entr = dt_gravkick = dt_hydrokick = (tend - tstart) * All.Timebase_interval;
                    dt_gravkick2 = dt_hydrokick2 = dt_entr2 = (tend - P[i].Ti_endstep) * All.Timebase_interval;
                }

                P[i].Ti_begstep = P[i].Ti_endstep;
                P[i].Ti_endstep = P[i].Ti_begstep + ti_step;

                /* do the kick */

                for(j = 0; j < 3; j++)
                {
                    dv[j] = P[i].GravAccel[j] * dt_gravkick;
                    P[i].Vel[j] += dv[j];
                }

                if(P[i].Type == 0)	/* SPH stuff */
                {
                    #ifdef MORRIS97VISC

                    soundspeed  = sqrt(GAMMA * SphP[i].Pressure / SphP[i].Density);
                    f_fac = fabs(SphP[i].DivVel) / (fabs(SphP[i].DivVel) + SphP[i].CurlVel +
                    0.0001 * soundspeed / SphP[i].Hsml);
                    tau = 2.5 * SphP[i].Hsml / soundspeed;
                    SphP[i].DAlphaDt = f_fac*dmax(-SphP[i].DivVel, 0) * (2.0 - SphP[i].Alpha) - (SphP[i].Alpha - .1)/tau;

                    /*printf("f_fac = %g, tau = %g, Alpha = %g, DivVel = %g, DAlphaDt = %g \n",
                     *		f_fac.toDouble(),
                     *		tau.toDouble(),
                     *		SphP[i].Alpha.toDouble(),
                     *		SphP[i].DivVel.toDouble(),
                     *		SphP[i].DAlphaDt.toDouble());*/
                    /* change this to input the maximum the viscosity can get to. */
                    /*  SphP[i].DAlphaDt = f_fac*dmax(-SphP[i].DivVel, 0) * (All.ArtBulkViscConst - SphP[i].Alpha) - (SphP[i].Alpha - .1)/tau;*/
                    #endif
                    for(j = 0; j < 3; j++)
                    {
                        dv[j] += SphP[i].HydroAccel[j] * dt_hydrokick;
                        P[i].Vel[j] += SphP[i].HydroAccel[j] * dt_hydrokick;

                        SphP[i].VelPred[j] =
                        P[i].Vel[j] - dt_gravkick2 * P[i].GravAccel[j] - dt_hydrokick2 * SphP[i].HydroAccel[j];
                    }

                    /* In case of cooling, we prevent that the entropy (and
                     *	         hence temperature decreases by more than a factor 0.5 */

                    if(SphP[i].DtEntropy * dt_entr > -0.5 * SphP[i].Entropy)
                        SphP[i].Entropy += SphP[i].DtEntropy * dt_entr;
                    else
                        SphP[i].Entropy *= 0.5;

                    if(All.MinEgySpec)
                    {
                        minentropy = All.MinEgySpec * GAMMA_MINUS1 / pow(SphP[i].Density * a3inv, GAMMA_MINUS1);
                        if(SphP[i].Entropy < minentropy)
                        {
                            SphP[i].Entropy = minentropy;
                            SphP[i].DtEntropy = 0;
                        }
                    }

                    /* In case the timestep increases in the new step, we
                     *	         make sure that we do not 'overcool' when deriving
                     *	         predicted temperatures. The maximum timespan over
                     *	         which prediction can occur is ti_step/2, i.e. from
                     *	         the middle to the end of the current step */

                    dt_entr = ti_step / 2 * All.Timebase_interval;
                    if(SphP[i].Entropy + SphP[i].DtEntropy * dt_entr < 0.5 * SphP[i].Entropy)
                        SphP[i].DtEntropy = -0.5 * SphP[i].Entropy / dt_entr;
                }


                /* if tree is not going to be reconstructed, kick parent nodes dynamically.
                 */
                if(All.NumForcesSinceLastDomainDecomp < All.TotNumPart * All.TreeDomainUpdateFrequency)
                {
                    no = Father[i];
                    while(no >= 0)
                    {
                        for(j = 0; j < 3; j++)
                            Extnodes[no].vs[j] += dv[j] * P[i].Mass / Nodes[no].u_d_mass;

                        no = Nodes[no].u.d.father;
                    }
                }
            }
        }


        t1 = second();
        timecomp += timediff(t0, t1);

        free(ndonelist);
        free(nsend);
        free(nsend_local);
        free(nbuffer);
        free(noffset);

}



/*! This function normally (for flag==0) returns the maximum allowed timestep
 *  of a particle, expressed in terms of the integer mapping that is used to
 *  represent the total simulated timespan. The physical acceleration is
 *  returned in `aphys'. The latter is used in conjunction with the
 *  PSEUDOSYMMETRIC integration option, which also makes of the second
 *  function of get_timestep. When it is called with a finite timestep for
 *  flag, it returns the physical acceleration that would lead to this
 *  timestep, assuming timestep criterion 0.
 */
int gadgetmp2::get_timestep(int p,         /*!< particle index */
                            my_float *aphys, /*!< acceleration (physical units) */
                            int flag       /*!< either 0 for normal operation, or finite timestep to get corresponding
                            aphys */ )
{
    my_float ax, ay, az, ac, csnd;
    my_float dt = 0, dt_courant = 0, dt_accel;
    int ti_step;

    #ifdef CONDUCTION
    my_float dt_cond;
    #endif

    if(flag == 0)
    {
        ax = fac1 * P[p].GravAccel[0] * 0;
        ay = fac1 * P[p].GravAccel[1] * 0;
        az = fac1 * P[p].GravAccel[2] * 0;

        if(P[p].Type == 0)
        {
            #ifndef TIMESTEP_UPDATE
            ax += fac2 * SphP[p].HydroAccel[0];
            ay += fac2 * SphP[p].HydroAccel[1];
            az += fac2 * SphP[p].HydroAccel[2];
            #else
            if(SphP[p].FeedbackFlag)
            {
                ax += fac2 * SphP[p].FeedAccel[0];
                ay += fac2 * SphP[p].FeedAccel[1];
                az += fac2 * SphP[p].FeedAccel[2];
            }
            else
            {
                ax += fac2 * SphP[p].HydroAccel[0];
                ay += fac2 * SphP[p].HydroAccel[1];
                az += fac2 * SphP[p].HydroAccel[2];
            }
            #endif
        }

        ac = sqrt(ax * ax + ay * ay + az * az);   /* this is now the physical acceleration */
        *aphys = ac;
    }
    else
        ac = *aphys;

    if(ac == 0)
        ac = 1.0e-30;

    switch (All.TypeOfTimestepCriterion)
    {
        case 0:
            if(flag > 0)
            {
                dt = flag * All.Timebase_interval;
                dt /= hubble_a;       /* convert dloga to physical timestep  */
                ac = 2 * All.ErrTolIntAccuracy * atime * All.SofteningTable[P[p].Type] / (dt * dt);
                *aphys = ac;
                return flag;
            }

            if(P[p].Type == 0)
                dt = dt_accel = sqrt(2 * All.ErrTolIntAccuracy * atime * dmin(SphP[p].Hsml, All.SofteningTable[P[p].Type]) / ac);
            else
                dt = dt_accel = sqrt(2 * All.ErrTolIntAccuracy * atime * All.SofteningTable[P[p].Type] / ac);

            #ifdef ADAPTIVE_GRAVSOFT_FORGAS
            if(P[p].Type == 0)
                dt = dt_accel = sqrt(2 * All.ErrTolIntAccuracy * atime * SphP[p].Hsml / 2.8 / ac);
            #endif
            break;
        default:
            endrun(888);
            break;
    }

    if(P[p].Type == 0)
    {
        csnd = sqrt(GAMMA * SphP[p].Pressure / SphP[p].Density);

        if(All.ComovingIntegrationOn)
            dt_courant = 2 * All.CourantFac * All.Time * SphP[p].Hsml / (fac3 * SphP[p].MaxSignalVel);
        else
            dt_courant = 2 * All.CourantFac * SphP[p].Hsml / SphP[p].MaxSignalVel;

        if(dt_courant < dt)
            dt = dt_courant;
    }

    dt *= hubble_a;

    if(dt >= All.MaxSizeTimestep)
        dt = All.MaxSizeTimestep;

    if(dt >= dt_displacement)
        dt = dt_displacement;


    if(dt < All.MinSizeTimestep)
    {
        #ifndef NOSTOP_WHEN_BELOW_MINTIMESTEP
        printf("warning: Timestep wants to be below the limit `MinSizeTimestep'\n");

        if(P[p].Type == 0)
        {
            printf
            ("Part-ID=%d  dt=%g dtc=%g ac=%g xyz=(%g|%g|%g)  hsml=%g  maxsignalvel=%g dt0=%g eps=%g\n",
             (int) P[p].ID, dt.toDouble(), (dt_courant * hubble_a).toDouble(), ac.toDouble(), P[p].Pos[0].toDouble(), P[p].Pos[1].toDouble(), P[p].Pos[2].toDouble(),
             SphP[p].Hsml.toDouble(), SphP[p].MaxSignalVel.toDouble(),
             (sqrt(2 * All.ErrTolIntAccuracy * atime * All.SofteningTable[P[p].Type] / ac) * hubble_a).toDouble(),
             All.SofteningTable[P[p].Type].toDouble());
        }
        else
        {
            printf("Part-ID=%d  dt=%g ac=%g xyz=(%g|%g|%g)\n", (int) P[p].ID, dt.toDouble(), ac.toDouble(), P[p].Pos[0].toDouble(), P[p].Pos[1].toDouble(),
                   P[p].Pos[2].toDouble());
        }
        fflush(stdout);
        endrun(888);
        #endif
        dt = All.MinSizeTimestep;
    }

    ti_step = (dt / All.Timebase_interval).toDouble();

    if(!(ti_step > 0 && ti_step < TIMEBASE) && !ZeroTimestepEncountered)
    {
        printf("\nError: A timestep of size zero was assigned on the integer timeline!\n"
        "We better stop.\n"
        "Task=%d Part-ID=%d dt=%g tibase=%g ti_step=%d ac=%g xyz=(%g|%g|%g) tree=(%g|%g|%g)\n\n",
               ThisTask, (int) P[p].ID, dt.toDouble(), All.Timebase_interval, ti_step, ac.toDouble(),
               P[p].Pos[0].toDouble(), P[p].Pos[1].toDouble(), P[p].Pos[2].toDouble(), P[p].GravAccel[0].toDouble(), P[p].GravAccel[1].toDouble(), P[p].GravAccel[2].toDouble());

        if(P[p].Type == 0)
            printf("hydro-frc=(%g|%g|%g)\n", SphP[p].HydroAccel[0].toDouble(), SphP[p].HydroAccel[1].toDouble(), SphP[p].HydroAccel[2].toDouble());

        #ifdef TIMESTEP_UPDATE
        if(P[p].Type == 0)
            printf("feedback-flag=%d feedback-frc=(%g|%g|%g)\n", SphP[p].FeedbackFlag,SphP[p].FeedAccel[0].toDouble(), SphP[p].FeedAccel[1].toDouble(), SphP[p].FeedAccel[2].toDouble());
        #endif

        fflush(stdout);
        //endrun(818);
        ZeroTimestepEncountered = 1; // Do not terminate the run, but let AMUSE handle it
        ti_step = 1;
    }
    #ifdef TIMESTEP_UPDATE
    if(P[p].Type == 0)
        if(SphP[p].FeedbackFlag<0 || SphP[p].FeedbackFlag>2)
        {
            printf("invalid feedback-flag, Task=%d p=%d Part-ID=%d\n", ThisTask, p, (int) P[p].ID);
        }
        #endif
        return ti_step;
}



/*! This function computes an upper limit ('dt_displacement') to the global
 *  timestep of the system based on the rms velocities of particles. For
 *  cosmological simulations, the criterion used is that the rms displacement
 *  should be at most a fraction MaxRMSDisplacementFac of the mean particle
 *  separation. Note that the latter is estimated using the assigned particle
 *  masses, separately for each particle type. If comoving integration is not
 *  used, the function imposes no constraint on the timestep.
 */
void gadgetmp2::find_dt_displacement_constraint(my_float hfac /*!<  should be  a^2*H(a)  */ )
{
    int i, j, type, *temp;
    int count[6];
    long long count_sum[6];
    my_float v[6], v_sum[6], mim[6], min_mass[6];
    my_float dt, dmean, asmth = 0;

    dt_displacement = All.MaxSizeTimestep;

    if(All.ComovingIntegrationOn)
    {
        for(type = 0; type < 6; type++)
        {
            count[type] = 0;
            v[type] = 0;
            mim[type] = 1.0e30;
        }

        for(i = 0; i < NumPart; i++)
        {
            v[P[i].Type] += P[i].Vel[0] * P[i].Vel[0] + P[i].Vel[1] * P[i].Vel[1] + P[i].Vel[2] * P[i].Vel[2];
            if(mim[P[i].Type] > P[i].Mass)
                mim[P[i].Type] = P[i].Mass;
            count[P[i].Type]++;
        }

        #ifndef NOMPI
        //MPI_Allreduce(v, v_sum, 6, MPI_DOUBLE, MPI_SUM, GADGET_WORLD);
        //MPI_Allreduce(mim, min_mass, 6, MPI_DOUBLE, MPI_MIN, GADGET_WORLD);
        for(i = 0; i < 6; i++){
            v_sum[i] = mpi_all_sum(v[i]);
            min_mass[i] = mpi_all_min(mim[i]);
        }
        #else
        for(i = 0; i < 6; i++){
            v_sum[i] = v[i];
            min_mass[i] = mim[i];
        }
        #endif
        temp = new int[NTask * 6];
        #ifndef NOMPI
        MPI_Allgather(count, 6, MPI_INT, temp, 6, MPI_INT, GADGET_WORLD);
        #else
        for(i = 0; i < 6; i++){
            temp[i] = count[i];
        }
        #endif
        for(i = 0; i < 6; i++)
        {
            count_sum[i] = 0;
            for(j = 0; j < NTask; j++)
                count_sum[i] += temp[j * 6 + i];
        }
        free(temp);

        for(type = 0; type < 6; type++)
        {
            if(count_sum[type] > 0)
            {
                if(type == 0)
                    dmean =
                    pow(min_mass[type] / (All.OmegaBaryon * 3 * All.Hubble * All.Hubble / (8 * M_PI * All.G)),
                        1.0 / 3);
                    else
                        dmean =
                        pow(min_mass[type] /
                        ((All.Omega0 - All.OmegaBaryon) * 3 * All.Hubble * All.Hubble / (8 * M_PI * All.G)),
                            1.0 / 3);

                        dt = All.MaxRMSDisplacementFac * hfac * dmean / sqrt(v_sum[type] / count_sum[type]);

                    if(ThisTask == 0)
                        printf("type=%d  dmean=%g asmth=%g minmass=%g a=%g  sqrt(<p^2>)=%g  dlogmax=%g\n",
                               type, dmean.toDouble(), asmth.toDouble(), min_mass[type].toDouble(), All.Time, sqrt(v_sum[type] / count_sum[type]).toDouble(), dt.toDouble());

                        if(dt < dt_displacement)
                            dt_displacement = dt;
            }
        }

        if(ThisTask == 0)
            printf("displacement time constraint: %g  (%g)\n", dt_displacement.toDouble(), All.MaxSizeTimestep);
    }
}

//#define TIMESTEP_LIMITER
#ifdef TIMESTEP_LIMITER
int gadgetmp2::time_limiter_evaluate(int target, int mode)
{
    int j, k, n, startnode, numngb_inbox;
    my_float h, h2;
    my_float r2, dx, dy, dz;

    my_float pos[3];

    int  sizestep_i, begstep_i;
    int  sizestep_j, begstep_j, endstep_j;
    int  sizestep_j_old, begstep_j_old;

    int  tstart, tend;

    int  CptShrink = 0;

    my_float dt_entr;
    my_float dt_gravkick;
    my_float dt_hydrokick;

    my_float dv[3];

    if(mode == 0)
    {
        //pos        =    P[target].Pos;
        pos[0]        =    P[target].Pos[0];
        pos[1]        =    P[target].Pos[1];
        pos[2]        =    P[target].Pos[2];
        h          = SphP[target].Hsml;
        sizestep_i =    P[target].Ti_sizestep;
        begstep_i  =    P[target].Ti_endstep;
    }
    else
    {
        //pos        =  TimeDataGet[target].Pos;
        pos[0]       =  TimeDataGet->read_re_init_Pos0(target);
        pos[1]       =  TimeDataGet->read_re_init_Pos1(target);
        pos[2]       =  TimeDataGet->read_re_init_Pos2(target);
        //h          =  TimeDataGet[target].Hsml;
        h            =  TimeDataGet->read_re_init_Hsml(target);
        //sizestep_i =  TimeDataGet[target].Size;
        sizestep_i   =  TimeDataGet->read_Size(target);
        //begstep_i  =  TimeDataGet[target].Begin;
        begstep_i    =  TimeDataGet->read_Begin(target);
    }

    h2   = h * h;

    startnode = All.MaxPart;

    do
    {
        numngb_inbox = ngb_treefind_variable(&pos[0], h, &startnode);

        for(n = 0; n < numngb_inbox; n++)
        {
            j = Ngblist[n];

            dx = P[j].Pos[0] - pos[0];
            dy = P[j].Pos[1] - pos[1];
            dz = P[j].Pos[2] - pos[2];

            #ifdef PERIODIC                 /*  now find the closest image in the given box size  */
            if(dx > boxHalf_X)
                dx -= boxSize_X;
            if(dx < -boxHalf_X)
                dx += boxSize_X;
            if(dy > boxHalf_Y)
                dy -= boxSize_Y;
            if(dy < -boxHalf_Y)
                dy += boxSize_Y;
            if(dz > boxHalf_Z)
                dz -= boxSize_Z;
            if(dz < -boxHalf_Z)
                dz += boxSize_Z;
            #endif
            r2 = dx * dx + dy * dy + dz * dz;
            if(r2 < h2)
            {
                /* ################################ */
                sizestep_j_old = P[j].Ti_sizestep;
                sizestep_j     = P[j].Ti_sizestep;
                begstep_j_old  = P[j].Ti_begstep;
                begstep_j      = P[j].Ti_begstep;
                endstep_j      = P[j].Ti_endstep;
                /* ################################ */

                if( (endstep_j - All.Ti_Current > TIMESTEP_LIMITER * sizestep_i)
                    || ( (sizestep_j > TIMESTEP_LIMITER * sizestep_i) && endstep_j == All.Ti_Current ) )
                {
                    if( endstep_j == All.Ti_Current ) /* the particle is synchronized */
                    {
                        /* just change the particle time-step. note: this may not work correctly since the
                         *                         particle is active at this time, and may not be able to communicate the new
                         *                         time-step to neighbours. in the end that is enforced in Saitoh's scheme, so we
                         *                         may do better here */

                        CptShrink++;

                        P[j].Ti_sizestep = TIMESTEP_LIMITER * sizestep_i;
                    }
                    else
                    {
                        CptShrink++;

                        #ifdef SYNCHRONIZATION
                        for(k = 0; begstep_j + k * TIMESTEP_LIMITER * sizestep_i <= All.Ti_Current; k++);

                        endstep_j  = begstep_j + k * TIMESTEP_LIMITER * sizestep_i;
                        sizestep_j = TIMESTEP_LIMITER * sizestep_i;
                        begstep_j  = endstep_j - sizestep_j;
                        #else
                        endstep_j  = begstep_i + TIMESTEP_LIMITER * sizestep_i;
                        sizestep_j = endstep_j - begstep_j;
                        #endif

                        /* ###################################################### */
                        /* it's needed to synchronize mid-step quantities to the
                         *			 new fictious time-step. */
                        tstart = begstep_j_old + sizestep_j_old / 2;
                        tend = (begstep_j + endstep_j) / 2;
                        /* ###################################################### */

                        if(All.ComovingIntegrationOn)
                        {
                            dt_entr = (tend - tstart) * All.Timebase_interval;
                            dt_gravkick = get_gravkick_factor(tstart, tend);
                            dt_hydrokick = get_hydrokick_factor(tstart, tend);
                        }
                        else
                            dt_entr = dt_gravkick = dt_hydrokick = (tend - tstart) * All.Timebase_interval;

                        for(k = 0; k < 3; k++)
                        {
                            dv[k] = P[j].GravAccel[k] * dt_gravkick;
                            P[j].Vel[k] += dv[k];
                        }

                        if(P[j].Type == 0)        /* SPH stuff */
                        {
                            /* note that VelPred is already at the right time */
                            for(k = 0; k < 3; k++)
                                P[j].Vel[k] += SphP[j].HydroAccel[k] * dt_hydrokick;

                            SphP[j].Entropy += SphP[j].DtEntropy * dt_entr;
                        }
                        /* ######################### */
                        P[j].Ti_begstep  = begstep_j;
                        P[j].Ti_endstep  = endstep_j;
                        P[j].Ti_sizestep = sizestep_j;
                        /* ######################### */
                    }
                }
            }
        }
    }
    while(startnode >= 0);

    return CptShrink;
}

/*! This is a comparison kernel for a sort routine, which is used to group
 *  particles that are going to be exported to the same CPU.
 */
int gadgetmp2::time_compare_key(const void *a, const void *b)
{
    if(((timedata_in *) a)->read_Task() < (((timedata_in *) b)->read_Task()))
        return -1;
    if(((timedata_in *) a)->read_Task() > (((timedata_in *) b)->read_Task()))
        return +1;
    return 0;
}






void gadgetmp2::make_it_active(int target)
{
    int k;
    int tstart, tend;

    my_float dt_entr;
    my_float dt_gravkick;
    my_float dt_hydrokick;

    my_float dv[3];


    tstart = (P[target].Ti_begstep + All.Ti_Current) / 2;       /* midpoint of shrinked step */
    tend = (P[target].Ti_begstep + P[target].Ti_endstep) / 2;   /* midpoint of old step */

    if(All.ComovingIntegrationOn)
    {
        dt_entr = (tend - tstart) * All.Timebase_interval;
        dt_gravkick = get_gravkick_factor(tstart, tend);
        dt_hydrokick = get_hydrokick_factor(tstart, tend);
    }
    else
        dt_entr = dt_gravkick = dt_hydrokick = (tend - tstart) * All.Timebase_interval;

    for(k = 0; k < 3; k++)
    {
        dv[k] = P[target].GravAccel[k] * dt_gravkick +
        SphP[target].HydroAccel[k] * dt_hydrokick;
        P[target].Vel[k] -= dv[k];
    }

    SphP[target].Entropy -= SphP[target].DtEntropy * dt_entr;

    P[target].Ti_endstep = All.Ti_Current;
    NumForceUpdate++;
}

#endif
