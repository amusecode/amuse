/*!
 * \copyright   This file is part of the public version of the AREPO code.
 * \copyright   Copyright (C) 2009-2019, Max-Planck Institute for Astrophysics
 * \copyright   Developed by Volker Springel (vspringel@MPA-Garching.MPG.DE) and
 *              contributing authors.
 * \copyright   Arepo is free software: you can redistribute it and/or modify
 *              it under the terms of the GNU General Public License as published by
 *              the Free Software Foundation, either version 3 of the License, or
 *              (at your option) any later version.
 *
 *              Arepo is distributed in the hope that it will be useful,
 *              but WITHOUT ANY WARRANTY; without even the implied warranty of
 *              MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *              GNU General Public License for more details.
 *
 *              A copy of the GNU General Public License is available under
 *              LICENSE as part of this program.  See also
 *              <https://www.gnu.org/licenses/>.
 *
 * \file        src/subfind/subfind_so.c
 * \date        05/2018
 * \brief       Spherical overdensity algorithm for subfind.
 * \details     contains functions:
 *                static void particle2in(data_in * in, int i, int firstnode)
 *                static void out2particle(data_out * out, int i, int mode)
 *                static void kernel_local(void)
 *                static void kernel_imported(void)
 *                double subfind_overdensity(void)
 *                static int subfind_overdensity_evaluate(int target, int mode,
 *                  int threadid)
 *
 *
 * \par Major modifications and contributions:
 *
 * - DD.MM.YYYY Description
 * - 14.05.2018 Prepared file for public release -- Rainer Weinberger
 */

#include <gsl/gsl_math.h>
#include <math.h>
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>
#include <sys/types.h>

#include "../main/allvars.h"
#include "../main/proto.h"

#ifdef SUBFIND

#include "../fof/fof.h"
#include "subfind.h"

static double *R200, *M200;

static char *Todo;
static MyFloat *Left, *Right;
static int mainstep;

static int subfind_overdensity_evaluate(int target, int mode, int threadid);

#ifdef SUBFIND_EXTENDED_PROPERTIES
/*! \brief Structure for angular momentum properties.
 */
static struct Angular_Momentum
{
  double Pmom[3];
  double MassType[NTYPES];
  double Jtot[3];
  double Jdm[3];
  double Jgas[3];
  double Jstars[3];
  int LenType[NTYPES];
  double CMFrac;
  double CMFracType[NTYPES];
  double Ekin;
  double Epot;
  double Ethr;
  double N200;
} * AngMom;
#endif /* #ifdef SUBFIND_EXTENDED_PROPERTIES */

/*! \brief Local data structure for collecting particle/cell data that is sent
 *         to other processors if needed. Type called data_in and static
 *         pointers DataIn and DataGet needed by generic_comm_helpers2.
 */
typedef struct
{
  MyDouble Pos[3];
  double R200;

#ifdef SUBFIND_EXTENDED_PROPERTIES
  double M200;
  int GrNr;
  int TaskOfGr;
  int LocGrIndex;
  struct Angular_Momentum AngMomIn;
#endif /* #ifdef SUBFIND_EXTENDED_PROPERTIES */

  int Firstnode;
} data_in;

static data_in *DataIn, *DataGet;

/*! \brief Routine that fills the relevant group data into the input
 *         structure defined above. Needed by generic_comm_helpers2.
 *
 *  \param[out] in Data structure to fill.
 *  \param[in] i Index of particle in group arrays.
 *  \param[in] firstnode First note of communication.
 *
 *  \return void
 */
static void particle2in(data_in *in, int i, int firstnode)
{
  in->Pos[0] = Group[i].Pos[0];
  in->Pos[1] = Group[i].Pos[1];
  in->Pos[2] = Group[i].Pos[2];
  in->R200   = R200[i];

#ifdef SUBFIND_EXTENDED_PROPERTIES
  in->GrNr       = Group[i].GrNr;
  in->TaskOfGr   = ThisTask;
  in->LocGrIndex = i;
  in->M200       = M200[i];
  in->AngMomIn   = AngMom[i];
#endif /* #ifdef SUBFIND_EXTENDED_PROPERTIES */

  in->Firstnode = firstnode;
}

/*! \brief Local data structure that holds results acquired on remote
 *         processors. Type called data_out and static pointers DataResult and
 *         DataOut needed by generic_comm_helpers2.
 */
typedef struct
{
  double Mass;

#ifdef SUBFIND_EXTENDED_PROPERTIES
  struct Angular_Momentum AngMomOut;
#endif /* #ifdef SUBFIND_EXTENDED_PROPERTIES */

} data_out;

static data_out *DataResult, *DataOut;

/*! \brief Routine to store or combine result data. Needed by
 *         generic_comm_helpers2.
 *
 *  \param[in] out Data to be moved to appropriate variables in global
 *  particle and group data arrays (AngMom,...)
 *  \param[in] i Index of particle in group arrays
 *  \param[in] mode Mode of function: local particles or information that was
 *  communicated from other tasks and has to be added locally?
 *
 *  \return void
 */
static void out2particle(data_out *out, int i, int mode)
{
  if(mode == MODE_LOCAL_PARTICLES) /* initial store */
    {
      if(mainstep == 0)
        M200[i] = out->Mass;
#ifdef SUBFIND_EXTENDED_PROPERTIES
      if(mainstep == 0)
        {
          for(int k = 0; k < 3; k++)
            AngMom[i].Pmom[k] = out->AngMomOut.Pmom[k];
          for(int k = 0; k < NTYPES; k++)
            {
              AngMom[i].MassType[k] = out->AngMomOut.MassType[k];
              AngMom[i].LenType[k]  = out->AngMomOut.LenType[k];
            }
          AngMom[i].N200 = out->AngMomOut.N200;
        }
      else if(mainstep == 1)
        {
          for(int k = 0; k < 3; k++)
            {
              AngMom[i].Jtot[k]   = out->AngMomOut.Jtot[k];
              AngMom[i].Jdm[k]    = out->AngMomOut.Jdm[k];
              AngMom[i].Jgas[k]   = out->AngMomOut.Jgas[k];
              AngMom[i].Jstars[k] = out->AngMomOut.Jstars[k];
            }
          AngMom[i].Ekin = out->AngMomOut.Ekin;
          AngMom[i].Ethr = out->AngMomOut.Ethr;
        }
      else if(mainstep == 2)
        {
          AngMom[i].CMFrac = out->AngMomOut.CMFrac;
          for(int k = 0; k < NTYPES; k++)
            AngMom[i].CMFracType[k] = out->AngMomOut.CMFracType[k];
        }
#endif /* #ifdef SUBFIND_EXTENDED_PROPERTIES */
    }
  else /* combine */
    {
      if(mainstep == 0)
        M200[i] += out->Mass;
#ifdef SUBFIND_EXTENDED_PROPERTIES
      if(mainstep == 0)
        {
          for(int k = 0; k < 3; k++)
            AngMom[i].Pmom[k] += out->AngMomOut.Pmom[k];
          for(int k = 0; k < NTYPES; k++)
            {
              AngMom[i].MassType[k] += out->AngMomOut.MassType[k];
              AngMom[i].LenType[k] += out->AngMomOut.LenType[k];
            }
          AngMom[i].N200 += out->AngMomOut.N200;
        }
      else if(mainstep == 1)
        {
          for(int k = 0; k < 3; k++)
            {
              AngMom[i].Jtot[k] += out->AngMomOut.Jtot[k];
              AngMom[i].Jdm[k] += out->AngMomOut.Jdm[k];
              AngMom[i].Jgas[k] += out->AngMomOut.Jgas[k];
              AngMom[i].Jstars[k] += out->AngMomOut.Jstars[k];
            }
          AngMom[i].Ekin += out->AngMomOut.Ekin;
          AngMom[i].Ethr += out->AngMomOut.Ethr;
        }
      else if(mainstep == 2)
        {
          AngMom[i].CMFrac += out->AngMomOut.CMFrac;
          for(int k = 0; k < NTYPES; k++)
            AngMom[i].CMFracType[k] += out->AngMomOut.CMFracType[k];
        }
#endif /* #ifdef SUBFIND_EXTENDED_PROPERTIES */
    }
}

#include "../utils/generic_comm_helpers2.h"

/*! \brief Routine that defines what to do with local particles.
 *
 *  Calls the *_evaluate function in MODE_LOCAL_PARTICLES.
 *
 *  \return void
 */
static void kernel_local(void)
{
  int i;

  {
    int threadid = get_thread_num();

    for(int j = 0; j < NTask; j++)
      Thread[threadid].Exportflag[j] = -1;

    while(1)
      {
        if(Thread[threadid].ExportSpace < MinSpace)
          break;

        i = NextParticle++;

        if(i >= Ngroups)
          break;

        if(Todo[i])
          {
            R200[i] = 0.5 * (Left[i] + Right[i]);
            subfind_overdensity_evaluate(i, MODE_LOCAL_PARTICLES, threadid);
          }
      }
  }
}

/*! \brief Routine that defines what to do with imported particles.
 *
 *  Calls the *_evaluate function in MODE_IMPORTED_PARTICLES.
 *
 *  \return void
 */
static void kernel_imported(void)
{
  /* now do the particles that were sent to us */
  int i, cnt = 0;

  {
    int threadid = get_thread_num();

    while(1)
      {
        i = cnt++;

        if(i >= Nimport)
          break;

        subfind_overdensity_evaluate(i, MODE_IMPORTED_PARTICLES, threadid);
      }
  }
}

/*! \brief Main routine executing the spherical overdensity algorithm.
 *
 *  \return Time needed for calculation.
 */
double subfind_overdensity(void)
{
  long long ntot;
  int i, npleft, rep, iter;
  double t0, t1, overdensity, Deltas[4], rhoback, z, omegaz, x, DeltaMean200, DeltaCrit200, DeltaCrit500, DeltaTopHat;
  double tstart = second();

  Left  = (MyFloat *)mymalloc("Left", sizeof(MyFloat) * Ngroups);
  Right = (MyFloat *)mymalloc("Right", sizeof(MyFloat) * Ngroups);
  R200  = (double *)mymalloc("R200", sizeof(double) * Ngroups);
  M200  = (double *)mymalloc("M200", sizeof(double) * Ngroups);
#ifdef SUBFIND_EXTENDED_PROPERTIES
  AngMom = (struct Angular_Momentum *)mymalloc("AngMom", sizeof(struct Angular_Momentum) * Ngroups);
  Paux   = (struct paux_data *)mymalloc("Paux", sizeof(struct paux_data) * NumPart);
#endif /* #ifdef SUBFIND_EXTENDED_PROPERTIES */

  Todo = mymalloc("Todo", sizeof(char) * Ngroups);

  if(All.ComovingIntegrationOn)
    z = 1 / All.Time - 1;
  else
    z = 0;

  rhoback = 3 * All.Omega0 * All.Hubble * All.Hubble / (8 * M_PI * All.G);

  omegaz =
      All.Omega0 * pow(1 + z, 3) / (All.Omega0 * pow(1 + z, 3) + (1 - All.Omega0 - All.OmegaLambda) * pow(1 + z, 2) + All.OmegaLambda);

  DeltaMean200 = 200.0;
  DeltaCrit200 = 200.0 / omegaz;
  DeltaCrit500 = 500.0 / omegaz;

  x           = omegaz - 1;
  DeltaTopHat = 18 * M_PI * M_PI + 82 * x - 39 * x * x;
  DeltaTopHat /= omegaz;

  Deltas[0] = DeltaMean200; /* standard fixed overdensity with respect to background */
  Deltas[1] = DeltaTopHat;  /* tophat overdensity with respect to background */
  Deltas[2] = DeltaCrit200; /* overdensity of 200 relative to critical, expressed relative to background density */
  Deltas[3] = DeltaCrit500; /* overdensity of 500 relative to critical, expressed relative to background density */

  generic_set_MaxNexport();

  for(rep = 0; rep < 4; rep++) /* repeat for all four overdensity values */
    {
#ifdef SUBFIND_EXTENDED_PROPERTIES
      int mainstepmax = 3;
#else  /* #ifdef SUBFIND_EXTENDED_PROPERTIES */
      int mainstepmax = 1;
#endif /* #ifdef SUBFIND_EXTENDED_PROPERTIES #else */
      for(mainstep = 0; mainstep < mainstepmax; mainstep++)
        {
          for(i = 0; i < Ngroups; i++)
            {
              if(Group[i].Nsubs > 0)
                {
                  if(mainstep == 0)
                    {
                      double rguess = pow(All.G * Group[i].Mass / (100 * All.Hubble * All.Hubble), 1.0 / 3);

                      Right[i] = 3 * rguess;
                      Left[i]  = 0;
                    }
                  Todo[i] = 1;
                }
              else
                {
                  Todo[i] = 0;
                }
            }

          iter = 0;

#ifdef SUBFIND_EXTENDED_PROPERTIES
          if(mainstep == 1)
            NumPaux = 0;
#endif /* #ifdef SUBFIND_EXTENDED_PROPERTIES */

          /* we will repeat the whole thing for those groups where we didn't converge to a SO radius yet */
          do
            {
              t0 = second();

              generic_comm_pattern(Ngroups, kernel_local, kernel_imported);

              if(mainstep == 0)
                {
                  /* do final operations on results */
                  for(i = 0, npleft = 0; i < Ngroups; i++)
                    {
                      if(Todo[i])
                        {
                          overdensity = M200[i] / (4.0 * M_PI / 3.0 * R200[i] * R200[i] * R200[i]) / rhoback;

                          if((Right[i] - Left[i]) > 1.0e-4 * Left[i])
                            {
                              /* need to redo this group */
                              npleft++;

                              if(overdensity > Deltas[rep])
                                Left[i] = R200[i];
                              else
                                Right[i] = R200[i];

                              if(iter >= MAXITER - 10)
                                {
                                  printf("gr=%d task=%d  R200=%g Left=%g Right=%g Menclosed=%g Right-Left=%g\n   pos=(%g|%g|%g)\n", i,
                                         ThisTask, R200[i], Left[i], Right[i], M200[i], Right[i] - Left[i], Group[i].Pos[0],
                                         Group[i].Pos[1], Group[i].Pos[2]);
                                  myflush(stdout);
                                }
                            }
                          else
                            Todo[i] = 0;
                        }
                    }
                }
              else
                for(i = 0, npleft = 0; i < Ngroups; i++)
                  Todo[i] = 0;

              sumup_large_ints(1, &npleft, &ntot);

              t1 = second();

              if(ntot > 0)
                {
                  iter++;

                  if(iter > 0)
                    mpi_printf("SUBFIND: SO iteration %2d: need to repeat for %12lld halo centers. (took %g sec)\n", iter, ntot,
                               timediff(t0, t1));

                  if(iter > MAXITER)
                    terminate("failed to converge in SO iteration");
                }
            }
          while(ntot > 0);
        } /* end of mainstep loop */

#ifdef SUBFIND_EXTENDED_PROPERTIES
      double *egypot = mymalloc("egypot", Ngroups * sizeof(double));

      subfind_so_potegy(egypot);

      for(i = 0; i < Ngroups; i++)
        {
          double rate;

          /* work out sampling rate */
          if(AngMom[i].N200 < SUBFIND_SO_POT_CALCULATION_PARTICLE_NUMBER)
            rate = 1.0;
          else
            rate = (SUBFIND_SO_POT_CALCULATION_PARTICLE_NUMBER / AngMom[i].N200);

          AngMom[i].Epot = egypot[i] / (rate * rate);
        }

      myfree(egypot);
#endif /* #ifdef SUBFIND_EXTENDED_PROPERTIES */

      for(i = 0; i < Ngroups; i++)
        {
          if(Group[i].Nsubs > 0)
            {
              overdensity = M200[i] / (4.0 * M_PI / 3.0 * R200[i] * R200[i] * R200[i]) / rhoback;

              if((overdensity - Deltas[rep]) > 0.1 * Deltas[rep])
                {
                  R200[i] = M200[i] = 0;
#ifdef SUBFIND_EXTENDED_PROPERTIES
                  memset(&AngMom[i], 0, sizeof(struct Angular_Momentum));
#endif /* #ifdef SUBFIND_EXTENDED_PROPERTIES */
                }
              else if(M200[i] < 5 * Group[i].Mass / Group[i].Len)
                {
                  R200[i] = M200[i] = 0;
#ifdef SUBFIND_EXTENDED_PROPERTIES
                  memset(&AngMom[i], 0, sizeof(struct Angular_Momentum));
#endif /* #ifdef SUBFIND_EXTENDED_PROPERTIES */
                }
            }
          else
            {
              R200[i] = M200[i] = 0;
#ifdef SUBFIND_EXTENDED_PROPERTIES
              memset(&AngMom[i], 0, sizeof(struct Angular_Momentum));
#endif /* #ifdef SUBFIND_EXTENDED_PROPERTIES */
            }

          switch(rep)
            {
              case 0:
                Group[i].M_Mean200 = M200[i];
                Group[i].R_Mean200 = R200[i];
#ifdef SUBFIND_EXTENDED_PROPERTIES
                Group[i].Ekin_Mean200   = AngMom[i].Ekin;
                Group[i].Ethr_Mean200   = AngMom[i].Ethr;
                Group[i].Epot_Mean200   = AngMom[i].Epot;
                Group[i].CMFrac_Mean200 = AngMom[i].CMFrac;
                for(int k = 0; k < NTYPES; k++)
                  {
                    Group[i].MassType_Mean200[k]   = AngMom[i].MassType[k];
                    Group[i].LenType_Mean200[k]    = AngMom[i].LenType[k];
                    Group[i].CMFracType_Mean200[k] = AngMom[i].CMFracType[k];
                  }
                for(int k = 0; k < 3; k++)
                  {
                    Group[i].J_Mean200[k]      = AngMom[i].Jtot[k];
                    Group[i].JDM_Mean200[k]    = AngMom[i].Jdm[k];
                    Group[i].JGas_Mean200[k]   = AngMom[i].Jgas[k];
                    Group[i].JStars_Mean200[k] = AngMom[i].Jstars[k];
                  }
#endif /* #ifdef SUBFIND_EXTENDED_PROPERTIES */
                break;
              case 1:
                Group[i].M_TopHat200 = M200[i];
                Group[i].R_TopHat200 = R200[i];
#ifdef SUBFIND_EXTENDED_PROPERTIES
                Group[i].Ekin_TopHat200   = AngMom[i].Ekin;
                Group[i].Ethr_TopHat200   = AngMom[i].Ethr;
                Group[i].Epot_TopHat200   = AngMom[i].Epot;
                Group[i].CMFrac_TopHat200 = AngMom[i].CMFrac;
                for(int k = 0; k < NTYPES; k++)
                  {
                    Group[i].MassType_TopHat200[k]   = AngMom[i].MassType[k];
                    Group[i].LenType_TopHat200[k]    = AngMom[i].LenType[k];
                    Group[i].CMFracType_TopHat200[k] = AngMom[i].CMFracType[k];
                  }
                for(int k = 0; k < 3; k++)
                  {
                    Group[i].J_TopHat200[k]      = AngMom[i].Jtot[k];
                    Group[i].JDM_TopHat200[k]    = AngMom[i].Jdm[k];
                    Group[i].JGas_TopHat200[k]   = AngMom[i].Jgas[k];
                    Group[i].JStars_TopHat200[k] = AngMom[i].Jstars[k];
                  }
#endif /* #ifdef SUBFIND_EXTENDED_PROPERTIES */
                break;
              case 2:
                Group[i].M_Crit200 = M200[i];
                Group[i].R_Crit200 = R200[i];
#ifdef SUBFIND_EXTENDED_PROPERTIES
                Group[i].Ekin_Crit200   = AngMom[i].Ekin;
                Group[i].Ethr_Crit200   = AngMom[i].Ethr;
                Group[i].Epot_Crit200   = AngMom[i].Epot;
                Group[i].CMFrac_Crit200 = AngMom[i].CMFrac;
                for(int k = 0; k < NTYPES; k++)
                  {
                    Group[i].MassType_Crit200[k]   = AngMom[i].MassType[k];
                    Group[i].LenType_Crit200[k]    = AngMom[i].LenType[k];
                    Group[i].CMFracType_Crit200[k] = AngMom[i].CMFracType[k];
                  }
                for(int k = 0; k < 3; k++)
                  {
                    Group[i].J_Crit200[k]      = AngMom[i].Jtot[k];
                    Group[i].JDM_Crit200[k]    = AngMom[i].Jdm[k];
                    Group[i].JGas_Crit200[k]   = AngMom[i].Jgas[k];
                    Group[i].JStars_Crit200[k] = AngMom[i].Jstars[k];
                  }
#endif /* #ifdef SUBFIND_EXTENDED_PROPERTIES */
                break;
              case 3:
                Group[i].M_Crit500 = M200[i];
                Group[i].R_Crit500 = R200[i];
#ifdef SUBFIND_EXTENDED_PROPERTIES
                Group[i].Ekin_Crit500   = AngMom[i].Ekin;
                Group[i].Ethr_Crit500   = AngMom[i].Ethr;
                Group[i].Epot_Crit500   = AngMom[i].Epot;
                Group[i].CMFrac_Crit500 = AngMom[i].CMFrac;
                for(int k = 0; k < NTYPES; k++)
                  {
                    Group[i].MassType_Crit500[k]   = AngMom[i].MassType[k];
                    Group[i].LenType_Crit500[k]    = AngMom[i].LenType[k];
                    Group[i].CMFracType_Crit500[k] = AngMom[i].CMFracType[k];
                  }
                for(int k = 0; k < 3; k++)
                  {
                    Group[i].J_Crit500[k]      = AngMom[i].Jtot[k];
                    Group[i].JDM_Crit500[k]    = AngMom[i].Jdm[k];
                    Group[i].JGas_Crit500[k]   = AngMom[i].Jgas[k];
                    Group[i].JStars_Crit500[k] = AngMom[i].Jstars[k];
                  }
#endif /* #ifdef SUBFIND_EXTENDED_PROPERTIES */
                break;
            }
        }
    }

  myfree(Todo);
#ifdef SUBFIND_EXTENDED_PROPERTIES
  myfree(Paux);
  myfree(AngMom);
#endif /* #ifdef SUBFIND_EXTENDED_PROPERTIES */
  myfree(M200);
  myfree(R200);
  myfree(Right);
  myfree(Left);

  double tend = second();
  return timediff(tstart, tend);
}

/*! \brief Evaluate function of subfind_overdensity.
 *
 *  \param[in] target Index of group.
 *  \param[in] mode Flag if it operates on local or imported data.
 *  \param[in] threadid ID of thread.
 *
 *  \return 0
 */
static int subfind_overdensity_evaluate(int target, int mode, int threadid)
{
  int k, p, no, numnodes, *firstnode;
  double hsml, mass;
  MyDouble *pos;
  struct NODE *current;
  MyDouble dx, dy, dz, dist, r2;
#define FACT2 0.86602540
  MyDouble xtmp, ytmp, ztmp;

  data_in local, *in;
  data_out out;

  if(mode == MODE_LOCAL_PARTICLES)
    {
      particle2in(&local, target, 0);
      in = &local;

      numnodes  = 1;
      firstnode = NULL;
    }
  else
    {
      in = &DataGet[target];

      generic_get_numnodes(target, &numnodes, &firstnode);
    }

  pos  = in->Pos;
  hsml = in->R200;
  mass = 0;

#ifdef SUBFIND_EXTENDED_PROPERTIES
  double Pmom[3], Mtot = 0, Jtot[3], Jdm[3], Jgas[3], Jstars[3], CMFrac = 0, N200 = 0;
  double ekin = 0, etherm = 0;
  double MassType[NTYPES], CMFracType[NTYPES];
  int LenType[NTYPES];

  for(int i = 0; i < 3; i++)
    {
      Pmom[i]   = 0;
      Jtot[i]   = 0;
      Jdm[i]    = 0;
      Jgas[i]   = 0;
      Jstars[i] = 0;
    }
  for(int i = 0; i < NTYPES; i++)
    {
      MassType[i]   = 0;
      LenType[i]    = 0;
      CMFracType[i] = 0;
    }

  if(mainstep == 1)
    {
      Mtot = in->M200;
      N200 = in->AngMomIn.N200;
      for(int i = 0; i < 3; i++)
        Pmom[i] = in->AngMomIn.Pmom[i];
    }
  else if(mainstep == 2)
    {
      Mtot = in->M200;
      for(int i = 0; i < 3; i++)
        {
          Pmom[i]   = in->AngMomIn.Pmom[i];
          Jtot[i]   = in->AngMomIn.Jtot[i];
          Jdm[i]    = in->AngMomIn.Jdm[i];
          Jgas[i]   = in->AngMomIn.Jgas[i];
          Jstars[i] = in->AngMomIn.Jstars[i];
        }
      for(int i = 0; i < NTYPES; i++)
        MassType[i] = in->AngMomIn.MassType[i];
    }
#endif /* #ifdef SUBFIND_EXTENDED_PROPERTIES */

  for(k = 0; k < numnodes; k++)
    {
      if(mode == MODE_LOCAL_PARTICLES)
        {
          no = Tree_MaxPart; /* root node */
        }
      else
        {
          no = firstnode[k];
          no = Nodes[no].u.d.nextnode; /* open it */
        }

      while(no >= 0)
        {
          if(no < Tree_MaxPart) /* single particle */
            {
              p  = no;
              no = Nextnode[no];

              dist = hsml;
              dx   = FOF_NEAREST_LONG_X(Tree_Pos_list[3 * p + 0] - pos[0]);
              if(dx > dist)
                continue;
              dy = FOF_NEAREST_LONG_Y(Tree_Pos_list[3 * p + 1] - pos[1]);
              if(dy > dist)
                continue;
              dz = FOF_NEAREST_LONG_Z(Tree_Pos_list[3 * p + 2] - pos[2]);
              if(dz > dist)
                continue;
              if(dx * dx + dy * dy + dz * dz > dist * dist)
                continue;

              if(mainstep == 0)
                mass += P[p].Mass;

#ifdef SUBFIND_EXTENDED_PROPERTIES
              if(mainstep == 0)
                {
                  for(int i = 0; i < 3; i++)
                    Pmom[i] += P[p].Mass * P[p].Vel[i] / All.cf_atime;  // units: 10^10 M_sol/h km/s

                  for(int i = 0; i < NTYPES; i++)
                    if(P[p].Type == i)
                      {
                        MassType[i] += P[p].Mass;

                        LenType[i]++;
                      }

                  N200 += 1.0;
                }
              else if(mainstep == 1)
                {
                  double rate;
                  /* work out sampling rate */
                  if(N200 < SUBFIND_SO_POT_CALCULATION_PARTICLE_NUMBER)
                    rate = 1.0;
                  else
                    rate = (SUBFIND_SO_POT_CALCULATION_PARTICLE_NUMBER / N200);

                  if(get_random_number_aux() < rate)
                    {
                      if(NumPaux >= NumPart)
                        terminate("NumPaux >= NumPart");

                      Paux[NumPaux].Pos[0]        = NEAREST_X(P[p].Pos[0] - pos[0]);
                      Paux[NumPaux].Pos[1]        = NEAREST_Y(P[p].Pos[1] - pos[1]);
                      Paux[NumPaux].Pos[2]        = NEAREST_Z(P[p].Pos[2] - pos[2]);
                      Paux[NumPaux].Mass          = P[p].Mass;
                      Paux[NumPaux].TaskOfGr      = in->TaskOfGr;
                      Paux[NumPaux].LocGrIndex    = in->LocGrIndex;
                      Paux[NumPaux].Type          = P[p].Type;
                      Paux[NumPaux].SofteningType = P[p].SofteningType;
                      NumPaux++;
                    }

                  int ptype = P[p].Type;

                  double Pos_pbc[3], Vel_centre[3], Vel_tot[3];
                  Pos_pbc[0] = NEAREST_X(P[p].Pos[0] - pos[0]) * All.cf_atime;
                  Pos_pbc[1] = NEAREST_Y(P[p].Pos[1] - pos[1]) * All.cf_atime;
                  Pos_pbc[2] = NEAREST_Z(P[p].Pos[2] - pos[2]) * All.cf_atime;

                  for(int i = 0; i < 3; i++)
                    Vel_centre[i] = (Pmom[i] / Mtot);  // units: km/s

                  for(int i = 0; i < 3; i++)
                    Vel_tot[i] = P[p].Vel[i] / All.cf_atime - Vel_centre[i] + All.cf_Hrate * Pos_pbc[i];

                  ekin += 0.5 * P[p].Mass * (Vel_tot[0] * Vel_tot[0] + Vel_tot[1] * Vel_tot[1] + Vel_tot[2] * Vel_tot[2]);

                  Jtot[0] += P[p].Mass * (Pos_pbc[1] * Vel_tot[2] - Pos_pbc[2] * Vel_tot[1]);
                  Jtot[1] += P[p].Mass * (Pos_pbc[2] * Vel_tot[0] - Pos_pbc[0] * Vel_tot[2]);
                  Jtot[2] += P[p].Mass * (Pos_pbc[0] * Vel_tot[1] - Pos_pbc[1] * Vel_tot[0]);

                  if(ptype == 1)  // dm illustris
                    {
                      Jdm[0] += P[p].Mass * (Pos_pbc[1] * Vel_tot[2] - Pos_pbc[2] * Vel_tot[1]);
                      Jdm[1] += P[p].Mass * (Pos_pbc[2] * Vel_tot[0] - Pos_pbc[0] * Vel_tot[2]);
                      Jdm[2] += P[p].Mass * (Pos_pbc[0] * Vel_tot[1] - Pos_pbc[1] * Vel_tot[0]);
                    }
                  if(ptype == 0)  // gas
                    {
                      etherm += P[p].Mass * PS[p].Utherm;

                      Jgas[0] += P[p].Mass * (Pos_pbc[1] * Vel_tot[2] - Pos_pbc[2] * Vel_tot[1]);
                      Jgas[1] += P[p].Mass * (Pos_pbc[2] * Vel_tot[0] - Pos_pbc[0] * Vel_tot[2]);
                      Jgas[2] += P[p].Mass * (Pos_pbc[0] * Vel_tot[1] - Pos_pbc[1] * Vel_tot[0]);
                    }
                  if(ptype == 4)  // stars
                    {
                      Jstars[0] += P[p].Mass * (Pos_pbc[1] * Vel_tot[2] - Pos_pbc[2] * Vel_tot[1]);
                      Jstars[1] += P[p].Mass * (Pos_pbc[2] * Vel_tot[0] - Pos_pbc[0] * Vel_tot[2]);
                      Jstars[2] += P[p].Mass * (Pos_pbc[0] * Vel_tot[1] - Pos_pbc[1] * Vel_tot[0]);
                    }
                }
              else if(mainstep == 2)
                {
                  int ptype = P[p].Type;

                  double Pos_pbc[3], Vel_centre[3], Vel_tot[3], jpart[3], Jtot[3];
                  Pos_pbc[0] = NEAREST_X(P[p].Pos[0] - pos[0]) * All.cf_atime;
                  Pos_pbc[1] = NEAREST_Y(P[p].Pos[1] - pos[1]) * All.cf_atime;
                  Pos_pbc[2] = NEAREST_Z(P[p].Pos[2] - pos[2]) * All.cf_atime;

                  for(int i = 0; i < 3; i++)
                    Vel_centre[i] = (Pmom[i] / Mtot);

                  for(int i = 0; i < 3; i++)
                    Vel_tot[i] = P[p].Vel[i] / All.cf_atime - Vel_centre[i] + All.cf_Hrate * Pos_pbc[i];

                  jpart[0] = P[p].Mass * (Pos_pbc[1] * Vel_tot[2] - Pos_pbc[2] * Vel_tot[1]);
                  jpart[1] = P[p].Mass * (Pos_pbc[2] * Vel_tot[0] - Pos_pbc[0] * Vel_tot[2]);
                  jpart[2] = P[p].Mass * (Pos_pbc[0] * Vel_tot[1] - Pos_pbc[1] * Vel_tot[0]);

                  if((Jtot[0] * jpart[0] + Jtot[1] * jpart[1] + Jtot[2] * jpart[2]) < 0.)
                    CMFrac += P[p].Mass / Mtot;

                  if(ptype == 1)  // dm
                    if((Jdm[0] * jpart[0] + Jdm[1] * jpart[1] + Jdm[2] * jpart[2]) < 0.)
                      CMFracType[1] += P[p].Mass / MassType[1];

                  if(ptype == 0)  // gas
                    if((Jgas[0] * jpart[0] + Jgas[1] * jpart[1] + Jgas[2] * jpart[2]) < 0.)
                      CMFracType[0] += P[p].Mass / MassType[0];

                  if(ptype == 4)  // stars
                    if((Jstars[0] * jpart[0] + Jstars[1] * jpart[1] + Jstars[2] * jpart[2]) < 0.)
                      CMFracType[4] += P[p].Mass / MassType[4];
                }
#endif /* #ifdef SUBFIND_EXTENDED_PROPERTIES */
            }
          else if(no < Tree_MaxPart + Tree_MaxNodes) /* internal node */
            {
              if(mode == MODE_IMPORTED_PARTICLES)
                {
                  if(no <
                     Tree_FirstNonTopLevelNode) /* we reached a top-level node again, which means that we are done with the branch */
                    break;
                }

              current = &Nodes[no];

              no = current->u.d.sibling; /* in case the node can be discarded */

              dist = hsml + 0.5 * current->len;
              dx   = FOF_NEAREST_LONG_X(current->center[0] - pos[0]);
              if(dx > dist)
                continue;
              dy = FOF_NEAREST_LONG_Y(current->center[1] - pos[1]);
              if(dy > dist)
                continue;
              dz = FOF_NEAREST_LONG_Z(current->center[2] - pos[2]);
              if(dz > dist)
                continue;
              /* now test against the minimal sphere enclosing everything */
              dist += FACT1 * current->len;
              if((r2 = (dx * dx + dy * dy + dz * dz)) > dist * dist)
                continue;

#ifndef SUBFIND_EXTENDED_PROPERTIES
              if(no >= Tree_FirstNonTopLevelNode) /* only do this for fully local nodes */
                {
                  /* test whether the node is contained within the sphere, which gives  short-cut if we only need the mass */
                  dist = hsml - FACT2 * current->len;
                  if(dist > 0)
                    if(r2 < dist * dist)
                      {
                        mass += current->u.d.mass;
                        continue;
                      }
                }
#endif /* #ifndef SUBFIND_EXTENDED_PROPERTIES */

              no = current->u.d.nextnode; /* ok, we need to open the node */
            }
          else if(no >= Tree_ImportedNodeOffset) /* point from imported nodelist */
            {
              int n = no - Tree_ImportedNodeOffset;
              no    = Nextnode[no - Tree_MaxNodes];

              dist = hsml;
              dx   = FOF_NEAREST_LONG_X(Tree_Points[n].Pos[0] - pos[0]);
              if(dx > dist)
                continue;
              dy = FOF_NEAREST_LONG_Y(Tree_Points[n].Pos[1] - pos[1]);
              if(dy > dist)
                continue;
              dz = FOF_NEAREST_LONG_Z(Tree_Points[n].Pos[2] - pos[2]);
              if(dz > dist)
                continue;
              if(dx * dx + dy * dy + dz * dz > dist * dist)
                continue;

              mass += Tree_Points[n].Mass;
            }
          else /* pseudo particle */
            {
              if(mode == MODE_IMPORTED_PARTICLES)
                terminate("mode == MODE_IMPORTED_PARTICLES");

              if(mode == MODE_LOCAL_PARTICLES)
                tree_treefind_export_node_threads(no, target, threadid);

              no = Nextnode[no - Tree_MaxNodes];
            }
        }
    }

  out.Mass = mass;

#ifdef SUBFIND_EXTENDED_PROPERTIES
  if(mainstep == 0)
    {
      for(int k = 0; k < 3; k++)
        out.AngMomOut.Pmom[k] = Pmom[k];
      for(int k = 0; k < NTYPES; k++)
        {
          out.AngMomOut.MassType[k] = MassType[k];
          out.AngMomOut.LenType[k]  = LenType[k];
        }

      out.AngMomOut.N200 = N200;
    }
  else if(mainstep == 1)
    {
      for(int k = 0; k < 3; k++)
        {
          out.AngMomOut.Jtot[k]   = Jtot[k];
          out.AngMomOut.Jdm[k]    = Jdm[k];
          out.AngMomOut.Jgas[k]   = Jgas[k];
          out.AngMomOut.Jstars[k] = Jstars[k];
        }

      out.AngMomOut.Ekin = ekin;
      out.AngMomOut.Ethr = etherm;
    }
  else if(mainstep == 2)
    {
      out.AngMomOut.CMFrac = CMFrac;
      for(int k = 0; k < NTYPES; k++)
        out.AngMomOut.CMFracType[k] = CMFracType[k];
    }
#endif /* #ifdef SUBFIND_EXTENDED_PROPERTIES */

  /* Now collect the result at the right place */
  if(mode == MODE_LOCAL_PARTICLES)
    out2particle(&out, target, MODE_LOCAL_PARTICLES);
  else
    DataResult[target] = out;

  return 0;
}

#endif /* #ifdef SUBFIND */
