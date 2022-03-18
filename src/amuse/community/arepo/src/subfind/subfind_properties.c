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
 * \file        src/subfind/subfind_properties.c
 * \date        05/2018
 * \brief       Calculation of the subgroup properties.
 * \details     contains functions:
 *                void subfind_determine_sub_halo_properties(struct
 *                  unbind_data *d, int num, struct subgroup_properties
 *                  *subgroup, int grnr, int subnr, int parallel_flag, int
 *                  nsubgroups_cat)
 *
 *
 * \par Major modifications and contributions:
 *
 * - DD.MM.YYYY Description
 * - 14.05.2018 Prepared file for public release -- Rainer Weinberger
 */

#include <math.h>
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "../main/allvars.h"
#include "../main/proto.h"

#ifdef SUBFIND
#include "../fof/fof.h"
#include "subfind.h"

/*! \brief Calculates subhalo properties.
 *
 *
 *  \param[in] d Unbind data.
 *  \param[in] num Length of d.
 *  \param[out] subgroup Data for subgroup properties.
 *  \param[in] grnr Index in GroupCat.
 *  \param[in] subnr Index of Subhalo in this group.
 *  \param[in] parallel_flag If set, the code calculates the properties for a
 *             subhalo distributed onto several processors.
 *  \param[in] nsubgroups_cat (unused)
 *
 *  \return void
 */
void subfind_determine_sub_halo_properties(struct unbind_data *d, int num, struct subgroup_properties *subgroup, int grnr, int subnr,
                                           int parallel_flag, int nsubgroups_cat)
{
  int i, j, p, len_type[NTYPES], len_type_loc[NTYPES], totlen;
  double s[3], v[3], pos[3], vel[3], spin[3], cm[3], veldisp, max, vel_to_phys, H_of_a, minpot;
#ifdef MHD
  double bfld_halo, bfld_disk, bfld_vol_halo, bfld_vol_disk;
#endif /* #ifdef MHD */
#ifdef SUBFIND_EXTENDED_PROPERTIES
  double Ekin = 0, Epot = 0, Ethr = 0, Jdm[3], Jgas[3], Jstars[3], CMFrac, CMFracType[NTYPES];
  double Jdm_inHalfRad[3], Jgas_inHalfRad[3], Jstars_inHalfRad[3], CMFrac_inHalfRad, CMFracType_inHalfRad[NTYPES];
  double Jdm_inRad[3], Jgas_inRad[3], Jstars_inRad[3], CMFrac_inRad, CMFracType_inRad[NTYPES];
  double jpart[3], Jtot[3], Jtot_inRad[3], Jtot_inHalfRad[3];
  double sinrad[3], sinhalfrad[3], vinrad[3], vinhalfrad[3];
#endif /* #ifdef SUBFIND_EXTENDED_PROPERTIES */
  double lx, ly, lz, dv[3], dx[3], disp, rr_tmp, disp_tmp, halfmassrad = 0, halfmassradtype[NTYPES];
  double boxsize, ddxx, vmax, vmaxrad, maxrad;
  double mass, massinrad, massinhalfrad, massinmaxrad;
  double mass_tab[NTYPES], massinrad_tab[NTYPES], massinhalfrad_tab[NTYPES], massinmaxrad_tab[NTYPES];
  double xtmp;

  sort_r2list *rr_list = 0;
  int minindex;
  MyIDType mostboundid;

#ifdef USE_SFR
  double sfr = 0, sfrinrad = 0, sfrinhalfrad = 0, sfrinmaxrad = 0, gasMassSfr = 0;
#endif /* #ifdef USE_SFR */

  boxsize = All.BoxSize;

  vel_to_phys = 1.0 / All.cf_atime;

  if(All.ComovingIntegrationOn)
    H_of_a = hubble_function(All.Time);
  else
    H_of_a = 0;

  mass = massinrad = massinhalfrad = massinmaxrad = 0;
  for(j = 0; j < NTYPES; j++)
    {
      len_type[j] = 0;
      mass_tab[j] = halfmassradtype[j] = massinrad_tab[j] = massinhalfrad_tab[j] = massinmaxrad_tab[j] = 0;
    }

  for(i = 0, minindex = -1, minpot = 1.0e30; i < num; i++)
    {
      p = d[i].index;
      if(PS[p].Potential < minpot || minindex == -1)
        {
          minpot   = PS[p].Potential;
          minindex = p;
        }

      len_type[P[p].Type]++;

#ifdef USE_SFR
      if(P[p].Type == 0)
        sfr += SphP[PS[p].OldIndex].Sfr; /* note: the SphP[] array has not been reordered */
#endif                                   /* #ifdef USE_SFR */
    }

  for(j = 0; j < NTYPES; j++)
    len_type_loc[j] = len_type[j];

  if(parallel_flag)
    {
      int len_typetot[NTYPES];
      MPI_Allreduce(len_type, len_typetot, NTYPES, MPI_INT, MPI_SUM, SubComm);
      for(j = 0; j < NTYPES; j++)
        len_type[j] = len_typetot[j];

      double *minpotlist = mymalloc("minpotlist", SubNTask * sizeof(double));
      MPI_Allgather(&minpot, 1, MPI_DOUBLE, minpotlist, 1, MPI_DOUBLE, SubComm);
      int mincpu;

      for(i = 0, mincpu = -1, minpot = 1.0e30; i < SubNTask; i++)
        if(minpotlist[i] < minpot)
          {
            mincpu = i;
            minpot = minpotlist[mincpu];
          }

      myfree(minpotlist);

      if(mincpu < 0)
        terminate("mincpu < 0");

      if(SubThisTask == mincpu)
        for(j = 0; j < 3; j++)
          {
#ifdef CELL_CENTER_GRAVITY
            if(P[minindex].Type == 0)
              pos[j] = SphP[PS[minindex].OldIndex].Center[j];
            else
#endif /* #ifdef CELL_CENTER_GRAVITY */
              pos[j] = P[minindex].Pos[j];
          }

      MPI_Bcast(pos, 3, MPI_DOUBLE, mincpu, SubComm);

#ifdef USE_SFR
      double sfrtot;
      MPI_Allreduce(&sfr, &sfrtot, 1, MPI_DOUBLE, MPI_SUM, SubComm);
      sfr = sfrtot;
#endif /* #ifdef USE_SFR */
    }
  else
    {
      if(minindex == -1)
        terminate("minindex == -1");

      for(j = 0; j < 3; j++)
        {
#ifdef CELL_CENTER_GRAVITY
          if(P[minindex].Type == 0)
            pos[j] = SphP[PS[minindex].OldIndex].Center[j];
          else
#endif /* #ifdef CELL_CENTER_GRAVITY */
            pos[j] = P[minindex].Pos[j];
        }
    }

  /* pos[] now holds the position of minimum potential */
  /* we'll take it that as the center */

  /* determine the particle ID with the smallest binding energy */
  for(i = 0, minindex = -1, minpot = 1.0e30; i < num; i++)
    {
      p = d[i].index;
      if(PS[p].BindingEnergy < minpot || minindex == -1)
        {
          minpot   = PS[p].BindingEnergy;
          minindex = p;
        }
    }

  if(parallel_flag)
    {
      double *minpotlist = mymalloc("minpotlist", SubNTask * sizeof(double));
      MPI_Allgather(&minpot, 1, MPI_DOUBLE, minpotlist, 1, MPI_DOUBLE, SubComm);
      int mincpu;

      for(i = 0, mincpu = -1, minpot = 1.0e30; i < SubNTask; i++)
        if(minpotlist[i] < minpot)
          {
            mincpu = i;
            minpot = minpotlist[mincpu];
          }

      myfree(minpotlist);

      if(mincpu < 0)
        terminate("mincpu < 0");

      if(SubThisTask == mincpu)
        {
          mostboundid = P[minindex].ID;
        }

      MPI_Bcast(&mostboundid, sizeof(mostboundid), MPI_BYTE, mincpu, SubComm);
    }
  else
    {
      if(minindex == -1)
        terminate("minindex == -1");

      mostboundid = P[minindex].ID;
    }

  /* let's get bulk velocity and the center-of-mass */
  /* here we still take all particles */

  for(j = 0; j < 3; j++)
    s[j] = v[j] = 0;

  for(i = 0; i < num; i++)
    {
      p = d[i].index;
      for(j = 0; j < 3; j++)
        {
          ddxx = GRAVITY_NEAREST_X(P[p].Pos[j] - pos[j]);
          s[j] += P[p].Mass * ddxx;
          v[j] += P[p].Mass * P[p].Vel[j];
        }
      mass += P[p].Mass;

      int ptype = P[p].Type;
      mass_tab[ptype] += P[p].Mass;
    }

  if(parallel_flag)
    {
      double stot[3], vtot[3], masstot, mass_tabtot[NTYPES];

      MPI_Allreduce(s, stot, 3, MPI_DOUBLE, MPI_SUM, SubComm);
      MPI_Allreduce(&mass, &masstot, 1, MPI_DOUBLE, MPI_SUM, SubComm);
      MPI_Allreduce(v, vtot, 3, MPI_DOUBLE, MPI_SUM, SubComm);
      MPI_Allreduce(mass_tab, mass_tabtot, NTYPES, MPI_DOUBLE, MPI_SUM, SubComm);

      mass = masstot;
      for(j = 0; j < 3; j++)
        {
          s[j] = stot[j];
          v[j] = vtot[j];
        }

      for(j = 0; j < NTYPES; j++)
        mass_tab[j] = mass_tabtot[j];
    }

  for(j = 0; j < 3; j++)
    {
      s[j] /= mass; /* center of mass */
      v[j] /= mass;
      vel[j] = vel_to_phys * v[j];
    }

  for(j = 0; j < 3; j++)
    {
      s[j] += pos[j];

      while(s[j] < 0)
        s[j] += boxsize;
      while(s[j] >= boxsize)
        s[j] -= boxsize;
      cm[j] = s[j];  // this is in comoving coordinates
    }

  disp = lx = ly = lz = 0;
#ifdef SUBFIND_EXTENDED_PROPERTIES
  Jtot[0] = Jtot[1] = Jtot[2] = 0;
  Jdm[0] = Jdm[1] = Jdm[2] = 0;
  Jgas[0] = Jgas[1] = Jgas[2] = 0;
  Jstars[0] = Jstars[1] = Jstars[2] = 0;
#endif /* #ifdef SUBFIND_EXTENDED_PROPERTIES */

  rr_list = mymalloc("rr_list", sizeof(sort_r2list) * (num + 1));

  for(i = 0; i < num; i++)
    {
      p = d[i].index;

      for(j = 0, rr_tmp = 0, disp_tmp = 0; j < 3; j++)
        {
          ddxx  = GRAVITY_NEAREST_X(P[p].Pos[j] - s[j]);
          dx[j] = All.cf_atime * ddxx;
          dv[j] = vel_to_phys * (P[p].Vel[j] - v[j]);
          dv[j] += H_of_a * dx[j];

          disp_tmp += P[p].Mass * dv[j] * dv[j];
          /* for rotation curve computation, take minimum of potential as center */
          ddxx = GRAVITY_NEAREST_X(P[p].Pos[j] - pos[j]);
          ddxx = All.cf_atime * ddxx;
          rr_tmp += ddxx * ddxx;
        }

      lx += P[p].Mass * (dx[1] * dv[2] - dx[2] * dv[1]);
      ly += P[p].Mass * (dx[2] * dv[0] - dx[0] * dv[2]);
      lz += P[p].Mass * (dx[0] * dv[1] - dx[1] * dv[0]);

#ifdef SUBFIND_EXTENDED_PROPERTIES
      for(j = 0; j < 3; j++)  // hubble drifts in velocity now with respect to pot min which we consider as the centre of rotation
        {
          ddxx  = GRAVITY_NEAREST_X(P[p].Pos[j] - pos[j]);
          dx[j] = All.cf_atime * ddxx;
          dv[j] = vel_to_phys * (P[p].Vel[j] - v[j]);
          dv[j] += H_of_a * dx[j];
        }

      int ptype = P[p].Type;

      Ekin += (P[p].Mass / 2) * (dv[0] * dv[0] + dv[1] * dv[1] + dv[2] * dv[2]);
      Epot += (P[p].Mass / 2) * PS[p].Potential;
      if(P[p].Type == 0)
        Ethr += P[p].Mass * SphP[PS[p].OldIndex].Utherm;

      Jtot[0] += P[p].Mass * (dx[1] * dv[2] - dx[2] * dv[1]);
      Jtot[1] += P[p].Mass * (dx[2] * dv[0] - dx[0] * dv[2]);
      Jtot[2] += P[p].Mass * (dx[0] * dv[1] - dx[1] * dv[0]);

      if(ptype == 1)  // dm illustris
        {
          Jdm[0] += P[p].Mass * (dx[1] * dv[2] - dx[2] * dv[1]);
          Jdm[1] += P[p].Mass * (dx[2] * dv[0] - dx[0] * dv[2]);
          Jdm[2] += P[p].Mass * (dx[0] * dv[1] - dx[1] * dv[0]);
        }
      if(ptype == 0)  // gas (incl. winds!)
        {
          Jgas[0] += P[p].Mass * (dx[1] * dv[2] - dx[2] * dv[1]);
          Jgas[1] += P[p].Mass * (dx[2] * dv[0] - dx[0] * dv[2]);
          Jgas[2] += P[p].Mass * (dx[0] * dv[1] - dx[1] * dv[0]);
        }
      if(ptype == 4)  // stars (previously: StarP[P[p].AuxDataID].BirthTime)
        {
          Jstars[0] += P[p].Mass * (dx[1] * dv[2] - dx[2] * dv[1]);
          Jstars[1] += P[p].Mass * (dx[2] * dv[0] - dx[0] * dv[2]);
          Jstars[2] += P[p].Mass * (dx[0] * dv[1] - dx[1] * dv[0]);
        }
#endif /* #ifdef SUBFIND_EXTENDED_PROPERTIES */

      rr_tmp = sqrt(rr_tmp);

      rr_list[i].mass = P[p].Mass;
      rr_list[i].r    = rr_tmp;
      disp += disp_tmp;
    }

  if(parallel_flag)
    {
      double spintot[3], disptot;
      spin[0] = lx;
      spin[1] = ly;
      spin[2] = lz;
      MPI_Allreduce(spin, spintot, 3, MPI_DOUBLE, MPI_SUM, SubComm);
      MPI_Allreduce(&disp, &disptot, 1, MPI_DOUBLE, MPI_SUM, SubComm);
      disp = disptot;
      lx   = spintot[0];
      ly   = spintot[1];
      lz   = spintot[2];
#ifdef SUBFIND_EXTENDED_PROPERTIES
      MPI_Allreduce(MPI_IN_PLACE, &Ekin, 1, MPI_DOUBLE, MPI_SUM, SubComm);
      MPI_Allreduce(MPI_IN_PLACE, &Epot, 1, MPI_DOUBLE, MPI_SUM, SubComm);
      MPI_Allreduce(MPI_IN_PLACE, &Ethr, 1, MPI_DOUBLE, MPI_SUM, SubComm);
      MPI_Allreduce(MPI_IN_PLACE, Jtot, 3, MPI_DOUBLE, MPI_SUM, SubComm);
      MPI_Allreduce(MPI_IN_PLACE, Jdm, 3, MPI_DOUBLE, MPI_SUM, SubComm);
      MPI_Allreduce(MPI_IN_PLACE, Jgas, 3, MPI_DOUBLE, MPI_SUM, SubComm);
      MPI_Allreduce(MPI_IN_PLACE, Jstars, 3, MPI_DOUBLE, MPI_SUM, SubComm);
#endif /* #ifdef SUBFIND_EXTENDED_PROPERTIES */
    }

  spin[0] = lx / mass;
  spin[1] = ly / mass;
  spin[2] = lz / mass;

  veldisp = sqrt(disp / (3 * mass)); /* convert to 1d velocity dispersion */

#ifdef SUBFIND_EXTENDED_PROPERTIES
  // counter rotating mass fractions
  CMFrac = 0;
  for(i = 0; i < NTYPES; i++)
    CMFracType[i] = 0;

  for(i = 0; i < num; i++)
    {
      /* identify particle type */
      p = d[i].index;

      /* calculate particle radius */
      for(j = 0; j < 3; j++)
        {
          ddxx  = GRAVITY_NEAREST_X(P[p].Pos[j] - pos[j]);  // counter-rotating mass calc with respect to pot min
          dx[j] = All.cf_atime * ddxx;
          dv[j] = vel_to_phys * (P[p].Vel[j] - v[j]);
          dv[j] += H_of_a * dx[j];
        }

      int ptype = P[p].Type;

      jpart[0] = P[p].Mass * (dx[1] * dv[2] - dx[2] * dv[1]);
      jpart[1] = P[p].Mass * (dx[2] * dv[0] - dx[0] * dv[2]);
      jpart[2] = P[p].Mass * (dx[0] * dv[1] - dx[1] * dv[0]);

      if((Jtot[0] * jpart[0] + Jtot[1] * jpart[1] + Jtot[2] * jpart[2]) < 0.)
        CMFrac += P[p].Mass / mass;

      if(ptype == 1)  // dm illustris
        if((Jdm[0] * jpart[0] + Jdm[1] * jpart[1] + Jdm[2] * jpart[2]) < 0.)
          CMFracType[1] += P[p].Mass / mass_tab[1];
      if(ptype == 0)  // gas (incl. winds!)
        if((Jgas[0] * jpart[0] + Jgas[1] * jpart[1] + Jgas[2] * jpart[2]) < 0.)
          CMFracType[0] += P[p].Mass / mass_tab[0];
      if(ptype == 4)  // stars
        if((Jstars[0] * jpart[0] + Jstars[1] * jpart[1] + Jstars[2] * jpart[2]) < 0.)
          CMFracType[4] += P[p].Mass / mass_tab[4];
    }

  if(parallel_flag)
    {
      MPI_Allreduce(MPI_IN_PLACE, &CMFrac, 1, MPI_DOUBLE, MPI_SUM, SubComm);
      MPI_Allreduce(MPI_IN_PLACE, CMFracType, NTYPES, MPI_DOUBLE, MPI_SUM, SubComm);
    }

#endif /* #ifdef SUBFIND_EXTENDED_PROPERTIES */

  if(parallel_flag)
    parallel_sort_comm(rr_list, num, sizeof(sort_r2list), subfind_compare_dist_rotcurve, SubComm);
  else
    mysort(rr_list, num, sizeof(sort_r2list), subfind_compare_dist_rotcurve);

  /* calculate cumulative mass */
  for(i = 1; i < num; i++)
    rr_list[i].mass += rr_list[i - 1].mass;

  if(parallel_flag)
    {
      double mass_part = 0;
      if(num)
        mass_part = rr_list[num - 1].mass;
      double *masslist = mymalloc("masslist", SubNTask * sizeof(double));
      MPI_Allgather(&mass_part, 1, MPI_DOUBLE, masslist, 1, MPI_DOUBLE, SubComm);

      double massbefore = 0;
      for(i = 0; i < SubThisTask; i++)
        massbefore += masslist[i];

      for(i = 0; i < num; i++)
        rr_list[i].mass += massbefore;

      myfree(masslist);

      /* now calculate rotation curve maximum and half mass radius */

      double halfmassrad_loc  = 0;
      sort_r2list *rr_lowlist = mymalloc("rr_lowlist", SubNTask * sizeof(sort_r2list));
      sort_r2list low_element;
      if(num > 0)
        low_element = rr_list[0];
      else
        {
          low_element.mass = 0;
          low_element.r    = 0;
        }
      MPI_Allgather(&low_element, sizeof(sort_r2list), MPI_BYTE, rr_lowlist, sizeof(sort_r2list), MPI_BYTE, SubComm);

      rr_list[num].mass = 0;
      rr_list[num].r    = 0;

      for(j = SubThisTask + 1; j < SubNTask; j++)
        if(rr_lowlist[j].mass > 0)
          {
            rr_list[num] = rr_lowlist[j];
            break;
          }

      myfree(rr_lowlist);

      int *numlist = mymalloc("numlist", SubNTask * sizeof(int));
      MPI_Allgather(&num, 1, MPI_INT, numlist, 1, MPI_INT, SubComm);

      int nbefore = 0;
      for(i = 0; i < SubThisTask; i++)
        nbefore += numlist[i];

      for(i = num - 1, max = 0, maxrad = 0; i >= 0; i--)
        {
          if((i + nbefore) > 5 && rr_list[i].mass > max * rr_list[i].r)
            {
              max    = rr_list[i].mass / rr_list[i].r;
              maxrad = rr_list[i].r;
            }

          if(rr_list[i].mass < 0.5 * mass && rr_list[i + 1].mass >= 0.5 * mass)
            halfmassrad_loc = 0.5 * (rr_list[i].r + rr_list[i + 1].r);
        }

      myfree(numlist);

      MPI_Allreduce(&halfmassrad_loc, &halfmassrad, 1, MPI_DOUBLE, MPI_MAX, SubComm);
      double *maxlist    = mymalloc("maxlist", SubNTask * sizeof(double));
      double *maxradlist = mymalloc("maxradlist", SubNTask * sizeof(double));
      MPI_Allgather(&max, 1, MPI_DOUBLE, maxlist, 1, MPI_DOUBLE, SubComm);
      MPI_Allgather(&maxrad, 1, MPI_DOUBLE, maxradlist, 1, MPI_DOUBLE, SubComm);
      for(i = 0, max = maxrad = 0; i < SubNTask; i++)
        {
          if(maxlist[i] > max)
            {
              max    = maxlist[i];
              maxrad = maxradlist[i];
            }
        }
      myfree(maxradlist);
      myfree(maxlist);
    }
  else
    {
      for(i = num - 1, max = 0, maxrad = 0; i >= 0; i--)
        {
          if(i > 5 && rr_list[i].mass > max * rr_list[i].r)
            {
              max    = rr_list[i].mass / rr_list[i].r;
              maxrad = rr_list[i].r;
            }

          if(i < num - 1)
            if(rr_list[i].mass < 0.5 * mass && rr_list[i + 1].mass >= 0.5 * mass)
              halfmassrad = 0.5 * (rr_list[i].r + rr_list[i + 1].r);
        }
    }

  halfmassrad /= All.cf_atime;
  vmax    = sqrt(All.G * max);
  vmaxrad = maxrad / All.cf_atime;

  myfree(rr_list);

  /* half mass radii for different types */
  /* need to recalculate len_type_loc first, because of special particle treatment in GFM */
  for(j = 0; j < NTYPES; j++)
    len_type_loc[j] = 0;

  for(i = 0; i < num; i++)
    {
      p         = d[i].index;
      int ptype = P[p].Type;

      len_type_loc[ptype]++;
    }

  int itmp, type;
  for(type = 0; type < NTYPES; type++)
    {
      rr_list = mymalloc("rr_list", sizeof(sort_r2list) * (len_type_loc[type] + 1));
      itmp    = 0;
      for(i = 0; i < num; i++)
        {
          p = d[i].index;

          int ptype = P[p].Type;

          if(ptype == type)
            {
              for(j = 0, rr_tmp = 0; j < 3; j++)
                {
                  ddxx = GRAVITY_NEAREST_X(P[p].Pos[j] - pos[j]);
                  rr_tmp += ddxx * ddxx;
                }

              rr_tmp = sqrt(rr_tmp);

              rr_list[itmp].mass = P[p].Mass;
              rr_list[itmp].r    = rr_tmp;
              itmp++;
            }
        }

      if(itmp != len_type_loc[type])
        terminate("should not occur: %d %d", itmp, len_type_loc[type]);

      if(parallel_flag)
        parallel_sort_comm(rr_list, len_type_loc[type], sizeof(sort_r2list), subfind_compare_dist_rotcurve, SubComm);
      else
        mysort(rr_list, len_type_loc[type], sizeof(sort_r2list), subfind_compare_dist_rotcurve);

      /* calculate cumulative mass */
      for(i = 1; i < len_type_loc[type]; i++)
        rr_list[i].mass = rr_list[i - 1].mass + rr_list[i].mass;

      if(parallel_flag)
        {
          double mass_part = 0;
          if(len_type_loc[type])
            mass_part = rr_list[len_type_loc[type] - 1].mass;
          double *masslist = mymalloc("masslist", SubNTask * sizeof(double));
          MPI_Allgather(&mass_part, 1, MPI_DOUBLE, masslist, 1, MPI_DOUBLE, SubComm);

          double massbefore = 0;
          for(i = 0; i < SubThisTask; i++)
            massbefore += masslist[i];

          for(i = 0; i < len_type_loc[type]; i++)
            rr_list[i].mass += massbefore;

          myfree(masslist);
        }

      /* now calculate half mass radii */
      if(parallel_flag)
        {
          double halfmassrad_loc  = 0;
          sort_r2list *rr_lowlist = mymalloc("rr_lowlist", SubNTask * sizeof(sort_r2list));
          sort_r2list low_element;
          if(len_type_loc[type] > 0)
            low_element = rr_list[0];
          else
            {
              low_element.mass = 0;
              low_element.r    = 0;
            }

          MPI_Allgather(&low_element, sizeof(sort_r2list), MPI_BYTE, rr_lowlist, sizeof(sort_r2list), MPI_BYTE, SubComm);

          rr_list[len_type_loc[type]].mass = 0;
          rr_list[len_type_loc[type]].r    = 0;
          for(j = SubThisTask + 1; j < SubNTask; j++)
            if(rr_lowlist[j].mass > 0)
              {
                rr_list[len_type_loc[type]] = rr_lowlist[j];
                break;
              }

          myfree(rr_lowlist);

          for(i = len_type_loc[type] - 1; i >= 0; i--)
            {
              if(rr_list[i].mass < 0.5 * mass_tab[type] && rr_list[i + 1].mass >= 0.5 * mass_tab[type])
                halfmassrad_loc = 0.5 * (rr_list[i].r + rr_list[i + 1].r);
            }

          MPI_Allreduce(&halfmassrad_loc, &halfmassradtype[type], 1, MPI_DOUBLE, MPI_MAX, SubComm);
        }
      else
        {
          for(i = len_type_loc[type] - 1; i >= 0; i--)
            {
              if(i < len_type_loc[type] - 1)
                if(rr_list[i].mass < 0.5 * mass_tab[type] && rr_list[i + 1].mass >= 0.5 * mass_tab[type])
                  halfmassradtype[type] = 0.5 * (rr_list[i].r + rr_list[i + 1].r);
            }
        }

      myfree(rr_list);
    }

    /* properties of 'central galaxies', defined in several ways as particles within some radius:
       either (stellar half mass radius) or SUBFIND_GAL_RADIUS_FAC*(stellar half mass radius) or (radius of Vmax) */
#ifdef SUBFIND_EXTENDED_PROPERTIES
  // centre of mass /velocity of particles in half/ stellar mass rad
  sinrad[0] = sinrad[1] = sinrad[2] = 0;
  sinhalfrad[0] = sinhalfrad[1] = sinhalfrad[2] = 0;
  vinrad[0] = vinrad[1] = vinrad[2] = 0;
  vinhalfrad[0] = vinhalfrad[1] = vinhalfrad[2] = 0;
#endif /* #ifdef SUBFIND_EXTENDED_PROPERTIES */

  for(i = 0; i < num; i++)
    {
      /* identify particle type */
      p         = d[i].index;
      int ptype = P[p].Type;

      /* calculate particle radius */
      for(j = 0, rr_tmp = 0; j < 3; j++)
        {
          ddxx = GRAVITY_NEAREST_X(P[p].Pos[j] - pos[j]);
          rr_tmp += ddxx * ddxx;
        }
      rr_tmp = sqrt(rr_tmp);

      /* properties inside SUBFIND_GAL_RADIUS_FAC*(stellar half mass radius) */
      if(rr_tmp < SUBFIND_GAL_RADIUS_FAC * halfmassradtype[4])
        {
          massinrad += P[p].Mass;
          massinrad_tab[ptype] += P[p].Mass;

#ifdef SUBFIND_EXTENDED_PROPERTIES
          for(j = 0; j < 3; j++)
            {
              ddxx = GRAVITY_NEAREST_X(P[p].Pos[j] - pos[j]);  // comoving (as it should be.)
              sinrad[j] += P[p].Mass * ddxx;
              vinrad[j] += P[p].Mass * P[p].Vel[j];
            }
#endif /* #ifdef SUBFIND_EXTENDED_PROPERTIES */

          if(ptype == 0)
            {
              if(P[p].Type == 0)
                {
#ifdef USE_SFR
                  sfrinrad += SphP[PS[p].OldIndex].Sfr; /* note: the SphP[] array has not been reordered */
#endif                                                  /* #ifdef USE_SFR */
                }
            }
        }

      /* properties inside (stellar half mass radius) */
      if(rr_tmp < 1.0 * halfmassradtype[4])
        {
          massinhalfrad += P[p].Mass;
          massinhalfrad_tab[ptype] += P[p].Mass;

#ifdef SUBFIND_EXTENDED_PROPERTIES
          for(j = 0; j < 3; j++)
            {
              ddxx = GRAVITY_NEAREST_X(P[p].Pos[j] - pos[j]);  // comoving (as it should be.)
              sinhalfrad[j] += P[p].Mass * ddxx;
              vinhalfrad[j] += P[p].Mass * P[p].Vel[j];
            }
#endif /* #ifdef SUBFIND_EXTENDED_PROPERTIES */

          if(ptype == 0)
            {
              if(P[p].Type == 0)
                {
#ifdef USE_SFR
                  sfrinhalfrad += SphP[PS[p].OldIndex].Sfr; /* note: the SphP[] array has not been reordered */
#endif                                                      /* #ifdef USE_SFR */
                }
            }
        }

      /* properties inside (radius of Vmax) */
      if(rr_tmp < 1.0 * vmaxrad)
        {
          massinmaxrad += P[p].Mass;
          massinmaxrad_tab[ptype] += P[p].Mass;

          if(ptype == 0)
            {
              if(P[p].Type == 0)
                {
#ifdef USE_SFR
                  sfrinmaxrad += SphP[PS[p].OldIndex].Sfr; /* note: the SphP[] array has not been reordered */
#endif                                                     /* #ifdef USE_SFR */
                }
            }
        }
    }

    /* properties of star forming gas */
#ifdef USE_SFR
  for(i = 0; i < num; i++)
    {
      p = d[i].index;

      if(P[p].Type == 0)
        {
          if(SphP[PS[p].OldIndex].Sfr > 0)
            {
              gasMassSfr += P[p].Mass;
            }
        }
    }
#endif /* #ifdef USE_SFR */

#ifdef MHD
  bfld_halo = bfld_disk = bfld_vol_halo = bfld_vol_disk = 0;

  for(i = 0; i < num; i++)
    {
      p = d[i].index;

      if(P[p].Type == 0)
        {
          double bfld2 = (SphP[PS[p].OldIndex].B[0] * SphP[PS[p].OldIndex].B[0]) +
                         (SphP[PS[p].OldIndex].B[1] * SphP[PS[p].OldIndex].B[1]) +
                         (SphP[PS[p].OldIndex].B[2] * SphP[PS[p].OldIndex].B[2]);
          double vol = SphP[PS[p].OldIndex].Volume;

          bfld_halo += bfld2 * vol;
          bfld_vol_halo += vol;

          /* calculate particle radius */
          for(j = 0, rr_tmp = 0; j < 3; j++)
            {
              ddxx = GRAVITY_NEAREST_X(P[p].Pos[j] - pos[j]);
              rr_tmp += ddxx * ddxx;
            }
          rr_tmp = sqrt(rr_tmp);

          if(rr_tmp < SUBFIND_GAL_RADIUS_FAC * halfmassradtype[4])
            {
              bfld_disk += bfld2 * vol;
              bfld_vol_disk += vol;
            }
        }
    }
#endif /* #ifdef MHD */

  if(parallel_flag)
    {
      double massinradtot, massinrad_tabtot[NTYPES];
      MPI_Allreduce(&massinrad, &massinradtot, 1, MPI_DOUBLE, MPI_SUM, SubComm);
      MPI_Allreduce(massinrad_tab, massinrad_tabtot, NTYPES, MPI_DOUBLE, MPI_SUM, SubComm);
      massinrad = massinradtot;
      for(j = 0; j < NTYPES; j++)
        massinrad_tab[j] = massinrad_tabtot[j];

      double massinhalfradtot, massinhalfrad_tabtot[NTYPES];
      MPI_Allreduce(&massinhalfrad, &massinhalfradtot, 1, MPI_DOUBLE, MPI_SUM, SubComm);
      MPI_Allreduce(massinhalfrad_tab, massinhalfrad_tabtot, NTYPES, MPI_DOUBLE, MPI_SUM, SubComm);
      massinhalfrad = massinhalfradtot;
      for(j = 0; j < NTYPES; j++)
        massinhalfrad_tab[j] = massinhalfrad_tabtot[j];

      double massinmaxradtot, massinmaxrad_tabtot[NTYPES];
      MPI_Allreduce(&massinmaxrad, &massinmaxradtot, 1, MPI_DOUBLE, MPI_SUM, SubComm);
      MPI_Allreduce(massinmaxrad_tab, massinmaxrad_tabtot, NTYPES, MPI_DOUBLE, MPI_SUM, SubComm);
      massinmaxrad = massinmaxradtot;
      for(j = 0; j < NTYPES; j++)
        massinmaxrad_tab[j] = massinmaxrad_tabtot[j];

#ifdef SUBFIND_EXTENDED_PROPERTIES
      MPI_Allreduce(MPI_IN_PLACE, sinrad, 3, MPI_DOUBLE, MPI_SUM, SubComm);
      MPI_Allreduce(MPI_IN_PLACE, vinrad, 3, MPI_DOUBLE, MPI_SUM, SubComm);
      MPI_Allreduce(MPI_IN_PLACE, sinhalfrad, 3, MPI_DOUBLE, MPI_SUM, SubComm);
      MPI_Allreduce(MPI_IN_PLACE, vinhalfrad, 3, MPI_DOUBLE, MPI_SUM, SubComm);
#endif /* #ifdef SUBFIND_EXTENDED_PROPERTIES */

#ifdef MHD
      double bfld_halo_tot, bfld_disk_tot, bfld_vol_halo_tot, bfld_vol_disk_tot;
      MPI_Allreduce(&bfld_halo, &bfld_halo_tot, 1, MPI_DOUBLE, MPI_SUM, SubComm);
      MPI_Allreduce(&bfld_vol_halo, &bfld_vol_halo_tot, 1, MPI_DOUBLE, MPI_SUM, SubComm);
      MPI_Allreduce(&bfld_disk, &bfld_disk_tot, 1, MPI_DOUBLE, MPI_SUM, SubComm);
      MPI_Allreduce(&bfld_vol_disk, &bfld_vol_disk_tot, 1, MPI_DOUBLE, MPI_SUM, SubComm);

      bfld_halo     = bfld_halo_tot;
      bfld_vol_halo = bfld_vol_halo_tot;
      bfld_disk     = bfld_disk_tot;
      bfld_vol_disk = bfld_vol_disk_tot;
#endif /* #ifdef MHD */

#ifdef USE_SFR
      double sfrinradtot;
      MPI_Allreduce(&sfrinrad, &sfrinradtot, 1, MPI_DOUBLE, MPI_SUM, SubComm);
      sfrinrad = sfrinradtot;

      double sfrinhalfradtot;
      MPI_Allreduce(&sfrinhalfrad, &sfrinhalfradtot, 1, MPI_DOUBLE, MPI_SUM, SubComm);
      sfrinhalfrad = sfrinhalfradtot;

      double sfrinmaxradtot;
      MPI_Allreduce(&sfrinmaxrad, &sfrinmaxradtot, 1, MPI_DOUBLE, MPI_SUM, SubComm);
      sfrinmaxrad = sfrinmaxradtot;

      double gasMassSfrtot;
      MPI_Allreduce(&gasMassSfr, &gasMassSfrtot, 1, MPI_DOUBLE, MPI_SUM, SubComm);
      gasMassSfr = gasMassSfrtot;
#endif /* #ifdef USE_SFR */
    }

  if(parallel_flag)
    MPI_Allreduce(&num, &totlen, 1, MPI_INT, MPI_SUM, SubComm);
  else
    totlen = num;

#ifdef MHD
  if(bfld_vol_halo > 0.)
    bfld_halo = sqrt(bfld_halo / bfld_vol_halo);
  if(bfld_vol_disk > 0.)
    bfld_disk = sqrt(bfld_disk / bfld_vol_disk);
#endif /* #ifdef MHD */

#ifdef SUBFIND_EXTENDED_PROPERTIES
  // finish centre of mass of spheres
  for(j = 0; j < 3; j++)
    {
      if(massinrad > 0)
        {
          sinrad[j] /= massinrad;
          sinrad[j] += pos[j];

          while(sinrad[j] < 0)
            sinrad[j] += boxsize;
          while(sinrad[j] >= boxsize)
            sinrad[j] -= boxsize;

          vinrad[j] /= massinrad;  // this is comoving (as it should be.)
        }

      if(massinhalfrad > 0)
        {
          sinhalfrad[j] /= massinhalfrad;
          sinhalfrad[j] += pos[j];

          while(sinhalfrad[j] < 0)
            sinhalfrad[j] += boxsize;
          while(sinhalfrad[j] >= boxsize)
            sinhalfrad[j] -= boxsize;

          vinhalfrad[j] /= massinhalfrad;
        }
    }

  Jtot_inHalfRad[0] = Jtot_inHalfRad[1] = Jtot_inHalfRad[2] = 0;
  Jdm_inHalfRad[0] = Jdm_inHalfRad[1] = Jdm_inHalfRad[2] = 0;
  Jgas_inHalfRad[0] = Jgas_inHalfRad[1] = Jgas_inHalfRad[2] = 0;
  Jstars_inHalfRad[0] = Jstars_inHalfRad[1] = Jstars_inHalfRad[2] = 0;
  Jtot_inRad[0] = Jtot_inRad[1] = Jtot_inRad[2] = 0;
  Jdm_inRad[0] = Jdm_inRad[1] = Jdm_inRad[2] = 0;
  Jgas_inRad[0] = Jgas_inRad[1] = Jgas_inRad[2] = 0;
  Jstars_inRad[0] = Jstars_inRad[1] = Jstars_inRad[2] = 0;

  for(i = 0; i < num; i++)
    {
      /* identify particle type */
      p = d[i].index;

      /* calculate particle radius */
      for(j = 0, rr_tmp = 0; j < 3; j++)
        {
          ddxx = GRAVITY_NEAREST_X(P[p].Pos[j] - pos[j]);
          rr_tmp += ddxx * ddxx;
        }
      rr_tmp = sqrt(rr_tmp);

      int ptype = P[p].Type;

      /* properties inside SUBFIND_GAL_RADIUS_FAC*(stellar half mass radius) */
      if((massinrad > 0) && (rr_tmp < SUBFIND_GAL_RADIUS_FAC * halfmassradtype[4]))
        {
          for(j = 0; j < 3; j++)
            {
              ddxx  = GRAVITY_NEAREST_X(P[p].Pos[j] - pos[j]);
              dx[j] = All.cf_atime * ddxx;
              dv[j] = vel_to_phys * (P[p].Vel[j] - vinrad[j]);
              dv[j] += H_of_a * dx[j];
            }

          Jtot_inRad[0] += P[p].Mass * (dx[1] * dv[2] - dx[2] * dv[1]);
          Jtot_inRad[1] += P[p].Mass * (dx[2] * dv[0] - dx[0] * dv[2]);
          Jtot_inRad[2] += P[p].Mass * (dx[0] * dv[1] - dx[1] * dv[0]);

          if(ptype == 1)  // dm illustris
            {
              Jdm_inRad[0] += P[p].Mass * (dx[1] * dv[2] - dx[2] * dv[1]);
              Jdm_inRad[1] += P[p].Mass * (dx[2] * dv[0] - dx[0] * dv[2]);
              Jdm_inRad[2] += P[p].Mass * (dx[0] * dv[1] - dx[1] * dv[0]);
            }
          if(ptype == 0)  // gas
            {
              Jgas_inRad[0] += P[p].Mass * (dx[1] * dv[2] - dx[2] * dv[1]);
              Jgas_inRad[1] += P[p].Mass * (dx[2] * dv[0] - dx[0] * dv[2]);
              Jgas_inRad[2] += P[p].Mass * (dx[0] * dv[1] - dx[1] * dv[0]);
            }
          if(ptype == 4)  // stars
            {
              Jstars_inRad[0] += P[p].Mass * (dx[1] * dv[2] - dx[2] * dv[1]);
              Jstars_inRad[1] += P[p].Mass * (dx[2] * dv[0] - dx[0] * dv[2]);
              Jstars_inRad[2] += P[p].Mass * (dx[0] * dv[1] - dx[1] * dv[0]);
            }
        }

      /* properties inside (stellar half mass radius) */
      if((massinhalfrad > 0) && (rr_tmp < 1.0 * halfmassradtype[4]))
        {
          for(j = 0; j < 3; j++)
            {
              ddxx  = GRAVITY_NEAREST_X(P[p].Pos[j] - pos[j]);
              dx[j] = All.cf_atime * ddxx;
              dv[j] = vel_to_phys * (P[p].Vel[j] - vinhalfrad[j]);
              dv[j] += H_of_a * dx[j];
            }

          Jtot_inHalfRad[0] += P[p].Mass * (dx[1] * dv[2] - dx[2] * dv[1]);
          Jtot_inHalfRad[1] += P[p].Mass * (dx[2] * dv[0] - dx[0] * dv[2]);
          Jtot_inHalfRad[2] += P[p].Mass * (dx[0] * dv[1] - dx[1] * dv[0]);

          if(ptype == 1)  // dm illustris
            {
              Jdm_inHalfRad[0] += P[p].Mass * (dx[1] * dv[2] - dx[2] * dv[1]);
              Jdm_inHalfRad[1] += P[p].Mass * (dx[2] * dv[0] - dx[0] * dv[2]);
              Jdm_inHalfRad[2] += P[p].Mass * (dx[0] * dv[1] - dx[1] * dv[0]);
            }
          if(ptype == 0)  // gas
            {
              Jgas_inHalfRad[0] += P[p].Mass * (dx[1] * dv[2] - dx[2] * dv[1]);
              Jgas_inHalfRad[1] += P[p].Mass * (dx[2] * dv[0] - dx[0] * dv[2]);
              Jgas_inHalfRad[2] += P[p].Mass * (dx[0] * dv[1] - dx[1] * dv[0]);
            }
          if(ptype == 4)  // stars
            {
              Jstars_inHalfRad[0] += P[p].Mass * (dx[1] * dv[2] - dx[2] * dv[1]);
              Jstars_inHalfRad[1] += P[p].Mass * (dx[2] * dv[0] - dx[0] * dv[2]);
              Jstars_inHalfRad[2] += P[p].Mass * (dx[0] * dv[1] - dx[1] * dv[0]);
            }
        }
    }

  if(parallel_flag)
    {
      MPI_Allreduce(MPI_IN_PLACE, Jtot_inRad, 3, MPI_DOUBLE, MPI_SUM, SubComm);
      MPI_Allreduce(MPI_IN_PLACE, Jdm_inRad, 3, MPI_DOUBLE, MPI_SUM, SubComm);
      MPI_Allreduce(MPI_IN_PLACE, Jgas_inRad, 3, MPI_DOUBLE, MPI_SUM, SubComm);
      MPI_Allreduce(MPI_IN_PLACE, Jstars_inRad, 3, MPI_DOUBLE, MPI_SUM, SubComm);
      MPI_Allreduce(MPI_IN_PLACE, Jtot_inHalfRad, 3, MPI_DOUBLE, MPI_SUM, SubComm);
      MPI_Allreduce(MPI_IN_PLACE, Jdm_inHalfRad, 3, MPI_DOUBLE, MPI_SUM, SubComm);
      MPI_Allreduce(MPI_IN_PLACE, Jgas_inHalfRad, 3, MPI_DOUBLE, MPI_SUM, SubComm);
      MPI_Allreduce(MPI_IN_PLACE, Jstars_inHalfRad, 3, MPI_DOUBLE, MPI_SUM, SubComm);
    }

  // counter rotating mass fractions
  CMFrac_inHalfRad = CMFrac_inRad = 0;
  for(i = 0; i < NTYPES; i++)
    CMFracType_inHalfRad[i] = CMFracType_inRad[i] = 0;

  for(i = 0; i < num; i++)
    {
      /* identify particle type */
      p = d[i].index;

      /* calculate particle radius */
      for(j = 0, rr_tmp = 0; j < 3; j++)
        {
          ddxx = GRAVITY_NEAREST_X(P[p].Pos[j] - pos[j]);  // counter-rotating mass calc with respect to pot min
          rr_tmp += ddxx * ddxx;
        }
      rr_tmp = sqrt(rr_tmp);

      int ptype = P[p].Type;

      /* properties inside SUBFIND_GAL_RADIUS_FAC*(stellar half mass radius) */
      if((massinrad > 0) && (rr_tmp < SUBFIND_GAL_RADIUS_FAC * halfmassradtype[4]))
        {
          for(j = 0; j < 3; j++)
            {
              ddxx  = GRAVITY_NEAREST_X(P[p].Pos[j] - pos[j]);
              dx[j] = All.cf_atime * ddxx;
              dv[j] = vel_to_phys * (P[p].Vel[j] - vinrad[j]);
              dv[j] += H_of_a * dx[j];
            }

          jpart[0] = P[p].Mass * (dx[1] * dv[2] - dx[2] * dv[1]);
          jpart[1] = P[p].Mass * (dx[2] * dv[0] - dx[0] * dv[2]);
          jpart[2] = P[p].Mass * (dx[0] * dv[1] - dx[1] * dv[0]);

          if((Jtot_inRad[0] * jpart[0] + Jtot_inRad[1] * jpart[1] + Jtot_inRad[2] * jpart[2]) < 0.)
            CMFrac_inRad += P[p].Mass / massinrad;

          if(ptype == 1)  // dm illustris
            if((Jdm_inRad[0] * jpart[0] + Jdm_inRad[1] * jpart[1] + Jdm_inRad[2] * jpart[2]) < 0.)
              CMFracType_inRad[1] += P[p].Mass / massinrad_tab[1];
          if(ptype == 0)  // gas (incl. winds!)
            if((Jgas_inRad[0] * jpart[0] + Jgas_inRad[1] * jpart[1] + Jgas_inRad[2] * jpart[2]) < 0.)
              CMFracType_inRad[0] += P[p].Mass / massinrad_tab[0];
          if(ptype == 4)  // stars
            if((Jstars_inRad[0] * jpart[0] + Jstars_inRad[1] * jpart[1] + Jstars_inRad[2] * jpart[2]) < 0.)
              CMFracType_inRad[4] += P[p].Mass / massinrad_tab[4];
        }

      /* properties inside (stellar half mass radius) */
      if((massinhalfrad > 0) && (rr_tmp < 1.0 * halfmassradtype[4]))
        {
          for(j = 0; j < 3; j++)
            {
              ddxx  = GRAVITY_NEAREST_X(P[p].Pos[j] - pos[j]);
              dx[j] = All.cf_atime * ddxx;
              dv[j] = vel_to_phys * (P[p].Vel[j] - vinhalfrad[j]);
              dv[j] += H_of_a * dx[j];
            }

          jpart[0] = P[p].Mass * (dx[1] * dv[2] - dx[2] * dv[1]);
          jpart[1] = P[p].Mass * (dx[2] * dv[0] - dx[0] * dv[2]);
          jpart[2] = P[p].Mass * (dx[0] * dv[1] - dx[1] * dv[0]);

          if((Jtot_inHalfRad[0] * jpart[0] + Jtot_inHalfRad[1] * jpart[1] + Jtot_inHalfRad[2] * jpart[2]) < 0.)
            CMFrac_inHalfRad += P[p].Mass / massinhalfrad;

          if(ptype == 1)  // dm illustris
            if((Jdm_inHalfRad[0] * jpart[0] + Jdm_inHalfRad[1] * jpart[1] + Jdm_inHalfRad[2] * jpart[2]) < 0.)
              CMFracType_inHalfRad[1] += P[p].Mass / massinhalfrad_tab[1];
          if(ptype == 0)  // gas (incl. winds!)
            if((Jgas_inHalfRad[0] * jpart[0] + Jgas_inHalfRad[1] * jpart[1] + Jgas_inHalfRad[2] * jpart[2]) < 0.)
              CMFracType_inHalfRad[0] += P[p].Mass / massinhalfrad_tab[0];
          if(ptype == 4)  // stars
            if((Jstars_inHalfRad[0] * jpart[0] + Jstars_inHalfRad[1] * jpart[1] + Jstars_inHalfRad[2] * jpart[2]) < 0.)
              CMFracType_inHalfRad[4] += P[p].Mass / massinhalfrad_tab[4];
        }
    }

  if(parallel_flag)
    {
      MPI_Allreduce(MPI_IN_PLACE, &CMFrac_inRad, 1, MPI_DOUBLE, MPI_SUM, SubComm);
      MPI_Allreduce(MPI_IN_PLACE, &CMFrac_inHalfRad, 1, MPI_DOUBLE, MPI_SUM, SubComm);
      MPI_Allreduce(MPI_IN_PLACE, CMFracType_inRad, NTYPES, MPI_DOUBLE, MPI_SUM, SubComm);
      MPI_Allreduce(MPI_IN_PLACE, CMFracType_inHalfRad, NTYPES, MPI_DOUBLE, MPI_SUM, SubComm);
    }
#endif /* #ifdef SUBFIND_EXTENDED_PROPERTIES */

  /* now store the calculated properties in the subgroup structure */
  if(parallel_flag == 0 || SubThisTask == 0)
    {
      subgroup->Len              = totlen;
      subgroup->Mass             = mass;
      subgroup->SubMassInRad     = massinrad;
      subgroup->SubMassInHalfRad = massinhalfrad;
      subgroup->SubMassInMaxRad  = massinmaxrad;
#ifdef SUBFIND_EXTENDED_PROPERTIES
      subgroup->Ekin             = Ekin;
      subgroup->Epot             = Epot;
      subgroup->Ethr             = Ethr;
      subgroup->CMFrac           = CMFrac;
      subgroup->CMFrac_inHalfRad = CMFrac_inHalfRad;
      subgroup->CMFrac_inRad     = CMFrac_inRad;
#endif /* #ifdef SUBFIND_EXTENDED_PROPERTIES */

#ifdef MHD
      subgroup->Bfld_Halo = bfld_halo;
      subgroup->Bfld_Disk = bfld_disk;
#endif /* #ifdef MHD */

      for(j = 0; j < 6; j++)
        {
          subgroup->MassType[j]             = mass_tab[j];
          subgroup->LenType[j]              = len_type[j];
          subgroup->SubHalfMassRadType[j]   = halfmassradtype[j];
          subgroup->SubMassInRadType[j]     = massinrad_tab[j];
          subgroup->SubMassInHalfRadType[j] = massinhalfrad_tab[j];
          subgroup->SubMassInMaxRadType[j]  = massinmaxrad_tab[j];
#ifdef SUBFIND_EXTENDED_PROPERTIES
          subgroup->CMFracType[j]           = CMFracType[j];
          subgroup->CMFracType_inHalfRad[j] = CMFracType_inHalfRad[j];
          subgroup->CMFracType_inRad[j]     = CMFracType_inRad[j];
#endif /* #ifdef SUBFIND_EXTENDED_PROPERTIES */
        }
      for(j = 0; j < 3; j++)
        {
          subgroup->Pos[j]  = pos[j];
          subgroup->Vel[j]  = vel[j];
          subgroup->CM[j]   = cm[j];
          subgroup->Spin[j] = spin[j];
#ifdef SUBFIND_EXTENDED_PROPERTIES
          subgroup->J[j]                = Jtot[j];
          subgroup->Jdm[j]              = Jdm[j];
          subgroup->Jgas[j]             = Jgas[j];
          subgroup->Jstars[j]           = Jstars[j];
          subgroup->J_inHalfRad[j]      = Jtot_inHalfRad[j];
          subgroup->Jdm_inHalfRad[j]    = Jdm_inHalfRad[j];
          subgroup->Jgas_inHalfRad[j]   = Jgas_inHalfRad[j];
          subgroup->Jstars_inHalfRad[j] = Jstars_inHalfRad[j];
          subgroup->J_inRad[j]          = Jtot_inRad[j];
          subgroup->Jdm_inRad[j]        = Jdm_inRad[j];
          subgroup->Jgas_inRad[j]       = Jgas_inRad[j];
          subgroup->Jstars_inRad[j]     = Jstars_inRad[j];
#endif /* #ifdef SUBFIND_EXTENDED_PROPERTIES */
        }

      subgroup->SubMostBoundID = mostboundid;
      subgroup->SubVelDisp     = veldisp;
      subgroup->SubVmax        = vmax;
      subgroup->SubVmaxRad     = vmaxrad;
      subgroup->SubHalfMassRad = halfmassrad;

#ifdef USE_SFR
      subgroup->Sfr          = sfr;
      subgroup->SfrInRad     = sfrinrad;
      subgroup->SfrInHalfRad = sfrinhalfrad;
      subgroup->SfrInMaxRad  = sfrinmaxrad;
      subgroup->GasMassSfr   = gasMassSfr;
#endif /* #ifdef USE_SFR */
    }
}

#endif /* #ifdef SUBFIND */
