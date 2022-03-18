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
 * \file        src/io/read_ic.c
 * \date        05/2018
 * \brief       Contains the routines needed to load initial conditions.
 * \details     contains functions:
 *                void read_ic(const char *fname, int readTypes)
 *                MyIDType determine_ids_offset(void)
 *                void empty_read_buffer(enum iofields blocknr, int offset,
 *                  int pc, int type)
 *                void share_particle_number_in_file(const char *fname, int
 *                  filenr, int readTask, int lastTask, int readTypes)
 *                void read_file(const char *fname, int filenr, int readTask,
 *                  int lastTask, int readTypes)
 *                int find_files(const char *fname)
 *                void distribute_file(int nfiles, int firstfile, int
 *                  firsttask, int lasttask, int *filenr, int *master, int
 *                  *last)
 *                herr_t hdf5_header_error_handler(void *unused)
 *                void read_header_attributes_in_hdf5(const char *fname)
 *                void read_header_attributes(FILE * fd)
 *                void swap_Nbyte(char *data, int n, int m)
 *                void swap_header()
 *                void tile_ics(void)
 *
 * \par Major modifications and contributions:
 *
 * - DD.MM.YYYY Description
 * - 08.05.2018 Prepared file for public release -- Rainer Weinberger
 */

#include <math.h>
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "../main/allvars.h"
#include "../main/proto.h"

#ifndef IDS_OFFSET
#ifdef LONGIDS
#define IDS_OFFSET 100000000000
#else /* #ifdef LONGIDS */
#define IDS_OFFSET 1000000000
#endif /* #ifdef LONGIDS #else */
#endif /* #ifndef IDS_OFFSET */

#define SKIP                                 \
  {                                          \
    my_fread(&blksize1, sizeof(int), 1, fd); \
  }
#define SKIP2                                \
  {                                          \
    my_fread(&blksize2, sizeof(int), 1, fd); \
  }

void read_header_attributes(FILE *fd);

#ifdef HAVE_HDF5
#include <hdf5.h>
void read_header_attributes_in_hdf5(const char *fname);
#endif /* #ifdef HAVE_HDF5 */

int num_files;

int swap_file = 8;

#if defined(RECOMPUTE_POTENTIAL_IN_SNAPSHOT)
/*! \brief Struct containing information about the number of particles per
 *         particle type.
 */
static struct ntypes_data
{
  int npart[NTYPES];
} * ntype_in_files;
#endif /* #if defined(RECOMPUTE_POTENTIAL_IN_SNAPSHOT) */

/*! \brief Reads initial conditions that are in one of the supported file
 *         formats.
 *
 *  Snapshot files can be used as input files. However, when a
 *  snapshot file is used as input, not all the information in the header is
 *  used: THE STARTING TIME NEEDS TO BE SET IN THE PARAMETERFILE.
 *  Alternatively, the code can be started with restartflag 2, then snapshots
 *  from the code can be used as initial conditions-files without having to
 *  change the parameter file. For gas particles, only the internal energy is
 *  read, the density and mean molecular weight will be recomputed by the code.
 *  When InitGasTemp>0 is given, the gas temperature will be initialized to
 *  this value assuming a mean molecular weight either corresponding to
 *  complete neutrality, or full ionization.
 *
 *  \param[in] fname File name of the ICs.
 *  \param[in] readTypes A bitfield that determines what particle types to
 *             read, only if the bit corresponding to a particle type is set,
 *             the corresponding data is loaded, otherwise its particle number
 *             is set to zero. (This is only implemented for HDF5 files.)
 *
 *  \return void
 */
void read_ic(const char *fname, int readTypes)
{
  int i, rep, rest_files, ngroups, gr, filenr, masterTask, lastTask, groupMaster;
  double u_init, molecular_weight;
  char buf[500];
  double t0, t1;

  if((All.ICFormat < 1) || (All.ICFormat > 4))
    {
      mpi_terminate("ICFormat=%d not supported.\n", All.ICFormat);
    }

  t0 = second();
  CPU_Step[CPU_MISC] += measure_time();

  num_files = find_files(fname);

#if defined(RECOMPUTE_POTENTIAL_IN_SNAPSHOT)
  ntype_in_files = mymalloc("ntype_in_files", num_files * sizeof(struct ntypes_data));
  memset(ntype_in_files, 0, num_files * sizeof(struct ntypes_data));
#endif /* #if defined(RECOMPUTE_POTENTIAL_IN_SNAPSHOT) */

  All.TotNumPart = 0;

  /* we repeat reading the headers of the files two times. In the first iteration, only the
   * particle numbers ending up on each processor are assembled, followed by memory allocation.
   * In the second iteration, the data is actually read in.
   */
  for(rep = 0; rep < 2; rep++)
    {
      NumPart = 0;
      NumGas  = 0;

#if defined(RECOMPUTE_POTENTIAL_IN_SNAPSHOT)
      if(rep == 1)
        MPI_Allreduce(MPI_IN_PLACE, ntype_in_files, num_files * NTYPES, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
#endif /* #if defined(RECOMPUTE_POTENTIAL_IN_SNAPSHOT) */

      rest_files = num_files;
      while(rest_files > NTask)
        {
          sprintf(buf, "%s.%d", fname, ThisTask + (rest_files - NTask));
          if(All.ICFormat == 3)
            sprintf(buf, "%s.%d.hdf5", fname, ThisTask + (rest_files - NTask));

          ngroups = NTask / All.NumFilesWrittenInParallel;
          if((NTask % All.NumFilesWrittenInParallel))
            ngroups++;
          groupMaster = (ThisTask / ngroups) * ngroups;

          for(gr = 0; gr < ngroups; gr++)
            {
              if(ThisTask == (groupMaster + gr)) /* ok, it's this processor's turn */
                {
                  if(rep == 0)
                    share_particle_number_in_file(buf, ThisTask + (rest_files - NTask), ThisTask, ThisTask, readTypes);
                  else
                    read_file(buf, ThisTask + (rest_files - NTask), ThisTask, ThisTask, readTypes);
                }
              MPI_Barrier(MPI_COMM_WORLD);
            }

          rest_files -= NTask;
        }

      if(rest_files > 0)
        {
          distribute_file(rest_files, 0, 0, NTask - 1, &filenr, &masterTask, &lastTask);

          if(num_files > 1)
            {
              sprintf(buf, "%s.%d", fname, filenr);
              if(All.ICFormat == 3)
                sprintf(buf, "%s.%d.hdf5", fname, filenr);
            }
          else
            {
              sprintf(buf, "%s", fname);
              if(All.ICFormat == 3)
                sprintf(buf, "%s.hdf5", fname);
            }

          ngroups = rest_files / All.NumFilesWrittenInParallel;
          if((rest_files % All.NumFilesWrittenInParallel))
            ngroups++;

          for(gr = 0; gr < ngroups; gr++)
            {
              if((filenr / All.NumFilesWrittenInParallel) == gr) /* ok, it's this processor's turn */
                {
                  if(rep == 0)
                    share_particle_number_in_file(buf, filenr, masterTask, lastTask, readTypes);
                  else
                    read_file(buf, filenr, masterTask, lastTask, readTypes);
                }
              MPI_Barrier(MPI_COMM_WORLD);
            }
        }

      /* now do the memory allocation */
      if(rep == 0)
        {
          int max_load, max_sphload;
          MPI_Allreduce(&NumPart, &max_load, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
          MPI_Allreduce(&NumGas, &max_sphload, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);

#ifdef GENERATE_GAS_IN_ICS
          if(max_sphload < max_load)
            max_sphload = max_load;
#endif /* #ifdef GENERATE_GAS_IN_ICS */

          All.MaxPart    = max_load / (1.0 - 2 * ALLOC_TOLERANCE);
          All.MaxPartSph = max_sphload / (1.0 - 2 * ALLOC_TOLERANCE);

#ifdef EXACT_GRAVITY_FOR_PARTICLE_TYPE
          if(All.TotPartSpecial != 0)
            All.MaxPartSpecial = (int)(All.TotPartSpecial);
          else
            terminate("Code compiled with option EXACT_GRAVITY_FOR_PARTICLE_TYPE but no particles of specified type found in ICs.");
#endif /* #ifdef EXACT_GRAVITY_FOR_PARTICLE_TYPE */
          allocate_memory();

          CommBuffer = mymalloc("CommBuffer", COMMBUFFERSIZE);
        }
    }

  myfree(CommBuffer);

#ifdef TILE_ICS
  tile_ics();
#endif /* #ifdef TILE_ICS */

  /* this makes sure that masses are initialized in the case that the mass-block
     is empty for this particle type */
  for(i = 0; i < NumPart; i++)
    {
      if(All.MassTable[P[i].Type] != 0)
        P[i].Mass = All.MassTable[P[i].Type];
    }

    /* If we are reading in Gadget2 ICs, we need to compute the material
       number from the ID  */
#ifdef READ_LEGACY_ICS
  if(header.flag_entropy_instead_u)
    {
      sprintf(buf, "\nProblem: Legacy ICs cannot contain entropy in the u field!\n");
      terminate(buf);
    }

  for(i = 0; i < NumGas; i++)
    {
      int j;

      double mat;

      modf(((double)(P[i].ID - EOS_ID_START)) / EOS_ID_SKIP, &mat); /* This stores the int part in variable mat and
                                                                       discards the remainder */
      int imat = mat;

      SphP[i].Composition[imat] = 1.0;
    }
#endif /* #ifdef READ_LEGACY_ICS */

#if defined(REFINEMENT) && defined(REFINEMENT_HIGH_RES_GAS)
  if(RestartFlag == 0) /* All gas that is already present in the ICs is allowed to be (de-)refined */
    {
      for(i = 0; i < NumGas; i++)
        {
          if(All.ReferenceGasPartMass == 0 || P[i].Mass < 1.2 * All.ReferenceGasPartMass)
            SphP[i].AllowRefinement = 1;
        }
    }
#endif /* #if defined (REFINEMENT) && defined (REFINEMENT_HIGH_RES_GAS) */

  for(i = 0; i < NumPart; i++)
    P[i].SofteningType = All.SofteningTypeOfPartType[P[i].Type];

#ifdef GENERATE_GAS_IN_ICS
  int count;
  double fac, d, a, b, rho;

  if(RestartFlag == 0)
    {
      header.flag_entropy_instead_u = 0;

      MyIDType ids_offset = determine_ids_offset();

      for(i = 0, count = 0; i < NumPart; i++)
#ifdef SPLIT_PARTICLE_TYPE
        if((1 << P[i].Type) & (SPLIT_PARTICLE_TYPE))
#else  /* #ifdef SPLIT_PARTICLE_TYPE */
        if(P[i].Type == 1)
#endif /* #ifdef SPLIT_PARTICLE_TYPE #else */
          count++;

      if(count)
        {
          domain_resize_storage(count, count, 0);

          memmove(P + count, P, sizeof(struct particle_data) * NumPart);

          NumPart += count;
          NumGas += count;

          if(NumGas > All.MaxPartSph)
            terminate("Task=%d ends up getting more SPH particles (%d) than allowed (%d)\n", ThisTask, NumGas, All.MaxPartSph);

#ifdef REFINEMENT_HIGH_RES_GAS
          for(i = 0; i < NumGas - count; i++) /* make sure that AllowRefinement is shifted with the particles */
            SphP[i + count].AllowRefinement = SphP[i].AllowRefinement;
          for(i = 0; i < count; i++) /* by default, new cells are not allowed to be refined */
            SphP[i].AllowRefinement = 0;
#endif /* #ifdef REFINEMENT_HIGH_RES_GAS */

          fac = All.OmegaBaryon / All.Omega0;
          rho = All.Omega0 * 3 * All.Hubble * All.Hubble / (8 * M_PI * All.G);

          int j;

          for(i = count, j = 0; i < NumPart; i++)
#ifdef SPLIT_PARTICLE_TYPE
            if((1 << P[i].Type) & (SPLIT_PARTICLE_TYPE))
#else  /* #ifdef SPLIT_PARTICLE_TYPE */
            if(P[i].Type == 1)
#endif /* #ifdef SPLIT_PARTICLE_TYPE #else */
              {
                d = pow(P[i].Mass / rho, 1.0 / 3);
                a = 0.5 * All.OmegaBaryon / All.Omega0 * d;
                b = 0.5 * (All.Omega0 - All.OmegaBaryon) / All.Omega0 * d;

                P[j] = P[i];

                P[j].Mass *= fac;
                P[i].Mass *= (1 - fac);
                P[j].Type = 0;
                P[j].ID += ids_offset;
                P[i].Pos[0] += a;
                P[i].Pos[1] += a;
                P[i].Pos[2] += a;
                P[j].Pos[0] -= b;
                P[j].Pos[1] -= b;
                P[j].Pos[2] -= b;

#ifdef REFINEMENT_HIGH_RES_GAS
                if(P[i].Type == 1) /* also allow gas which is produced by splitting a high res DM particle to be (de-) refined */
                  SphP[j].AllowRefinement = 2;
#endif /* #ifdef REFINEMENT_HIGH_RES_GAS */

                j++;
              }

          All.MassTable[0] = 0;

#ifdef SPLIT_PARTICLE_TYPE
          for(i = 1; i < NTYPES; i++)
            if((1 << i) & (SPLIT_PARTICLE_TYPE))
              All.MassTable[i] *= (1 - fac);
#else  /* #ifdef SPLIT_PARTICLE_TYPE */
          All.MassTable[1] *= (1 - fac);
#endif /* #ifdef SPLIT_PARTICLE_TYPE #else */
        }
    }
#endif /* #ifdef GENERATE_GAS_IN_ICS */

#ifdef READ_DM_AS_GAS
  {
    domain_resize_storage(0, NumPart, 0);

    if(NumGas > All.MaxPartSph)
      terminate("Task=%d ends up getting more SPH particles (%d) than allowed (%d)\n", ThisTask, NumGas, All.MaxPartSph);

    for(i = 0; i < NumPart; i++)
      {
        P[i].Type      = 0;
        SphP[i].Utherm = 1.0;
      }

    All.MassTable[0] = 0;

    header.npartTotal[0]         = header.npartTotal[1];
    header.npartTotalHighWord[0] = header.npartTotalHighWord[1];
    header.npart[0]              = header.npart[1];
    header.npartTotal[1]         = 0;
    header.npartTotalHighWord[1] = 0;
    header.npart[1]              = 0;
    NumGas                       = NumPart;
    All.TotNumGas                = All.TotNumPart;
    mpi_printf("READ_DM_AS_GAS: generated %lld gas particles from type %d\n",
               header.npartTotal[0] + (((long long)header.npartTotalHighWord[0]) << 32), 0);
  }
#endif /* #ifdef READ_DM_AS_GAS */

#ifdef USE_SFR
  if(RestartFlag == 0)
    {
      if(All.MassTable[4] == 0 && All.MassTable[0] > 0)
        {
          All.MassTable[0] = 0;
          All.MassTable[4] = 0;
        }
    }
#endif

  u_init = (1.0 / GAMMA_MINUS1) * (BOLTZMANN / PROTONMASS) * All.InitGasTemp;
  u_init *= All.UnitMass_in_g / All.UnitEnergy_in_cgs; /* unit conversion */

  if(All.InitGasTemp > 1.0e4) /* assuming FULL ionization */
    molecular_weight = 4 / (8 - 5 * (1 - HYDROGEN_MASSFRAC));
  else /* assuming NEUTRAL GAS */
    molecular_weight = 4 / (1 + 3 * HYDROGEN_MASSFRAC);

  u_init /= molecular_weight;

  All.InitGasU = u_init;

  header.mass[0]   = 0; /* to make sure that the variable masses are stored in output file */
  All.MassTable[0] = 0;

  if(RestartFlag == 0)
    {
#if defined(REFINEMENT_HIGH_RES_GAS)
      for(i = 0; i < NumGas; i++)
        if(SphP[i].AllowRefinement)
          SphP[i].HighResMass = P[i].Mass;
        else
          SphP[i].HighResMass = 0;
#endif /* #if defined(REFINEMENT_HIGH_RES_GAS) */

      if(All.InitGasTemp > 0)
        {
          for(i = 0; i < NumGas; i++)
            {
              if(ThisTask == 0 && i == 0 && SphP[i].Utherm == 0)
                printf("READIC: Initializing u from InitGasTemp!\n");

              if(SphP[i].Utherm == 0)
                SphP[i].Utherm = All.InitGasU;
              /* Note: the coversion to entropy will be done in the function init(),
                 after the densities have been computed */
            }
        }
    }

  for(i = 0; i < NumGas; i++)
    {
      SphP[i].Utherm = dmax(All.MinEgySpec, SphP[i].Utherm);
      if(SphP[i].Density > 0)
        SphP[i].Volume = P[i].Mass / SphP[i].Density;
    }

  MPI_Barrier(MPI_COMM_WORLD);

  t1 = second();
  mpi_printf("READIC: reading done (took %g sec).\n", timediff(t0, t1));

  /* verify number of particles */
  int num = 0;
  long long glob_num;
  for(i = 0; i < NumPart; i++)
    num += 1;
  sumup_large_ints(1, &num, &glob_num);
  if(glob_num != All.TotNumPart)
    terminate("glob_num (=%lld) != All.TotNumPart (=%lld)", glob_num, All.TotNumPart);

  mpi_printf("READIC: Total number of particles :  %lld\n\n", All.TotNumPart);

  CPU_Step[CPU_SNAPSHOT] += measure_time();
}

/*! \brief This function computes a suitable offset for the particle IDs in
 *         case gas should be generated in the ICs.
 *
 *  If the macro OFFSET_FOR_NON_CONTIGUOUS_IDS is not defined the code reverts
 *  to a fixed offset defined at the beginning of the file.
 *
 *  \return Offset for the gas particles to be generated.
 */
MyIDType determine_ids_offset(void)
{
#ifndef OFFSET_FOR_NON_CONTIGUOUS_IDS
  MyIDType ids_offset = IDS_OFFSET;
#else /* #ifndef OFFSET_FOR_NON_CONTIGUOUS_IDS */
  if(All.MaxID == 0) /* MaxID not calculated yet */
    calculate_maxid();

  int bits_used       = 1;
  int bits_available  = CHAR_BIT * sizeof(MyIDType);
  MyIDType ids_offset = 1;

  while(ids_offset <= All.MaxID && ids_offset > 0)
    {
      ids_offset <<= 1;
      bits_used++;
    }

  All.MaxID = 0; /* reset to allow recomputing */

  if(ids_offset <= 0)
    terminate("not enough memory to generate id offsets. Used %d bits out of %d\n", bits_used, bits_available);

#ifdef LONGIDS
  mpi_printf("GENERATE_GAS_IN_ICS: determined id offset as %llu. Used %d bits out of %d\n", ids_offset, bits_used, bits_available);
#else  /* #ifdef LONGIDS */
  mpi_printf("GENERATE_GAS_IN_ICS: determined id offset as %u. Used %d bits out of %d\n", ids_offset, bits_used, bits_available);
#endif /* #ifdef LONGIDS #else */

#endif /* #ifndef OFFSET_FOR_NON_CONTIGUOUS_IDS */
  return ids_offset;
}

/*! \brief Reads out the io buffer that was filled with particle data.
 *
 *  The data in the io buffer is put in the appropriate places of the particle
 *  structures.
 *
 * \param[in] blocknr Data block present in io buffer.
 * \param[in] offset Particle corresponding to the first element in io buffer.
 * \param[in] pc Number of elements in the io buffer.
 * \param[in] type If blocknr=IO_POS P[n].Type is set to type.
 *
 * \return void
 */
void empty_read_buffer(enum iofields blocknr, int offset, int pc, int type)
{
  int n, k;
  MyInputFloat *fp;
  double *doublep;
  MyIDType *ip;
  int *intp;
  float *floatp;

  int vt, vpb;
  char *cp;

  fp      = (MyInputFloat *)CommBuffer;
  doublep = (double *)CommBuffer;
  ip      = (MyIDType *)CommBuffer;
  intp    = (int *)CommBuffer;
  floatp  = (float *)CommBuffer;

  cp  = (char *)CommBuffer;
  vt  = get_datatype_in_block(blocknr, 1);
  vpb = get_values_per_blockelement(blocknr);
  if(vt == 2)
    swap_Nbyte(cp, pc * vpb, 8);
  else
    {
#ifdef INPUT_IN_DOUBLEPRECISION
      if(vt == 1)
        swap_Nbyte(cp, pc * vpb, 8);
      else
#endif /* #ifdef INPUT_IN_DOUBLEPRECISION */
        swap_Nbyte(cp, pc * vpb, 4);
    }

  int field = -1;
  int f;
  for(f = 0; f < N_IO_Fields; f++)
    {
      if(IO_Fields[f].field == blocknr)
        {
          field = f;
          break;
        }
    }

  if(field < 0)
    terminate("error: field not found");

  for(n = 0; n < pc; n++)
    {
      if(IO_Fields[field].io_func)
        {
          int particle;
          switch(IO_Fields[field].array)
            {
              case A_NONE:
              case A_SPHP:
              case A_P:
                particle = offset + n;
                break;
              case A_PS:
                terminate("Not good, trying to read into PS[]?\n");
                break;
              default:
                terminate("ERROR in empty_read_buffer: Array not found!\n");
                break;
            }

          switch(IO_Fields[field].type_in_file_input)
            {
              case FILE_NONE:
                terminate("error");
                break;
              case FILE_INT:
                IO_Fields[field].io_func(particle, IO_Fields[field].values_per_block, intp, 1);
                intp += IO_Fields[field].values_per_block;
                break;
              case FILE_MY_ID_TYPE:
                IO_Fields[field].io_func(particle, IO_Fields[field].values_per_block, ip, 1);
                ip += IO_Fields[field].values_per_block;
                break;
              case FILE_MY_IO_FLOAT:
                IO_Fields[field].io_func(particle, IO_Fields[field].values_per_block, fp, 1);
                fp += IO_Fields[field].values_per_block;
                break;
              case FILE_DOUBLE:
                IO_Fields[field].io_func(particle, IO_Fields[field].values_per_block, doublep, 1);
                doublep += IO_Fields[field].values_per_block;
                break;
              case FILE_FLOAT:
                IO_Fields[field].io_func(particle, IO_Fields[field].values_per_block, floatp, 1);
                floatp += IO_Fields[field].values_per_block;
                break;
            }
        }
      else
        {
          void *array_pos;
          switch(IO_Fields[field].array)
            {
              case A_NONE:
                array_pos = 0;
                break;
              case A_SPHP:
                array_pos = SphP + offset + n;
                break;
              case A_P:
                array_pos = P + offset + n;
                break;
              case A_PS:
                terminate("Not good, trying to read into PS[]?\n");
                break;
              default:
                terminate("ERROR in empty_read_buffer: Array not found!\n");
                break;
            }

          for(k = 0; k < IO_Fields[field].values_per_block; k++)
            {
              double value = 0;
              switch(IO_Fields[field].type_in_file_input)
                {
                  case FILE_MY_IO_FLOAT:
                    value = *fp;
                    fp++;
                    break;
                  case FILE_DOUBLE:
                    value = *doublep;
                    doublep++;
                    break;
                  case FILE_FLOAT:
                    value = *floatp;
                    floatp++;
                    break;
                  default:
                    break;
                }

              switch(IO_Fields[field].type_in_memory)
                {
                  case MEM_INT:
                    *((int *)((size_t)array_pos + IO_Fields[field].offset + k * sizeof(int))) = *intp;
                    intp++;
                    break;
                  case MEM_MY_ID_TYPE:
                    *((MyIDType *)((size_t)array_pos + IO_Fields[field].offset + k * sizeof(MyIDType))) = *ip;
                    ip++;
                    break;
                  case MEM_FLOAT:
                    *((float *)((size_t)array_pos + IO_Fields[field].offset + k * sizeof(float))) = value;
                    break;

                  case MEM_DOUBLE:
                    *((double *)((size_t)array_pos + IO_Fields[field].offset + k * sizeof(double))) = value;
                    break;

                  case MEM_MY_SINGLE:
                    *((MySingle *)((size_t)array_pos + IO_Fields[field].offset + k * sizeof(MySingle))) = value;
                    break;

                  case MEM_MY_FLOAT:
                    *((MyFloat *)((size_t)array_pos + IO_Fields[field].offset + k * sizeof(MyFloat))) = value;
                    break;

                  case MEM_MY_DOUBLE:
                    *((MyDouble *)((size_t)array_pos + IO_Fields[field].offset + k * sizeof(MyDouble))) = value;
                    break;

                  default:
                    terminate("ERROR in empty_read_buffer: Type not found!\n");
                    break;
                }
            }
        }
    }

  if(blocknr == IO_VEL)
    {
      for(n = 0; n < pc; n++)
        P[offset + n].Type = type; /* initialize type here as well */
    }
}

/*! \brief Distributes the particle numbers in the file fname
 *         to tasks 'readTask' to 'lastTask', and calculates the number of
 *         particles each task gets.
 *
 *  \param[in] fname Filename to be read.
 *  \param[in] readTask Task responsible for reading the file fname.
 *  \param[in] lastTask Last task which gets data contained in the file.
 *  \param[in] readTypes A bitfield that determines what particle  types to
 *             read, only if the bit corresponding to a particle type is set,
 *             the corresponding data is loaded, otherwise its particle number
 *             is set to zero. (This is only implemented for HDF5 files.)
 *
 *  \return void
 */
void share_particle_number_in_file(const char *fname, int filenr, int readTask, int lastTask, int readTypes)
{
  int i, n_in_file, n_for_this_task, ntask, task;
  int blksize1, blksize2;
  MPI_Status status;
  FILE *fd = 0;
  int type;
  char label[4], buf[500];
  int nextblock;
#ifdef HAVE_HDF5
  hid_t hdf5_file = 0, hdf5_grp[NTYPES];
#endif /* #ifdef HAVE_HDF5 */

  if(ThisTask == readTask)
    {
      if(All.ICFormat == 1 || All.ICFormat == 2)
        {
          if(!(fd = fopen(fname, "r")))
            {
              sprintf(buf, "can't open file `%s' for reading initial conditions.\n", fname);
              terminate(buf);
            }

          if(All.ICFormat == 2)
            {
              SKIP;
              swap_file = blksize1;
              my_fread(&label, sizeof(char), 4, fd);
              my_fread(&nextblock, sizeof(int), 1, fd);
              swap_Nbyte((char *)&nextblock, 1, 4);
              printf("Reading header => '%c%c%c%c' (%d byte)\n", label[0], label[1], label[2], label[3], nextblock);
              SKIP2;
            }

          SKIP;
          if(All.ICFormat == 1)
            {
              if(blksize1 != 256)
                swap_file = 1;
            }
          read_header_attributes(fd);
          SKIP2;
          swap_Nbyte((char *)&blksize1, 1, 4);
          swap_Nbyte((char *)&blksize2, 1, 4);

          if(blksize1 != 256 || blksize2 != 256)
            terminate("incorrect header format blocksize %d, %d\n", blksize1, blksize2);

          swap_header();

#ifdef COMBINETYPES
          header.npartTotal[3] += header.npartTotal[4] + header.npartTotal[5];
          header.npart[3] += header.npart[4] + header.npart[5];
          header.npartTotal[4] = 0;
          header.npartTotal[5] = 0;
          header.npart[4]      = 0;
          header.npart[5]      = 0;
#endif /* #ifdef COMBINETYPES */
        }

#ifdef HAVE_HDF5
      if(All.ICFormat == 3)
        {
          read_header_attributes_in_hdf5(fname);

          hdf5_file = my_H5Fopen(fname, H5F_ACC_RDONLY, H5P_DEFAULT);
          if(hdf5_file < 0)
            terminate("cannot read initial conditions file %s", fname);

          for(type = 0; type < NTYPES; type++)
            {
              if(header.npart[type] > 0 && (readTypes & (1 << type)))
                {
                  sprintf(buf, "/PartType%d", type);
                  hdf5_grp[type] = my_H5Gopen(hdf5_file, buf);
                }
              if(!(readTypes & (1 << type)))
                {
                  // Override particle number in file. If we don't
                  // read the type, both npart and npartTotal will be 0
                  header.npartTotal[type]         = 0;
                  header.npart[type]              = 0;
                  header.npartTotalHighWord[type] = 0;
                  header.mass[type]               = 0;
                }
            }
        }
#endif /* #ifdef HAVE_HDF5 */

      for(task = readTask + 1; task <= lastTask; task++)
        {
          MPI_Ssend(&header, sizeof(header), MPI_BYTE, task, TAG_HEADER, MPI_COMM_WORLD);
          MPI_Ssend(&swap_file, sizeof(swap_file), MPI_BYTE, task, TAG_KEY, MPI_COMM_WORLD);
        }
    }
  else
    {
      MPI_Recv(&header, sizeof(header), MPI_BYTE, readTask, TAG_HEADER, MPI_COMM_WORLD, &status);
      MPI_Recv(&swap_file, sizeof(swap_file), MPI_BYTE, readTask, TAG_KEY, MPI_COMM_WORLD, &status);
    }

  if(header.num_files != num_files)
    warn("header.num_files=%d != num_files=%d", header.num_files, num_files);

  if(All.TotNumPart == 0)
    {
      if(num_files == 1)
        for(type = 0; type < NTYPES; type++)
          {
            if(header.npartTotal[type] != header.npart[type])
              {
                warn("header.npartTotal[%d]=%d != header.npart[%d]=%d, setting header.npartTotal[%d] = header.npart[%d]\n", type,
                     header.npartTotal[type], type, header.npart[type], type, type);
                header.npartTotal[type] = header.npart[type];
              }
#ifdef USE_SFR
            header.npartTotalHighWord[type] = 0;
#endif
          }

      All.TotNumGas = header.npartTotal[0] + (((long long)header.npartTotalHighWord[0]) << 32);
#ifdef EXACT_GRAVITY_FOR_PARTICLE_TYPE
      All.TotPartSpecial = header.npartTotal[EXACT_GRAVITY_FOR_PARTICLE_TYPE] +
                           (((long long)header.npartTotalHighWord[EXACT_GRAVITY_FOR_PARTICLE_TYPE]) << 32);
      mpi_printf("Tot Special %d %d %d %d\n", All.TotPartSpecial, EXACT_GRAVITY_FOR_PARTICLE_TYPE, header.npart[4],
                 header.npartTotal[4]);
#endif /* #ifdef EXACT_GRAVITY_FOR_PARTICLE_TYPE */

      for(type = 0, All.TotNumPart = 0; type < NTYPES; type++)
        {
          All.TotNumPart += header.npartTotal[type];
          All.TotNumPart += (((long long)header.npartTotalHighWord[type]) << 32);
        }

#ifdef GENERATE_GAS_IN_ICS
      if(RestartFlag == 0)
        {
          if(All.TotNumGas > 0)
            terminate("You specified GENERATE_GAS_IN_ICS but your ICs already contain gas! (namely %lld gas cells)\n", All.TotNumGas);

#ifdef SPLIT_PARTICLE_TYPE
          for(i = 0; i < NTYPES; i++)
            if((1 << i) & (SPLIT_PARTICLE_TYPE))
              {
                All.TotNumGas += header.npartTotal[i] + (((long long)header.npartTotalHighWord[i]) << 32);
                All.TotNumPart += header.npartTotal[i] + (((long long)header.npartTotalHighWord[i]) << 32);
                mpi_printf("GENERATE_GAS_IN_ICS: generated %lld gas particles from type %d\n",
                           header.npartTotal[i] + (((long long)header.npartTotalHighWord[i]) << 32), i);
              }
#else  /* #ifdef SPLIT_PARTICLE_TYPE */
          All.TotNumGas += header.npartTotal[1] + (((long long)header.npartTotalHighWord[1]) << 32);
          All.TotNumPart += header.npartTotal[1] + (((long long)header.npartTotalHighWord[1]) << 32);
          mpi_printf("GENERATE_GAS_IN_ICS: generated %lld gas particles from type 1\n",
                     header.npartTotal[1] + (((long long)header.npartTotalHighWord[1]) << 32));
#endif /* #ifdef SPLIT_PARTICLE_TYPE #else */
        }
#endif /* #ifdef GENERATE_GAS_IN_ICS */

#ifdef TILE_ICS
      All.TotNumPart *= All.TileICsFactor * All.TileICsFactor * All.TileICsFactor;
      All.TotNumGas *= All.TileICsFactor * All.TileICsFactor * All.TileICsFactor;
#ifdef EXACT_GRAVITY_FOR_PARTICLE_TYPE
      All.TotPartSpecial *= All.TileICsFactor * All.TileICsFactor * All.TileICsFactor;
#endif /* #ifdef EXACT_GRAVITY_FOR_PARTICLE_TYPE */
#endif /* #ifdef TILE_ICS */

      for(i = 0; i < NTYPES; i++)
        All.MassTable[i] = header.mass[i];

      if(RestartFlag >= 2)
        All.Time = All.TimeBegin = header.time;
      else
        All.Time = All.TimeBegin;

      set_cosmo_factors_for_current_time();
    }

  if(ThisTask == readTask)
    {
      for(type = 0, n_in_file = 0; type < NTYPES; type++)
        n_in_file += header.npart[type];

      printf("READIC: Reading file `%s' on task=%d and distribute it to %d to %d (contains %d particles).\n", fname, ThisTask,
             readTask, lastTask, n_in_file);

      myflush(stdout);
    }

  for(type = 0; type < NTYPES; type++)
    {
      n_in_file       = header.npart[type];
      ntask           = lastTask - readTask + 1;
      n_for_this_task = n_in_file / ntask;
      if((ThisTask - readTask) < (n_in_file % ntask))
        n_for_this_task++;

      NumPart += n_for_this_task;

      if(type == 0)
        NumGas += n_for_this_task;
    }

  if(ThisTask == readTask)
    {
      if(All.ICFormat == 1 || All.ICFormat == 2)
        fclose(fd);
#ifdef HAVE_HDF5
      if(All.ICFormat == 3)
        {
          for(type = NTYPES - 1; type >= 0; type--)
            if(header.npart[type] > 0)
              {
                sprintf(buf, "/PartType%d", type);
                my_H5Gclose(hdf5_grp[type], buf);
              }
          my_H5Fclose(hdf5_file, fname);
        }
#endif /* #ifdef HAVE_HDF5 */

#if defined(RECOMPUTE_POTENTIAL_IN_SNAPSHOT)
      for(int type = 0; type < NTYPES; type++)
        ntype_in_files[filenr].npart[type] = header.npart[type];
#endif /* #if defined(RECOMPUTE_POTENTIAL_IN_SNAPSHOT) */
    }
}

/*! \brief Reads a single snapshot file.
 *
 *  This routine reads a single file. The data it contains is
 *  distributed to tasks 'readTask' to 'lastTask'.
 *
 *  \param[in] fname Filename to be read.
 *  \param[in] readTask Task responsible for reading the file fname
 *  \param[in] lastTask Last task which gets data contained in the file
 *  \param[in] readTypes readTypes is a bitfield that determines what particle
 *             types to read, only if the bit corresponding to a particle type
 *             is set, the corresponding data is loaded, otherwise its particle
 *             number is set to zero. (This is only implemented for HDF5
 *             files.)
 *
 *  \return void
 */
void read_file(const char *fname, int filenr, int readTask, int lastTask, int readTypes)
{
  int blockmaxlen;
  int n_in_file, n_for_this_task, ntask, pc, offset = 0, task;
  int blksize1, blksize2;
  MPI_Status status;
  FILE *fd = 0;
  int nall;
  int type, bnr;
  char label[4], expected_label[4], buf[500];
  int nstart, bytes_per_blockelement, npart, nextblock, typelist[NTYPES];
  enum iofields blocknr;

#ifdef HAVE_HDF5
  int rank, pcsum;
  hid_t hdf5_file     = 0, hdf5_grp[NTYPES], hdf5_dataspace_in_file;
  hid_t hdf5_datatype = 0, hdf5_dataspace_in_memory, hdf5_dataset;
  hsize_t dims[2], count[2], start[2];
#endif /* #ifdef HAVE_HDF5 */

  if(ThisTask == readTask)
    {
      if(All.ICFormat == 1 || All.ICFormat == 2)
        {
          if(!(fd = fopen(fname, "r")))
            {
              sprintf(buf, "can't open file `%s' for reading initial conditions.\n", fname);
              terminate(buf);
            }

          if(All.ICFormat == 2)
            {
              SKIP;
              swap_file = blksize1;
              my_fread(&label, sizeof(char), 4, fd);
              my_fread(&nextblock, sizeof(int), 1, fd);
              swap_Nbyte((char *)&nextblock, 1, 4);
              SKIP2;
            }

          SKIP;
          if(All.ICFormat == 1)
            {
              if(blksize1 != 256)
                swap_file = 1;
            }
          read_header_attributes(fd);
          SKIP2;
          swap_Nbyte((char *)&blksize1, 1, 4);
          swap_Nbyte((char *)&blksize2, 1, 4);

          swap_header();

#ifdef COMBINETYPES
          header.npartTotal[3] += header.npartTotal[4] + header.npartTotal[5];
          header.npart[3] += header.npart[4] + header.npart[5];
          header.npartTotal[4] = 0;
          header.npartTotal[5] = 0;
          header.npart[4]      = 0;
          header.npart[5]      = 0;
#endif /* #ifdef COMBINETYPES */
        }

#ifdef HAVE_HDF5
      if(All.ICFormat == 3)
        {
          read_header_attributes_in_hdf5(fname);

          hdf5_file = my_H5Fopen(fname, H5F_ACC_RDONLY, H5P_DEFAULT);
          if(hdf5_file < 0)
            terminate("cannot read initial conditions file %s", fname);

          for(type = 0; type < NTYPES; type++)
            {
              if(header.npart[type] > 0 && (readTypes & (1 << type)))
                {
                  sprintf(buf, "/PartType%d", type);
                  hdf5_grp[type] = my_H5Gopen(hdf5_file, buf);
                }
              if(!(readTypes & (1 << type)))
                {
                  // Override particle number in file. If we don't
                  // read the type, both npart and npartTotal will be 0
                  header.npartTotal[type]         = 0;
                  header.npart[type]              = 0;
                  header.npartTotalHighWord[type] = 0;
                  header.mass[type]               = 0;
                }
            }
        }
#endif /* #ifdef HAVE_HDF5 */

      for(task = readTask + 1; task <= lastTask; task++)
        MPI_Ssend(&header, sizeof(header), MPI_BYTE, task, TAG_HEADER, MPI_COMM_WORLD);
    }
  else
    MPI_Recv(&header, sizeof(header), MPI_BYTE, readTask, TAG_HEADER, MPI_COMM_WORLD, &status);

#ifdef INPUT_IN_DOUBLEPRECISION
  if(header.flag_doubleprecision == 0)
    {
      sprintf(buf, "\nProblem: Code compiled with INPUT_IN_DOUBLEPRECISION, but input files are in single precision!\n");
      terminate(buf);
    }
#else  /* #ifdef INPUT_IN_DOUBLEPRECISION */
  if(header.flag_doubleprecision)
    {
      sprintf(buf, "\nProblem: Code not compiled with INPUT_IN_DOUBLEPRECISION, but input files are in double precision!\n");
      terminate(buf);
    }
#endif /* #ifdef INPUT_IN_DOUBLEPRECISION #else */

  if(ThisTask == readTask)
    {
      if(filenr == 0)
        mpi_printf(
            "\nREADIC: filenr=%d, '%s' contains:\n"
            "READIC: Type 0 (gas):   %8d  (tot=%15lld) masstab= %g\n"
            "READIC: Type 1 (halo):  %8d  (tot=%15lld) masstab= %g\n"
            "READIC: Type 2 (disk):  %8d  (tot=%15lld) masstab= %g\n"
            "READIC: Type 3 (bulge): %8d  (tot=%15lld) masstab= %g\n"
            "READIC: Type 4 (stars): %8d  (tot=%15lld) masstab= %g\n"
            "READIC: Type 5 (bndry): %8d  (tot=%15lld) masstab= %g\n\n",
            filenr, fname, header.npart[0], header.npartTotal[0] + (((long long)header.npartTotalHighWord[0]) << 32), All.MassTable[0],
            header.npart[1], header.npartTotal[1] + (((long long)header.npartTotalHighWord[1]) << 32), All.MassTable[1],
            header.npart[2], header.npartTotal[2] + (((long long)header.npartTotalHighWord[2]) << 32), All.MassTable[2],
            header.npart[3], header.npartTotal[3] + (((long long)header.npartTotalHighWord[3]) << 32), All.MassTable[3],
            header.npart[4], header.npartTotal[4] + (((long long)header.npartTotalHighWord[4]) << 32), All.MassTable[4],
            header.npart[5], header.npartTotal[5] + (((long long)header.npartTotalHighWord[5]) << 32), All.MassTable[5]);
    }

  /* to collect the gas particles all at the beginning (in case several
     snapshot files are read on the current CPU) we move the collisionless
     particles such that a gap of the right size is created */

  for(type = 0, nall = 0; type < NTYPES; type++)
    {
      n_in_file       = header.npart[type];
      ntask           = lastTask - readTask + 1;
      n_for_this_task = n_in_file / ntask;
      if((ThisTask - readTask) < (n_in_file % ntask))
        n_for_this_task++;

      nall += n_for_this_task;
    }

  memmove(&P[NumGas + nall], &P[NumGas], (NumPart - NumGas) * sizeof(struct particle_data));
  nstart = NumGas;

  for(bnr = 0; bnr < 1000; bnr++)
    {
      blocknr = (enum iofields)bnr;

      if(blocknr == IO_LASTENTRY)
        {
#if defined(RECOMPUTE_POTENTIAL_IN_SNAPSHOT)
          int pc = nstart;

          for(int type = 0; type < NTYPES; type++)
            {
              int n_in_file = header.npart[type];

              long long nprevious = 0;
              for(int t = 0; t < type; t++)
                nprevious += header.npartTotal[t] + (((long long)header.npartTotalHighWord[t]) << 32);

              for(int nr = 0; nr < filenr; nr++)
                nprevious += ntype_in_files[nr].npart[type];

              for(int task = readTask; task <= lastTask; task++)
                {
                  int n_for_this_task = n_in_file / ntask;
                  if((task - readTask) < (n_in_file % ntask))
                    n_for_this_task++;

                  if(ThisTask == task)
                    {
                      for(int i = 0; i < n_for_this_task; i++)
                        P[pc++].FileOrder = nprevious++;
                    }
                  else
                    nprevious += n_for_this_task;
                }
            }
#endif /* #if defined(RECOMPUTE_POTENTIAL_IN_SNAPSHOT) */
          break;
        }

      /* proceed reading this field only if we are expecting it */
      if(blockpresent(blocknr, 0))
        {
          if(ThisTask == readTask)
            {
              get_dataset_name(blocknr, buf);
              if(filenr == 0)
                mpi_printf("READIC: reading block %d (%s)...\n", blocknr, buf);
              myflush(stdout);
            }

          bytes_per_blockelement = get_bytes_per_blockelement(blocknr, 1);

          blockmaxlen = (int)(COMMBUFFERSIZE / bytes_per_blockelement);

          npart = get_particles_in_block(blocknr, &typelist[0]);

          if(npart > 0)
            {
              if(ThisTask == readTask)
                {
                  if(All.ICFormat == 2)
                    {
                      SKIP;
                      my_fread(&label, sizeof(char), 4, fd);
                      my_fread(&nextblock, sizeof(int), 1, fd);
                      swap_Nbyte((char *)&nextblock, 1, 4);
                      printf("Reading header => '%c%c%c%c' (%d byte)\n", label[0], label[1], label[2], label[3], nextblock);
                      SKIP2;

                      get_Tab_IO_Label(blocknr, expected_label);
                      if(strncmp(label, expected_label, 4) != 0)
                        {
                          sprintf(buf, "incorrect block-structure!\nexpected '%c%c%c%c' but found '%c%c%c%c'\n", expected_label[0],
                                  expected_label[1], expected_label[2], expected_label[3], label[0], label[1], label[2], label[3]);
                          terminate(buf);
                        }
                    }

                  if(All.ICFormat == 1 || All.ICFormat == 2)
                    SKIP;
                }

              for(type = 0, offset = 0; type < NTYPES; type++)
                {
                  n_in_file = header.npart[type];
#ifdef HAVE_HDF5
                  pcsum = 0;
#endif /* #ifdef HAVE_HDF5 */
                  if(typelist[type] == 0)
                    {
                      /* we are expecting (npart>0) this block, but not for this particle type */
                      n_for_this_task = n_in_file / ntask;
                      if((ThisTask - readTask) < (n_in_file % ntask))
                        n_for_this_task++;

                      offset += n_for_this_task;
                    }
                  else
                    {
                      /* we are expecting (npart>0) this block for this particle type, read or recv */
                      for(task = readTask; task <= lastTask; task++)
                        {
                          n_for_this_task = n_in_file / ntask;
                          if((task - readTask) < (n_in_file % ntask))
                            n_for_this_task++;

                          if(task == ThisTask)
                            if(NumPart + n_for_this_task > All.MaxPart)
                              terminate("too many particles. %d %d %d\n", NumPart, n_for_this_task, All.MaxPart);

                          /* blocked load to fit in finite size of CommBuffer */
                          do
                            {
                              pc = n_for_this_task;

                              if(pc > blockmaxlen)
                                pc = blockmaxlen;

                              if(ThisTask == readTask)
                                {
                                  if(All.ICFormat == 1 || All.ICFormat == 2)
                                    my_fread(CommBuffer, bytes_per_blockelement, pc, fd);
#ifdef HAVE_HDF5
                                  if(All.ICFormat == 3 && pc > 0)
                                    {
                                      /* configure HDF5 dataspaces and hyperslab selection */
                                      dims[0] = header.npart[type];
                                      dims[1] = get_values_per_blockelement(blocknr);
                                      if(dims[1] == 1)
                                        rank = 1;
                                      else
                                        rank = 2;

                                      hdf5_dataspace_in_file = my_H5Screate_simple(rank, dims, NULL);

                                      dims[0]                  = pc;
                                      hdf5_dataspace_in_memory = my_H5Screate_simple(rank, dims, NULL);

                                      start[0] = pcsum;
                                      start[1] = 0;

                                      count[0] = pc;
                                      count[1] = get_values_per_blockelement(blocknr);
                                      pcsum += pc;

                                      my_H5Sselect_hyperslab(hdf5_dataspace_in_file, H5S_SELECT_SET, start, NULL, count, NULL);

                                      switch(get_datatype_in_block(blocknr, 1))
                                        {
                                          case FILE_INT:
                                            hdf5_datatype = my_H5Tcopy(H5T_NATIVE_UINT);
                                            break;
                                          case FILE_MY_IO_FLOAT:
#ifdef INPUT_IN_DOUBLEPRECISION
                                            hdf5_datatype = my_H5Tcopy(H5T_NATIVE_DOUBLE);
#else  /* #ifdef INPUT_IN_DOUBLEPRECISION */
                                            hdf5_datatype = my_H5Tcopy(H5T_NATIVE_FLOAT);
#endif /* #ifdef INPUT_IN_DOUBLEPRECISION #else */
                                            break;
                                          case FILE_MY_ID_TYPE:
#ifdef LONGIDS
                                            hdf5_datatype = my_H5Tcopy(H5T_NATIVE_UINT64);
#else  /* #ifdef LONGIDS */
                                            hdf5_datatype = my_H5Tcopy(H5T_NATIVE_UINT32);
#endif /* #ifdef LONGIDS #else */
                                            break;
                                          case FILE_DOUBLE:
                                            hdf5_datatype = my_H5Tcopy(H5T_NATIVE_DOUBLE);
                                            break;
                                          case FILE_FLOAT:
                                            hdf5_datatype = my_H5Tcopy(H5T_NATIVE_FLOAT);
                                            break;
                                          default:
                                            terminate("can't process this input type");
                                            break;
                                        }

                                      /* test if HDF5 dataset is actually present */
                                      get_dataset_name(blocknr, buf);

                                      hdf5_dataset = my_H5Dopen_if_existing(hdf5_grp[type], buf);

                                      if(hdf5_dataset < 0)
                                        {
                                          // no, pad with zeros
                                          if((ThisTask == readTask) && (task == ThisTask))
                                            mpi_printf("\tDataset %s not present for particle type %d, using zero.\n", buf, type);
                                          memset(CommBuffer, 0, dims[0] * dims[1] * my_H5Tget_size(hdf5_datatype));
                                        }
                                      else
                                        {
                                          // yes, read into CommBuffer
                                          my_H5Dread(hdf5_dataset, hdf5_datatype, hdf5_dataspace_in_memory, hdf5_dataspace_in_file,
                                                     H5P_DEFAULT, CommBuffer, buf);
                                          my_H5Dclose(hdf5_dataset, buf);
                                        }
                                      my_H5Tclose(hdf5_datatype);
                                      my_H5Sclose(hdf5_dataspace_in_memory, H5S_SIMPLE);
                                      my_H5Sclose(hdf5_dataspace_in_file, H5S_SIMPLE);

                                    } /* All.ICFormat == 3 */
#endif                                /* #ifdef HAVE_HDF5 */
                                }

                              if(ThisTask == readTask && task != readTask && pc > 0)
                                MPI_Ssend(CommBuffer, bytes_per_blockelement * pc, MPI_BYTE, task, TAG_PDATA, MPI_COMM_WORLD);

                              if(ThisTask != readTask && task == ThisTask && pc > 0)
                                MPI_Recv(CommBuffer, bytes_per_blockelement * pc, MPI_BYTE, readTask, TAG_PDATA, MPI_COMM_WORLD,
                                         &status);

                              /* copy CommBuffer contents into actual particle data structs */
                              if(ThisTask == task)
                                {
                                  empty_read_buffer(blocknr, nstart + offset, pc, type);

                                  offset += pc;
                                }

                              n_for_this_task -= pc;
                            } /* do */
                          while(n_for_this_task > 0);

                        } /* task loop */
                    }     /* typelist[type] > 0 */
                }         /* type loop */

              if(ThisTask == readTask)
                {
                  if(All.ICFormat == 1 || All.ICFormat == 2)
                    {
                      SKIP2;
                      swap_Nbyte((char *)&blksize1, 1, 4);
                      swap_Nbyte((char *)&blksize2, 1, 4);
                      if(blksize1 != blksize2)
                        {
                          sprintf(buf, "incorrect block-sizes detected!\n Task=%d   blocknr=%d  blksize1=%d  blksize2=%d\n", ThisTask,
                                  blocknr, blksize1, blksize2);
                          if(blocknr == IO_ID)
                            {
                              strcat(buf, "Possible mismatch of 32bit and 64bit ID's in IC file and AREPO compilation !\n");
                            }
                          terminate(buf);
                        }
                    }
                }

            } /* npart > 0 */
        }     /* blockpresent */
    }         /* blocknr loop */

  for(type = 0; type < NTYPES; type++)
    {
      n_in_file = header.npart[type];

      n_for_this_task = n_in_file / ntask;
      if((ThisTask - readTask) < (n_in_file % ntask))
        n_for_this_task++;

      NumPart += n_for_this_task;

      if(type == 0)
        NumGas += n_for_this_task;
    }

  if(ThisTask == readTask)
    {
      if(All.ICFormat == 1 || All.ICFormat == 2)
        fclose(fd);
#ifdef HAVE_HDF5
      if(All.ICFormat == 3)
        {
          for(type = NTYPES - 1; type >= 0; type--)
            if(header.npart[type] > 0)
              {
                sprintf(buf, "/PartType%d", type);
                my_H5Gclose(hdf5_grp[type], buf);
              }
          my_H5Fclose(hdf5_file, fname);
        }
#endif /* #ifdef HAVE_HDF5 */
    }
}

/*! \brief Determines on how many files a given snapshot is distributed.
 *
 *  \param[in] fname File name of the snapshot as given in the parameter file.
 *
 *  \return Number of files; -1: could not find files.
 */
int find_files(const char *fname)
{
  FILE *fd;
  char buf[200], buf1[200];
  int dummy;

  sprintf(buf, "%s.%d", fname, 0);
  sprintf(buf1, "%s", fname);

  if(All.ICFormat == 3)
    {
      sprintf(buf, "%s.%d.hdf5", fname, 0);
      sprintf(buf1, "%s.hdf5", fname);
    }

#ifndef HAVE_HDF5
  if(All.ICFormat == 3)
    {
      mpi_terminate("Code wasn't compiled with HDF5 support enabled!\n");
    }
#endif /* #ifndef HAVE_HDF5 */

  header.num_files = 0;

  if(ThisTask == 0)
    {
      if((fd = fopen(buf, "r")))
        {
          if(All.ICFormat == 1 || All.ICFormat == 2)
            {
              if(All.ICFormat == 2)
                {
                  my_fread(&dummy, sizeof(dummy), 1, fd);
                  swap_file = dummy;
                  my_fread(&dummy, sizeof(dummy), 1, fd);
                  my_fread(&dummy, sizeof(dummy), 1, fd);
                  my_fread(&dummy, sizeof(dummy), 1, fd);
                }

              my_fread(&dummy, sizeof(dummy), 1, fd);
              if(All.ICFormat == 1)
                {
                  if(dummy == 256)
                    swap_file = 8;
                  else
                    swap_file = dummy;
                }
              read_header_attributes(fd);

              swap_header();

#ifdef COMBINETYPES
              header.npartTotal[3] += header.npartTotal[4] + header.npartTotal[5];
              header.npart[3] += header.npart[4] + header.npart[5];
              header.npartTotal[4] = 0;
              header.npartTotal[5] = 0;
              header.npart[4]      = 0;
              header.npart[5]      = 0;
#endif /* #ifdef COMBINETYPES */

              my_fread(&dummy, sizeof(dummy), 1, fd);
            }
          fclose(fd);

#ifdef HAVE_HDF5
          if(All.ICFormat == 3)
            read_header_attributes_in_hdf5(buf);
#endif /* #ifdef HAVE_HDF5 */
        }
    }

  MPI_Bcast(&swap_file, sizeof(swap_file), MPI_BYTE, 0, MPI_COMM_WORLD);
  MPI_Bcast(&header, sizeof(header), MPI_BYTE, 0, MPI_COMM_WORLD);

  if(header.num_files < 0)
    terminate("header.num_files < 0");
  if(header.num_files > 100000)
    terminate("header.num_files=%d read from %s does not make sense - header possibly corrupt.", header.num_files, buf);
  if(header.num_files > 0)
    return header.num_files;

  if(ThisTask == 0)
    {
      if((fd = fopen(buf1, "r")))
        {
          if(All.ICFormat == 1 || All.ICFormat == 2)
            {
              if(All.ICFormat == 2)
                {
                  my_fread(&dummy, sizeof(dummy), 1, fd);
                  swap_file = dummy;
                  my_fread(&dummy, sizeof(dummy), 1, fd);
                  my_fread(&dummy, sizeof(dummy), 1, fd);
                  my_fread(&dummy, sizeof(dummy), 1, fd);
                }

              my_fread(&dummy, sizeof(dummy), 1, fd);
              if(All.ICFormat == 1)
                {
                  if(dummy == 256)
                    swap_file = 8;
                  else
                    swap_file = dummy;
                }
              read_header_attributes(fd);
              swap_header();

#ifdef COMBINETYPES
              header.npartTotal[3] += header.npartTotal[4] + header.npartTotal[5];
              header.npart[3] += header.npart[4] + header.npart[5];
              header.npartTotal[4] = 0;
              header.npartTotal[5] = 0;
              header.npart[4]      = 0;
              header.npart[5]      = 0;
#endif /* #ifdef COMBINETYPES */

              my_fread(&dummy, sizeof(dummy), 1, fd);
            }
          fclose(fd);

#ifdef HAVE_HDF5
          if(All.ICFormat == 3)
            read_header_attributes_in_hdf5(buf1);
#endif /* #ifdef HAVE_HDF5 */

          header.num_files = 1;
        }
    }

  MPI_Bcast(&swap_file, sizeof(swap_file), MPI_BYTE, 0, MPI_COMM_WORLD);
  MPI_Bcast(&header, sizeof(header), MPI_BYTE, 0, MPI_COMM_WORLD);

  if(header.num_files > 0)
    return header.num_files;

  mpi_terminate("\nCan't find initial conditions file, neither as '%s'\nnor as '%s'\n", buf, buf1);
  return -1;
}

/*! \brief This function assigns a certain number of tasks to each file.
 *
 *  These tasks are containing the content of that file after the ICs have been
 *  read. The number of tasks per file is as homogeneous as possible.
 *  The number of files may at most be equal to the number of tasks.
 *
 *  \param[in] nfiles Number of files of which the snapshot is distributed.
 *  \param[in] filenr Contains the file number to which this task belongs.
 *  \param[in] master The number of the task responsible to read the file.
 *  \param[in] last Number of the last task belonging to the same file as this
 *             task.
 *
 *  \return void
 */
void distribute_file(int nfiles, int firstfile, int firsttask, int lasttask, int *filenr, int *master, int *last)
{
  int i, group;
  int tasks_per_file = NTask / nfiles;
  int tasks_left     = NTask % nfiles;

  if(tasks_left == 0)
    {
      group   = ThisTask / tasks_per_file;
      *master = group * tasks_per_file;
      *last   = (group + 1) * tasks_per_file - 1;
      *filenr = group;
      return;
    }

  double tpf = ((double)NTask) / nfiles;

  for(i = 0, *last = -1; i < nfiles; i++)
    {
      *master = *last + 1;
      *last   = (i + 1) * tpf;
      if(*last >= NTask)
        *last = *last - 1;
      if(*last < *master)
        terminate("last < master");
      *filenr = i;

      if(i == nfiles - 1)
        *last = NTask - 1;

      if(ThisTask >= *master && ThisTask <= *last)
        return;
    }
}

#ifdef HAVE_HDF5
/*! \brief The error handler used during the loading of the hdf5 header.
 *
 *  \param[in] unused The parameter is not used, but it is necessary for
 *             compatibility with the HDF5 library.
 *  \return 1 if the write error is tolerated, otherwise the run is terminated.
 */
herr_t hdf5_header_error_handler(void *unused)
{
#ifdef TOLERATE_WRITE_ERROR
  write_error(3, 0, 0);
  return 1;
#else
  terminate("Failed to read HDF5 header attribute. Probably your file is corrupt.\n");
  return 0;
#endif
}

/*! \brief This function reads the snapshot header in case of hdf5 files
 *        (i.e. format 3).
 *
 *  \param[in] fname File name of the snapshot as given in the parameter file.
 *
 *  \return void
 */
void read_header_attributes_in_hdf5(const char *fname)
{
  hid_t hdf5_file, hdf5_headergrp, hdf5_attribute;
  hssize_t scalar_attr_dim = 1;
  hssize_t vector_attr_dim = NTYPES;

  hdf5_file      = my_H5Fopen(fname, H5F_ACC_RDONLY, H5P_DEFAULT);
  hdf5_headergrp = my_H5Gopen(hdf5_file, "/Header");

  hdf5_attribute = my_H5Aopen_name(hdf5_headergrp, "NumPart_ThisFile");
  my_H5Aread(hdf5_attribute, H5T_NATIVE_INT, header.npart, "NumPart_ThisFile", vector_attr_dim);
  my_H5Aclose(hdf5_attribute, "NumPart_ThisFile");

  hdf5_attribute = my_H5Aopen_name(hdf5_headergrp, "NumPart_Total");
  my_H5Aread(hdf5_attribute, H5T_NATIVE_UINT, header.npartTotal, "NumPart_Total", vector_attr_dim);
  my_H5Aclose(hdf5_attribute, "NumPart_Total");

  hdf5_attribute = my_H5Aopen_name(hdf5_headergrp, "NumPart_Total_HighWord");
  my_H5Aread(hdf5_attribute, H5T_NATIVE_UINT, header.npartTotalHighWord, "NumPart_Total_HighWord", vector_attr_dim);
  my_H5Aclose(hdf5_attribute, "NumPart_Total_HighWord");

  hdf5_attribute = my_H5Aopen_name(hdf5_headergrp, "MassTable");
  my_H5Aread(hdf5_attribute, H5T_NATIVE_DOUBLE, header.mass, "MassTable", vector_attr_dim);
  my_H5Aclose(hdf5_attribute, "MassTable");

  hdf5_attribute = my_H5Aopen_name(hdf5_headergrp, "Time");
  my_H5Aread(hdf5_attribute, H5T_NATIVE_DOUBLE, &header.time, "Time", scalar_attr_dim);
  my_H5Aclose(hdf5_attribute, "Time");

  hdf5_attribute = my_H5Aopen_name(hdf5_headergrp, "Redshift");
  my_H5Aread(hdf5_attribute, H5T_NATIVE_DOUBLE, &header.redshift, "Redshift", scalar_attr_dim);
  my_H5Aclose(hdf5_attribute, "Redshift");

  hdf5_attribute = my_H5Aopen_name(hdf5_headergrp, "BoxSize");
  my_H5Aread(hdf5_attribute, H5T_NATIVE_DOUBLE, &header.BoxSize, "BoxSize", scalar_attr_dim);
  my_H5Aclose(hdf5_attribute, "BoxSize");

  hdf5_attribute = my_H5Aopen_name(hdf5_headergrp, "NumFilesPerSnapshot");
  my_H5Aread(hdf5_attribute, H5T_NATIVE_INT, &header.num_files, "NumFilesPerSnapshot", scalar_attr_dim);
  my_H5Aclose(hdf5_attribute, "NumFilesPerSnapshot");

  hdf5_attribute = my_H5Aopen_name(hdf5_headergrp, "Omega0");
  my_H5Aread(hdf5_attribute, H5T_NATIVE_DOUBLE, &header.Omega0, "Omega0", scalar_attr_dim);
  my_H5Aclose(hdf5_attribute, "Omega0");

  hdf5_attribute = my_H5Aopen_name(hdf5_headergrp, "OmegaLambda");
  my_H5Aread(hdf5_attribute, H5T_NATIVE_DOUBLE, &header.OmegaLambda, "OmegaLambda", scalar_attr_dim);
  my_H5Aclose(hdf5_attribute, "OmegaLambda");

  hdf5_attribute = my_H5Aopen_name(hdf5_headergrp, "HubbleParam");
  my_H5Aread(hdf5_attribute, H5T_NATIVE_DOUBLE, &header.HubbleParam, "HubbleParam", scalar_attr_dim);
  my_H5Aclose(hdf5_attribute, "HubbleParam");

  hdf5_attribute = my_H5Aopen_name(hdf5_headergrp, "Flag_Sfr");
  my_H5Aread(hdf5_attribute, H5T_NATIVE_INT, &header.flag_sfr, "Flag_Sfr", scalar_attr_dim);
  my_H5Aclose(hdf5_attribute, "Flag_Sfr");

  hdf5_attribute = my_H5Aopen_name(hdf5_headergrp, "Flag_Cooling");
  my_H5Aread(hdf5_attribute, H5T_NATIVE_INT, &header.flag_cooling, "Flag_Cooling", scalar_attr_dim);
  my_H5Aclose(hdf5_attribute, "Flag_Cooling");

  hdf5_attribute = my_H5Aopen_name(hdf5_headergrp, "Flag_StellarAge");
  my_H5Aread(hdf5_attribute, H5T_NATIVE_INT, &header.flag_stellarage, "Flag_StellarAge", scalar_attr_dim);
  my_H5Aclose(hdf5_attribute, "Flag_StellarAge");

  hdf5_attribute = my_H5Aopen_name(hdf5_headergrp, "Flag_Metals");
  my_H5Aread(hdf5_attribute, H5T_NATIVE_INT, &header.flag_metals, "Flag_Metals", scalar_attr_dim);
  my_H5Aclose(hdf5_attribute, "Flag_Metals");

  hdf5_attribute = my_H5Aopen_name(hdf5_headergrp, "Flag_Feedback");
  my_H5Aread(hdf5_attribute, H5T_NATIVE_INT, &header.flag_feedback, "Flag_Feedback", scalar_attr_dim);
  my_H5Aclose(hdf5_attribute, "Flag_Feedback");

  hdf5_attribute = my_H5Aopen_name(hdf5_headergrp, "Flag_DoublePrecision");
  my_H5Aread(hdf5_attribute, H5T_NATIVE_INT, &header.flag_doubleprecision, "Flag_DoublePrecision", scalar_attr_dim);
  my_H5Aclose(hdf5_attribute, "Flag_DoublePrecision");

  my_H5Gclose(hdf5_headergrp, "/Header");
  my_H5Fclose(hdf5_file, fname);
}
#endif /* #ifdef HAVE_HDF5 */

/*! \brief Reads the snapshot header in case of non-hdf5 files (i.e. formats 1
 *         and 2).
 *
 * \param[in] *fd Pointer to snapshot file.
 *
 * \return void
 */
void read_header_attributes(FILE *fd)
{
#ifdef NTYPES_ICS
  int type;
  if(RestartFlag == 0)
    {
      my_fread(&header_ICs, sizeof(header_ICs), 1, fd);

      for(type = 0; type < NTYPES_ICS; type++)
        {
          header.npart[type]              = header_ICs.npart[type];
          header.mass[type]               = header_ICs.mass[type];
          header.npartTotal[type]         = header_ICs.npartTotal[type];
          header.npartTotalHighWord[type] = header_ICs.npartTotalHighWord[type];
        }
      for(type = NTYPES_ICS; type < NTYPES; type++)
        {
          header.npart[type]              = 0;
          header.mass[type]               = 0;
          header.npartTotal[type]         = 0;
          header.npartTotalHighWord[type] = 0;
        }

      header.time                      = header_ICs.time;
      header.redshift                  = header_ICs.redshift;
      header.flag_sfr                  = header_ICs.flag_sfr;
      header.flag_feedback             = header_ICs.flag_feedback;
      header.flag_cooling              = header_ICs.flag_cooling;
      header.num_files                 = header_ICs.num_files;
      header.BoxSize                   = header_ICs.BoxSize;
      header.Omega0                    = header_ICs.Omega0;
      header.OmegaLambda               = header_ICs.OmegaLambda;
      header.HubbleParam               = header_ICs.HubbleParam;
      header.flag_stellarage           = header_ICs.flag_stellarage;
      header.flag_metals               = header_ICs.flag_metals;
      header.flag_entropy_instead_u    = header_ICs.flag_entropy_instead_u;
      header.flag_doubleprecision      = header_ICs.flag_doubleprecision;
      header.flag_lpt_ics              = header_ICs.flag_lpt_ics;
      header.lpt_scalingfactor         = header_ICs.lpt_scalingfactor;
      header.flag_tracer_field         = header_ICs.flag_tracer_field;
      header.composition_vector_length = header_ICs.composition_vector_length;
    }
  else
    my_fread(&header, sizeof(header), 1, fd);
#else  /* #ifdef NTYPES_ICS */
  my_fread(&header, sizeof(header), 1, fd);
#endif /* #ifdef NTYPES_ICS #else */
}

/*! \brief Swaps endiannes of data.
 *
 * \param[in, out] data Pointer to the data.
 * \param[in] n Number of elements to swap.
 * \param[in] m Size of single element to swap: int, float = 4; double = 8.
 *
 * \return void
 */
void swap_Nbyte(char *data, int n, int m)
{
  int i, j;
  char old_data[16];

  if(swap_file != 8)
    {
      for(j = 0; j < n; j++)
        {
          memcpy(&old_data[0], &data[j * m], m);
          for(i = 0; i < m; i++)
            {
              data[j * m + i] = old_data[m - i - 1];
            }
        }
    }
}

/*! \brief Swaps the endianness of the snapshot header.
 *
 *  \return void
 */
void swap_header()
{
  swap_Nbyte((char *)&header.npart, NTYPES, 4);
  swap_Nbyte((char *)&header.mass, NTYPES, 8);
  swap_Nbyte((char *)&header.time, 1, 8);
  swap_Nbyte((char *)&header.redshift, 1, 8);
  swap_Nbyte((char *)&header.flag_sfr, 1, 4);
  swap_Nbyte((char *)&header.flag_feedback, 1, 4);
  swap_Nbyte((char *)&header.npartTotal, NTYPES, 4);
  swap_Nbyte((char *)&header.flag_cooling, 1, 4);
  swap_Nbyte((char *)&header.num_files, 1, 4);
  swap_Nbyte((char *)&header.BoxSize, 1, 8);
  swap_Nbyte((char *)&header.Omega0, 1, 8);
  swap_Nbyte((char *)&header.OmegaLambda, 1, 8);
  swap_Nbyte((char *)&header.HubbleParam, 1, 8);
  swap_Nbyte((char *)&header.flag_stellarage, 1, 4);
  swap_Nbyte((char *)&header.flag_metals, 1, 4);
  swap_Nbyte((char *)&header.npartTotalHighWord, NTYPES, 4);
  swap_Nbyte((char *)&header.flag_entropy_instead_u, 1, 4);
  swap_Nbyte((char *)&header.flag_doubleprecision, 1, 4);
  swap_Nbyte((char *)&header.flag_lpt_ics, 1, 4);
  swap_Nbyte((char *)&header.lpt_scalingfactor, 1, 4);
  swap_Nbyte((char *)&header.flag_tracer_field, 1, 4);
  swap_Nbyte((char *)&header.composition_vector_length, 1, 4);
}

#ifdef TILE_ICS
/*! \brief Duplicates ICs and lines TileICsFactor of them up in each dimension.
 *
 *  \return void
 */
void tile_ics(void)
{
  mpi_printf("TILE_ICS: tiling by a factor of %d...\n", All.TileICsFactor);

  /* allocate memory for new particles */
  domain_resize_storage(NumPart * (All.TileICsFactor * All.TileICsFactor * All.TileICsFactor - 1),
                        NumGas * (All.TileICsFactor * All.TileICsFactor * All.TileICsFactor - 1), 0);

  /* tile gas particles at the beginning of P[] */
  int N_others = NumPart - NumGas;
  memmove(&P[NumGas * All.TileICsFactor * All.TileICsFactor * All.TileICsFactor], &P[NumGas], N_others * sizeof(struct particle_data));
  int i, j, ix, iy = 0, iz = 0;
  for(i = 0; i < NumGas; i++)
    {
      for(ix = 0; ix < All.TileICsFactor; ix++)
        {
#ifndef ONEDIMS
          for(iy = 0; iy < All.TileICsFactor; iy++)
#endif /* #ifndef ONEDIMS */
            {
#if !defined(TWODIMS) && !defined(ONEDIMS)
              for(iz = 0; iz < All.TileICsFactor; iz++)
#endif /* #if !defined(TWODIMS) && !defined(ONEDIMS) */
                {
                  if(ix == 0 && iy == 0 && iz == 0)
                    continue;
                  j       = i + NumGas * ix + NumGas * All.TileICsFactor * iy + NumGas * All.TileICsFactor * All.TileICsFactor * iz;
                  P[j]    = P[i];
                  P[j].ID = P[i].ID + IDS_OFFSET * ix + IDS_OFFSET * All.TileICsFactor * iy +
                            IDS_OFFSET * All.TileICsFactor * All.TileICsFactor * iz;
                  P[j].Pos[0] += All.BoxSize / All.TileICsFactor * ix;
                  P[j].Pos[1] += All.BoxSize / All.TileICsFactor * iy;
                  P[j].Pos[2] += All.BoxSize / All.TileICsFactor * iz;
                  SphP[j] = SphP[i];
                }
            }
        }
    }
  /* tile the other particle types */
  iy = 0;
  iz = 0;
  for(i = NumGas * All.TileICsFactor * All.TileICsFactor * All.TileICsFactor;
      i < NumGas * All.TileICsFactor * All.TileICsFactor * All.TileICsFactor + N_others; i++)
    {
      for(ix = 0; ix < All.TileICsFactor; ix++)
        {
#ifndef ONEDIMS
          for(iy = 0; iy < All.TileICsFactor; iy++)
#endif /* #ifndef ONEDIMS */
            {
#if !defined(TWODIMS) && !defined(ONEDIMS)
              for(iz = 0; iz < All.TileICsFactor; iz++)
#endif /* #if !defined(TWODIMS) && !defined(ONEDIMS) */
                {
                  if(ix == 0 && iy == 0 && iz == 0)
                    continue;
                  j    = i + N_others * ix + N_others * All.TileICsFactor * iy + N_others * All.TileICsFactor * All.TileICsFactor * iz;
                  P[j] = P[i];
                  P[j].ID = P[i].ID + IDS_OFFSET * ix + IDS_OFFSET * All.TileICsFactor * iy +
                            IDS_OFFSET * All.TileICsFactor * All.TileICsFactor * iz;
                  P[j].Pos[0] += All.BoxSize / All.TileICsFactor * ix;
                  P[j].Pos[1] += All.BoxSize / All.TileICsFactor * iy;
                  P[j].Pos[2] += All.BoxSize / All.TileICsFactor * iz;
                }
            }
        }
    }

  NumGas *= All.TileICsFactor * All.TileICsFactor * All.TileICsFactor;
  NumPart *= All.TileICsFactor * All.TileICsFactor * All.TileICsFactor;
}
#endif /* #ifdef TILE_ICS */
