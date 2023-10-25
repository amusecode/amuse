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
 * \file        src/io/parameters.c
 * \date        05/2018
 * \brief       Parses the parameter file.
 * \details     This file contains the routine to parse the parameter file.
 *              Additionally the output list is also parsed.
 *              contains functions:
 *                void read_parameter_file(char *fname)
 *                void check_parameters()
 *                int read_outputlist(char *fname)
 *
 * \par Major modifications and contributions:
 *
 * - DD.MM.YYYY Description
 * - 06.05.2018 Prepared file for public release -- Rainer Weinberger
 */

#include <gsl/gsl_rng.h>
#include <math.h>
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <unistd.h>

#include "../main/allvars.h"
#include "../main/proto.h"

/*! \brief This function parses the parameter file.
 *
 *  Each parameter is defined by a keyword (`tag'), and can be either
 *  of type douple, int, or character string. Three arrays containing the name,
 *  type and address of the parameter are filled first. The routine then parses
 *  the parameter file and fills the referenced variables. The routine makes
 *  sure that each parameter appears exactly once in the parameter file,
 *  otherwise error messages are produced that complain about the missing
 *  parameters.
 *
 *  \param[in] fname The file name of the parameter file
 *
 *  \return void
 */
void read_parameter_file(char *fname)
{
#define REAL 1
#define STRING 2
#define INT 3

  FILE *fd, *fdout;
  char buf[MAXLEN_PARAM_TAG + MAXLEN_PARAM_VALUE + 200], buf1[MAXLEN_PARAM_TAG + 200], buf2[MAXLEN_PARAM_VALUE + 200],
      buf3[MAXLEN_PARAM_TAG + MAXLEN_PARAM_VALUE + 400];
  int i, j, nt;
  int id[MAX_PARAMETERS];
  void *addr[MAX_PARAMETERS];
  char tag[MAX_PARAMETERS][MAXLEN_PARAM_TAG];
  int param_handled[MAX_PARAMETERS];
  int errorFlag = 0;

  All.StarformationOn = 0; /* defaults */

  for(i = 0; i < MAX_PARAMETERS; i++)
    {
      param_handled[i] = 0;
    }

  if(sizeof(long long) != 8)
    {
      mpi_terminate("\nType `long long' is not 64 bit on this platform. Stopping.\n\n");
    }

  if(sizeof(int) != 4)
    {
      mpi_terminate("\nType `int' is not 32 bit on this platform. Stopping.\n\n");
    }

  if(sizeof(float) != 4)
    {
      mpi_terminate("\nType `float' is not 32 bit on this platform. Stopping.\n\n");
    }

  if(sizeof(double) != 8)
    {
      mpi_terminate("\nType `double' is not 64 bit on this platform. Stopping.\n\n");
    }

  if(ThisTask == 0) /* read parameter file on process 0 */
    {
      nt = 0;

      strcpy(tag[nt], "InitCondFile");
      addr[nt] = All.InitCondFile;
      id[nt++] = STRING;

      strcpy(tag[nt], "OutputDir");
      addr[nt] = All.OutputDir;
      id[nt++] = STRING;

#ifdef TOLERATE_WRITE_ERROR
      strcpy(tag[nt], "AlternativeOutputDir");
      addr[nt] = AlternativeOutputDir;
      id[nt++] = STRING;
#endif /* #ifdef TOLERATE_WRITE_ERROR */

      strcpy(tag[nt], "SnapshotFileBase");
      addr[nt] = All.SnapshotFileBase;
      id[nt++] = STRING;

      strcpy(tag[nt], "ResubmitCommand");
      addr[nt] = All.ResubmitCommand;
      id[nt++] = STRING;

      strcpy(tag[nt], "OutputListFilename");
      addr[nt] = All.OutputListFilename;
      id[nt++] = STRING;

      strcpy(tag[nt], "OutputListOn");
      addr[nt] = &All.OutputListOn;
      id[nt++] = INT;

      strcpy(tag[nt], "Omega0");
      addr[nt] = &All.Omega0;
      id[nt++] = REAL;

      strcpy(tag[nt], "OmegaBaryon");
      addr[nt] = &All.OmegaBaryon;
      id[nt++] = REAL;

      strcpy(tag[nt], "OmegaLambda");
      addr[nt] = &All.OmegaLambda;
      id[nt++] = REAL;

      strcpy(tag[nt], "HubbleParam");
      addr[nt] = &All.HubbleParam;
      id[nt++] = REAL;

      strcpy(tag[nt], "BoxSize");
      addr[nt] = &All.BoxSize;
      id[nt++] = REAL;

      strcpy(tag[nt], "PeriodicBoundariesOn");
      addr[nt] = &All.PeriodicBoundariesOn;
      id[nt++] = INT;

      strcpy(tag[nt], "MaxMemSize");
      addr[nt] = &All.MaxMemSize;
      id[nt++] = INT;

      strcpy(tag[nt], "TimeOfFirstSnapshot");
      addr[nt] = &All.TimeOfFirstSnapshot;
      id[nt++] = REAL;

      strcpy(tag[nt], "CpuTimeBetRestartFile");
      addr[nt] = &All.CpuTimeBetRestartFile;
      id[nt++] = REAL;

#ifdef REDUCE_FLUSH
      strcpy(tag[nt], "FlushCpuTimeDiff");
      addr[nt] = &All.FlushCpuTimeDiff;
      id[nt++] = REAL;
#endif /* #ifdef REDUCE_FLUSH */

      strcpy(tag[nt], "TimeBetStatistics");
      addr[nt] = &All.TimeBetStatistics;
      id[nt++] = REAL;

      strcpy(tag[nt], "TimeBegin");
      addr[nt] = &All.TimeBegin;
      id[nt++] = REAL;

      strcpy(tag[nt], "TimeMax");
      addr[nt] = &All.TimeMax;
      id[nt++] = REAL;

      strcpy(tag[nt], "TimeBetSnapshot");
      addr[nt] = &All.TimeBetSnapshot;
      id[nt++] = REAL;

      strcpy(tag[nt], "UnitVelocity_in_cm_per_s");
      addr[nt] = &All.UnitVelocity_in_cm_per_s;
      id[nt++] = REAL;

      strcpy(tag[nt], "UnitLength_in_cm");
      addr[nt] = &All.UnitLength_in_cm;
      id[nt++] = REAL;

      strcpy(tag[nt], "UnitMass_in_g");
      addr[nt] = &All.UnitMass_in_g;
      id[nt++] = REAL;

      strcpy(tag[nt], "ErrTolIntAccuracy");
      addr[nt] = &All.ErrTolIntAccuracy;
      id[nt++] = REAL;

      strcpy(tag[nt], "ErrTolTheta");
      addr[nt] = &All.ErrTolTheta;
      id[nt++] = REAL;

      strcpy(tag[nt], "ErrTolForceAcc");
      addr[nt] = &All.ErrTolForceAcc;
      id[nt++] = REAL;

      strcpy(tag[nt], "MaxSizeTimestep");
      addr[nt] = &All.MaxSizeTimestep;
      id[nt++] = REAL;

      strcpy(tag[nt], "MinSizeTimestep");
      addr[nt] = &All.MinSizeTimestep;
      id[nt++] = REAL;

      strcpy(tag[nt], "CourantFac");
      addr[nt] = &All.CourantFac;
      id[nt++] = REAL;

      strcpy(tag[nt], "LimitUBelowThisDensity");
      addr[nt] = &All.LimitUBelowThisDensity;
      id[nt++] = REAL;

      strcpy(tag[nt], "LimitUBelowCertainDensityToThisValue");
      addr[nt] = &All.LimitUBelowCertainDensityToThisValue;
      id[nt++] = REAL;

      strcpy(tag[nt], "DesNumNgb");
      addr[nt] = &All.DesNumNgb;
      id[nt++] = INT;

      strcpy(tag[nt], "MultipleDomains");
      addr[nt] = &All.MultipleDomains;
      id[nt++] = INT;

      strcpy(tag[nt], "TopNodeFactor");
      addr[nt] = &All.TopNodeFactor;
      id[nt++] = REAL;

      strcpy(tag[nt], "ActivePartFracForNewDomainDecomp");
      addr[nt] = &All.ActivePartFracForNewDomainDecomp;
      id[nt++] = REAL;

#ifdef SUBFIND
      strcpy(tag[nt], "DesLinkNgb");
      addr[nt] = &All.DesLinkNgb;
      id[nt++] = INT;

      strcpy(tag[nt], "ErrTolThetaSubfind");
      addr[nt] = &All.ErrTolThetaSubfind;
      id[nt++] = REAL;
#endif /* #ifdef SUBFIND */

#if defined(ISOTHERM_EQS)
      strcpy(tag[nt], "IsoSoundSpeed");
      addr[nt] = &All.IsoSoundSpeed;
      id[nt++] = REAL;
#endif /* #if defined(ISOTHERM_EQS) */

      strcpy(tag[nt], "MaxNumNgbDeviation");
      addr[nt] = &All.MaxNumNgbDeviation;
      id[nt++] = REAL;

      strcpy(tag[nt], "ComovingIntegrationOn");
      addr[nt] = &All.ComovingIntegrationOn;
      id[nt++] = INT;

      strcpy(tag[nt], "ICFormat");
      addr[nt] = &All.ICFormat;
      id[nt++] = INT;

      strcpy(tag[nt], "SnapFormat");
      addr[nt] = &All.SnapFormat;
      id[nt++] = INT;

      strcpy(tag[nt], "NumFilesPerSnapshot");
      addr[nt] = &All.NumFilesPerSnapshot;
      id[nt++] = INT;

      strcpy(tag[nt], "NumFilesWrittenInParallel");
      addr[nt] = &All.NumFilesWrittenInParallel;
      id[nt++] = INT;

      strcpy(tag[nt], "ResubmitOn");
      addr[nt] = &All.ResubmitOn;
      id[nt++] = INT;

      strcpy(tag[nt], "CoolingOn");
      addr[nt] = &All.CoolingOn;
      id[nt++] = INT;

      strcpy(tag[nt], "StarformationOn");
      addr[nt] = &All.StarformationOn;
      id[nt++] = INT;

      strcpy(tag[nt], "TypeOfTimestepCriterion");
      addr[nt] = &All.TypeOfTimestepCriterion;
      id[nt++] = INT;

      strcpy(tag[nt], "TypeOfOpeningCriterion");
      addr[nt] = &All.TypeOfOpeningCriterion;
      id[nt++] = INT;

      strcpy(tag[nt], "TimeLimitCPU");
      addr[nt] = &All.TimeLimitCPU;
      id[nt++] = REAL;

      strcpy(tag[nt], "GasSoftFactor");
      addr[nt] = &All.GasSoftFactor;
      id[nt++] = REAL;

      for(i = 0; i < NSOFTTYPES; i++)
        {
          char buf[100];
          sprintf(buf, "SofteningComovingType%d", i);
          strcpy(tag[nt], buf);
          addr[nt] = &All.SofteningComoving[i];
          id[nt++] = REAL;
        }

      for(i = 0; i < NSOFTTYPES; i++)
        {
          char buf[100];
          sprintf(buf, "SofteningMaxPhysType%d", i);
          strcpy(tag[nt], buf);
          addr[nt] = &All.SofteningMaxPhys[i];
          id[nt++] = REAL;
        }

      for(i = 0; i < NTYPES; i++)
        {
          char buf[100];
          sprintf(buf, "SofteningTypeOfPartType%d", i);
          strcpy(tag[nt], buf);
          addr[nt] = &All.SofteningTypeOfPartType[i];
          id[nt++] = INT;
        }

#ifdef ADAPTIVE_HYDRO_SOFTENING
      strcpy(tag[nt], "MinimumComovingHydroSoftening");
      addr[nt] = &All.MinimumComovingHydroSoftening;
      id[nt++] = REAL;

      strcpy(tag[nt], "AdaptiveHydroSofteningSpacing");
      addr[nt] = &All.AdaptiveHydroSofteningSpacing;
      id[nt++] = REAL;
#endif /* #ifdef ADAPTIVE_HYDRO_SOFTENING */

      strcpy(tag[nt], "GravityConstantInternal");
      addr[nt] = &All.GravityConstantInternal;
      id[nt++] = REAL;

      strcpy(tag[nt], "InitGasTemp");
      addr[nt] = &All.InitGasTemp;
      id[nt++] = REAL;

      strcpy(tag[nt], "MinGasTemp");
      addr[nt] = &All.MinGasTemp;
      id[nt++] = REAL;

      strcpy(tag[nt], "MinEgySpec");
      addr[nt] = &All.MinEgySpec;
      id[nt++] = REAL;

      strcpy(tag[nt], "MinimumDensityOnStartUp");
      addr[nt] = &All.MinimumDensityOnStartUp;
      id[nt++] = REAL;

#ifdef NODEREFINE_BACKGROUND_GRID
      strcpy(tag[nt], "MeanVolume");
      addr[nt] = &All.MeanVolume;
      id[nt++] = REAL;
#endif /* #ifdef NODEREFINE_BACKGROUND_GRID */

#ifndef VORONOI_STATIC_MESH
#ifdef REGULARIZE_MESH_FACE_ANGLE
      strcpy(tag[nt], "CellMaxAngleFactor");
      addr[nt] = &All.CellMaxAngleFactor;
      id[nt++] = REAL;
#else  /* #ifdef REGULARIZE_MESH_FACE_ANGLE */
      strcpy(tag[nt], "CellShapingFactor");
      addr[nt] = &All.CellShapingFactor;
      id[nt++] = REAL;
#endif /* #ifdef REGULARIZE_MESH_FACE_ANGLE #else */

      strcpy(tag[nt], "CellShapingSpeed");
      addr[nt] = &All.CellShapingSpeed;
      id[nt++] = REAL;
#endif /* #ifndef VORONOI_STATIC_MESH */

#if defined(COOLING)
      strcpy(tag[nt], "TreecoolFile");
      addr[nt] = &All.TreecoolFile;
      id[nt++] = STRING;
#endif /* #if defined(COOLING) */

#if defined(REFINEMENT)
      strcpy(tag[nt], "ReferenceGasPartMass");
      addr[nt] = &All.ReferenceGasPartMass;
      id[nt++] = REAL;

      strcpy(tag[nt], "TargetGasMassFactor");
      addr[nt] = &All.TargetGasMassFactor;
      id[nt++] = REAL;

      strcpy(tag[nt], "RefinementCriterion");
      addr[nt] = &All.RefinementCriterion;
      id[nt++] = INT;

      strcpy(tag[nt], "DerefinementCriterion");
      addr[nt] = &All.DerefinementCriterion;
      id[nt++] = INT;
#endif /* #if defined(REFINEMENT) */

#ifdef USE_SFR
      strcpy(tag[nt], "CritOverDensity");
      addr[nt] = &All.CritOverDensity;
      id[nt++] = REAL;

      strcpy(tag[nt], "TemperatureThresh");
      addr[nt] = &All.TemperatureThresh;
      id[nt++] = REAL;

      strcpy(tag[nt], "CritPhysDensity");
      addr[nt] = &All.CritPhysDensity;
      id[nt++] = REAL;

      strcpy(tag[nt], "FactorSN");
      addr[nt] = &All.FactorSN;
      id[nt++] = REAL;

      strcpy(tag[nt], "FactorEVP");
      addr[nt] = &All.FactorEVP;
      id[nt++] = REAL;

      strcpy(tag[nt], "TempSupernova");
      addr[nt] = &All.TempSupernova;
      id[nt++] = REAL;

      strcpy(tag[nt], "TempClouds");
      addr[nt] = &All.TempClouds;
      id[nt++] = REAL;

      strcpy(tag[nt], "MaxSfrTimescale");
      addr[nt] = &All.MaxSfrTimescale;
      id[nt++] = REAL;
#endif /* #ifdef USE_SFR */

#ifdef MHD_SEEDFIELD
      strcpy(tag[nt], "MHDSeedDir");
      addr[nt] = &All.B_dir;
      id[nt++] = INT;

      strcpy(tag[nt], "MHDSeedValue");
      addr[nt] = &All.B_value;
      id[nt++] = REAL;
#endif /* #ifdef MHD_SEEDFIELD */

#ifdef REFINEMENT_VOLUME_LIMIT
      strcpy(tag[nt], "MaxVolumeDiff");
      addr[nt] = &All.MaxVolumeDiff;
      id[nt++] = REAL;

      strcpy(tag[nt], "MinVolume");
      addr[nt] = &All.MinVolume;
      id[nt++] = REAL;

      strcpy(tag[nt], "MaxVolume");
      addr[nt] = &All.MaxVolume;
      id[nt++] = REAL;
#endif /* #ifdef REFINEMENT_VOLUME_LIMIT */

#ifdef TILE_ICS
      strcpy(tag[nt], "TileICsFactor");
      addr[nt] = &All.TileICsFactor;
      id[nt++] = INT;
#endif /* #ifdef TILE_ICS */

#ifdef ADDBACKGROUNDGRID
      strcpy(tag[nt], "GridSize");
      addr[nt] = &All.GridSize;
      id[nt++] = INT;
#endif /* #ifdef ADDBACKGROUNDGRID */

#ifdef ONEDIMS_SPHERICAL
      strcpy(tag[nt], "CoreRadius");
      addr[nt] = &All.CoreRadius;
      id[nt++] = REAL;

      strcpy(tag[nt], "CoreMass");
      addr[nt] = &All.CoreMass;
      id[nt++] = REAL;
#endif /* #ifdef ONEDIMS_SPHERICAL */

      if((fd = fopen(fname, "r")))
        {
          sprintf(buf, "%s%s", fname, "-usedvalues");
          if(!(fdout = fopen(buf, "w")))
            {
              printf("error opening file '%s' \n", buf);
              errorFlag = 1;
            }
          else
            {
              printf("Obtaining parameters from file '%s':\n\n", fname);
              while(!feof(fd))
                {
                  *buf = 0;
                  fgets(buf, MAXLEN_PARAM_TAG + MAXLEN_PARAM_VALUE + 200, fd);
                  if(sscanf(buf, "%s%s%s", buf1, buf2, buf3) < 2)
                    continue;

                  if(buf1[0] == '%')
                    continue;

                  for(i = 0, j = -1; i < nt; i++)
                    if(strcmp(buf1, tag[i]) == 0)
                      {
                        if(param_handled[i] == 0)
                          {
                            j                = i;
                            param_handled[i] = 1;
                            break;
                          }
                        else
                          {
                            j = -2;
                            break;
                          }
                      }

                  if(j >= 0)
                    {
                      switch(id[j])
                        {
                          case REAL:
                            *((double *)addr[j]) = atof(buf2);
                            sprintf(buf3, "%%-%ds%%g\n", MAXLEN_PARAM_TAG);
                            fprintf(fdout, buf3, buf1, *((double *)addr[j]));
                            fprintf(stdout, "        ");
                            fprintf(stdout, buf3, buf1, *((double *)addr[j]));
                            break;
                          case STRING:
                            strcpy((char *)addr[j], buf2);
                            sprintf(buf3, "%%-%ds%%s\n", MAXLEN_PARAM_TAG);
                            fprintf(fdout, buf3, buf1, buf2);
                            fprintf(stdout, "        ");
                            fprintf(stdout, buf3, buf1, buf2);
                            break;
                          case INT:
                            *((int *)addr[j]) = atoi(buf2);
                            sprintf(buf3, "%%-%ds%%d\n", MAXLEN_PARAM_TAG);
                            fprintf(fdout, buf3, buf1, *((int *)addr[j]));
                            fprintf(stdout, "        ");
                            fprintf(stdout, buf3, buf1, *((int *)addr[j]));
                            break;
                        }
                    }
                  else if(j == -2)
                    {
#ifdef ALLOWEXTRAPARAMS
                      warn("Tag '%s' ignored from file %s !", buf1, fname);
#else  /* #ifdef ALLOWEXTRAPARAMS */
                      fprintf(stdout, "Error in file %s:   Tag '%s' multiply defined.\n", fname, buf1);
                      errorFlag = 1;
#endif /* #ifdef ALLOWEXTRAPARAMS #else */
                    }
                  else
                    {
#ifdef ALLOWEXTRAPARAMS
                      warn("Tag '%s' ignored from file %s !", buf1, fname);
#else  /* #ifdef ALLOWEXTRAPARAMS */
                      fprintf(stdout, "Error in file %s:   Tag '%s' not allowed\n", fname, buf1);
                      errorFlag = 1;
#endif /* #ifdef ALLOWEXTRAPARAMS #else */
                    }
                }
              fclose(fd);
              fclose(fdout);
              printf("\n");

              i = strlen(All.OutputDir);
              if(i > 0)
                if(All.OutputDir[i - 1] != '/')
                  strcat(All.OutputDir, "/");

              mkdir(All.OutputDir, 02755);
              sprintf(buf1, "%s%s", fname, "-usedvalues");
              sprintf(buf2, "%s%s", All.OutputDir, "parameters-usedvalues");
              sprintf(buf3, "cp %s %s", buf1, buf2);
#ifndef NOCALLSOFSYSTEM
              if(errorFlag == 0)
                system(buf3);
#endif /* #ifndef NOCALLSOFSYSTEM */
            }
        }
      else
        {
          printf("Parameter file %s not found.\n", fname);
          errorFlag = 1;
        }

      for(i = 0; i < nt; i++)
        {
          if(param_handled[i] != 1)
            {
              printf("Error. I miss a value for tag '%s' in parameter file '%s'.\n", tag[i], fname);
              errorFlag = 1;
            }
        }

      if(All.OutputListOn && errorFlag == 0)
        errorFlag += read_outputlist(All.OutputListFilename);
      else
        All.OutputListLength = 0;
    }

  MPI_Bcast(&errorFlag, 1, MPI_INT, 0, MPI_COMM_WORLD);

  if(errorFlag)
    {
      MPI_Finalize();
      exit(errorFlag);
    }

  All.NParameters = nt;

  /* now communicate the relevant parameters to the other processes */
  MPI_Bcast(&All, sizeof(struct global_data_all_processes), MPI_BYTE, 0, MPI_COMM_WORLD);

#ifdef TOLERATE_WRITE_ERROR
  MPI_Bcast(AlternativeOutputDir, MAXLEN_PATH, MPI_BYTE, 0, MPI_COMM_WORLD);
#endif /* #ifdef TOLERATE_WRITE_ERROR */

#ifdef HOST_MEMORY_REPORTING
  check_maxmemsize_setting();
#endif /* #ifdef HOST_MEMORY_REPORTING */

  mymalloc_init();

  Parameters      = (char(*)[MAXLEN_PARAM_TAG])mymalloc("Parameters", All.NParameters * MAXLEN_PARAM_TAG * sizeof(char));
  ParametersValue = (char(*)[MAXLEN_PARAM_VALUE])mymalloc("ParametersValue", All.NParameters * MAXLEN_PARAM_VALUE * sizeof(char));
  ParametersType  = mymalloc("ParamtersType", All.NParameters * sizeof(char));

  if(ThisTask == 0)
    {
      for(i = 0; i < All.NParameters; i++)
        {
          strncpy(Parameters[i], tag[i], MAXLEN_PARAM_TAG);
          ParametersType[i] = id[i];
          void *tmp         = ParametersValue[i];
          switch(id[i])
            {
              case REAL:
                *((double *)tmp) = *((double *)addr[i]);
                break;
              case STRING:
                strncpy(tmp, addr[i], MAXLEN_PARAM_VALUE);
                break;
              case INT:
                tmp           = ParametersValue[i];
                *((int *)tmp) = *((int *)addr[i]);
                break;
            }
        }
    }

  MPI_Bcast(Parameters, sizeof(char) * All.NParameters * MAXLEN_PARAM_TAG, MPI_BYTE, 0, MPI_COMM_WORLD);
  MPI_Bcast(ParametersValue, sizeof(char) * All.NParameters * MAXLEN_PARAM_VALUE, MPI_BYTE, 0, MPI_COMM_WORLD);
  MPI_Bcast(ParametersType, sizeof(char) * All.NParameters, MPI_BYTE, 0, MPI_COMM_WORLD);

#undef REAL
#undef STRING
#undef INT
}

/*! \brief This function checks the consistency of the input parameters.
 *
 *  If you encounter some possible misuse and a corresponding error message
 *  that is hard to interpret, a check should be placed in this function with
 *  a terminate statement and a clear explanation why this does not work.
 *
 *  \return void
 */
void check_parameters()
{
  int i, errorFlag = 0;

  /* check whether time max is larger than max timestep */
  if(All.TimeMax - All.TimeBegin <= All.MaxSizeTimestep)
    {
      printf("PARAMETERS: check_parameters: TimeBegin = %g, TimeMax = %g, MaxSizeTimestep = %g \n", All.TimeBegin, All.TimeMax,
             All.MaxSizeTimestep);
      terminate(
          "check_parameters: Your total runtime is smaller than the maximum allowed timestep! Choose an appropriate value for "
          "MaxSizeTimestep < TimeMax-TimeBegin! \n");
    }

  /* check softening types */
  for(i = 0; i < NTYPES; i++)
    {
      if(All.SofteningTypeOfPartType[i] >= NSOFTTYPES || All.SofteningTypeOfPartType[i] < 0)
        {
          mpi_printf("SofteningTypeOfPartType%  invalid (NSOFTTYPES=%d)\n", i, NSOFTTYPES);
          errorFlag = 1;
        }
    }

  if(errorFlag)
    mpi_terminate("Softening invalid!");

  if(All.NumFilesWrittenInParallel > NTask)
    {
      if(ThisTask == 0)
        warn("NOTICE: Reducing requested NumFilesWrittenInParallel=%d to %d\n", All.NumFilesWrittenInParallel, NTask);
      All.NumFilesWrittenInParallel = NTask;
    }

  if(All.NumFilesWrittenInParallel == 0)
    {
      mpi_printf("NOTICE: All.NumFilesWrittenInParallel has been set to be equal to the number of processors\n");
      All.NumFilesWrittenInParallel = NTask;
    }

#ifndef GRAVITY_NOT_PERIODIC
  if(All.PeriodicBoundariesOn == 0)
    {
      mpi_terminate(
          "Code was compiled with gravity periodic boundary conditions switched on.\nYou must set `PeriodicBoundariesOn=1', or "
          "recompile the code.\n");
    }
#else  /* #ifndef GRAVITY_NOT_PERIODIC */
  if(All.PeriodicBoundariesOn == 1)
    {
      mpi_terminate(
          "Code was compiled with gravity periodic boundary conditions switched off.\nYou must set `PeriodicBoundariesOn=0', or "
          "recompile the code.\n");
    }
#endif /* #ifndef GRAVITY_NOT_PERIODIC #else */

#ifdef COOLING
  if(All.CoolingOn == 0)
    {
      mpi_terminate("Code was compiled with cooling switched on.\nYou must set `CoolingOn=1', or recompile the code.\n");
    }
#else  /* #ifdef COOLING */
  if(All.CoolingOn == 1)
    {
      mpi_terminate("Code was compiled with cooling switched off.\nYou must set `CoolingOn=0', or recompile the code.\n");
    }
#endif /* #ifdef COOLING #else */

  if(All.TypeOfTimestepCriterion >= 3)
    {
      mpi_terminate("The specified timestep criterion\nis not valid\n");
    }

#if(NTYPES < 6)
  mpi_terminate("NTYPES < 6 is not allowed.\n");
#endif /* #if (NTYPES < 6) */

#if(NTYPES > 15)
  mpi_terminate("NTYPES > 15 is not supported yet.\n");
#endif /* #if (NTYPES > 15) */

#if(NTYPES > 8)
  if(All.ICFormat == 1 || All.ICFormat == 2)
    {
      mpi_terminate("NTYPES>8 is not allowed with ICFormat=%d, since the header block is limited to 256 bytes.\n", All.ICFormat);
    }
#endif /* #if (NTYPES > 8) */

#ifdef USE_SFR
  if(All.StarformationOn == 0)
    {
      mpi_terminate("Code was compiled with star formation switched on.\nYou must set `StarformationOn=1', or recompile the code.\n");
    }
  if(All.CoolingOn == 0)
    {
      mpi_terminate(
          "You try to use the code with star formation enabled,\nbut you did not switch on cooling.\nThis mode is not supported.\n");
    }
#else  /* #ifdef USE_SFR */
  if(All.StarformationOn == 1)
    {
      mpi_terminate("Code was compiled with star formation switched off.\nYou must set `StarformationOn=0', or recompile the code.\n");
    }
#endif /* #ifdef USE_SFR #else */

#if defined(ENFORCE_JEANS_STABILITY_OF_CELLS) && defined(USE_SFR)
  if(ThisTask == 0)
    warn("Code was compiled with ENFORCE_JEANS_STABILITY_OF_CELLS together with another EOS. Please make sure you really want this.");
#endif /* #if defined(ENFORCE_JEANS_STABILITY_OF_CELLS) && (defined(ISOTHERM_EQS) || (defined(USE_SFR) && !defined(FM_SFR))) */
}

/*! \brief This function reads a table with a list of desired output times.
 *
 *  The table does not have to be ordered in any way, but may not contain more
 *  than MAXLEN_OUTPUTLIST entries.
 *
 *  \param[in] fname The file name of the outputlist.
 *
 *  \return 0: success  1: unable to open file.
 */
int read_outputlist(char *fname)
{
  FILE *fd;
  int count, flag;
  char buf[512], msg[512];

  if(!(fd = fopen(fname, "r")))
    {
      printf("can't read output list in file '%s'\n", fname);
      return 1;
    }

  All.OutputListLength = 0;

  while(1)
    {
      if(fgets(buf, 500, fd) != buf)
        break;

      count = sscanf(buf, " %lg %d ", &All.OutputListTimes[All.OutputListLength], &flag);

      if(count == 1)
        flag = 1;

      if(count == 1 || count == 2)
        {
          if(All.OutputListLength >= MAXLEN_OUTPUTLIST)
            {
              sprintf(msg, "\ntoo many entries in output-list. You should increase MAXLEN_OUTPUTLIST=%d.\n", (int)MAXLEN_OUTPUTLIST);
              terminate(msg);
            }

          All.OutputListFlag[All.OutputListLength] = flag;
          All.OutputListLength++;
        }
    }

  fclose(fd);

  printf("\nBEGRUN: found %d times in output-list.\n", All.OutputListLength);

  return 0;
}
