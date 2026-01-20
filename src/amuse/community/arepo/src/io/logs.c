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
 * \file        src/io/logs.c
 * \date        05/2018
 * \brief       Log-files handling.
 * \details     contains functions:
 *                void open_logfiles(void)
 *                void close_logfiles(void)
 *                void output_log_messages(void)
 *                void init_cpu_log(void)
 *                void write_cpu_log(void)
 *                void put_symbol(char *string, double t0, double t1, char c)
 *
 *
 * \par Major modifications and contributions:
 *
 * - DD.MM.YYYY Description
 * - 07.05.2018 Prepared file for public release -- Rainer Weinberger
 */

#include <ctype.h>
#include <math.h>
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>
#include <unistd.h>

#include "../main/allvars.h"
#include "../main/proto.h"

#include "../mesh/voronoi/voronoi.h"

#define CPU_STRING_LEN 120

/*! \brief Contains informations about the used CPU timers like it's name,
 * symbols etc.
 */
struct timer_d Timer_data[CPU_LAST + 1];

enum timers TimerStack[TIMER_STACK_DEPTH];
int TimerStackPos = 0;

/*! \brief Opens files for logging.
 *
 *   This function opens various log-files that report on the status and
 *   performance of the simulation. Upon restart, the code will append to
 *   these files.
 *
 *   \return void
 */
void open_logfiles(void)
{
  char mode[2], buf[1000], msg[1000];

  if(RestartFlag == 0)
    strcpy(mode, "w");
  else
    strcpy(mode, "a");

  if(ThisTask == 0)
    mkdir(All.OutputDir, 02755);

  MPI_Barrier(MPI_COMM_WORLD);

#ifdef DETAILEDTIMINGS
  sprintf(buf, "%stimings_detailed_%d.txt", All.OutputDir, ThisTask);
  if(!(FdDetailed = fopen(buf, mode)))
    terminate("error in opening file '%s'\n", buf);
#endif /* #ifdef DETAILEDTIMINGS */

  if(ThisTask != 0) /* only the root processors writes to the log files */
    return;

  sprintf(buf, "%s%s", All.OutputDir, "cpu.txt");
  if(!(FdCPU = fopen(buf, mode)))
    {
      sprintf(msg, "error in opening file '%s'\n", buf);
      terminate(msg);
    }

  sprintf(buf, "%s%s", All.OutputDir, "info.txt");
  if(!(FdInfo = fopen(buf, mode)))
    {
      sprintf(msg, "error in opening file '%s'\n", buf);
      terminate(msg);
    }

  sprintf(buf, "%s%s", All.OutputDir, "energy.txt");
  if(!(FdEnergy = fopen(buf, mode)))
    {
      sprintf(msg, "error in opening file '%s'\n", buf);
      terminate(msg);
    }

  sprintf(buf, "%s%s", All.OutputDir, "timings.txt");
  if(!(FdTimings = fopen(buf, mode)))
    {
      sprintf(msg, "error in opening file '%s'\n", buf);
      terminate(msg);
    }

  sprintf(buf, "%s%s", All.OutputDir, "balance.txt");
  if(!(FdBalance = fopen(buf, mode)))
    {
      sprintf(msg, "error in opening file '%s'\n", buf);
      terminate(msg);
    }

  sprintf(buf, "%s%s", All.OutputDir, "timebins.txt");
  if(!(FdTimebin = fopen(buf, mode)))
    {
      sprintf(msg, "error in opening file '%s'\n", buf);
      terminate(msg);
    }

  sprintf(buf, "%s%s", All.OutputDir, "domain.txt");
  if(!(FdDomain = fopen(buf, mode)))
    {
      sprintf(msg, "error in opening file '%s'\n", buf);
      terminate(msg);
    }

  sprintf(buf, "%s%s", All.OutputDir, "memory.txt");
  if(!(FdMemory = fopen(buf, mode)))
    {
      sprintf(msg, "error in opening file '%s'\n", buf);
      terminate(msg);
    }

#ifdef FORCETEST
  sprintf(buf, "%s%s", All.OutputDir, "forcetest.txt");
  if(!(FdForceTest = fopen(buf, mode)))
    {
      sprintf(msg, "error in opening file '%s'\n", buf);
      terminate(msg);
    }
  fclose(FdForceTest);
#endif /* #ifdef FORCETEST */

#ifdef RESTART_DEBUG
  sprintf(buf, "%s%s", All.OutputDir, "restartdebug.txt");
  if(!(FdRestartTest = fopen(buf, mode)))
    {
      sprintf(msg, "error in opening file '%s'\n", buf);
      terminate(msg);
    }
#endif /* #ifdef RESTART_DEBUG */

#ifdef OUTPUT_CPU_CSV
  sprintf(buf, "%s%s", All.OutputDir, "cpu.csv");
  if(!(FdCPUCSV = fopen(buf, mode)))
    {
      sprintf(msg, "error in opening file '%s'\n", buf);
      terminate(msg);
    }
#endif /* #ifdef OUTPUT_CPU_CSV */

#ifdef USE_SFR
  sprintf(buf, "%s%s", All.OutputDir, "sfr.txt");
  if(!(FdSfr = fopen(buf, mode)))
    {
      sprintf(msg, "error in opening file '%s'\n", buf);
      terminate(msg);
    }
#endif /* #ifdef USE_SFR */

  int i = 0;
  fprintf(FdBalance, "\n");

#ifdef OUTPUT_CPU_CSV
  fprintf(FdCPUCSV, "STEP, TIME, CPUS, MULTIPLEDOMAIN, HIGHESTTIMEBIN, ");
#endif /* #ifdef OUTPUT_CPU_CSV */
  for(; i < CPU_LAST; i++)
    {
      if(Timer_data[i].symb != 0 && Timer_data[i].symbImbal != 0)
        {
          fprintf(FdBalance, "%-20s = '%c' / '%c'\n", Timer_data[i].longname, Timer_data[i].symb, Timer_data[i].symbImbal);
        }
#ifdef OUTPUT_CPU_CSV
      fprintf(FdCPUCSV, "%s1, %s2, %s3, ", Timer_data[i].shortname, Timer_data[i].shortname, Timer_data[i].shortname);
#endif /* #ifdef OUTPUT_CPU_CSV */
    }
  fprintf(FdBalance, "\n");

#ifdef OUTPUT_CPU_CSV
  fprintf(FdCPUCSV, "\n");
#endif /* #ifdef OUTPUT_CPU_CSV */
}

/*! \brief Closes the global log-files.
 *
 *  \return void
 */
void close_logfiles(void)
{
  if(ThisTask != 0) /* only the root processors writes to the log files */
    return;

  fclose(FdCPU);
  fclose(FdInfo);
  fclose(FdEnergy);
  fclose(FdTimings);
  fclose(FdBalance);
  fclose(FdTimebin);

#ifdef OUTPUT_CPU_CSV
  fclose(FdCPUCSV);
#endif /* #ifdef OUTPUT_CPU_CSV */

#ifdef USE_SFR
  fclose(FdSfr);
#endif /* #ifdef USE_SFR */
}

/*! \brief Writes log messages in log-files.
 *
 *  At each time step this function writes on to two log-files.
 *  In FdInfo, it just lists the timesteps that have been done, while in
 *  FdTimeBin it outputs information about the active and occupied time-bins.
 *  Additionally, reports to memory log-files are written.
 *
 *  \return void
 */
void output_log_messages(void)
{
  double z;
  int i, j, write_logs = 1;
  double sum, avg_CPU_TimeBin[TIMEBINS], frac_CPU_TimeBin[TIMEBINS];
  int weight, corr_weight;
  long long tot_cumulative_grav[TIMEBINS], tot_cumulative_sph[TIMEBINS];
  long long tot_grav, tot_sph;

  TIMER_START(CPU_LOGS);

  if(write_logs)
    report_detailed_memory_usage_of_largest_task();

  long long count[4 * TIMEBINS], tot_count[4 * TIMEBINS];
  long long *tot_count_grav = &tot_count[0], *tot_count_sph = &tot_count[TIMEBINS];
  int nelem = 2 * TIMEBINS;

  for(int i = 0; i < TIMEBINS; i++)
    count[i] = TimeBinsGravity.TimeBinCount[i];

  for(int i = 0; i < TIMEBINS; i++)
    count[i + TIMEBINS] = TimeBinsHydro.TimeBinCount[i];

  MPI_Reduce(count, tot_count, nelem, MPI_LONG_LONG_INT, MPI_SUM, 0, MPI_COMM_WORLD);

  if(ThisTask == 0)
    {
      if(All.ComovingIntegrationOn)
        {
          z = 1.0 / (All.Time) - 1;

          if(write_logs)
            fprintf(FdInfo,
                    "\nSync-Point %d, TimeBin=%d, Time: %g, Redshift: %g, Systemstep: %g, Dloga: %g, Nsync-grv: %10llu, Nsync-hyd: "
                    "%10llu\n",
                    All.NumCurrentTiStep, All.HighestActiveTimeBin, All.Time, z, All.TimeStep,
                    log(All.Time) - log(All.Time - All.TimeStep), All.GlobalNSynchronizedGravity, All.GlobalNSynchronizedHydro);

          printf("\n\nSync-Point %d, Time: %g, Redshift: %g, Systemstep: %g, Dloga: %g, Nsync-grv: %10llu, Nsync-hyd: %10llu\n",
                 All.NumCurrentTiStep, All.Time, z, All.TimeStep, log(All.Time) - log(All.Time - All.TimeStep),
                 All.GlobalNSynchronizedGravity, All.GlobalNSynchronizedHydro);

          if(write_logs)
            fprintf(FdTimebin, "\nSync-Point %d, Time: %g, Redshift: %g, Systemstep: %g, Dloga: %g\n", All.NumCurrentTiStep, All.Time,
                    z, All.TimeStep, log(All.Time) - log(All.Time - All.TimeStep));

          myflush(FdInfo);
        }
      else
        {
          if(write_logs)
            fprintf(FdInfo, "\nSync-Point %d, TimeBin=%d, Time: %g, Systemstep: %g, Nsync-grv: %10llu, Nsync-hyd: %10llu\n",
                    All.NumCurrentTiStep, All.HighestActiveTimeBin, All.Time, All.TimeStep, All.GlobalNSynchronizedGravity,
                    All.GlobalNSynchronizedHydro);

          printf("\n\nSync-Point %d, Time: %g, Systemstep: %g, Nsync-grv: %10llu, Nsync-hyd: %10llu\n", All.NumCurrentTiStep, All.Time,
                 All.TimeStep, All.GlobalNSynchronizedGravity, All.GlobalNSynchronizedHydro);

          if(write_logs)
            fprintf(FdTimebin, "\nSync-Point %d, Time: %g, Systemstep: %g\n", All.NumCurrentTiStep, All.Time, All.TimeStep);

          myflush(FdInfo);
        }

      for(i = 1, tot_cumulative_grav[0] = tot_count_grav[0], tot_cumulative_sph[0] = tot_count_sph[0]; i < TIMEBINS; i++)
        {
          tot_cumulative_grav[i] = tot_count_grav[i] + tot_cumulative_grav[i - 1];
          tot_cumulative_sph[i]  = tot_count_sph[i] + tot_cumulative_sph[i - 1];
        }

      for(i = 0; i < TIMEBINS; i++)
        {
          for(j = 0, sum = 0; j < All.CPU_TimeBinCountMeasurements[i]; j++)
            sum += All.CPU_TimeBinMeasurements[i][j];
          if(All.CPU_TimeBinCountMeasurements[i])
            avg_CPU_TimeBin[i] = sum / All.CPU_TimeBinCountMeasurements[i];
          else
            avg_CPU_TimeBin[i] = 0;
        }

      for(i = All.HighestOccupiedTimeBin, weight = 1, sum = 0; i >= 0 && tot_count_grav[i] > 0; i--, weight *= 2)
        {
          if(weight > 1)
            corr_weight = weight / 2;
          else
            corr_weight = weight;

          frac_CPU_TimeBin[i] = corr_weight * avg_CPU_TimeBin[i];
          sum += frac_CPU_TimeBin[i];
        }

      for(i = All.HighestOccupiedTimeBin; i >= 0 && tot_count_grav[i] > 0; i--)
        {
          if(sum)
            frac_CPU_TimeBin[i] /= sum;
        }

      char tracerString[13];

      sprintf(tracerString, "%s", "");

      char dustString[13];
      sprintf(dustString, "%s", "");
      if(write_logs)
        fprintf(FdTimebin,
                "Occupied timebins: gravity      hydro     %s     %s     dt              cumul-grav   cumul-sph A D    avg-time  "
                "cpu-frac\n",
                tracerString, dustString);

      for(i = TIMEBINS - 1, tot_grav = tot_sph = 0; i >= 0; i--)
        {
          int binUsed = 0;

#if(defined(SELFGRAVITY) || defined(EXTERNALGRAVITY) || defined(EXACT_GRAVITY_FOR_PARTICLE_TYPE)) && !defined(MESHRELAX)
          if(tot_count_grav[i] > 0)
            binUsed = 1;
#endif /* #if (defined(SELFGRAVITY) || defined(EXTERNALGRAVITY) || defined(EXACT_GRAVITY_FOR_PARTICLE_TYPE)) && !defined(MESHRELAX) \
        */

          if(tot_count_sph[i] > 0)
            binUsed = 1;

          sprintf(tracerString, "%s", "");

          if(binUsed)
            {
              if(write_logs)
                fprintf(FdTimebin, " %c  bin=%2d      %10llu  %10llu  %s  %s  %16.12f       %10llu  %10llu %c %c  %10.2f    %5.1f%%\n",
                        TimeBinSynchronized[i] ? 'X' : ' ', i, tot_count_grav[i], tot_count_sph[i], tracerString, dustString,
                        i > 0 ? (((integertime)1) << i) * All.Timebase_interval : 0.0, tot_cumulative_grav[i], tot_cumulative_sph[i],
                        (i == All.HighestActiveTimeBin) ? '<' : ' ',
                        (All.HighestActiveTimeBin >= All.SmallestTimeBinWithDomainDecomposition && i == All.HighestActiveTimeBin)
                            ? '*'
                            : ' ',
                        avg_CPU_TimeBin[i], 100.0 * frac_CPU_TimeBin[i]);

              if(TimeBinSynchronized[i])
                {
                  tot_grav += tot_count_grav[i];
                  tot_sph += tot_count_sph[i];
                }
            }
        }

      if(write_logs)
        {
          fprintf(FdTimebin, "               ------------------------\n");
        }

      sprintf(tracerString, "%s", "");
      sprintf(dustString, "%s", "");

      if(write_logs)
        {
#ifdef PMGRID
          if(All.PM_Ti_endstep == All.Ti_Current)
            {
              fprintf(FdTimebin, "PM-Step. Total: %10llu  %10llu  %s  %s\n", tot_grav, tot_sph, tracerString, dustString);
            }
          else
#endif /* #ifdef PMGRID */
            {
              fprintf(FdTimebin, "Total active:   %10llu  %10llu  %s  %s\n", tot_grav, tot_sph, tracerString, dustString);
            }

          fprintf(FdTimebin, "\n");
        }

      myflush(FdTimebin);
    }

#ifdef RESTART_DEBUG
  log_restart_debug();
#endif /* #ifdef RESTART_DEBUG */

  TIMER_STOP(CPU_LOGS);
}

/*! \brief Initializes cpu log file.
 *
 *  \return void
 */
void init_cpu_log(void)
{
  int i = 0;

#define TIMER_STRUCT
#include "../utils/timer.h"

  for(i = 0; i < CPU_LAST; i++)
    {
      if(Timer_data[i].parent >= 0)
        Timer_data[i].depth = Timer_data[Timer_data[i].parent].depth + 1;
      else
        Timer_data[i].depth = 0;
    }

  for(i = 0; i < CPU_LAST; i++)
    {
      All.CPU_Sum[i] = 0.;
      CPU_Step[i]    = 0.;
    }

  TimerStackPos = 0;
  TimerStack[0] = CPU_MISC;

  CPUThisRun = 0.;

  WallclockTime = second();
  StartOfRun    = second();
}

/*! \brief Write the FdBalance and FdCPU files.
 *
 *  At each time step this function writes on to two log-files.
 *  In FdBalance, it outputs in a graphical way the amount of
 *  time spent in the various parts of the code, while
 *  in FdCPU it writes information about the cpu-time consumption
 *  of the various modules.
 *
 * \return void
 */
void write_cpu_log(void)
{
  int write_logs = 1;
  double max_CPU_Step[CPU_LAST], avg_CPU_Step[CPU_LAST], summed_CPU_Step[CPU_LAST];
  double t0, t1, tsum;
  double avg_total   = 0;
  double local_total = 0;
  double max_total   = 0;
  int i;

  TIMER_START(CPU_LOGS);

  for(i = 0; i < CPU_LAST; i++)
    {
      local_total += CPU_Step[i];
    }

  MPI_Reduce(CPU_Step, max_CPU_Step, CPU_LAST, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
  MPI_Reduce(&local_total, &max_total, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
  MPI_Reduce(CPU_Step, avg_CPU_Step, CPU_LAST, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

  if(ThisTask == 0)
    {
      /* sum up cpu items into groups */
      for(i = 0; i < CPU_LAST; i++)
        {
          summed_CPU_Step[i] = avg_CPU_Step[i];
        }
      for(i = CPU_LAST - 1; i > CPU_ALL; i--)
        {
          if(Timer_data[i].parent >= 0)
            {
              summed_CPU_Step[Timer_data[i].parent] += summed_CPU_Step[i];
            }
        }

      /* calc averages, update All.CPU_Sum */
      for(i = 0; i < CPU_LAST; i++)
        {
          avg_CPU_Step[i] /= NTask;
          avg_total += avg_CPU_Step[i];

          summed_CPU_Step[i] /= NTask;
          All.CPU_Sum[i] += summed_CPU_Step[i];
        }

      /* create balance.txt string */
      char cpu_String[CPU_STRING_LEN + 1];
      put_symbol(cpu_String, 0., 1.0, '-');

      for(i = 1, tsum = 0.0; i < CPU_LAST; i++)
        {
          if(max_CPU_Step[i] > 0 && Timer_data[i].symb != 0 && Timer_data[i].symbImbal != 0)
            {
              t0 = tsum;
              t1 = tsum + avg_CPU_Step[i] * (avg_CPU_Step[i] / max_CPU_Step[i]);
              put_symbol(cpu_String, t0 / avg_total, t1 / avg_total, Timer_data[i].symb);
              tsum += t1 - t0;

              t0 = tsum;
              t1 = tsum + avg_CPU_Step[i] * ((max_CPU_Step[i] - avg_CPU_Step[i]) / max_CPU_Step[i]);
              put_symbol(cpu_String, t0 / avg_total, t1 / avg_total, Timer_data[i].symbImbal);
              tsum += t1 - t0;
            }
        }

      if(write_logs)
        {
          fprintf(FdBalance, "Step=%7d  sec=%10.3f Nsync-grv=%10llu Nsync-hyd=%10llu  %s\n", All.NumCurrentTiStep, max_total,
                  All.GlobalNSynchronizedGravity, All.GlobalNSynchronizedHydro, cpu_String);
        }

      myflush(FdBalance);

      if(All.CPU_TimeBinCountMeasurements[All.HighestActiveTimeBin] == NUMBER_OF_MEASUREMENTS_TO_RECORD)
        {
          All.CPU_TimeBinCountMeasurements[All.HighestActiveTimeBin]--;
          memmove(&All.CPU_TimeBinMeasurements[All.HighestActiveTimeBin][0], &All.CPU_TimeBinMeasurements[All.HighestActiveTimeBin][1],
                  (NUMBER_OF_MEASUREMENTS_TO_RECORD - 1) * sizeof(double));
        }

      All.CPU_TimeBinMeasurements[All.HighestActiveTimeBin][All.CPU_TimeBinCountMeasurements[All.HighestActiveTimeBin]++] = max_total;

      if(write_logs)
        {
#ifdef OUTPUT_CPU_CSV
          fprintf(FdCPUCSV, "%d, %g, %d, %d, %d, ", All.NumCurrentTiStep, All.Time, NTask, All.MultipleDomains,
                  All.HighestActiveTimeBin);
#endif /* #ifdef OUTPUT_CPU_CSV */
          fprintf(FdCPU, "Step %d, Time: %g, CPUs: %d, MultiDomains: %d, HighestActiveTimeBin: %d\n", All.NumCurrentTiStep, All.Time,
                  NTask, All.MultipleDomains, All.HighestActiveTimeBin);

          fprintf(FdCPU, "                          diff               cumulative\n");

          for(i = 0; i < CPU_LAST; i++)
            {
              fprintf(FdCPU, "%*s%*s%10.2f  %5.1f%% %10.2f  %*s%5.1f%%\n", 2 * Timer_data[i].depth, "", -20 + 2 * Timer_data[i].depth,
                      Timer_data[i].longname, summed_CPU_Step[i], summed_CPU_Step[i] / summed_CPU_Step[CPU_ALL] * 100., All.CPU_Sum[i],
                      5 * Timer_data[i].depth, "", All.CPU_Sum[i] / All.CPU_Sum[CPU_ALL] * 100.);

#ifdef OUTPUT_CPU_CSV
              fprintf(FdCPUCSV, "%f, %f, %f, ", summed_CPU_Step[i], All.CPU_Sum[i], All.CPU_Sum[i] / All.CPU_Sum[CPU_ALL] * 100.);
#endif /* #ifdef OUTPUT_CPU_CSV */
            }

          fprintf(FdCPU, "\n");
        }

      myflush(FdCPU);

#ifdef OUTPUT_CPU_CSV
      if(write_logs)
        fprintf(FdCPUCSV, "\n");

      myflush(FdCPUCSV);
#endif /* #ifdef OUTPUT_CPU_CSV */
    }

  for(i = 0; i < CPU_LAST; i++)
    CPU_Step[i] = 0.;

  CPUThisRun = timediff(StartOfRun, second());

  TIMER_STOP(CPU_LOGS);
}

/*! \brief Fill the cpu balance string representing the cpu usage in a
 *         graphical way.
 *
 *  This function fills a fraction, specified by the parameters t0 and t1,
 *  of the array string with the debug symbol given by c.
 *
 *  \param[out] string String to fill.
 *  \param[in] t0 Initial position of the symbol in the array as a fraction of
 *             its maximum dimension.
 *  \param[in] t1 Final position of the symbol in the array as a fraction of
 *             its maximum dimension.
 *  \param[in] c Symbol to be put on string.
 *
 *  \return void
 */
void put_symbol(char *string, double t0, double t1, char c)
{
  int i, j;

  i = (int)(t0 * CPU_STRING_LEN + 0.5);
  j = (int)(t1 * CPU_STRING_LEN);

  if(i < 0)
    i = 0;
  if(j >= CPU_STRING_LEN)
    j = CPU_STRING_LEN;

  while(i <= j)
    string[i++] = c;

  string[CPU_STRING_LEN] = 0;
}
