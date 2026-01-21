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
 * \file        src/utils/timer.h
 * \date        05/2018
 * \brief       Timer macros for Arepo.
 * \details
 *
 * \par Major modifications and contributions:
 *
 * - DD.MM.YYYY Description
 * - 28.05.2018 Prepared file for public release -- Rainer Weinberger
 */

#if !defined(TIMER_H) || defined(TIMER_STRUCT)
#define TIMER_H

#define DETAILED_TIMING_GRAVWALK 0
#define DETAILED_TIMING_STELLARDENSITY 1

#define TIMER_INSTRUMENT_START(counter)
#define TIMER_INSTRUMENT_STOP(counter)
#define TIMER_INSTRUMENT_CREATE(name, descr) ;

#ifdef TIMER_STRUCT
#undef TIMER_CREATE
/*! \def TIMER_CREATE(name,desc, par, symba, symbb )
 *  \brief creates a new CPU timer
 *
 *  \param name name used in the code to reference this timer
 *  \param desc description string used in output files
 *  \param parent parent of this timer to build a tree-like hierarchy of timers
 *  \param symba character used for active time in balance.txt
 *  \param symbb character used for imbalance in balance.txt
 *
 */
#define TIMER_CREATE(name, desc, par, symba, symbb) \
  Timer_data[name].parent = par;                    \
  strncpy(Timer_data[name].shortname, #name, 40);   \
  strncpy(Timer_data[name].longname, (desc), 40);   \
  Timer_data[name].symb      = (symba);             \
  Timer_data[name].symbImbal = (symbb);             \
  TIMER_INSTRUMENT_CREATE(name, desc)

#else /* #ifdef TIMER_STRUCT */

#define TIMER_STACK_DEPTH 30
#define TIMER_CREATE(name, desc, parent, symba, symbb) name,

/*! \def  TIMER_START(counter)
 *  \brief Starts the timer counter.
 *
 *  Use this macro instead of directly accessing the CPU_Step array,
 *  so manual  instrumentation APIs can be attached.
 *
 *  \param[in] counter Name of the timer to start.
 */
#define TIMER_START_INTERNAL(counter)                                             \
  {                                                                               \
    TIMER_INSTRUMENT_START(counter);                                              \
    CPU_Step[TimerStack[TimerStackPos]] += measure_time();                        \
    int itimer;                                                                   \
    for(itimer = 0; itimer <= TimerStackPos; itimer++)                            \
      if(counter == TimerStack[itimer])                                           \
        {                                                                         \
          printf("Try to start timer %d, but it is already running.\n", counter); \
          terminate("fail")                                                       \
        };                                                                        \
    if(++TimerStackPos >= TIMER_STACK_DEPTH)                                      \
      {                                                                           \
        terminate("Run out of timer stack space, increase TIMER_STACK_DEPTH");    \
      }                                                                           \
    else                                                                          \
      {                                                                           \
        TimerStack[TimerStackPos] = (counter);                                    \
      }                                                                           \
  }

#define TIMER_START(counter) TIMER_START_INTERNAL(counter)

/*! \def TIMER_STOP(counter)
 *  \brief Stops the timer counter
 *
 *  Use this macro instead of directly accessing the CPU_Step array,
 *  so manual instrumentation APIs can be attached.
 *
 *  \param[in] counter Name of the timer to stop.
 */
#define TIMER_STOP_INTERNAL(counter)                                                \
  {                                                                                 \
    if(TimerStack[TimerStackPos] != (counter))                                      \
      {                                                                             \
        terminate("Wrong use of TIMER_STOP, you must stop the timer started last"); \
      }                                                                             \
    CPU_Step[TimerStack[TimerStackPos--]] += measure_time();                        \
    if(TimerStackPos < 0)                                                           \
      {                                                                             \
        terminate("Do not stop the out CPU_MISC timer");                            \
      }                                                                             \
    TIMER_INSTRUMENT_STOP(counter);                                                 \
  }

#define TIMER_STOP(counter) TIMER_STOP_INTERNAL(counter)

/*! \def TIMER_STOPSTART(stop, start)
 *  \brief Stops the timer 'stop' and starts the timer 'start'
 *
 *  Use this macro instead of directly accessing the CPU_Step array,
 *  so manual instrumentation APIs can be attached.
 *
 *  \param[in] stop Name of the timer to stop
 *  \param[in] start Name of the timer to start
 */
#define TIMER_STOPSTART(stop, start) \
  {                                  \
    TIMER_STOP_INTERNAL(stop);       \
    TIMER_START_INTERNAL(start);     \
  }

/*! \def TIMER_ADD(counter, amount)
 *  \brief Adds amount to the timer counter.

 *  \param[in] counter Name of the timer to add to.
 *  \param[in] amount Amount to add to timer counter.
 */
#define TIMER_ADD(counter, amount) CPU_Step[counter] += (amount);

/*! \def TIMER_DIFF(counter)
 *  \brief Returns amount elapsed for the timer since last save with
 *         TIMER_STORE.
 *
 *  \param[in] counter Name of the timer to add to.
 */
#define TIMER_DIFF(counter) (CPU_Step[counter] - CPU_Step_Stored[counter])

/*! \def TIMER_STORE
 *  \brief Copies the current value of CPU times to a stored variable, such
 *         that differences with respect to this reference can be calculated.
 */
#define TIMER_STORE memcpy(CPU_Step_Stored, CPU_Step, sizeof(CPU_Step));

enum timers
{
  CPU_NONE = -2,                /*!< used for counters without a parent */
  CPU_ROOT = -1,                /*!< root node of the tree */
#endif /* #ifdef TIMER_STRUCT #else */

/* possible characters to use for marking the parts:
 *
 *   abdefghijklmnopqrstuvABCDEFGHHIJKLMNOPQRSTUV
 *   0123456789
 *   -:.*=[]^&;~/_$()?+"<>@#!|\
 */

/*add your counter here, they must appear in the right order*/

TIMER_CREATE(CPU_ALL, "total", CPU_ROOT, '-', '-') /*!< root timer, everything should be below this timer */
TIMER_CREATE(CPU_TREE, "treegrav", CPU_ALL, 'a', ')')
TIMER_CREATE(CPU_TREEBUILD, "treebuild", CPU_TREE, 'b', '(')
TIMER_CREATE(CPU_TREEBUILD_INSERT, "insert", CPU_TREEBUILD, 'c', '*')
TIMER_CREATE(CPU_TREEBUILD_BRANCHES, "branches", CPU_TREEBUILD, 'd', '&')
TIMER_CREATE(CPU_TREEBUILD_TOPLEVEL, "toplevel", CPU_TREEBUILD, 'e', '^')
TIMER_CREATE(CPU_TREECOSTMEASURE, "treecostm", CPU_TREE, 'f', '%')
TIMER_CREATE(CPU_TREEWALK, "treewalk", CPU_TREE, 'g', '$')
TIMER_CREATE(CPU_TREEWALK1, "treewalk1", CPU_TREEWALK, 'h', '#')
TIMER_CREATE(CPU_TREEWALK2, "treewalk2", CPU_TREEWALK, 'i', '@')
TIMER_CREATE(CPU_TREEBALSNDRCV, "treebalsndrcv", CPU_TREE, 'j', '!')
TIMER_CREATE(CPU_TREESENDBACK, "treeback", CPU_TREE, 'm', '7')
TIMER_CREATE(CPU_TREEDIRECT, "treedirect", CPU_TREE, 'r', '2')
#ifdef PMGRID
TIMER_CREATE(CPU_PM_GRAVITY, "pm_grav", CPU_ALL, 's', '1')
#endif /* #ifdef PMGRID */
TIMER_CREATE(CPU_NGBTREEBUILD, "ngbtreebuild", CPU_ALL, 't', 'Z')
TIMER_CREATE(CPU_NGBTREEUPDATEVEL, "ngbtreevelupdate", CPU_ALL, 'u', 'Y')
TIMER_CREATE(CPU_MESH, "voronoi", CPU_ALL, 'v', 'X')
TIMER_CREATE(CPU_MESH_INSERT, "insert", CPU_MESH, 'w', 'W')
TIMER_CREATE(CPU_MESH_FIND_DP, "findpoints", CPU_MESH, 'x', 'V')
TIMER_CREATE(CPU_MESH_CELLCHECK, "cellcheck", CPU_MESH, 'y', 'U')
TIMER_CREATE(CPU_MESH_GEOMETRY, "geometry", CPU_MESH, 'z', 'T')
TIMER_CREATE(CPU_MESH_EXCHANGE, "exchange", CPU_MESH, 'A', 'S')
TIMER_CREATE(CPU_MESH_DYNAMIC, "dynamic", CPU_MESH, 'B', 'R')
TIMER_CREATE(CPU_HYDRO, "hydro", CPU_ALL, 'C', 'Q')
TIMER_CREATE(CPU_GRADIENTS, "gradients", CPU_HYDRO, 'D', 'P')
TIMER_CREATE(CPU_FLUXES, "fluxes", CPU_HYDRO, 'F', 'N')
TIMER_CREATE(CPU_FLUXES_COMM, "fluxcomm", CPU_HYDRO, 'H', 'L')
TIMER_CREATE(CPU_CELL_UPDATES, "updates", CPU_HYDRO, 'J', 'j')
TIMER_CREATE(CPU_SET_VERTEXVELS, "vertex vel", CPU_HYDRO, 'K', 'I')
TIMER_CREATE(CPU_MHD, "mhd", CPU_HYDRO, '4', 'p')
TIMER_CREATE(CPU_DOMAIN, "domain", CPU_ALL, 'U', 'y')
TIMER_CREATE(CPU_PEANO, "peano", CPU_ALL, 'V', 'x')
TIMER_CREATE(CPU_DRIFTS, "drift/kicks", CPU_ALL, 'W', 'w')
TIMER_CREATE(CPU_TIMELINE, "timeline", CPU_ALL, 'X', 'v')
#ifdef TREE_BASED_TIMESTEPS
TIMER_CREATE(CPU_TREE_TIMESTEPS, "treetimesteps", CPU_ALL, 'Y', 'u')
#endif /* #ifdef TREE_BASED_TIMESTEPS */
TIMER_CREATE(CPU_SNAPSHOT, "i/o", CPU_ALL, 'Z', 't')
TIMER_CREATE(CPU_LOGS, "logs", CPU_ALL, '1', 's')
TIMER_CREATE(CPU_COOLINGSFR, "sfrcool", CPU_ALL, '2', 'r')
#ifdef FOF
TIMER_CREATE(CPU_FOF, "fof", CPU_ALL, '#', 'h')
#endif /* #ifdef FOF */
#ifdef SUBFIND
TIMER_CREATE(CPU_SUBFIND, "subfind", CPU_ALL, '$', 'g')
#endif /* #ifdef SUBFIND */
TIMER_CREATE(CPU_REFINE, "refine", CPU_ALL, '%', 'f')
TIMER_CREATE(CPU_DEREFINE, "mesh_derefine", CPU_ALL, '^', 'e')
TIMER_CREATE(CPU_MAKEIMAGES, "images", CPU_ALL, '&', 'd')
TIMER_CREATE(CPU_INIT, "initializ.", CPU_ALL, '*', 'c')
TIMER_CREATE(CPU_RESTART, "restart", CPU_ALL, '(', 'b')
TIMER_CREATE(CPU_MISC, "misc", CPU_ALL, ')', 'a')
TIMER_CREATE(CPU_LAST, "LAST", CPU_NONE, ' ', ' ') /*!<last item, do not use! */
#ifndef TIMER_STRUCT
}
;

extern enum timers TimerStack[TIMER_STACK_DEPTH];
extern int TimerStackPos;

/*! \brief struct containing the information of a CPU timer
 *
 */
struct timer_d
{
  int parent;         /*!< id of the parent timer */
  char shortname[40]; /*!< string containing the internal name of the timer */
  char longname[40];  /*!< name of the timer */
  char symb;          /*!< symbol used in balance.txt for the active part */
  char symbImbal;     /*!< symbol used in balance.txt for imbalance */
  char depth;         /*!< depth in the tree-like structure of this timer */
};
extern struct timer_d Timer_data[CPU_LAST + 1];
#else /* #ifndef TIMER_STRUCT */
#undef TIMER_STRUCT
#endif /* #ifndef TIMER_STRUCT #else */
#endif /* #if !defined(TIMER_H) || defined(TIMER_STRUCT) */
