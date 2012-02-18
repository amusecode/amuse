/*
 *  time.C: Interface to some standard UNIX time functions.
 *
 *          WARNING! -- the "standard" functions are anything but standard,.
 *                      They may be both non-portable and subject to change!
 *
 *                      The code below is tested only in SunOS 4.1.3 and HPUX.
 *
 *.............................................................................
 *    version 1:  Mar 1994   Steve McMillan
 *    version 2:  
 *.............................................................................
 *  non-local functions: 
 *    cpu_time
 *.............................................................................
 */

#include "stdinc.h"

#ifndef NO_CPU_TIME
#   include <sys/times.h>
    struct tms buffer;
#   include <unistd.h>
    static long ticks_per_sec = 0;
#endif

#ifndef NO_CPU_TIME
  static real initial_cpu = 0;
#endif

//-----------------------------------------------------------------------------
//  cpu_init: initialize the CPU timer.
//-----------------------------------------------------------------------------

void cpu_init()
{
#ifndef NO_CPU_TIME
    times(&buffer);

    // Use both system and user time because of ambiguities
    // with Linux multiprocessing...

    initial_cpu = (real) (buffer.tms_utime + buffer.tms_stime);

    ticks_per_sec = sysconf(_SC_CLK_TCK);	// Clock ticks per second
#endif
}

//-----------------------------------------------------------------------------
//  cpu_time: return total processor time (in seconds) used since the timer
//            was last initialized by a call to cpu_init.
//-----------------------------------------------------------------------------

real cpu_time()
{
#ifndef NO_CPU_TIME

    if (!ticks_per_sec)
	cpu_init();

    times(&buffer);
    return ((real) (buffer.tms_utime + buffer.tms_stime - initial_cpu))
      			/ ticks_per_sec;
#else
    return 0;
#endif
}

// STARLAB_WAIT:	Wait a specified number of microseconds.

#ifdef HAS_NANOSLEEP
#include <time.h>
#endif

void starlab_wait(int iwait)
{
#ifdef HAS_USLEEP
    usleep(iwait);
#else
#ifdef HAS_NANOSLEEP
    timespec wait_time;
    twait_ime.tv_nsec = 1000*iwait;
    nanosleep(wait_time);
#endif
#endif
}
