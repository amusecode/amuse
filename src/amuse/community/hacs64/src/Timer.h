#ifndef __TIMER_H__
#define __TIMER_H__

#include <cstdio>

#include <sys/time.h>
inline double get_wtime() {
  struct timeval Tvalue;
  struct timezone dummy;
  
  gettimeofday(&Tvalue,&dummy);
  return ((double) Tvalue.tv_sec +1.e-6*((double) Tvalue.tv_usec));
}

struct Timer{
	const char *name;
	FILE       *fp;
	const char *format;
	double tstart;

	Timer(
			const char *_name,
			FILE       *_fp     = stdout,
			const char *_format = " %-10s : %f sec\n")
	   : name(_name), fp(_fp), format(_format)
	{
		tstart = get_wtime();
	}
	~Timer(){
		double tend = get_wtime();
		fprintf(fp, format, name, tend - tstart);
		fflush(fp);
	}
};

#endif // __TIMER_H__
