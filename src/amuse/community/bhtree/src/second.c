#define RETURN_CPU
/*
 * SECOND.c  : returns the CPU time
 *
 * Jun Makino 89/09/20 Ver. 1.00 (created)
 *            89/10/02 Ver. 1.01 debugged on UNIX system.
 *
 * functions :
 *             void timer_init() : reset the timer
 *             void second(&double): returns the cpu seconds
 *             void cpumin(&double): returns the cpu minutes
 */
 
#ifdef MSC
#include <time.h>

static long start_time;

void timer_init()
{
	time(&start_time);
}

second(dtime)
double * dtime;
{
	long int timenow;
	time(&timenow);
	*dtime = timenow - start_time;
}

cpumin(dmin)
double *dmin ;
{
	double sec;
	second(&sec);
	*dmin = sec/60.0;
}
#else

/*
 *  SECONDS:  returns number of CPU seconds spent so far
 *		This version can be called from fortran routines
 *		in double precision mode
 */

#include <stdio.h>
#include <sys/types.h>
#include <sys/time.h>
#include <sys/times.h>
#include <sys/resource.h>
#include <limits.h>
static double tstart;
static struct timeval timearg;
static struct timezone zonearg;

#define RETURN_CPU    

void timer_init()
{
#ifndef RETUEN_CPU
    if(gettimeofday(&timearg,&zonearg)){
	fprintf(stderr,"Time initialization failed\n");
    }
    tstart = timearg.tv_sec + timearg.tv_usec*1e-6;
#else
#ifndef SOLARIS
    struct rusage usage;
    if(getrusage(RUSAGE_SELF,&usage)){
	fprintf(stderr,"getrusage failed\n");
    }
    tstart =  usage.ru_utime.tv_sec + usage.ru_utime.tv_usec*1e-6;
#else
    struct tms buffer;
    times(&buffer);
    tstart = (buffer.tms_utime + 0.0)/CLK_TCK;
#endif
#endif
}

void tminit_()
{
    timer_init();
}

cpumin(t)
double * t;
{
	xcpumin_(t);
}

double cpusec()
{
    double t;
    second(&t);
    return t;
}
second(t)
double * t;
{
	second_(t);
}
xcpumin_(t)
    double*t;
{
    double sec;
    second_(&sec);
    *t=sec/60.0;
}

second_(t)			
double *t;
{
#ifdef RETURN_CPU    
#ifdef SOLARIS
    struct tms buffer;
        
        if (times(&buffer) == -1) {
	        printf("times() call failed\n");
	        exit(1);
	}
	*t =  (buffer.tms_utime / (CLK_TCK+0.0));
#else
    struct rusage usage;
    if(getrusage(RUSAGE_SELF,&usage)){
	fprintf(stderr,"getrusage failed\n");
    }
    *t =  usage.ru_utime.tv_sec + usage.ru_utime.tv_usec*1e-6;
#endif      
#else
#ifndef xxxx
    if(gettimeofday(&timearg,&zonearg)){
	fprintf(stderr,"Timer failed\n");
    }
    *t =  timearg.tv_sec + timearg.tv_usec*1e-6 - tstart;
#else
    struct tms tbuf;
    int it0;
    it0 = times(&tbuf);
    if(it0 == -1){
	fprintf(stderr,"Time initialization failed\n");
    }
    *t  = it0*0.01 - tstart;
#endif
#endif    
}
#endif

#ifdef TEST
main()
{
	double t;
	timer_init();
	for(;;){
		second(&t);
		printf("time = %f\n", t);
	}
}
#endif
#ifdef SHORT_REAL
#  define COMPTYPE float
#else
#  define COMPTYPE double
#endif
typedef union double_and_int{
    unsigned int iwork;
    COMPTYPE fwork;
} DOUBLE_INT;
void c_assignbody_(pos,cpos,bsub)
register COMPTYPE  pos[];
register COMPTYPE  cpos[];
int * bsub;
{
    register int k,l,k1,k2;
    register DOUBLE_INT tmpx,tmpy,tmpz ;
#if 1
    k = k1 = k2 =0;
    if(pos[0] >= cpos[0]) k=1;
    if(pos[1] >= cpos[1]) k1=2;
    if(pos[2] >= cpos[2]) k2=4;
    *bsub = 1 + k + k1 + k2;
#endif
#if 0
    k += signbit(pos[0]-cpos[0]);
    k += signbit(pos[1]-cpos[1])*2;
    k += signbit(pos[2]-cpos[2])*4;
    k = (k + 7)>>1;
    *bsub += k;
#endif
#if 0
    tmpx.fwork = cpos[0]-pos[0];
    tmpy.fwork = cpos[1]-pos[1];
    tmpz.fwork = cpos[2]-pos[2];
    *bsub = 1+(tmpx.iwork >> 31)
       + ((tmpy.iwork >> 30) & 2   )
       + ((tmpz.iwork >> 29) & 4   );
#if 0
    if(pos[0] >= cpos[0]) k++;
    if(pos[1] >= cpos[1]) k+=2;
    if(pos[2] >= cpos[2]) k+=4;
    printf("%f %f %f %x %x %x %x %x %x %x %x\n",
	   tmpx.fwork,tmpy.fwork,tmpz.fwork,
	   tmpx.iwork,tmpy.iwork,tmpz.iwork,
	   tmpx.iwork>>31,tmpy.iwork>>30,tmpz.iwork>>29,
	   *bsub,k);
#endif    
#endif
}

