
// General functions and random numbers.

#include "stdinc.h"

void warning(const char *s)  {cout << s << endl << flush;}
void err_exit(const char *s) {warning(s); exit(1);}

#include <sys/time.h>
#ifdef _WIN32
#include <windows.h>
#else
#include <sys/resource.h>
#endif

// Note: the first call to a timer function sets its zero point.

static bool etset = false;
static timeval et0;
real get_elapsed_time()
{
    if (!etset) {
	gettimeofday(&et0, NULL);
	etset = true;
	cout << "initialized elapsed time counter" << endl << flush;
	return 0;
    } else {
	timeval et;			// et.tv_set = time(NULL), note
	gettimeofday(&et, NULL);
	return et.tv_sec - et0.tv_sec + 1.e-6*(et.tv_usec - et0.tv_usec);
    }
}

static bool ctset = false;
static real user0, sys0;
void get_cpu_time(real& user_time, real& system_time)
{
#ifdef _WIN32
    HANDLE hProcess = GetCurrentProcess();
    FILETIME ftCreation, ftExit, ftUser, ftKernel;
    ULONGLONG user;
    GetProcessTimes (hProcess, &ftCreation, &ftExit, &ftKernel, &ftUser); 
    user =  (((ULONGLONG) ftUser.dwHighDateTime) << 32) + ftUser.dwLowDateTime;
    user_time = user / 10000000U; 
    user =  (((ULONGLONG) ftKernel.dwHighDateTime) << 32) + ftKernel.dwLowDateTime;
    system_time = user / 10000000U; 
    
#else
    struct rusage tim;
    getrusage(RUSAGE_SELF, &tim);

    if (!ctset) {
	user0 = tim.ru_utime.tv_sec + 1.e-6*tim.ru_utime.tv_usec;
	sys0 = tim.ru_stime.tv_sec + 1.e-6*tim.ru_stime.tv_usec;
	user_time = system_time = 0;
	ctset = true;
	cout << "initialized CPU time counter" << endl << flush;
    } else {
	user_time = tim.ru_utime.tv_sec + 1.e-6*tim.ru_utime.tv_usec - user0;
	system_time = tim.ru_stime.tv_sec + 1.e-6*tim.ru_stime.tv_usec - sys0;
    }
#endif
}

static real scale = 1.0/(((unsigned long)1<<31)-1.0);
real randinter(real a, real b)
{
#ifdef _WIN32
    unsigned long r = rand();
#else
    unsigned long r = random();
#endif
    real x = scale*r;
    return a + (b-a)*x;

}
