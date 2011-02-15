
#include "stdinc.h"
#include "jdata.h"

// Global diagnostics:

real TDEBUG = -1; // 3.8593;

void printq(int j, real2 q, const char *s)
{
    cout << "    " << s << ": " << flush;
    real q2 = 0;
    for (int k = 0; k < 3; k++) {
	q2 += pow(q[j][k],2);
	cout << q[j][k] << " ";
    }
    cout << "   " << sqrt(q2) << endl << flush;
}

void print_particle(int j, jdata &jd)
{
    cout << "j = " << j << " (" << jd.id[j] << ")" << endl << flush;
    printq(j, jd.pos, "pos");
    printq(j, jd.vel, "vel");
    printq(j, jd.acc, "acc");
    printq(j, jd.jerk, "jerk");
    printq(j, jd.pred_pos, "ppos");
    printq(j, jd.pred_vel, "pvel");
}

void print_list(jdata &jd)
{
    int jlist[4] = {0,1,709,982}, nlist = 4;

    int p = cout.precision(10);
    cout << "print_list at time " << jd.system_time << endl;
    for (int jj = 0; jj < nlist; jj++) print_particle(jlist[jj], jd);
    cout.precision(p);
}

bool twiddles(real q, real v, real tol = 1.e-3)
{
    if (fabs(v) < 10*tol)
	return (fabs(q-v) < tol);
    else
	return (fabs(q-v) < tol*fabs(v));
}

