#include "debug.h"

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

void print_list(int jlist[], int nlist, jdata &jd)
{
    int p = cout.precision(10);
    cout << "print_list at time " << jd.system_time << endl;
    for (int jj = 0; jj < nlist; jj++) print_particle(jlist[jj], jd);
    cout.precision(p);
}

void print_list(jdata &jd)
{
    int jlist[4] = {0,1,709,982}, nlist = 4;
    print_list(jlist, nlist, jd);
}

bool twiddles(real q, real v, real tol)
{
    if (fabs(v) < 10*tol)
	return (fabs(q-v) < tol);
    else
	return (fabs(q-v) < tol*fabs(v));
}

real total_energy(int jlist[], int n, jdata& jd)
{
    // Compute the total self-energy of the listed particles.

    real ke = 0, pe = 0;
    for (int i = 0; i < n; i++) {
	int j = jlist[i];
	real v2 = 0;
	for (int k = 0; k < 3; k++) v2 += pow(jd.vel[j][k], 2);
	ke += 0.5*jd.mass[j]*v2;
	real ppe = 0;
	for (int i1 = i+1; i1 < n; i1++) {
	    int j1 = jlist[i1];
	    real r2 = jd.eps2;
	    for (int k = 0; k < 3; k++)
		r2 += pow(jd.pos[j][k] - jd.pos[j1][k], 2);
	    ppe -= jd.mass[j1]/sqrt(r2);
	}
	pe += jd.mass[j]*ppe;
    }

    return pe+ke;
}
