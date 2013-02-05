#include <functions.h>
#include <types.h>
#include <cmath>

using namespace std;


HostError Corrector(double *GTIME, double *ATIME, double *local_time, double *step, int *next, unsigned long nextsize, double4 *pos_CH, double4 *vel_CH, double4 *a_H1, double4 *a_H0, double4 *a3_H, double4 *p_v_a3, double ETA6, double ETA4, double DTMAX, double DTMIN, unsigned int N){

   for(unsigned int i = 0; i < nextsize; i++){

      double dt;
      int who = next[i];
      int who1 = who + N;
      int who2 = who1 + N;

      double amod = a_H0[who].x * a_H0[who].x +
                    a_H0[who].y * a_H0[who].y +
                    a_H0[who].z * a_H0[who].z ;

      double adotmod = a_H0[who1].x * a_H0[who1].x +
                       a_H0[who1].y * a_H0[who1].y +
                       a_H0[who1].z * a_H0[who1].z ;

      double a2dotsmod =  a_H0[who2].x * a_H0[who2].x +
                          a_H0[who2].y * a_H0[who2].y +
                          a_H0[who2].z * a_H0[who2].z ;

      double h = *GTIME-local_time[who];
      local_time[who] = *GTIME;


      double h1 = 0.5*h;
      double h2 = h*h*0.1;
      double h3 = h*h*h / 120.0;


   pos_CH[who].x = pos_CH[who].x +
     h1 * vel_CH[who].x -
     h2 * (a_H0[who].x - a_H1[who].x) +
     h3 * (a_H1[who1].x + a_H0[who1].x);
   pos_CH[who].y = pos_CH[who].y +
     h1 * vel_CH[who].y -
     h2 * (a_H0[who].y - a_H1[who].y) +
     h3 * (a_H1[who1].y + a_H0[who1].y);
	pos_CH[who].z = pos_CH[who].z +
     h1 * vel_CH[who].z -
     h2 * (a_H0[who].z - a_H1[who].z) +
     h3 * (a_H1[who1].z + a_H0[who1].z);
   vel_CH[who].x = vel_CH[who].x +
     h1 * (a_H0[who].x + a_H1[who].x) -
     h2 * (a_H0[who1].x - a_H1[who1].x) +
     h3 * (a_H0[who2].x + a_H1[who2].x);
   vel_CH[who].y = vel_CH[who].y +
     h1 * (a_H0[who].y  + a_H1[who].y) -
     h2 * (a_H0[who1].y - a_H1[who1].y) +
     h3 * (a_H0[who2].y + a_H1[who2].y);
   vel_CH[who].z = vel_CH[who].z +
     h1 * (a_H0[who].z  + a_H1[who].z) -
     h2 * (a_H0[who1].z - a_H1[who1].z) +
     h3 * (a_H0[who2].z + a_H1[who2].z);


   pos_CH[who].x += h1*vel_CH[who].x;
   pos_CH[who].y += h1*vel_CH[who].y;
   pos_CH[who].z += h1*vel_CH[who].z;

   p_v_a3[i] = pos_CH[who];

   p_v_a3[i+nextsize].x = vel_CH[who].x;
   p_v_a3[i+nextsize].y = vel_CH[who].y;
   p_v_a3[i+nextsize].z = vel_CH[who].z;

   h2 = h1*h1;
   h3 = 1.0/(h2*h1);
   double h4 = h3/h1;
   double h5 = h4/h1;
	
	h3 *= 0.75;
   h4 *= 1.5;
   h5 *= 7.5;

    double Amin = a_H0[who].x - a_H1[who].x;
    double Jmin = h1 * (a_H0[who1].x - a_H1[who1].x);
    double Jplu = h1 * (a_H0[who1].x + a_H1[who1].x);
    double Smin = h2 * (a_H0[who2].x - a_H1[who2].x);
    double Splu = h2 * (a_H0[who2].x + a_H1[who2].x);

    a3_H[who].x = h3*(-5.0*Amin + 5.0*Jplu - Smin);
    double a4halfx = h4*(-Jmin + Splu);
    double a5halfx = h5*(3.0*Amin - 3.0*Jplu + Smin);
    a3_H[who].x += h1*a4halfx + h2/2.0*a5halfx;
    a4halfx += h1*a5halfx;

    Amin = a_H0[who].y - a_H1[who].y;
    Jmin = h1 * (a_H0[who1].y - a_H1[who1].y);
    Jplu = h1 * (a_H0[who1].y + a_H1[who1].y);
    Smin = h2 * (a_H0[who2].y - a_H1[who2].y);
    Splu = h2 * (a_H0[who2].y + a_H1[who2].y);

    a3_H[who].y = h3*(-5.0*Amin + 5.0*Jplu - Smin);
    double a4halfy = h4*(-Jmin + Splu);
    double a5halfy = h5*(3.0*Amin - 3.0*Jplu + Smin);

    a3_H[who].y += h1*a4halfy + h2/2.0*a5halfy;
    a4halfy += h1*a5halfy;

    Amin = a_H0[who].z - a_H1[who].z;
    Jmin = h1 * (a_H0[who1].z - a_H1[who1].z);
    Jplu = h1 * (a_H0[who1].z + a_H1[who1].z);
    Smin = h2 * (a_H0[who2].z - a_H1[who2].z);
    Splu = h2 * (a_H0[who2].z + a_H1[who2].z);

    a3_H[who].z = h3*(-5.0*Amin + 5.0*Jplu - Smin);
    double a4halfz = h4*(-Jmin + Splu);
    double a5halfz = h5*(3.0*Amin - 3.0*Jplu + Smin);

    a3_H[who].z += h1*a4halfz + h2/2.0*a5halfz;
    a4halfz += h1*a5halfz;

	double a3dotsmod = sqrt(a3_H[who].x*a3_H[who].x + a3_H[who].y*a3_H[who].y + a3_H[who].z*a3_H[who].z);

    double a4mod = sqrt(a4halfx*a4halfx + a4halfy*a4halfy + a4halfz*a4halfz);
    double a5mod = sqrt(a5halfx*a5halfx + a5halfy*a5halfy + a5halfz*a5halfz);

    double    dt6 = (sqrt(amod*a2dotsmod) + adotmod) / (a5mod*a3dotsmod + a4mod*a4mod);
    dt6 = ETA6 * pow(dt6,1.0/6.0);

     double stp = h;
      double overh3 = 1.0/(stp*stp*stp);
      double overh2 = 1.0/(stp*stp);

    double a2dx = overh2 * (-6.0 * (a_H1[who].x - a_H0[who].x) -
             stp * (4.0 * a_H0[who1].x + 2.0 * a_H1[who1].x));
    double a2dy = overh2 * (-6.0 * (a_H1[who].y - a_H0[who].y) -
             stp * (4.0 * a_H0[who1].y + 2.0 * a_H1[who1].y));
    double a2dz = overh2 * (-6.0 * (a_H1[who].z - a_H0[who].z) -
             stp * (4.0 * a_H0[who1].z + 2.0 * a_H1[who1].z));

    double a3dx = overh3 * (12.0 * (a_H1[who].x - a_H0[who].x) +
             6.0 * stp * (a_H0[who1].x + a_H1[who1].x));
    double a3dy = overh3 * (12.0 * (a_H1[who].y - a_H0[who].y) +
             6.0 * stp * (a_H0[who1].y + a_H1[who1].y));
    double a3dz = overh3 * (12.0 * (a_H1[who].z - a_H0[who].z) +
             6.0 * stp * (a_H0[who1].z + a_H1[who1].z));

    a2dx += h*a3dx;
    a2dy += h*a3dy;
    a2dx += h*a3dz;

    a2dotsmod =  a2dx*a2dx + a2dy*a2dy + a2dz*a2dz;
    a3dotsmod = a3dx*a3dx + a3dy*a3dy + a3dz*a3dz;

    double dt4 = sqrt(ETA4*(sqrt(amod*a2dotsmod) + adotmod) / (sqrt(adotmod*a3dotsmod) + a2dotsmod));


    dt = 0.5*dt4+0.5*dt6;

	 double rest = *GTIME / (2.0 * step[who]);
   rest = (double)((int)(rest)) - rest;

   if(dt > 2.0*step[who] && rest == 0.0 && 2.0*step[who] <= DTMAX)
     step[who] *= 2.0;
   else if (dt < 0.5*step[who])
     step[who] *= 0.25;
   else if (dt < step[who])
      step[who]*=0.5;

   if(step[who] < DTMIN)
      step[who] = DTMIN;

   p_v_a3[i+2*nextsize].x = a3_H[who].x;
   p_v_a3[i+2*nextsize].y = a3_H[who].y;
   p_v_a3[i+2*nextsize].z = a3_H[who].z;

   *ATIME = min (local_time[who] + step[who], *ATIME);

   a_H1[who].x = a_H0[who].x;
   a_H1[who].y = a_H0[who].y;
   a_H1[who].z = a_H0[who].z;
   a_H1[who1].x = a_H0[who1].x;
   a_H1[who1].y = a_H0[who1].y;
   a_H1[who1].z = a_H0[who1].z;
   a_H1[who2].x = a_H0[who2].x;
   a_H1[who2].y = a_H0[who2].y;
   a_H1[who2].z = a_H0[who2].z;


}


   return HNoError;

}

