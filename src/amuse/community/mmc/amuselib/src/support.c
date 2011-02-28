/* 
For info on the many points code (idea taken from this ref):
============================================================

Saff, E.B., Kuijlaars, A.B.J., "Distributing many points on a
sphere" The Mathematical Intelligencer, 19, 1, 5-11, 1997 

For info on the rnd points code (idea taken from this ref):
============================================================

*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define ABS(x) (x < 0 ? -(x) : (x))

typedef struct {
  double x,y,z;
} vector3;

void asl_normalize(vector3 *,double);
double asl_distance(vector3, vector3);

int asl_add(double x, double y, double *z)
{
  *z = x + y;
  return 0;
}

int asl_polar_to_cart(double r, double theta, double phi, 
		      double *x, double *y, double *z) 
{
  *x = r * sin(theta) * cos(phi);
  *y = r * sin(theta) * sin(phi);
  *z = r * cos(theta);
  return 0;
}

int asl_cart_to_polar()
{
  return 0;
}

double asl_real_modulo(double x, double n)
{
  return x - n * floor(x/n);
}

void asl_normalize(vector3 *p, double r)
{
   double l;

   l = r / sqrt(p->x*p->x + p->y*p->y + p->z*p->z);
   p->x *= l;
   p->y *= l;
   p->z *= l;
}

double asl_distance(vector3 p1, vector3 p2)
{
   vector3 p;

   p.x = p1.x - p2.x;
   p.y = p1.y - p2.y;
   p.z = p1.z - p2.z;
   return(sqrt(p.x*p.x + p.y*p.y + p.z*p.z));
}

int asl_rnd_points_on_sphere(double *x_req, double *y_req, double *z_req, 
			      int N)
{
  //Paul Bourke, July 1996

  int i,j,n;
  int counter = 0,countmax = 100;
  int minp1,minp2;
  double r,d,mind,maxd;
  vector3 p1,p2;
  vector3 *p;

  p = (vector3 *) malloc(N*sizeof(vector3));
  n = N;
  r = 1.0;
  countmax = 100;

  /* Create the initial random cloud */
  for (i=0;i<n;i++) {
    p[i].x = (rand()%1000)-500;
    p[i].y = (rand()%1000)-500;
    p[i].z = (rand()%1000)-500;
    asl_normalize(&p[i],r);
  }
  
  while (counter < countmax) {
    
    /* Find the closest two points */
    minp1 = 0;
    minp2 = 1;
    mind = asl_distance(p[minp1],p[minp2]);
    maxd = mind;
    for (i=0;i<n-1;i++) {
      for (j=i+1;j<n;j++) {
	if ((d = asl_distance(p[i],p[j])) < mind) {
	  mind = d;
	  minp1 = i;
	  minp2 = j;
	}
	if (d > maxd)
	  maxd = d;
      }
    }
    
    /*
      Move the two minimal points apart, in this case by 1%
      but should really vary this for refinement
    */
    p1 = p[minp1];
    p2 = p[minp2];
    p[minp2].x = p1.x + 1.01 * (p2.x - p1.x);
    p[minp2].y = p1.y + 1.01 * (p2.y - p1.y);
    p[minp2].z = p1.z + 1.01 * (p2.z - p1.z);
    p[minp1].x = p1.x - 0.01 * (p2.x - p1.x);
    p[minp1].y = p1.y - 0.01 * (p2.y - p1.y);
    p[minp1].z = p1.z - 0.01 * (p2.z - p1.z);
    asl_normalize(&p[minp1],r);
    asl_normalize(&p[minp2],r);
    
    counter++;
  }
  
  /* Write out the points in your favorite format */
  //CLEAN UP! CELLO
  // 
  //\/\/
  
  for (i=0; i<n; i++) {
    x_req[i] = p[i].x;
    y_req[i] = p[i].y;
    z_req[i] = p[i].z;
  }
  
  free(p);
  return 0;
}

int asl_many_points_on_sphere(double *x_req, double *y_req, double *z_req, 
			      int N)
{
  double h;
  double x, y, z;
  int k;

  double *phi ;
  double *theta ;

  theta = (double *)malloc(N * sizeof(double));
  phi = (double *)malloc(N * sizeof(double));

  for (k=1; k<N; k++) {
    h = -1.0 + 2.0*(k-1)/(N-1);
    theta[k] = acos(h);
    if ((k==1) | (k==N)) {
      phi[k] = 0.0;
    } else {
      phi[k] = asl_real_modulo(phi[k-1] + 3.6/sqrt(N*(1-h*h)), (2.0*M_PI));
    }

    asl_polar_to_cart(1.0, theta[k], phi[k], &x, &y, &z);

    x_req[k-1] = x;
    y_req[k-1] = y;
    z_req[k-1] = z;

    //fprintf(stdout, "%f %f %f\n", x, y, z);
  }
  //fflush(stdout);
  free(theta);
  free(phi);
  return 0;
}

