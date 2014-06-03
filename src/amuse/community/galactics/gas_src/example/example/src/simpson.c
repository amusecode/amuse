/* Integrate the function 'fcn' using Simpson's rule */
/*  x0, xn - limits of integrand
 		 n - number of divisions used in Simpson's rule */

float simpson(fcn,x0,xn,n)
float (*fcn)();
float x0, xn;
int n;
{
	int i;
	float dx, x, sum;

	dx = (xn-x0)/n;
	sum = (*fcn)(x0) + 4*(*fcn)(x0+dx) + (*fcn)(xn);
	for(i=2; i<n; i+=2) {
		x = x0 + i*dx;
		sum += 2*(*fcn)(x) + 4*(*fcn)(x+dx);
	}
	return(sum*dx/3.0);
}
