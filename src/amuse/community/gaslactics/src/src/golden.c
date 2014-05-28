#include <math.h>

#define R 0.61803399
#define C (1.0-R)
#define SHFT(a,b,c,d) (a)=(b);(b)=(c);(c)=(d);

float golden(ax,bx,cx,f,tol,xmin)
float ax,bx,cx,tol,*xmin;
float (*f)();	/* ANSI: float (*f)(float); */
{
	float f0,f1,f2,f3,x0,x1,x2,x3;

	x0=ax;
	x3=cx;
	if (fabs(cx-bx) > fabs(bx-ax)) {
		x1=bx;
		x2=bx+C*(cx-bx);
	} else {
		x2=bx;
		x1=bx-C*(bx-ax);
	}
	f1=(*f)(x1);
	f2=(*f)(x2);
	while (fabs(x3-x0) > tol*(fabs(x1)+fabs(x2))) {
		if (f2 < f1) {
			SHFT(x0,x1,x2,R*x1+C*x3)
			SHFT(f0,f1,f2,(*f)(x2))
		} else {
			SHFT(x3,x2,x1,R*x2+C*x0)
			SHFT(f3,f2,f1,(*f)(x1))
		}
	}
	if (f1 < f2) {
		*xmin=x1;
		return f1;
	} else {
		*xmin=x2;
		return f2;
	}
}

#undef C
#undef R
