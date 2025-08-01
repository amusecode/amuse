/* nstab.f -- translated by f2c (version 20100827).
   You must link the resulting object file with libf2c:
	on Microsoft Windows system, link with libf2c.lib;
	on Linux or Unix systems, link with .../path/to/libf2c.a -lm
	or, if you install libf2c.a in a standard place, with -lf2c -lm
	-- in that order, at the end of the command line, as in
		cc *.o -lf2c -lm
	Source for libf2c is in /netlib/f2c/libf2c.zip, e.g.,

		http://www.netlib.org/f2c/libf2c.zip
*/

#include "f2c.h"

/* Common Block Declarations */

union {
    struct {
	doublereal mm1, mm2, mm3;
    } _1;
    struct {
	doublereal m1, m2, m3;
    } _2;
} params_;

#define params_1 (params_._1)
#define params_2 (params_._2)

/* Table of constant values */

static doublereal c_b2 = .66666666666666663;
static doublereal c_b3 = 1.3333333333333333;
static integer c__2 = 2;
static doublereal c_b6 = .666;
static integer c__0 = 0;
static doublereal c_b30 = 1.5;
static real c_b32 = 2.f;
static doublereal c_b39 = .3333;

/* 	general three-body stability algorithm */

/* 	system is unstable if nstab=1 is returned */
/* 	system is stable if nstab=0 is returned */

/* 	Rosemary Mardling */
/* 	School of Mathematical Sciences, Monash University */

/* 	version as of 16-11-07 */
/* 	email mardling@sci.monash.edu.au to be added to updates email list */
/* 	preprint on astro-ph by New Year :-) */

/* 	sigma=period ratio (outer/inner) **should be > 1** */
/* 	ei0=initial inner eccentricity */
/* 	eo=outer eccentricity */
/* 	relinc=relative inclination (radians) */
/* 	m1, m2, m3=masses (any units; m3=outer body) */

/*       valid for all inclinations */

/* 	MASS RATIO CONDITIONS */
/* 	valid for systems with at least one of m_2/m_1>0.05 OR */
/* 	m_3/m_1>0.05 (so that one could have, for example, m_2/m_1=0 and */
/* 	m_3/m_1=0.1) OR BOTH m2/m1>0.01 AND m3/m1>0.01 */
/* 	**future version will include other resonances to cover smaller */
/* 	mass ratios */

/* 	assumes resonance angle phi=0 because resonance overlap */
/* 	criterion doesn't recognize instability outside separatrix. */

integer nstab_(doublereal *sigma, doublereal *ei0, doublereal *eo, doublereal 
	*relinc, doublereal *m1, doublereal *m2, doublereal *m3)
{
    /* System generated locals */
    integer ret_val, i__1;
    doublereal d__1, d__2, d__3, d__4;

    /* Builtin functions */
    double atan(doublereal), pow_dd(doublereal *, doublereal *), sqrt(
	    doublereal), cos(doublereal), sin(doublereal);

    /* Local variables */
    static doublereal a, e;
    static integer n;
    static doublereal z__, c20, c22, c31, c33, ei, an, m12, ek, en, pi, m123, 
	    s201, mi2, mi3, s221, mo2, mo3, del, f20n, f22n, phi, win;
    extern doublereal ein_induced__(doublereal *, doublereal *, doublereal *, 
	    doublereal *);
    static doublereal enp1, s22m1;
    extern doublereal eoct_(doublereal *, doublereal *, doublereal *), flmn_(
	    integer *, integer *, integer *, doublereal *);
    static doublereal gam200, gam220, gam222, cosik, sinik, gam22m2, eoctmax;

    pi = atan(1.) * 4.;
    params_1.mm1 = *m1;
    params_1.mm2 = *m2;
    params_1.mm3 = *m3;
    m12 = *m1 + *m2;
    m123 = m12 + *m3;
    mi2 = *m3 / m123;
/* Computing 2nd power */
    d__1 = m12;
    d__2 = m12 / m123;
    mo2 = *m1 * *m2 / (d__1 * d__1) * pow_dd(&d__2, &c_b2);
    d__1 = m12 / m123;
    mi3 = *m3 / m12 * pow_dd(&d__1, &c_b3) * (*m1 - *m2) / m12;
/* Computing 2nd power */
    d__1 = m12;
    mo3 = *m1 * *m2 / (d__1 * d__1) * (m12 / m123) * (*m1 - *m2) / m12;
    c22 = .375f;
    c20 = .25f;
    c31 = sqrt(3.f) / 4.f;
    c33 = -sqrt(5.f) / 4.f;
    e = *eo;
/* 	inclination coefficients */
    win = 0.;
/* Computing 2nd power */
    d__1 = *ei0;
    a = sqrt(1 - d__1 * d__1) * cos(*relinc);
/* Computing 2nd power */
    d__1 = *ei0;
/* Computing 2nd power */
    d__2 = sin(*relinc);
/* Computing 2nd power */
    d__3 = *ei0;
/* Computing 2nd power */
    d__4 = sin(win) * sin(*relinc);
    z__ = (1 - d__1 * d__1) * (d__2 * d__2 + 1) + d__3 * d__3 * 5 * (d__4 * 
	    d__4);
/* Computing 2nd power */
    d__1 = z__;
/* Computing 4th power */
    d__2 = a, d__2 *= d__2;
/* Computing 2nd power */
    d__3 = a;
/* Computing 2nd power */
    d__4 = a;
    del = d__1 * d__1 + 25 + d__2 * d__2 * 16 - z__ * 10 - d__3 * d__3 * 20 - 
	    d__4 * d__4 * 8 * z__;
/* Computing 2nd power */
    d__2 = a;
    ek = sqrt((d__1 = (z__ + 1 - d__2 * d__2 * 4 + sqrt(del)) / 6.f, abs(d__1)
	    ));
/* Computing 2nd power */
    d__1 = ek;
    cosik = a / sqrt(1 - d__1 * d__1);
/* Computing 2nd power */
    d__1 = cosik;
    sinik = sqrt(1 - d__1 * d__1);
/* Computing 2nd power */
    d__1 = cosik + 1;
    gam222 = d__1 * d__1 * .25f;
/* Computing 2nd power */
    d__1 = 1 - cosik;
    gam22m2 = d__1 * d__1 * .25f;
/* Computing 2nd power */
    d__1 = sinik;
    gam220 = sqrt(1.5f) * .5f * (d__1 * d__1);
/* Computing 2nd power */
    d__1 = cosik;
    gam200 = (d__1 * d__1 * 3 - 1) * .5f;
/* 	induced inner eccentricity */
    ei = ein_induced__(sigma, ei0, &e, relinc);
/* 	octopole emax */
    if (*m1 != *m2) {
	eoctmax = eoct_(sigma, ei0, &e);
	ei = max(eoctmax,ei);
    }
    ei = max(ek,ei);
    ei = min(ei,1.);
    n = (integer) (*sigma);
    ret_val = 0;
/* 	[n:1](222) resonance */
/* Computing 3rd power */
    d__1 = ei;
/* Computing 5th power */
    d__2 = ei, d__3 = d__2, d__2 *= d__2;
    s221 = ei * -3 + d__1 * (d__1 * d__1) * 1.625f + d__3 * (d__2 * d__2) * 
	    .026041666666666668f;
/* Computing 3rd power */
    d__1 = 1 - e;
    f22n = flmn_(&c__2, &c__2, &n, &e) / (d__1 * (d__1 * d__1));
    an = (d__1 = c22 * 6 * s221 * f22n * (mi2 + mo2 * pow_dd(sigma, &c_b6)) * 
	    gam222, abs(d__1));
    phi = 0.;
/* Computing 2nd power */
    d__1 = *sigma - n;
    en = d__1 * d__1 * .5f - an * (cos(phi) + 1);
/* 	[n+1:1](222) resonance */
    i__1 = n + 1;
/* Computing 3rd power */
    d__1 = 1 - e;
    f22n = flmn_(&c__2, &c__2, &i__1, &e) / (d__1 * (d__1 * d__1));
    an = (d__1 = c22 * 6 * s221 * f22n * (mi2 + mo2 * pow_dd(sigma, &c_b6)) * 
	    gam222, abs(d__1));
/* Computing 2nd power */
    d__1 = *sigma - (n + 1);
    enp1 = d__1 * d__1 * .5f - an * (cos(phi) + 1);
    if (en < 0. && enp1 < 0.) {
	ret_val = 1;
    }
/* 	[n:1](22-2) resonance */
/* Computing 3rd power */
    d__1 = ei;
/* Computing 2nd power */
    d__2 = ei;
/* Computing 4th power */
    d__3 = ei, d__3 *= d__3;
    s22m1 = -(d__1 * (d__1 * d__1) * (d__2 * d__2 * 1880 + 4480 + d__3 * d__3 
	    * 1091)) / 15360.f;
/* Computing 3rd power */
    d__1 = 1 - e;
    f22n = flmn_(&c__2, &c__2, &n, &e) / (d__1 * (d__1 * d__1));
    an = (d__1 = c22 * 6 * s22m1 * f22n * (mi2 + mo2 * pow_dd(sigma, &c_b6)) *
	     gam22m2, abs(d__1));
    phi = 0.;
/* Computing 2nd power */
    d__1 = *sigma - n;
    en = d__1 * d__1 * .5f - an * (cos(phi) + 1);
/* 	[n+1:1](22-2) resonance */
    i__1 = n + 1;
/* Computing 3rd power */
    d__1 = 1 - e;
    f22n = flmn_(&c__2, &c__2, &i__1, &e) / (d__1 * (d__1 * d__1));
    an = (d__1 = c22 * 6 * s22m1 * f22n * (mi2 + mo2 * pow_dd(sigma, &c_b6)) *
	     gam22m2, abs(d__1));
/* Computing 2nd power */
    d__1 = *sigma - (n + 1);
    enp1 = d__1 * d__1 * .5f - an * (cos(phi) + 1);
    if (en < 0. && enp1 < 0.) {
	ret_val = 1;
    }
/* 	[n:1](202) resonance */
/* Computing 2nd power */
    d__1 = ei;
/* Computing 4th power */
    d__2 = ei, d__2 *= d__2;
/* Computing 6th power */
    d__3 = ei, d__3 *= d__3;
    s201 = ei * (d__1 * d__1 * 1152 - 9216 - d__2 * d__2 * 48 + d__3 * (d__3 *
	     d__3)) / 9216.f;
/* Computing 3rd power */
    d__1 = 1 - e;
    f22n = flmn_(&c__2, &c__2, &n, &e) / (d__1 * (d__1 * d__1));
    an = (d__1 = sqrt(c20 * c22) * 6 * s201 * f22n * (mi2 + mo2 * pow_dd(
	    sigma, &c_b6)) * gam220, abs(d__1));
    phi = 0.;
/* Computing 2nd power */
    d__1 = *sigma - n;
    en = d__1 * d__1 * .5f - an * (cos(phi) + 1);
/* 	[n+1:1](202) resonance */
    i__1 = n + 1;
/* Computing 3rd power */
    d__1 = 1 - e;
    f22n = flmn_(&c__2, &c__2, &i__1, &e) / (d__1 * (d__1 * d__1));
    an = (d__1 = sqrt(c20 * c22) * 6 * s201 * f22n * (mi2 + mo2 * pow_dd(
	    sigma, &c_b6)) * gam220, abs(d__1));
/* Computing 2nd power */
    d__1 = *sigma - (n + 1);
    enp1 = d__1 * d__1 * .5f - an * (cos(phi) + 1);
    if (en < 0. && enp1 < 0.) {
	ret_val = 1;
    }
/* 	[n:1](002) resonance */
/* Computing 2nd power */
    d__1 = ei;
/* Computing 4th power */
    d__2 = ei, d__2 *= d__2;
/* Computing 6th power */
    d__3 = ei, d__3 *= d__3;
    s201 = ei * (d__1 * d__1 * 1152 - 9216 - d__2 * d__2 * 48 + d__3 * (d__3 *
	     d__3)) / 9216.f;
/* Computing 3rd power */
    d__1 = 1 - e;
    f20n = flmn_(&c__2, &c__0, &n, &e) / (d__1 * (d__1 * d__1));
    an = (d__1 = c20 * 3 * s201 * f20n * (mi2 + mo2 * pow_dd(sigma, &c_b6)) * 
	    gam200, abs(d__1));
    phi = 0.;
/* Computing 2nd power */
    d__1 = *sigma - n;
    en = d__1 * d__1 * .5f - an * (cos(phi) + 1);
/* 	[n+1:1](002) resonance */
    i__1 = n + 1;
/* Computing 3rd power */
    d__1 = 1 - e;
    f20n = flmn_(&c__2, &c__0, &i__1, &e) / (d__1 * (d__1 * d__1));
    an = (d__1 = c20 * 3 * s201 * f20n * (mi2 + mo2 * pow_dd(sigma, &c_b6)) * 
	    gam200, abs(d__1));
/* Computing 2nd power */
    d__1 = *sigma - (n + 1);
    enp1 = d__1 * d__1 * .5f - an * (cos(phi) + 1);
    if (en < 0. && enp1 < 0.) {
	ret_val = 1;
    }
/* L8: */
    return ret_val;
} /* nstab_ */

/* --------------------------------------------------------------------------- */
/* 	Asymptotic expression for f^(lm)_n(e) for all e<1 and n. */

doublereal flmn_(integer *l, integer *m, integer *n, doublereal *e)
{
    /* System generated locals */
    doublereal ret_val, d__1, d__2, d__3;

    /* Builtin functions */
    double atan(doublereal), pow_dd(doublereal *, doublereal *), sqrt(
	    doublereal), exp(doublereal), pow_ri(real *, integer *), pow_di(
	    doublereal *, integer *);

    /* Local variables */
    static doublereal pi, xi, rho;
    extern doublereal acosh_(doublereal *), facfac_(integer *, integer *);

    if (*e < .005f) {
	if (*m == *n) {
	    ret_val = 1.;
	} else {
	    ret_val = 0.;
	}
	return ret_val;
    }
    pi = atan(1.) * 4.;
    d__1 = 1 - *e;
    rho = *n * pow_dd(&d__1, &c_b30);
    d__1 = 1 / *e;
/* Computing 2nd power */
    d__2 = *e;
    d__3 = 1 - *e;
    xi = (acosh_(&d__1) - sqrt(1 - d__2 * d__2)) / pow_dd(&d__3, &c_b30);
    d__1 = *e + 1;
    d__2 = (doublereal) ((real) (*m * 3 - *l - 1) / 4.f);
    d__3 = (doublereal) ((real) (*l + *m + 1) / 2.f);
    ret_val = 1 / (pi * 2 * *n) * pow_ri(&c_b32, m) * (sqrt(pi * 2) / facfac_(
	    l, m)) * (pow_dd(&d__1, &d__2) / pow_di(e, m)) * pow_dd(&rho, &
	    d__3) * exp(-rho * xi);
    return ret_val;
} /* flmn_ */

/* --------------------------------------------------------------------------- */
doublereal ein_induced__(doublereal *sigma, doublereal *ei0, doublereal *e, 
	doublereal *relinc)
{
    /* System generated locals */
    doublereal ret_val, d__1, d__2;

    /* Builtin functions */
    double atan(doublereal), cos(doublereal), sqrt(doublereal), sin(
	    doublereal);

    /* Local variables */
    static doublereal a;
    static integer n;
    static doublereal pi, m123, f20n, f22n;
    extern doublereal flmn_(integer *, integer *, integer *, doublereal *);
    static doublereal prod, gam200, gam220, gam222, prod200, prod220, prod222;

    pi = atan(1.) * 4.;
    m123 = params_2.m1 + params_2.m2 + params_2.m3;
    n = (integer) (*sigma);
/* Computing 2nd power */
    d__1 = cos(*relinc) + 1;
    gam222 = d__1 * d__1 * .25f;
/* Computing 2nd power */
    d__1 = sin(*relinc);
    gam220 = sqrt(1.5f) * .5f * (d__1 * d__1);
/* Computing 2nd power */
    d__1 = cos(*relinc);
    gam200 = (d__1 * d__1 * 3 - 1) * .5f;
/* Computing 3rd power */
    d__1 = 1 - *e;
    f22n = flmn_(&c__2, &c__2, &n, e) / (d__1 * (d__1 * d__1));
/* Computing 3rd power */
    d__1 = 1 - *e;
    f20n = flmn_(&c__2, &c__0, &n, e) / (d__1 * (d__1 * d__1));
    prod222 = f22n * gam222;
    prod220 = f22n * gam220;
    prod200 = f20n * gam200;
/* Computing MAX */
    d__1 = max(prod222,prod220);
    prod = max(d__1,prod200);
/* Computing 2nd power */
    d__1 = *sigma;
    a = params_2.m3 / m123 * 4.5f * (pi * 2 * n) * prod / (d__1 * d__1);
/* Computing 2nd power */
    d__1 = *ei0;
/* Computing 2nd power */
    d__2 = a;
    ret_val = sqrt(d__1 * d__1 + d__2 * d__2);
    return ret_val;
} /* ein_induced__ */

/* ------------------------------------------------------------------------------ */
/* 	eoct.f */

/* 	calculates maximum eccentricity for arbitrary coplanar system */
/* 	using Mardling (2007) MNRAS in press */

doublereal eoct_(doublereal *sigma, doublereal *ei0, doublereal *eo)
{
    /* System generated locals */
    doublereal ret_val, d__1, d__2;

    /* Builtin functions */
    double atan(doublereal), pow_dd(doublereal *, doublereal *), sqrt(
	    doublereal);

    /* Local variables */
    static doublereal aa, al, m12, pi, m123, eeq, aoai, epso;

    pi = atan(1.) * 4.;
    m12 = params_2.m1 + params_2.m2;
    m123 = m12 + params_2.m3;
/* Computing 2nd power */
    d__2 = *sigma;
    d__1 = m123 / m12 * (d__2 * d__2);
    aoai = pow_dd(&d__1, &c_b39);
    al = 1 / aoai;
/* Computing 2nd power */
    d__1 = *eo;
    epso = sqrt(1 - d__1 * d__1);
/* Computing 2nd power */
    d__2 = epso;
    eeq = al * 1.25f * *eo / (d__2 * d__2) / (d__1 = 1 - sqrt(al) * (
	    params_2.m2 / params_2.m3) / epso, abs(d__1));
    aa = (d__1 = 1 - *ei0 / eeq, abs(d__1));
    if (aa < 1.) {
	ret_val = (aa + 1) * eeq;
    } else {
	ret_val = *ei0 + eeq * 2;
    }
    return ret_val;
} /* eoct_ */

/* --------------------------------------------------------------------------- */
doublereal acosh_(doublereal *x)
{
    /* System generated locals */
    doublereal ret_val, d__1;

    /* Builtin functions */
    double sqrt(doublereal), log(doublereal);

/* Computing 2nd power */
    d__1 = *x;
    ret_val = log(*x + sqrt(d__1 * d__1 - 1.));
    return ret_val;
} /* acosh_ */

/* --------------------------------------------------------------------------- */
doublereal sgn_(doublereal *x)
{
    /* System generated locals */
    doublereal ret_val;

    if (*x < 0.) {
	ret_val = -1.;
    } else {
	ret_val = 1.;
    }
    return ret_val;
} /* sgn_ */

/* --------------------------------------------------------------------------- */
doublereal facfac_(integer *l, integer *m)
{
    /* System generated locals */
    integer i__1;
    doublereal ret_val;

    /* Local variables */
    static integer i__, n;
    static doublereal prod;

    prod = 1.;
    n = *l + *m - 1;
    i__1 = n;
    for (i__ = 1; i__ <= i__1; i__ += 2) {
	prod *= i__;
    }
    ret_val = prod;
    return ret_val;
} /* facfac_ */

