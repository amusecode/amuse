/* CAUTION: This is the ANSI C (only) version of the Numerical Recipes
   utility file nr.h.  Do not confuse this file with the same-named
   file nr.h that may be supplied in a 'misc' subdirectory.
   *That* file is the one from the book, and contains both ANSI and
   traditional K&R versions, along with #ifdef macros to select the
   correct version.  *This* file contains only ANSI C.               */

#ifndef _NR_H_
#define _NR_H_

#ifndef _FCOMPLEX_DECLARE_T_
typedef struct FCOMPLEX {double r,i;} fcomplex;
#define _FCOMPLEX_DECLARE_T_
#endif /* _FCOMPLEX_DECLARE_T_ */

#ifndef _ARITHCODE_DECLARE_T_
typedef struct {
	unsigned long *ilob,*iupb,*ncumfq,jdif,nc,minint,nch,ncum,nrad;
} arithcode;
#define _ARITHCODE_DECLARE_T_
#endif /* _ARITHCODE_DECLARE_T_ */

#ifndef _HUFFCODE_DECLARE_T_
typedef struct {
	unsigned long *icod,*ncod,*left,*right,nch,nodemax;
} huffcode;
#define _HUFFCODE_DECLARE_T_
#endif /* _HUFFCODE_DECLARE_T_ */

#include <stdio.h>

void addint(double **uf, double **uc, double **res, int nf);
void airy(double x, double *ai, double *bi, double *aip, double *bip);
void amebsa(double **p, double y[], int ndim, double pb[],	double *yb,
	double ftol, double (*funk)(double []), int *iter, double temptr);
void amoeba(double **p, double y[], int ndim, double ftol,
	double (*funk)(double []), int *iter);
double amotry(double **p, double y[], double psum[], int ndim,
	double (*funk)(double []), int ihi, double fac);
double amotsa(double **p, double y[], double psum[], int ndim, double pb[],
	double *yb, double (*funk)(double []), int ihi, double *yhi, double fac);
void anneal(double x[], double y[], int iorder[], int ncity);
double anorm2(double **a, int n);
void arcmak(unsigned long nfreq[], unsigned long nchh, unsigned long nradd,
	arithcode *acode);
void arcode(unsigned long *ich, unsigned char **codep, unsigned long *lcode,
	unsigned long *lcd, int isign, arithcode *acode);
void arcsum(unsigned long iin[], unsigned long iout[], unsigned long ja,
	int nwk, unsigned long nrad, unsigned long nc);
void asolve(unsigned long n, double b[], double x[], int itrnsp);
void atimes(unsigned long n, double x[], double r[], int itrnsp);
void avevar(double data[], unsigned long n, double *ave, double *var);
void balanc(double **a, int n);
void banbks(double **a, unsigned long n, int m1, int m2, double **al,
	unsigned long indx[], double b[]);
void bandec(double **a, unsigned long n, int m1, int m2, double **al,
	unsigned long indx[], double *d);
void banmul(double **a, unsigned long n, int m1, int m2, double x[], double b[]);
void bcucof(double y[], double y1[], double y2[], double y12[], double d1,
	double d2, double **c);
void bcuint(double y[], double y1[], double y2[], double y12[],
	double x1l, double x1u, double x2l, double x2u, double x1,
	double x2, double *ansy, double *ansy1, double *ansy2);
void beschb(double x, double *gam1, double *gam2, double *gampl,
	double *gammi);
double bessi(int n, double x);
double bessi0(double x);
double bessi1(double x);
void bessik(double x, double xnu, double *ri, double *rk, double *rip,
	double *rkp);
double bessj(int n, double x);
double bessj0(double x);
double bessj1(double x);
void bessjy(double x, double xnu, double *rj, double *ry, double *rjp,
	double *ryp);
double bessk(int n, double x);
double bessk0(double x);
double bessk1(double x);
double bessy(int n, double x);
double bessy0(double x);
double bessy1(double x);
double beta(double z, double w);
double betacf(double a, double b, double x);
double betai(double a, double b, double x);
double bico(int n, int k);
void bksub(int ne, int nb, int jf, int k1, int k2, double ***c);
double bnldev(double pp, int n, long *idum);
double brent(double ax, double bx, double cx,
	double (*f)(double), double tol, double *xmin);
void broydn(double x[], int n, int *check,
	void (*vecfunc)(int, double [], double []));
void bsstep(double y[], double dydx[], int nv, double *xx, double htry,
	double eps, double yscal[], double *hdid, double *hnext,
	void (*derivs)(double, double [], double []));
void caldat(long julian, int *mm, int *id, int *iyyy);
void chder(double a, double b, double c[], double cder[], int n);
double chebev(double a, double b, double c[], int m, double x);
void chebft(double a, double b, double c[], int n, double (*func)(double));
void chebpc(double c[], double d[], int n);
void chint(double a, double b, double c[], double cint[], int n);
double chixy(double bang);
void choldc(double **a, int n, double p[]);
void cholsl(double **a, int n, double p[], double b[], double x[]);
void chsone(double bins[], double ebins[], int nbins, int knstrn,
	double *df, double *chsq, double *prob);
void chstwo(double bins1[], double bins2[], int nbins, int knstrn,
	double *df, double *chsq, double *prob);
void cisi(double x, double *ci, double *si);
void cntab1(int **nn, int ni, int nj, double *chisq,
	double *df, double *prob, double *cramrv, double *ccc);
void cntab2(int **nn, int ni, int nj, double *h, double *hx, double *hy,
	double *hygx, double *hxgy, double *uygx, double *uxgy, double *uxy);
void convlv(double data[], unsigned long n, double respns[], unsigned long m,
	int isign, double ans[]);
void copy(double **aout, double **ain, int n);
void correl(double data1[], double data2[], unsigned long n, double ans[]);
void cosft(double y[], int n, int isign);
void cosft1(double y[], int n);
void cosft2(double y[], int n, int isign);
void covsrt(double **covar, int ma, int ia[], int mfit);
void crank(unsigned long n, double w[], double *s);
void cyclic(double a[], double b[], double c[], double alpha, double beta,
	double r[], double x[], unsigned long n);
void daub4(double a[], unsigned long n, int isign);
double dawson(double x);
double dbrent(double ax, double bx, double cx,
	double (*f)(double), double (*df)(double), double tol, double *xmin);
void ddpoly(double c[], int nc, double x, double pd[], int nd);
int decchk(char string[], int n, char *ch);
void derivs(double x, double y[], double dydx[]);
double df1dim(double x);
void dfour1(double data[], unsigned long nn, int isign);
void dfpmin(double p[], int n, double gtol, int *iter, double *fret,
	double (*func)(double []), void (*dfunc)(double [], double []));
double dfridr(double (*func)(double), double x, double h, double *err);
void dftcor(double w, double delta, double a, double b, double endpts[],
	double *corre, double *corim, double *corfac);
void dftint(double (*func)(double), double a, double b, double w,
	double *cosint, double *sinint);
void difeq(int k, int k1, int k2, int jsf, int is1, int isf,
	int indexv[], int ne, double **s, double **y);
void dlinmin(double p[], double xi[], int n, double *fret,
	double (*func)(double []), void (*dfunc)(double [], double[]));
double dpythag(double a, double b);
void drealft(double data[], unsigned long n, int isign);
void dsprsax(double sa[], unsigned long ija[], double x[], double b[],
	unsigned long n);
void dsprstx(double sa[], unsigned long ija[], double x[], double b[],
	unsigned long n);
void dsvbksb(double **u, double w[], double **v, int m, int n, double b[],
	double x[]);
void dsvdcmp(double **a, int m, int n, double w[], double **v);
void eclass(int nf[], int n, int lista[], int listb[], int m);
void eclazz(int nf[], int n, int (*equiv)(int, int));
double ei(double x);
void eigsrt(double d[], double **v, int n);
double elle(double phi, double ak);
double ellf(double phi, double ak);
double ellpi(double phi, double en, double ak);
void elmhes(double **a, int n);
double erfcc(double x);
/*double erff(double x);*/
double erffc(double x);
void eulsum(double *sum, double term, int jterm, double wksp[]);
double evlmem(double fdt, double d[], int m, double xms);
double expdev(long *idum);
double expint(int n, double x);
double f1(double x);
double f1dim(double x);
double f2(double y);
double f3(double z);
double factln(int n);
double factrl(int n);
void fasper(double x[], double y[], unsigned long n, double ofac, double hifac,
	double wk1[], double wk2[], unsigned long nwk, unsigned long *nout,
	unsigned long *jmax, double *prob);
void fdjac(int n, double x[], double fvec[], double **df,
	void (*vecfunc)(int, double [], double []));
void fgauss(double x, double a[], double *y, double dyda[], int na);
void fill0(double **u, int n);
void fit(double x[], double y[], int ndata, double sig[], int mwt,
	double *a, double *b, double *siga, double *sigb, double *chi2, double *q);
void fitexy(double x[], double y[], int ndat, double sigx[], double sigy[],
	double *a, double *b, double *siga, double *sigb, double *chi2, double *q);
void fixrts(double d[], int m);
void fleg(double x, double pl[], int nl);
void flmoon(int n, int nph, long *jd, double *frac);
/* double fmin(double x[]); */
void four1(double data[], unsigned long nn, int isign);
void fourew(FILE *file[5], int *na, int *nb, int *nc, int *nd);
void fourfs(FILE *file[5], unsigned long nn[], int ndim, int isign);
void fourn(double data[], unsigned long nn[], int ndim, int isign);
void fpoly(double x, double p[], int np);
void fred2(int n, double a, double b, double t[], double f[], double w[],
	double (*g)(double), double (*ak)(double, double));
double fredin(double x, int n, double a, double b, double t[], double f[], double w[],
	double (*g)(double), double (*ak)(double, double));
void frenel(double x, double *s, double *c);
void frprmn(double p[], int n, double ftol, int *iter, double *fret,
	double (*func)(double []), void (*dfunc)(double [], double []));
void ftest(double data1[], unsigned long n1, double data2[], unsigned long n2,
	double *f, double *prob);
double gamdev(int ia, long *idum);
double gammln(double xx);
double gammp(double a, double x);
double gammq(double a, double x);
double gasdev(long *idum);
void gaucof(int n, double a[], double b[], double amu0, double x[], double w[]);
void gauher(double x[], double w[], int n);
void gaujac(double x[], double w[], int n, double alf, double bet);
void gaulag(double x[], double w[], int n, double alf);
void gauleg(double x1, double x2, double x[], double w[], int n);
void gaussj(double **a, int n, double **b, int m);
void gcf(double *gammcf, double a, double x, double *gln);
double golden(double ax, double bx, double cx, double (*f)(double), double tol,
	double *xmin);
void gser(double *gamser, double a, double x, double *gln);
void hpsel(unsigned long m, unsigned long n, double arr[], double heap[]);
void hpsort(unsigned long n, double ra[]);
void hqr(double **a, int n, double wr[], double wi[]);
void hufapp(unsigned long index[], unsigned long nprob[], unsigned long n,
	unsigned long i);
void hufdec(unsigned long *ich, unsigned char *code, unsigned long lcode,
	unsigned long *nb, huffcode *hcode);
void hufenc(unsigned long ich, unsigned char **codep, unsigned long *lcode,
	unsigned long *nb, huffcode *hcode);
void hufmak(unsigned long nfreq[], unsigned long nchin, unsigned long *ilong,
	unsigned long *nlong, huffcode *hcode);
void hunt(double xx[], unsigned long n, double x, unsigned long *jlo);
void hypdrv(double s, double yy[], double dyyds[]);
fcomplex hypgeo(fcomplex a, fcomplex b, fcomplex c, fcomplex z);
void hypser(fcomplex a, fcomplex b, fcomplex c, fcomplex z,
	fcomplex *series, fcomplex *deriv);
unsigned short icrc(unsigned short crc, unsigned char *bufptr,
	unsigned long len, short jinit, int jrev);
unsigned short icrc1(unsigned short crc, unsigned char onech);
unsigned long igray(unsigned long n, int is);
void iindexx(unsigned long n, long arr[], unsigned long indx[]);
void indexx(unsigned long n, double arr[], unsigned long indx[]);
void interp(double **uf, double **uc, int nf);
int irbit1(unsigned long *iseed);
int irbit2(unsigned long *iseed);
void jacobi(double **a, int n, double d[], double **v, int *nrot);
void jacobn(double x, double y[], double dfdx[], double **dfdy, int n);
long julday(int mm, int id, int iyyy);
void kendl1(double data1[], double data2[], unsigned long n, double *tau, double *z,
	double *prob);
void kendl2(double **tab, int i, int j, double *tau, double *z, double *prob);
void kermom(double w[], double y, int m);
void ks2d1s(double x1[], double y1[], unsigned long n1,
	void (*quadvl)(double, double, double *, double *, double *, double *),
	double *d1, double *prob);
void ks2d2s(double x1[], double y1[], unsigned long n1, double x2[], double y2[],
	unsigned long n2, double *d, double *prob);
void ksone(double data[], unsigned long n, double (*func)(double), double *d,
	double *prob);
void kstwo(double data1[], unsigned long n1, double data2[], unsigned long n2,
	double *d, double *prob);
void laguer(fcomplex a[], int m, fcomplex *x, int *its);
void lfit(double x[], double y[], double sig[], int ndat, double a[], int ia[],
	int ma, double **covar, double *chisq, void (*funcs)(double, double [], int));
void linbcg(unsigned long n, double b[], double x[], int itol, double tol,
	 int itmax, int *iter, double *err);
void linmin(double p[], double xi[], int n, double *fret,
	double (*func)(double []));
void lnsrch(int n, double xold[], double fold, double g[], double p[], double x[],
	 double *f, double stpmax, int *check, double (*func)(double []));
void load(double x1, double v[], double y[]);
void load1(double x1, double v1[], double y[]);
void load2(double x2, double v2[], double y[]);
void locate(double xx[], unsigned long n, double x, unsigned long *j);
void lop(double **out, double **u, int n);
void lubksb(double **a, int n, int *indx, double b[]);
void ludcmp(double **a, int n, int *indx, double *d);
void machar(int *ibeta, int *it, int *irnd, int *ngrd,
	int *machep, int *negep, int *iexp, int *minexp, int *maxexp,
	double *eps, double *epsneg, double *xmin, double *xmax);
void matadd(double **a, double **b, double **c, int n);
void matsub(double **a, double **b, double **c, int n);
void medfit(double x[], double y[], int ndata, double *a, double *b, double *abdev);
void memcof(double data[], int n, int m, double *xms, double d[]);
int metrop(double de, double t);
void mgfas(double **u, int n, int maxcyc);
void mglin(double **u, int n, int ncycle);
double midexp(double (*funk)(double), double aa, double bb, int n);
double midinf(double (*funk)(double), double aa, double bb, int n);
double midpnt(double (*func)(double), double a, double b, int n);
double midsql(double (*funk)(double), double aa, double bb, int n);
double midsqu(double (*funk)(double), double aa, double bb, int n);
void miser(double (*func)(double []), double regn[], int ndim, unsigned long npts,
	double dith, double *ave, double *var);
void mmid(double y[], double dydx[], int nvar, double xs, double htot,
	int nstep, double yout[], void (*derivs)(double, double[], double[]));
void mnbrak(double *ax, double *bx, double *cx, double *fa, double *fb,
	double *fc, double (*func)(double));
void mnewt(int ntrial, double x[], int n, double tolx, double tolf);
void moment(double data[], int n, double *ave, double *adev, double *sdev,
	double *var, double *skew, double *curt);
void mp2dfr(unsigned char a[], unsigned char s[], int n, int *m);
void mpadd(unsigned char w[], unsigned char u[], unsigned char v[], int n);
void mpdiv(unsigned char q[], unsigned char r[], unsigned char u[],
	unsigned char v[], int n, int m);
void mpinv(unsigned char u[], unsigned char v[], int n, int m);
void mplsh(unsigned char u[], int n);
void mpmov(unsigned char u[], unsigned char v[], int n);
void mpmul(unsigned char w[], unsigned char u[], unsigned char v[], int n,
	int m);
void mpneg(unsigned char u[], int n);
void mppi(int n);
void mprove(double **a, double **alud, int n, int indx[], double b[],
	double x[]);
void mpsad(unsigned char w[], unsigned char u[], int n, int iv);
void mpsdv(unsigned char w[], unsigned char u[], int n, int iv, int *ir);
void mpsmu(unsigned char w[], unsigned char u[], int n, int iv);
void mpsqrt(unsigned char w[], unsigned char u[], unsigned char v[], int n,
	int m);
void mpsub(int *is, unsigned char w[], unsigned char u[], unsigned char v[],
	int n);
void mrqcof(double x[], double y[], double sig[], int ndata, double a[],
	int ia[], int ma, double **alpha, double beta[], double *chisq,
	void (*funcs)(double, double [], double *, double [], int));
void mrqmin(double x[], double y[], double sig[], int ndata, double a[],
	int ia[], int ma, double **covar, double **alpha, double *chisq,
	void (*funcs)(double, double [], double *, double [], int), double *alamda);
void newt(double x[], int n, int *check,
	void (*vecfunc)(int, double [], double []));
void odeint(double ystart[], int nvar, double x1, double x2,
	double eps, double h1, double hmin, int *nok, int *nbad,
	void (*derivs)(double, double [], double []),
	void (*rkqs)(double [], double [], int, double *, double, double,
	double [], double *, double *, void (*)(double, double [], double [])));
void orthog(int n, double anu[], double alpha[], double beta[], double a[],
	double b[]);
void pade(double cof[], int n, double *resid);
void pccheb(double d[], double c[], int n);
void pcshft(double a, double b, double d[], int n);
void pearsn(double x[], double y[], unsigned long n, double *r, double *prob,
	double *z);
void period(double x[], double y[], int n, double ofac, double hifac,
	double px[], double py[], int np, int *nout, int *jmax, double *prob);
void piksr2(int n, double arr[], double brr[]);
void piksrt(int n, double arr[]);
void pinvs(int ie1, int ie2, int je1, int jsf, int jc1, int k,
	double ***c, double **s);
double plgndr(int l, int m, double x);
double poidev(double xm, long *idum);
void polcoe(double x[], double y[], int n, double cof[]);
void polcof(double xa[], double ya[], int n, double cof[]);
void poldiv(double u[], int n, double v[], int nv, double q[], double r[]);
void polin2(double x1a[], double x2a[], double **ya, int m, int n,
	double x1, double x2, double *y, double *dy);
void polint(double xa[], double ya[], int n, double x, double *y, double *dy);
void powell(double p[], double **xi, int n, double ftol, int *iter, double *fret,
	double (*func)(double []));
void predic(double data[], int ndata, double d[], int m, double future[], int nfut);
double probks(double alam);
void psdes(unsigned long *lword, unsigned long *irword);
void pwt(double a[], unsigned long n, int isign);
void pwtset(int n);
double pythag(double a, double b);
void pzextr(int iest, double xest, double yest[], double yz[], double dy[],
	int nv);
double qgaus(double (*func)(double), double a, double b);
void qrdcmp(double **a, int n, double *c, double *d, int *sing);
double qromb(double (*func)(double), double a, double b);
double qromo(double (*func)(double), double a, double b,
	double (*choose)(double (*)(double), double, double, int));
void qroot(double p[], int n, double *b, double *c, double eps);
void qrsolv(double **a, int n, double c[], double d[], double b[]);
void qrupdt(double **r, double **qt, int n, double u[], double v[]);
double qsimp(double (*func)(double), double a, double b);
double qtrap(double (*func)(double), double a, double b);
double quad3d(double (*func)(double, double, double), double x1, double x2);
void quadct(double x, double y, double xx[], double yy[], unsigned long nn,
	double *fa, double *fb, double *fc, double *fd);
void quadmx(double **a, int n);
void quadvl(double x, double y, double *fa, double *fb, double *fc, double *fd);
double ran0(long *idum);
double ran1(long *idum);
double ran2(long *idum);
double ran3(long *idum);
double ran4(long *idum);
void rank(unsigned long n, unsigned long indx[], unsigned long irank[]);
void ranpt(double pt[], double regn[], int n);
void ratint(double xa[], double ya[], int n, double x, double *y, double *dy);
void ratlsq(double (*fn)(double), double a, double b, int mm, int kk,
	double cof[], double *dev);
double ratval(double x, double cof[], int mm, int kk);
double rc(double x, double y);
double rd(double x, double y, double z);
void realft(double data[], unsigned long n, int isign);
void rebin(double rc, int nd, double r[], double xin[], double xi[]);
void red(int iz1, int iz2, int jz1, int jz2, int jm1, int jm2, int jmf,
	int ic1, int jc1, int jcf, int kc, double ***c, double **s);
void relax(double **u, double **rhs, int n);
void relax2(double **u, double **rhs, int n);
void resid(double **res, double **u, double **rhs, int n);
double revcst(double x[], double y[], int iorder[], int ncity, int n[]);
void reverse(int iorder[], int ncity, int n[]);
double rf(double x, double y, double z);
double rj(double x, double y, double z, double p);
void rk4(double y[], double dydx[], int n, double x, double h, double yout[],
	void (*derivs)(double, double [], double []));
void rkck(double y[], double dydx[], int n, double x, double h,
	double yout[], double yerr[], void (*derivs)(double, double [], double []));
void rkdumb(double vstart[], int nvar, double x1, double x2, int nstep,
	void (*derivs)(double, double [], double []));
void rkqs(double y[], double dydx[], int n, double *x,
	double htry, double eps, double yscal[], double *hdid, double *hnext,
	void (*derivs)(double, double [], double []));
void rlft3(double ***data, double **speq, unsigned long nn1,
	unsigned long nn2, unsigned long nn3, int isign);
double rofunc(double b);
void rotate(double **r, double **qt, int n, int i, double a, double b);
void rsolv(double **a, int n, double d[], double b[]);
void rstrct(double **uc, double **uf, int nc);
double rtbis(double (*func)(double), double x1, double x2, double xacc);
double rtflsp(double (*func)(double), double x1, double x2, double xacc);
double rtnewt(void (*funcd)(double, double *, double *), double x1, double x2,
	double xacc);
double rtsafe(void (*funcd)(double, double *, double *), double x1, double x2,
	double xacc);
double rtsec(double (*func)(double), double x1, double x2, double xacc);
void rzextr(int iest, double xest, double yest[], double yz[], double dy[], int nv);
void savgol(double c[], int np, int nl, int nr, int ld, int m);
void score(double xf, double y[], double f[]);
void scrsho(double (*fx)(double));
/* double select(unsigned long k, unsigned long n, double arr[]); */
double selip(unsigned long k, unsigned long n, double arr[]);
void shell(unsigned long n, double a[]);
void shoot(int n, double v[], double f[]);
void shootf(int n, double v[], double f[]);
void simp1(double **a, int mm, int ll[], int nll, int iabf, int *kp,
	double *bmax);
void simp2(double **a, int n, int l2[], int nl2, int *ip, int kp, double *q1);
void simp3(double **a, int i1, int k1, int ip, int kp);
void simplx(double **a, int m, int n, int m1, int m2, int m3, int *icase,
	int izrov[], int iposv[]);
void simpr(double y[], double dydx[], double dfdx[], double **dfdy,
	int n, double xs, double htot, int nstep, double yout[],
	void (*derivs)(double, double [], double []));
void sinft(double y[], int n);
void slvsm2(double **u, double **rhs);
void slvsml(double **u, double **rhs);
void sncndn(double uu, double emmc, double *sn, double *cn, double *dn);
double snrm(unsigned long n, double sx[], int itol);
void sobseq(int *n, double x[]);
void solvde(int itmax, double conv, double slowc, double scalv[],
	int indexv[], int ne, int nb, int m, double **y, double ***c, double **s);
void sor(double **a, double **b, double **c, double **d, double **e,
	double **f, double **u, int jmax, double rjac);
void sort(unsigned long n, double arr[]);
void sort2(unsigned long n, double arr[], double brr[]);
void sort3(unsigned long n, double ra[], double rb[], double rc[]);
void spctrm(FILE *fp, double p[], int m, int k, int ovrlap);
void spear(double data1[], double data2[], unsigned long n, double *d, double *zd,
	double *probd, double *rs, double *probrs);
void sphbes(int n, double x, double *sj, double *sy, double *sjp, double *syp);
void splie2(double x1a[], double x2a[], double **ya, int m, int n, double **y2a);
void splin2(double x1a[], double x2a[], double **ya, double **y2a, int m, int n,
	double x1, double x2, double *y);
void spline(double x[], double y[], int n, double yp1, double ypn, double y2[]);
void splint(double xa[], double ya[], double y2a[], int n, double x, double *y);
void spread(double y, double yy[], unsigned long n, double x, int m);
void sprsax(double sa[], unsigned long ija[], double x[], double b[],
	unsigned long n);
void sprsin(double **a, int n, double thresh, unsigned long nmax, double sa[],
	unsigned long ija[]);
void sprspm(double sa[], unsigned long ija[], double sb[], unsigned long ijb[],
	double sc[], unsigned long ijc[]);
void sprstm(double sa[], unsigned long ija[], double sb[], unsigned long ijb[],
	double thresh, unsigned long nmax, double sc[], unsigned long ijc[]);
void sprstp(double sa[], unsigned long ija[], double sb[], unsigned long ijb[]);
void sprstx(double sa[], unsigned long ija[], double x[], double b[],
	unsigned long n);
void stifbs(double y[], double dydx[], int nv, double *xx,
	double htry, double eps, double yscal[], double *hdid, double *hnext,
	void (*derivs)(double, double [], double []));
void stiff(double y[], double dydx[], int n, double *x,
	double htry, double eps, double yscal[], double *hdid, double *hnext,
	void (*derivs)(double, double [], double []));
void stoerm(double y[], double d2y[], int nv, double xs,
	double htot, int nstep, double yout[],
	void (*derivs)(double, double [], double []));
void svbksb(double **u, double w[], double **v, int m, int n, double b[],
	double x[]);
void svdcmp(double **a, int m, int n, double w[], double **v);
void svdfit(double x[], double y[], double sig[], int ndata, double a[],
	int ma, double **u, double **v, double w[], double *chisq,
	void (*funcs)(double, double [], int));
void svdvar(double **v, int ma, double w[], double **cvm);
void toeplz(double r[], double x[], double y[], int n);
void tptest(double data1[], double data2[], unsigned long n, double *t, double *prob);
void tqli(double d[], double e[], int n, double **z);
double trapzd(double (*func)(double), double a, double b, int n);
void tred2(double **a, int n, double d[], double e[]);
void tridag(double a[], double b[], double c[], double r[], double u[],
	unsigned long n);
double trncst(double x[], double y[], int iorder[], int ncity, int n[]);
void trnspt(int iorder[], int ncity, int n[]);
void ttest(double data1[], unsigned long n1, double data2[], unsigned long n2,
	double *t, double *prob);
void tutest(double data1[], unsigned long n1, double data2[], unsigned long n2,
	double *t, double *prob);
void twofft(double data1[], double data2[], double fft1[], double fft2[],
	unsigned long n);
void vander(double x[], double w[], double q[], int n);
void vegas(double regn[], int ndim, double (*fxn)(double [], double), int init,
	unsigned long ncall, int itmx, int nprn, double *tgral, double *sd,
	double *chi2a);
void voltra(int n, int m, double t0, double h, double *t, double **f,
	double (*g)(int, double), double (*ak)(int, int, double, double));
void wt1(double a[], unsigned long n, int isign,
	void (*wtstep)(double [], unsigned long, int));
void wtn(double a[], unsigned long nn[], int ndim, int isign,
	void (*wtstep)(double [], unsigned long, int));
void wwghts(double wghts[], int n, double h,
	void (*kermom)(double [], double ,int));
int zbrac(double (*func)(double), double *x1, double *x2);
void zbrak(double (*fx)(double), double x1, double x2, int n, double xb1[],
	double xb2[], int *nb);
double zbrent(double (*func)(double), double x1, double x2, double tol);
void zrhqr(double a[], int m, double rtr[], double rti[]);
double zriddr(double (*func)(double), double x1, double x2, double xacc);
void zroots(fcomplex a[], int m, fcomplex roots[], int polish);



void rlft3dbl(double ***data, double **speq, unsigned long nn1,
	unsigned long nn2, unsigned long nn3, int isign);
void fourndbl(double data[], unsigned long nn[], int ndim, int isign);


void linmindbl(double p[], double xi[], int n, double *fret, double (*func)(double []));
void frprmndbl(double p[], int n, double ftol, int *iter, double *fret,
	double (*func)(double []), void (*dfunc)(double [], double []));
void mnbrakdbl(double *ax, double *bx, double *cx, double *fa, double *fb, double *fc,
	double (*func)(double));
double f1dimdbl(double x);
double brentdbl(double ax, double bx, double cx, double (*f)(double), double tol,
	double *xmin);



#endif /* _NR_H_ */
