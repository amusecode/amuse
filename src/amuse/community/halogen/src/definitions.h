/* 
** definitions.h
**
** Some definitions for halogen4muse
*/

#ifdef _WIN32
double drand48(void);
void srand48(long seed);
#endif

#define NGRIDR 2001
#define NGRIDDF 101
#define NINTMIN 5
#define NINTMAX 28
#define NINTMINDF 8
#define NINTMAXDF 28
#define TOL 1e-10
#define TOLDF 1e-3
#define TOLLININT 1e-10
#define DFFAILUREMAX 1e20
#define FACTORRINNER 1e-6
#define FACTORROUTER 1e20
#define SBI 1e100
#define CutoffFac 0.3
#define G 1
#define STRINGSIZE 50
#define INT int
#define FLOAT float
#define DOUBLE double
#define CHAR char
#define OFI1 "%d"
#define OFI2 "%-3d"
#define OFI3 "%-10d"
#define OFD1 "%g"
#define OFD2 "%8.7e"
#define OFD3 "%8.7e"
#define OFD4 "%+8.7e"
#define OFD5 "%16.15e"
#define OFD6 "%+16.15e"

typedef struct gridr {

    DOUBLE r[NGRIDR];
    DOUBLE logr[NGRIDR];
    DOUBLE rho[NGRIDR];
    DOUBLE logrho[NGRIDR];
    DOUBLE rhoHalo[NGRIDR];
    DOUBLE logrhoHalo[NGRIDR];
    DOUBLE rhoenc[NGRIDR];
    DOUBLE logrhoenc[NGRIDR];
    DOUBLE rhoencHalo[NGRIDR];
    DOUBLE logrhoencHalo[NGRIDR];
    DOUBLE Menc[NGRIDR];
    DOUBLE logMenc[NGRIDR];
    DOUBLE MencHalo[NGRIDR];
    DOUBLE logMencHalo[NGRIDR];
    DOUBLE Pot[NGRIDR];
    DOUBLE logPot[NGRIDR];
    DOUBLE Potoutr[NGRIDR];
    DOUBLE eqrvcmax[NGRIDR];
    } GRIDR;

typedef struct griddf {

    DOUBLE r[NGRIDDF];
    DOUBLE logr[NGRIDDF];
    DOUBLE E[NGRIDDF];
    DOUBLE logE[NGRIDDF];
    DOUBLE fE[NGRIDDF];
    DOUBLE logfE[NGRIDDF];
    } GRIDDF;

typedef struct systemparameters {

    DOUBLE alpha;
    DOUBLE beta;
    DOUBLE gamma;
    DOUBLE delta;
    DOUBLE M;
    DOUBLE rho0;
    DOUBLE rs;
    DOUBLE rvir;
    DOUBLE rcutoff;
    DOUBLE rdecay;
    DOUBLE rhalf;
    } SP;
    
typedef struct particle {

    DOUBLE r[4];
    DOUBLE v[4];
    DOUBLE mass;
    } PARTICLE;

typedef struct stuff {

    INT N;
    DOUBLE Mp;
    DOUBLE Ekin;
    DOUBLE Epot;
    DOUBLE Etot;
    DOUBLE Cr[4];
    DOUBLE Cv[4];
    } STUFF;

typedef struct systeminfo {

    INT N0;
    INT N;
    DOUBLE Mmin;
    DOUBLE Mmax;
    DOUBLE mass;
    DOUBLE rimp;
    DOUBLE r1;
    DOUBLE r100;
    DOUBLE logr[NGRIDR];
    DOUBLE logMenc[NGRIDR];
    DOUBLE logrhoenc[NGRIDR];
    SP *sp;
    PARTICLE *p;
    GRIDDF *griddf;
    CHAR systemname[STRINGSIZE];
    } SI;

typedef struct generalinfo {

    INT dorvirexact;
    DOUBLE rinner;
    DOUBLE router;
    STUFF *stuff;
    GRIDR *gridr;
    } GI;

