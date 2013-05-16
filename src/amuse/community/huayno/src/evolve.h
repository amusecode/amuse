#define FLOAT double
#define CLFLOAT cl_float
#define CLFLOAT4 cl_float4
#define DOUBLE long double
#define INT int
#define UINT unsigned int
#define LONG long
#define ULONG unsigned long

#define RVTIMESTEP
#define RATIMESTEP
#define RARVRATIO   1.

#define MPWORKLIMIT 1000
#define CLWORKLIMIT 100000

#define MAXLEVEL  64

#define COMPENSATED_SUMMP
//#define COMPENSATED_SUMMV  

struct particle
{
  UINT id;
  FLOAT mass;
  FLOAT radius; /*only used for stopping conditions*/
  DOUBLE pos[3];
  DOUBLE vel[3];
#ifdef COMPENSATED_SUMMP
  DOUBLE pos_e[3];
#endif
#ifdef COMPENSATED_SUMMV
  DOUBLE vel_e[3];
#endif
  DOUBLE pot;
  DOUBLE postime;
  FLOAT timestep;
};

struct jparticle
{
  FLOAT mass;
  FLOAT pos[3];
  FLOAT vel[3];
};

struct sys
{
  UINT n; 
  struct particle *part;
  struct particle *last;
  struct sys *next_cc; // used in the CC split only
};

enum intopt
{
  CONSTANT,   // 0
  SHARED2,    // 1
  PASS,       // 2
  HOLD,       // 3
  BRIDGE,     // 4
  NAIVE,      // 5
  VARIABLE,   // 6
  PASS_DKD,   // 7
  HOLD_DKD,   // 8
  PPASS_DKD,  // 9
  BRIDGE_DKD, // 10
  CC,         // 11
  CC_KEPLER,  // 12
  OK,         // 13
  KEPLER,     // 14
  SHARED4,    // 15
  FOURTH_M4,  // 16 
  FOURTH_M5,  // 17 
  SHARED6,    // 18
  SHARED8,    // 19
  SHARED10,   // 20
  SHAREDBS,    // 21
  CCC,         // 22
  CCC_KEPLER,  // 23
  CC_BS,       // 24
  CCC_BS,       // 25
  BS_CC_KEPLER,    // 26
  CC_BSA,       // 27
  CCC_BSA,       // 28
  SHARED6_COLLISIONS, // 29
};

extern FLOAT eps2;
extern FLOAT dt_param;
#pragma omp threadprivate(dt_param)

extern struct sys zerosys;

extern int fixed_j;
extern DOUBLE bs_target_error;

/* diagnostics */
struct diagnostics {
  DOUBLE simtime;
  DOUBLE timetrack;
  unsigned long deepsteps;
  unsigned long tcount[MAXLEVEL],kcount[MAXLEVEL],dcount[MAXLEVEL];
  unsigned long tstep[MAXLEVEL],kstep[MAXLEVEL],dstep[MAXLEVEL];
  unsigned long cefail[MAXLEVEL],cecount[MAXLEVEL]; // call/fail counts of the Kepler solver
  unsigned long bsstep[MAXLEVEL],jcount[MAXLEVEL]; /* count + jcount of BS evolve */
#ifdef EVOLVE_OPENCL
  unsigned long cpu_step,cl_step,cpu_count,cl_count;
#endif
  int ntasks[MAXLEVEL],taskcount[MAXLEVEL];
  unsigned long taskdrift,taskkick;
};

extern struct diagnostics global_diag;
extern struct diagnostics *diag;
#pragma omp threadprivate(diag)


void init_code();
void stop_code();
void init_evolve(struct sys s, int inttype);
void do_evolve(struct sys s, double dt, int inttype);

void system_center_of_mass(struct sys s, DOUBLE *cmpos, DOUBLE *cmvel);
void move_system(struct sys s, DOUBLE dpos[3],DOUBLE dvel[3],int dir);
FLOAT system_potential_energy(struct sys s);
FLOAT system_kinetic_energy(struct sys s);

void drift(int clevel,struct sys s, DOUBLE etime, DOUBLE dt); /* drift sys */
void kick(int clevel,struct sys s1, struct sys s2, DOUBLE dt); /* =kick sys1 for interactions with sys2  */

void kdk(int clevel,struct sys s1,struct sys s2, DOUBLE stime, DOUBLE etime, DOUBLE dt);
void dkd(int clevel,struct sys s1,struct sys s2, DOUBLE stime, DOUBLE etime, DOUBLE dt);

void timestep(int clevel,struct sys s1, struct sys s2,int dir);
FLOAT timestep_ij(struct particle *i, struct particle *j,int dir);
FLOAT global_timestep(struct sys s);

void potential(struct sys s1, struct sys s2);

struct sys join(struct sys s1,struct sys s2);

#define SWAP(a,b,c) {c t;t=(a);(a)=(b);(b)=t;}

#define ABS(X) (((X) >= 0) ? (X) : -(X))
#define SIGN(X)   ((X>0)-(X<0))

#define LOG(fmt, ...) {\
  printf("%s:%d\t", __FILE__, __LINE__);\
  printf(fmt, ## __VA_ARGS__);\
}

#define ENDRUN(fmt, ...) { \
  printf("ENDRUN at %s:%d ", __FILE__, __LINE__);\
  printf(fmt, ## __VA_ARGS__);\
  fflush(stdout);\
  exit(-1);\
}

#ifdef COMPENSATED_SUMMP
#define COMPSUMP(sum,err,delta) \
  { \
    DOUBLE a; \
    a=sum; \
    err=err+delta; \
    sum=a+err; \
    err=err+(a-sum); \
  }
#else
#define COMPSUMP(sum,err,delta)  {sum+=delta;}
#endif

#ifdef COMPENSATED_SUMMV
#define COMPSUMV(sum,err,delta) \
  { \
    DOUBLE a; \
    a=sum; \
    err=err+delta; \
    sum=a+err; \
    err=err+(a-sum); \
  }
#else
#define COMPSUMV(sum,err,delta)  {sum+=delta;}
#endif

#define COMPSUM(sum,err,delta) \
  { \
    DOUBLE a; \
    a=sum; \
    err=err+delta; \
    sum=a+err; \
    err=err+(a-sum); \
  }

#define COMPSUM1(sum,err,delta) \
  { \
    DOUBLE t,y; \
    y=(delta)-err; \
    t=sum+y; \
    err=(t-sum)-y; \
    sum=t; \
  }
