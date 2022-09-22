#define FLOAT double
#define CLFLOAT cl_double
#define CLFLOAT4 cl_double4
#define DOUBLE long double
#define INT int
#define UINT unsigned int
#define LONG long
#define ULONG unsigned long

#define RVTIMESTEP
#define RATIMESTEP
#define RARVRATIO   1.

#define MPWORKLIMIT 1000
#define CLWORKLIMIT 40000

#define MAXLEVEL  64

#define COMPENSATED_SUMMP
//~ #define COMPENSATED_SUMMV  

//~ #define CONSISTENCY_CHECKS // perform (time-consuming, but thorough) sanity checks

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
  UINT n, nzero; // n=total particles, nzero=# zero mass particles
  struct particle *part; // start of particles, NULL iff n==0 
  struct particle *zeropart; // start of zero mass particles nzero>0, otherwise NULL
};

#define GETPART(s, i)   ((i)<(s).n-(s).nzero ? (s).part+(i) : (s).zeropart+((i)-((s).n-(s).nzero)))    
#define LAST(s)   ((s).part==NULL || (s).n-(s).nzero==0 ? NULL : (s).part+((s).n-(s).nzero)-1)
#define LASTZERO(s)   ((s).zeropart==NULL || (s).nzero==0 ? NULL : (s).zeropart+(s).nzero-1)

extern struct sys debugsys; // for monitoring purposes 

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
  SHARED2_COLLISIONS, // 29
  SHARED4_COLLISIONS, // 30
  SHARED6_COLLISIONS, // 31
  SHARED8_COLLISIONS, // 32
  SHARED10_COLLISIONS, // 33
  CONSTANT2,   // 34
  CONSTANT4,   // 35
  CONSTANT6,   // 36
  CONSTANT8,   // 37
  CONSTANT10,   // 38
  ERROR_CONTROL, // 39
  CC_SHARED10, // 40
  CCC_SHARED10 // 41
};

extern int verbosity;
extern FLOAT eps2;
extern FLOAT dt_param;
#pragma omp threadprivate(dt_param)
extern int accel_zero_mass;

extern struct sys zerosys;

extern int fixed_j;
extern DOUBLE bs_target_error;
extern int opencl_device_type;

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
DOUBLE system_potential_energy(struct sys s);
DOUBLE system_kinetic_energy(struct sys s);

void drift(int clevel,struct sys s, DOUBLE etime, DOUBLE dt); /* drift sys */
void kick(int clevel,struct sys s1, struct sys s2, DOUBLE dt); /* =kick sys1 for interactions with sys2  */

void kdk(int clevel,struct sys s1,struct sys s2, DOUBLE stime, DOUBLE etime, DOUBLE dt);
void dkd(int clevel,struct sys s1,struct sys s2, DOUBLE stime, DOUBLE etime, DOUBLE dt);

void timestep(int clevel,struct sys s1, struct sys s2,int dir);
FLOAT timestep_ij(struct particle *i, struct particle *j,int dir);
FLOAT global_timestep(struct sys s);
FLOAT max_global_timestep(struct sys s);

void potential(struct sys s1, struct sys s2);

struct sys join(struct sys s1,struct sys s2);
void split_zeromass(struct sys *s);
void verify_split_zeromass(struct sys s);

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

#define CHECK_TIMESTEP(etime,stime,dt,clevel) \
  if(sizeof(dt)==sizeof(long double)) { \
  if(etime == stime ||  dt==0 || clevel>=MAXLEVEL) \
    ENDRUN("timestep too small: etime=%Le stime=%Le dt=%Le clevel=%u\n", etime, stime, dt, clevel); \
  } else { \
  if(etime == stime ||  dt==0 || clevel>=MAXLEVEL) \
    ENDRUN("timestep too small: etime=%le stime=%le dt=%le clevel=%u\n", (double) etime, (double) stime, (double) dt, clevel); \
  }  
