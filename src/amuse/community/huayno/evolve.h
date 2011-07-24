#define FLOAT double
#define CLFLOAT cl_double
#define CLFLOAT4 cl_double4
#define DOUBLE long double
#define INT int
#define UINT unsigned int
#define LONG long
#define ULONG unsigned long

#define SWAP(a,b,c) {c t;t=(a);(a)=(b);(b)=t;}

#define ABS(X) (((X) >= 0) ? (X) : -(X))

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

struct particle
{
  UINT id;
  FLOAT mass;
  FLOAT radius; /*not used*/
  DOUBLE pos[3];
  DOUBLE vel[3];
  DOUBLE pot;
  DOUBLE postime;
  FLOAT timestep;
  INT level; 
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
  UNSPLIT,    // 1
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
};

extern FLOAT eps2;
extern FLOAT dt_param;

void init_code();
void stop_code();
void init_evolve(struct sys s);
void do_evolve(struct sys s, double dt, int inttype);
FLOAT system_potential_energy(struct sys s);
FLOAT system_kinetic_energy(struct sys s);

#define RVTIMESTEP
#define RATIMESTEP
#define RARVRATIO   0.25

#define MPWORKLIMIT 100
#define CLWORKLIMIT 10000

#define MAXLEVEL  64

extern FLOAT eps2;
extern FLOAT dt_param;
extern struct sys zerosys;

/* diagnostics */
extern DOUBLE simtime;
extern int clevel;
extern unsigned long tcount[MAXLEVEL],kcount[MAXLEVEL],dcount[MAXLEVEL],deepsteps;
extern unsigned long tstep[MAXLEVEL],kstep[MAXLEVEL],dstep[MAXLEVEL];
extern unsigned long cefail[MAXLEVEL],cecount[MAXLEVEL]; // call/fail counts of the Kepler solver
#ifdef EVOLVE_OPENCL
extern unsigned long cpu_step,cl_step,cpu_count,cl_count;
#endif

/*static*/ void drift(struct sys s, DOUBLE etime, DOUBLE dt); /* drift sys */
/*static*/ void kick(struct sys s1, struct sys s2, DOUBLE dt); /* =kick sys1 for interactions with sys2  */

void evolve_cc2(struct sys s, DOUBLE stime, DOUBLE etime, DOUBLE dt);


// Optimal Kick header
struct force {
    struct particle *parti;
    struct particle *partj;
    FLOAT timestep;
};
struct forces {
    UINT n;
    struct force *forc;
    struct force *last;
};
extern struct forces zeroforces;

void evolve_ok_init(struct sys s);
void evolve_ok_stop();
void evolve_ok2(struct sys s, struct forces f, DOUBLE stime, DOUBLE etime, DOUBLE dt, int calc_timestep);
