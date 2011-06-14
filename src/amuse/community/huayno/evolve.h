#define FLOAT double
#define CLFLOAT cl_double
#define CLFLOAT4 cl_double4
#define DOUBLE long double
#define INT int
#define UINT unsigned int
#define LONG long
#define ULONG unsigned long

#define SWAP(a,b,c) {c t;t=(a);(a)=(b);(b)=t;}

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
};

enum intopt
{
  CONSTANT,  // 0
  UNSPLIT,   // 1
  PASS,      // 2
  HOLD,      // 3
  BRIDGE,    // 4
  NAIVE,     // 5
  VARIABLE,  // 6
  PASS_DKD,  // 7
  HOLD_DKD,  // 8
  PPASS_DKD, // 9
  BRIDGE_DKD // 10
};

extern FLOAT eps2;
extern FLOAT dt_param;

void init_code();
void stop_code();
void init_evolve(struct sys s);
void do_evolve(struct sys s, double dt, int inttype);
FLOAT system_potential_energy(struct sys s);
FLOAT system_kinetic_energy(struct sys s);
