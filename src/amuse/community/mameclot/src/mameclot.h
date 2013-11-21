
typedef struct input{
  int N;
  int N2;
  int Ncl;
  int model;
  int imftype; // 0=equal mass; 1=Kroupa (2001)
  int spin;
  float frot;
  double ra;
  char name[15];
  int gamma;
  unsigned int seed;
  double q;
  double rcut;
  double rbar;
  double eta; // r_2/r_1 = q^eta
  double d;
  double Ehat;  double Lhat;
}INPUT;

typedef struct star{
  double mass;
  double kin, phi;
  double pos[3];
  double vel[3];
} STAR;

typedef struct cluster{
  int N;
  int model; // 0=Cored gamma-model; 1=Hernquist; 2=Jaffe; 3=Isochrone; 4=Plummer
  int id;
  char name[15];
  int gamma;
  int imftype;
  int spin;
  float frot;
  float ra;
  float rh_over_r0;
  float rv_over_r0;
  float rh_over_rv;
  float rmax_over_r0;
  float E;
  float W;
  float K;
  float Lz;
  float Krot;
  float lambda;
  double vrms; 
  double rvir; 
  double rcut; 
  double M; 
  double mmean; //Msun
  double trh; 
  double compos[3];
  double comvel[3];
  STAR *stars;
} CLUSTER;

typedef struct system{
  int N;
  int Ncl;
  int seed;
  double q;
  double eta; // m-r relation
  double d;
  double b;
  double dx;
  double tstar;
  double vstar;
  double mstar;
  double rstar; 
  double vrel;
  double Lz;
  double lambda;
  double lambda_orb;
  double tenc;
  double mu;
  double Ehat;
  double Lhat;
  double mrsig;
  double msig2;
  double rfac;
  double vfac;
  CLUSTER *clusters;
} SYSTEM;

double myrand();
void get_args(int argc, char* argv[], INPUT *parameters);
void parameter_use();
void parameter_check(INPUT *parameters);
void allocate_memory(SYSTEM **system, INPUT *parameters);
void initialize(SYSTEM **system,  INPUT parameters);
void create(SYSTEM **system);
void set_scalings(CLUSTER *cluster);
void imf(CLUSTER *cluster);
void get_pos_vel(CLUSTER *cluster);
double get_r(CLUSTER *cluster);
double get_v(double r, int model);
void get_osipkov_merrit_v_eta(double r, double ra,int model, double *v, double *eta);
void get_q_eta_pair(double p, double *q, double *eta);
void scale(CLUSTER *cluster);
void scale_system(SYSTEM *system);
double Lz(SYSTEM *system);
extern void calculate_potential();
extern void calculate_cross_potential();
void twobody_orbit(SYSTEM *system);
double dawson(double x);
void shell(double a[], int n);
void output(SYSTEM *system);
void free_memory(SYSTEM *system);
