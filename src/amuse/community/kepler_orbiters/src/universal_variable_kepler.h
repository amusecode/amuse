#define DOUBLE double

#define COMPARE(x, y) (((x) > (y)) - ((x) < (y)))
#define SIGN(x) COMPARE(x, 0)

int universal_variable_kepler_solver(DOUBLE dt,DOUBLE mu,DOUBLE pos0[3], 
             DOUBLE vel0[3],DOUBLE pos[3], DOUBLE vel[3]);
