#define CLCONTEXT   stdcpu
#define NTHREAD  64
#define BLOCKSIZE 64

void init_cl();
void close_cl();
void kick_cl(struct sys s1, struct sys s2, DOUBLE dt);
void potential_cl(struct sys s1, struct sys s2);
void timestep_cl(struct sys s1, struct sys s2);

