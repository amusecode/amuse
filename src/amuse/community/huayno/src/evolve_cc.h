//~ #define CC2_SPLIT_SHORTCUTS
void evolve_cc2(int clevel,struct sys s, DOUBLE stime, DOUBLE etime, DOUBLE dt, int inttype,int recenter);

#ifdef CC2_SPLIT_SHORTCUTS
void evolve_cc2_shortcut(int clevel,struct sys s, DOUBLE stime, DOUBLE etime, DOUBLE dt, int inttype, int recenter, FLOAT dtsys);
#endif
