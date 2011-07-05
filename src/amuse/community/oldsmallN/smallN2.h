#include "hdyn.h"

// in integrate.cc

bool integrate(hdyn *b,
	       real t_end, real dt_struc,
	       real dt_dyn = 1,
	       real eps2 = 0,
	       bool verbose = true,
	       real break_r2 = VERY_LARGE_NUMBER,
	       real dt_log = 0,
	       real dt_energy = 0,
	       real dt_snap = 0);

// in analyze.cc

dyn *flat_copy_tree(dyn *bin);
bool check_structure(hdyn *bin, real eps2 = 0, bool verbose = true);
dyn* get_tree(hdyn *bin, real eps2 = 0, bool verbose = false);
void my_sys_stats(dyn *b, int verbose = 1);
