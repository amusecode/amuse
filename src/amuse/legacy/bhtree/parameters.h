/* Time-stamp: <21 June 2007 14:15:10 CEST bon@science.uva.nl>

   Swig needs this information to make these globals available to python so
   that MUSE can access them as, e.g., md.timestep = 0.015625
*/


extern real timestep;
extern real eps2_for_gravity;
extern real theta_for_tree;
extern int  use_self_gravity;
extern int  ncrit_for_tree;
extern real dt_dia;
extern BHTC_SYSTEM bhtcs;
