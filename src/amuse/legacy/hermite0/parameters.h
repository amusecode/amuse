/* Time-stamp: <6 September 2006 16:32:45 CEST bon@science.uva.nl>

   We duplicate the typing information from muse_dynamics.cc here.
   (It should probably be in the header file anyway.)
   Swig needs this information to make these globals available to python so
   that MUSE can access them as, e.g., md.dt_param = 4.0
*/

typedef double  real;

// List here any global variables we wish to make accessible to python
// via SWIG.  See muse_dynamics.py for more details on how this occurs.

extern real t;
extern real dt_param;
extern real dt_dia;
extern real eps2;
extern bool flag_collision;
