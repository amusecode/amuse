#include<fstream>
#include<vector>

#include "integrator.h"


typedef double  real;

//-------------------------------------------------------------------------
//
/// Structure defining particle dynamics data.

//  Use components to avoid possible SWIG problems with vectors.

typedef struct {
    int id;                                             /// identifier
    double mass;                                        /// mass
    double radius;                                      /// radius
    double x, y, z;                                     /// position
    double vx, vy, vz;                                  /// velocity
} dynamics_state;



#define PR(x)  cerr << #x << " = " << x << " "
#define PRC(x) cerr << #x << " = " << x << ",  "
#define PRL(x) cerr << #x << " = " << x << "\n"

int get_index_from_identity(int id);
