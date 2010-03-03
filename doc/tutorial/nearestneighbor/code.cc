#include "code.h"
#include <math.h>

/**
 * return the distance between point0 and point1
 */
double distance_between_points(
    double x0, double y0, double z0, 
    double x1, double y1, double z1)
{
    double dx = x1 - x0;
    double dy = y1 - y0;
    double dz = z1 - z0;
    return sqrt(dx * dx  + dy * dy + dz * dz);
}

class DistanceAndIndex{
    
public:
    double r;
    int index;
    
    DistanceAndIndex():r(-1),index(-1) {}
    DistanceAndIndex(DistanceAndIndex & original):r(original.r),index(original.index){}
};

/**
 * Find the nearest neigbors of all the points (specified with
 * x, y and z).
 * 
 * Fills the n1, n2 and n3 arrays with the closest, second clostest
 * and third closests points.
 * 
 **/
int find_nearest_neighbors(int npoints, 
    double * x, double * y, double * z,
    int * n0, int * n1, int * n2)
{
    for(int i = 0; i < npoints; i++)
    {
        DistanceAndIndex found[3];
        
        for(int j = 0; j < npoints; j++)
        {
            if(j == i) {
                continue;
            }
            
            double x0 = x[i];
            double y0 = x[i];
            double z0 = x[i];
            
            double x1 = x[j];
            double y1 = x[j];
            double z1 = x[j];
            
            double r = distance_between_points(x0, y0, z0, x1, y1 , z1);
            
            for(int k = 0; k < 3; k++) {
                if(found[k].index == -1) {
                    found[k].index = j;
                    found[k].r = r;
                    break;
                } else {
                    if(r < found[k].r) {
                        for(int l = 2;l > k; l--) {
                            found[l] = found[l-1];
                        }
                        found[k].index = j;
                        found[k].r = r;
                        break;
                    }
                }
            }
            
        }
        
        n0[i] = found[0].index;
        n1[i] = found[1].index;
        n2[i] = found[2].index;
    }
    
    
    return 0;
}
