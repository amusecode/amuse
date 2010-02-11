#ifndef _G6LIB_
#define _G6LIB_

#ifdef __cplusplus
extern "C" {
#endif

int g6_open_(int *id);
int g6_close_(int *id);
int g6_npipes_();
int g6_set_tunit_(int*);
int g6_set_xunit_(int*);
int g6_set_ti_(int *id, double *ti);
int g6_set_j_particle_(int *cluster_id,
         int *address,
         int *index,
         double *tj, double *dtj,
         double *mass,
         double k18[3], double j6[3],
         double a2[3], double v[3], double x[3]);
void g6calc_firsthalf_(int *cluster_id,
         int *nj, int *ni,
         int index[], 
         double xi[][3], double vi[][3],
         double aold[][3], double j6old[][3],
         double phiold[3], 
         double *eps2, double h2[]);
int g6calc_lasthalf_(int *cluster_id,
           int *nj, int *ni,
           int index[], 
           double xi[][3], double vi[][3],
           double *eps2, double h2[],
           double acc[][3], double jerk[][3], double pot[]);
int g6calc_lasthalf2_(int *cluster_id,
        int *nj, int *ni,
        int index[], 
        double xi[][3], double vi[][3],
        double *eps2, double h2[],
        double acc[][3], double jerk[][3], double pot[],
        int *inn);

int g6_initialize_jp_buffer_(int* cluster_id, int* buf_size);
int g6_flush_jp_buffer_(int* cluster_id);
int g6_reset_(int* cluster_id);
int g6_reset_fofpga_(int* cluster_id);

int g6_read_neighbour_list_(int* cluster_id);

int g6_get_neighbour_list_(int *cluster_id,
             int *ipipe,
             int *maxlength,
             int *n_neighbours,
             int neighbour_list[]);
             
             
             
             
int g6_open(int clusterid);
int g6_close(int clusterid);

int g6_npipes();

int g6_set_tunit(int newtunit);
int g6_set_xunit(int newxunit);

//void g6_set_njp(int clusterid, int njp);

int g6_set_ti(int clusterid, double ti);
int g6_set_j_particle(int clusterid, int address,
         int index,
         double tj, /* particle time */
         double dtj, /* particle time */
         double mass,
         double a2by18[3], /* a2dot divided by 18 */
         double a1by6[3], /* a1dot divided by 6 */
         double aby2[3], /* a divided by 2 */
         double v[3], /* velocity */
         double x[3] /* position */);
         
int g6_read_neighbour_list(int clusterid);
void g6_set_neighbour_list_sort_mode(int mode);;
int g6_get_neighbour_list_sort_mode();
int g6_get_neighbour_list(int clusterid,
			   int ipipe,
			   int maxlength,
			   int *nblen,
			   int nbl[]);

int g6_get_number_of_pipelines();
void g6_reset(int devid);

void g6calc_firsthalf(int clusterid, 
                      int nj,  
                      int ni,  
                      int index[],  
                      double xi[][3],  
                      double vi[][3],  
                      double fold[][3],
                      double jold[][3],  
                      double phiold[],  
                      double eps2,   
                      double h2[]);   
int g6calc_lasthalf(int clusterid,
                     int nj,
                     int ni,
                     int index[],
                     double xi[][3],
                     double vi[][3],
                     double eps2,
                     double h2[], 
                     double acc[][3],
                     double jerk[][3],
                     double pot[]);   
int g6calc_lasthalf2(int clusterid,
                     int nj,
                     int ni,
                     int index[],
                     double xi[][3],
                     double vi[][3],
                     double eps2,
                     double h2[], 
                     double acc[][3],
                     double jerk[][3],
                     double pot[],     
                      int nnbindex[]);
void g6_reinitialize(int clusterid);
int g6_initialize_jp_buffer(int clusterid, int size);
int g6_flush_jp_buffer(int  clusterid);
#ifdef __cplusplus
 }
#endif 

#endif
