#include "Main.h"
#include "configfile.h"

class AMUSE_SimpleX : public SimpleX {

  friend class SimpleX;
  
 public:
  
 AMUSE_SimpleX(string output_path): SimpleX(output_path) {
    total_time=0.;
    syncflag=-1;
  };
  
  ~AMUSE_SimpleX(){};
  
  int add_vertex(long *id, double x,double y,double z,double rho,
		 double flux,double xion, double);
  int add_site(long *id, double x,double y,double z,double rho,
	       double flux,double xion,double);
  int remove_site(long id);
  int setup_simplex();

  void convert_units();
   
  int initialize(double current_time);
  int reinitialize();
  int evolve(double dt, int sync);
  int get_site(int id, double *x,double *y,double *z,double *rho,
	       double *flux,double *xion, double *uInt);
  int get_position(int id, double *x, double *y, double *z);
  int get_density(int id, double *rho);
  int get_flux(int id, double *flux);
  int get_ionisation(int id, double *xion);
  int get_internalEnergy(int id, double *uInt);
  int get_dinternalEnergydt(int id, double *dudt);
  int set_site(int id, double x, double y, double z, double rho,
	       double flux, double xion, double uInt);
  int set_position(int id, double x, double y, double z);
  int set_density(int id, double rho);
  int set_flux(int id, double flux);
  int set_ionisation(int id, double xion);
  int set_internalEnergy(int id, double uInt);
  int set_dinternalEnergydt(int id, double dudt);
  int cleanup();
  int get_syncflag(){return syncflag;}
  int set_UNIT_T(double ts){ UNIT_T_MYR=ts;return 0;}
  int get_UNIT_T(double *ts){ *ts=UNIT_T_MYR;return 0;}
  int set_sizeBox(double bs){ sizeBox=bs;return 0;}
  int get_sizeBox(double *bs){ *bs=sizeBox;return 0;}
  int set_hilbert_order(int ho){ hilbert_order=ho;return 0;}
  int get_hilbert_order(int *ho){ *ho=hilbert_order;return 0;}

  int set_sourceTeff(double ts){ sourceTeff=ts;return 0;}
  int get_sourceTeff(double *ts){ *ts=sourceTeff;return 0;}
  int set_numFreq(int ts){ numFreq=ts;return 0;}
  int get_numFreq(int *ts){ *ts=numFreq;return 0;}
  int set_heat_cool(int ts){ heat_cool=ts;return 0;}
  int get_heat_cool(int *ts){ *ts=heat_cool;return 0;}
  int set_metal_cooling(int ts){ metal_cooling=ts;return 0;}
  int get_metal_cooling(int *ts){ *ts=metal_cooling;return 0;}
  int set_coll_ion(int ts){ coll_ion=ts;return 0;}
  int get_coll_ion(int *ts){ *ts=heat_cool;return 0;}
  int set_blackBody(int ts){ blackBody=ts;return 0;}
  int get_blackBody(int *ts){ *ts=blackBody;return 0;}



 private:

  void read_parameters();
  double total_time;     
  int syncflag; 
};
