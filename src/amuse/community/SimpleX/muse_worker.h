extern "C" {
  
  
  int evolve(double time_end, int synchronize);
  
  
  int remove_particle(int id);
  
  
  int cleanup_module();
  
  
  int initialize(double current_time);
  
  
  void get_state(int id, double * x, double * y, double * z, double * rho, double * flux, double * xion);
  
  
  int reinitialize();
  
  
  int add_particle(int * id, double x, double y, double z, double rho, double flux, double xion);
  
  
  int setup_module();
  
}
