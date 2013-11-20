

double set_halo_positions(void);
double set_disk_positions(void);
double set_bulge_positions(void);

double set_halo_velocities(void);
double set_disk_velocities(void);
double set_bulge_velocities(void);




void save_particles(char *fname);




double comp_Dphi_z_halo(double R,double z);
double comp_Dphi_R_halo(double R,double z);
double comp_rho_halo(double R,double z);

double comp_Dphi_R_disk_razorthin(double RR,double zz);

double comp_Dphi_z_bulge(double R,double z);
double comp_Dphi_R_bulge(double R,double z);
double comp_rho_bulge(double R,double z);


double comp_Dphi_z_disk(double R,double z);
double comp_Dphi_R_disk(double R,double z);
double comp_rho_disk(double R,double z);


double comp_Dphi_z(double R,double z);
double comp_Dphi_R(double R,double z);

double epicyclic_kappa2(double R);


double splint_zs_y_D2y(double z);
double splint_xl_yl_D2yl(double t);


double   drand48(void);



double mass_cumulative_disk(double);
double mass_cumulative_bulge(double);
