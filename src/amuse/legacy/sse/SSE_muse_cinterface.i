%module SSE_muse_interface

%include typemaps.i

%{
int initialize( double z_in, double neta_in, double bwind_in,
                double hewind_in, double sigma_in, int ifflag_in,
                int wdflag_in, int bhflag_in, int nsflag_in, double mxns_in,
                double pts1_in, double pts2_in, double pts3_in );

void evolve(int kw, double mass, double mt, double r, double lum,
            double mc, double rc, double menv, double renv,
            double ospin, double epoch, double tm, double tphys,
            double tphysf,
            int *kw1, double *mass1, double *mt1, double *r1, double *lum1,
            double *mc1, double *rc1, double *menv1, double *renv1,
            double *ospin1, double *epoch1, double *tm1, double *tphys1,
            double *tphysf1);
double get_time_step(int kw, double mass, double age, double mt, double tm, double epoch);
%}

int initialize( double z_in, double neta_in, double bwind_in,
                double hewind_in, double sigma_in, int ifflag_in,
                int wdflag_in, int bhflag_in, int nsflag_in, double mxns_in,
                double pts1_in, double pts2_in, double pts3_in );

void evolve(int kw, double mass, double mt, double r, double lum,
            double mc, double rc, double menv, double renv,
            double ospin, double epoch, double tm, double tphys,
            double tphysf,
            int *OUTPUT, double *OUTPUT, double *OUTPUT, double *OUTPUT,
            double *OUTPUT, double *OUTPUT, double *OUTPUT, double *OUTPUT,
            double *OUTPUT, double *OUTPUT, double *OUTPUT, double *OUTPUT,
            double *OUTPUT, double *OUTPUT);

double get_time_step(int kw, double mass, double age, double mt, double tm, double epoch);

