
#include "hacs64.h"

#ifndef __APPLE__
#define __LINUX__
#endif

#ifdef __APPLE__
#include <Accelerate/Accelerate.h>
#include <xmmintrin.h>
inline void fpe_catch() {
	_mm_setcsr( _MM_MASK_MASK &~
			(_MM_MASK_OVERFLOW|_MM_MASK_INVALID|_MM_MASK_DIV_ZERO) );
}
#elif defined __LINUX__
#include <fenv.h>
void fpe_catch(void) {
	/* Enable some exceptions. At startup all exceptions are masked. */
	feenableexcept(FE_INVALID|FE_DIVBYZERO|FE_OVERFLOW);
}
#else
void fpe_catch(void) {}
#endif

std::vector<particle> node::ptcl;
std::vector<node>     node::node_heap;
std::vector<std::pair<node*, node*> > node::pair_list;


void read_dumbp(
		FILE *fin,
		std::vector<dvec3>  &pos,
		std::vector<dvec3>  &vel,
		std::vector<double> &mass)
{
	const int linesz=256;
	char line[linesz];

	pos.clear();
	vel.clear();
	mass.clear();

	pos.reserve(128);
	vel.reserve(128);
	mass.reserve(128);

	fgets(line, linesz, fin);
	while (line[0] == ';') 
		fgets(line, linesz,fin);

	while(!feof(fin)) {
		dvec3 ipos, ivel;
		double imass;
		int dummy;
		sscanf(line, "%d %lf %lf %lf %lf %lf %lf %lf",
				&dummy, 
				&imass,
				&ipos.x, &ipos.y, &ipos.z,
				&ivel.x, &ivel.y, &ivel.z);
		pos. push_back(ipos);
		vel. push_back(ivel);
		mass.push_back(imass);
		fgets(line, linesz, fin);
	}
	const int nbodies = pos.size() - 1;
	pos.resize(nbodies);
	vel.resize(nbodies);
	mass.resize(nbodies);
	fprintf(stderr, " Read %d bodies \n", nbodies);
}


int main(int argc, char * argv[]) 
{
	fpe_catch();

	std::vector<dvec3 > pos, vel;
	std::vector<double> mass;

	read_dumbp(stdin, pos, vel, mass);

	hacs64::Nbody nbody;

	const int nbodies = pos.size();

	double eps     = 0.0/nbodies;
  double eps2    = eps*eps;
	double hmax    = std::sqrt(0.5);
	double eta_irr = 0.6;
	double eta_reg = 0.1;
#if 1
  eta_irr = 0.8;
  eta_reg = 0.14;
//  eta_reg = 0.04;
#endif
	double dt_max  = 1.0/(1 << 4);


	fprintf(stderr, " Initializing system ... \n");
	nbody.initialize(
			nbodies,
			&mass[0],
			&pos [0],
			&vel [0],
			eps2,
			hmax,
			eta_irr,
			eta_reg,
			dt_max);
	fprintf(stderr, " ... done ... \n");

	nbody.potential();

	hacs64::Energy E0(nbody.ptcl, nbody.eps2);
  double e_prev = E0.e;
	E0.print(stderr, "Initial E : ");
	E0.print_mom(stderr, "Initial Momentum : ");

	double dt_log = dt_max;
	//	double dt_out = 1.0;
	double  t_log = 0.0;
	//	double  t_out = 0.0;
	double de_max = 0.0;

  double  t_out = 0.0;
  double dt_out = 100.0;
	double  t_end = 1.0;

	nbody.scheduler.debug_dump();

	int niter_max = 1000000000;
  const double t_beg = get_wtime();
	while (nbody.t_global < t_end && nbody.iteration < niter_max) 
  {

		const double tcur = nbody.t_global;
#if 0
		{
			Timer time("  integrate for dt_max : ", stderr);
			while (nbody.t_global - tcur != dt_max) {
				nbody.iterate();
				if (nbody.iteration % 100 == 0) 
					nbody.scheduler.debug_dump();
				fprintf(stderr, " t= %g, dt= %g, eta=%g iteration= %d , nsteps= %lld : nsteps_irr= %lld  nsteps_reg= %lld  \n",
						nbody.t_global,
						nbody.dt_global,
						nbody.eta_reg,
						nbody.iteration, 
						nbody.nsteps,
						nbody.nsteps_irr,
						nbody.nsteps_reg);
			}
		}
#endif
    const double t0_dtmax = get_wtime();
		{
      if (nbody.t_global >= t_out)
      {
        char fn[256];
        sprintf(fn, "%s/iter_%.5d.dumbp", "out", int(t_out/dt_out));
        fprintf(stderr, " ------ writing snapshout @ t=  %g to  %s ------- \n", nbody.t_global, fn);
        nbody.write_dumbp(fn);
        t_out += dt_out;
      }
      while (nbody.t_global - tcur != dt_max) {
        nbody.iterate();
//        if (nbody.iteration % 100 == 0) 
      if (0)
        {
          fprintf(stderr, " t= %g, dt= %g, eta=(%g, %g) iteration= %d , nsteps= %lld : nsteps_irr= %lld  nsteps_reg= %lld  Tcpu= %g h\n",
              nbody.t_global,
              nbody.dt_global,
              nbody.eta_irr,
              nbody.eta_reg,
              nbody.iteration, 
              nbody.nsteps,
              nbody.nsteps_irr,
              nbody.nsteps_reg,
              (get_wtime() - t_beg)/3600.0);
        }
     //   if (nbody.iteration % 1000 == 0) 
    //      nbody.scheduler.debug_dump();
      }

    }
    fprintf(stderr," >>> integrate for dt_max: %g  sec  :: all= %g  reg= %g irr= %g  ngb= %g  corr= [%g %g %g := %g] diff= %g\n",
        get_wtime() - t0_dtmax, 
        nbody.dt_all,
        nbody.dt_reg,
        nbody.dt_irr,
        nbody.dt_ngb,
        nbody.dt_corr_reg,
        nbody.dt_corr_irr,
        nbody.dt_corr_mix,
        nbody.dt_corr_reg+nbody.dt_corr_irr+nbody.dt_corr_mix,
        nbody.dt_all - (nbody.dt_reg + nbody.dt_irr + nbody.dt_ngb + nbody.dt_corr_reg + nbody.dt_corr_irr + nbody.dt_corr_mix));
    nbody.dt_all = nbody.dt_reg = nbody.dt_irr = nbody.dt_ngb = 0;
    nbody.dt_corr_reg = nbody.dt_corr_irr = nbody.dt_corr_mix = 0;

    if (nbody.t_global == t_log) 
    {
      nbody.potential();
      E0 = hacs64::Energy(nbody.ptcl, nbody.eps2);
    }

    if (nbody.t_global >= t_log) 
    {
      nbody.potential();
      const hacs64::Energy E1(nbody.ptcl, nbody.eps2);
      de_max = std::max(de_max, std::abs((E1.e - E0.e)/E0.e));

      fprintf(stderr,  " t= %g  dt= %g  Tcpu= %g h: ninter= %lld %lld [ %g ] %lld  ",
          nbody.t_global, nbody.dt_global,
          (get_wtime() - t_beg)/3600.0,
          nbody.ninter_irr,
          nbody.ninter_reg,
          nbody.ninter_irr*1.0/nbody.ninter_reg,
          nbody.ninter_irr + nbody.ninter_reg);
      fprintf(stderr, " d(de) = %g  ", E1.e - e_prev);
      (E1 - E0).print(stderr, "Delta E : ");
      e_prev = E1.e;
      (E1 - E0).print_mom(stderr, "Momentum  : ");
      fprintf(stderr,  " >>> de_max= %g <<< \n", de_max);
      t_log += dt_log;
    }

  }

  std::cerr << "end-of-program" << std::endl;
}
