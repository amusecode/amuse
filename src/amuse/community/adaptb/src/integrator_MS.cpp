/////////////////////////////////////////////////////////////////////////////
// integrator_MS is an N-Body code with the following features: 
//
// Integrator	: Bulirsch-Stoer w. leapfrog
// Time step	: Shared, adaptive w. dt ~ { sqrt(r_ij/da_ij) }_min
// Data type	: mpreal
// Cores	: Single
//
// Compile	: see makefile
// Run		: ./integrator_MS.exe file_in file_out file_log t_sim dt_print dt_max dt_factor epsilon numBits softening t_lim
//
// Tjarda Boekholt
// 17-12-2011
/////////////////////////////////////////////////////////////////////////////

#include <iostream>
using namespace std;

#include "Bs_integrator.h"
#include "Clock.h"

int main(int argc, char* argv[])
{
  ////////////////////////////////////////////////////////
  // Set precision
  ////////////////////////////////////////////////////////
  int numBits = atoi(argv[9]);  
  int Lw = numBits/4;
  mpreal::set_default_prec(numBits);  
  cout.precision(16);
  ofstream odata;
  odata.precision(Lw);

  ////////////////////////////////////////////////////////
  // Declare user variables
  ////////////////////////////////////////////////////////
  string file_in, file_out, file_log;
  mpreal t_sim, dt_print, dt_max, dt_factor;  
  mpreal epsilon;
  mpreal softening;
  mpreal t_lim;

  int sim_state = 0;

  ////////////////////////////////////////////////////////
  // Set user variables
  ////////////////////////////////////////////////////////
  file_in = argv[1];
  file_out = argv[2];
  file_log = argv[3];

  t_sim = argv[4];
  dt_print = argv[5];
  dt_max = argv[6];
  dt_factor = argv[7];

  epsilon = argv[8];

  softening = argv[10];

  t_lim = argv[11];

  ////////////////////////////////////////////////////////
  // Declare program variables
  ////////////////////////////////////////////////////////
  int n_max, k_max;
  mpreal E0, E, dE;

  ////////////////////////////////////////////////////////
  // Set program variables
  ////////////////////////////////////////////////////////
  n_max = 64;
  k_max = 64;

  ////////////////////////////////////////////////////////
  // Configure Objects
  ////////////////////////////////////////////////////////
  Force Force(softening);
  Cluster cluster(file_in, Force);
  Bs_integrator bs(epsilon, n_max, k_max);
  Clock clock(cluster.get_t(), t_sim, dt_print, dt_max, dt_factor);

  ////////////////////////////////////////////////////////
  // Run Simulation
  ////////////////////////////////////////////////////////
  clock.Start_timer();

  odata.open( file_out.c_str() );
  if( !odata ) {
    cerr << "Could not open " << file_out << "!" << endl;
    sim_state = 1;
    exit(1);
  }
  else {
    cluster.print(odata);

    cerr << endl;
    cerr << clock.get_progress() << "%" << endl;

    // Initial calculations
    E0 = cluster.get_E();

    while( !clock.alarm() )
    {
      cluster.calc_a_dt();
      clock.calc_dt( cluster.get_dt() );

      mpreal dt = clock.get_dt();
      bs.integrate(cluster, dt);
      clock.set_dt(dt);	

      if( !bs.converged() ) {
        cerr << "No Convergence Reached, Simulation Aborted!" << endl;
        sim_state = 2;
        clock.abort();
      }
      else {
        clock.tick();
        cluster.set_t( clock.get_t() );

        if( clock.to_print() ) {
	  Cluster cl_exp = cluster;

          cl_exp.calc_a();
	  mpreal dt_exp = clock.get_t_print() - clock.get_t();

	  bs.integrate(cl_exp, dt_exp);

	  cl_exp.set_t( clock.get_t_print() );
          cl_exp.print(odata);
        }

	if( clock.read() > t_lim) {
          sim_state = 3;
          clock.abort();	  
	}

        cerr << clock.get_progress() << "%" << endl;
      }

    }

  }
  odata.close();  

  t_sim = cluster.get_t();
  E = cluster.get_E();
  dE = log10(abs((E-E0)/E0));
  clock.stop_timer();

  odata.open( file_log.c_str() );
  if( !odata )
  {
    cerr << "Could not open " << file_log << "!" << endl;
    exit(1);
  }
  else
  {
    odata << "sim_state = " << sim_state << endl;

    odata << "N         = " << cluster.get_N() << endl;

    odata << "t_sim     = " << t_sim << endl;
    odata << "dt_print  = " << dt_print << endl;
    odata << "dt_max    = " << dt_max << endl;
    odata << "dt_factor = " << dt_factor << endl;

    odata << "epsilon   = " << epsilon << endl;
    odata << "numBits   = " << numBits << endl;

    odata << "softening = " << softening << endl;

    odata << "t_cpu     = " << clock.get_timer() << endl;
    odata << "dE        = " << dE << endl;
  }  

  ////////////////////////////////////////////////////////
  // Finish Program
  ////////////////////////////////////////////////////////
  cout << "t_cpu dE: " << clock.get_timer() << " " << dE << endl;

  return 0;
}


