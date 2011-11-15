/*************************************************************************
description:
This file contains the main routine to do radiative transfer calculations
with SimpleX.

Copyright Jan-Pieter Paardekooper and Chael Kruip October 2011

This file is part of SimpleX.

SimpleX is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

SimpleX is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with SimpleX.  If not, see <http://www.gnu.org/licenses/>.
**************************************************************************/

#include "Main.h"

using namespace std;

int main(int argc, char** argv) {  

  cerr.precision(5);

  MPI::Init( argc, argv );

  unsigned int COMM_RANK = MPI::COMM_WORLD.Get_rank();  /* Get my rank          */    
  //int COMM_SIZE = MPI::COMM_WORLD.Get_size();  /* Get the total number of processors */

  double t0 = MPI::Wtime();

  if (argc != 2) {
    
    if(COMM_RANK == 0){
      cerr << "Use:" << endl
           << "mpirun < mpi param > " << argv[0] << " < input file >" << endl;
    }
    
    MPI::Finalize();
    return EXIT_SUCCESS;
  }

  //initialise the simulation
  SimpleX SimpleXGrid;

  if(COMM_RANK==0){

    cerr << endl << endl
	 << "**********************************************************************" << endl
	 << "*                            SimpleX 2.0                             *" << endl
	 << "**********************************************************************" << endl
	 << "*                                                                    *" << endl
	 << "*                       Jan-Pieter Paardekooper                      *" << endl
	 << "*                     jppaarde@strw.leidenuniv.nl                    *" << endl
	 << "*                                                                    *" << endl    
	 << "*                           Chael Kruip                              *" << endl
	 << "*                       kruip@strw.leidenuniv.nl                     *" << endl
	 << "*                                                                    *" << endl
	 << "*                        Version: May 2010                           *" << endl
	 << "**********************************************************************" << endl << endl;
    
  }

  //===============   Triangulation Section  ==========================//

  if(COMM_RANK == 0){
    cerr << " (" << COMM_RANK << ") Computing triangulation" << endl;
  }

  unsigned int start_run = 0;
  SimpleXGrid.init_triangulation(argv[1]);
  
  unsigned int numRuns = SimpleXGrid.get_numRuns();
  
  /****** MAIN LOOP ******/
  for( unsigned int run=start_run; run<numRuns; run++ ){

    if(COMM_RANK == 0){

      cerr << endl << endl
        << "**********************************************************************" << endl
        << "                                Run: " << run+1 << "                              " << endl
        << "**********************************************************************" << endl << endl;

      SimpleXGrid.simpleXlog << endl << endl
        << "**********************************************************************" << endl
        << "                                Run: " << run+1 << "                              " << endl
        << "**********************************************************************" << endl << endl;


    }

    //TESTING
    if(run>0)
      SimpleXGrid.reinit_triangulation();

    // For the first run, compute the properties of the triangulation
    if(run == 0){

      if(COMM_RANK == 0){
        cerr << " (" << COMM_RANK << ") Computing triangulation properties" << endl;
      }

      SimpleXGrid.compute_site_properties();

      //================  Physics Section     =============================//

      MPI::COMM_WORLD.Barrier();

      if(COMM_RANK == 0)
        cerr << " (" << COMM_RANK << ") Computing physical properties" << endl;

      SimpleXGrid.compute_physics( run );
      SimpleXGrid.remove_border_simplices();

    }

    SimpleXGrid.parameter_output( run );

    //================  Radiative Transport Section     =============================//

    if(COMM_RANK == 0){
      cerr << " (" << COMM_RANK << ") Entering Radiative Transport Routine" << endl;
      SimpleXGrid.simpleXlog << endl << "  Entering Radiative Transport Routine" << endl;
    }

    int stop = 0;
    SimpleXGrid.radiation_transport( run );
    
    if(stop){
      break;
    }

    //================  Output Section     =============================//

    SimpleXGrid.generate_output( run );

    if(COMM_RANK==0){
      cerr << endl << endl << endl << " END OF RUN " << endl << endl;
    }

  }

  //SimpleXGrid.clear_temporary();

  double t1 = MPI::Wtime();
  if(COMM_RANK==0){
    cerr << endl << " *****  Total Simulation time: " << t1-t0 << " seconds *****" << endl;
    SimpleXGrid.simpleXlog << endl << " *****  Total Simulation time: " << t1-t0 << " seconds *****" << endl;
  }



  MPI::Finalize();
  return EXIT_SUCCESS;

}
