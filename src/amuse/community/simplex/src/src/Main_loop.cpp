/*************************************************************************
file:         Main.cpp
author:       Jan-Pieter Paardekooper
mail:         jppaarde@strw.leidenuniv.nl
version:      0.1
last change:  04.04.2008
---------------------------------------------------------------------
description:
This file contains the main routine to do radiative transfer calculations
with SimpleX.
**************************************************************************/
/*
 * Date: Name
 * Put additional comments here
 *
 * 03.04.08 Jan-Pieter Paardekooper
 * Put comments suitable for doxygen in the relevant places
 *
*/

/***** To Do *******
 *
 *******************/

#include "Main_loop.h"

using namespace std;

//extern int main_loop(int argc, char** argv, SimpleX *SimpleXGrid, unsigned int COMM_RANK);

int main_loop(int argc, char** argv, SimpleX *SimpleXGrid) {  

  unsigned int COMM_RANK = MPI::COMM_WORLD.Get_rank();  /* Get my rank          */    
  //int COMM_SIZE = MPI::COMM_WORLD.Get_size();  /* Get the total number of processors */


  double t0 = MPI::Wtime();

  if(COMM_RANK==0){
    cerr << endl << endl
	 << "**********************************************************************" << endl
	 << "*                            SimpleX 2.2                             *" << endl
	 << "**********************************************************************" << endl
	 << "*                                                                    *" << endl
	 << "*                       Jan-Pieter Paardekooper                      *" << endl
	 << "*                     jppaarde@strw.leidenuniv.nl                    *" << endl
	 << "*                                                                    *" << endl    
	 << "*                           Chael Kruip                              *" << endl
	 << "*                       kruip@strw.leidenuniv.nl                     *" << endl
	 << "*                                                                    *" << endl    
	 << "*                         Jelle Ritzerveld                           *" << endl
	 << "*                        ritzervj@gmail.com                          *" << endl
	 << "*                                                                    *" << endl
	 << "*                        Version: Dec 2008                           *" << endl
	 << "**********************************************************************" << endl << endl;
  }


  //===============   Triangulation Section  ==========================//

  if(COMM_RANK == 0){
    cerr << " (" << COMM_RANK << ") Computing triangulation" << endl;
  }

  (*SimpleXGrid).init_triangulation(argv[1]);

  unsigned int numRuns = (*SimpleXGrid).get_numRuns();
  bool dynamic_updates = (*SimpleXGrid).get_updates();
  bool vertex_movement = (*SimpleXGrid).get_vertex_movement();
  bool virtual_movement = (*SimpleXGrid).get_virtual_movement();
  bool dynamic_grid = 0;
  if( dynamic_updates || vertex_movement ){
    dynamic_grid=1;
  }

  /****** MAIN LOOP ******/
  for( unsigned int run=0; run<numRuns; run++ ){

    if(COMM_RANK == 0){

      cerr << endl << endl
	   << "**********************************************************************" << endl
	   << "                                Run: " << run+1 << "                              " << endl
	   << "**********************************************************************" << endl << endl;

      (*SimpleXGrid).simpleXlog << endl << endl
		 << "**********************************************************************" << endl
		 << "                                Run: " << run+1 << "                              " << endl
		 << "**********************************************************************" << endl << endl;


    }

    //bool to determine whether triangulation has changed since previous run
    bool calcGridProps = 1;
    
    //for first run always true, so check 
    if(run > 0){
      //dynamic grid: vertex movement and/or site deletion
      if(dynamic_grid){

	if(COMM_RANK == 0){
	  cerr << " (" << COMM_RANK << ") Computing grid updates" << endl;
	}

        calcGridProps = (*SimpleXGrid).grid_dynamics();

      }else{

	if( virtual_movement ){
	  (*SimpleXGrid).move_vertices_virtually();
	}

	calcGridProps=0;

      }
    }

    if(calcGridProps){

      if(COMM_RANK == 0){
	cerr << " (" << COMM_RANK << ") Computing triangulation properties" << endl;
      }

      (*SimpleXGrid).compute_site_properties();

      //================  Physics Section     =============================//

      MPI::COMM_WORLD.Barrier();

      if(COMM_RANK == 0)
	cerr << " (" << COMM_RANK << ") Computing physical properties" << endl;

      (*SimpleXGrid).compute_physics( run );

      (*SimpleXGrid).remove_border_simplices();

    }

    (*SimpleXGrid).parameter_output( run );

    //================  Radiative Transport Section     =============================//

    if(COMM_RANK == 0){
      cerr << " (" << COMM_RANK << ") Entering Radiative Transport Routine" << endl;
      (*SimpleXGrid).simpleXlog << endl << "  Entering Radiative Transport Routine" << endl;
    }

    (*SimpleXGrid).radiation_transport( run );

//     MPI::COMM_WORLD.Barrier();
//     //test, wait 5 secdonds between outputs
//     double t_begin = MPI::Wtime();
//     bool stop = 0;
//     while(!stop){
//       double t_end = MPI::Wtime();
//       if( t_end - t_begin > 50){
// 	stop = 1;
//       }      
//     }


    //================  Output Section     =============================//

    (*SimpleXGrid).generate_output( run );


  }

  //(*SimpleXGrid).clear_temporary();

  double t1 = MPI::Wtime();
  if(COMM_RANK==0){
    cerr << endl << " *****  Total Simulation time: " << t1-t0 << " seconds *****" << endl;
    (*SimpleXGrid).simpleXlog << endl << " *****  Total Simulation time: " << t1-t0 << " seconds *****" << endl;
  }

  return 0;

}
