/*************************************************************************
file:         SimpleX.cpp
author:       Jan-Pieter Paardekooper
mail:         jppaarde@strw.leidenuniv.nl
version:      0.1
last change:  15.01.2009
---------------------------------------------------------------------
description:
This file contains the implementation of the SimpleX class that does
the radiative transfer calculations
**************************************************************************/
/*
* Date: Name
* Put additional comments here
*
* 22.03.2009: JPP
* Added opportunity to preserve directions during grid updates (although memory
* consuming). Grid updates are now possible. Fixed bug in redistribute_photons(), 
* intensity is now thrown away when cell is flash ionised.
*
* 18.03.2009: JPP
* Changed the SimpleX algorithm so that directions are preserved. A tesselation
* of the unit sphere is superimposed on every site in which the intensities are
* stored. Currently the tesselation of the sphere is given by HEALPix.  
* 
*/

/***** Bugs ********
*
* In parallel hdf5 writing, chunk_size has to be as big as total number of simplices
*
* 
*******************/

/***** To Do *******
*
* In send_sites_marked_for_removal() sites need only be send to master,
* instead of to all procs. Same for send_sites_to_remove().
*
* In update routine
*  * If more than half of points flagged need to be removed, turn around criterium (faster)
*  * Let number of bins be a function of number of points to be deleted
*  
*
*******************/


#include "SimpleX.h"

#include "Map16.h"
#include "Map21.h"
#include "Map32.h"
#include "Map42.h"
#include "Map64.h"
#include "Map84.h"

#include <algorithm>

using namespace std;

// constructor
SimpleX::SimpleX(){

  //rank of this processor and total number of processors
  COMM_RANK = MPI::COMM_WORLD.Get_rank();    
  COMM_SIZE = MPI::COMM_WORLD.Get_size();    

  //cross section and fraction of photons that is capable of ionising
  cross_H = 0.0;
  cross_dust = 0.0;

  //conversion factors between numerical and physical units
  UNIT_T = 0.0;
  UNIT_L = 0.0;

  //total volume, total number of atoms and maximal resolution of/in simulation domain
  totalVolume = 0.0;
  totalAtoms = 0.0;
  maxRes = 0;

  time_conversion = 1.0;

  dust_to_gas_ratio = 0.0;
  cross_dust = 0.0;

  //euler angles
  euler_phi = 0.0;
  euler_theta = 0.0;
  euler_psi = 0.0;

  //number of pixels of HEAL_PIX sphere
  numPixels = 0;

  //orientation of the unit sphere tesselation
  orientation_index = 0;

  //only 3D at the moment, this gotta change at one point!
  dimension = 3;

  // Random number generator
  ran = gsl_rng_alloc (gsl_rng_taus);

  //output to log file
  simpleXlog.open( "SimpleX.log" );

}

SimpleX::~SimpleX(){

  // free the random number generator 
  gsl_rng_free ( ran );

  //close output files
  simpleXlog.close();
  IFront_output.close();

  //clear data structures
  vertices.clear();
  clear_temporary();
  sites.clear();
  simplices.clear();


}



/****************************************************************************************/
/*                       Grid Initialisation Functions                                  */
/****************************************************************************************/


/****  Initialize the simulation and compute the triangulation  ****/
void SimpleX::init_triangulation(char* inputName){

  //keep track of cpu time
  double t0 = MPI::Wtime();

  //read in the user specified parameter file
  read_parameters( inputName );

  //keep track of the number of points the simulation started with
  origNumSites = numSites;

  //this might have to change when the number of sites is possible
  //to get bigger than the original number during the simulation
  vertex_id_max = numSites+borderSites;

  //calculate the total number of sweeps from the total simulation time
  //and the time step
  unsigned int total_numSweeps = simTime/UNIT_T;

  if( updates ){
    //number of runs is the minumum of the output
    //and the grid dynamics
    numRuns = max( numGridDyn, numOutputs);
  }else{
    numRuns = numOutputs;
  }

  //convert time for nice output
  time_conversion = 1.0;
  while( time_conversion*simTime/numOutputs < 1.0 ){
    time_conversion *= 1e3;
  }

  //calculate the number of sweeps per run
  numSweeps = total_numSweeps/numRuns;

  // initialise random number generator with random seed
  gsl_rng_set ( ran , randomSeed );

  //set the correct directions headers
  set_direction_bins();

  //read the input file with vertex positions 
  if(fillChoice != AUTOMATIC){
    read_vertex_list();
  }

  //only the master proc creates the point distribution (in case of automatic filling) and
  //the boundary around the domain (always)
  if(COMM_RANK == 0){

    //create homogeneous Poisson distribution of points
    if(fillChoice == AUTOMATIC){
      poisson_square();
    }

    //set boundary around unity domain
    if(periodic){
      create_periodic_boundary();
    }else{
      create_boundary();
    }

    //create octree on master proc
    create_vertex_tree();

    //decompose the domain
    decompose_domain();

    //assign process to vertices
    assign_process();

  }

  if( COMM_SIZE > 1 ){
    //send the domain decomposition to other procs
    send_dom_dec();
    //send the vertex positions round to all procs, so every proc (temporary!) has its own copy
    send_vertices();
  }

  //now that all procs have a list of vertices, create octree
  if( COMM_RANK != 0 ){
    create_vertex_tree();
  }

  MPI::COMM_WORLD.Barrier();

  // triangulate the subboxes
  compute_triangulation();

  //create a site at each vertex from the list of simplices that was obtained 
  //from the triangulation functions, containing the physical parameters
  create_sites();

  if(periodic){
    //set_periodic_sites_to_send();
    remove_periodic_sites();
  }
  
  //the list of sites was obtained from list of simplices, and therefore 
  //have an order which might lead to problems when using the dynamic update routines
  //shuffle_sites();
  
  assign_site_ids();  

  //determine which sites are ballistic and which direction conserving
  initiate_ballistic_sites();

  //write the time it took to triangulate the points to log file
  double t1 = MPI::Wtime();
  if( COMM_RANK == 0 ){
    simpleXlog << endl << "  Calculating triangulation took " << t1-t0 << " seconds" << endl << endl;
  }

}


/**** read in parameter file  ****/
void SimpleX::read_parameters( char* inputName ){

  //create the keyvalue variable
  string strFilename = inputName;
  ConfigFile KeyValueFile(strFilename);

  //strings needed for reading in strings
  string fill, bbChoice, recChoice, diffChoice, upd, 
    mov, movVirt, keep_dir, cycl, hp, phot_cons, dust;

  //read in fill choice
  KeyValueFile.readInto(fill, "fillChoice", string("AUTOMATIC"));
  if(fill == "AUTOMATIC"){
    fillChoice = AUTOMATIC;
  }else if(fill == "READ"){
    fillChoice = READ;
  }else{
    fillChoice = AUTOMATIC;
  }
  //read in input files
  KeyValueFile.readInto(inputDensityFile, "inputName", string("./points.pnt") );
  //random seed to be used
  KeyValueFile.readInto(randomSeed, "randomSeed");
  //dimension
  KeyValueFile.readInto(dimension, "dimension", (short) 3 );
  //Number of points in case of automatic filling
  KeyValueFile.readInto(numSites, "numSites", (unsigned long long int) 262144 );
  //Periodic boundary conditions?
  KeyValueFile.readInto(periodic, "periodic", (bool) 0 );
  //Number of points in boundary
  KeyValueFile.readInto(borderSites, "borderPoints", (unsigned) 25000 );
  //Buffer around unity domain in which extra points are placed
  KeyValueFile.readInto(borderBox, "padding", (float) 0.1 );
  //hilbert order to determine domain decomposition
  KeyValueFile.readInto(hilbert_order, "hilbert_order", (unsigned int) 2 );
  //Buffer around subbox in which extra points are triangulated
  KeyValueFile.readInto(padding_subbox, "padding_subbox", (float) 0.25 );

  //Number of sweeps per run
  KeyValueFile.readInto( UNIT_T, "time_step", (double) 0.05 );
  //Number density in homogeneous case
  KeyValueFile.readInto(homDens, "homDens", (float) 0.01 );
  //Size of simulation domain
  KeyValueFile.readInto(sizeBox, "sizeBox", (double) 13200.0 );
  //Simulation time
  KeyValueFile.readInto(simTime, "simTime", (double) 500.0 );
  //Black body spectrum for source?
  KeyValueFile.readInto( blackBody, "bbSpectrum", (bool) 0);
  //Effective temperature for source
  KeyValueFile.readInto(sourceTeff, "sourceTeff", (double) 0.0 );
  //Temperature of ionized gas
  KeyValueFile.readInto(gasIonTemp, "gasIonTemp", (double) 1.0e4 );

  //Include dust?
  KeyValueFile.readInto(dust, "dust_model", string("NO_DUST") );
  if(dust == "SMC"){
    dust_model = SMC;
  }else if(dust == "LMC"){
    dust_model = LMC;
  }else if(dust == "MW"){
    dust_model = MW;
  }else{
    dust_model = NO_DUST;
  }

  //metallicity
  KeyValueFile.readInto(metallicity, "metallicity", (double) 0.2 );

  KeyValueFile.readInto(dust_sublimation, "dust_sublimation", (bool) 0 );

  //Source strength (number of ionising photons)
  KeyValueFile.readInto(source_strength, "sourceStrength", (float) 5.0 );
  //Units of source
  KeyValueFile.readInto(UNIT_I, "sourceUnits", (double) 1.0e48 );
  //source radius
  KeyValueFile.readInto(source_radius, "sourceRadius", (float) 0.02 );
  //number of points in source
  KeyValueFile.readInto(source_points, "sourcePoints", (unsigned) 25 );
  //x-position
  KeyValueFile.readInto(source_x, "sourceX", (float) 0.5 );
  //y-position
  KeyValueFile.readInto(source_y, "sourceY", (float) 0.5 );
  //z-position
  KeyValueFile.readInto(source_z, "sourceZ", (float) 0.5 );

  //source resides in cell or is entire cell
  KeyValueFile.readInto(source_inside_cell, "source_inside_cell", (bool) 0);

  //Include recombination?
  KeyValueFile.readInto( recombination, "recombination", (bool) 1 );

  //Only ballistic transport?
  KeyValueFile.readInto( ballisticTransport, "onlyBallistic", (bool) 0 );
  //Only direction_conserving transport?
  KeyValueFile.readInto( dirConsTransport, "onlyDirCons", (bool) 0 );
  //Combined transport?
  KeyValueFile.readInto( combinedTransport, "combinedTransport", (bool) 1 );

  //Number of directions
  KeyValueFile.readInto( number_of_directions, "number_of_directions", (unsigned) 42 );
  //Number of orientations
  KeyValueFile.readInto( number_of_orientations, "number_of_orientations", (unsigned) 100 );

  //Number of outputs
  KeyValueFile.readInto(numOutputs, "outputs", (unsigned) 50 );

  //output IFront for source in centre?
  KeyValueFile.readInto(give_IFront, "IFront", (bool) 0 );
  //calculate escape fraction?
  KeyValueFile.readInto(give_escape_fraction, "escape_fraction", (bool) 0 );


  //Chunk size for hdf5 writing
  KeyValueFile.readInto( chunk_size, "chunkSize", (unsigned) 100000 );

  //Maximum number of messages to send in MPI routines
  KeyValueFile.readInto( max_msg_to_send, "max_msg_to_send", (unsigned) 100000 );

  //number of reference pixels in the case HEAL_PIX is used
  KeyValueFile.readInto( num_ref_pix_HP, "number_of_reference_pixels_HealPix", 5);

  //Use dynamical updates?
  KeyValueFile.readInto( updates, "remove_vertices", (bool) 0 );
  KeyValueFile.readInto( numGridDyn, "number_of_grid_dyn", (unsigned) 100 );

  // Maximal angle between straightforward direction and Delaunay lines (second and third)
  KeyValueFile.readInto( straightAngle, "straight_angle", (float) 90.0 );
  straightAngle *= M_PI/180.0;

  KeyValueFile.readInto(straight_from_tess, "straight_from_tess", (bool) 1 );

  // Do cyclic check?
  KeyValueFile.readInto( cyclic_check, "cyclic_check", (bool) 0 );

  //Vertex deletion criterium
  KeyValueFile.readInto(switchTau, "maxRatioDiff");
  //Maximum number of neighbours that can be deleted
  KeyValueFile.readInto(minResolution, "minResolution");
  //number of bins in update routine
  KeyValueFile.readInto(nbins, "nbins", (unsigned int) 100 );

  //Use temporal photon conservation?
  KeyValueFile.readInto( photon_conservation, "temporal_photon_conservation", (bool) 1 );

  //write all information to logfile
  if(COMM_RANK == 0){
    simpleXlog << endl <<" *****  This is SimpleX version 2.4  ***** " << endl << endl;
    simpleXlog << "  Number of processors: " << COMM_SIZE << endl;
    if(fillChoice == AUTOMATIC){
      simpleXlog << "  Fill Choice: AUTOMATIC" << endl;
      simpleXlog << "  Point process : Poisson" << endl;
      simpleXlog << "  Homogeneous density: " << homDens << " cm^-1" << endl;
      simpleXlog << "  Size of the simulation domain: " << sizeBox << " pc" << endl;
      simpleXlog << "  Source strength: " << source_strength*UNIT_I << endl;
      simpleXlog << "  Number of source points: " << source_points << endl;
      simpleXlog << "  Source radius: " << source_radius << endl;
      simpleXlog << "  Source position: (" << source_x << ", " << source_y << ", " 
        << source_z << ")" << endl;
    } else{ 
      simpleXlog << "  Fill Choice: Read file" << endl;
      simpleXlog << "  Input file: " << inputDensityFile << endl;
    }
    if(periodic){
      simpleXlog << "  Boundaries are periodic" << endl;
    } 
    simpleXlog << "  Number of boundary points    : " << borderSites << endl;
    simpleXlog << "  Buffer around unity domain in which boundary points are placed: " << borderBox << endl;
    simpleXlog << "  Hilbert order                : " << hilbert_order << endl;
    simpleXlog << "  Time step                    : " << UNIT_T << endl;
    simpleXlog << "  Number of outputs            : " << numOutputs << endl;
    simpleXlog << "  Chunk size for hdf5 writing  : " << chunk_size << endl;     
    simpleXlog << "  Total simulation time        : " << simTime << " Myr" << endl;
    if(blackBody){ 
      simpleXlog << "  Source has black body spectrum of " << sourceTeff << " K" << endl 
        << "  Resulting temperature of the ionized gas is " << gasIonTemp << " K" << endl;
    }
    if(source_inside_cell){
      simpleXlog << " Source is inside the cell and ionises the gas in its own cell " << endl;
    }else{
      simpleXlog << " Source is made up of entire cell and does not ionise its own cell " << endl;
    }

    if(dust_model == SMC){
      simpleXlog << " SMC dust model";
      if(dust_sublimation){
        simpleXlog << " with total dust sublimation  " << endl;
      }else{
        simpleXlog << ", no dust sublimation " << endl;
      }
    }else if(dust_model == LMC){
      simpleXlog << " LMC dust model";
      if(dust_sublimation){
        simpleXlog << " with total dust sublimation  " << endl;
      }else{
        simpleXlog << ", no dust sublimation " << endl;
      }
    }else if(dust_model == MW){
      simpleXlog << " Milky Way dust model";
      if(dust_sublimation){
        simpleXlog << " with total dust sublimation  " << endl;
      }else{
        simpleXlog << ", no dust sublimation " << endl;
      }
    }else{
      simpleXlog << "  No dust" << endl;
    }


    if(recombination){
      simpleXlog << "  Recombination is included" << endl;
    }else{
      simpleXlog << "  Recombination is not included" << endl;
    }

    if(ballisticTransport){
      simpleXlog << endl << "  Mode of transport: ballistic transport " << endl << endl;
    }
    if(dirConsTransport){
      simpleXlog << endl << "  Mode of transport: direction conserving transport " << endl;
    }
    if(combinedTransport){
      simpleXlog << endl << "  Mode of transport: combined transport " << endl;
    }

    if(dirConsTransport || combinedTransport){
      simpleXlog << "  Number of direction bins used        : " << number_of_directions << endl;
      simpleXlog << "  Number of orientations in header file: " << number_of_orientations << endl << endl;
    }

    if(updates){
      simpleXlog << "  Vertex removal is included  " << endl;
      simpleXlog << "  Vertices removed " << numGridDyn << " times" << endl;
    }

#ifdef HEALPIX
    simpleXlog << "  Healpix sphere with " << num_ref_pix_HP 
      << " pixels is used for quick referencing in compute_solid_angles()" << endl;
#endif

    simpleXlog << "  Switch between ballistic and dirCons transport: " << switchTau << endl;
    simpleXlog << "  Minimum resolution                            : " << minResolution << endl;

    simpleXlog << "  Maximal angle between forward direction and real Delaunay direction: " 
      << straightAngle << " degrees. " << endl;

    simpleXlog << "  Check for (and repair) cyclic connections in the grid? " << cyclic_check << endl << endl;
  }


  if(give_IFront){
    //output to IFront file
    IFront_output.open( "IFront.dat" );
  }
  if(give_escape_fraction){
    //output to escape fraction file
    f_esc_output.open( "escape_fraction.dat" );
  }


}

/**** Set the  direction bins with the correct number of directions  ****/
void SimpleX::set_direction_bins(){

  //first create orient array of correct size
  orient = new float**[number_of_orientations];
  for(unsigned int i=0; i<number_of_orientations; i++){
    orient[i] = new float*[number_of_directions];
    for(unsigned int j=0; j<number_of_directions; j++){
      orient[i][j] = new float[dimension];
    }
  }

  //fill the orientations from the headers
  if(number_of_directions == 16){
    for(unsigned int i=0; i<number_of_orientations; i++){
      for(unsigned int j=0; j<number_of_directions; j++){
        for(short int k=0; k<dimension; k++){
          orient[i][j][k] = orient_16[i][j][k];
        }
      }
    }
  }else if(number_of_directions == 21){
    for(unsigned int i=0; i<number_of_orientations; i++){
      for(unsigned int j=0; j<number_of_directions; j++){
        for(short int k=0; k<dimension; k++){
          orient[i][j][k] = orient_21[i][j][k];
        }
      }
    }
  }else if(number_of_directions == 32){
    for(unsigned int i=0; i<number_of_orientations; i++){
      for(unsigned int j=0; j<number_of_directions; j++){
        for(short int k=0; k<dimension; k++){
          orient[i][j][k] = orient_32[i][j][k];
        }
      }
    }
  }else if(number_of_directions == 42){
    for(unsigned int i=0; i<number_of_orientations; i++){
      for(unsigned int j=0; j<number_of_directions; j++){
        for(short int k=0; k<dimension; k++){
          orient[i][j][k] = orient_42[i][j][k];
        }
      }
    }
  }else if(number_of_directions == 64){
    for(unsigned int i=0; i<number_of_orientations; i++){
      for(unsigned int j=0; j<number_of_directions; j++){
        for(short int k=0; k<dimension; k++){
          orient[i][j][k] = orient_64[i][j][k];
        }
      }
    }
  }else if(number_of_directions == 84){
    for(unsigned int i=0; i<number_of_orientations; i++){
      for(unsigned int j=0; j<number_of_directions; j++){
        for(short int k=0; k<dimension; k++){
          orient[i][j][k] = orient_84[i][j][k];
        }
      }
    }
  }else{
    cerr << " (" << COMM_RANK << ") Incorrect header, orientations not found " << endl;
    MPI::COMM_WORLD.Abort(-1);
  }


  //Now create mapping of correct size
  maps = new unsigned int**[number_of_orientations];
  for(unsigned int i=0; i<number_of_orientations; i++){
    maps[i] = new unsigned int*[number_of_orientations];
    for(unsigned int j=0; j<number_of_orientations; j++){
      maps[i][j] = new unsigned int[number_of_directions];
    }
  }

  //fill the mapping
  if(number_of_directions == 16){
    for(unsigned int i=0; i<number_of_orientations; i++){
      for(unsigned int j=0; j<number_of_orientations; j++){
        for(unsigned int k=0; k<number_of_directions; k++){
          maps[i][j][k] = maps_16[i][j][k];
        }
      }
    }
  }else if(number_of_directions == 21){
    for(unsigned int i=0; i<number_of_orientations; i++){
      for(unsigned int j=0; j<number_of_orientations; j++){
        for(unsigned int k=0; k<number_of_directions; k++){
          maps[i][j][k] = maps_21[i][j][k];
        }
      }
    }
  }else if(number_of_directions == 32){
    for(unsigned int i=0; i<number_of_orientations; i++){
      for(unsigned int j=0; j<number_of_orientations; j++){
        for(unsigned int k=0; k<number_of_directions; k++){
          maps[i][j][k] = maps_32[i][j][k];
        }
      }
    }
  }else if(number_of_directions == 42){
    for(unsigned int i=0; i<number_of_orientations; i++){
      for(unsigned int j=0; j<number_of_orientations; j++){
        for(unsigned int k=0; k<number_of_directions; k++){
          maps[i][j][k] = maps_42[i][j][k];
        }
      }
    }
  }else if(number_of_directions == 64){
    for(unsigned int i=0; i<number_of_orientations; i++){
      for(unsigned int j=0; j<number_of_orientations; j++){
        for(unsigned int k=0; k<number_of_directions; k++){
          maps[i][j][k] = maps_64[i][j][k];
        }
      }
    }
  }else if(number_of_directions == 84){
    for(unsigned int i=0; i<number_of_orientations; i++){
      for(unsigned int j=0; j<number_of_orientations; j++){
        for(unsigned int k=0; k<number_of_directions; k++){
          maps[i][j][k] = maps_84[i][j][k];
        }
      }
    }
  }else{
    cerr << " (" << COMM_RANK << ") Incorrect header, mappings not found " << endl;
    MPI::COMM_WORLD.Abort(-1);
  }

  // double total = sizeof(orient_16) + sizeof(maps_16) +
  //   sizeof(orient_21) + sizeof(maps_21) +
  //   sizeof(orient_32) + sizeof(maps_32) +
  //   sizeof(orient_42) + sizeof(maps_42) +
  //   sizeof(orient_64) + sizeof(maps_64) +
  //   sizeof(orient_84) + sizeof(maps_84);

  // cerr << " Total spurious memory consumption is " << total/1.e6 << " Mb " << endl;

  // //clear the mappings
  // for(unsigned int i=0; i<number_of_orientations; i++){

  //   for(unsigned int j=0; j<16; j++){
  //     delete [] orient_16[i][j];
  //   }
  //   for(unsigned int j=0; j<21; j++){
  //     delete [] orient_21[i][j];
  //   }
  //   for(unsigned int j=0; j<32; j++){
  //     delete [] orient_32[i][j];
  //   }
  //   for(unsigned int j=0; j<42; j++){
  //     delete [] orient_42[i][j];
  //   }
  //   for(unsigned int j=0; j<64; j++){
  //     delete [] orient_64[i][j];
  //   }
  //   for(unsigned int j=0; j<84; j++){
  //     delete [] orient_84[i][j];
  //   }

  //   delete [] orient_16[i];
  //   delete [] orient_21[i];
  //   delete [] orient_32[i];
  //   delete [] orient_42[i];
  //   delete [] orient_64[i];
  //   delete [] orient_84[i];

  // }

  // delete [] orient_16;
  // delete [] orient_21;
  // delete [] orient_32;
  // delete [] orient_42;
  // delete [] orient_64;
  // delete [] orient_84;

}

/**** Create a simple homogeneous point distribution using the gsl random number generator  ****/
void SimpleX::poisson_square() {

  //create temporary structure to store vertex
  Vertex tempVert; 
  //make sure the vertices vector is completely empty
  vertices.clear();
  vector< Vertex >().swap( vertices );

  //first create a source at the user specified position
  //first source point is at the exact source position
  tempVert.set_x( (float) source_x );
  tempVert.set_y( (float) source_y );
  tempVert.set_z( (float) source_z );
  tempVert.set_vertex_id( (unsigned long long int) 0 );
  tempVert.set_border( 0 );
  tempVert.set_process( 0 );

  vertices.push_back( tempVert );

  //in the case of more than one source point, the rest is placed in sphere around this point,
  //at user-specified radius
  for(unsigned int i=1;i<source_points;i++){

    //draw random positions
    double x = gsl_rng_uniform( ran );
    double y = gsl_rng_uniform( ran );
    double z = gsl_rng_uniform( ran );

    //check if this position is inside the specified radius, if so, at to vertex list
    double radius = sqrt( pow( x - source_x , 2 ) + pow( y - source_y, 2 ) + pow( z - source_z, 2 ) );  
    if( radius <= source_radius ){

      tempVert.set_x( (float) x );
      tempVert.set_y( (float) y );
      tempVert.set_z( (float) z );
      tempVert.set_vertex_id( (unsigned long long int) i );
      tempVert.set_border( 0 );
      tempVert.set_process( 0 );

      vertices.push_back( tempVert );

    }else{
      i--;
    }

  }

  //now continue with the 'non-source' points
  for( unsigned long long int i=source_points; i<numSites; i++ ){

    //positions are random between 0 and 1
    double x = gsl_rng_uniform( ran );
    double y = gsl_rng_uniform( ran );
    double z = gsl_rng_uniform( ran );

    tempVert.set_x( (float) x );
    tempVert.set_y( (float) y );
    tempVert.set_z( (float) z );
    tempVert.set_vertex_id( (unsigned long long int) i );
    tempVert.set_border( 0 );
    tempVert.set_process( 0 );

    vertices.push_back( tempVert );

  }

  // ofstream vertex_output;
  // vertex_output.open( "vertices.txt" );
  // vertex_output << "id    x    y   z   n_H   flux   x_ion" << endl;
  // vector<Vertex>::iterator it=vertices.begin();
  // vertex_output << it->get_vertex_id() << " " << it->get_x() << " " << it->get_y() << " " << it->get_z() << " " << 1.0 << " " << 1.0 << " " << 0.0 << endl;
  // it++; 
  // while( it!=vertices.end() ){
  //   vertex_output << it->get_vertex_id() << " " << it->get_x() << " " << it->get_y() << " " << it->get_z() << " " << 1.0 << " " << 0.0 << " " << 0.0 << endl;
  //   it++;
  // }


  if( COMM_RANK == 0 ){
    simpleXlog << "  Homogeneous point distribution created " << endl;
  }

}

/****  Read in vertices from hdf5 file  ****/
void SimpleX::read_vertex_list(){

  //open hdf5 file
  char fileName[200];
  sprintf(fileName, inputDensityFile.c_str());
  h5w file(fileName,'o');        // open file for reading

  //read in total number of sites
  unsigned int temp;
  file.read_attr("/Header","number_of_sites", &temp); 
  numSites = (unsigned long long int) temp;  

  //read in conversion factors for dimensionless units 
  file.read_attr("/Header","UNIT_D",&UNIT_D);
  file.read_attr("/Header","UNIT_I",&UNIT_I);
  file.read_attr("/Header","UNIT_L",&UNIT_L);
  file.read_attr("/Header","UNIT_M",&UNIT_M);
  UNIT_V = pow( UNIT_L, 3.0 );


  //read in the size of the simulation box
  file.read_attr("/Header","box_size", &sizeBox);

  //create vector containing the vertices, and the temporary lists containing fluxes and masses
  //the last two are necessary because they are needed in check_undersampled()
  vertices.resize(numSites);
  temp_n_HI_list.resize(numSites);
  temp_flux_list.resize(numSites);
  temp_n_HII_list.resize(numSites);

  //structures needed for reading in the values from hdf5
  arr_1D<float> double_arr;
  arr_1D<unsigned int> int_arr;

  //arrays holding the dimensions of the data and the offset in case the data is read in in chunks
  int dims[2], offset[2];

  unsigned int chunk_size_read = chunk_size;  
  if (chunk_size_read > numSites){
    chunk_size_read = numSites;
  }

  // do the reading !!!
  offset[1] = 0;
  for( unsigned long long int i=0; i<numSites; i+=chunk_size_read ){

    offset[0] = i;
    if (i+chunk_size_read >= numSites) // make sure not to write outside of data range
      dims[0] = numSites - i;
    else
      dims[0] = chunk_size_read;

    // coordinates
    dims[1] = 3; 
    double_arr.reinit(2,dims);
    file.read_data("/Vertices/coordinates",offset, &double_arr);
    for(int j=0; j<dims[0]; j++){
      vertices[ j + i ].set_x( double_arr(j,0) ); // coordinates
      vertices[ j + i ].set_y( double_arr(j,1) ); 
      vertices[ j + i ].set_z( double_arr(j,2) );   
    }

    dims[1]=1;
    //number density
    double_arr.reinit(1,dims);
    file.read_data("Vertices/n_HI",offset, &double_arr);
    for( int j=0; j<dims[0]; j++ ){
      temp_n_HI_list[ j + i ] = double_arr(j);
    }

    //ionised fraction
    double_arr.reinit(1,dims);
    file.read_data("Vertices/n_HII",offset, &double_arr);
    for( int j=0; j<dims[0]; j++ ){
      temp_n_HII_list[ j + i ] = double_arr(j);
    }

    //flux
    double_arr.reinit(1,dims);
    file.read_data("Vertices/flux",offset, &double_arr);
    for( int j=0; j<dims[0]; j++ ){
      temp_flux_list[ j + i ] = double_arr(j);
    }


  }//for all sites to read

  //set the vertex id
  //perhaps in the future this should be read in as well?
  unsigned long long int i=0;
  for( VERTEX_ITERATOR it=vertices.begin(); it!=vertices.end(); it++, i++ ){
    it->set_vertex_id( i );
  }

  origNumSites = vertices.size();
  vertex_id_max = vertices.size() + borderSites;

  file.close();

  if(COMM_RANK == 0){
    simpleXlog << "  Read " << numSites << " sites from file " << endl;
    simpleXlog << "  Size of the simulation domain: " << sizeBox << " pc" << endl;
  }

}


/**** Create homogenous boundary around simulation domain ****/
// Width of the boundary is variable borderBox
// Number of points in the boundary is variable borderSites
void SimpleX::create_boundary(){

  //coordinates of border point
  double x1,x2,x3;

  //add the number of border sites to be added to the number of sites
  numSites+=borderSites;

  //loop over the extra sites
  bool stop;
  for( unsigned long long int i=0; i<borderSites; i++) {
    stop=0;
    //loop until one point in border is found
    while(!stop) {

      //random numbers between 0 and 1
      x1 = gsl_rng_uniform( ran );
      x2 = gsl_rng_uniform( ran );
      x3 = gsl_rng_uniform( ran );

      //convert from random number between 0 and 1 to between -border and 1.0+border
      x1 = x1 * (1.0 + 2 * borderBox) - borderBox;
      x2 = x2 * (1.0 + 2 * borderBox) - borderBox;
      x3 = x3 * (1.0 + 2 * borderBox) - borderBox;


      //if the point is not in de domain itself (so in border), stop
      //and use this point
      if(x1<0.0||x1>1.0||x2<0.0||x2>1.0||x3<0.0||x3>1.0) {
        stop=1;
      }
    }

    //add the thus acquired point to list of vertices
    Vertex tempVert; 

    tempVert.set_x( x1 );
    tempVert.set_y( x2 );
    tempVert.set_z( x3 );
    tempVert.set_vertex_id( (unsigned long long int ) origNumSites+i );
    tempVert.set_border(1);

    vertices.push_back( tempVert );

  }//for all border sites

  if( COMM_RANK == 0 )
    simpleXlog << "  Boundary points placed " << endl;

}

/**** Create periodic boundary around simulation domain ****/
// Width of the boundary in which periodicity
// is determined is variable borderBox
void SimpleX::create_periodic_boundary(){

  unsigned int count = 0;

  //loop over permutations of -1 0 and 1 in 3 dimensions
  for( short int i=-1; i<=1; i++ ) {
    for( short int j=-1; j<=1; j++ ) {
      for( short int k=-1; k<=1; k++ ) {
        for( unsigned int v = 0; v < numSites; v++ ) {
          //exclude the case i=j=k=0, that's the simulation domain itself
          if( (i||j||k) &&
            ( vertices[v].get_x() + 1.0*i ) >= -borderBox && ( vertices[v].get_x() + 1.0*i ) <= ( 1.0 + borderBox ) &&
            ( vertices[v].get_y() + 1.0*j ) >= -borderBox && ( vertices[v].get_y() + 1.0*j ) <= ( 1.0 + borderBox ) &&
            ( vertices[v].get_z() + 1.0*k ) >= -borderBox && ( vertices[v].get_z() + 1.0*k ) <= ( 1.0 + borderBox ) ){

            Vertex tempVert;
            tempVert = vertices[v];

            tempVert.set_x( vertices[v].get_x() + 1.0*i );       
            tempVert.set_y( vertices[v].get_y() + 1.0*j );       
            tempVert.set_z( vertices[v].get_z() + 1.0*k );       

            //make sure the vertex id is not 0!
            tempVert.set_border( vertices[v].get_vertex_id() + 1 );

            tempVert.set_vertex_id( numSites + count );

            vertices.push_back(tempVert);

            count++;
          }     

        }
      }
    }
  }

  //set the correct vertex_id_max
  vertex_id_max = numSites + count;
  numSites += count;

  if( COMM_RANK == 0 ){
    simpleXlog << "  Created periodic boundary around unity domain" << endl;
  }

}


/****  Create octree of vertices  ****/
//size of octree is 4 * #subboxes, to ensure fast localisation of vertices
//in subbox and subbox boundaries
void SimpleX::create_vertex_tree(){

  //make sure tree is empty
  vert_tree.delete_octree();

  //number of subboxes in one dimension
  unsigned int hilbert_resolution = pow( 2, hilbert_order );
  //tree size is 4*number of subboxes (in one dimension) to make search for boundary vertices most efficient
  unsigned int tree_size = 4*hilbert_resolution;

  //initilise the tree with this computed size
  double temp_borderBox = (double) borderBox;
  vert_tree.init_octree( tree_size, temp_borderBox );

  //put the vertices in the tree
  for( unsigned long long int i=0; i<vertices.size(); i++ ){
    Vertex temp = vertices[i];
    //it is essential that the vertex_id of the vertex in teh tree points
    //to the place in the vertices vector, which is not the case after
    //one run since unused boundary points are deleted
    temp.set_vertex_id(i);
    vert_tree.insert_vertex( temp.get_x(), temp.get_y(), temp.get_z(), temp );
  }

  if( COMM_RANK == 0 ){
    simpleXlog << "  Vertex tree created " << endl;
  }

}


/****  Decompose domain  ****/
// Count the number of sites contained in one subbox and add those until 
// number of points is approximately equal per processor
void SimpleX::decompose_domain(){

  //number of subboxes is number of cells in hilbert curve
  unsigned int number_of_subboxes = pow( 2, hilbert_order*dimension );

  //hilbert resolution is the number of cells in one direction
  unsigned int hilbert_resolution = pow( 2, hilbert_order );

  //width of the subbox is 1/number of cells in one dimension
  double subbox_width = 1.0/hilbert_resolution;

  //make sure dom_dec is empty
  dom_dec.clear();

  //array that will hold the process number for every subbox
  dom_dec.resize(number_of_subboxes);

  //work load is approximately equal for every site, so divide approx equal number of points over procs
  unsigned int numSites_per_proc = unsigned( ceil(  (double) vertices.size()/(COMM_SIZE)  ));

  //current proc
  unsigned int proc = 0;
  //number of points already on current proc
  unsigned int total_points_on_proc = 0;

  //loop over all subboxes and count the number of vertices in each box 
  //to determine the domain decomposition
  for( unsigned int i=0; i<number_of_subboxes; i++ ){

    //coordinates of the hilbert cell
    unsigned long long int coord[dimension];
    //hilbert number (place on hilbert curve)
    unsigned long long r = (unsigned long long) i;

    //find coordinates of subbox number r
    hilbert_i2c( dimension, hilbert_order, r, coord );

    //calculate minimum and maximum coordinates of subbox and boundary around it
    //minimum coordinates of subbox should not be bigger than boundary around unity domain
    double x_min = ( coord[0] == 0 ) ? (0.0 - borderBox) : (double) coord[0]/hilbert_resolution;
    double x_max = ( (double) coord[0]/hilbert_resolution + subbox_width  >= 1.0  ) ? 
      (1.0 + borderBox) : (double) coord[0]/hilbert_resolution + subbox_width;

    double y_min = ( coord[1] == 0 ) ? (0.0 - borderBox) : (double) coord[1]/hilbert_resolution;
    double y_max = ( (double) coord[1]/hilbert_resolution + subbox_width  >= 1.0  ) ? 
      (1.0 + borderBox) : (double) coord[1]/hilbert_resolution + subbox_width;

    double z_min = ( coord[2] == 0 ) ? (0.0 - borderBox) : (double) coord[2]/hilbert_resolution;
    double z_max = ( (double) coord[2]/hilbert_resolution + subbox_width  >= 1.0  ) ? 
      (1.0 + borderBox) : (double) coord[2]/hilbert_resolution + subbox_width;

    //find id's of vertices that lie in this subbox
    vector< unsigned long long int > in_box = vert_tree.find_vertices_in_box( x_min, x_max, y_min, y_max, z_min, z_max );

    //add the number of vertices to the total on this proc
    total_points_on_proc += in_box.size();

    //assign this subbox to current proc
    dom_dec[i] = proc;

    //if the total number of points in this subbox exceeds 
    //the number of points assigned to this proc, continue
    //with next proc
    if( total_points_on_proc >= ( numSites_per_proc - (numSites_per_proc/number_of_subboxes) ) ){
      total_points_on_proc = 0;
      proc++;
      //if this is already the last proc, stay on this proc
      if( proc >= COMM_SIZE ){
        proc = COMM_SIZE - 1;
      }
    }  

  }//for all subboxes

  if( dom_dec.size() < COMM_SIZE ){
    cerr << " (" << COMM_RANK << ") Error: not enough subboxes for domain decomposition " << endl;
    MPI::COMM_WORLD.Abort( -1 );
  }

  bool correct = 0;
  unsigned int i = 1;
  while( !correct ){

    if( dom_dec[ dom_dec.size() - i ] < (COMM_SIZE - i ) ){
      dom_dec[ dom_dec.size() - i ] = COMM_SIZE - i;
    }else{
      correct = 1;
    }
    i++;

  }

  if( COMM_RANK == 0 ){
    simpleXlog << "  Domain decomposed " << endl;
    simpleXlog << "      Domain decomposition: " << endl;
    for( unsigned int i=0; i<number_of_subboxes; i++ ){
      simpleXlog << "        Subbox: " << i << " proc: " << dom_dec[i] << endl; 
    } 
    simpleXlog << endl;
  }

}

void SimpleX::assign_process(){


  //hilbert resolution is the number of cells in one direction
  unsigned int hilbert_resolution = pow( 2, hilbert_order );

  //width of the subbox is 1/number of cells in one dimension
  double subbox_width = 1.0/hilbert_resolution;

  //loop over all subboxes
  for( unsigned int i=0; i<dom_dec.size(); i++ ){

    //coordinates of the hilbert cell
    unsigned long long int coord[dimension];
    //hilbert number (place on hilbert curve)
    unsigned long long r = (unsigned long long) i;

    //find coordinates of subbox number r
    hilbert_i2c( dimension, hilbert_order, r, coord );

    //calculate minimum and maximum coordinates of subbox and boundary around it
    //minimum coordinates of subbox should not be bigger than boundary around unity domain
    double x_min = ( coord[0] == 0 ) ? (0.0 - borderBox) : (double) coord[0]/hilbert_resolution;
    double x_max = ( (double) coord[0]/hilbert_resolution + subbox_width  >= 1.0  ) ? 
      (1.0 + borderBox) : (double) coord[0]/hilbert_resolution + subbox_width;

    double y_min = ( coord[1] == 0 ) ? (0.0 - borderBox) : (double) coord[1]/hilbert_resolution;
    double y_max = ( (double) coord[1]/hilbert_resolution + subbox_width  >= 1.0  ) ? 
      (1.0 + borderBox) : (double) coord[1]/hilbert_resolution + subbox_width;

    double z_min = ( coord[2] == 0 ) ? (0.0 - borderBox) : (double) coord[2]/hilbert_resolution;
    double z_max = ( (double) coord[2]/hilbert_resolution + subbox_width  >= 1.0  ) ? 
      (1.0 + borderBox) : (double) coord[2]/hilbert_resolution + subbox_width;

    //find id's of vertices that lie in this subbox
    vector< unsigned long long int > in_box = vert_tree.find_vertices_in_box( x_min, x_max, y_min, y_max, z_min, z_max );

    //assign all vertices in this subbox to this proc
    for( vector< unsigned long long int >::iterator it=in_box.begin(); it!=in_box.end(); it++ ){
      if( *it > vertices.size() ){
        cerr << " in_box gives: " << *it << " number of vertices: " << vertices.size() << endl;
        MPI::COMM_WORLD.Abort( -1 );
      }else{
        vertices[ *it ].set_process( dom_dec[i] ); 
      }
    }

  }//for all subboxes

}

//check if vertex is in domain
bool SimpleX::inDomain(const float& x, const float& y, const float& z, const unsigned int& rank){
//bool SimpleX::inDomain(float x, float y, float z, unsigned int rank){
  
  bool in_domain = 0;
  unsigned int hilbert_resolution = pow( 2, hilbert_order );


  //determine subbox the vertex is in

  //first find the hilbert coordinates of the vertex
  int x_hilbert = (int) floor(x*hilbert_resolution);
  if( x_hilbert < 0 ){
    x_hilbert = 0;
  }else if( x_hilbert >= (int) hilbert_resolution ){
    x_hilbert = hilbert_resolution - 1;
  }
  int y_hilbert = (int) floor(y*hilbert_resolution); 
  if( y_hilbert < 0 ){
    y_hilbert = 0;
  }else if( y_hilbert >= (int) hilbert_resolution ){
    y_hilbert = hilbert_resolution - 1;
  }
  int z_hilbert = (int) floor(z*hilbert_resolution);
  if( z_hilbert < 0 ){
    z_hilbert = 0;
  }else if( z_hilbert >= (int) hilbert_resolution ){
    z_hilbert = hilbert_resolution - 1;
  }

  unsigned long long int coord_vertex[dimension];
  coord_vertex[0] = (unsigned long long) x_hilbert;
  coord_vertex[1] = (unsigned long long) y_hilbert;
  coord_vertex[2] = (unsigned long long) z_hilbert;

  //determine the hilbert number of the subbox the vertex is in
  unsigned long long r_vertex = hilbert_c2i( dimension, hilbert_order, coord_vertex );

  if( dom_dec[ r_vertex ] == rank ){
    in_domain = 1;
  }

  return in_domain;

}

/****  Compute triangulation  ****/
//compute triangulation of every subbox separately and check whether the
//triangulation is valid, i.e. whether the boundaries around the subbox
//are sufficiently large
void SimpleX::compute_triangulation(){


  // ---------- define subboxes  ---------------//


  //number of subboxes is number of cells in hilbert curve
  unsigned int number_of_subboxes = pow( 2, hilbert_order*dimension );

  //hilbert resolution is the number of cells in one direction
  unsigned int hilbert_resolution = pow( 2, hilbert_order );

  //width of the subbox is 1/number of cells in one dimension
  double subbox_width = 1.0/hilbert_resolution;
  //border around subbox is padding_subbox*subbox width to begin with
  double border_subbox = padding_subbox*subbox_width;

  //number of subboxes that is assigned to this proc
  unsigned int subboxes_on_proc = 0;
  //first subbox on this proc on Hilbert curve
  unsigned int start_number = 0;

  for(unsigned int i=0; i<number_of_subboxes; i++){
    if(dom_dec[i] == COMM_RANK){
      start_number = i;
      subboxes_on_proc++;
    }
  }

  //compute correct starting subbox by substracting number of subboxes-1
  start_number -= subboxes_on_proc-1;

  //hilbert number of subbox
  unsigned long long this_subbox = (unsigned long long) start_number;

  simplices.clear();

  //loop over number of subboxes on this proc and triangulate and check for correct triangulation per subbox
  for( unsigned int i=0; i<subboxes_on_proc; i++ ){

    //vector to hold the simplices on this proc
    vector< Simpl > simplices_subbox;

    //flag that determines whether triangulation can be trusted
    bool correct = 0;

    //coordinates of the hilbert cell
    unsigned long long int coord[dimension];

    //find coordinates of this subbox
    hilbert_i2c( dimension, hilbert_order, this_subbox, coord );

    //minimum and maximum coordinates of subbox and boundary around it
    //minimum coordinates of subbox should not be bigger than boundary around unity domain
    double temp_x_min = (double) coord[0]/hilbert_resolution - border_subbox;
    //double x_min = (  temp_x_min < (0.0 - borderBox) ) ? (0.0 - borderBox) : temp_x_min;
    double x_min = (  temp_x_min < 0.0 ) ? (0.0 - borderBox) : temp_x_min;

    double temp_x_max = (double) coord[0]/hilbert_resolution + subbox_width + border_subbox;
    //double x_max = (  temp_x_max > (1.0 + borderBox) ) ? (1.0 + borderBox) : temp_x_max;
    double x_max = (  temp_x_max > 1.0 ) ? (1.0 + borderBox) : temp_x_max;

    double temp_y_min = (double) coord[1]/hilbert_resolution - border_subbox;
    //double y_min = (  temp_y_min < (0.0 - borderBox) ) ? (0.0 - borderBox) : temp_y_min;
    double y_min = (  temp_y_min < 0.0 ) ? (0.0 - borderBox) : temp_y_min;

    double temp_y_max = (double) coord[1]/hilbert_resolution + subbox_width + border_subbox;
    //double y_max = (  temp_y_max > (1.0 + borderBox) ) ? (1.0 + borderBox) : temp_y_max;
    double y_max = (  temp_y_max > 1.0 ) ? (1.0 + borderBox) : temp_y_max;

    double temp_z_min = (double) coord[2]/hilbert_resolution - border_subbox;
    //double z_min = (  temp_z_min < (0.0 - borderBox) ) ? (0.0 - borderBox) : temp_z_min;
    double z_min = (  temp_z_min < 0.0 ) ? (0.0 - borderBox) : temp_z_min;

    double temp_z_max = (double) coord[2]/hilbert_resolution + subbox_width + border_subbox;
    //double z_max = (  temp_z_max > (1.0 + borderBox) ) ? (1.0 + borderBox) : temp_z_max;
    double z_max = (  temp_z_max > 1.0 ) ? (1.0 + borderBox) : temp_z_max;

    //if the triangulation is not correct the first time, do it again
    //with new boundaries that are correct
    while(!correct){

      //empty the vector, it might have been filled in the previous iteration
      simplices_subbox.clear();

      //------ find vertices in subbox --------//

      //find id's of all vertices in subbox + border around subbox
      vector< unsigned long long int > in_box = vert_tree.find_vertices_in_box( x_min, x_max, y_min, y_max, z_min, z_max );

      //check whether there's enough points in subbox to do triangulation.
      unsigned int in_subbox = 0;
      for( unsigned int q=0; q<in_box.size(); q++){
        if( !vertices[ in_box[q] ].get_border() ){
          in_subbox++;
        }
      }

      // cerr << " (" << COMM_RANK << ") Subbox " << i << " contains " << in_subbox << " points" << endl
      // 	   << " Coordinates: (" << x_min << "," << x_max << ") (" << y_min << "," << y_max << ") (" << z_min << "," << z_max << ")";

      //change the abort to extension of the borders!
      //if( in_box.size() <= (unsigned int) dimension ){ 
      if( in_subbox <= (unsigned int) dimension ){ 
        cerr << "Too few points in subbox to do tessellation." << endl;
        MPI::COMM_WORLD.Abort( -1 );
      }

      //minimum and maximum coordinates of subbox without the boundary
      //minimum x-coordinate of the subbox for check
      double x_min_subbox = (double) coord[0]/hilbert_resolution;
      double x_min_subbox_check = (double) coord[0]/hilbert_resolution;
      //if x_min is equal to or smaller than 0.0, take into account outer boundary
      if( x_min_subbox <= 0.0 ){
        x_min_subbox = x_min ; 
        x_min_subbox_check = 0.0;
      }
      //maximum x-coordinate of subbox for check
      double x_max_subbox = (double) coord[0]/hilbert_resolution + subbox_width;
      double x_max_subbox_check = (double) coord[0]/hilbert_resolution + subbox_width;
      //if x_max is larger than or equal to 1.0, it's in the outer boundary, so take into account
      if( x_max_subbox >= 1.0 ){
        x_max_subbox = x_max;
        x_max_subbox_check = 1.0;
      }

      //minimum y-coordinate of the subbox for check
      double y_min_subbox = (double) coord[1]/hilbert_resolution;
      double y_min_subbox_check = (double) coord[1]/hilbert_resolution;
      //if y_min is equal to or smaller than 0.0, take into account outer boundary
      if( y_min_subbox <= 0.0 ){
        y_min_subbox = y_min;
        y_min_subbox_check = 0.0;;
      }
      //maximum y-coordinate of subbox for check
      double y_max_subbox = (double) coord[1]/hilbert_resolution + subbox_width;
      double y_max_subbox_check = (double) coord[1]/hilbert_resolution + subbox_width;
      //if y_max is larger than or equal to 1.0, it's in the outer boundary, so take into account
      if( y_max_subbox >= 1.0 ){
        y_max_subbox = y_max;
        y_max_subbox_check = 1.0;
      }

      //minimum z-coordinate of the subbox for check
      //if z_min is equal to or smaller than 0.0, take into account outer boundary
      double z_min_subbox = (double) coord[2]/hilbert_resolution;
      double z_min_subbox_check = (double) coord[2]/hilbert_resolution;
      if( z_min_subbox <= 0.0 ){
        z_min_subbox = z_min;
        z_min_subbox_check = 0.0;
      }
      //maximum z-coordinate of subbox for check
      double z_max_subbox = (double) coord[2]/hilbert_resolution + subbox_width;
      double z_max_subbox_check = (double) coord[2]/hilbert_resolution + subbox_width;
     //if z_max is larger than or equal to 1.0, it's in the outer boundary, so take into account
      if( z_max_subbox >= 1.0 ){
        z_max_subbox = z_max;
        z_max_subbox_check = 1.0;
      }

      // ---------- triangulate this subbox  ---------------//

      //variables needed for qhull call
      FILE *errfile=stderr;
      boolT ismalloc=False;
      char flags[250];
      sprintf(flags,"qhull d Qbb T0");
      facetT *facet;
      vertexT *vertex, **vertexp;
      int curlong, totlong;

      //arrays needed for qhull call
      int* idList=NULL;
      double* pt_array=NULL;
      pt_array=new double[dimension*in_box.size()];      
      idList=new int[dimension*in_box.size()];
      int ids[4];

      //fill arrays with coordinates and id's
      for( unsigned int q=0; q<in_box.size(); q++ ){
        //id's
        idList[q]=q;
        //coordinates
        pt_array[q*dimension]   = vertices[ in_box[q] ].get_x();
        pt_array[q*dimension+1] = vertices[ in_box[q] ].get_y();
        pt_array[q*dimension+2] = vertices[ in_box[q] ].get_z();
      }

      //assume for now the triangulation will be correct, check will follow
      correct = 1;      

      //call to qhull to do triangulation
      if (!qh_new_qhull(dimension, in_box.size(), pt_array, ismalloc, flags, NULL, errfile)) {
        //loop over all facets
        FORALLfacets {
          if (!facet->upperdelaunay) {

            //store the simplex in simplices vector
            unsigned int r=0;
            Simpl tempSimpl;
            tempSimpl.set_volume( qh_facetarea(facet) );

            FOREACHvertex_ (facet->vertices) {
              //store the indices of each simplex
              ids[r++]=idList[qh_pointid(vertex->point)];
            }//for each facet

            tempSimpl.set_id1( in_box[ ids[0] ] );
            tempSimpl.set_id2( in_box[ ids[1] ] );
            tempSimpl.set_id3( in_box[ ids[2] ] );
            tempSimpl.set_id4( in_box[ ids[3] ] );


      // ---------- check if simplex is in this subbox  ---------------//


            //coordinates of tetrahedron, to be used to calculate circumsphere
            double xTet[4][3];

            //store the coordinates of the 4 vertices that make 
            //up this simplex in xTet
            for(short int p=0; p<dimension+1; p++) {
              xTet[p][0] = (double) vertices[ in_box[ids[p]] ].get_x();
              xTet[p][1] = (double) vertices[ in_box[ids[p]] ].get_y();
              xTet[p][2] = (double) vertices[ in_box[ids[p]] ].get_z();
            }

            //coordinates of the centre of the circumsphere
            double xCC, yCC, zCC;

            //calculate the centers of the circumsphere
            CalcCircumCenter(xTet, xCC, yCC, zCC);

      //check if the centre of the circumsphere is inside this subbox
      //in that case, the simplex belongs to this subbox and will be
      //added to the list of simplices. The boundary around the unity domain 
      //has to be take into account, while the subbox boundaries inside the 
      //domain have to be discarded

            bool xOK=0;
            bool yOK=0;
            bool zOK=0;


           //check whether centre of circumsphere is inside the boundaries defined above
            if( xCC >= (x_min_subbox) && xCC < (x_max_subbox) ){
              xOK = 1;
            }
            if( yCC >= (y_min_subbox) && yCC < (y_max_subbox) ){
              yOK = 1;
            }
            if( zCC >= (z_min_subbox) && zCC < (z_max_subbox) ){
              zOK = 1;
            }

            //if all coordinates are inside subbox, add simplex to list if it's not in boundary
            if( xOK && yOK && zOK){

              //exclude simplices that are entirely inside the boundary
              if( !vertices[ in_box[ids[0]] ].get_border() || !vertices[ in_box[ids[1]] ].get_border() ||
              !vertices[ in_box[ids[2]] ].get_border() || !vertices[ in_box[ids[3]] ].get_border() ){

                //add simplex to list
                simplices_subbox.push_back(tempSimpl);

              } //if one vertex not in border

            }else{ 

        //if simplex is not in subbox it might have to be taken into account if
        //the simplex is on another proc, but only if it has not already been added
        //by another subbox on this proc


        //determine subbox the simplex is in, by calculating the hilbert number 
        //of the subbox the circumcentre of this simplex is in

              //first find the hilbert coordinates of the circumcentre of the simplex
              int x_hilbert = (int) floor(xCC*hilbert_resolution);
              if( x_hilbert < 0 ){
                x_hilbert = 0;
              }else if( x_hilbert >= (int) hilbert_resolution ){
                x_hilbert = hilbert_resolution - 1;
              }
              int y_hilbert = (int) floor(yCC*hilbert_resolution); 
              if( y_hilbert < 0 ){
                y_hilbert = 0;
              }else if( y_hilbert >= (int) hilbert_resolution ){
                y_hilbert = hilbert_resolution - 1;
              }
              int z_hilbert = (int) floor(zCC*hilbert_resolution);
              if( z_hilbert < 0 ){
                z_hilbert = 0;
              }else if( z_hilbert >= (int) hilbert_resolution ){
                z_hilbert = hilbert_resolution - 1;
              }

              unsigned long long int coord_CC[dimension];
              coord_CC[0] = (unsigned long long) x_hilbert;
              coord_CC[1] = (unsigned long long) y_hilbert;
              coord_CC[2] = (unsigned long long) z_hilbert;

              //determine the hilbert number of the subbox the circumcentre is in
              unsigned long long r_CC = hilbert_c2i( dimension, hilbert_order, coord_CC );

              //if this simplex belongs to other proc, determine the lowest hilbert number
              //of the vertices that are on this proc
              if( dom_dec[ r_CC ] != COMM_RANK ){

                short int hilbert_number[dimension+1];
    //calculate the hilbert number of subboxes the vertices that live on this proc are in
                for(short int p=0; p<dimension+1; p++) {
      //the vertex should be on this proc
                  if( vertices[ in_box[ids[p]] ].get_process() == COMM_RANK ){

        //calculate coordinates in units of the hilbert resolution
                    int x_vert = (int) floor( vertices[ in_box[ids[p]] ].get_x()*hilbert_resolution);
                    if( x_vert < 0 ){
                      x_vert = 0;
                    }else if( x_vert >= (int) hilbert_resolution ){
                      x_vert = hilbert_resolution - 1;
                    }
                    int y_vert = (int) floor( vertices[ in_box[ids[p]] ].get_y()*hilbert_resolution);
                    if( y_vert < 0 ){
                      y_vert = 0;
                    }else if( y_vert >= (int) hilbert_resolution ){
                      y_vert = hilbert_resolution - 1;
                    }
                    int z_vert = (int) floor( vertices[ in_box[ids[p]] ].get_z()*hilbert_resolution);
                    if( z_vert < 0 ){
                      z_vert = 0;
                    }else if( z_vert >= (int) hilbert_resolution ){
                      z_vert = hilbert_resolution - 1;
                    }

                    unsigned long long int coord_vert[dimension];
                    coord_vert[0] = (unsigned long long) x_vert;
                    coord_vert[1] = (unsigned long long) y_vert;
                    coord_vert[2] = (unsigned long long) z_vert;

        //calculate the hilbert number of the subbox the vertex is in
                    unsigned long long r_vert = hilbert_c2i( dimension, hilbert_order, coord_vert );

                    //store the hilbert number
                    hilbert_number[p] = r_vert;

                  }else{
                    //if the vertex is not on this proc, it shouldn't be the minimum
                    hilbert_number[p] = number_of_subboxes+1;
                  }

                }//for all vertices in this simplex

                //determine the minimum hilbert number of vertices on this proc
                unsigned int min = *min_element( hilbert_number, hilbert_number + 4);


                //if this subbox has the lowest hilbert number that contains a vertex on this
                //proc, include the simplex 
                if( min == this_subbox ){
                 //exclude simplices that are entirely inside the boundary
                  if( !vertices[ in_box[ids[0]] ].get_border() || !vertices[ in_box[ids[1]] ].get_border() ||
                  !vertices[ in_box[ids[2]] ].get_border() || !vertices[ in_box[ids[3]] ].get_border() ){

                    simplices_subbox.push_back(tempSimpl);
                  }
                }

              }//if simplex not on this proc

            }//if simplex in this subbox

      //check if boundary is sufficiently large
      //if at least one of the vertices is in subbox, check the boundary
            bool vertex_in_subbox = 0;
            for(short int p=0; !vertex_in_subbox && p<dimension+1; p++) {

              bool x_in = 0;
              bool y_in = 0;
              bool z_in = 0;

        //check whether vertex coordinates are inside the subbox
              if( xTet[p][0] >= (x_min_subbox_check) && xTet[p][0] < (x_max_subbox_check) ){
                x_in = 1;
              }
              if( xTet[p][1] >= (y_min_subbox_check) && xTet[p][1] < (y_max_subbox_check) ){
                y_in = 1;
              }
              if( xTet[p][2] >= (z_min_subbox_check) && xTet[p][2] < (z_max_subbox_check) ){
                z_in = 1;
              }

              if( x_in && y_in && z_in ){
                vertex_in_subbox = 1;
              }
            }

      //if one of the vertices is inside the subbox, check the boundary
            if(vertex_in_subbox){
    // ---------- check whether the boundaries are sufficiently large ---------------//

    //calculate the radius of the circumsphere
              double radiusCC = sqrt( pow( xCC - (double) vertices[ in_box[ ids[0] ] ].get_x(), 2) +
                pow( yCC - (double) vertices[ in_box[ ids[0] ] ].get_y(), 2) + 
                pow( zCC - (double) vertices[ in_box[ ids[0] ] ].get_z(), 2) );

    //the radius of the circumsphere should be inside the boundaries on all sides
    //if not, extend x_min and x_max
              if( (xCC - radiusCC) < x_min ){
                x_min = xCC - radiusCC;
                correct = 0;
              }
              if( (xCC + radiusCC) > x_max ){
                x_max = xCC + radiusCC;
                correct = 0;
              }
              if( (yCC - radiusCC) < y_min ){
                y_min = yCC - radiusCC;
                correct = 0;
              }
              if( (yCC + radiusCC) > y_max ){
                y_max = yCC + radiusCC;
                correct = 0;
              }
              if( (zCC - radiusCC) < z_min ){
                z_min = zCC - radiusCC;
                correct = 0;
              }
              if( (zCC + radiusCC) > z_max ){
                z_max = zCC + radiusCC;
                correct = 0;
              }

              if( x_min < (0.0 - borderBox) || x_max > (1.0 + borderBox) || 
                y_min < (0.0 - borderBox) || y_max > (1.0 + borderBox) || 
              z_min < (0.0 - borderBox) || z_max > (1.0 + borderBox) ){

      // cerr << " (" << COMM_RANK << ") Simplex: (" 
      //      << vertices[ in_box[ids[0]] ].get_x() << "," 
      //      << vertices[ in_box[ids[0]] ].get_y() << "," 
      //      << vertices[ in_box[ids[0]] ].get_z() << ")" << endl << " (" 
      //      << vertices[ in_box[ids[1]] ].get_x() << "," 
      //      << vertices[ in_box[ids[1]] ].get_y() << "," 
      //      << vertices[ in_box[ids[1]] ].get_z() << ")" << endl << " (" 
      //      << vertices[ in_box[ids[2]] ].get_x() << "," 
      //      << vertices[ in_box[ids[2]] ].get_y() << "," 
      //      << vertices[ in_box[ids[2]] ].get_z() << ")" << endl << " (" 
      //      << vertices[ in_box[ids[3]] ].get_x() << "," 
      //      << vertices[ in_box[ids[3]] ].get_y() << "," 
      //      << vertices[ in_box[ids[3]] ].get_z() << ")" << endl 
      //      << " r_cc: " << radiusCC << endl;

      //MPI::COMM_WORLD.Abort(-1);
              }

            }//if at least one vertex is in subbox

          }//if upper delaunay
        }//for all facets
      }//if no qhull error    


      //check if boundaries aren't outside domain
      if( x_min < (0.0 - borderBox) || x_max > (1.0 + borderBox) || 
        y_min < (0.0 - borderBox) || y_max > (1.0 + borderBox) || 
      z_min < (0.0 - borderBox) || z_max > (1.0 + borderBox) ){

        cerr << endl << "  (" << COMM_RANK << ") Boundary chosen around unity domain does not contain enough points, exiting" << endl;
        cerr << "  (" << COMM_RANK << ")   box: (" << x_min << ", " << x_max << ")  (" 
          << y_min << ", " << y_max << ")  (" << z_min << ", " << z_max << ")" << endl;

        MPI::COMM_WORLD.Abort( -1 );

      }

      //free long memory
      qh_freeqhull(!qh_ALL);

      //free short memory and memory allocator
      qh_memfreeshort (&curlong, &totlong);
      if (curlong || totlong) {
        fprintf(errfile, "QHull: did not free %d bytes of long memory (%d pieces)", totlong, curlong);
      }

      //free dynamically allocated arrays
      if(pt_array) {
        delete [] pt_array;
        pt_array=NULL;
      }

      if(idList) {
        delete [] idList;
        idList=NULL;
      }

      if(!correct){
        cerr << " (" << COMM_RANK << ") Subbox " << i << " boundaries too small, retriangulating with subbox: (" << x_min << "," << x_max << ") (" 
          << y_min << "," << y_max << ") (" << z_min << "," << z_max << ")" << endl;
      }

    }//while triangulation is not correct

    simplices.insert( simplices.end(), simplices_subbox.begin(), simplices_subbox.end() );
    simplices_subbox.clear();

    //if all went well, increase subbox number and triangulate the next
    this_subbox++;
  }//for all subboxes on this proc


  if( COMM_RANK == 0 ){
    simpleXlog << "  Triangulation computed " << endl;
  }

}


/****  Create the sites array  ****/
//create a list of sites from the list of simplices that was created 
//during the triangulation stage
void SimpleX::create_sites(){

  //make sure the sites vector is empty before filling it
  sites.clear();

  //create a vector holding all the vertex indices
  //this vector will be used to check whether a vertex has already been added
  //this is necessary to avoid duplices because vertices belong to more than one simplex
  vector< unsigned int > indices( vertices.size(), vertices.size()+1 );

  //temporary Site structure
  Site tempSite;

    //loop over all simplices to extract the sites
  for( vector< Simpl >::iterator it=simplices.begin(); it!=simplices.end(); it++ ){

    //id's of the vertices, place in vertices vector
    unsigned int id1 = it->get_id1();
    unsigned int id2 = it->get_id2();
    unsigned int id3 = it->get_id3();
    unsigned int id4 = it->get_id4();

    //if site hasn't been added to list yet, the entry in indices is higher
    //than vertices.size()
    if( indices[id1] > vertices.size() ){
      //set indices to sites.size() to show site has been added, and to be 
      //used to correct the indices of the simplices
      indices[id1] = sites.size();
      tempSite = vertices[id1];
      //add to sites vector
      sites.push_back( tempSite );
    }

    if( indices[id2] > vertices.size() ){
      //set indices to sites.size() to show site has been added, and to be 
      //used to correct the indices of the simplices
      indices[id2] = sites.size();
      tempSite = vertices[id2];
      //add to sites vector
      sites.push_back( tempSite );
    }

    if( indices[id3] > vertices.size() ){
      //set indices to sites.size() to show site has been added, and to be 
      //used to correct the indices of the simplices
      indices[id3] = sites.size();
      tempSite = vertices[id3];
      //add to sites vector
      sites.push_back( tempSite );
    }

    if( indices[id4] > vertices.size() ){
      //set indices to sites.size() to show site has been added, and to be 
      //used to correct the indices of the simplices
      indices[id4] = sites.size();
      tempSite = vertices[id4];
      //add to sites vector
      sites.push_back( tempSite );
    }

    //the list of simplices should be up to date with
    //the new site indices
    it->set_id1( indices[id1] );
    it->set_id2( indices[id2] );
    it->set_id3( indices[id3] );
    it->set_id4( indices[id4] );

  }//for all simplices

  //It could be that not all border sites are taken into account, so 
  //the numSites might have changed
  unsigned int local_numSites = 0;
  for( SITE_ITERATOR it=sites.begin(); it!=sites.end(); it++){
    if( it->get_process() == COMM_RANK ){
      local_numSites++;
    }
  }

  MPI::COMM_WORLD.Allreduce( &local_numSites, &numSites, 1, MPI::UNSIGNED, MPI::SUM );

  if( COMM_RANK == 0 ){
    cerr << " (" << COMM_RANK << ") Number of sites in triangulation is " << numSites << endl;
    simpleXlog << "  Final triangulation contains " << numSites << " sites " << endl;
  }

//   if( COMM_RANK == 0 ){
//     cerr << " (" << COMM_RANK << ") Number of simplices in triangulation is " << simplices.size() << endl;
//   }

  //if periodic boundaries are applied, keep track of the site id of the corresponding vertex
  //do this only for sites on the same proc, others are kept for communication
  // if(periodic){
  //   for(SITE_ITERATOR it=sites.begin(); it!=sites.end(); it++){
  //     if(it->get_border()){
  // 	unsigned int vert_id = it->get_border() - 1;
  // 	if( vert_id < sites.size() ){
  // 	  if(it->get_process() == sites[ indices[ vert_id ] ].get_process()){
  // 	    it->set_border( indices[ vert_id ] + 1 );
  // 	  }
  // 	}
  //     }
  //   }
  // }

  //vertices vector no longer needed
  vertices.clear();
  vector< Vertex >().swap(vertices);

  //clear indices vector
  indices.clear();
  vector< unsigned int >().swap(indices);

  if( COMM_RANK == 0 ){
    simpleXlog << "  Sites array created" << endl;
  }

}

/****  Remove periodic boundary sites  ****/
void SimpleX::remove_periodic_sites(){
  
  //if multiple procs are used, use periodic sites to communicate between procs
  if(COMM_SIZE > 1){

    //store the boundaries of the domains
    vector< vector<double> > dom_bounds( dom_dec.size(), vector<double>(6,0.0) );

    //hilbert resolution is the number of cells in one direction
    unsigned int hilbert_resolution = pow( 2, hilbert_order );

    //width of the subbox is 1/number of cells in one dimension
    double subbox_width = 1.0/hilbert_resolution;

    //loop over all subboxes
    for( unsigned int i=0; i<dom_dec.size(); i++ ){

      //coordinates of the hilbert cell
      unsigned long long int coord[dimension];
      //hilbert number (place on hilbert curve)
      unsigned long long r = (unsigned long long) i;

      //find coordinates of subbox number r
      hilbert_i2c( dimension, hilbert_order, r, coord );

      //calculate minimum and maximum coordinates of subbox and boundary around it
      //minimum coordinates of subbox should not be bigger than boundary around unity domain
      double x_min = ( coord[0] == 0 ) ? 0.0  : (double) coord[0]/hilbert_resolution;
      double x_max = ( (double) coord[0]/hilbert_resolution + subbox_width  >= 1.0  ) ? 
        1.0  : (double) coord[0]/hilbert_resolution + subbox_width;

      double y_min = ( coord[1] == 0 ) ? 0.0 : (double) coord[1]/hilbert_resolution;
      double y_max = ( (double) coord[1]/hilbert_resolution + subbox_width  >= 1.0  ) ? 
        1.0  : (double) coord[1]/hilbert_resolution + subbox_width;

      double z_min = ( coord[2] == 0 ) ? 0.0  : (double) coord[2]/hilbert_resolution;
      double z_max = ( (double) coord[2]/hilbert_resolution + subbox_width  >= 1.0  ) ? 
        1.0  : (double) coord[2]/hilbert_resolution + subbox_width;

      dom_bounds[i][0] = x_min;
      dom_bounds[i][1] = x_max;
      dom_bounds[i][2] = y_min;
      dom_bounds[i][3] = y_max;
      dom_bounds[i][4] = z_min;
      dom_bounds[i][5] = z_max;

    }

    //point sites on other proc towards correct process
    unsigned int border_point = 0;
    unsigned int total_border_points = 0;
    for( SITE_ITERATOR it=sites.begin(); it!=sites.end(); it++ ){
      if( it->get_border() ){
      //if( it->get_border() && it->get_process() == COMM_RANK){
        //do reverse permutations
        for( short int i=-1; i<=1; i++ ) {
          for( short int j=-1; j<=1; j++ ) {
            for( short int k=-1; k<=1; k++ ) {
              if( (i||j||k) &&
                ( it->get_x() + 1.0*i ) >= 0.0 && ( it->get_x() + 1.0*i ) <= 1.0  &&
                ( it->get_y() + 1.0*j ) >= 0.0 && ( it->get_y() + 1.0*j ) <= 1.0  &&
                ( it->get_z() + 1.0*k ) >= 0.0 && ( it->get_z() + 1.0*k ) <= 1.0 ){

                double x = it->get_x() + 1.0*i;
                double y = it->get_y() + 1.0*j;
                double z = it->get_z() + 1.0*k;

                //determine process
                unsigned int process = COMM_SIZE;
                int found = 0;
                for(unsigned int p=0; p<dom_dec.size(); p++){
                  if( x >= dom_bounds[p][0] && x < dom_bounds[p][1] &&
                      y >= dom_bounds[p][2] && y < dom_bounds[p][3] &&
                      z >= dom_bounds[p][4] && z < dom_bounds[p][5] ){

                    process = dom_dec[p];
                    found++;
                  }
                }

                //check if found
                if(process == COMM_SIZE){
                  cerr << " (" << COMM_RANK << ") Error in assignment of periodic sites to processes!" << endl;
                  MPI::COMM_WORLD.Abort(-1);
                }

                if(found > 1){
                  cerr << " (" << COMM_RANK << ") Error in assignment of periodic sites to processes, found multiple processes!" << endl;
                  MPI::COMM_WORLD.Abort(-1);
                }

                //if(process != it->get_process() ){
                if(process != COMM_RANK){
                  it->set_vertex_id( it->get_border() - 1 );
                  it->set_border( 0 );
                  it->set_process( process );

                  border_point++;
                }

                total_border_points++;

              } //if    
            }//k
          }//j
        }//i

      }//if in border
    }//for all sites

    cerr << " (" << COMM_RANK << ") number of border points that changed process: " << border_point << " out of " << total_border_points << endl;

    dom_bounds.clear();
  }
  
  //it could be that sites from other proc are twice in the list, so remove double sites
  
  //set current site id
  unsigned int i=0;
  for(SITE_ITERATOR it=sites.begin(); it!=sites.end(); it++, i++){
    it->set_site_id(i);
  }
    
  //sort on vertex id
  //sort the sites on vertex id
  sort( sites.begin(), sites.end(), compare_vertex_id_site );

  //vector to keep track of the indices
  vector<unsigned long long int> indices( sites.size(), 0 );

  //store the new place of the sites in the indices array at the place of their 
  //former index
  for( unsigned long long int i=0; i<sites.size(); i++ ){
    indices[ sites[i].get_site_id() ] = i;
  }
  
  vector<Site> temp_sites = sites;
  sites.clear();
  
  //make sure the first site is accepted
  unsigned int vert_id = sites[0].get_vertex_id()+1;
  //number of rejected sites
  unsigned int rejected = 0;
  
  //vector to keep track of the indices
  vector<unsigned long long int> indices2( temp_sites.size(), 0 );
  
  i=0;
  unsigned int j=0;
  for(SITE_ITERATOR it=temp_sites.begin(); it!=temp_sites.end(); it++, i++){
    if(it->get_vertex_id() != vert_id){
      sites.push_back(*it);
      vert_id = it->get_vertex_id();
      j++;
    }else{
      rejected++;
    } 
    //place of temp_site in new sites vector
    indices2[i]=j-1;   
  }

  cerr << " (" << COMM_RANK << ") rejected " << rejected << " double sites" << endl;
  
  
  // if(COMM_RANK == 0){
  //   unsigned int k=100;
  //   cerr << " (" << COMM_RANK << ") " << sites[k].get_vertex_id() << " " << sites[indices2[ indices[sites[k].get_site_id()] ]].get_vertex_id() << endl;
  // }
  
  // if(COMM_RANK == 0){
  //   unsigned int k = 0;
  //   cerr << " (" << COMM_RANK << ") simplex " << k << ": " << compare_sites[simplices[k].get_id1()].get_vertex_id() << " " << compare_sites[simplices[k].get_id2()].get_vertex_id() << " "
  //        << compare_sites[simplices[k].get_id3()].get_vertex_id() << " " << compare_sites[simplices[k].get_id4()].get_vertex_id() << " " << endl;    
  // }
  
  //update the simplices
  for( SIMPL_ITERATOR it=simplices.begin(); it!=simplices.end(); it++ ){

     //former place of the sites in sites array
    unsigned int id1 = it->get_id1();
    unsigned int id2 = it->get_id2();
    unsigned int id3 = it->get_id3();
    unsigned int id4 = it->get_id4();

     //set the id to the new id of the sites
    it->set_id1(  indices2[ indices[id1] ] );
    it->set_id2(  indices2[ indices[id2] ] );
    it->set_id3(  indices2[ indices[id3] ] );
    it->set_id4(  indices2[ indices[id4] ] );

  }

  // if(COMM_RANK == 0){
  //   unsigned int k = 0;
  //   cerr << " (" << COMM_RANK << ") simplex " << k << ": " << sites[simplices[k].get_id1()].get_vertex_id() << " " << sites[simplices[k].get_id2()].get_vertex_id() << " "
  //        << sites[simplices[k].get_id3()].get_vertex_id() << " " << sites[simplices[k].get_id4()].get_vertex_id() << " " << endl;    
  // }

  //remove the indices vector
  indices.clear();
  indices2.clear();
  
  //set correct site id
  i=0;
  for( SITE_ITERATOR it=sites.begin(); it!=sites.end(); it++, i++ ){
    it->set_site_id(i);
  }  
  
  
  //-----------------
  //the only sites that are now in the border are double sites on this proc
  
  //find all the correct site ids
  for( SITE_ITERATOR it=sites.begin(); it!=sites.end(); it++ ){
    if(it->get_border() ){
      for( unsigned int i=0; i<sites.size(); i++ ){
        if( (it->get_border()-1) == sites[i].get_vertex_id() ){
          it->set_border( sites[i].get_site_id() + 1 );
        }
      }
    }
  }

  // //update the simplices to point to the correct sites
  for( SIMPL_ITERATOR it=simplices.begin(); it!=simplices.end(); it++ ){
  
     //place of the sites in sites array
     unsigned long long int id1 = it->get_id1();
     //if the site is in the border, point towards correct site
     if(sites[id1].get_border()){
       it->set_id1( sites[id1].get_border() - 1 );
     }
  
     //place of the sites in sites array
     unsigned long long int id2 = it->get_id2();
     //if the site is in the border, point towards correct site
     if(sites[id2].get_border()){
       it->set_id2( sites[id2].get_border() - 1 );
     }
     
     //place of the sites in sites array
     unsigned long long int id3 = it->get_id3();
     //if the site is in the border, point towards correct site
     if(sites[id3].get_border()){
       it->set_id3( sites[id3].get_border() - 1 );
     }
     
     //place of the sites in sites array
     unsigned long long int id4 = it->get_id4();
     //if the site is in the border, point towards correct site
     if(sites[id4].get_border()){
       it->set_id4( sites[id4].get_border() - 1 );
     }
  
   }
  
   //-- remove remaining border sites --//
  
     
    //vector to keep track of the indices
    indices.resize( sites.size(), 0 );
      
    //remove border sites
    temp_sites = sites;
    sites.clear();
    for( SITE_ITERATOR it=temp_sites.begin(); it!=temp_sites.end(); it++ ){
      if(!it->get_border()){
        sites.push_back( *it );
      }
    }
      
    temp_sites.clear();
      
    //store the new place of the sites in the indices array at the place of their 
    //former index
    for( unsigned long long int i=0; i<sites.size(); i++ ){
      indices[ sites[i].get_site_id() ] = i;
    }   
      
    //correct the indices of the simplices
    for( SIMPL_ITERATOR it=simplices.begin(); it!=simplices.end(); it++ ){
      
      //former place of the sites in sites array
      unsigned long long int id1 = it->get_id1();
      unsigned long long int id2 = it->get_id2();
      unsigned long long int id3 = it->get_id3();
      unsigned long long int id4 = it->get_id4();
      
      //set the id to the new id of the sites
      it->set_id1( (unsigned long long int) indices[id1] );
      it->set_id2( (unsigned long long int) indices[id2] );
      it->set_id3( (unsigned long long int) indices[id3] );
      it->set_id4( (unsigned long long int) indices[id4] );
      
    }
    //remove the indices vector
    indices.clear();
      
    //set correct site id
    i=0;
    for( SITE_ITERATOR it=sites.begin(); it!=sites.end(); it++, i++ ){
      it->set_site_id(i);
    }
      
    //It could be that not all border sites are taken into account, so 
    //the numSites might have changed
    unsigned int local_numSites = 0;
    for( SITE_ITERATOR it=sites.begin(); it!=sites.end(); it++){
      if( it->get_process() == COMM_RANK ){
        local_numSites++;
      }
    }
      
    MPI::COMM_WORLD.Allreduce( &local_numSites, &numSites, 1, MPI::UNSIGNED, MPI::SUM );
      
    if( COMM_RANK == 0 ){
      cerr << " (" << COMM_RANK << ") Number of sites in triangulation after periodic border has been removed is " << numSites << endl;
      simpleXlog << "  Final triangulation after periodic border has been removed contains " << numSites << " sites " << endl;
    }
   
  
  
  
}

/****  Assign id's to sites  ****/
void SimpleX::assign_site_ids(){

  //if periodic boundaries are applied, give sites the correct vertex index again
  //this makes sure that the sites are in the correct order
  // if(periodic){
  //   for(SITE_ITERATOR it=sites.begin(); it!=sites.end(); it++){
  //     if(it->get_border()){
  //       unsigned int real_vert_id = it->get_border() - 1;
  //       unsigned int temp_vert_id = it->get_vertex_id();
  //       it->set_vertex_id(real_vert_id);
  //       it->set_border(temp_vert_id);
  //     }
  //   }
  // }

  //set current site id
  unsigned long long int i=0;
  for( SITE_ITERATOR it=sites.begin(); it!=sites.end(); it++, i++ ){
    it->set_site_id(i);
  }

  //sort the sites on vertex id
  sort( sites.begin(), sites.end(), compare_vertex_id_site );

  //change back
  // if(periodic){
  //   for(SITE_ITERATOR it=sites.begin(); it!=sites.end(); it++){
  //     if(it->get_border()){
  //       unsigned int real_vert_id = it->get_vertex_id();
  //       unsigned int temp_vert_id = it->get_border();
  //       it->set_vertex_id(temp_vert_id);
  //       it->set_border(real_vert_id+1);
  //     }
  //   }
  // }


  //vector to keep track of the indices
  vector<unsigned long long int> indices( sites.size(), 0 );

  //store the new place of the sites in the indices array at the place of their 
  //former index
  for( unsigned long long int i=0; i<sites.size(); i++ ){
    indices[ sites[i].get_site_id() ] = i;
  }   

  //loop over all simplices
  for( SIMPL_ITERATOR it=simplices.begin(); it!=simplices.end(); it++ ){

    //former place of the sites in sites array
    unsigned long long int id1 = it->get_id1();
    unsigned long long int id2 = it->get_id2();
    unsigned long long int id3 = it->get_id3();
    unsigned long long int id4 = it->get_id4();

    //set the id to the new id of the sites
    it->set_id1( (unsigned long long int) indices[id1] );
    it->set_id2( (unsigned long long int) indices[id2] );
    it->set_id3( (unsigned long long int) indices[id3] );
    it->set_id4( (unsigned long long int) indices[id4] );

  }

  //in case of periodic boundaries update site id of corresponding vertex
  //do this only for sites on same proc
  // if(periodic){
  //   for( SITE_ITERATOR it=sites.begin(); it!=sites.end(); it++ ){
  //     if(it->get_border() ){
  // 	unsigned int site_id = it->get_border() - 1;
  // 	if(it->get_process() == sites[ indices[site_id] ].get_process()){
  // 	  it->set_border( indices[site_id] + 1 );
  // 	}
  //     }
  //   }
  // }

  //remove the indices vector
  indices.clear();

  //set correct site id
  i=0;
  for( SITE_ITERATOR it=sites.begin(); it!=sites.end(); it++, i++ ){
    it->set_site_id(i);
  }

}


/****  Initiate which site use ballistic transport  ****/
//Determine which sites are to use ballistic transport 
//and which direction conserving
void SimpleX::initiate_ballistic_sites(){


  if( ballisticTransport ){
     //if only ballistic transport is done, set ballistic to one at every site
    for( SITE_ITERATOR it=sites.begin();it!=sites.end();it++ ){
      it->set_ballistic( 1 );
    }
  }else if( dirConsTransport ){
     //direction conserving transport
    for( SITE_ITERATOR it=sites.begin();it!=sites.end();it++ ){
      it->set_ballistic( 0 );
    }
  }else if ( combinedTransport ){
     //in case of combined transport, all sites start ballistic, except for
     //the source
    if(fillChoice == AUTOMATIC){
       //loop over all sites
      for(SITE_ITERATOR it=sites.begin(); it!=sites.end(); it++){
   //sources are put first in vertex list so have
   //lowest id's
        if( it->get_vertex_id() < source_points ){
     //source is always direction conserving
          it->set_ballistic(0);
        }else{
     //non-source points start out ballistic
          it->set_ballistic(1);
        }//if source
      }//for all sites

    }else{
       //in case of properties read from file, fluxes are in
       //flux_list
       //loop over all sites
      for(SITE_ITERATOR it=sites.begin(); it!=sites.end(); it++){
   //non-source points start out ballistic
        it->set_ballistic(1);

   //if flux is bigger than 0.0, we have a source

   //exclude boundary points
        if( it->get_vertex_id() < temp_flux_list.size() ){
          if( temp_flux_list[ it->get_vertex_id() ] > 0.0){
       //source is always direction conserving
            it->set_ballistic(0);
          } //if flux
        }//if valid id
      }//for all sites
    }//if automatic filling
  }//if combined transport

}


/****  Shuffle the list of sites  ****/
//this is necessary since the sites list was created from the list of
//simplices and is thus ordered according to the simplices they belong to
void SimpleX::shuffle_sites(){

  //loop over all sites
  for( unsigned int i=0; i<sites.size(); i++ ){

    //draw a random number and convert it to a integer in the range 0-sites.size() 
    unsigned int swapIndex = (unsigned int) floor( sites.size() * gsl_rng_uniform( ran ) );

    //index should be smaller than sites.size()
    if( swapIndex == sites.size() ){
      swapIndex--;
    }

    //temporarily store the site itself and the site at the random drawn index
    Site tempSite1 = sites[ i ];
    Site tempSite2 = sites[ swapIndex ];

    //swap this site with the random drawn site
    sites[ i ] = tempSite2;
    sites[ swapIndex ] = tempSite1;

  }

  //now the sites are shuffled, so the simplices need to know the new place
  //of their sites in the sites vector

  //vector to keep track of the indices
  vector<unsigned int> indices( sites.size(), 0 );

  //store the new place of the sites in the indices array at the place of their 
  //former index
  for( unsigned int i=0; i<sites.size(); i++ ){
    indices[ sites[i].get_site_id() ] = i;
  }   

  //loop over all simplices
  for( SIMPL_ITERATOR it=simplices.begin(); it!=simplices.end(); it++ ){

    //former place of the sites in sites array
    unsigned int id1 = it->get_id1();
    unsigned int id2 = it->get_id2();
    unsigned int id3 = it->get_id3();
    unsigned int id4 = it->get_id4();

    //set the id to the new id of the sites
    it->set_id1( (unsigned int) indices[id1] );
    it->set_id2( (unsigned int) indices[id2] );
    it->set_id3( (unsigned int) indices[id3] );
    it->set_id4( (unsigned int) indices[id4] );

  }

  //remove the indices vector
  indices.clear();
  vector<unsigned int>().swap(indices);

  //now that the simplices are updates, we can assign the new 
  //indices to the sites themselves
  assign_site_ids();  


  if( COMM_RANK == 0 ){
    simpleXlog << "  Sites shuffled " << endl;
  }

}




/****************************************************************************************/
/*                           Site Properties Functions                                  */
/****************************************************************************************/


/****  Compute the properties of the sites  ****/
void SimpleX::compute_site_properties(){

  //keep track of cpu time
  double t0 = MPI::Wtime();

  //compute the site id's of the neighbours and the (numerical) volume of a site
  compute_neighbours();
  compute_volumes();

  //calculate the mean Delaunay length of every site
  calculate_line_lengths();

  //periodic sites that are shared over procs are needed for communication
  if(periodic){
    remove_periodic_simplices();
  }

  //create a list of sites existing on this proc but with a virtual copy on another proc
  fill_send_list();
  //send the properties of sites with a virtual copy to other procs
  send_site_properties();

  //in the case of ballistic transport, calculate the d (=dimension) 
  //most straight forward neighbours and the isotropic directions 
  //in which every outgoing line lies
  calculate_straight();

  //match the local neighbour id's of ballistic sites
  match_neighbours();
  
  //compute the solid angles of the direction conserving sites
  //the argument 1 means that the sphere will be rotated
  if(!ballisticTransport){
    compute_solid_angles(1);
  }

  //check for cyclic connections in the grid
  if( cyclic_check ){
    check_cyclic_connections();
  }
     
  //send the properties of neighbours needed for ballistic transport to other procs
  send_neighbour_properties();

  // if(periodic){
  //   remove_periodic_sites_2();
  // }

  //create the intensity arrays, in the case of diffuse transport just one value, in the case of 
  //ballistic transport the number of HEALPIX directions
  for( SITE_ITERATOR it=sites.begin();it!=sites.end();it++ ){

    //in case of ballistic transport, intensity has size of number of neighbours;
    //in case of direction conserving transport, intensity has 
    //the size of the tesselation of the unit sphere
    numPixels = ( it->get_ballistic() ) ? it->get_numNeigh() : number_of_directions;

    //create the intensity arrays
    it->create_intensityIn( numPixels );
    it->create_intensityOut( numPixels );

  }

  //write the time it took to triangulate the points to log file
  double t1 = MPI::Wtime();

  if( COMM_RANK == 0 ){
    simpleXlog << endl << "  Calculating sites properties took " << t1-t0 << " seconds" << endl << endl;
  }

}


/****  Compute the neighbours of every site ****/
// Generate the neighbour vectors for each vertex from the simplex vector
// Care should be taken so that no vertex is put twice in the neighbour list
void SimpleX::compute_neighbours(){

  //2D vector to hold for every site all neighbour id's
  vector< unsigned int >* neighVec;
  neighVec = new vector< unsigned int >[sites.size()];
  //iterator to loop through this vector
  vector< unsigned int >::iterator ItrT;

  //fill neighbour vector by using the list of simplices
  for( vector< Simpl >::iterator it=simplices.begin(); it!=simplices.end(); it++ ){

    //check if first vertex of this simplex is in list of the first vertex of the simplex, 
    //and add it if not
    bool found1 = 0;
    bool found2 = 0;
    bool found3 = 0;

    //loop over already existing neighbours and check whether site is already there
    for( ItrT = neighVec[ it->get_id1() ].begin(); ItrT!=neighVec[ it->get_id1() ].end(); ItrT++){
      //check if the second vertex of this simplex is already in neighbour list
      if( *ItrT == it->get_id2() ){
        found1 = 1;
      }
      //check if the third vertex of this simplex is already in neighbour list
      if( *ItrT == it->get_id3() ){
        found2 = 1;
      }
      //check if the third vertex of this simplex is already in neighbour list
      if( *ItrT == it->get_id4() ){
        found3 = 1; 
      }
    }

    //if one of the other vertices in this simplex was not yet in the neighbour 
    //list of this vertex, add it
    if(!found1){
      neighVec[ it->get_id1() ].push_back( it->get_id2() );
    }
    if(!found2){
      neighVec[ it->get_id1() ].push_back( it->get_id3() );
    }
    if(!found3){
      neighVec[ it->get_id1() ].push_back( it->get_id4() );
    }

    //check if vertex is in list of the second vertex of the simplex, and add it if not
    found1 = 0;
    found2 = 0;
    found3 = 0;
    for(ItrT = neighVec[ it->get_id2() ].begin(); ItrT!=neighVec[ it->get_id2() ].end(); ItrT++){

      if( *ItrT == it->get_id1() ){
        found1 = 1;
      }
      if( *ItrT == it->get_id3() ){
        found2 = 1;
      }
      if( *ItrT == it->get_id4() ){
        found3 = 1;
      }

    }
    //if not found, put it in the neighbour list
    if(!found1){
      neighVec[ it->get_id2() ].push_back( it->get_id1() );
    }
    if(!found2){
      neighVec[ it->get_id2() ].push_back( it->get_id3() );
    }
    if(!found3){
      neighVec[ it->get_id2() ].push_back( it->get_id4() ); 
    }


    //check if the vertex is in list of the third vertex of the simplex, and add it if not
    found1 = 0;
    found2 = 0;
    found3 = 0;
    for( ItrT = neighVec[ it->get_id3() ].begin(); ItrT!=neighVec[ it->get_id3() ].end(); ItrT++){

      if( *ItrT == it->get_id1() ){
        found1 = 1;
      }
      if( *ItrT == it->get_id2() ){
        found2 = 1;
      }
      if( *ItrT == it->get_id4() ){
        found3 = 1;
      }

    }

    //if not found, put it in the neighbour list
    if(!found1){
      neighVec[ it->get_id3() ].push_back( it->get_id1() );
    }
    if(!found2){
      neighVec[ it->get_id3() ].push_back( it->get_id2() );
    }
    if(!found3){
      neighVec[ it->get_id3() ].push_back( it->get_id4() ); 
    }


    //check if vertex is in list of the fourth vertex of the simplex, and add it if not
    found1 = 0;
    found2 = 0;
    found3 = 0;
    for(ItrT = neighVec[ it->get_id4() ].begin(); ItrT!=neighVec[ it->get_id4() ].end(); ItrT++){

      if( *ItrT == it->get_id1() ){
        found1=1;
      }
      if( *ItrT == it->get_id2() ){
        found2=1;
      }
      if( *ItrT == it->get_id3() ){
        found3=1;
      }

    }

    //if not found, put it in the neighbour list
    if(!found1){
      neighVec[ it->get_id4() ].push_back( it->get_id1() );
    }
    if(!found2){
      neighVec[ it->get_id4() ].push_back( it->get_id2() );
    }
    if(!found3){
      neighVec[ it->get_id4() ].push_back( it->get_id3() ); 
    }

  }//for all simplices

  //create neighbour array in Sites
  unsigned int i=0;
  for( vector< Site >::iterator it=sites.begin(); it!=sites.end(); it++, i++){
    it->set_numNeigh( neighVec[i].size() );
    it->create_neighId();
    unsigned int j=0;
    for( ItrT = neighVec[i].begin(); ItrT!=neighVec[i].end(); ItrT++,j++){
      it->set_neighId( j, *ItrT );
    }
    neighVec[i].clear();
  }

  if(neighVec){
    delete [] neighVec;
    neighVec = NULL;
  }

  if( COMM_RANK == 0 ){
    simpleXlog << "  Neighbour array created " << endl;
  }

}

/****  Compute site volumes  ****/
// Compute the volume of the voronoi cell around each vertex from the volumes 
// of the simplices that the vertex is part of
void SimpleX::compute_volumes(){

  //farction of the simplex volume that belongs to one site
  const double volFrac = ONE_FOURTH;

  //loop over all simplices and add volume to each site belonging to it
  for( SIMPL_ITERATOR it=simplices.begin(); it!=simplices.end(); it++){

    //volume to be added is 0.25*simplex volume
    double volSimpl = (double) it->get_volume();
    volSimpl *= volFrac;

    //add volume to first site
    double volSite1 = (double) sites[ it->get_id1() ].get_volume();
    double vol1 = volSite1+volSimpl;
    sites[ it->get_id1() ].set_volume( (float) vol1 );

    //add volume to second site
    double volSite2 = (double) sites[ it->get_id2() ].get_volume();
    double vol2 = volSite2+volSimpl;
    sites[ it->get_id2() ].set_volume( (float) vol2 );

    //add volume to third site
    double volSite3 = (double) sites[ it->get_id3() ].get_volume();
    double vol3 = volSite3+volSimpl;
    sites[ it->get_id3() ].set_volume( (float) vol3 );

    //add volume to fourth site
    double volSite4 = (double) sites[ it->get_id4() ].get_volume();
    double vol4 = volSite4+volSimpl;
    sites[ it->get_id4() ].set_volume( (float) vol4 );

  }//for all simplices

  double volume=0.0;
  for(unsigned int i=0; i<sites.size(); i++){
    if( sites[i].get_process() == COMM_RANK && !sites[i].get_border() )
      volume+=sites[i].get_volume();
  }

  MPI::COMM_WORLD.Reduce(&volume,&totalVolume,1,MPI::DOUBLE,MPI::SUM,0);

  if( COMM_RANK == 0 ){
    simpleXlog << "  Volumes computed " << endl;
  }  

}

void SimpleX::remove_periodic_simplices(){
    
  //create simplex list with only real simplices for plotting
  vector< Simpl > temp_simplices;
  Simpl tempSimpl;
  temp_simplices = simplices;
  simplices.clear();
  for( SIMPL_ITERATOR it = temp_simplices.begin(); it!=temp_simplices.end(); it++ ) {
    
    //calculate radius of circumcircle
    unsigned int simpl[4];
    simpl[0] = it->get_id1();
    simpl[1] = it->get_id2();
    simpl[2] = it->get_id3();
    simpl[3] = it->get_id4();

    //coordinates of tetrahedron, to be used to calculate circumsphere
    double xTet[4][3];

    //store the coordinates of the 4 vertices that make 
    //up this simplex in xTet
    for(short int p=0; p<dimension+1; p++) {
      xTet[p][0] = (double) sites[ simpl[p]  ].get_x();
      xTet[p][1] = (double) sites[ simpl[p]  ].get_y();
      xTet[p][2] = (double) sites[ simpl[p]  ].get_z();
    }

    // //coordinates of the centre of the circumsphere
    // double xCC, yCC, zCC;
    // 
    // //calculate the centers of the circumsphere
    // CalcCircumCenter(xTet, xCC, yCC, zCC);
    // 
    // //radius of circle
    // double r_CC = sqrt( pow(xCC - xTet[0][0],2) + pow(yCC - xTet[0][1],2) + pow(zCC - xTet[0][2],2) );
    
    //hilbert resolution is the number of cells in one direction
    unsigned int hilbert_resolution = pow( 2, hilbert_order );

    //width of the subbox is 1/number of cells in one dimension
    double subbox_width = 1.0/hilbert_resolution;
    
    vector<double> dist(6,0.0);
    dist[0] = sqrt( pow(xTet[0][0] - xTet[1][0],2) + pow(xTet[0][1] - xTet[1][1],2) + pow(xTet[0][2] - xTet[1][2],2) );
    dist[1] = sqrt( pow(xTet[0][0] - xTet[2][0],2) + pow(xTet[0][1] - xTet[2][1],2) + pow(xTet[0][2] - xTet[2][2],2) );
    dist[2] = sqrt( pow(xTet[0][0] - xTet[3][0],2) + pow(xTet[0][1] - xTet[3][1],2) + pow(xTet[0][2] - xTet[3][2],2) );
    dist[3] = sqrt( pow(xTet[1][0] - xTet[2][0],2) + pow(xTet[1][1] - xTet[2][1],2) + pow(xTet[1][2] - xTet[2][2],2) );
    dist[4] = sqrt( pow(xTet[1][0] - xTet[3][0],2) + pow(xTet[1][1] - xTet[3][1],2) + pow(xTet[1][2] - xTet[3][2],2) );
    dist[5] = sqrt( pow(xTet[2][0] - xTet[3][0],2) + pow(xTet[2][1] - xTet[3][1],2) + pow(xTet[2][2] - xTet[3][2],2) );
    
    double max=0.0;
    for(unsigned int k=0;k<6;k++){
      if(dist[k] > max){
        max = dist[k];
      }
    }
    
    if( max < 0.5*subbox_width){

      simplices.push_back( *it );
    }
  }
  temp_simplices.clear();
  
}

void SimpleX::set_periodic_sites_to_send(){

  // if(COMM_RANK == 1){
  //   for(SITE_ITERATOR it=sites.begin();it!=sites.end();it++){
  //     if(it->get_vertex_id() == 47){
  //       cerr << "  (" << COMM_RANK << ") site " << it->get_vertex_id() << ": (" << it->get_process() << ": " << it->get_x() << "," << it->get_y() << "," << it->get_z() << ") " 
  //            << it->get_border() << endl;
  //     }
  //   }
  // }

  //if multiple procs are used, use periodic sites to communicate between procs
  if(COMM_SIZE > 1){

    //store the boundaries of the domains
    vector< vector<double> > dom_bounds( dom_dec.size(), vector<double>(6,0.0) );

    //hilbert resolution is the number of cells in one direction
    unsigned int hilbert_resolution = pow( 2, hilbert_order );

    //width of the subbox is 1/number of cells in one dimension
    double subbox_width = 1.0/hilbert_resolution;

    //loop over all subboxes
    for( unsigned int i=0; i<dom_dec.size(); i++ ){

      //coordinates of the hilbert cell
      unsigned long long int coord[dimension];
      //hilbert number (place on hilbert curve)
      unsigned long long r = (unsigned long long) i;

      //find coordinates of subbox number r
      hilbert_i2c( dimension, hilbert_order, r, coord );

      //calculate minimum and maximum coordinates of subbox and boundary around it
      //minimum coordinates of subbox should not be bigger than boundary around unity domain
      double x_min = ( coord[0] == 0 ) ? 0.0  : (double) coord[0]/hilbert_resolution;
      double x_max = ( (double) coord[0]/hilbert_resolution + subbox_width  >= 1.0  ) ? 
        1.0  : (double) coord[0]/hilbert_resolution + subbox_width;

      double y_min = ( coord[1] == 0 ) ? 0.0 : (double) coord[1]/hilbert_resolution;
      double y_max = ( (double) coord[1]/hilbert_resolution + subbox_width  >= 1.0  ) ? 
        1.0  : (double) coord[1]/hilbert_resolution + subbox_width;

      double z_min = ( coord[2] == 0 ) ? 0.0  : (double) coord[2]/hilbert_resolution;
      double z_max = ( (double) coord[2]/hilbert_resolution + subbox_width  >= 1.0  ) ? 
        1.0  : (double) coord[2]/hilbert_resolution + subbox_width;

      dom_bounds[i][0] = x_min;
      dom_bounds[i][1] = x_max;
      dom_bounds[i][2] = y_min;
      dom_bounds[i][3] = y_max;
      dom_bounds[i][4] = z_min;
      dom_bounds[i][5] = z_max;

    }

    //point sites on other proc towards correct process
    unsigned int border_point = 0;
    unsigned int total_border_points = 0;
    for( SITE_ITERATOR it=sites.begin(); it!=sites.end(); it++ ){
      if( it->get_border() ){
      //if( it->get_border() && it->get_process() == COMM_RANK){
        //do reverse permutations
        for( short int i=-1; i<=1; i++ ) {
          for( short int j=-1; j<=1; j++ ) {
            for( short int k=-1; k<=1; k++ ) {
              if( (i||j||k) &&
                ( it->get_x() + 1.0*i ) >= 0.0 && ( it->get_x() + 1.0*i ) <= 1.0  &&
                ( it->get_y() + 1.0*j ) >= 0.0 && ( it->get_y() + 1.0*j ) <= 1.0  &&
                ( it->get_z() + 1.0*k ) >= 0.0 && ( it->get_z() + 1.0*k ) <= 1.0 ){

                double x = it->get_x() + 1.0*i;
                double y = it->get_y() + 1.0*j;
                double z = it->get_z() + 1.0*k;

                //determine process
                unsigned int process = COMM_SIZE;
                int found = 0;
                for(unsigned int p=0; p<dom_dec.size(); p++){
                  if( x >= dom_bounds[p][0] && x < dom_bounds[p][1] &&
                      y >= dom_bounds[p][2] && y < dom_bounds[p][3] &&
                      z >= dom_bounds[p][4] && z < dom_bounds[p][5] ){

                    process = dom_dec[p];
                    found++;
                  }
                }

                //check if found
                if(process == COMM_SIZE){
                  cerr << " (" << COMM_RANK << ") Error in assignment of periodic sites to processes!" << endl;
                  MPI::COMM_WORLD.Abort(-1);
                }

                if(found > 1){
                  cerr << " (" << COMM_RANK << ") Error in assignment of periodic sites to processes, found multiple processes!" << endl;
                  MPI::COMM_WORLD.Abort(-1);
                }

                //if(process != it->get_process() ){
                if(process != COMM_RANK){
                  it->set_vertex_id( it->get_border() - 1 );
                  it->set_border( 0 );
                  it->set_process( process );

                  border_point++;
                }

                total_border_points++;

              } //if    
            }//k
          }//j
        }//i

      }//if in border
    }//for all sites

    cerr << " (" << COMM_RANK << ") number of border points that changed process: " << border_point << " out of " << total_border_points << endl;

    // if(COMM_RANK == 0){
    //   cerr << " Domain decomposition: " << endl;
    //   for(unsigned int p=0; p<dom_dec.size(); p++){
    // 	cerr << " (" << dom_dec[p] << ") : (" << dom_bounds[p][0] << "," << dom_bounds[p][1] << ") (" 
    // 	     << dom_bounds[p][2] << "," << dom_bounds[p][3] << ") (" 
    // 	     << dom_bounds[p][4] << "," << dom_bounds[p][5] << ")" << endl;
    //   }
    // }

    dom_bounds.clear();
  }

  // if(COMM_RANK == 1){
  //   cerr << "---" << endl;
  //   for(SITE_ITERATOR it=sites.begin();it!=sites.end();it++){
  //     if(it->get_vertex_id() == 47){
  //       cerr << "  (" << COMM_RANK << ") site " << it->get_vertex_id() << ": (" << it->get_process() << ": " << it->get_x() << "," << it->get_y() << "," << it->get_z() << ") " 
  //            << it->get_border() << endl;
  //     }
  //   }
  // }

}


/****  Calculate most straight forward neighbours  ****/
// Calculate the d most straight paths for every incoming direction,
// but exclude neighbours that deviate straightAngle degrees (input
// variable) from going straight
void SimpleX::calculate_straight(){

  //keep track of straight neighbours that are not included
  unsigned int backPointing = 0;
  //keep track of the total number of straight neighbours that are included
  unsigned int outgoingCount = 0;
  //calculate the cosine of the straight angle
  double cosStraightAngle = cos( straightAngle );
  //keep track of the total sum of cosines for the correction factor
  double cosineSum = 0.0;

  //loop over all sites
  for( SITE_ITERATOR it=sites.begin(); it!=sites.end(); it++ ){
    //exclude sites on another proc and border sites
    if( it->get_process() == COMM_RANK && !it->get_border() ){

      //create the straight vector in this site
      it->create_straight();

      //loop over all neighbours to calculate the d most straightforward
      for( unsigned int j=0; j<it->get_numNeigh(); j++ ){

  // Calculate the inproduct of all Delaunay-line combo's

  //site id of this neighbour
        unsigned int neighId = it->get_neighId(j);

  //arrays to hold the neighbour vector 
  //and the maximum inner product
        double max[3], vector1[3];
  //array to hold the indices of the most straight forward neighbours
        unsigned int index[3];

  //initialise the max and index arrays
        for( short int q=0; q<dimension; q++ ){
          max[q] = -FLT_MAX; 
          index[q] = -1;
        }//for q

  //fill vector1 with neighbour vector
        vector1[0] = (double) it->get_x() - sites[ neighId ].get_x();
        vector1[1] = (double) it->get_y() - sites[ neighId ].get_y();
        vector1[2] = (double) it->get_z() - sites[ neighId ].get_z();

  //loop over all neighbours except the jth neighbour to calculate 
  //the most straight forward ones with respect to neighbour j
        for( unsigned int k=0; k<it->get_numNeigh(); k++ ){

    //site id of this neighbour
          neighId = it->get_neighId(k);

    //if the neighbour is not the neighbour we are looking at
          if( k != j ) {

      //arrays to hold the neighbour vector of this neighbour
            double vector2[3];

      //fill vector2 with neighbour vector
            vector2[0] = (double) sites[ neighId ].get_x() - it->get_x();
            vector2[1] = (double) sites[ neighId ].get_y() - it->get_y();
            vector2[2] = (double) sites[ neighId ].get_z() - it->get_z();

      //calculate inner product of both vectors
            double inprod = inproduct(vector1, vector2, dimension);

      //store d largest inner products in max array
      //store the according neighbour indices in index
            if( inprod > max[0] ) {
              max[2] = max[1];
              max[1] = max[0]; 
              max[0] = inprod;
              index[2] = index[1];
              index[1] = index[0];
              index[0] = k;
            } else if ( inprod > max[1]) {
              max[2] = max[1];
              max[1] = inprod;
              index[2] = index[1];
              index[1] = k;
        // Third most straigthforward only needed in 3D case
            } else if ( dimension == 3 && inprod > max[2] ) { 
              max[2] = inprod;
              index[2] = k;
            }
          }
        }//for k 


  //------ Check if one of the neighbours lies outside cone specified by straight_angle ------//
      //------ Calculate the deviation from a straight line for correction factor ------//

  // Always accept the first line
        it->add_straight( j, index[0] );
  // Second and third depending on angle
        for(int l=1; l < dimension; l++) {
    //add the neighbour if the maximum inner product
    //is larger than the cosine of the largest angle
          if( max[l] > cosStraightAngle ){

            it->add_straight( j, index[l] );

      // Count total number of outgoing lines
            outgoingCount++;

      // Sum of the cosine of the angles
            cosineSum += max[l];
          }
          else{
      //sum of all straight neighbours not included
            backPointing++;
          }
        }//for l

      }//for all neighbours
    }//if on this proc
  }//for all sites

  //the correction factor to simulate straight lines on the grid
  straight_correction_factor = cosineSum/double(outgoingCount);

  //add up number backpointing neighbours on all procs
  unsigned int total_backPointing;
  MPI::COMM_WORLD.Allreduce(&backPointing, &total_backPointing, 1, MPI::UNSIGNED, MPI::SUM);

  if( COMM_RANK == 0 ){
    simpleXlog << "  Most straight forward neighbours computed " << endl; 
    simpleXlog << "    straight_correction_factor = " << straight_correction_factor << endl;
    simpleXlog << "    number of straight neighbours that were outside cone of " 
      << straightAngle*180/M_PI << " degrees: " << total_backPointing << endl;
  }

}

/****  determine the place of this in neighbour array of neighbour site  ****/
void SimpleX::match_neighbours(){

  //keep track of neighbours that were not found
  int tellerNotFound=0;

  //loop over all ballistic sites and create outgoing vector
  for( SITE_ITERATOR it=sites.begin(); it!=sites.end(); it++){
    //only consider vertices on this proc and inside simulation domain
    if( it->get_process() == COMM_RANK && !it->get_border() && it->get_ballistic() ){

      //make sure outgoing is empty before creating it
      it->delete_outgoing();
      it->create_outgoing( it->get_numNeigh() );

      //loop over all neighbours of this site
      for(unsigned int j=0; j<it->get_numNeigh(); j++){
        //check if site has bee nfound
        bool found=0;
        //neighbour id of this site in neighbour array of neighbour
        unsigned int localID;
        //site id of this neighbour
        unsigned int neigh = it->get_neighId(j);
        //loop over all neighbours of neighbouring site
        for(unsigned int k=0; !found && k<sites[neigh].get_numNeigh(); k++) {
          //check if current site is in neighbour array
          if( sites[neigh].get_neighId(k) == it->get_site_id() ) {
            //if found, store the place in neighbour array
            found=1;
            localID=k;
          }	  
        }//for k

        if(!found) {
          tellerNotFound++;
        } else {
    //store the place of this site in neighbour 
    //array of neighbour site in outgoing
          it->set_outgoing( j, localID );
        }
      }//for all neighbours
    }
  }

  //if some neighbours not matched, issue warning
  if(tellerNotFound) {
    cerr << endl << endl
      << "      WARNING: " << tellerNotFound << " neighbours not matched!" << endl
      << endl;
  }

  if( COMM_RANK == 0 ){
    simpleXlog << "  Neighbours of ballistic sites matched " << endl;
  }

}


#ifdef HEAL_PIX

/****  compute the solid angles and the outgoing vector  ****/
//Associate every outgoing direction with the Delaunay lines of
//every site using healpix for quick referencing
void SimpleX::compute_solid_angles( const bool& rotate ){

  //directions are stored in header
  float **refVector;

  //get number of pixels from the value in the header of the unit sphere tesselation
  numPixels = number_of_directions;
  unsigned int numOrient = number_of_orientations;

  //orientation index has to the same on all procs, so only master computes it
  if( rotate && COMM_RANK == 0 ){
    //draw random number between 0 and number_of_orientations - 1
    orientation_index = (short int) floor( (numOrient-1)*gsl_rng_uniform( ran ) + 0.5 );
  }

  //communicate orientation index to all procs
  MPI::COMM_WORLD.Bcast( &orientation_index, 1, MPI::SHORT, 0 );

  //if the straight neighbours used are from the Delauany edge, associate every line with direction bins
  //if(straight_from_tess){
    // Assign memory to the refVector
  refVector = new float*[numPixels];
  for( unsigned int m=0; m<numPixels; m++ ){
    refVector[m]=new float[3];
  }

    //assign the first orientation to refVector
  for( unsigned int m=0; m<numPixels; m++ ){
    refVector[m][0] = (float) orient[orientation_index][m][0];
    refVector[m][1] = (float) orient[orientation_index][m][1];
    refVector[m][2] = (float) orient[orientation_index][m][2];
  }
  
  // Determine for every site which outgoing Delaunay line 
  //must be associated with a solid angle

  //create high res healpix sphere for quick referencing
  Healpix_Base healpix_base_ref(num_ref_pix_HP,RING);
  unsigned int numPixels_ref = healpix_base_ref.Npix();

    //store all directions of the healpix sphere in refVector_ref
  float **refVector_ref;
  refVector_ref = new float*[numPixels_ref];
  for( unsigned int m=0; m<numPixels_ref; m++ ){
    refVector_ref[m]=new float[3];
  }
  
  //associate this high res healpix sphere with the directions 
  //of the intensities using inner products

  //store the mapping between the directions vector and the 
  //healpix sphere
  unsigned int mapping_pixels[numPixels][numPixels_ref];
    //loop over all directions of the intensities
  for( unsigned int m=0; m<numPixels; m++ ){

    //storage for the maximum inner product and 
    //associated index in refVector_ref
    vector< float > max( numPixels_ref, -1.0 );
    vector< int > index( numPixels_ref, -1 );


    //fill max and index with inner products

    //loop over all pixels of the high res healpix sphere
    for( unsigned int p=0; p<numPixels_ref; p++ ){

      //use healpix function to find x-,y- and z-direction
      //of reference sphere
      vec3 pixpt = healpix_base_ref.pix2vec( p );
      refVector_ref[p][0] = pixpt.x;
      refVector_ref[p][1] = pixpt.y;
      refVector_ref[p][2] = pixpt.z;

      //take inner product with the intensity directions
      float tempInprod=inproduct(refVector[m], refVector_ref[p], 3);
      
      //store the inner products and indices
      max[p] = tempInprod ;
      index[p] = p;

    }    

      //sort the inner products and according indices 
    quickSortPerm( max, index, 0, numPixels_ref-1 );

    //loop over all pixels of healpix sphere
    for( unsigned int p=0; p<numPixels_ref; p++ ){
      //make sure the index value is in correct interval
      if( index[p] >=0 && index[p] < (int) numPixels_ref){
        //fill the mapping
        mapping_pixels[m][p] = index[p];
      }else{
        cerr << " (" << COMM_RANK << ") ERROR: in compute_solid_angles, mapping_pixels contains entry < 0 or > numPixels_ref, exiting" << endl;
        cerr << " (" << COMM_RANK << ") index[p]: " << index[p] << endl;
        MPI::COMM_WORLD.Abort( -1 );
      }
    }

    //clear the max and index vector
    max.clear();
    index.clear();

  }//for numPixels
 
  // Determine for every site which outgoing direction must be associated with a solid angle
  // using the calculated mapping

  //loop over all sites
  for( SITE_ITERATOR it=sites.begin(); it!=sites.end(); it++){
    //only consider vertices on this proc and inside simulation domain
    if( it->get_process() == COMM_RANK && !it->get_border() && !it->get_ballistic() ){

      //make sure outgoin is empty before creating it
      it->delete_outgoing();
      it->create_outgoing(numPixels);

      //fill the outgoing vector with temporary values
      for(unsigned int m = 0; m<numPixels; m++){
        it->set_outgoing( m, numPixels_ref+1);
      }

      // Create the vector that points from the site to the neighbour
      vec3 pixpt;
      int reference_sphere[ numPixels_ref ];
      //fill the vector with temporary values
      for( unsigned int p=0; p<numPixels_ref; p++ ){
        reference_sphere[p]=-1;
      }

      //loop over all neighbours of the site
      for(unsigned int j=0; j<it->get_numNeigh(); j++) { 

        //calculate the vector of this Delaunay line
        pixpt.x = sites[it->get_neighId(j)].get_x()-it->get_x();
        pixpt.y = sites[it->get_neighId(j)].get_y()-it->get_y();
        pixpt.z = sites[it->get_neighId(j)].get_z()-it->get_z();

        //calculate pixel on high res HEALPix sphere in which this vector resides
        int pixel_ref = healpix_base_ref.vec2pix( pixpt );
        //using previously calculated mapping, associate this with outgoing direction
        reference_sphere[ pixel_ref ] = j;
      }

      //loop over all intensity directions
      for( unsigned int m = 0; m<numPixels; m++ ){
        bool found = 0;
        //loop over all pixels of reference healpix sphere
        for( unsigned int p=0; !found && p<numPixels_ref; p++ ){
          //if the value in the mapping is fiducial
          if( reference_sphere[ mapping_pixels[m][p] ] >= 0 ){
            // associate the closest neighbour with this reference direction
            it->set_outgoing( m, reference_sphere[ mapping_pixels[m][p] ] );
            found = 1;
          }
        }
        //check if every direction is found
        if(!found){
          cerr << " Outgoing direction not found! " << endl;
          MPI::COMM_WORLD.Abort( -1 );

        }

      }//for all intensity directions

    }//if inside domain on this proc
  }// for sites 

  // free memory of the refVector_ref
  for(unsigned int m=0; m<numPixels_ref; m++){
    delete [] refVector_ref[m];
  }
  delete [] refVector_ref;

    // free memory of the refVector
  for(unsigned int m=0; m<numPixels; m++){
    delete [] refVector[m];
  }
  delete [] refVector;

    //}//if straight_from_tess
}


#else


/****  compute the solid angles and the outgoing vector  ****/
//Associate every outgoing direction with the Delaunay lines of
//every site using the inner product
void SimpleX::compute_solid_angles( const bool& rotate ){

  //directions are stored in header
  float **refVector;

  //get number of pixels from the value in the header of the unit sphere tesselation
  numPixels = number_of_directions;
  unsigned int numOrient = number_of_orientations;

  //orientation index has to the same on all procs, so only master computes it
  if( rotate && COMM_RANK == 0 ){
    //draw random number between 0 and number_of_orientations - 1
    orientation_index = (short int) floor( (numOrient-1)*gsl_rng_uniform( ran ) + 0.5 );
  }

  //communicate orientation index to all procs
  MPI::COMM_WORLD.Bcast( &orientation_index, 1, MPI::SHORT, 0 );

  //if the straight neighbours used are from the Delauany edge, associate every line with direction bins
  //if(straight_from_tess){

    // Assign memory to the refVector
  refVector = new float*[numPixels];
  for( unsigned int m=0; m<numPixels; m++ ){
    refVector[m]=new float[3];
  }

    //assign the first orientation to refVector
  for( unsigned int m=0; m<numPixels; m++ ){
    refVector[m][0] = (float) orient[orientation_index][m][0];
    refVector[m][1] = (float) orient[orientation_index][m][1];
    refVector[m][2] = (float) orient[orientation_index][m][2];
  }

    // Determine for every site which outgoing Delaunay line 
    //must be associated with a solid angle
  for( SITE_ITERATOR it=sites.begin(); it!=sites.end(); it++){
      //only consider sites on this proc and inside simulation domain
    if( it->get_process() == COMM_RANK && !it->get_border() && !it->get_ballistic() ){

  //make sure outgoing is empty
      it->delete_outgoing();
  //create the outgoing vector with numPixels entries
      it->create_outgoing(numPixels);

  //array to store neighbour vectors
      float **vectorNeigh;
  //array to store the length of the neighbour vectors
      float *vectorNeighLength;

  //assign memory to vectorNeigh and vectorNeighLength
      vectorNeigh = new float*[it->get_numNeigh()];
      for(unsigned int j=0; j<it->get_numNeigh(); j++){
        vectorNeigh[j]= new float[3];
      }
      vectorNeighLength = new float[it->get_numNeigh()];

  // Create the vector that points from the site to the neighbour
      for( unsigned int j=0; j<it->get_numNeigh(); j++ ){ 
        vectorNeigh[j][0] = sites[ it->get_neighId(j) ].get_x() - it->get_x();
        vectorNeigh[j][1] = sites[ it->get_neighId(j) ].get_y() - it->get_y();
        vectorNeigh[j][2] = sites[ it->get_neighId(j) ].get_z() - it->get_z();
      }

  // Find the neighVec that is closest (in angle) to the refVec
      for(unsigned int m=0; m<numPixels; m++) {
        float highestInprod=-FLT_MAX;
        int highestInprodNumber=-1;// to check later if a neighbour is found
        for( unsigned int j=0; j<it->get_numNeigh(); j++ ){
      //calculate inner product between the refVector and the current neighbour
          float tempInprod = inproduct(vectorNeigh[j], refVector[m], 3);
      //store the highest inner product and its id
          if( tempInprod > highestInprod ){
            highestInprod = tempInprod;
            highestInprodNumber = j;
          }
        }
    //check if inner products have been calculated correctly
        if( highestInprodNumber == -1 ){
          cerr << "Number of neighbours of site " << it->get_site_id() << " is: " 
            << int(it->get_numNeigh()) << endl;
          cerr << "In routine compute_solid_angles: there is no highest inproduct." << endl;
          MPI::COMM_WORLD.Abort( -1 );

        }
    // associate the closest neighbour with this reference direction
        it->set_outgoing(m, highestInprodNumber);
      }  //for all pixels

  //free memory of neighbour vectors
      for (unsigned int j=0;j<it->get_numNeigh();j++) {
        delete [] vectorNeigh[j]; 
      }
      delete [] vectorNeigh;
      delete [] vectorNeighLength;

    }
  }// for sites 

    // free memory of the refVector
  for(unsigned int m=0; m<numPixels; m++){
    delete [] refVector[m];
  }
  delete [] refVector;
    //}
}

#endif

//check whether there are cyclic connections in the grid and remove them if possible
void SimpleX::check_cyclic_connections(){


  // Check if cyclic connections exist and remove them
  unsigned int cyclic=0, noncyclic=0;
  for( SITE_ITERATOR it=sites.begin(); it!=sites.end(); it++){
    if( it->get_process() == COMM_RANK && !it->get_border() ){
      for(unsigned int m=0; m<numPixels; m++) {

  // For every MSN (most straightforward neighbour) check if it points to the original site
        unsigned int out = (unsigned) it->get_outgoing( m );
        for(short int k=0; k < it->get_numStraight( it->get_outgoing( m ) ); k++){

          unsigned int straight = it->get_straight( out, k );
          unsigned int msn = it->get_neighId( straight );

          unsigned int out2 = (unsigned) sites[msn].get_outgoing( m );
          for(short int l=0; l< sites[msn].get_numStraight(sites[msn].get_outgoing( m )); l++){

            unsigned int straight2 = sites[msn].get_straight( out2, l );
            unsigned int receiving = sites[msn].get_neighId( straight2 );

            if( receiving == it->get_site_id() ){ 
              cyclic++;

        // Remove the connection
              sites[msn].remove_straight( sites[msn].get_outgoing( m ), l );

        // Check if this direction has at least one outgoing edge left
              bool cycFlag = 1;
              if( !sites[msn].get_numStraight(sites[msn].get_outgoing( m )) ){

    // If not: Send radiation to a neighbouring site (other than i)
                unsigned int j=0;
                while( ( j<sites[msn].get_numNeigh() ) && cycFlag ){

                  if(j != it->get_site_id() ){

        // id of this site
                    unsigned int neighId = sites[msn].get_neighId( j );

        // Check if the connection is non-cyclic
                    cycFlag = 0;
                    unsigned int out3 = (unsigned) sites[neighId].get_outgoing( m );
                    for(short int k=0; k< sites[neighId].get_numStraight( sites[neighId].get_outgoing( m ) ); k++){

                      unsigned int straight3 = sites[neighId].get_straight( out3, l );
                      unsigned int receiving3 = sites[neighId].get_neighId( straight3 );

                      if(receiving3==msn){ 

      // Cyclic connection
                        cycFlag = 1;
                      }
                    }
                  }
                  j++;
                }

    // Add as connection
                sites[msn].add_straight( sites[msn].get_outgoing( m ), j );
              }

        // Check if a new connection has been found
              if( cycFlag ){ cerr << "WARNING: No noncyclic connection could be established for direction " << m <<  " of site " << msn << "." << endl;}

            }
            else{
              noncyclic++;
            }
          }
        }
      }
    }
  }

  if( COMM_RANK == 0 ){
    simpleXlog << "  (" << COMM_RANK << ") Removed " << cyclic << " cyclic connections, which is " << 100*double(cyclic)/double(noncyclic) << "% of total." << endl;
  }

}


void SimpleX::calculate_line_lengths(){

  vector<double> neigh_dist;
  //  double total_min_source_dist = 0.0;
  //unsigned int count = 0;
  for( SITE_ITERATOR it=sites.begin(); it!=sites.end(); it++ ){

    double totalLength=0.0;
    //double min_source_dist=FLT_MAX;
    //unsigned int cnt = 0;
    for( unsigned int j=0; j<it->get_numNeigh(); j++ ){    

      double length=0.0;
      double tmp = 0.0;

      double x1 = (double) it->get_x();
      double x2 = (double) sites[ it->get_neighId(j) ].get_x();
      tmp = x1-x2;
      if(tmp < 0.8){
        length += pow( x1 - x2, 2);
      }else{
        length += pow( x1 - x2 - 1.0, 2);
      }
      //length += pow( x1 - x2, 2);
      
      double y1 = (double) it->get_y();
      double y2 = (double) sites[ it->get_neighId(j) ].get_y();
      length += pow( y1 - y2, 2);

      double z1 = (double) it->get_z();
      double z2 = (double) sites[ it->get_neighId(j) ].get_z();
      length += pow( z1 - z2, 2);

      totalLength+=sqrt(length);

      // double source_dist=sqrt(length);
      // if(source_dist < min_source_dist && !sites[it->get_neighId(j)].get_border() ){
      // 	min_source_dist=source_dist;
      // 	cnt++;
      // }

    }

    it->set_neigh_dist( totalLength/it->get_numNeigh() );
    //it->set_neigh_dist( 1.0/32. );

    //to calculate a mean smallest neighbour distance
    if( !it->get_border() && it->get_process() == COMM_RANK ){// &&  cnt > 0 ){
      neigh_dist.push_back(it->get_neigh_dist());
    //   total_min_source_dist += min_source_dist;
    //   count++;
    }


  }//for all vertices


  //sort the distance
  sort( neigh_dist.begin(), neigh_dist.end());

  //Send the first N_mean vertices to master
  unsigned int N_mean = 0.005*origNumSites;

  vector<double> minDist_proc(N_mean*COMM_SIZE, 0);
  for(unsigned int i=0; i<N_mean; i++){
    minDist_proc[COMM_RANK*N_mean + i] = neigh_dist[i];
  }

  vector<double> minDist_total(COMM_SIZE*N_mean, 0.0);

  MPI::COMM_WORLD.Reduce(&minDist_proc[0], &minDist_total[0], N_mean*COMM_SIZE, MPI::DOUBLE, MPI::SUM, 0);

  if(COMM_RANK == 0){

    sort(minDist_total.begin(), minDist_total.end() );

    double minDist_mean = 0.0;
    for(unsigned int i=0; i<N_mean; i++){
      minDist_mean+=minDist_total[i];
    }

    minDist_mean /= (double) N_mean;

    //not yet known here...
     // UNIT_L = sizeBox * parsecToCm; 
     // cerr << " Mean minimum neighbour distance averaged over " << N_mean << " vertices is: " << minDist_mean << endl
     // 	  << " Corresponding maximum resolution: " << (int) ceil( 1.0/minDist_mean ) << endl
     // 	  << " Minimum physical scale resolved : " << minDist_mean*UNIT_L/parsecToCm << " pc" << endl
     // 	  << " Averge minimum distance between particles: " << total_min_source_dist/(double) count << endl;

    maxRes = (int) ceil( 1.0/minDist_mean );
  }

  if( COMM_RANK == 0 ){
    simpleXlog << "  Delaunay line lengths computed " << endl;
  }

}


//remove periodic sites
void SimpleX::remove_periodic_sites_2(){


  // if(COMM_RANK == 1){
  //   for( SITE_ITERATOR it=sites.begin(); it!=sites.end(); it++ ){
  //     if(it->get_vertex_id() == 47){
  //       cerr << " (" << COMM_RANK << ") site " << it->get_vertex_id() << " (" << it->get_process() 
  //            << ":" << it->get_x() << "," << it->get_y() << "," << it->get_z() << ") " << it->get_border() << endl;
  //     }
  //   }
  // }

  //find all the correct site ids
  for( SITE_ITERATOR it=sites.begin(); it!=sites.end(); it++ ){
    if(it->get_border() && it->get_process() == COMM_RANK){
    //if(it->get_border() ){
      for( unsigned int i=0; i<sites.size(); i++ ){
        if( (it->get_border()-1) == sites[i].get_vertex_id() ){
          it->set_border( sites[i].get_site_id() + 1 );
        }
      }
    }
  }

  //update sites and all neighbours on this proc
  for( SITE_ITERATOR it=sites.begin(); it!=sites.end(); it++ ){
    for(unsigned int j=0; j<it->get_numNeigh(); j++){
      if(sites[ it->get_neighId(j)].get_border() && sites[ it->get_neighId(j)].get_process() == COMM_RANK ){
      //if(sites[ it->get_neighId(j)].get_border() ){
        it->set_neighId( j, sites[ it->get_neighId(j)].get_border() - 1 );
      }
    }    
  }


  //-- remove remaining border sites --//

  //this only works if the vertex ids of the border sites are larger than those of the other vertices
  //otherwise neighbours need to matched again

  //vector to keep track of the indices
  vector<unsigned long long int> indices( sites.size(), 0 );

  //remove border sites
  vector<Site> temp_sites = sites;
  sites.clear();
  for( SITE_ITERATOR it=temp_sites.begin(); it!=temp_sites.end(); it++ ){
    if(!it->get_border()){
      sites.push_back( *it );
    }
  }

  temp_sites.clear();

  //store the new place of the sites in the indices array at the place of their 
  //former index
  for( unsigned long long int i=0; i<sites.size(); i++ ){
    indices[ sites[i].get_site_id() ] = i;
  }   

  //correct the indices of the simplices
  for( SIMPL_ITERATOR it=simplices.begin(); it!=simplices.end(); it++ ){

    //former place of the sites in sites array
    unsigned long long int id1 = it->get_id1();
    unsigned long long int id2 = it->get_id2();
    unsigned long long int id3 = it->get_id3();
    unsigned long long int id4 = it->get_id4();

    //set the id to the new id of the sites
    it->set_id1( (unsigned long long int) indices[id1] );
    it->set_id2( (unsigned long long int) indices[id2] );
    it->set_id3( (unsigned long long int) indices[id3] );
    it->set_id4( (unsigned long long int) indices[id4] );

  }
  //remove the indices vector
  indices.clear();

  //set correct site id
  unsigned int i=0;
  for( SITE_ITERATOR it=sites.begin(); it!=sites.end(); it++, i++ ){
    it->set_site_id(i);
  }

  //It could be that not all border sites are taken into account, so 
  //the numSites might have changed
  unsigned int local_numSites = 0;
  for( SITE_ITERATOR it=sites.begin(); it!=sites.end(); it++){
    if( it->get_process() == COMM_RANK ){
      local_numSites++;
    }
  }

  MPI::COMM_WORLD.Allreduce( &local_numSites, &numSites, 1, MPI::UNSIGNED, MPI::SUM );

  if( COMM_RANK == 0 ){
    cerr << " (" << COMM_RANK << ") Number of sites in triangulation after periodic border has been removed is " << numSites << endl;
    simpleXlog << "  Final triangulation after periodic border has been removed contains " << numSites << " sites " << endl;
  }

}


//create list of only local simplices
void SimpleX::remove_border_simplices(){

  //first create simplex list with only local simplices
  vector< Simpl > temp_simplices;
  Simpl tempSimpl;
  temp_simplices = simplices;
  simplices.clear();
  for( SIMPL_ITERATOR it = temp_simplices.begin(); it!=temp_simplices.end(); it++ ) {

    unsigned int simpl[4];
    simpl[0] = it->get_id1();
    simpl[1] = it->get_id2();
    simpl[2] = it->get_id3();
    simpl[3] = it->get_id4();


    unsigned int one   = it->get_id1();
    unsigned int two   = it->get_id2();
    unsigned int three = it->get_id3();
    unsigned int four  = it->get_id4();

    // ---------- check if simplex is on this proc  ---------------//


    //coordinates of tetrahedron, to be used to calculate circumsphere
    double xTet[4][3];

    //store the coordinates of the 4 vertices that make 
    //up this simplex in xTet
    for(short int p=0; p<dimension+1; p++) {
      xTet[p][0] = (double) sites[ simpl[p]  ].get_x();
      xTet[p][1] = (double) sites[ simpl[p]  ].get_y();
      xTet[p][2] = (double) sites[ simpl[p]  ].get_z();
    }

    //coordinates of the centre of the circumsphere
    double xCC, yCC, zCC;

    //calculate the centers of the circumsphere
    CalcCircumCenter(xTet, xCC, yCC, zCC);
    //determine subbox the simplex is in, by calculating the hilbert number 
    //of the subbox the circumcentre of this simplex is in
    unsigned int hilbert_resolution = pow( 2, hilbert_order );

    //first find the hilbert coordinates of the circumcentre of the simplex
    int x_hilbert = (int) floor(xCC*hilbert_resolution);
    if( x_hilbert < 0 ){
      x_hilbert = 0;
    }else if( x_hilbert >= (int) hilbert_resolution ){
      x_hilbert = hilbert_resolution - 1;
    }
    int y_hilbert = (int) floor(yCC*hilbert_resolution); 
    if( y_hilbert < 0 ){
      y_hilbert = 0;
    }else if( y_hilbert >= (int) hilbert_resolution ){
      y_hilbert = hilbert_resolution - 1;
    }
    int z_hilbert = (int) floor(zCC*hilbert_resolution);
    if( z_hilbert < 0 ){
      z_hilbert = 0;
    }else if( z_hilbert >= (int) hilbert_resolution ){
      z_hilbert = hilbert_resolution - 1;
    }

    unsigned long long int coord_CC[dimension];
    coord_CC[0] = (unsigned long long) x_hilbert;
    coord_CC[1] = (unsigned long long) y_hilbert;
    coord_CC[2] = (unsigned long long) z_hilbert;

    //determine the hilbert number of the subbox the circumcentre is in
    unsigned long long r_CC = hilbert_c2i( dimension, hilbert_order, coord_CC );

    //if this simplex belongs to this proc, we can add it
    if( dom_dec[ r_CC ] == COMM_RANK ){

      // if( !sites[ one ].get_border() && !sites[ two ].get_border() &&
      // 	  !sites[ three ].get_border() && !sites[ four ].get_border() ){

      if( !sites[ one ].get_border() || !sites[ two ].get_border() ||
      !sites[ three ].get_border() || !sites[ four ].get_border() ){

        tempSimpl.set_id1( sites[ one ].get_vertex_id() );
        tempSimpl.set_id2( sites[ two ].get_vertex_id() );
        tempSimpl.set_id3( sites[ three ].get_vertex_id() );
        tempSimpl.set_id4( sites[ four ].get_vertex_id() );
        tempSimpl.set_volume( it->get_volume() );

        simplices.push_back( tempSimpl );

      }
    }//if on this proc
  }//for all cells


  temp_simplices.clear();
  vector< Simpl >().swap(temp_simplices);

//   unsigned int numSimplLocal = simplices.size();
//   unsigned int numSimpl = 0;

//   MPI::COMM_WORLD.Allreduce(&numSimplLocal, &numSimpl, 1, MPI::UNSIGNED, MPI::SUM);

  if( COMM_RANK == 0 ){
    simpleXlog << "  Created vector with only local simplices " << endl;
  }

}


/****************************************************************************************/
/*                         Grid Update Routines                                         */
/****************************************************************************************/

/****  Rotate the solid angles of the unit sphere  ****/
void SimpleX::rotate_solid_angles(){


  //Rotation only needed in case of direction conserving transport
  if( dirConsTransport || combinedTransport ){

    //store the intensities, because the outgoing vector changes

    //make sure that photons have left all ghost vertices 
    //before storing intensities
    send_intensities();
    //make sure that the vectors that will be filled are empty 
    site_intensities.clear();
    intens_ids.clear();
    //now store the intensities
    store_intensities();

    //remember the orientation index with which the intensities were stored
    orientation_index_old = orientation_index;

    //for combined transport, calculate whether there are ballistic sites
    //that should be changed to direction conserving sites
    if(combinedTransport){
      check_ballistic_sites();
    }

    // Associate the outgoing directions with the new positions
    // Also rotate the unit sphere directions
    compute_solid_angles(1);

    //return the stored intensities to the sites
    //use the orientation index with which the intensities were stored
    return_intensities();

    //empty vectors that are no longer needed
    site_intensities.clear();
    intens_ids.clear();

    //to completely empty site_intensities, swap it with an empty vector
    vector< float >().swap( site_intensities );
    vector< unsigned  long long int >().swap( intens_ids );

  }//if rotation needed

}


//calculate normal distributed random numbers using the Box-Muller transformation
void SimpleX::randGauss( float& y1, float& y2 ){

  double w = 66.6;
  double x1, x2;

  while( w >= 1.0){
    x1 = 2.0*gsl_rng_uniform( ran ) - 1.0;
    x2 = 2.0*gsl_rng_uniform( ran ) - 1.0;
    w = x1*x1 + x2*x2;
  }

  w = sqrt( (-2.0 * log( w ) ) / w );
  y1 = x1*w;
  y2 = x2*w;

}


/**** store the intensities in big array ****/
//big array can be send by MPI
void SimpleX::store_intensities(){

  //loop over all sites
  for( SITE_ITERATOR it=sites.begin(); it!=sites.end(); it++ ){
    //only include sites in simulation domain and on this proc 
    if( it->get_process() == COMM_RANK && !it->get_border() && !it->get_ballistic() ){

      numPixels = number_of_directions;

      //store at the place in site_intensities vector
      //the vertex id
      intens_ids.push_back( it->get_vertex_id() );
      //resize the site_intensities vector
      site_intensities.insert( site_intensities.end(), numPixels, 0.0 );

      //loop over directions
      for( unsigned int j=0; j<numPixels; j++ ){
  //position of this direction in array
        unsigned int pos = j + site_intensities.size() - numPixels;

  //add all intensity to site_intensities in correct place
        site_intensities[pos] += it->get_intensityIn(j) + it->get_intensityOut(j);

  //now that they are stored, set them to zero
        it->set_intensityOut(j,0.0);
        it->set_intensityIn(j,0.0);

      }//for all pixels
    }//if
  }//for all sites

}

const vector<unsigned long long int> SimpleX::get_ballistic_sites_to_store(){

  vector<unsigned long long int> sites_to_store;

  //loop over all sites
  for( SITE_ITERATOR it=sites.begin(); it!=sites.end(); it++){
    //only consider vertices on this proc and inside simulation domain
    if( it->get_process() == COMM_RANK && !it->get_border() && it->get_ballistic() ){
      //loop over all neighbours to see if they have intensities
      bool use = 0;
      for( unsigned int j=0; !use && j<it->get_numNeigh(); j++ ){
        if( it->get_intensityIn(j) > 0.0 || it->get_intensityOut(j) > 0.0 ){
          use = 1;
        }//if there is intensity in this site
      }//for all neighbours
      if(use){
        sites_to_store.push_back( it->get_site_id() );
      }
    }//if
  }//for all sites

  return sites_to_store;

}

#ifdef HEAL_PIX
/****  Store the ballistic intensities  ****/
//Ballistic intensities are interpolated to the fixed unit sphere
//tesselation, the same that the direction conserving transport uses
//Use faster healpix referencing
void SimpleX::store_ballistic_intensities( const vector<unsigned long long int>& sites_to_store ){

  //directions are stored in header
  float **refVector;

  //get number of pixels from the value in the header of the unit sphere tesselation
  numPixels = number_of_directions;

  // Assign memory to the refVector
  refVector = new float*[numPixels];
  for( unsigned int m=0; m<numPixels; m++ ){
    refVector[m]=new float[3];
  }

  //assign the first orientation to refVector
  for( unsigned int m=0; m<numPixels; m++ ){
    refVector[m][0] = (float) orient[orientation_index][m][0];
    refVector[m][1] = (float) orient[orientation_index][m][1];
    refVector[m][2] = (float) orient[orientation_index][m][2];
  }

  //Determine the mapping between unit sphere tesselation and heal
  //pix sphere

  //create high res healpix sphere for quick referencing
  Healpix_Base healpix_base_ref(num_ref_pix_HP,RING);
  unsigned int numPixels_ref = healpix_base_ref.Npix();

  //store all directions of the healpix sphere in refVector_ref
  float **refVector_ref;
  refVector_ref = new float*[numPixels_ref];
  for( unsigned int m=0; m<numPixels_ref; m++ ){
    refVector_ref[m]=new float[3];
  }

  //associate this high res healpix sphere with the directions 
  //of the intensities using inner products

  //store the mapping between the directions vector and the 
  //healpix sphere
  unsigned int mapping_pixels[numPixels][numPixels_ref];
  //loop over all directions of the intensities
  for( unsigned int m=0; m<numPixels; m++ ){

    //storage for the maximum inner product and 
    //associated index in refVector_ref
    vector< float > max( numPixels_ref, -1.0 );
    vector< int > index( numPixels_ref, -1 );


    //fill max and index with inner products

    //loop over all pixels of the high res healpix sphere
    for( unsigned int p=0; p<numPixels_ref; p++ ){

      //use healpix function to find x-,y- and z-direction
      //of reference sphere
      vec3 pixpt = healpix_base_ref.pix2vec( p );
      refVector_ref[p][0] = pixpt.x;
      refVector_ref[p][1] = pixpt.y;
      refVector_ref[p][2] = pixpt.z;

      //take inner product with the intensity directions
      float tempInprod=inproduct(refVector[m], refVector_ref[p], 3);
      //store the inner products and indices
      max[p] = tempInprod ;
      index[p] = p;

    }    

    //sort the inner products and according indices 
    quickSortPerm( max, index, 0, numPixels_ref-1 );

    //loop over all pixels of healpix sphere
    for( unsigned int p=0; p<numPixels_ref; p++ ){
      //make sure the index value is in correct interval
      if( index[p] >=0 && index[p] < (int) numPixels_ref){
  //fill the mapping
        mapping_pixels[m][p] = index[p];
      }else{
        cerr << " (" << COMM_RANK << ") ERROR: in store_ballisitc_intensities, mapping_pixels contains entry < 0 or > numPixels_ref, exiting" << endl;
        cerr << " (" << COMM_RANK << ") index[p]: " << index[p] << endl;
        MPI::COMM_WORLD.Abort( -1 );

      }
    }

    //clear the max and index vector
    max.clear();
    index.clear();

  }//for numPixels

  // Determine for every site which outgoing direction must be associated with a solid angle
  // using the calculated mapping

  //loop over all sites to store
  for( unsigned long long int i=0; i<sites_to_store.size(); i++ ){

    unsigned long long int site_id = sites_to_store[i];

    // Create the vector that points from the site to the neighbour
    vec3 pixpt;
    int reference_sphere[ numPixels_ref ];
    //fill the vector with temporary values
    for( unsigned int p=0; p<numPixels_ref; p++ ){
      reference_sphere[p]=-1;
    }

    //loop over all neighbours of the site
    for(unsigned int j=0; j<sites[ site_id ].get_numNeigh(); j++) { 

      //calculate the vector of this Delaunay line
      pixpt.x = sites[sites[ site_id ].get_neighId(j)].get_x()-sites[ site_id ].get_x();
      pixpt.y = sites[sites[ site_id ].get_neighId(j)].get_y()-sites[ site_id ].get_y();
      pixpt.z = sites[sites[ site_id ].get_neighId(j)].get_z()-sites[ site_id ].get_z();

      //calculate pixel on high res HEALPix sphere in which this vector resides
      int pixel_ref = healpix_base_ref.vec2pix( pixpt );
      //using previously calculated mapping, associate this with outgoing direction
      reference_sphere[ pixel_ref ] = j;
    }

    //loop over all intensity directions
    vector<unsigned int> directions(numPixels, sites[ site_id ].get_numNeigh()+1);;
    for( unsigned int m = 0; m<numPixels; m++ ){
      bool found = 0;
      //loop over all pixels of reference healpix sphere
      for( unsigned int p=0; !found && p<numPixels_ref; p++ ){
  //if the value in the mapping is fiducial
        if( reference_sphere[ mapping_pixels[m][p] ] >= 0 ){
    // associate the closest neighbour with this reference direction
          directions[m] = reference_sphere[ mapping_pixels[m][p] ];
          found = 1;
        }
      }
      //check if every direction is found
      if(!found){
        cerr << " Outgoing direction not found! " << endl;
        MPI::COMM_WORLD.Abort( -1 );
      }

    }//for all intensity directions

    //store at the place in sites vector the place 
    //in site_intensities vector
    intens_ids.push_back( sites[ site_id ].get_vertex_id() );
    //resize the site_intensities vector
    site_intensities.insert( site_intensities.end(), numPixels, 0.0 );

    //loop over all neighbours
    for( unsigned int j=0; j<sites[ site_id ].get_numNeigh(); j++ ){

      //see how many directions are associated with this neighbour
      unsigned int count = 0;
      for(unsigned int m=0; m<numPixels; m++) {
        if( directions[m] == j ){
          count++;
        }
      }
      //check if this neighbour belongs to pixels in unit sphere tesselation
      if( count > 0 ){
  //give intensities to big array
        for(unsigned int m=0; m<numPixels; m++) {
          if( directions[m] == j ){
            unsigned int pos = m + site_intensities.size() - numPixels;

            double inten = ( (double) sites[ site_id ].get_intensityIn(j) + (double) sites[ site_id ].get_intensityOut(j) )/count;
            site_intensities[pos] += (float) inten;

          }
        }
  //now that they are stored, set them to zero
        sites[ site_id ].set_intensityOut(j,0.0);
        sites[ site_id ].set_intensityIn(j,0.0);
      }else{
  //if not, calculate the closest refVector entry

        float vectorNeigh[3];
        vectorNeigh[0] = sites[ sites[ site_id ].get_neighId(j) ].get_x() - sites[ site_id ].get_x();
        vectorNeigh[1] = sites[ sites[ site_id ].get_neighId(j) ].get_y() - sites[ site_id ].get_y();
        vectorNeigh[2] = sites[ sites[ site_id ].get_neighId(j) ].get_z() - sites[ site_id ].get_z();

        float highestInprod=-FLT_MAX;
        int highestInprodNumber=-1;
        for(unsigned int m=0; m<numPixels; m++) {
          float tempInprod = inproduct(vectorNeigh, refVector[m], 3);
          if( tempInprod > highestInprod ){
            highestInprod = tempInprod;
            highestInprodNumber = m;
          }
        }
        unsigned long long int pos = highestInprodNumber + site_intensities.size() - numPixels;

        double inten = ( (double) sites[ site_id ].get_intensityIn(j) + (double) sites[ site_id ].get_intensityOut(j) );
        site_intensities[pos] += (float) inten;

  //now that they are stored, set them to zero
        sites[ site_id ].set_intensityOut(j,0.0);
        sites[ site_id ].set_intensityIn(j,0.0);

      }//if count

    }//for all neighbours
  }// for sites to use

  // free memory of the refVector_ref
  for(unsigned int m=0; m<numPixels_ref; m++){
    delete [] refVector_ref[m];
  }
  delete [] refVector_ref;

  // free memory of the refVector
  for(unsigned int m=0; m<numPixels; m++){
    delete [] refVector[m];
  }
  delete [] refVector;

}

#else
/****  Store the ballistic intensities  ****/
//Ballistic intensities are interpolated to the fixed unit sphere
//tesselation, the same that the direction conserving transport uses
//use inner products for referencing
//This routine may seem to do redundant calculations, but the reason for this 
//is that we want to fill every pixel in the unit sphere tesselation
void SimpleX::store_ballistic_intensities( const vector<unsigned long long int>& sites_to_store ){


   //directions are stored in header
  float **refVector;

  //get number of pixels from the value in the header of the unit sphere tesselation
  numPixels = number_of_directions;

  // Assign memory to the refVector
  refVector = new float*[numPixels];
  for( unsigned int m=0; m<numPixels; m++ ){
    refVector[m]=new float[3];
  }

  //assign the first orientation to refVector
  for( unsigned int m=0; m<numPixels; m++ ){
    refVector[m][0] = (float) orient[orientation_index][m][0];
    refVector[m][1] = (float) orient[orientation_index][m][1];
    refVector[m][2] = (float) orient[orientation_index][m][2];
  }

  //loop over all sites to store
  for( unsigned long long int i=0; i<sites_to_store.size(); i++ ){

    //place in sites vector of this site
    unsigned long long int site_id = sites_to_store[i];


    //first associate every direction with a Delaunay line

    //array to store neighbour vectors
    float **vectorNeigh;

    //assign memory to vectorNeigh
    vectorNeigh = new float*[ sites[ site_id ].get_numNeigh()];
    for(unsigned int j=0; j<sites[ site_id ].get_numNeigh(); j++){
      vectorNeigh[j]= new float[3];
    }

    // Create the vector that points from the site to the neighbour
    for( unsigned int j=0; j<sites[ site_id ].get_numNeigh(); j++ ){ 
      vectorNeigh[j][0] = sites[ sites[ site_id ].get_neighId(j) ].get_x() - sites[ site_id ].get_x();
      vectorNeigh[j][1] = sites[ sites[ site_id ].get_neighId(j) ].get_y() - sites[ site_id ].get_y();
      vectorNeigh[j][2] = sites[ sites[ site_id ].get_neighId(j) ].get_z() - sites[ site_id ].get_z();
    }

    // Find the neighVec that is closest (in angle) to the refVec
    vector<unsigned int> directions(numPixels, sites[ site_id ].get_numNeigh()+1);
    for(unsigned int m=0; m<numPixels; m++) {
      float highestInprod=-FLT_MAX;
      int highestInprodNumber=-1;// to check later if a neighbour is found
      for( unsigned int j=0; j<sites[ site_id ].get_numNeigh(); j++ ){
  //calculate inner product between the refVector and the current neighbour
        float tempInprod = inproduct(vectorNeigh[j], refVector[m], 3);
  //store the highest inner product and its id
        if( tempInprod > highestInprod ){
          highestInprod = tempInprod;
          highestInprodNumber = j;
        }
      }
      //check if inner products have been calculated correctly
      if( highestInprodNumber == -1 ){
        cerr << "Number of neighbours of site " << sites[ site_id ].get_site_id() << " is: " 
          << int(sites[ site_id ].get_numNeigh()) << endl;
        cerr << "In routine compute_solid_angles: there is no highest inproduct." << endl;
        MPI::COMM_WORLD.Abort( -1 );
      }

      // associate the closest neighbour with this reference direction
      directions[m] = highestInprodNumber;

    }  //for all pixels

    //store at the place in sites vector the place 
    //in site_intensities vector
    intens_ids.push_back( sites[ site_id ].get_vertex_id() );
    //resize the site_intensities vector
    site_intensities.insert( site_intensities.end(), numPixels, 0.0 );

    //loop over all neighbours
    for( unsigned int j=0; j<sites[ site_id ].get_numNeigh(); j++ ){

      //see how many directions are associated with this neighbour
      unsigned int count = 0;
      for(unsigned int m=0; m<numPixels; m++) {
        if( directions[m] == j ){
          count++;
        }
      }

      //check if this neighbour is associated with any refVector entries
      if( count > 0 ){
  //give intensities to big array
        for(unsigned int m=0; m<numPixels; m++) {
          if( directions[m] == j ){
            unsigned long long int pos = m + site_intensities.size() - numPixels;

            double inten = ( (double) sites[ site_id ].get_intensityIn(j) + (double) sites[ site_id ].get_intensityOut(j) )/(double) count;
            site_intensities[pos] += (float) inten;

          }
        }
  //now that they are stored, set them to zero
        sites[ site_id ].set_intensityOut(j,0.0);
        sites[ site_id ].set_intensityIn(j,0.0);

      }else{
  //if not, calculate the closest refVector entry
        float highestInprod=-FLT_MAX;
        int highestInprodNumber=-1;
        for(unsigned int m=0; m<numPixels; m++) {
          float tempInprod = inproduct(vectorNeigh[j], refVector[m], 3);
          if( tempInprod > highestInprod ){
            highestInprod = tempInprod;
            highestInprodNumber = m;
          }
        }
        unsigned long long int pos = highestInprodNumber + site_intensities.size() - numPixels;

        double inten = ( (double) sites[ site_id ].get_intensityIn(j) + (double) sites[ site_id ].get_intensityOut(j) );
        site_intensities[pos] += (float) inten;

  //now that they are stored, set them to zero
        sites[ site_id ].set_intensityOut(j,0.0);
        sites[ site_id ].set_intensityIn(j,0.0);

      }

    }//for all neighbours

    //free memory of neighbour vectors
    for (unsigned int j=0;j<sites[ site_id ].get_numNeigh();j++) {
      delete [] vectorNeigh[j]; 
    }
    delete [] vectorNeigh;

  }//for all sites

  // free memory of the refVector
  for(unsigned int m=0; m<numPixels; m++){
    delete [] refVector[m];
  }
  delete [] refVector;

}

#endif


/**** return the stored intensities to sites ****/
void SimpleX::return_intensities(){

  //set correct number of pixels
  numPixels = number_of_directions;

  //return the intensities in chunks
  unsigned int num_chunks = 1;
  if( vertex_id_max > max_msg_to_send ){
    num_chunks = (unsigned int) ceil( (double) vertex_id_max/max_msg_to_send );
  }

  //since sites are ordered according to vertex id, we can keep the iterator
  //out of the loop, so we don't have to iterate over all sites
  SITE_ITERATOR it = sites.begin();
  for( unsigned int q=0; q<num_chunks; q++ ){

    //size of the chunk to send
    unsigned int this_chunk_size = max_msg_to_send;
    //start id of current chunk
    unsigned int start_vertex_id = q*max_msg_to_send; 
    //make sure the final chunk has correct size
    if( start_vertex_id + this_chunk_size >= vertex_id_max ){
      this_chunk_size = vertex_id_max - start_vertex_id;
    }

    //mapping from vertex id to site id
    vector< unsigned long long int > mapping_sites(this_chunk_size,vertex_id_max+1);
    bool stop = 0;
    while( it != sites.end() && !stop ){ 
      //only include direction conserving sites on this proc not in border
      if( it->get_process() == COMM_RANK && !it->get_border() && !it->get_ballistic() ){  
  //associate every vertex id with its place in the array 
        if( ( (long long int) it->get_vertex_id() - (long long int) start_vertex_id) >= 0 && 
        (it->get_vertex_id() - start_vertex_id) < this_chunk_size ){
          mapping_sites[ it->get_vertex_id() - start_vertex_id ] = it->get_site_id();
        }
      }           
      it++;
      if(it->get_vertex_id() >= (start_vertex_id + this_chunk_size) ){
        stop = 1;
      }
    }  

    //loop over site_intensities to return the intensities to sites in this chunk
    for( unsigned int i=0; i<intens_ids.size(); i++ ){
      //only include vertex ids in current chunk
      if(intens_ids[i] >= start_vertex_id && intens_ids[i] < (start_vertex_id+this_chunk_size) ){
  //get site id from vertex id using pre-calculated mapping
        unsigned int site_id = mapping_sites[ intens_ids[i] - start_vertex_id ];
  //only take into account if site id is valid, otherwise the
  //intensities should not be given back at this point
        if( site_id < sites.size() ){

    //loop over all directions
    //be careful, j is the position of the direction in the orientation of the PREVIOUS run!!!
          for( unsigned int j=0; j<numPixels; j++ ){

      //get the intensity that belongs to this site
            double inten = (double) site_intensities[ numPixels*i + j ];

      //get the closest associate with this neighbour from the mapping 
      //stored in the header file
            unsigned int pos = maps[orientation_index][orientation_index_old][j];

      //in case there is already intensity in this direction
      //obsolete with the new header file
            inten += (double) sites[site_id].get_intensityOut(pos);

      //assign the intensity to the site
            sites[site_id].set_intensityOut( pos, (float) inten );
            sites[site_id].set_intensityIn( pos, 0.0 );

          }//for all directions
        }//if valid id
      }
    } //for all intens_ids

    mapping_sites.clear();

  }

}

#ifdef HEAL_PIX

/****  Return the ballistic intensities to sites  ****/
//Interpolate from unit sphere tesselation to Delaunay lines
//using heal_pix referencing
void SimpleX::return_ballistic_intensities(){

  //create reference vector with the direction of the intensities 
  //that were stored

  //directions of stored intensities are stored in header
  float **refVector;

  //get number of pixels from the value in the header of the unit sphere tesselation
  numPixels = number_of_directions;

  // Assign memory to the refVector
  refVector = new float*[numPixels];
  for( unsigned int m=0; m<numPixels; m++ ){
    refVector[m]=new float[3];
  }

  //assign the orientation with which the intensities were stored to refVector
  for( unsigned int m=0; m<numPixels; m++ ){
    refVector[m][0] = (float) orient[orientation_index_old][m][0];
    refVector[m][1] = (float) orient[orientation_index_old][m][1];
    refVector[m][2] = (float) orient[orientation_index_old][m][2];
  }

  //Determine the mapping between unit sphere tesselation and heal
  //pix sphere

  //create high res healpix sphere for quick referencing
  Healpix_Base healpix_base_ref(num_ref_pix_HP,RING);
  unsigned int numPixels_ref = healpix_base_ref.Npix();

  //store all directions of the healpix sphere in refVector_ref
  float **refVector_ref;
  refVector_ref = new float*[numPixels_ref];
  for( unsigned int m=0; m<numPixels_ref; m++ ){
    refVector_ref[m]=new float[3];
  }

  //associate this high res healpix sphere with the directions 
  //of the intensities using inner products

  //store the mapping between the directions vector and the 
  //healpix sphere
  unsigned int mapping_pixels[numPixels][numPixels_ref];
  //loop over all directions of the intensities
  for( unsigned int m=0; m<numPixels; m++ ){

    //storage for the maximum inner product and 
    //associated index in refVector_ref
    vector< float > max( numPixels_ref, -1.0 );
    vector< int > index( numPixels_ref, -1 );


    //fill max and index with inner products

    //loop over all pixels of the high res healpix sphere
    for( unsigned int p=0; p<numPixels_ref; p++ ){

      //use healpix function to find x-,y- and z-direction
      //of reference sphere
      vec3 pixpt = healpix_base_ref.pix2vec( p );
      refVector_ref[p][0] = pixpt.x;
      refVector_ref[p][1] = pixpt.y;
      refVector_ref[p][2] = pixpt.z;

      //take inner product with the intensity directions
      float tempInprod=inproduct(refVector[m], refVector_ref[p], 3);
      //store the inner products and indices
      max[p] = tempInprod ;
      index[p] = p;

    }    

    //sort the inner products and according indices 
    quickSortPerm( max, index, 0, numPixels_ref-1 );

    //loop over all pixels of healpix sphere
    for( unsigned int p=0; p<numPixels_ref; p++ ){
      //make sure the index value is in correct interval
      if( index[p] >=0 && index[p] < (int) numPixels_ref){
  //fill the mapping
        mapping_pixels[m][p] = index[p];
      }else{
        cerr << " (" << COMM_RANK << ") ERROR: in return_ballisitc_intensities, mapping_pixels contains entry < 0 or > numPixels_ref, exiting" << endl;
        cerr << " (" << COMM_RANK << ") index[p]: " << index[p] << endl;
        MPI::COMM_WORLD.Abort( -1 );
      }
    }

    //clear the max and index vector
    max.clear();
    index.clear();

  }//for numPixels

  //return the intensities in chunks
  unsigned int num_chunks = 1;
  if( vertex_id_max > max_msg_to_send ){
    num_chunks = (unsigned int) ceil( (double) vertex_id_max/max_msg_to_send );
  }

  //since sites are ordered according to vertex id, we can keep the iterator
  //out of the loop, so we don't have to iterate over all sites
  SITE_ITERATOR it = sites.begin();
  for( unsigned int q=0; q<num_chunks; q++ ){

    //size of the chunk to send
    unsigned int this_chunk_size = max_msg_to_send;
    //start id of current chunk
    unsigned long long int start_vertex_id = q*max_msg_to_send; 
    //make sure the final chunk has correct size
    if( start_vertex_id + this_chunk_size >= vertex_id_max ){
      this_chunk_size = vertex_id_max - start_vertex_id;
    }

   //mapping from vertex id to site id
    vector< unsigned long long int > mapping_sites(this_chunk_size,vertex_id_max+1);
    bool stop = 0;
    while( it != sites.end() && !stop ){ 
      //only include direction conserving sites on this proc not in border                                                                                      
      if( it->get_process() == COMM_RANK && !it->get_border() && it->get_ballistic() ){  
  //associate every vertex id with its place in the array                                                                                                 
        if( ( (long long int) it->get_vertex_id() - (long long int) start_vertex_id ) >= 0 && 
        (it->get_vertex_id() - start_vertex_id) < this_chunk_size ){
          mapping_sites[ it->get_vertex_id() - start_vertex_id ] = it->get_site_id();
        }
      }           
      it++;
      if(it->get_vertex_id() >= (start_vertex_id + this_chunk_size) ){
        stop = 1;
      }
    }  

    //loop over site_intensities to return the intensities to sites in this chunk
    for( unsigned long long int i=0; i<intens_ids.size(); i++ ){
      //only include vertex ids in current chunk
      if(intens_ids[i] >= start_vertex_id && intens_ids[i] < (start_vertex_id+this_chunk_size) ){
  //get site id from vertex id using pre-calculated mapping
        unsigned long long int site_id = mapping_sites[ intens_ids[i] - start_vertex_id ];
  //only take into account if site id is valid, otherwise the
  //intensities should not be given back at this point
        if( site_id < sites.size() ){

    // Create the vector that points from the site to the neighbour
          vec3 pixpt;
          int reference_sphere[ numPixels_ref ];
    //fill the vector with temporary values
          for( unsigned int p=0; p<numPixels_ref; p++ ){
            reference_sphere[p]=-1;
          }

    //loop over all neighbours of the site
          for( unsigned int j=0; j<sites[site_id].get_numNeigh(); j++ ){ 

      //calculate the vector of this Delaunay line
            pixpt.x = sites[ sites[site_id].get_neighId(j) ].get_x() - sites[site_id].get_x();
            pixpt.y = sites[ sites[site_id].get_neighId(j) ].get_y() - sites[site_id].get_y();
            pixpt.z = sites[ sites[site_id].get_neighId(j) ].get_z() - sites[site_id].get_z();

      //calculate pixel on high res HEALPix sphere in which this vector resides
            int pixel_ref = healpix_base_ref.vec2pix( pixpt );
      //using previously calculated mapping, associate this with outgoing direction
            reference_sphere[ pixel_ref ] = j;
          }

    //loop over all intensity directions
          for( unsigned int m = 0; m<numPixels; m++ ){
            bool found = 0;
      //loop over all pixels of reference healpix sphere
            for( unsigned int p=0; !found && p<numPixels_ref; p++ ){
        //if the value in the mapping is fiducial
              if( reference_sphere[ mapping_pixels[m][p] ] >= 0 ){


    // associate the closest neighbour with this reference direction
                unsigned int j = reference_sphere[ mapping_pixels[m][p] ];
                found = 1;

    //give this neighbour the intensity in site_intensities
                double inten = site_intensities[ i*numPixels + m ];
                inten += (double) sites[site_id].get_intensityOut(j);

                sites[site_id].set_intensityOut(j,inten);
                sites[site_id].set_intensityIn(j,0.0);

              }
            }
      //check if every direction is found
            if(!found){
              cerr << " Outgoing direction not found! " << endl;
              MPI::COMM_WORLD.Abort( -1 );
            }

          }//for all intensity directions



        }//if valid id
      }//if vertex id in range
    }//for all intens_ids
    mapping_sites.clear();
  }//for all chunks

  // free memory of the refVector_ref
  for(unsigned int m=0; m<numPixels_ref; m++){
    delete [] refVector_ref[m];
  }
  delete [] refVector_ref;

  // free memory of the refVector
  for(unsigned int m=0; m<numPixels; m++){
    delete [] refVector[m];
  }
  delete [] refVector;

}

#else

/****  Return the ballistic intensities to sites  ****/
//Interpolate from unit sphere tesselation to Delaunay lines
//using inner products
void SimpleX::return_ballistic_intensities(){


  //create reference vector with the direction of the intensities 
  //that were stored

  //directions of stored intensities are stored in header
  float **refVector;

  //get number of pixels from the value in the header of the unit sphere tesselation
  numPixels = number_of_directions;

  // Assign memory to the refVector
  refVector = new float*[numPixels];
  for( unsigned int m=0; m<numPixels; m++ ){
    refVector[m]=new float[3];
  }

  //assign the orientation with which the intensities were stored to refVector
  for( unsigned int m=0; m<numPixels; m++ ){
    refVector[m][0] = (float) orient[orientation_index_old][m][0];
    refVector[m][1] = (float) orient[orientation_index_old][m][1];
    refVector[m][2] = (float) orient[orientation_index_old][m][2];
  }

  //return the intensities in chunks
  unsigned int num_chunks = 1;
  if( vertex_id_max > max_msg_to_send ){
    num_chunks = (unsigned int) ceil( (double) vertex_id_max/max_msg_to_send );
  }

  //since sites are ordered according to vertex id, we can keep the iterator
  //out of the loop, so we don't have to iterate over all sites
  SITE_ITERATOR it = sites.begin();
  for( unsigned int q=0; q<num_chunks; q++ ){

    //size of the chunk to send
    unsigned int this_chunk_size = max_msg_to_send;
    //start id of current chunk
    unsigned long long int start_vertex_id = q*max_msg_to_send; 
    //make sure the final chunk has correct size
    if( start_vertex_id + this_chunk_size >= vertex_id_max ){
      this_chunk_size = vertex_id_max - start_vertex_id;
    }

   //mapping from vertex id to site id
    vector< unsigned long long int > mapping_sites(this_chunk_size,vertex_id_max+1);
    bool stop = 0;
    while( it != sites.end() && !stop ){ 
      //only include direction conserving sites on this proc not in border                                                                                      
      if( it->get_process() == COMM_RANK && !it->get_border() && it->get_ballistic() ){  
  //associate every vertex id with its place in the array                                                                                                 
        if( (it->get_vertex_id() - start_vertex_id) >= 0 && (it->get_vertex_id() - start_vertex_id) < this_chunk_size ){
          mapping_sites[ it->get_vertex_id() - start_vertex_id ] = it->get_site_id();
        }
      }           
      it++;
      if(it->get_vertex_id() >= (start_vertex_id + this_chunk_size) ){
        stop = 1;
      }
    }  

    //loop over site_intensities to return the intensities to sites in this chunk
    for( unsigned long long int i=0; i<intens_ids.size(); i++ ){
      //only include vertex ids in current chunk
      if(intens_ids[i] >= start_vertex_id && intens_ids[i] < (start_vertex_id+this_chunk_size) ){
  //get site id from vertex id using pre-calculated mapping
        unsigned long long int site_id = mapping_sites[ intens_ids[i] - start_vertex_id ];
  //only take into account if site id is valid, otherwise the
  //intensities should not be given back at this point
        if( site_id < sites.size() ){

    //array to store neighbour vectors
          float **vectorNeigh;

    //assign memory to vectorNeigh and vectorNeighLength
          vectorNeigh = new float*[sites[site_id].get_numNeigh()];
          for(unsigned int j=0; j<sites[site_id].get_numNeigh(); j++){
            vectorNeigh[j]= new float[3];
          }
    // Create the vector that points from the site to the neighbour
          for( unsigned int j=0; j<sites[site_id].get_numNeigh(); j++ ){ 
            vectorNeigh[j][0] = sites[ sites[site_id].get_neighId(j) ].get_x() - sites[site_id].get_x();
            vectorNeigh[j][1] = sites[ sites[site_id].get_neighId(j) ].get_y() - sites[site_id].get_y();
            vectorNeigh[j][2] = sites[ sites[site_id].get_neighId(j) ].get_z() - sites[site_id].get_z();
          }

    // Find the neighVec that is closest (in angle) to the refVec
          for(unsigned int m=0; m<numPixels; m++) {
            float highestInprod=-FLT_MAX;
            int highestInprodNumber=-1;// to check later if a neighbour is found
            for( unsigned int j=0; j<sites[site_id].get_numNeigh(); j++ ){
        //calculate inner product between the refVector and the current neighbour
              float tempInprod = inproduct(vectorNeigh[j], refVector[m], 3);
        //store the highest inner product and its id
              if( tempInprod > highestInprod ){
                highestInprod = tempInprod;
                highestInprodNumber = j;
              }
            }

      //check if inner products have been calculated correctly
            if( highestInprodNumber == -1 ){
              cerr << "Number of neighbours of site " << site_id << " is: " 
                << int(sites[site_id].get_numNeigh()) << endl;
              cerr << "In routine return_ballistic_intensities: there is no highest inproduct." << endl;
              MPI::COMM_WORLD.Abort( -1 );
            }

            unsigned int j = highestInprodNumber;

      //give intensity of this pixel to Delaunay line with highest inner product
            double inten = site_intensities[ i*numPixels + m ];

            inten += (double) sites[site_id].get_intensityOut(j);

            sites[site_id].set_intensityOut(j,inten);
            sites[site_id].set_intensityIn(j,0.0);

            total_inten += site_intensities[ numPixels*i + m ];

          }  //for all pixels


    //free memory of neighbour vectors
          for (unsigned int j=0;j<sites[site_id].get_numNeigh();j++) {
            delete [] vectorNeigh[j]; 
          }
          delete [] vectorNeigh;



        }//if valid id
      }//if vertex id in range
    }//for all intens_ids
    mapping_sites.clear();
  }//for all chunks

}

#endif


/****  Check which sites are ballistic and which not  ****/
//Only used in case of combined transport
void SimpleX::check_ballistic_sites(){

  //vector to hold ballistic sites that need to be stored
  vector< unsigned long long int > sites_to_store;
  //vector to hold all sites in internal boundary that have
  //been changed
  vector< unsigned long long int > sites_to_send;

  //loop over sites
  for( SITE_ITERATOR it=sites.begin(); it!=sites.end(); it++ ){
    //only include sites on this proc inside simulation domain
    if( it->get_process() == COMM_RANK && !it->get_border()  ){

      //site is considered for DCT only when ionising radiation has come through
      //or the density is 0.0
      bool is_ionised = 0;
      if( it->get_n_HII() > 0.0 || ( it->get_n_HI() + it->get_n_HII() ) == 0.0 ){
        is_ionised = 1;
      }
      //calculate the optical depth
      double tau = it->get_n_HI() * UNIT_D * it->get_neigh_dist() * UNIT_L * cross_H * straight_correction_factor;
      //if optical depth is smaller than the switch set by user, switch from
      //ballistic to direction conserving if the site has seen ionising radiation
      if( tau < switchTau && is_ionised  ){
  //check if site already was direction conserving
  //if not, change
        if( it->get_ballistic() == 1 ){

          it->set_ballistic( 0 );

    //store the id of this site in sites_to store so 
    //intensity of this site will be stored in site_intensities
          sites_to_store.push_back( it->get_site_id() );

    //if this site has incarnation on other proc, 
    //put it in the list to send
          bool flag = 0;
    //loop over all neighbours to check whether this site has at least
    //one neighbour on another proc
          for( unsigned int j=0; !flag && j<it->get_numNeigh(); j++ ) {
            unsigned long long int neigh = it->get_neighId(j);
      //check if neighbour belongs to other proc
            if( sites[neigh].get_process() != COMM_RANK ){
        //if so, set flag
              flag = 1;	
              sites_to_send.push_back( it->get_site_id() );
            }
          }//for all neighbours
        }//if switch from ballistic to dirCons

      }else{
  //check if site already was ballistic
  //if not, change
        if( it->get_ballistic() == 0 && it->get_flux() == 0.0 ){

          it->set_ballistic( 1 );

    //if this site has incarnation on other proc, 
    //put it in the list to send
          bool flag = 0;
    //loop over all neighbours to check whether this site has at least
    //one neighbour on another proc
          for( unsigned int j=0; !flag && j<it->get_numNeigh(); j++ ) {
            unsigned long long int neigh = it->get_neighId(j);
      //check if neighbour belongs to other proc
            if( sites[neigh].get_process() != COMM_RANK ){
        //if so, set flag
              flag = 1;	
              sites_to_send.push_back( it->get_site_id() );
            }
          }//for all neighbours

    //put the intensities in site_intensities
          numPixels = number_of_directions;

    //store at the place in site_intensities vector
    //the vertex id
          intens_ids.push_back( it->get_vertex_id() );
    //resize the site_intensities vector
          site_intensities.insert( site_intensities.end(), numPixels, 0.0 );
    //loop over directions
          for( unsigned int j=0; j<numPixels; j++ ){
      //position of this direction in array
            unsigned long long int pos = j + site_intensities.size() - numPixels;

      //add all intensity to site_intensities in correct place
            site_intensities[pos] += it->get_intensityIn(j) + it->get_intensityOut(j);

      //now that they are stored, set them to zero
            it->set_intensityOut(j,0.0);
            it->set_intensityIn(j,0.0);

          }//for all pixels

    //delete the intensity arrays
          it->delete_intensityOut();
          it->delete_intensityIn();

    //create new intensity arrays with correct size
          it->create_intensityIn( it->get_numNeigh() );
          it->create_intensityOut( it->get_numNeigh() );


    //match the neighbours of this site

    //make sure outgoing is empty before creating it
          it->delete_outgoing();
          it->create_outgoing( it->get_numNeigh() );

    //loop over all neighbours of this site
          for(unsigned int j=0; j<it->get_numNeigh(); j++){
      //check if site has bee nfound
            bool found=0;
      //neighbour id of this site in neighbour array of neighbour
            unsigned int localID = INT_MAX;
      //site id of this neighbour
            unsigned long long int neigh = it->get_neighId(j);
      //loop over all neighbours of neighbouring site
            for(unsigned int k=0; !found && k<sites[neigh].get_numNeigh(); k++) {
        //check if current site is in neighbour array
              if( sites[neigh].get_neighId(k) == it->get_site_id() ) {
    //if found, store the place in neighbour array
                found=1;
                localID=k;
              }	  
            }//for k

      //store the place of this site in neighbour 
      //array of neighbour site in outgoing
            if( localID < INT_MAX ){
              it->set_outgoing( j, localID );
            }else{
              cerr << " (" << COMM_RANK << ") Error in check_ballistic_sites(): Neighbour not matched " << endl;
              MPI::COMM_WORLD.Abort( -1 );
            }

          }//for all neighbours

        }//if switch from dirCons to ballistic
      }//if tau smaller than switch

    }//if
  }//for all sites

  //store the intensities of sites that were switched from ballistic
  //to direction conserving
  store_ballistic_intensities( sites_to_store );

  //give the new properties to incarnation on other proc
  send_site_ballistics( sites_to_send );

  //give the sites that were switched from ballistic
  //to direction conserving correct intensity arrays
  for(unsigned long long int i=0; i<sites_to_store.size(); i++ ){

    unsigned long long int site_id = sites_to_store[i];

    sites[ site_id ].delete_intensityOut();
    sites[ site_id ].delete_intensityIn();

    sites[ site_id ].create_intensityOut( number_of_directions );
    sites[ site_id ].create_intensityIn( number_of_directions );

  }

  sites_to_store.clear();
  sites_to_send.clear();

  //return the intensities of the sites that were switched from 
  //direction conserving to ballistic
  return_ballistic_intensities();

}



/****  Do grid dynamics  ****/
//Vertices can be removed according to the new situation in the simulation,
//or vertices can be moved around to avoid Poisson noise
bool SimpleX::grid_dynamics(){

  double t0 = MPI::Wtime();

  //bool to check whether grid has changed
  bool didUpdate = 0;

  //number of sites that was removed from simulation
  unsigned int numSites_to_remove = 0;

  //removal of points in low opacity regions
  if(updates){
    numSites_to_remove = mark_sites_for_removal();
    if( numSites_to_remove > 0 ){
      didUpdate = 1;
    }
  }

  //only recalculate grid if something has changed
  if(didUpdate){


    //make sure that photons have left all ghost vertices 
    //before calculating local intensities
    send_intensities();

     //if dynamic updates are done, sites need to be removed
    if( numSites_to_remove > 0 ){

      //compute relevant update properties and send them to master
      //total_sites_marked_for_removal vector is filled here for master
      send_sites_marked_for_removal();

      //master proc calculates which sites will be removed
      if( COMM_RANK == 0 ){
        get_sites_to_remove(numSites_to_remove);
      }

      //send the properties of sites to be removed to other procs
      //total_sites_marked_for_removal vector is filled 
      //here for other procs
      send_sites_to_remove();

      //redistribute the properties of sites that are removed
      //to their neighbours
      redistribute_site_properties();

    }

   //make sure that the vectors that will be filled are empty 
    site_intensities.clear();
    intens_ids.clear();

    //store the intensities in big array
    store_intensities();

    //in this case, also store the ballistic intensities
    vector< unsigned long long int > sites_to_store = get_ballistic_sites_to_store();

    store_ballistic_intensities( sites_to_store );
    sites_to_store.clear();

    //remember the orientation index with which the intensities were stored
    orientation_index_old = orientation_index;


    //store the relevant properties of the sites to be used 
    //in the coming run
    store_site_properties();

    //create a list of vertices to be triangulated
    create_new_vertex_list();

//     for( SITE_ITERATOR it=sites.begin(); it!=sites.end(); it++ ){
//       it->delete_intensityIn();
//     }


    //clear temporary structures in the sites 
    //and completely clear the sites vector
    clear_temporary();
    sites.clear();
    vector< Site >().swap(sites);

    simplices.clear();
    vector< Simpl >().swap(simplices);

    //send the list to master proc
    send_new_vertex_list();


    if( COMM_RANK == 0 ){
      cerr << " (" << COMM_RANK << ") Computing triangulation" << endl;
    }

    if(COMM_RANK == 0){

      //set boundary around unity domain
      if(periodic){
        create_periodic_boundary();
      }else{
        create_boundary();
      }

      //create octree
      create_vertex_tree();


      //for now, don't do a new domain decomposition
      //decompose the domain
      //decompose_domain();

      //assign process to vertices
      //this is necessary to include boundary sites
      assign_process();

    }else{
      //clear vertices on other procs
      vertices.clear();
      vector< Vertex >().swap(vertices);
    }

    //send the vertex positions round to all procs, so every proc (temporary!) has its own copy
    send_vertices();

    //for now, don't do a new domain decomposition
    //send the domain decomposition to other procs
    //send_dom_dec();

    //now that all procs have a list of vertices, create octree
    if( COMM_RANK != 0 ){
      create_vertex_tree();
    }

    //compute the triangulation
    compute_triangulation();

    //create the sites vector from the vertex list
    create_sites();

    //the list of sites was obtained from list of simplices, and therefore have an order which might lead 
    //to problems when using the dynamic update routines
    //shuffle_sites();

    //assign the correct site ids to the sites
    assign_site_ids();  

    //since no new domain decomposition is done, no need to send site_properties
    //send the properties of the sites to other procs
    //send_site_physics();


    //return the physical properties to the sites
    return_physics();

  }

  double t1 = MPI::Wtime();
  if( COMM_RANK == 0 ){
    simpleXlog << endl << "  Calculating updates took " << t1-t0 << " seconds" << endl << endl;
  }

  return didUpdate;

}


/****  Store the relevant properties of the sites  ****/
// stored are:  vertex_id
//              site_id
//              n_HI
//              n_HII  
//              flux
//              ballistic          
void SimpleX::store_site_properties(){

  //make sure site properties vector is empty
  site_properties.clear();

  //first determine if there are changes from 
  //ballistic to dirCons sites or back

  if(combinedTransport){
    //loop over sites
    for( SITE_ITERATOR it=sites.begin(); it!=sites.end(); it++ ){
      //check if site is on this proc and if tau is bigger than zero.
      //if neighbour distance<0.0, site is flagged as removed from simulation
      if( it->get_process() == COMM_RANK && !it->get_border() && it->get_neigh_dist() >= 0.0 ){

        bool is_ionised = 0;
        if( it->get_n_HII() > 0.0 || (it->get_n_HII() + it->get_n_HI()) == 0.0 ){
          is_ionised = 1;
        }

  //calculate optical depth

  //use the volume to find average neighbour distance, since the 
  //triangulation hasn't been updated yet after sites were removed
        double aver_neigh_dist = 3*pow( (double) it->get_volume(), 1.0/3.0 )/(4.0*M_PI);
        double tau = it->get_n_HI() * UNIT_D * aver_neigh_dist * UNIT_L * cross_H * straight_correction_factor;
  //if optical depth is smaller than the switch set by user, switch from
  //ballistic to direction conserving
        if( tau < switchTau && is_ionised){
          it->set_ballistic( 0 );
        }else{
          it->set_ballistic( 1 );
        }
  //sources are always direction conserving
        if( it->get_flux() > 0.0 ){
          it->set_ballistic( 0 );
        }
      }
    }
  }

  //loop over sites
  for( SITE_ITERATOR it=sites.begin(); it!=sites.end(); it++ ){
    //check if site is on this proc and if tau is bigger than zero.
    //if tau<0.0, site is flagged as removed from simulation
    if( it->get_process() == COMM_RANK && !it->get_border() && it->get_neigh_dist() >= 0.0 ){
      //create structure that stores site properties
      Site_Update temp;
      //assign the relevant properties
      temp = *it;

      //test whether it's better to send atoms instead 
      //double N_H = (double) it->get_number_density() * (double) it->get_volume();
      //temp.set_number_density( (float) N_H );

      //put the properties in vector
      site_properties.push_back( temp );
    }
  }

  if( COMM_RANK == 0 ){
    simpleXlog << " New sites array filled" << endl;
  }

}

/****  Create new list of vertices to be triangulated  ****/
void SimpleX::create_new_vertex_list(){

  //temporary structure to hold vertex information
  Vertex tempVert; 

  //make sure the vertices vector is empty
  vertices.clear();

  //loop over all sites
  for( SITE_ITERATOR it=sites.begin(); it!=sites.end(); it++ ){
    //only sites on this proc. Sites with tau<0.0 are flagged
    //to be removed from the simulation, so don't use them
    if( it->get_process() == COMM_RANK && !it->get_border() && it->get_neigh_dist() >= 0.0 ){

      //assign the relevant properties to the vertex
      tempVert.set_x( it->get_x() );
      tempVert.set_y( it->get_y() );
      tempVert.set_z( it->get_z() );
      tempVert.set_border( it->get_border() );
      tempVert.set_vertex_id( it->get_vertex_id() );
      tempVert.set_process( it->get_process() );

      //put them in the vertices array
      vertices.push_back( tempVert );
    }
  }

  cerr << " (" << COMM_RANK << ") number of vertices: " << vertices.size() << endl;

  if( COMM_RANK == 0 ){
    simpleXlog << "  Created new vertex list" << endl;
  }

}


/**** Calculate the number of sites to be removed from simulation  ****/
//Number of sites that can be removed from simulation depends on optical 
//depth of sites and the minimum resolution set by the user
unsigned int SimpleX::mark_sites_for_removal(){

  //keep track of the volume that sites marked for removal represent
  double local_volume = 0.0;

  //vector to hold teh indices of sites marked for removal must be empty
  sites_marked_for_removal.clear();

  //loop over all sites
  for( SITE_ITERATOR it=sites.begin(); it!=sites.end(); it++ ){
    //only consider sites on this proc inside simulation domain
    if( it->get_process() == COMM_RANK && !it->get_border() ){

      bool is_ionised = 0;
      if( it->get_n_HII() > 0.0 || (it->get_n_HII() + it->get_n_HI()) == 0.0 ){
        is_ionised = 1;
      }

      //if the neutral fraction of the cell has changed, 
      //the minimum resolution is not exceeded and the 
      //site is not a source, mark it for possible removal
      double tau = it->get_n_HI() * UNIT_D * (double) it->get_neigh_dist() * UNIT_L * cross_H * straight_correction_factor;
      if( is_ionised && 
        tau < switchTau && 
        it->get_neigh_dist() < 1.0/minResolution && 
      it->get_flux() == 0.0 ){

  //put the site id of the site in vector
        sites_marked_for_removal.push_back( it->get_site_id() );
  //add the volume of the site to the total volume
  //occupied by the sites marked for removal
        local_volume += it->get_volume();

      }//if update criterium is satisfied
    }//if inside simulation domain
  }//for all sites

  //add up the number of sites marked for removal of all procs
  unsigned int local_numSites_to_remove = sites_marked_for_removal.size();
  unsigned int numSites_to_remove;
  MPI::COMM_WORLD.Allreduce( &local_numSites_to_remove, &numSites_to_remove, 1, MPI::UNSIGNED, MPI::SUM );

  //add up the occupied volume of these sites of all procs
  double totalVolume;
  MPI::COMM_WORLD.Allreduce( &local_volume, &totalVolume, 1, MPI::DOUBLE, MPI::SUM );


  if( COMM_RANK == 0){
    simpleXlog << "  Number of candidates to remove: " << numSites_to_remove << endl;
  }

  //maximum number of sites to avoid numerical diffusion
  //Assume the ionised region is spherical and the source sits in the centre.
  //maximum number of steps to I-Front is minResolution
  //Number of steps relates to number of points in volume as:
  // (ksi * S)^3 = (3 * N)/(4 * M_PI)
  //with ksi = 1.237 this gives N = 7.93 * S^3
  //  double max_numSteps = minResolution;
  //unsigned int max_numSites = 8 * (unsigned int)pow( max_numSteps, 3 );

  //minimum number of sites to retain minimum resolution
  unsigned int min_numSites = (unsigned int) ceil( totalVolume*( (double) pow( minResolution, 3 ) ) );

  //calculate the number of sites that need to be removed if the minimum 
  //number of sites in the optically thin volume is min_numSites
  if( numSites_to_remove > min_numSites ){
    numSites_to_remove -= min_numSites;
  }else{
    numSites_to_remove = 0;
    sites_marked_for_removal.clear();
    vector< unsigned long long int >().swap(sites_marked_for_removal);
  }

  if( COMM_RANK == 0){
    simpleXlog << "  Total number of sites to remove: " << numSites_to_remove << endl;
  }

  return numSites_to_remove;

}


/****  Calculate the sites that are to be removed from simulation  ****/
void SimpleX::get_sites_to_remove( const unsigned int& numSites_to_remove ){

  // Create a vector of vectors of ints to be used as a histogram.
  //number of bins is user specified variable
  vector< vector<Site_Remove> > sites_marked_for_removal_hist;
  sites_marked_for_removal_hist.assign( nbins, vector<Site_Remove>() );

  //shuffle the total_sites_marked_for_removal so there is 
  //no preference for the processor that communicated last
  for( unsigned int i=0; i<total_sites_marked_for_removal.size(); i++ ){

    unsigned int swapIndex = (int)floor( total_sites_marked_for_removal.size() * gsl_rng_uniform( ran ) );
    if( swapIndex == total_sites_marked_for_removal.size() ) swapIndex--;

    Site_Remove tempSite1 = total_sites_marked_for_removal[ i ];
    Site_Remove tempSite2 = total_sites_marked_for_removal[ swapIndex ];

    total_sites_marked_for_removal[ i ] = tempSite2;
    total_sites_marked_for_removal[ swapIndex ] = tempSite1;

  }

  // Create the 'histogram' based on the min and max of the point density found above.
  // A given bin of the histogram contains integers that correspond to candidates with 
  // point densities appropriate for that bin.
  for( unsigned int i=0; i<total_sites_marked_for_removal.size(); i++ ){
    double delta_y = total_sites_marked_for_removal[i].get_update_property();
    sites_marked_for_removal_hist[ unsigned( floor( delta_y * double( nbins - 1.0 )  )) ].push_back( total_sites_marked_for_removal[i] );
  }

  //total_sites_to_update can now be cleared for later use
  total_sites_marked_for_removal.clear();
  vector<Site_Remove>().swap( total_sites_marked_for_removal );

  // Create a cumulative distribution from the cumulativeDist to be used as the
  // integral up to a given density of the probability distribution.
  vector<unsigned int> cumulative_dist( nbins, 0 );
  cumulative_dist[0] = sites_marked_for_removal_hist[0].size();

  for( unsigned int i=1; i<nbins; i++ ){
    cumulative_dist[i] = sites_marked_for_removal_hist[i].size() + cumulative_dist[i-1];
  }

  // Make an array based on the cumulative distribution containing the
  // iDs of the flagged points in the original cumulative_dist.  Note
  // that the iDs of several bins of the original cumulative_dist can
  // be grouped in one bin of the new array.

  // The array has the same number of entries as the distribution but this
  // need not be.
  vector< vector<Site_Remove> >  draw_array;
  draw_array.assign( nbins, vector<Site_Remove>() );

  // Use the min and max for mapping
  unsigned int min_dist = cumulative_dist[0];
  unsigned int max_dist = cumulative_dist[ cumulative_dist.size() - 1 ];
  unsigned int delta_dist = max_dist - min_dist;

  //fill draw_array with sites marked for removal according to cumulative distribution
  for( unsigned int i=0; i<nbins; i++ ){
    unsigned int pos = (unsigned int) floor( ( ( cumulative_dist[i] - min_dist ) / double( delta_dist) ) * double( nbins - 1 ) );
    draw_array[pos].insert( draw_array[pos].end(), sites_marked_for_removal_hist[i].begin(), sites_marked_for_removal_hist[i].end() );
  }


  // From this array, the points are drawn at random.
  for( unsigned int i=0; i<numSites_to_remove; i++ ){

    unsigned int pos = (unsigned int) floor( gsl_rng_uniform( ran ) * double( nbins - 1 ) ) ;

    // If an empty bin is struck, the next non-empty bin higher up is used.
    // This makes sure that points at high point density are removed preferentially.

    // One could think of a smarter way of doing this...
    while( draw_array[pos].empty() ){
      pos++;
      if( pos > nbins )
        break;
    }

     //total_sites_to_update contains now the sites drawn from the histogram
    if( pos < nbins ){ //make sure that if all bins above are empty, point is ignored
      total_sites_marked_for_removal.push_back( draw_array[pos].back() );
      draw_array[pos].pop_back();// This works only is the sites have no order, otherwise, another random number should be used
    }else{
      i--;
    }

  }//for all sites to remove

  //free memory
  draw_array.clear();
  cumulative_dist.clear();
  sites_marked_for_removal_hist.clear();

  if( COMM_RANK == 0 ){
    simpleXlog << "  Total number of sites selected for removal: " << total_sites_marked_for_removal.size() << endl;
  }

}

/****  Redistribute the porperties of the removed sites  ****/
//Give the relevant properties of sites that are to be removed from
//the simulation to the still existing neighbours of the site
//If sites are deleted, tau is set to zero
void SimpleX::redistribute_site_properties(){

  //determine which of the sites is on this proc
  vector<Site_Remove> sites_to_remove;
  for( unsigned long long int i=0; i<total_sites_marked_for_removal.size(); i++ ){
    if( total_sites_marked_for_removal[i].get_process() == COMM_RANK ){
      sites_to_remove.push_back( total_sites_marked_for_removal[i] );
    }
  }

  cerr << " (" << COMM_RANK << ") Total number of sites to remove: " << sites_to_remove.size() 
    << " out of " << total_sites_marked_for_removal.size() << endl;

  //clear the vector for use in next run
  total_sites_marked_for_removal.clear();
  vector<Site_Remove>().swap(total_sites_marked_for_removal);

  //make sure numPixels has the correct value
  numPixels = number_of_directions;

  //keep track of sites for which the redistribution failed
  unsigned int give_up=0;
  //now give  the properties of deleted vertices to existing neighbours on this proc
  for( unsigned long long int i=0; i<sites_to_remove.size(); i++ ){

    //site id of the site
    unsigned long long int id = sites_to_remove[i].get_site_id();

    //check if the id is valid
    if( id > sites.size() ){
      cerr << " (" << COMM_RANK << ") id bigger than sites.size(): " << id << " " << sites.size() << endl;
    } 

    //site properties are only distributed to neighbours that 
    //live on this site, to avoid communications
    unsigned int on_proc = 0;
    //first make sure there are neighbours on this proc
    for(unsigned int j=0; j<sites[ id ].get_numNeigh(); j++){

      //site id of the neighbour
      unsigned int neighId = sites[ id ].get_neighId(j);

      //check if the neighbour is on this proc, not in border and
      //not already deleted
      if( sites[ neighId ].get_process() == COMM_RANK && 
        !sites[ neighId ].get_border() &&               
      sites[ neighId ].get_neigh_dist() >= 0.0 ){

        on_proc++;

      }
    }//for all neighbours

    //if there are neighbours that are allowed to receive the information
    //redistribute the properties
    if( on_proc > 0 ){
      //loop over all neighbours
      for(unsigned int j=0; j<sites[ id ].get_numNeigh(); j++){

  //site id of the neighbour
        unsigned int neighId = sites[ id ].get_neighId(j);

  //check if the neighbour is on this proc, not in border and
  //not already deleted
        if( sites[ neighId ].get_process() == COMM_RANK && 
          !sites[ neighId ].get_border() &&               
        sites[ neighId ].get_neigh_dist() >= 0.0 ){

    //total volume of the new site
          double volume = (double) sites[neighId].get_volume() + (double) sites[id].get_volume()/on_proc;

    //number of neutral atoms in combined cell
          double N_HI = (double) sites[id].get_n_HI() * (double) sites[id].get_volume();
          N_HI /= on_proc;
          N_HI += (double) sites[neighId].get_n_HI() * (double) sites[neighId].get_volume();

    //number of ionised atoms in combined cell
          double N_HII = (double) sites[id].get_n_HII() * (double) sites[id].get_volume();
          N_HII /= on_proc;
          N_HII += (double) sites[neighId].get_n_HII() * (double) sites[neighId].get_volume();

    //number density and ionised fraction in combined cell
          double n_HI = N_HI/volume;
          double n_HII = N_HII/volume;

    //assign properties of combined cell to neighbour
          sites[neighId].set_volume( (float) volume );
          sites[neighId].set_n_HI( (float) n_HI );
          sites[neighId].set_n_HII( (float) n_HII );

          unsigned int num = sites[id].get_ballistic() ? sites[id].get_numNeigh() : number_of_directions;
          for(unsigned int q=0;q<num;q++){

            double inten = ( sites[id].get_intensityIn(q) + sites[id].get_intensityOut(q) )/on_proc; 

      //assign intensities
            if( sites[id].get_ballistic()){
              if( sites[ neighId ].get_ballistic() ){
                unsigned int neighIdLoc = (unsigned) sites[id].get_outgoing(j);
                sites[ neighId ].addRadiationDiffOut( neighIdLoc, (float) inten );
              }else{
                unsigned int neighIdLoc = (unsigned) sites[id].get_outgoing(j);

    //if the neighbour is not ballistic, find out in which direction 
    //bins the photons should go
                vector< unsigned int > dir_to_use;
    //this trick only works for neighbours on this proc,
    //otherwise outgoing array does not exist
                for( unsigned int n=0; n<number_of_directions; n++ ){
                  if( sites[ neighId ].get_outgoing(n) == neighIdLoc ){
                    dir_to_use.push_back( n );
                  }
                }//for all directions

    //if no direction is associated with this neighbour, find one
                if(dir_to_use.size() == 0){

      //find closest associated direction for this Delaunay line

      //directions are stored in header
                  float **refVector;
      // Assign memory to the refVector
                  refVector = new float*[number_of_directions];
                  for( unsigned int n=0; n<number_of_directions; n++ ){
                    refVector[n]=new float[3];
                  }
      //assign the first orientation to refVector
                  for( unsigned int n=0; n<number_of_directions; n++ ){
                    refVector[n][0] = (float) orient[orientation_index][n][0];
                    refVector[n][1] = (float) orient[orientation_index][n][1];
                    refVector[n][2] = (float) orient[orientation_index][n][2];
                  }

                  float vectorNeigh[3];
                  vectorNeigh[0] = sites[id].get_x() - sites[ neighId ].get_x();
                  vectorNeigh[1] = sites[id].get_y() - sites[ neighId ].get_y();
                  vectorNeigh[2] = sites[id].get_z() - sites[ neighId ].get_z();

      // Find the neighVec that is closest (in angle) to the refVec
                  float highestInprod=-FLT_MAX;
                  int highestInprodNumber=-1;
                  for(unsigned int n=0; n<number_of_directions; n++) {
                    float tempInprod = inproduct(vectorNeigh, refVector[n], 3);
                    if( tempInprod > highestInprod ){
                      highestInprod = tempInprod;
                      highestInprodNumber = n;
                    }
                  }
                  dir_to_use.push_back( highestInprodNumber );

      // free memory of the refVector
                  for(unsigned int n=0; n<number_of_directions; n++){
                    delete [] refVector[n];
                  }
                  delete [] refVector;

                }//if no direction associated with this linw

                for(unsigned int n=0; n<dir_to_use.size() ; n++ ){
                  sites[ neighId ].addRadiationDiffOut( dir_to_use[n], (float) inten/dir_to_use.size() );
                }//for all neighbours to use
                dir_to_use.clear();


              }//if neighbour is ballistic

            }else{ //if site itself is not ballistic

              if(sites[ neighId ].get_ballistic() ){
    //if the neighbour is ballistic, find out which 
    //Delaunay line is associated with this direction
                bool found = 0;
                for( unsigned int n=0; !found && n<sites[ neighId ].get_numNeigh(); n++ ){
                  if( sites[neighId].get_neighId(n) == sites[id].get_site_id() ){
                    found = 1;
                    sites[neighId].addRadiationDiffOut( n, (float) inten );
                  } 
                }
                if(!found){
                  cerr << " Error in redistribute_site_properties(): neighbour not found! " << endl;
                }
              }else{
    //if not ballistic,
    //send the photons to the neighbour in the same direction bin
                sites[ neighId ].addRadiationDiffOut( q , (float) inten );

              }//if neighbour is ballistic

            }//if site is ballistic
          }//for all neighbours/pixels

    //	  for(unsigned int j=0; j<numPixels; j++ ){

// 	    double intensityOut = 0.0; 

// 	    //get the intensities of this site from site_intensities
// 	    unsigned int intens_id_site = mapping_intens[ sites[id].get_vertex_id() ];
// 	    //check if there were intensities stored for this site
// 	    if( intens_id_site < intens_ids.size() ){
// 	      intensityOut = (double) site_intensities[ intens_id_site + j ]/on_proc;
// 	    }

// 	    //only add the intensity if it is bigger than 0.0
// 	    if( intensityOut > 0.0 ){

// 	      //get the position in the site_intensities array
// 	      unsigned int intens_id_neigh = mapping_intens[ sites[neighId].get_vertex_id() ];
// 	      //it is possible that this site didn't have any intensities yet and so is not stored
// 	      //in that case, create an entry
// 	      if( intens_id_neigh >= intens_ids.size() ){
// 		//keep mapping up to date
// 		mapping_intens[ sites[neighId].get_vertex_id() ] = intens_ids.size();
// 		//add to intens_ids
// 		intens_ids.push_back( sites[neighId].get_vertex_id() );
// 		//give site_intensities extra entries for the new intensities
// 		site_intensities.insert( site_intensities.end(), numPixels, 0.0 );
// 	      }

// 	      //now the entry in the mapping surely exists
// 	      intens_id_neigh = mapping_intens[ sites[neighId].get_vertex_id() ];

// 	      site_intensities[ intens_id_neigh + j ] += (float) intensityOut;

//	    }

    //}//for all pixels


        }//if
      }//for all neighbours

      //set the properties of the site to zero after they
      //have been assigned to neighbour
      //a negative neigh_dist means the site has been removed
      sites[id].set_neigh_dist( -1.0 );
      sites[id].set_volume( 0.0 );
      sites[id].set_n_HI( 0.0 );
      sites[id].set_n_HII( 0.0 );
      unsigned int num = sites[id].get_ballistic() ? sites[id].get_numNeigh() : number_of_directions;
      for(unsigned int j=0; j<num; j++){
        sites[id].set_intensityIn(j,0.0);
        sites[id].set_intensityOut(j,0.0);
      }

    }else{

       //if all neighbours are already deleted, or on another proc, 
       //try neighbours neighbours, otherwise give up

       //loop over all neighbours
      for(unsigned int j=0; j<sites[ id ].get_numNeigh(); j++){

  //site id of the neighbour
        unsigned int neighId = sites[ id ].get_neighId(j);
  //check if site id is valid
        if( neighId > sites.size() ){
          cerr << " (" << COMM_RANK << ") neighId bigger than sites.size(): " << neighId << " " << sites.size() << endl;
        }

  //loop over neighbours neighbours
        for( unsigned int k=0; k<sites[neighId].get_numNeigh(); k++ ){
    //site id of neighbour
          unsigned int neighNeighId = sites[ neighId ].get_neighId(k);

    //check if the neighbour's neighbour is not the site itself
          if( neighNeighId != id){
      //check if the neighbour's neighbour is on this proc, not in border and
      //not already deleted
            if( sites[ neighNeighId ].get_process() == COMM_RANK && 
              !sites[ neighNeighId ].get_border() &&               
              sites[ neighNeighId ].get_neigh_dist() >= 0.0 )             
            {
              on_proc++;
            }
          }
        }//for all neighbours neighbours
      } //for all neighbours

      //if there are valid neighbour's neighbours
      if(on_proc){
  //loop over all neighbours
        for(unsigned int j=0; j<sites[ id ].get_numNeigh(); j++){

    //site id of the neighbour
          unsigned int neighId = sites[ id ].get_neighId(j);
    //check if site id is valid
          if( neighId > sites.size() ){
            cerr << " (" << COMM_RANK << ") neighId bigger than sites.size(): " << neighId << " " << sites.size() << endl;
          }

    //loop over neighbours neighbours
          for( unsigned int k=0; k<sites[neighId].get_numNeigh(); k++ ){
      //site id of neighbour
            unsigned int neighNeighId = sites[ neighId ].get_neighId(k);

      //check if the neighbour's neighbour is not the site itself
            if( neighNeighId != id){
        //check if the neighbour's neighbour is on this proc, not in border and
        //not already deleted
              if( sites[ neighNeighId ].get_process() == COMM_RANK && 
                !sites[ neighNeighId ].get_border() &&               
              sites[ neighNeighId ].get_neigh_dist() >= 0.0 ){

    //total volume of the new site
                double volume = (double) sites[neighNeighId].get_volume() + (double) sites[id].get_volume()/on_proc;

    //number of neutral atoms in combined cell
                double N_HI = (double) sites[id].get_n_HI() * (double) sites[id].get_volume();
                N_HI /= on_proc;
                N_HI += (double) sites[neighNeighId].get_n_HI() * (double) sites[neighNeighId].get_volume();

    //number of ionised atoms in combined cell
                double N_HII = (double) sites[id].get_n_HII() * (double) sites[id].get_volume();
                N_HII /= on_proc;
                N_HII += (double) sites[neighNeighId].get_n_HII() * (double) sites[neighNeighId].get_volume();

    //number density and ionised fraction in combined cell
                double n_HI = N_HI/volume;
                double n_HII = N_HII/volume;

    //assign properties of combined cell to neighbour
                sites[neighNeighId].set_volume( (float) volume );
                sites[neighNeighId].set_n_HI( (float) n_HI );
                sites[neighNeighId].set_n_HII( (float) n_HII );

                unsigned int num = sites[id].get_ballistic() ? sites[id].get_numNeigh() : number_of_directions;
                for(unsigned int q=0;q<num;q++){

                  double inten = ( sites[id].get_intensityIn(q) + sites[id].get_intensityOut(q) )/on_proc; 

      //assign intensities
                  if( sites[id].get_ballistic()){
                    if( sites[ neighNeighId ].get_ballistic() ){
                      bool found = 0;
                      for( unsigned int n=0; !found && n<sites[ neighNeighId ].get_numNeigh(); n++ ){
                        if( sites[neighNeighId].get_neighId(n) == sites[id].get_site_id() ){
                          found = 1;
                          sites[ neighNeighId ].addRadiationDiffOut( n, (float) inten );
                        }
                      }
          //unsigned int neighIdLoc = (unsigned) sites[neighId].get_outgoing(k);
          //sites[ neighNeighId ].addRadiationDiffOut( neighIdLoc, (float) inten );
                    }else{
                      unsigned int neighIdLoc = (unsigned) sites[neighId].get_outgoing(k);
          //if the neighbour is not ballistic, find out in which direction 
          //bins the photons should go
                      vector< unsigned int > dir_to_use;
          //this trick only works for neighbours on this proc,
          //otherwise outgoing array does not exist
                      for( unsigned int n=0; n<number_of_directions; n++ ){
                        if( sites[ neighNeighId ].get_outgoing(n) == neighIdLoc ){
                          dir_to_use.push_back( n );
                        }
                      }//for all directions

          //if no direction is associated with this neighbour, find one
                      if(dir_to_use.size() == 0){

      //find closest associated direction for this Delaunay line

      //directions are stored in header
                        float **refVector;
      // Assign memory to the refVector
                        refVector = new float*[number_of_directions];
                        for( unsigned int n=0; n<number_of_directions; n++ ){
                          refVector[n]=new float[3];
                        }
      //assign the first orientation to refVector
                        for( unsigned int n=0; n<number_of_directions; n++ ){
                          refVector[n][0] = (float) orient[orientation_index][n][0];
                          refVector[n][1] = (float) orient[orientation_index][n][1];
                          refVector[n][2] = (float) orient[orientation_index][n][2];
                        }

                        float vectorNeigh[3];
                        vectorNeigh[0] = sites[id].get_x() - sites[ neighNeighId ].get_x();
                        vectorNeigh[1] = sites[id].get_y() - sites[ neighNeighId ].get_y();
                        vectorNeigh[2] = sites[id].get_z() - sites[ neighNeighId ].get_z();

      // Find the neighVec that is closest (in angle) to the refVec
                        float highestInprod=-FLT_MAX;
                        int highestInprodNumber=-1;
                        for(unsigned int n=0; n<number_of_directions; n++) {
                          float tempInprod = inproduct(vectorNeigh, refVector[n], 3);
                          if( tempInprod > highestInprod ){
                            highestInprod = tempInprod;
                            highestInprodNumber = n;
                          }
                        }
                        dir_to_use.push_back( highestInprodNumber );

      // free memory of the refVector
                        for(unsigned int n=0; n<number_of_directions; n++){
                          delete [] refVector[n];
                        }
                        delete [] refVector;

                      }//if no direction associated with this linw

                      for(unsigned int n=0; n<dir_to_use.size() ; n++ ){
                        sites[ neighNeighId ].addRadiationDiffOut( dir_to_use[n], (float) inten/dir_to_use.size() );
                      }//for all neighbours to use
                      dir_to_use.clear();


                    }//if neighbour is ballistic

                  }else{
                    if(sites[ neighNeighId ].get_ballistic() ){
          //if the neighbour is ballistic, find out which 
          //Delaunay line is associated with this direction
                      bool found = 0;
                      for( unsigned int n=0; !found && n<sites[ neighNeighId ].get_numNeigh(); n++ ){
                        if( sites[neighNeighId].get_neighId(n) == sites[neighId].get_site_id() ){
                          found = 1;
                          sites[neighNeighId].addRadiationDiffOut( n, (float) inten );
                        } 
                      }
                      if(!found){
                        cerr << " Error in redistribute_site_properties(): neighbour not found! " << endl;
                      }
                    }else{
          //if not ballistic,
          //send the photons to the neighbour in the same direction bin
                      sites[ neighNeighId ].addRadiationDiffOut( q , (float) inten );

                    }//if neighbour is ballistic

                  }//if site is ballistic
                }//for all neighbours/pixels



// 		//assign intensities
// 		for(unsigned int j=0; j<numPixels; j++ ){

// 		  double intensityOut = 0.0; 

// 		  //get the intensities of this site from site_intensities
// 		  unsigned int intens_id_site = mapping_intens[ sites[id].get_vertex_id() ];
// 		  //check if there were intensities stored for this site
// 		  if( intens_id_site < intens_ids.size() ){
// 		    intensityOut = (double) site_intensities[ intens_id_site + j ]/on_proc;
// 		  }

// 		  //only add the intensity if it is bigger than 0.0
// 		  if( intensityOut > 0.0 ){
// 		    //get the position in the site_intensities array
// 		    unsigned int intens_id_neigh = mapping_intens[ sites[neighNeighId].get_vertex_id() ];
// 		    //it is possible that this site didn't have any intensities yet and so is not stored
// 		    //in that case, create an entry
// 		    if( intens_id_neigh >= intens_ids.size() ){
// 		      //keep mapping up to date
// 		      mapping_intens[ sites[neighNeighId].get_vertex_id() ] = intens_ids.size();
// 		      //add to intens_ids
// 		      intens_ids.push_back( sites[neighNeighId].get_vertex_id() );
// 		      //give site_intensities extra entries for the new intensities
// 		      site_intensities.insert( site_intensities.end(), numPixels, 0.0 );
// 		    }

// 		    //now the entry in the mapping surely exists
// 		    intens_id_neigh = mapping_intens[ sites[neighNeighId].get_vertex_id() ];

// 		    site_intensities[ intens_id_neigh + j ] += (float) intensityOut;

// 		  }

//		}//for all pixels

              }//if valid site
            }//if not this site
          }//for all neighbour's neighbours
        }//for all neighbours

  //set the properties of teh site to zero after they
  //have been assigned to neighbour
  //a negative tau means the site has been removed
        sites[id].set_neigh_dist( -1.0 );
        sites[id].set_volume( 0.0 );
        sites[id].set_n_HI( 0.0 );
        sites[id].set_n_HII( 0.0 );
        unsigned int num = sites[id].get_ballistic() ? sites[id].get_numNeigh() : number_of_directions;
        for(unsigned int j=0; j<num; j++){
          sites[id].set_intensityIn(j,0.0);
          sites[id].set_intensityOut(j,0.0);
        }

      }else{
  //if not neighbour's neighbours are valid, give up
  //and do not remove site
        give_up++;
      }      


    }//if on_proc 
  }//for all sites to delete

  //free memory
  sites_to_remove.clear();

  if(give_up){
    cerr << "  (" << COMM_RANK << ") Gave up on updating " << give_up << " vertices" << endl;
  }

  if( COMM_RANK == 0 ){
    simpleXlog << "  Sites to be removed redistributed " << endl;
  }

}


/****************************************************************************************/
/*                            Send Routines                                             */
/****************************************************************************************/

/****  Send the domain decomposition to all procs  ****/
void SimpleX::send_dom_dec(){

  //number of subboxes is number of cells in hilbert curve
  unsigned int number_of_subboxes = pow( 2, hilbert_order*dimension );

  //if there's a lot of subboxes, send in chunks
  unsigned int num_chunks = 1;
  if( number_of_subboxes > max_msg_to_send ){
    num_chunks = (unsigned int) ceil( (double) number_of_subboxes/max_msg_to_send);
  }

  if(COMM_RANK != 0){

    //make sure dom_dec is empty
    dom_dec.clear();
    vector<unsigned int>().swap(dom_dec);
    //resize to the number of subboxes
    //dom_dec.resize(number_of_subboxes);

  }

  //send the domain decomposition to other procs in chunks
  for(unsigned int i=0; i<num_chunks; i++ ){

    //temporary vector to hold the vertices
    vector< unsigned int > temp;
    //size of the chunk to send
    unsigned int this_chunk_size = max_msg_to_send;
    //start id of current chunk
    unsigned int start_id = i*max_msg_to_send; 
    //make sure the final chunk has correct size
    if( start_id + this_chunk_size >= number_of_subboxes ){
      this_chunk_size = number_of_subboxes - start_id;
    }

    //master fills temp
    if(COMM_RANK == 0){
      for(unsigned int j = start_id; j<(start_id+this_chunk_size); j++){
        temp.push_back( dom_dec[j] );
      }
      //other procs resize to correct sizw  
    }else{
      temp.resize( this_chunk_size );  
    }

    //broadcast this chunk
    MPI::COMM_WORLD.Bcast(&temp[0],this_chunk_size,MPI::UNSIGNED,0);

    //procs not master fill dom_dec vector from temp vector
    if( COMM_RANK != 0 ){
      dom_dec.insert( dom_dec.end(), temp.begin(), temp.end() );
    }

    //clear temp for next use
    temp.clear();

  }

  //if not every proc has received vertices, exit
  if( dom_dec[ dom_dec.size() - 1 ] != (COMM_SIZE - 1) ){
    if( COMM_RANK == 0 ){
      cerr << " (" << COMM_RANK << ") Error, domain decomposition and number of processes don't match, exiting" << endl;
      simpleXlog << "  Error, domain decomposition and number of processes don't match, exiting" << endl;
    }
    MPI::COMM_WORLD.Abort( -1 );
  }

  if( COMM_RANK == 0 ){
    simpleXlog << "  Domain decomposition sent to all procs " << endl;
  }

}


/****  Send all vertices to all processors  ****/
void SimpleX::send_vertices(){

  //make sure every proc knows the number of sites
  MPI::COMM_WORLD.Bcast(&numSites,1,MPI::UNSIGNED,0);

  //if there's a lot of vertices, send in chunks
  unsigned int num_chunks = 1;
  if( numSites > max_msg_to_send ){
    num_chunks = (unsigned int) ceil( (double) numSites/max_msg_to_send);
  }

  //make sure vertices vector is empty on all procs except master
  //and give correct size
  if(COMM_RANK != 0){
    vertices.clear();
    vector<Vertex>().swap(vertices);
  }

  //send the vertices to other procs in chunks
  for(unsigned int i=0; i<num_chunks; i++ ){

    //size of the chunk to send
    unsigned int this_chunk_size = max_msg_to_send;
    //start id of current chunk
    unsigned int start_vertex_id = i*max_msg_to_send; 
    //make sure the final chunk has correct size
    if( start_vertex_id + this_chunk_size >= numSites ){
      this_chunk_size = numSites - start_vertex_id;
    }

    //temporary vector to hold the vertices
    vector< Vertex > temp(this_chunk_size);
    //temp.resize( this_chunk_size );    
    //Vertex.construct_datatype();

    //master fills temp
    if(COMM_RANK == 0){
      for(unsigned int j = start_vertex_id; j<(start_vertex_id+this_chunk_size); j++){
        temp[ j - start_vertex_id ] = vertices[j];
      }
    }

    //broadcast this chunk
    MPI::COMM_WORLD.Bcast(&temp[0],this_chunk_size*sizeof(Vertex),MPI::BYTE,0);
    //MPI::COMM_WORLD.Bcast(&temp[0],this_chunk_size,Vertex::MPI_Type,0);


    //procs not master fill vertices vector from temp vector
    if( COMM_RANK != 0 ){
      vertices.insert( vertices.end(), temp.begin(), temp.end() );
    }

    //clear temp for next use
    temp.clear();

  }

  if( COMM_RANK == 0 ){
    simpleXlog << "  Vertices sent to all procs " << endl;
  }

}

/****  Create list of sites that belong to other proc  ****/
void SimpleX::fill_send_list(){

  //make sure list is empty before filling it
  send_list.clear();

  //loop over all sites
  for( SITE_ITERATOR it=sites.begin(); it!=sites.end(); it++ ){
    //if site is not on this proc and not in border add to list
    if( it->get_process() != COMM_RANK && !it->get_border() )
      send_list.push_back( it->get_site_id() );    
  }

  if( COMM_RANK == 0 ){
    simpleXlog << "  Send list created. " << endl;
    simpleXlog << "  Number of duplicate sites needed to properly send photons among procs: " << send_list.size() << endl;
  }

}

/****  Send properties of the sites  ****/
//send the properties of sites that belong to this proc 
//but also exist on another proc to other procs
void SimpleX::send_site_properties(){


  //temporary structure to hold site information
  Send_Site tempSend;
  //vector to hold the site properties that need to be send
  vector< Send_Site > sites_to_send;
  //vector to hold received sites
  vector< Send_Site > sites_to_update;

  bool total_sites_to_send = 1;

  //loop while there's sites to send 
  SITE_ITERATOR it=sites.begin(); 
  while( total_sites_to_send ){

    //vector to hold the sites that should be send 
    //from this proc to other proc
    vector<unsigned int> nsend_local( COMM_SIZE, 0 );

    //while the maximum number of sites to send is not yet reached, 
    //search for sites that need to be send
    //note that if a site should be send to more than one proc,
    //the chunk will be slightly bigger than the maximum size,
    //but never more than COMM_SIZE
    while( sites_to_send.size() < max_msg_to_send && it!=sites.end() ){

      //only include sites on this proc not in boundary
      if( it->get_process() == COMM_RANK && !it->get_border() ){

  //keep track of the procs the site has already been sent to
        vector< bool > proc(COMM_SIZE,0);
  //loop over all neighbours to check whether 
  //this site has neighbours on another proc
        for( unsigned int j=0; j<it->get_numNeigh(); j++ ){

    //place in sites array of this neighbour
          unsigned int neigh = it->get_neighId(j);
    //check if neighbour belongs to other proc and if the site is not already sent there
          if( sites[neigh].get_process() != COMM_RANK && !proc[sites[neigh].get_process()] ){

      //put neighbour in list to send
            tempSend.set_vertex_id( it->get_vertex_id() );
      //watch out, this is the process the site will be send to,
      //not the process the site belongs to!
            tempSend.set_process( sites[neigh].get_process() );
            tempSend.set_site_id( it->get_site_id() );
            tempSend.set_ballistic( (unsigned int) it->get_ballistic() );

      //store properties in vector
            sites_to_send.push_back( tempSend );
      //keep track of the number of sites to send to
      //which proc
            nsend_local[ sites[neigh].get_process() ]++;
      //set bool at the proc the site will be send to
            proc[ sites[neigh].get_process() ] = 1;

          }//if neighbour is on other proc
        }//for all neighbours
        proc.clear();
      }//if on this proc in domain
      it++;
    }//while sites to send is not too big

    //sort the sites to send 
    sort( sites_to_send.begin(), sites_to_send.end(), compare_process_send_site);

    //define the offset of different procs
    vector< unsigned int > offset( COMM_SIZE, 0 );
    for( unsigned int p=1; p<COMM_SIZE; p++ ){
      offset[p] = offset[p - 1] + nsend_local[p - 1];
    }

    //gather all the local numbers of sites to send to all procs
    vector< unsigned int > nsend( COMM_SIZE * COMM_SIZE, 0 );
    MPI::COMM_WORLD.Allgather(&nsend_local[0], COMM_SIZE, MPI::UNSIGNED, &nsend[0], COMM_SIZE, MPI::UNSIGNED );

    //calculate number of sites to receive
    unsigned int num_sites_to_receive = 0;
    for( unsigned int p=0; p<COMM_SIZE; p++){
      if( p != COMM_RANK ){
        num_sites_to_receive += nsend[ p*COMM_SIZE + COMM_RANK ];
      }
    }

    //give the vector correct size
    sites_to_update.resize( num_sites_to_receive );

    int PTask,recvTask;
    int NTask = (int) COMM_SIZE;
    int ThisTask = (int) COMM_RANK;
    // calculate PTask,  the smallest integer that satisfies COMM_SIZE<=2^PTask
    for(PTask = 0; NTask > (1 << PTask); PTask++) ;

    //loop over all tasks
    for(int level = 1; level < (1 << PTask); level++){
      //vector to hold buffer with number of sites
      //already send to that proc
      vector< unsigned int > nbuffer( COMM_SIZE, 0 );

      int ngrp;
      for(ngrp = level; ngrp < (1 << PTask); ngrp++){

        recvTask = ThisTask ^ ngrp;

        if(recvTask < NTask){
    //check if there's sites to receive or to send
          if(nsend[ ThisTask * NTask + recvTask ] > 0 || nsend[ recvTask * NTask + ThisTask ] > 0){

      //do the sending
            MPI::COMM_WORLD.Sendrecv(&sites_to_send[ offset[recvTask] ], nsend_local[recvTask] * sizeof(Send_Site),
              MPI::BYTE, recvTask, 0, 
              &sites_to_update[ nbuffer[ThisTask] ], nsend[recvTask * NTask + ThisTask] * sizeof(Send_Site),
              MPI::BYTE, recvTask, 0);

          }//if there's sites to send or receive

        }//if receiving task is smaller than COMM_SIZE

  //update the buffer size
        for( int p = 0; p < NTask; p++){
          if( (p ^ ngrp) < NTask ){
            nbuffer[p] += nsend[(p ^ ngrp) * NTask + p];
          }
        }

      }
      //why is this? It has to be there...
      level = ngrp - 1;
      nbuffer.clear();
    } //for all levels

    sites_to_send.clear();
    nsend_local.clear();
    nsend.clear();
    offset.clear();

    //now all information that is needed has been received,
    //so update the sites that are needed on this proc but
    //belong to other proc. The id's of these sites are listed 
    //in the send_list

    //iterator to loop through sites that need their information updated 
    vector<Send_Site>::iterator sit;
    //loop over all sites that this proc has send
    for( sit=sites_to_update.begin(); sit!=sites_to_update.end(); sit++ ){
      bool found = 0;
      //as long as the site has not been found on this proc, 
      //loop through the send_list
      for( unsigned int i=0; !found && i<send_list.size(); i++ ){ 
  //index of the site in the sites array
        unsigned int index = send_list[i];
  //if the vertex id of this site agrees with that of the site
  //on the other proc, assign relevant information
        if( sites[index].get_vertex_id() == sit->get_vertex_id() ){
          found = 1;
    //the site now knows the place of this site in the 
    //sites list on the proc the site belongs to!
          sites[index].set_site_id( sit->get_site_id() );
          sites[index].set_ballistic( (bool) sit->get_ballistic() );
        }

      }//for all sites to update
    }//for all updates

    //clear the vector
    sites_to_update.clear();

    //check if a proc still has sites to send left
    int local_sites_to_send = ( it == sites.end() ) ? 0 : 1;
    int sum;
    //MPI::LOR doesn't work on paracluster :S
    MPI::COMM_WORLD.Allreduce(&local_sites_to_send, &sum, 1, MPI::INT, MPI::SUM );
    if(sum == 0){
      total_sites_to_send = 0;
    }


  }//while still sites to send

  if( COMM_RANK == 0 ){
    simpleXlog << "  Site properties sent " << endl;
  }

}


/**** Send neighbour information  ****/
//send the properties of the neighbours of sites that belong to this proc
//but also exist on another proc to other procs
void SimpleX::send_neighbour_properties(){

  if( COMM_RANK == 0 ){
    simpleXlog << "  Sending neighbour properties" << endl;
  }

  if(COMM_RANK == 0){
    for(SITE_ITERATOR it=sites.begin();it!=sites.end();it++){
      SITE_ITERATOR tmp = it;
      tmp++;
      if( tmp->get_vertex_id() == it->get_vertex_id() ){
        cerr << " (" << COMM_RANK << ") double vertex: (" << tmp->get_process() << "|" << tmp->get_vertex_id() << " and (" << it->get_process() << "|" << it->get_vertex_id() << ")" << endl; 
      }      
    }
  }
  MPI::COMM_WORLD.Barrier();
  if(COMM_RANK == 1){
    for(SITE_ITERATOR it=sites.begin();it!=sites.end();it++){
      SITE_ITERATOR tmp = it;
      tmp++;
      if(tmp->get_vertex_id() == it->get_vertex_id()){
        cerr << " (" << COMM_RANK << ") double vertex: (" << tmp->get_process() << "|" << tmp->get_vertex_id() << " and (" << it->get_process() << "|" << it->get_vertex_id() << ")" << endl; 
      }
    }
  }
  
  
// if(COMM_RANK == 0){
//   for(SITE_ITERATOR it=sites.begin();it!=sites.end();it++){
//     if(it->get_vertex_id() == 47){
//       cerr << "  (" << COMM_RANK << ") site " << it->get_vertex_id() << ": (" << it->get_process() << ": " << it->get_x() << "," << it->get_y() << "," << it->get_z() << ") " 
//            << it->get_border() << endl;
//          for(unsigned int j=0; j<it->get_numNeigh();j++){
//            cerr << sites[it->get_neighId(j)].get_vertex_id() << " (" << sites[it->get_neighId(j)].get_process() << ": " 
//            << sites[it->get_neighId(j)].get_x() << "," << sites[it->get_neighId(j)].get_y() << "," << sites[it->get_neighId(j)].get_z() << ") " << endl;
//          }
//          cerr << endl;
//            
//     }
//     if(it->get_vertex_id() == 13119){
//       cerr << "  (" << COMM_RANK << ") site " << it->get_vertex_id() << ": (" << it->get_process() << ": " << it->get_x() << "," << it->get_y() << "," << it->get_z() << ") " 
//            << it->get_border() << endl;
//     }
//   }
// }

MPI::COMM_WORLD.Barrier();

// if(COMM_RANK == 1){
//   for(SITE_ITERATOR it=sites.begin();it!=sites.end();it++){
//     if(it->get_vertex_id() == 47){
//       cerr << "  (" << COMM_RANK << ") site " << it->get_vertex_id() << ": (" << it->get_process() << ": " << it->get_x() << "," << it->get_y() << "," << it->get_z() << ") " 
//            << it->get_border() << endl;
//            for(unsigned int j=0; j<it->get_numNeigh();j++){
//              cerr << sites[it->get_neighId(j)].get_vertex_id() << " (" << sites[it->get_neighId(j)].get_process() 
//                << ": " << sites[it->get_neighId(j)].get_x() << "," << sites[it->get_neighId(j)].get_y() << "," << sites[it->get_neighId(j)].get_z() << ")" << endl;
//            }
//            cerr << endl;
//     }
//     if(it->get_vertex_id() == 13119){
//       cerr << "  (" << COMM_RANK << ") site " << it->get_vertex_id() << ": (" << it->get_process() << ": " << it->get_x() << "," << it->get_y() << "," << it->get_z() << ") " 
//            << it->get_border() << endl;
//     }
//   }
// }

  //first create the straight array for all sites on other proc
  //that are ballistic
  //loop over all sites that need to get information
  for( unsigned int i=0; i<send_list.size(); i++ ){ 
    //site id
    unsigned int index = send_list[i];
    //only create the straight array in case the site is ballistic
    if( sites[index].get_ballistic() ){
      sites[ index ].create_straight();
    }
  }

  //temporary structure to hold neighbour information
  Send_Neigh tempNeigh;
  //vector to hold neighbour info that needs to be send
  vector< Send_Neigh > neigh_to_send;
  //vector to hold received neighbours
  vector<Send_Neigh> neigh_recv;

  //sets if one of the procs still has sites to send
  bool total_neigh_to_send = 1;

  //loop over all sites
  SITE_ITERATOR it=sites.begin(); 
  while( total_neigh_to_send ){

    //vector to hold the sites that should be send 
    //from this proc to every other proc
    vector<unsigned int> nsend_local( COMM_SIZE, 0 );

    //while the maximum number of sites to send is not yet reached, 
    //search for sites that need to be send
    //note that if a site should be send to more than one proc,
    //the chunk will be slightly bigger than the maximum size,
    //but never more than COMM_SIZE
    while( neigh_to_send.size() < max_msg_to_send && it!=sites.end() ){

      //only include sites on this proc not in boundary
      if( it->get_process() == COMM_RANK && !it->get_border() && it->get_ballistic() ){

        //keep track of the procs the site has already been sent to
        vector< bool > proc(COMM_SIZE,0);

        //loop over all neighbours to check whether 
        //this site has neighbours on another proc
        for( unsigned int j=0; j<it->get_numNeigh(); j++ ){
          //place in sites array of this neighbour
          unsigned long long int neigh = it->get_neighId(j);
          //check if neighbour belongs to other proc and if the site is not already sent there
          if( sites[neigh].get_process() != COMM_RANK && !proc[sites[neigh].get_process()] ){

            //in this case we have to send all neighbours that belong to other proc
            for( unsigned int k=0; k<it->get_numNeigh(); k++ ){
              //check if this neighbour belongs to other proc as well
              if( sites[ it->get_neighId(k) ].get_process() == sites[neigh].get_process() ){
              //store the needed properties of the neighbours of this site

                //watch out, this is the process the site will be send to,
                //not the process the site belongs to!
                tempNeigh.set_process( sites[neigh].get_process() );

                //global id of this site
                tempNeigh.set_vertex_id( it->get_vertex_id() );
                //process of this site
                tempNeigh.set_ballistic( (unsigned int) it->get_ballistic() );
                //global id of neighbour
                tempNeigh.set_neighId( sites[ it->get_neighId(k) ].get_vertex_id() );
                //local place of neighbour in neighbour array of this site
                tempNeigh.set_neighIdLoc( k );

                //put the neighbour information in the structure to be send
                neigh_to_send.push_back( tempNeigh );

                //keep track of the number of sites to send to
                //which proc
                nsend_local[ sites[neigh].get_process() ]++;

                //set bool at the proc the site will be send to
                proc[ sites[neigh].get_process() ] = 1;
              }
            }//for all neighbours
          }//if neighbour not on this proc
        }//for all neighbours
      }//if ballistic on this proc and not in border
      it++;
    }//while neigh to send is not too big

    //sort the sites to send 
    sort( neigh_to_send.begin(), neigh_to_send.end(), compare_process_send_neigh );

    //define the offset of different procs
    vector< unsigned int > offset( COMM_SIZE, 0 );
    for( unsigned int p=1; p<COMM_SIZE; p++ ){
      offset[p] = offset[p - 1] + nsend_local[p - 1];
    }

    //gather all the local numbers of sites to send to all procs
    vector< unsigned int > nsend( COMM_SIZE * COMM_SIZE, 0 );
    MPI::COMM_WORLD.Allgather(&nsend_local[0], COMM_SIZE, MPI::UNSIGNED, &nsend[0], COMM_SIZE, MPI::UNSIGNED );

    //calculate number of sites to receive
    unsigned int num_neigh_to_receive = 0;
    for( unsigned int p=0; p<COMM_SIZE; p++){
      if( p != COMM_RANK ){
        num_neigh_to_receive += nsend[ p*COMM_SIZE + COMM_RANK ];
      }
    }

    //give the vector correct size
    neigh_recv.resize( num_neigh_to_receive );

    int PTask,recvTask;
    int NTask = (int) COMM_SIZE;
    int ThisTask = (int) COMM_RANK;
    // calculate PTask,  the smallest integer that satisfies COMM_SIZE<=2^PTask
    for(PTask = 0; NTask > (1 << PTask); PTask++) ;

    //loop over all tasks
    for(int level = 1; level < (1 << PTask); level++){
      //vector to hold buffer with number of sites
      //already send to that proc
      vector< unsigned int > nbuffer( COMM_SIZE, 0 );

      int ngrp;
      for(ngrp = level; ngrp < (1 << PTask); ngrp++){

        recvTask = ThisTask ^ ngrp;

        if(recvTask < NTask){
    //check if there's sites to receive or to send
          if(nsend[ ThisTask * NTask + recvTask ] > 0 || nsend[ recvTask * NTask + ThisTask ] > 0){

      //do the sending
            MPI::COMM_WORLD.Sendrecv(&neigh_to_send[ offset[recvTask] ], nsend_local[recvTask] * sizeof(Send_Neigh),
              MPI::BYTE, recvTask, 0, 
              &neigh_recv[ nbuffer[ThisTask] ], nsend[recvTask * NTask + ThisTask] * sizeof(Send_Neigh),
              MPI::BYTE, recvTask, 0);

          }//if there's sites to send or receive

        }//if receiving task is smaller than COMM_SIZE

  //update the buffer size
        for( int p = 0; p < NTask; p++){
          if( (p ^ ngrp) < NTask ){
            nbuffer[p] += nsend[(p ^ ngrp) * NTask + p];
          }
        }

      }
      //why is this? It has to be there...
      level = ngrp - 1;
      nbuffer.clear();
    } //for all levels

    neigh_to_send.clear();
    nsend_local.clear();
    nsend.clear();
    offset.clear();

    //sort the received neighbours on vertex id
    sort( neigh_recv.begin(), neigh_recv.end(), compare_vertex_id_send_neigh );

    //loop over all received neighbour information
    //this loop is only possible if sites vector and neigh_recv vector
    //are sorted on vertex id
    SITE_ITERATOR site_it = sites.begin();
    vector<Send_Neigh>::iterator sit = neigh_recv.begin();
    while( sit != neigh_recv.end() ){
      bool found = 0;
      if( site_it == sites.end()){
        cerr << " (" << COMM_RANK << ") Error in send_neighbour_properties: vertex " << sit->get_vertex_id() << " not found " << endl;
        MPI::COMM_WORLD.Abort( -1 );
      }
      //check if we have the correct site
      bool found_site = 0;
      if( site_it->get_vertex_id() == sit->get_vertex_id() ){
        found_site = 1;
        for( unsigned int j=0; !found && j<site_it->get_numNeigh(); j++ ){
    //check if we have the correct neighbour
          if( sit->get_neighId() == sites[ site_it->get_neighId( j ) ].get_vertex_id() ){
      //add the local id of this neighbour on proc it belongs to to the straight array
      //of the site on this proc
            site_it->add_straight( j, sit->get_neighIdLoc() );
            found = 1;
          }
        }
      }

      if(found_site && !found){

        if(sit->get_neighId() < origNumSites){
          cerr << " (" << COMM_RANK << ") Error in send_neighbour_properties: found site " << site_it->get_vertex_id() 
               //<< " (" << site_it->get_x() << "," << site_it->get_y() << "," << site_it->get_z() << ")"
            << " but not neighbour " << sit->get_neighId() << endl;
          bool present = 0;
          for( vector<Send_Neigh>::iterator sit2 = neigh_recv.begin(); sit2!=neigh_recv.end(); sit2++ ){
            if (sit->get_neighId() == sit2->get_neighId() ) present = 1;
            break;
          }
          if(present)
            cerr << " (" << COMM_RANK << ") neighbour is in list but was not found"  << endl;
          else
            cerr << " (" << COMM_RANK << ") neighbour is not in list"  << endl;

          MPI::COMM_WORLD.Abort( -1 );
        }else{
          cerr << " (" << COMM_RANK << ") Warning, in send_neighbour_properties: found site " << site_it->get_vertex_id() 
            << " but not neighbour " << sit->get_neighId() << ". This neighbour is in boundary. " << endl;
        }
      }

      //if the neighbour has been found, go the next
      //else, the neighbour belongs to other site
      if(found_site){
        sit++;
      }else{
        site_it++;
      }

    }//for neighbours received

    //clear memory
    neigh_recv.clear();

    //check if a proc still has sites to send left

    //bool local_neigh_to_send = ( it == sites.end() ) ? 0 : 1;
    //MPI::COMM_WORLD.Allreduce(&local_neigh_to_send, &total_neigh_to_send, 1, MPI::BOOL, MPI::LOR );

    //MPI::LOR doesn't work on paracluster...
    int local_neigh_to_send = ( it == sites.end() ) ? 0 : 1;
    int sum;
    //MPI::LOR doesn't work on paracluster :S
    MPI::COMM_WORLD.Allreduce(&local_neigh_to_send, &sum, 1, MPI::INT, MPI::SUM );
    if(sum == 0){
      total_neigh_to_send = 0;
    }


  }//while still neighbours to send

  if( COMM_RANK == 0 ){
    simpleXlog << "  Neighbour properties sent " << endl;
  }

}


/****  Send the intensities among procs between sweeps  ****/
void SimpleX::send_intensities(){

  //temporary structure to hold intensity information
  Send_Intensity temp;
  //vector to hold the intensities that need to be send
  vector<Send_Intensity> sites_to_send;
  //vector to hold the received intensities
  vector<Send_Intensity> sites_recv;


  bool total_sites_to_send = 1;

  //loop while there's still sites to send
  unsigned int i = 0;
  while( total_sites_to_send ){

    //vector to hold the sites that should be send 
    //from this proc to every other proc
    vector<unsigned int> nsend_local( COMM_SIZE, 0 );

    //while the maximum number of sites to send is not yet reached, 
    //search for sites that need to be send
    //note that if a site should be send to more than one proc,
    //the chunk will be slightly bigger than the maximum size,
    //but never more than COMM_SIZE
    while( sites_to_send.size() < max_msg_to_send && i < send_list.size() ){

      //index in the sites vector
      unsigned int index = send_list[i];

      //loop over all directions

      //in case of ballistic transport, intensity has size of number of neighbours;
      //in case of direction conserving transport, intensity has 
      //the size of the tesselation of the unit sphere
      numPixels = ( sites[index].get_ballistic() ) ? sites[index].get_numNeigh() : number_of_directions;


      for( unsigned int j=0; j<numPixels; j++ ){

  //if the site holds intensity, it needs to be send 
        if( sites[index].get_intensityIn( j ) > 0.0 || sites[index].get_intensityOut( j ) > 0.0  ){

    //assign the relevant properties that need to be send 
    //to ensure fast communication

    //the process the information is meant for
          temp.set_process( sites[index].get_process() );
    //the neighbour that the intensity is coming from
          unsigned int neighId = 0;
          if( sites[index].get_ballistic() ){
      //in case of ballistic transport, the local position
      //of the neighbour on teh other proc is stored in
      //the straight vector, that is not used for sites
      //that belong to other proc
            short int l=0;
            neighId = sites[index].get_straight(j,l);
          }else{
            neighId = j;
          }

          temp.set_neighId( neighId );
    //the place of the site in the sites vector 
    //(note that this is the place in the sites vector 
    //on the proc the site belongs to!)
          temp.set_id( sites[index].get_site_id() );
    //incoming intensity
          temp.set_intensityIn( sites[index].get_intensityIn( j ) );
    //ougoing intensity
          temp.set_intensityOut( sites[index].get_intensityOut( j ) );

    //store the information in vector
          sites_to_send.push_back(temp);

    //photons will be on other proc, so remove them here 
    //to conserve photons
          sites[index].set_intensityIn( j, 0.0 ); 
          sites[index].set_intensityOut( j, 0.0 );

    //keep track of the number of sites to send to
    //which proc
          nsend_local[ sites[index].get_process() ]++;

        }//if site holds intensity

      }//for all pixels/neighbours


      i++;
    }//while there's sites to send and maximum msg size has not been reached


    //sort the sites to send 
    sort( sites_to_send.begin(), sites_to_send.end(), compare_process_send_intensity );
    //define the offset of different procs
    vector< unsigned int > offset( COMM_SIZE, 0 );
    for( unsigned int p=1; p<COMM_SIZE; p++ ){
      offset[p] = offset[p - 1] + nsend_local[p - 1];
    }

    //gather all the local numbers of sites to send to all procs
    vector< unsigned int > nsend( COMM_SIZE * COMM_SIZE, 0 );
    MPI::COMM_WORLD.Allgather(&nsend_local[0], COMM_SIZE, MPI::UNSIGNED, &nsend[0], COMM_SIZE, MPI::UNSIGNED );

    //calculate number of sites to receive
    unsigned int num_sites_to_receive = 0;
    for( unsigned int p=0; p<COMM_SIZE; p++){
      if( p != COMM_RANK ){
        num_sites_to_receive += nsend[ p*COMM_SIZE + COMM_RANK ];
      }
    }

    //give the vector correct size
    sites_recv.resize( num_sites_to_receive );

    //now that we have stored the information to be send, we can actually send it

    //first determine how many sites every proc will send
    int PTask, recvTask;
    int ThisTask = (int) COMM_RANK;
    int NTask = (int) COMM_SIZE;

    //determine PTask,the smallest integer that satisfies NTask<=2^PTask
    for(PTask = 0; NTask > (1 << PTask); PTask++) ;  

    //loop over all tasks
    for(int level = 1; level < (1 << PTask); level++){

      //vector to hold buffer with number of sites
      //already send to that proc
      vector< unsigned int > nbuffer( COMM_SIZE, 0 );

      int ngrp;
      for(ngrp = level; ngrp < (1 << PTask); ngrp++){   

        recvTask = ThisTask ^ ngrp;                       

        if(recvTask < NTask){                            
    //check if there's sites to receive or to send
          if(nsend[ ThisTask * NTask + recvTask ] > 0 || nsend[ recvTask * NTask + ThisTask ] > 0){

      //send the intensity information and receive it from other procs in sites_recv
            MPI::COMM_WORLD.Sendrecv(&sites_to_send[ offset[recvTask] ], nsend_local[recvTask] * sizeof(Send_Intensity),
              MPI::BYTE, recvTask, 0, 
              &sites_recv[ nbuffer[ThisTask] ], nsend[recvTask * NTask + ThisTask] * sizeof(Send_Intensity),
              MPI::BYTE, recvTask, 0);

          }//if there's sites to send or receive
        }//if receiving task is smaller than COMM_SIZE

  //update the buffer size
        for( int p = 0; p < NTask; p++){
          if( (p ^ ngrp) < NTask ){
            nbuffer[p] += nsend[(p ^ ngrp) * NTask + p];
          }
        }
      }
      //why is this? It has to be there...
      level = ngrp - 1;
      //clear buffer
      nbuffer.clear();

    }//for all levels

    //clear vectors
    sites_to_send.clear();
    nsend_local.clear();
    nsend.clear();
    offset.clear();


    //All relevant information was received, now give the photons 
    //to the correct sites

    //iterator to loop over the vector holding the intensity information
    vector<Send_Intensity>::iterator it;
    //loop over received intensity information
    for( it=sites_recv.begin(); it!=sites_recv.end(); it++ ){

      //incoming intensities
      if( it->get_intensityIn() ){
        sites[ it->get_id() ].addRadiationDiffIn( it->get_neighId(), it->get_intensityIn() );
      }

      //outgoing intensities
      if( it->get_intensityOut() ){
        sites[ it->get_id() ].addRadiationDiffOut( it->get_neighId(), it->get_intensityOut() );
      }

    }//for all updates

    sites_recv.clear();

    //check if a proc still has sites to send left
    // bool local_sites_to_send = ( i >= send_list.size() ) ? 0 : 1;
    //MPI::COMM_WORLD.Allreduce(&local_sites_to_send, &total_sites_to_send, 1, MPI::BOOL, MPI::LOR );
    int local_sites_to_send = ( i >= send_list.size() ) ? 0 : 1;
    int sum;
    //MPI::LOR doesn't work on paracluster :S
    MPI::COMM_WORLD.Allreduce(&local_sites_to_send, &sum, 1, MPI::INT, MPI::SUM );
    if(sum == 0){
      total_sites_to_send = 0;
    }


  }//while there's intensity to send  

}

/**** Send the sites marked for removal to master ****/
//To do: send only to master
void SimpleX::send_sites_marked_for_removal(){

  //first determine the relevant properties of the sites
  //marked for removal. This is done by every proc for the 
  //sites that live on that proc

  // the minimum and maximum line length
  double local_min_length = FLT_MAX, local_max_length = -FLT_MAX; 
  //the minimum and maximum gradient of the ionised fraction
  double local_min_frac_grad = FLT_MAX, local_max_frac_grad = -FLT_MAX;

  //loop over all sites marked for removal
  for( unsigned int i=0; i<sites_marked_for_removal.size(); i++ ){

    //site index of the site marked for removal 
    unsigned int index = sites_marked_for_removal[i];

    //-- line lengths --//
    // Determine the maximum and minimum cell size (approx. by line length) on this proc. 
    double line_length = sites[index].get_neigh_dist();
    if( line_length > local_max_length ){
      local_max_length = line_length;
    }
    if( line_length < local_min_length ){
      local_min_length = line_length;
    }

    // -- gradient of ionised fractions --//
    // Determine the maximum difference in log(ionised fraction) between this cell and its neighbours
    double min_frac_neigh = FLT_MAX, max_frac_neigh = -FLT_MAX;
    //loop over all neighbours 
    for( unsigned int j=0; j< sites[index].get_numNeigh(); j++ ){

      unsigned int neigh = sites[index].get_neighId( j );
      double ion_frac_neigh = (double) sites[ neigh ].get_n_HII()/( (double) sites[ neigh ].get_n_HI() + (double) sites[ neigh ].get_n_HII() );

      //store the maximum and minimum ionised fraction of all neighbours
      if( ion_frac_neigh > max_frac_neigh )
        max_frac_neigh = ion_frac_neigh;
      if( ion_frac_neigh < min_frac_neigh )
        min_frac_neigh = ion_frac_neigh;
    }

    //ionised fraction of this site
    double ion_frac = (double) sites[ index ].get_n_HII()/( (double) sites[ index ].get_n_HI() + (double) sites[ index ].get_n_HII() );
    //maximum gradient in ionised fraction between this site and its neighbours
    double frac_grad = max( log10( max_frac_neigh ) - log10( ion_frac ), log10( ion_frac ) - log10( min_frac_neigh ) );

    //store the maximum gradient of all sites on this proc
    if( frac_grad > local_max_frac_grad ){
      local_max_frac_grad = frac_grad;
    }
    if( frac_grad < local_min_frac_grad ){
      local_min_frac_grad = frac_grad;
    }

  }//for all sites marked for removal

  //determine minimum and maximum line length among processors
  double min_length, max_length;
  MPI::COMM_WORLD.Allreduce(&local_min_length,&min_length,1,MPI::DOUBLE,MPI::MIN);
  MPI::COMM_WORLD.Allreduce(&local_max_length,&max_length,1,MPI::DOUBLE,MPI::MAX);

  //determine the minimum and maximum gradient of ionised fraction between processors
  double min_frac_grad, max_frac_grad;
  MPI::COMM_WORLD.Allreduce(&local_min_frac_grad,&min_frac_grad,1,MPI::DOUBLE,MPI::MIN);
  MPI::COMM_WORLD.Allreduce(&local_max_frac_grad,&max_frac_grad,1,MPI::DOUBLE,MPI::MAX);

  //vector to hold the sites that are to be removed
  vector< Site_Remove > sites_marked_for_removal_to_send;
  Site_Remove tempRemove;

  //structure that will receive the sites (only filled on master)
  total_sites_marked_for_removal.clear();                                                                                                                    

  //temporary structure to which the sites are sent
  vector< Site_Remove > sites_recv;

  bool total_sites_to_send = 1;

  //loop while there's still sites to send
  unsigned int i = 0;
  while( total_sites_to_send ){

    //vector to hold the sites that should be send 
    //from this proc to every other proc
    vector<unsigned int> nsend_local( COMM_SIZE, 0 );

    //while the maximum number of sites to send is not yet reached, 
    //search for sites that need to be send
    //note that if a site should be send to more than one proc,
    //the chunk will be slightly bigger than the maximum size,
    //but never more than COMM_SIZE
    while( sites_marked_for_removal_to_send.size() < max_msg_to_send && i < sites_marked_for_removal.size() ){

      //site index of the site marked for removal 
      unsigned int index = sites_marked_for_removal[i];

      //difference between size of this cell and the maximum and minumum cell size
      double line_length = sites[index].get_neigh_dist();
      double delta_x = (max_length - line_length)/(max_length - min_length);

      //difference between gradient in ionised fraction of this cell 
      //and the maximum and minumum of all cells
      //first determine the gradient of this cell
      double min_frac_neigh = FLT_MAX, max_frac_neigh = -FLT_MAX; 
      //loop over all neighbours
      for( unsigned int j=0; j< sites[index].get_numNeigh(); j++ ){

        unsigned int neigh = sites[index].get_neighId( j );
        double ion_frac_neigh = (double) sites[ neigh ].get_n_HII()/( (double) sites[ neigh ].get_n_HI() + (double) sites[ neigh ].get_n_HII() );

  //store the maximum and minimum ionised fraction of all neighbours
        if( ion_frac_neigh > max_frac_neigh ){
          max_frac_neigh = ion_frac_neigh;
        }
        if( ion_frac_neigh < min_frac_neigh ){
          min_frac_neigh = ion_frac_neigh;
        }

      }//for all neighbours

      //uncomment if removal criterion is gradient in ionised fraction instead of cell size

      //double ion_frac = sites[ index ].get_ionised_fraction();
      //double frac_grad = max( log10( max_frac_neigh ) - log10( ion_frac ), log10( ion_frac ) - log10( min_frac_neigh ) );//log10( max_frac_neigh ) - log10( min_frac_neigh );

      //update according to gradient of fractions in cell and size of cell. The line length cancels out of the equations
      //double delta_y = ( max_frac_grad - frac_grad )/delta_frac_grad;

      //put the relevant properties of sites marked for removal in vector
      tempRemove.set_process( sites[index].get_process() );
      tempRemove.set_vertex_id( sites[index].get_vertex_id() );
      tempRemove.set_site_id( sites[index].get_site_id() );
      tempRemove.set_update_property( delta_x );


      if(COMM_RANK == 0){
        total_sites_marked_for_removal.push_back( tempRemove );
      }else{

        sites_marked_for_removal_to_send.push_back( tempRemove );
  //keep track of the number of sites to send to
  //which proc
        nsend_local[ 0 ]++;
      }

      i++;

    }//while there's sites to send and maximum msg size has not been reached

    //sort the sites to send 
    //sort( sites_marked_for_removal_to_send.begin(), sites_marked_for_removal_to_send.end());
    //define the offset of different procs
    vector< unsigned int > offset( COMM_SIZE, 0 );
    for( unsigned int p=1; p<COMM_SIZE; p++ ){
      offset[p] = offset[p - 1] + nsend_local[p - 1];
    }

    //gather all the local numbers of sites to send to all procs
    vector< unsigned int > nsend( COMM_SIZE * COMM_SIZE, 0 );
    MPI::COMM_WORLD.Allgather(&nsend_local[0], COMM_SIZE, MPI::UNSIGNED, &nsend[0], COMM_SIZE, MPI::UNSIGNED );

    //calculate number of sites to receive
    unsigned int num_sites_to_receive = 0;
    for( unsigned int p=0; p<COMM_SIZE; p++){
      if( p != COMM_RANK ){
        num_sites_to_receive += nsend[ p*COMM_SIZE + COMM_RANK ];
      }
    }

    //give the vector correct size
    sites_recv.resize( num_sites_to_receive );


    //================ Send the properties of sites marked for removal to master =====================//

    int PTask, recvTask;
    int ThisTask = COMM_RANK;
    int NTask = COMM_SIZE;

    //Determine PTask, the smallest integer that satisfies NTask <=2^PTask
    for(PTask = 0; NTask > (1 << PTask); PTask++) ; 

    //loop over all tasks
    for(int level = 1; level < (1 << PTask); level++){

      //vector to hold buffer with number of sites
      //already send to that proc
      vector< unsigned int > nbuffer( COMM_SIZE, 0 );

      int ngrp;
      for(ngrp = level; ngrp < (1 << PTask); ngrp++){   

        recvTask = ThisTask ^ ngrp;                       

        if(recvTask < NTask){                            
    //check if there's sites to receive or to send
          if(nsend[ ThisTask * NTask + recvTask ] > 0 || nsend[ recvTask * NTask + ThisTask ] > 0){

      //send the intensity information and receive it from other procs in sites_recv
            MPI::COMM_WORLD.Sendrecv(&sites_marked_for_removal_to_send[ offset[recvTask] ], nsend_local[recvTask] * sizeof(Site_Remove),
              MPI::BYTE, recvTask, 0, 
              &sites_recv[ nbuffer[ThisTask] ], nsend[recvTask * NTask + ThisTask] * sizeof(Site_Remove),
              MPI::BYTE, recvTask, 0);

          }//if there's sites to send or receive
        }//if receiving task is smaller than COMM_SIZE

  //update the buffer size
        for( int p = 0; p < NTask; p++){
          if( (p ^ ngrp) < NTask ){
            nbuffer[p] += nsend[(p ^ ngrp) * NTask + p];
          }
        }
      }
      //why is this? It has to be there...
      level = ngrp - 1;
      //clear buffer
      nbuffer.clear();

    }//for all levels

    sites_marked_for_removal_to_send.clear();
    nsend_local.clear();
    nsend.clear();
    offset.clear();

    //only the master selects the vertices to remove                                                                                                           
    if(COMM_RANK == 0){                                                                                                                                        
      //first add sites on this proc to list                                                                                                                   
      //total_sites_marked_for_removal.insert( total_sites_marked_for_removal.end(), 
      //					     sites_marked_for_removal_to_send.begin(),  sites_marked_for_removal_to_send.end() );
      //now add the sites that were send here 
      total_sites_marked_for_removal.insert( total_sites_marked_for_removal.end(), sites_recv.begin(), sites_recv.end() ); 

    }  

    //clear vector
    sites_recv.clear();

    //check if a proc still has sites to send left
    //bool local_sites_to_send = ( i >= sites_marked_for_removal.size() ) ? 0 : 1;
    //MPI::COMM_WORLD.Allreduce(&local_sites_to_send, &total_sites_to_send, 1, MPI::BOOL, MPI::LOR );

    int local_sites_to_send = ( i >= sites_marked_for_removal.size() ) ? 0 : 1;
    int sum;
    //MPI::LOR doesn't work on paracluster :S
    MPI::COMM_WORLD.Allreduce(&local_sites_to_send, &sum, 1, MPI::INT, MPI::SUM );
    if(sum == 0){
      total_sites_to_send = 0;
    }

  }//while there's sites to send  


  if( COMM_RANK == 0 ){
    simpleXlog << "  Total number of sites marked for removal: " << total_sites_marked_for_removal.size() << endl;
  }

}


/****  Send the sites to be removed to all procs  ****/
void SimpleX::send_sites_to_remove(){

  //master proc has sites that need to be removed
  unsigned long long int numSites_marked_for_removal = total_sites_marked_for_removal.size();

  //broadcast it to other procs
  MPI::COMM_WORLD.Bcast(&numSites_marked_for_removal,1,MPI::UNSIGNED_LONG_LONG,0);

  //if there's a lot of sites, send in chunks
  unsigned int num_chunks = 1;
  if( numSites_marked_for_removal > max_msg_to_send ){
    num_chunks = (unsigned int) ceil( (double) numSites_marked_for_removal/max_msg_to_send);
  }

  //make sure receiving vector is empty on all procs except master
  //and give correct size
  if(COMM_RANK != 0){
    total_sites_marked_for_removal.clear();
    vector<Site_Remove>().swap(total_sites_marked_for_removal);
  }

  //send to other procs in chunks
  for(unsigned int i=0; i<num_chunks; i++ ){

    //temporary vector to hold the vertices
    vector< Site_Remove > temp;
    //size of the chunk to send
    unsigned int this_chunk_size = max_msg_to_send;
    //start id of current chunk
    unsigned int start_id = i*max_msg_to_send; 
    //make sure the final chunk has correct size
    if( start_id + this_chunk_size >= numSites_marked_for_removal ){
      this_chunk_size = numSites_marked_for_removal - start_id;
    }

    //master fills temp
    if(COMM_RANK == 0){
      for(unsigned int j = start_id; j<(start_id+this_chunk_size); j++){
        temp.push_back( total_sites_marked_for_removal[j] );
      }
      //other procs resize to correct sizw  
    }else{
      temp.resize( this_chunk_size );  
    }

    //broadcast this chunk
    MPI::COMM_WORLD.Bcast(&temp[0],this_chunk_size*sizeof(Site_Remove),MPI::BYTE,0);

    //procs not master fill vector from temp vector
    if( COMM_RANK != 0 ){
      total_sites_marked_for_removal.insert( total_sites_marked_for_removal.end(), temp.begin(), temp.end() );
    }

    //clear temp for next use
    temp.clear();

  }


  if( COMM_RANK == 0 ){
    simpleXlog << "  Sites to be removed sent to all procs " << endl;
  }

}


/****  Send the updates list of vertices to master  ****/
//vertices no longer exist only on master proc as when the grid
//was first created, so every proc has to send its own list to
//master
//to do: send only to master in a better way
void SimpleX::send_new_vertex_list(){

  //vector that holds the vertices that will be send
  vector< Vertex > vertices_to_send;

  //vector that holds the vertices that will be received
  vector< Vertex > vertices_recv;

  //temporary vector that holds vertices currently on this proc
  vector< Vertex > temp_vertices = vertices;
  vertices.clear();

  bool total_vertices_to_send = 1;

  vector< Vertex >::iterator it=temp_vertices.begin();
  while( total_vertices_to_send ){

    //vector to hold the sites that should be send 
    //from this proc to other proc
    vector<unsigned int> nsend_local( COMM_SIZE, 0 );

    while( vertices_to_send.size() < max_msg_to_send && it!=temp_vertices.end() ){

      if(COMM_RANK != 0){
        vertices_to_send.push_back( *it );
  //keep track of the number of sites to send to
  //which proc
  //send only to master in this case
        nsend_local[ 0 ]++;
      }else{
        vertices.push_back( *it );
      }

      it++;
    }//while there's vertices to send and maximum msg size has not been reached

    //define the offset of different procs
    vector< unsigned int > offset( COMM_SIZE, 0 );
    for( unsigned int p=1; p<COMM_SIZE; p++ ){
      offset[p] = offset[p - 1] + nsend_local[p - 1];
    }

    //gather all the local numbers of sites to send to all procs
    vector< unsigned int > nsend( COMM_SIZE * COMM_SIZE, 0 );
    MPI::COMM_WORLD.Allgather(&nsend_local[0], COMM_SIZE, MPI::UNSIGNED, &nsend[0], COMM_SIZE, MPI::UNSIGNED );

    //calculate number of sites to receive
    unsigned int num_vertices_to_receive = 0;
    for( unsigned int p=0; p<COMM_SIZE; p++){
      if( p != COMM_RANK ){
        num_vertices_to_receive += nsend[ p*COMM_SIZE + COMM_RANK ];
      }
    }

    //give the vector correct size
    vertices_recv.resize( num_vertices_to_receive );

   //================ Send the vertices to master =====================//

    int PTask, recvTask;
    int ThisTask = COMM_RANK;
    int NTask = COMM_SIZE;

    //Determine PTask, the smallest integer that satisfies NTask <=2^PTask
    for(PTask = 0; NTask > (1 << PTask); PTask++) ; 

    //loop over all tasks
    for(int level = 1; level < (1 << PTask); level++){

      //vector to hold buffer with number of sites
      //already send to that proc
      vector< unsigned int > nbuffer( COMM_SIZE, 0 );

      int ngrp;
      for(ngrp = level; ngrp < (1 << PTask); ngrp++){   

        recvTask = ThisTask ^ ngrp;                       

        if(recvTask < NTask){                            
    //check if there's sites to receive or to send
          if(nsend[ ThisTask * NTask + recvTask ] > 0 || nsend[ recvTask * NTask + ThisTask ] > 0){

      //send the intensity information and receive it from other procs in sites_recv
            MPI::COMM_WORLD.Sendrecv(&vertices_to_send[ offset[recvTask] ], nsend_local[recvTask] * sizeof(Vertex),
              MPI::BYTE, recvTask, 0, 
              &vertices_recv[ nbuffer[ThisTask] ], nsend[recvTask * NTask + ThisTask] * sizeof(Vertex),
              MPI::BYTE, recvTask, 0);

          }//if there's sites to send or receive
        }//if receiving task is smaller than COMM_SIZE

  //update the buffer size
        for( int p = 0; p < NTask; p++){
          if( (p ^ ngrp) < NTask ){
            nbuffer[p] += nsend[(p ^ ngrp) * NTask + p];
          }
        }
      }
      //why is this? It has to be there...
      level = ngrp - 1;
      //clear buffer
      nbuffer.clear();

    }//for all levels

    vertices_to_send.clear();
    nsend_local.clear();
    nsend.clear();
    offset.clear();

    //only the master gets the vertices
    if(COMM_RANK == 0){                                                                                                                                        
      vertices.insert( vertices.end(), vertices_recv.begin(), vertices_recv.end() ); 
    }

    vertices_recv.clear();

    int local_sites_to_send = ( it == temp_vertices.end() ) ? 0 : 1;
    int sum;
    //MPI::LOR doesn't work on paracluster :S
    MPI::COMM_WORLD.Allreduce(&local_sites_to_send, &sum, 1, MPI::INT, MPI::SUM );
    if(sum == 0){
      total_vertices_to_send = 0;
    }

  }

  temp_vertices.clear();

  //set the correct number of sites
  numSites = vertices.size();

  if( COMM_RANK == 0 ){
    simpleXlog << "  Created new list of vertices " << endl;
    simpleXlog << "  New number of sites: " << numSites << endl;
  }

}


/****  Send the physical params of the sites to procs that need them ****/
//During vertex movement sites might have been migrated to 
//another proc, so the properties of these sites have
//to be send
void SimpleX::send_site_physics(){

//   //temporary structure to hold site information
//   Site_Update tempUpdate;
//   //vector to hold the site properties that need to be send
//   vector< Site_Update > sites_to_send;
//   //vector to temporarily store the sites that were on this proc
//   vector< Site_Update > temp_site_properties = site_properties;

//   //clear site_properties to later fill it with sites that belong to this proc
//   site_properties.clear();

//   bool total_sites_to_send = 1;

//   //loop while there's sites to send 
//   vector< Site_Update >::iterator it=temp_site_properties.begin(); 
//   while( total_sites_to_send ){

//     //vector to hold the sites that should be send 
//     //from this proc to other proc
//     vector<unsigned int> nsend_local( COMM_SIZE, 0 );

//     //while the maximum number of sites to send is not yet reached, 
//     //search for sites that need to be send
//     //note that if a site should be send to more than one proc,
//     //the chunk will be slightly bigger than the maximum size,
//     //but never more than COMM_SIZE
//     while( sites_to_send.size() < max_msg_to_send && it!=temp_site_properties.end() ){


//     }//while message not too big

//   }//while total_site_to_send

//   //determine which sites have moved from this proc, and need
//   //to be send

//   //create mapping from vertex id to place in sites vector
//   vector< unsigned int > mapping_sites( vertex_id_max, vertex_id_max+1);
//   for( SITE_ITERATOR it=sites.begin(); it!=sites.end(); it++ ){
//     if( it->get_process() == COMM_RANK && !it->get_border() ){
//       mapping_sites[ it->get_vertex_id() ] = it->get_site_id();
//     }
//   }

//   //fill the sites_to_send vector with the properties the other proc needs
//   vector<Site_Update> sites_to_send;
//   Site_Update temp;

//   //loop over all sites
//   for( unsigned int i=0; i<site_properties.size(); i++ ){

//     //get place in sites array from mapping
//     unsigned int site_id = mapping_sites[ site_properties[i].get_vertex_id() ];
//     //if site not on this proc, site_id is larger than numSites 
//     //and properties have to be send

//    if( site_id > vertex_id_max ){

//       //assign site_properties to temp
//       temp = site_properties[i];
//       //put it in array to be send
//       sites_to_send.push_back(temp);

//     }//if no valid id ont his proc
//   }//for all site_properties


//   //first send the number of sites that will be sent
//   unsigned int number_of_sites_send = sites_to_send.size();
//   vector< unsigned int > number_of_sites_recv( COMM_SIZE, 0);

//   int PTask, recvTask;
//   int ThisTask = COMM_RANK;
//   int NTask = COMM_SIZE;

//   //determine PTask, the smallest integer that satisfies NTask <=2^PTask
//   for(PTask = 0; NTask > (1 << PTask); PTask++) ;

//   //loop over all tasks
//   for(int level = 1; level < (1 << PTask); level++){
//     for(int ngrp = level; ngrp < (1 << PTask); ngrp++){  

//       recvTask = ThisTask ^ ngrp;                       

//       if(recvTask < NTask){                             

// 	//send and receive number of sites to be send
// 	MPI::COMM_WORLD.Sendrecv(&number_of_sites_send, 1,
// 				 MPI::UNSIGNED, recvTask, 0, 
// 				 &number_of_sites_recv[recvTask], 1,
// 				 MPI::UNSIGNED, recvTask ,0);

//       }
//     }
//   }


//   //Now site information can be send, since every proc know its size

//   //vector to hold the intensities
//   vector<Site_Update>* properties_recv;

//   //give vector correct size
//   properties_recv = new vector<Site_Update>[COMM_SIZE];
//   for( unsigned int i=0; i<COMM_SIZE; i++ ){
//     properties_recv[i].resize( number_of_sites_recv[i] );
//   }

//   //loop over all tasks
//   for(int level = 1; level < (1 << PTask); level++){
//     for(int ngrp = level; ngrp < (1 << PTask); ngrp++){   

//       recvTask = ThisTask ^ ngrp;                       

//       if(recvTask < NTask){                            

// 	//send and receive the site intensities
// 	MPI::COMM_WORLD.Sendrecv(&sites_to_send[0], number_of_sites_send*sizeof(Site_Update),
// 				 MPI::BYTE, recvTask, 0, 
// 				 &properties_recv[recvTask][0], number_of_sites_recv[recvTask]*sizeof(Site_Update),
// 				 MPI::BYTE, recvTask ,0);

//       }
//     }
//   }

//   //free memory
//   number_of_sites_recv.clear();
//   sites_to_send.clear();

//   //Update the sites with the received intensities


//   //loop over all received properties and add properties
//   //to the vector on this proc if site is here
//   vector<Site_Update>::iterator it;
//   //loop over all procs
//   for( unsigned int p=0; p<COMM_SIZE; p++){
//     //exclude this proc
//     if(p!=COMM_RANK){
//       //loop over all intensities that were send to this proc
//       for( it=properties_recv[p].begin(); it!=properties_recv[p].end(); it++ ){
// 	unsigned int site_id = mapping_sites[ it->get_vertex_id() ];
// 	//if site with this vertex id is on this proc, site_id 
// 	//is smaller than number of sites
// 	if( site_id < sites.size() ){

// 	  site_properties.push_back( *it );


// 	}//if site on this proc
//       }//for all received props from p
//       properties_recv[p].clear();
//     }//if not this proc
//   }//for all procs

//   //free memory
//   if(properties_recv){
//     delete [] properties_recv;
//     properties_recv = NULL;
//   }

//   mapping_sites.clear();

  if( COMM_RANK == 0 ){
    simpleXlog << "  Site physics sent " << endl;
  }

}


/****  Send the stored intensities to procs that need them  ****/
//During vertex movement sites might have been migrated to 
//another proc, so the intensities of these sites have
//to be send
void SimpleX::send_stored_intensities(){

  if( COMM_RANK == 0 ){
    simpleXlog << "  Sending stored intensities " << endl;
  }

  //set correct size of numPixels
  numPixels = number_of_directions;

  //determine which sites have moved to this proc, and need intensities

  //create a mapping from vertex id to place in sites array
  vector< unsigned int > mapping_sites( vertex_id_max, vertex_id_max+1);
  for( SITE_ITERATOR it=sites.begin(); it!=sites.end(); it++ ){
    if( it->get_process() == COMM_RANK && !it->get_border() ){
      mapping_sites[ it->get_vertex_id() ] = it->get_site_id();
    }
  }

  //fill the sites_to_send vector with the intensity information the other proc needs
  // neighId is the place in the reference pixel array
  // id is the vertex id of the site
  // intensityOut is the intensity belonging to this pixel of the relevant site
  vector<Send_Intensity> sites_to_send;
  Send_Intensity temp;
  //if entry in intens_ids points to a vertex not on this processor, it should be send
  for( unsigned int i=0; i<intens_ids.size(); i++ ){
    //place in sites array obtained from vertex id
    unsigned int site_id = mapping_sites[ intens_ids[i] ];
    //if site not on this proc, site_id is larger than numSites 
    //and intensities have to be send
    if( site_id > vertex_id_max ){
      //first entry of the site intensities in site_intensities vector
      unsigned int index = number_of_directions*i;
      //loop over all directions
      for( unsigned int j=0; j<numPixels; j++ ){
  //give relevant properties to temp
        temp.set_process( COMM_RANK );
        temp.set_neighId( j );
        temp.set_id( intens_ids[i] );
        temp.set_intensityIn( 0.0 );
        temp.set_intensityOut( site_intensities[ index + j ] );

  //put it in array to be send
        sites_to_send.push_back(temp);

      }//for all directions
    }//if site has moved from this proc

  }//for all sites formerly on this proc

  //first send the number of sites that will be sent
  unsigned int number_of_sites_send = sites_to_send.size();
  vector< unsigned int > number_of_sites_recv( COMM_SIZE, 0);

  int PTask, recvTask;
  int ThisTask = COMM_RANK;
  int NTask = COMM_SIZE;


  //determine PTask, the smallest integer that satisfies NTask <=2^PTask
  for(PTask = 0; NTask > (1 << PTask); PTask++) ;

  //loop over all tasks
  for(int level = 1; level < (1 << PTask); level++){
    for(int ngrp = level; ngrp < (1 << PTask); ngrp++){  

      recvTask = ThisTask ^ ngrp;                       

      if(recvTask < NTask){                             

  //send and receive number of sites to be send
        MPI::COMM_WORLD.Sendrecv(&number_of_sites_send, 1,
          MPI::UNSIGNED, recvTask, 0, 
          &number_of_sites_recv[recvTask], 1,
          MPI::UNSIGNED, recvTask ,0);

      }
    }
  }


  //Now site information can be send, since every proc know its size

  //vector to hold the intensities
  vector<Send_Intensity>* intensities_recv;

  //give vector correct size
  intensities_recv = new vector<Send_Intensity>[COMM_SIZE];
  for( unsigned int i=0; i<COMM_SIZE; i++ ){
    intensities_recv[i].resize( number_of_sites_recv[i] );
  }

  //loop over all tasks
  for(int level = 1; level < (1 << PTask); level++){
    for(int ngrp = level; ngrp < (1 << PTask); ngrp++){   

      recvTask = ThisTask ^ ngrp;                       

      if(recvTask < NTask){                            

  //send and receive the site intensities
        MPI::COMM_WORLD.Sendrecv(&sites_to_send[0], number_of_sites_send*sizeof(Send_Intensity),
          MPI::BYTE, recvTask, 0, 
          &intensities_recv[recvTask][0], number_of_sites_recv[recvTask]*sizeof(Send_Intensity),
          MPI::BYTE, recvTask ,0);

      }
    }
  }

  //free memory
  number_of_sites_recv.clear();
  sites_to_send.clear();

  //Update the sites with the received intensities

  //create vector to map vertex ids to place in intens_ids vector
  vector<unsigned int> mapping_intens(vertex_id_max,vertex_id_max+1);
  for(unsigned int i=0; i<intens_ids.size(); i++){
    mapping_intens[ intens_ids[i] ] = i;
  }

  //loop over all received intensities and add intensities
  //to the site_intensities on this proc 
  vector<Send_Intensity>::iterator it;
  //loop over all procs
  for( unsigned int p=0; p<COMM_SIZE; p++){
    //exclude this proc
    if(p!=COMM_RANK){
      //loop over all intensities that were send to this proc
      for( it=intensities_recv[p].begin(); it!=intensities_recv[p].end(); it++ ){
        unsigned int site_id = mapping_sites[ it->get_id() ];
  //if site with this vertex id is on this proc, site_id 
  //is smaller than number of sites
        if( site_id < sites.size() ){

    //check if the site intensities are already created
          if( mapping_intens[ it->get_id() ] > intens_ids.size() ){
      //if not, add the site before filling

      //keep mapping up to date
            mapping_intens[ it->get_id() ] = intens_ids.size();
            intens_ids.push_back( it->get_id() );
      //give site_intensities extra entries for the new intensities
            site_intensities.insert( site_intensities.end(), numPixels, 0.0 );
          }
    //put the intensity in correct place

    //neighId is place in intensityOut/In array
          unsigned int pos = it->get_neighId();
    //get correct place in site_intensities array
          pos += mapping_intens[it->get_id()]*numPixels;
    //add intensity to vector
          site_intensities[pos] += it->get_intensityOut();

        }

      }//for intensities received from proc p

      //clear memory
      intensities_recv[p].clear();

    }//if not htis proc
  }//for all procs

  //free memory
  if(intensities_recv){
    delete [] intensities_recv;
    intensities_recv = NULL;
  }

  mapping_sites.clear();
  mapping_intens.clear();

}

/****  Send the ballistic properties of the changed sites to other proc  ****/
void SimpleX::send_site_ballistics( const vector< unsigned long long int >& sites_to_send ){

   //send the properties of the vertex on this proc to other proc
  Send_Neigh tempNeigh;
  vector< Send_Neigh > neigh_to_send;
  //vector to hold received neighbours
  vector<Send_Neigh> neigh_recv;

  //sets if one of the procs still has sites to send
  bool total_neigh_to_send = 1;

  //loop over all sites
  unsigned long long int i=0;
  while( total_neigh_to_send ){

    //vector to hold the sites that should be send 
    //from this proc to every other proc
    vector<unsigned int> nsend_local( COMM_SIZE, 0 );

    //while the maximum number of sites to send is not yet reached, 
    //search for sites that need to be send
    //note that if a site should be send to more than one proc,
    //the chunk will be slightly bigger than the maximum size,
    //but never more than COMM_SIZE
    while( neigh_to_send.size() < max_msg_to_send && i<sites_to_send.size() ){

      //site id of the site to be send
      unsigned long long int site_id = sites_to_send[i];

      //keep track of the procs the site has already been sent to
      vector< bool > proc(COMM_SIZE,0);
      //loop over all neighbours to check to which 
      //procs this site has to be send
      for( unsigned int j=0; j<sites[site_id].get_numNeigh(); j++ ){

  //place in sites array of this neighbour
        unsigned long long int neigh = sites[site_id].get_neighId(j);
  //check if neighbour belongs to other proc and if the site is not already sent there
        if( sites[neigh].get_process() != COMM_RANK && !proc[sites[neigh].get_process()] ){

    //check if the site to be send is ballistic
          if( sites[site_id].get_ballistic() ){
      //if the site is now ballistic, all neighbours need to be send

      //loop over all neighbours
            for( unsigned int k=0; k<sites[site_id].get_numNeigh(); k++ ){
              if( sites[ sites[site_id].get_neighId(k) ].get_process() == sites[neigh].get_process() ){
    //store the needed properties of the neighbours of this site

    //watch out, this is the process the site will be send to,
    //not the process the site belongs to!
                tempNeigh.set_process( sites[neigh].get_process() );

    //global id of this site
                tempNeigh.set_vertex_id( sites[site_id].get_vertex_id() );
    //process of this site
                tempNeigh.set_ballistic( 1 );
    //global id of neighbour
                tempNeigh.set_neighId( sites[ sites[site_id].get_neighId(k) ].get_vertex_id() );
    //local place of neighbour in neighbour array of this site
                tempNeigh.set_neighIdLoc( k );

    //put the neighbour information in the structure to be send
                neigh_to_send.push_back( tempNeigh );

    //keep track of the number of sites to send to
    //which proc
                nsend_local[ sites[neigh].get_process() ]++;

    //set bool at the proc the site will be send to
                proc[ sites[neigh].get_process() ] = 1;
              }
            }//for all neighbours

          }else{

      //if the site to be send is not ballistic, only ballistic() needs
      //to be send
            tempNeigh.set_ballistic( 0 );
      //global id of this site
            tempNeigh.set_vertex_id( sites[site_id].get_vertex_id() );

      //watch out, this is the process the site will be send to,
      //not the process the site belongs to!
            tempNeigh.set_process( sites[neigh].get_process() );

      //global id of neighbour
            tempNeigh.set_neighId( sites[ neigh ].get_vertex_id() );
      //local place of neighbour in neighbour array of this site
            tempNeigh.set_neighIdLoc( 0 );

      //put the neighbour information in the structure to be send
            neigh_to_send.push_back( tempNeigh );

      //keep track of the number of sites to send to
      //which proc
            nsend_local[ sites[neigh].get_process() ]++;

      //set bool at the proc the site will be send to
            proc[ sites[neigh].get_process() ] = 1;

          }//if ballistic
        }//if neighbour is on other proc and site is not send there yet
      }//for all neighbours

      i++;
    }//while message not too big

    //sort the sites to send 
    sort( neigh_to_send.begin(), neigh_to_send.end(), compare_process_send_neigh );

    //define the offset of different procs
    vector< unsigned int > offset( COMM_SIZE, 0 );
    for( unsigned int p=1; p<COMM_SIZE; p++ ){
      offset[p] = offset[p - 1] + nsend_local[p - 1];
    }

    //gather all the local numbers of sites to send to all procs
    vector< unsigned int > nsend( COMM_SIZE * COMM_SIZE, 0 );
    MPI::COMM_WORLD.Allgather(&nsend_local[0], COMM_SIZE, MPI::UNSIGNED, &nsend[0], COMM_SIZE, MPI::UNSIGNED );

    //calculate number of sites to receive
    unsigned int num_neigh_to_receive = 0;
    for( unsigned int p=0; p<COMM_SIZE; p++){
      if( p != COMM_RANK ){
        num_neigh_to_receive += nsend[ p*COMM_SIZE + COMM_RANK ];
      }
    }

    //give the vector correct size
    neigh_recv.resize( num_neigh_to_receive );

    int PTask,recvTask;
    int NTask = (int) COMM_SIZE;
    int ThisTask = (int) COMM_RANK;

    // calculate PTask,  the smallest integer that satisfies COMM_SIZE<=2^PTask
    for(PTask = 0; NTask > (1 << PTask); PTask++) ;

     //loop over all tasks
    for(int level = 1; level < (1 << PTask); level++){
      //vector to hold buffer with number of sites
      //already send to that proc
      vector< unsigned int > nbuffer( COMM_SIZE, 0 );

      int ngrp;
      for(ngrp = level; ngrp < (1 << PTask); ngrp++){

        recvTask = ThisTask ^ ngrp;

        if(recvTask < NTask){
    //check if there's sites to receive or to send
          if(nsend[ ThisTask * NTask + recvTask ] > 0 || nsend[ recvTask * NTask + ThisTask ] > 0){

      //do the sending
            MPI::COMM_WORLD.Sendrecv(&neigh_to_send[ offset[recvTask] ], nsend_local[recvTask] * sizeof(Send_Neigh),
              MPI::BYTE, recvTask, 0, 
              &neigh_recv[ nbuffer[ThisTask] ], nsend[recvTask * NTask + ThisTask] * sizeof(Send_Neigh),
              MPI::BYTE, recvTask, 0);

          }//if there's sites to send or receive

        }//if receiving task is smaller than COMM_SIZE

  //update the buffer size
        for( int p = 0; p < NTask; p++){
          if( (p ^ ngrp) < NTask ){
            nbuffer[p] += nsend[(p ^ ngrp) * NTask + p];
          }
        }

      }
      //why is this? It has to be there...
      level = ngrp - 1;
      nbuffer.clear();
    } //for all levels

    neigh_to_send.clear();
    nsend_local.clear();
    nsend.clear();
    offset.clear();

    //sort the received neighbours on vertex id
    sort( neigh_recv.begin(), neigh_recv.end(), compare_vertex_id_send_neigh );
    //loop over all received neighbour information
    //this loop is only possible if sites vector and neigh_recv vector
    //are sorted on vertex id
    SITE_ITERATOR site_it = sites.begin();
    vector<Send_Neigh>::iterator sit = neigh_recv.begin();
    while( sit != neigh_recv.end() ){
      bool found = 0;
      //check if we have the correct site
      if( site_it->get_vertex_id() == sit->get_vertex_id() ){

  //check if received site is ballistic
        if( sit->get_ballistic() ){

    //if site is not yet set to ballistic, this means it has
    //not been visited yet, so create arrays
          if( site_it->get_ballistic() == 0 ){

            site_it->set_ballistic(1);

      //delete the intensity arrays
            site_it->delete_intensityOut();
            site_it->delete_intensityIn();

      //create new intensity arrays with correct size
            site_it->create_intensityIn( site_it->get_numNeigh() );
            site_it->create_intensityOut( site_it->get_numNeigh() );

      //create new straight array
            site_it->delete_straight();
            site_it->create_straight();


          }//if site not visited yet

    //check if neighbour has been found
          for( unsigned int j=0; !found && j<site_it->get_numNeigh(); j++ ){
      //check if we have the correct neighbour
            if( sit->get_neighId() == sites[ site_it->get_neighId( j ) ].get_vertex_id() ){
              found = 1;
        //add the local id of this neighbour on proc it belongs to to the straight array
        //of the site on this proc
              site_it->add_straight( j, sit->get_neighIdLoc() );

            }//if right neighbour
          }//for all neighbours


        }else{

    //if received site is not ballistic, change size of intensity vector
          site_it->set_ballistic( 0 );

    //delete the intensity arrays
          site_it->delete_intensityOut();
          site_it->delete_intensityIn();

    //create new intensity arrays with correct size
          site_it->create_intensityIn( number_of_directions );
          site_it->create_intensityOut( number_of_directions );

    //straight array no longer needed
    //site_it->delete_straight();

          found = 1;

        }//if ballistic
      }//if vertex ids match

      //if the neighbour has been found, go the next
      //else, the neighbour belongs to other site
      if(found){
        sit++;
      }else{
        site_it++;
      }

    }//for neighbours received

    neigh_recv.clear();

    //check if a proc still has sites to send left
//     bool local_neigh_to_send = ( i >= sites_to_send.size() ) ? 0 : 1;
//     MPI::COMM_WORLD.Allreduce(&local_neigh_to_send, &total_neigh_to_send, 1, MPI::BOOL, MPI::LOR );

    //check if a proc still has sites to send left
    int local_neigh_to_send = ( i >= sites_to_send.size() ) ? 0 : 1;
    int sum;
    //MPI::LOR doesn't work on paracluster :S
    MPI::COMM_WORLD.Allreduce(&local_neigh_to_send, &sum, 1, MPI::INT, MPI::SUM );
    if(sum == 0){
      total_neigh_to_send = 0;
    }

  }//while still neighbours to send

}



/****************************************************************************************/
/*                            Physics Functions                                         */
/****************************************************************************************/

void SimpleX::compute_physics( const unsigned int& run ){

  double t0 = MPI::Wtime();

  if( run == 0 ){

    initialise_physics_params();

    if(fillChoice==AUTOMATIC){
      set_homogeneous_number_density( homDens );
    }else{ //fillChoice == READ
      assign_read_properties();
    }

  }else{

    //no new domain decomposition is done, so this is not necessary
    //send the intensities to other sites that need them
    //send_stored_intensities();

    //return the stored intensities to the sites
    //use the orientation index with which the intensities were stored
    return_intensities();
    //also return intensities to the ballistic sites
    return_ballistic_intensities();

    //empty vectors that are no longer needed
    site_intensities.clear();
    intens_ids.clear();

    //to completely empty site_intensities, swap it with an empty vector
    vector< float >().swap( site_intensities );
    vector< unsigned long long int >().swap( intens_ids );

  }

  output_number_of_atoms();

  //parameter_output( run );

  double t1 = MPI::Wtime();
  if( COMM_RANK == 0 )
    simpleXlog << endl << "  Calculating physics parameters took " << t1-t0 << " seconds" << endl << endl;

}

//Initialize physical parameters
void SimpleX::initialise_physics_params() {  

  //physical units
  simTime *= secondsPerMyr;  
  UNIT_T *= secondsPerMyr;  
  UNIT_L = sizeBox * parsecToCm; 
  UNIT_V = pow(UNIT_L, 3.0);

  if(blackBody) {   // -> Blackbody spectrum
    black_body_source( sourceTeff );
  } else {             // -> Monochromatic spectrum
    cross_H = 6.3e-18;
  }

  if( COMM_RANK == 0 )
    simpleXlog << "  Essential physical parameters calculated " << endl;

}

//Set the homogeneous number density.
//Only used in case of automatic filling
void SimpleX::set_homogeneous_number_density(const float& nH){

  //unit mass is average number density * unit volume
  UNIT_D = nH; 

  UNIT_M = UNIT_D * UNIT_V;

  for( SITE_ITERATOR it=sites.begin(); it!=sites.end(); it++ ){
    if( it->get_process() == COMM_RANK){
      it->set_n_HI( 1.0 );
      it->set_n_HII( 0.0 );
    }
  }//for all sites

  if( COMM_RANK == 0 ){
    simpleXlog << "  Homogeneous number density created " << endl;
  }

}

/****  Assign the number densities and fluxes from read in to sites  ****/
void SimpleX::assign_read_properties(){

  //unit conversion were done when creating inout file, so UNIT_M and
  //UNIT_I are already known

  //loop over all sites
  for( SITE_ITERATOR it=sites.begin(); it!=sites.end(); it++ ){
    if( it->get_process() == COMM_RANK){
      if( !it->get_border() ){
  //get the number density and flux from list that was 
  //created when reading the input
        it->set_n_HI( temp_n_HI_list[ it->get_vertex_id() ] );
  //set ionised fraction
        it->set_n_HII( temp_n_HII_list[ it->get_vertex_id() ] );
  //change if there's more than one freq!
        it->set_flux( temp_flux_list[ it->get_vertex_id() ] );
      }else{
        it->set_n_HI( 0.0 );
        it->set_n_HII( 0.0 );
        it->set_flux( 0.0 );
      }
    }
  }

  //clear the lists
  temp_n_HI_list.clear();
  vector< float >().swap( temp_n_HI_list );
  temp_n_HII_list.clear();
  vector< float >().swap( temp_n_HII_list );
  temp_flux_list.clear();
  vector< float >().swap( temp_flux_list );

  if( COMM_RANK == 0 ){
    simpleXlog << "  Assigned number density and flux to sites " << endl;
  }

}

/****  Return the physical params to the sites  ****/
//this function only works if site_properties contains all sites on this proc
//and no others
void SimpleX::return_physics(){

  //sort the site_properties
  sort( site_properties.begin(), site_properties.end(), compare_vertex_id_site_update );

  //loop over all sites
  unsigned long long int i=0;
  for( SITE_ITERATOR it=sites.begin(); it!=sites.end(); it++ ){
    //only include sites on this proc that are not in border
    if( it->get_process() == COMM_RANK ){ 
      if( !it->get_border() ) {
  //check if the vertex ids match
        if(it->get_vertex_id() == site_properties[i].get_vertex_id() ){

    //assign properties
          it->set_n_HI( site_properties[ i ].get_n_HI() );
          it->set_n_HII( site_properties[ i ].get_n_HII() ); 
          it->set_flux( site_properties[ i ].get_flux() );
          it->set_ballistic( site_properties[ i ].get_ballistic() );

          i++;
        }else{
          cerr << " (" << COMM_RANK << ") Error in return_physics(): site not mapped" << endl;
          cerr << it->get_vertex_id() << " " << i << endl;
          MPI::COMM_WORLD.Abort( -1 );
        }
      }else{
  //in case of direction conserving transport, set boundary to 
  //direction conserving
        if( dirConsTransport ){
          it->set_ballistic( 0 );
        }else{
    //in all other cases, set boundary to ballistic
          it->set_ballistic( 1 );
        }

      }//if not in border
    }//if on this proc 
  }//for all sites


  if( COMM_RANK == 0 ){
    simpleXlog << " Returned physics to proper sites" << endl;
  }

}


//calculate effective cross section and spectrum of black body source
void SimpleX::black_body_source(const double& tempSource) {

//   double frac = 0.0;
//   double redTemp = tempSource/1e4;

//   //BBeval values are located in Commmon.cpp
//   if( redTemp > (double) BBevalNum ) 
//     redTemp = (double) BBevalNum;
//   else if( redTemp < 1.0) 
//     redTemp = 1.0;

//   int indTemp = (int) floor( redTemp );
//   double ratio = redTemp - (double) indTemp;

//   frac = ratio * BBevalEffAbs[ indTemp-2 ] + (1.0 - ratio ) * BBevalEffAbs[ indTemp-1 ];

//   cross = frac * A0;

  //Determine the frequency integrated cross section of hydrogen for a black bosy source of temperature tempSource
  //Units are scaled to close to 1

  //upper and lower bounds of integration in terms of the ionisation frequency for hydrogen
  const double upperBound = 1e2;
  const double lowerBound = 1.0;

  //integral of hydrogen cross section times Planck curve
  double intHI = qromb(fHI, lowerBound, upperBound, tempSource ) ;
  //integral over Planck curve divided by nu
  double bbOverNu = qromb(PlanckDivNu, lowerBound, upperBound, tempSource);
  //integral over Planck curve
  double bb = qromb(Planck, lowerBound, upperBound, tempSource);

  //divide the two to get the desired effective cross section
  intHI /= bbOverNu;


  //back to physical units
  cross_H = intHI*crossFactor;

  //In case dust is included, calculate effective cross section
  if(dust_model == SMC){

    double intDust = qromb(f_dust_SMC,lowerBound,upperBound,tempSource);
    intDust /= bb;

    cross_dust = intDust * sigma_dust_SMC;

    //cerr << " Old cross section: " << intDust*sigma_dust_SMC*bb/bbOverNu << "; updated cross section: " << cross_dust << endl;

  }else if(dust_model == LMC){

    double intDust = qromb(f_dust_LMC,lowerBound,upperBound,tempSource);
    intDust /= bb;

    cross_dust = intDust * sigma_dust_LMC;


  }else if(dust_model == MW){

    double intDust = qromb(f_dust_MW,lowerBound,upperBound,tempSource);
    intDust /= bb;

    cross_dust = intDust * sigma_dust_MW;

  }else{
    cross_dust = 0.0;
  }


  //Determine dust to gas ratio
  //reference metallicity
  double Z_ref=1.0;
  switch(dust_model){

    //values for metallicity MCs from Welty et al. 1997 and 1999
    case SMC:
      //assuming overall metallicity -0.6 dex relative to solar
    Z_ref = 0.25;
    break;
    case LMC:
      //assuming overall metallicity -0.3 dex relative to solar
    Z_ref = 0.5;
    break;
    case MW:
      //assume solar metallicity
    Z_ref = 1.0;
    break;
  }

  dust_to_gas_ratio = metallicity/Z_ref;

  if( COMM_RANK == 0 ){
    simpleXlog << "  Calculated effective cross sections for black body source" << endl;
  }

}


//compute the total number of atoms on this proc, used for output
void SimpleX::output_number_of_atoms() {  

  double atoms=0.0;  
  totalAtoms=0.0;
  for( SITE_ITERATOR it=sites.begin(); it!=sites.end(); it++ ){
    if( !it->get_border() && it->get_process() == COMM_RANK ){ 
      double N_HI  = (double) it->get_n_HI() * UNIT_D * (double) it->get_volume() * UNIT_V;
      double N_HII = (double) it->get_n_HII() * UNIT_D * (double) it->get_volume() * UNIT_V;
      atoms += N_HI + N_HII;  
    }
  }  

  MPI::COMM_WORLD.Reduce(&atoms,&totalAtoms,1,MPI::DOUBLE,MPI::SUM,0);

  if( COMM_RANK == 0 )
    simpleXlog << "  Total number of atoms on grid calculated and sent " << endl;

}  

//compute the mean optical depth on this proc, used for output
double SimpleX::output_optical_depth() {  

  double tau=0.0; 
  unsigned int count = 0; 
  for( SITE_ITERATOR it=sites.begin(); it!=sites.end(); it++ ){
    if( !it->get_border() && it->get_process() == COMM_RANK ){ 
      tau += it->get_n_HI() * UNIT_D * (double) it->get_neigh_dist() * UNIT_L * cross_H * straight_correction_factor;
      count++;  
    }
  }  

  double totalTau;
  MPI::COMM_WORLD.Allreduce(&tau,&totalTau,1,MPI::DOUBLE,MPI::SUM);
  unsigned int totalCount;
  MPI::COMM_WORLD.Allreduce(&count,&totalCount,1,MPI::UNSIGNED,MPI::SUM);

  double mean_tau = 0.0;
  if(totalCount > 0){
    mean_tau = totalTau/( (double) totalCount );
  }

  return mean_tau;

} 

//output relevant physical quantities to screen
void SimpleX::parameter_output( const unsigned int& run ){

  // int maxRes;

  // MPI::COMM_WORLD.Reduce(&localMaxRes,&maxRes,1,MPI::INT,MPI::MIN,0);

  double mean_tau = output_optical_depth();

  if(COMM_RANK == 0){

    cerr << endl << endl
      << "  Total volume: " << totalVolume << " (should be close to " << 1.0 << ")." << endl
      << "  Max resolution: " << maxRes << "^3, smallest scale resolved is " << UNIT_L/(maxRes*parsecToCm) << " pc" << endl
      << endl         
      << "  Effective Photo-Ionization Cross Section     : " << cross_H << " cm^2" << endl;
    if( dust_model != NO_DUST){
      cerr << "  Effective dust cross section                 : " << cross_dust << " cm^2 "<< endl;
    }  
    cerr << endl
      << "  Size of Simulation Domain                    : " << sizeBox << " pc" << endl  
      << "  Time elapsed                                 : " << run*simTime/( numRuns*secondsPerMyr ) << " Myr" << endl
      << "  This run                                     : " << simTime/( numRuns*secondsPerMyr ) << " Myr" << endl
      << "  Total Simulation Time                        : " << simTime/secondsPerMyr << " Myr" << endl  
      << endl  
      << "  Average density SimpleX grid                 : " << totalAtoms/pow( UNIT_L, 3 ) << " cm^-3" << endl  
      << "  Total number of atoms within SimpleX grid    : " << totalAtoms << endl  
      << endl
      << "  Average optical depth                        : " << mean_tau << endl
      << endl;

    simpleXlog << endl << endl
      << "  Total volume: " << totalVolume << " (should be close to " << 1.0 << ")." << endl
      << "  Max resolution: " << maxRes << "^3, smallest scale resolved is " << UNIT_L/(maxRes*parsecToCm) << " pc" << endl
      << endl         
      << "  Effective Photo-Ionization Cross Section     : " << cross_H << " cm^2" << endl;
    if( dust_model != NO_DUST){
      simpleXlog << "  Effective dust cross section                 : " << cross_dust << " cm^2 "<< endl;
    }  
    simpleXlog << endl  
      << "  Average density SimpleX grid                 : " << totalAtoms/pow( UNIT_L, 3 ) << " cm^-3" << endl  
      << "  Total number of atoms within SimpleX grid    : " << totalAtoms << endl  
      << endl
      << "  Average optical depth                        : " << mean_tau << endl
      << endl;

  }


}

//Put source in the centre of domain
void SimpleX::set_source_in_centre(){


  if( COMM_RANK == 0 )
    cerr << " (" << COMM_RANK << ") Placing a source with flux " << source_strength*UNIT_I << " in centre of domain, ";  

  double temp;  
  int tellerLocal=0;
  int tellerGlobal;

  for( SITE_ITERATOR it=sites.begin(); it!=sites.end(); it++ ){
    if( it->get_process() == COMM_RANK && !it->get_border() ){
      temp = sqrt( pow( it->get_x() - source_x , 2 ) + pow( it->get_y() - source_y, 2 ) + pow( it->get_z() - source_z, 2 ) );  
      if( temp<=source_radius ) tellerLocal++;  
    }
  }

  MPI::COMM_WORLD.Allreduce(&tellerLocal, &tellerGlobal, 1, MPI::INT, MPI::SUM );

  for( SITE_ITERATOR it=sites.begin(); it!=sites.end(); it++ ){
    if( it->get_process() == COMM_RANK  && !it->get_border() ){
      temp = sqrt( pow( it->get_x() - source_x , 2 ) + pow( it->get_y() - source_y, 2 ) + pow( it->get_z() - source_z, 2 ) );  
      if( temp<=source_radius ){
        double flux = (double) source_strength/tellerGlobal;
        it->set_flux( (float) flux );  
      }
    } 
  }  

  if( COMM_RANK == 0 ){
    cerr << "consisting of " << tellerGlobal; 
    simpleXlog << "  Source in centre contains " << tellerGlobal;
    if(tellerGlobal == 1){
      cerr << " point." << endl;
      simpleXlog << " point." << endl;
    }else{
      cerr << " points." << endl;
      simpleXlog << " points." << endl;
    }
  }

}


//return recombination coefficient
double SimpleX::recombCoeff( const double& tempGas) {  

  //double recomb_coefficient_old = alpha_B * pow( tempIonGas/tempGas, recombCoefficient);  

  //recombination coefficient from Hui&Gnedin, 1998
  double lambda_HI = 2.0 * 1.57807e5/tempGas;
  double recomb_coefficient = 2.753e-14 * pow( lambda_HI, 1.5 ) / pow( ( 1.0 + pow( 0.364963503 * lambda_HI, 0.407 )), 2.242 );  

  return recomb_coefficient;

}

//calculate recombination
double SimpleX::recombine(){

  double totalRecombs=0.0;  

  for( SITE_ITERATOR it=sites.begin(); it!=sites.end(); it++ ){
    if( it->get_process() == COMM_RANK && !it->get_border() ){

      // ( number of recombinations ) = ( free electrons ) * ( free protons ) * ( recombination coefficient )
      //t_rec_inv = ion_frac * dens * recombCoeff( gasIonTemp );  //ion_frac*dens = number of available protons in the volume

      double N_HI = (double) it->get_n_HI() * UNIT_D * (double) it->get_volume() * UNIT_V;
      double N_HII = (double) it->get_n_HII() * UNIT_D * (double) it->get_volume() * UNIT_V;
      double t_rec_inv = (double) it->get_n_HII() * UNIT_D * recombCoeff( gasIonTemp );  //ion_frac*dens = number of available protons in the volume
      double dxBYdt = N_HII * t_rec_inv;  //number of recombinations per second
      double dx = dxBYdt * UNIT_T;  

      if( dx > N_HII ) {
        dx = N_HII;  
      }

      N_HI += dx;
      N_HII -= dx;
      double tmp = N_HI/(UNIT_D * (double) it->get_volume() * UNIT_V);
      it->set_n_HI( (float) tmp );
      tmp = N_HII/(UNIT_D * (double) it->get_volume() * UNIT_V);
      it->set_n_HII( (float) tmp );

      totalRecombs += dx;  

    }
  } //for all vertices 

  return totalRecombs;

}

/****************************************************************************************/
/*                       Radiative Transport Functions                                  */
/****************************************************************************************/



//Calculate the ionisations and recombination
double SimpleX::solve_rate_equation( Site& site ){

  //calculate ionised fraction, neutral fraction, number density, optical depth 
  double n_H = ( (double) site.get_n_HI() + (double) site.get_n_HII() );

  double initial_ionised_fraction = (double) site.get_n_HII()/n_H;
  //double initial_neutral_fraction = (double) site.get_n_HI()/n_H;

  double initial_tau = (double) site.get_n_HI() * UNIT_D * (double) site.get_neigh_dist() * UNIT_L * cross_H * straight_correction_factor;

  //dirty hack to avoid nans
  if(initial_tau == 0.0 && n_H > 0.0){
    initial_tau = FLT_MIN;//1.e-15;
  }

  double tau_dust = 0.0;
  if(dust_model != NO_DUST){
    if(dust_sublimation){
      //dust scales with only neutral gas column density, dust is destroyed by radiation
      tau_dust = dust_to_gas_ratio * (double) site.get_n_HI() * UNIT_D * (double) site.get_neigh_dist() * UNIT_L * cross_dust * straight_correction_factor;
    }else{
      //dust scales with total gas column density, no dust is destroyed
      tau_dust = dust_to_gas_ratio * n_H * UNIT_D * (double) site.get_neigh_dist() * UNIT_L * cross_dust * straight_correction_factor;
    }
  }

  double tau_eff = initial_tau + tau_dust;

  //total number of incoming and outgoing photons
  double N_in_total = 0.0;
  double N_out_total = 0.0;

  //number of neutral and ionised atoms in the cell
  double initial_N_HI = (double) site.get_n_HI() * UNIT_D * (double) site.get_volume() * UNIT_V;
  double initial_N_HII = (double) site.get_n_HII() * UNIT_D * (double) site.get_volume() * UNIT_V;


  //in case of ballistic transport, intensity has size of number of neighbours;
  //in case of direction conserving transport, intensity has 
  //the size of the tesselation of the unit sphere
  numPixels = ( site.get_ballistic() ) ? site.get_numNeigh() : number_of_directions;


  for( unsigned int j=0; j < numPixels; j++) {
    if( site.get_intensityOut( j ) > 0.0 ) { 
      N_in_total += (double) site.get_intensityOut( j );
    } //if intensityOut
  } //for all neighbours

  //if rho is zero, no need to go into loop
  if(n_H == 0.0){
    N_in_total *= UNIT_I;
    N_out_total = N_in_total;
  }


  if(n_H > 0.0){

    //transform the number of photons to physical units
    N_in_total *= UNIT_I;
  
    
    //calculate the number of absorptions and ionisations if N_HI and N_HII were constant during RT time step
    //double initial_N_retained_total = (1.0 - exp(-initial_tau)) * N_in_total;
    //double initial_numIonised = ( initial_N_retained_total >= N_HI ) ? N_HI : initial_N_retained_total;

    //absorptions by gas and dust
    double initial_N_retained_total = (1.0 - exp(-tau_eff)) * N_in_total;
    //fraction available for H
    double initial_N_retained_H = initial_N_retained_total*initial_tau/tau_eff;

    double initial_numIonised = ( initial_N_retained_H >= initial_N_HI ) ? initial_N_HI : initial_N_retained_H;

    //use temporal photon condervation or not?
    if(photon_conservation){
      //calculate relevant time scales to see whether assumption of constant N_HI and N_HII is correct
      //double t_ion = UNIT_T * N_HI/initial_N_retained_total; //ionisation time scale
      //double t_ion = UNIT_T * initial_N_HI/initial_N_retained_H; //ionisation time scale
      //double t_rec = ( site.get_volume() * pow( UNIT_L, 3.0 ) )/( ( initial_N_HII + initial_numIonised ) * recombCoeff( gasIonTemp ) ); //recombination time scale
      //double t_eq = t_ion*t_rec/(t_ion + t_rec); //'equality' time scale


      //ionisation time scale, if no ionisations it is taken to be much larger than time step
      double t_ion = ( initial_N_retained_H > 0.0) ? UNIT_T * initial_N_HI/initial_N_retained_H : 1.e5*UNIT_T; 

      //total #ions available for recombinations
      double tot_ions = initial_N_HII + initial_numIonised;
      //recombination time scale, if no recombinations it is taken to be much larger than time step
      double t_rec = (tot_ions > 0.0 ) ? ( site.get_volume() * pow( UNIT_L, 3.0 ) )/
        ( tot_ions * recombCoeff( gasIonTemp ) ) : 1.e5*UNIT_T;

      double t_eq = t_ion*t_rec/(t_ion + t_rec); //'equality' time scale


      if(t_eq != t_eq){
        simpleXlog << endl << " (" << COMM_RANK << ") NAN detected! t_eq: " << t_eq << " " << t_ion << " " << t_rec << endl;
        simpleXlog << "    initial_N_retained_total: " << initial_N_retained_total 
          << "  initial_tau: " << initial_tau << " initial_ionised_fraction: " << initial_ionised_fraction 
          << " n_H: " << n_H << endl;
        exit(-1);
      }

      // photoionisation rate

      //double Gamma = initial_N_retained_total/(UNIT_T * N_HI); 
      //double Gamma = initial_N_retained_H/(UNIT_T * initial_N_HI); 


      //subcycling is necessary when the RT time step is larger than 
      //some factor times ionisation time scale or recombination time scale 
      unsigned long long int Nsteps = 1;
      double factor = 0.01;
      if( UNIT_T > factor*t_eq ){
        Nsteps = (unsigned long long int) ceil( UNIT_T/(factor*t_eq));
      }

      // if(Nsteps > 5){
      // 	cerr << " (" << COMM_RANK << ") Nsteps: " << Nsteps << endl;
      // }

//       unsigned long long int max_Nsteps = 100000;
//       if( Nsteps > max_Nsteps ){
// 	Nsteps = max_Nsteps;
//       } 


//       cerr << " t_eq: " << t_eq << " t_ion: " << t_ion << " t_rec: " << t_rec 
// 	   << " Nsteps: " << Nsteps << " N_ret: " << initial_N_retained_total 
// 	   << " N_HI: " << N_HI << endl;//" (" << site.get_x() << "," << site.get_y() << "," << site.get_z() << ")" << endl;

      //keep track of number of recombinations and ionisations
      double number_of_ionisations = 0.0;
      double number_of_recombinations= 0.0;

      //neutral fraction and ionised fraction during subcycle step
      //double neutral_fraction_step = initial_neutral_fraction;
      //double ionised_fraction_step = initial_ionised_fraction;

      double N_HI_step = initial_N_HI;
      double N_HII_step = initial_N_HII;

      //time step during subcycle
      double dt = UNIT_T/Nsteps;

      //speed up computation time:

      //effective time step in calculations of ionised atoms
      double dt_eff = dt/UNIT_T;
      //one over the initial neuatral atoms
      double one_over_initial_N_HI = (initial_N_HI > 0.0 ) ? 1.0/initial_N_HI : 1.0;
      //one over the physical volume
      double one_over_volume = 1.0/( site.get_volume() * UNIT_V );

      //-------------------------  Loop over subcycle steps  ---------------------------------//
      for( unsigned int i=0; i<Nsteps; i++ ){

  //time step, optical depth and photoionisation rate during subcycle step
  //double tau_step = (initial_tau/tau_eff)*( neutral_fraction_step )/( initial_neutral_fraction );
  //double Gamma_step = Gamma * ( 1.0 - exp(-tau_step) ) * ( initial_neutral_fraction )/( ( 1.0 - exp(-initial_tau) ) * ( neutral_fraction_step ) );

  //double tau_step = tau_eff*( neutral_fraction_step )/( initial_neutral_fraction );
  //double Gamma_step = Gamma * ( 1.0 - exp(-tau_step) ) * ( initial_neutral_fraction )/( ( 1.0 - exp(-tau_eff) ) * ( neutral_fraction_step ) );

  //number of ionisations during subcylce step
  //double numIonised = Gamma_step*dt*N_HI;
  //double numIonised = Gamma_step*dt*N_HI_step*initial_tau/tau_eff;

        double tau_step = tau_eff*N_HI_step*one_over_initial_N_HI;
        double numIonised = N_in_total * ( 1.0 - exp(-tau_step) ) * dt_eff;


  //check for unphysical values
        if( N_HI_step < 0 ||  numIonised < 0){
          simpleXlog << " (" << COMM_RANK << ") ERROR: negative quantity; numIonised: " <<  numIonised << ", N_HI: " << N_HI_step << " " << i << endl;
          simpleXlog << site.get_vertex_id();
          MPI::COMM_WORLD.Abort(-1);
        }

  //make sure number of ionisations is less than number of neutral atoms. Only necessary if subcycle step is too large
        if( numIonised > N_HI_step){
          simpleXlog << " Warning, subcycle step too large: numIonised: " << numIonised << " N_HI_step: " << N_HI_step << endl;
          numIonised = N_HI_step;
        }

  //keep track of ionisations
        number_of_ionisations += numIonised;

  //adjust number of neutrals and ioniseds accordingly
        N_HI_step  -= numIonised;
        N_HII_step += numIonised;

  //calculate new ionised fraction
  //ionised_fraction_step = N_HII_step/(N_HI_step + N_HII_step);

  //calculate number of recombinations  
        double t_rec_inv = N_HII_step * recombCoeff( gasIonTemp )*one_over_volume ;  //ion_frac*dens = number of available protons in the volume
        double dxBYdt = N_HII_step * t_rec_inv;  //number of recombinations per second
        double dx = dxBYdt * dt;  

  //if flag for recombinations is not set, set dx to zero
        if(!recombination){
          dx = 0;
        }

  //actual recombinations
        if( dx > N_HII_step){
          dx = N_HII_step;
          cerr << " Warning, subcycle step too large: dx: " << dx << " N_HII_step: " << N_HII_step << endl;
        }

        number_of_recombinations += dx;

        N_HI_step += dx;
        N_HII_step -= dx;

  //if photoionisation equilibrium is reached, no need to recalculate these
        if( fabs(numIonised - dx)/dx < 1e-5 && (int) i < (int) (Nsteps-2) ){

    //number of steps left to take
          unsigned int Nleft = Nsteps - i - 1;

          number_of_ionisations += numIonised*Nleft;
          number_of_recombinations += dx*Nleft;

          i = Nsteps-1;

        }//if photoionsation equilibrium

      }//for all steps

      // if(site.get_vertex_id() == 6466){
      // 	cerr << " (" << COMM_RANK << ") Number of recombinations: " << number_of_recombinations 
      // 	     << " number of ionisations: " << number_of_ionisations << " ratio: " << number_of_ionisations/number_of_recombinations << endl;
      // }

      //set the new number_densities
      double tmp = N_HI_step/( UNIT_D * site.get_volume() * UNIT_V );
      site.set_n_HI( (float) tmp );
      tmp = N_HII_step/( UNIT_D * site.get_volume() * UNIT_V );
      site.set_n_HII( (float) tmp );

      if( site.get_n_HI() != site.get_n_HI() ){
        simpleXlog << " Nan detected! " << endl;
        simpleXlog << N_HI_step << " " << UNIT_D << " " << site.get_volume() << " " << UNIT_V << endl;
        exit(-1);
        
      }


      //calculate the number of outgoing photons after the subcycling
      N_out_total = (N_in_total - number_of_ionisations);


      if(N_out_total != N_out_total){
        simpleXlog << endl << " (" << COMM_RANK << ") NAN detected! N_out_total: " << N_out_total << " " << N_in_total << " " << number_of_ionisations << endl;
        exit(-1);
      }

    }else{ //don't use temporal photon conservation scheme

      //make sure number of ionisations is less than number of neutral atoms. Only necessary if time step is too large
      // if( initial_numIonised > initial_N_HI){
      // 	initial_numIonised = initial_N_HI;
      // }

      //adjust number of neutrals and ioniseds accordingly
    initial_N_HI  -= initial_numIonised;
    initial_N_HII += initial_numIonised;

    double tmp = initial_N_HI/( UNIT_D * site.get_volume() * UNIT_V );
    site.set_n_HI( (float) tmp );
    tmp = initial_N_HII/( UNIT_D * site.get_volume() * UNIT_V );
    site.set_n_HII( (float) tmp );

    N_out_total = N_in_total - initial_numIonised;

  }

}

  //convert back to numerical units
N_in_total /= UNIT_I;
N_out_total /= UNIT_I;

return N_out_total;

}



void SimpleX::source_transport( Site& site ){

  //flux to send
  double flux_to_send = (double) site.get_flux() * UNIT_T;


  if(source_inside_cell){

    //source ionises its own cell, so put the radiation in the insitensity bins belonging to this site
    //technically this is not entirely correct because if the source is in the centre of the cell radiation
    //would see half the optical depth it does now, but only the source radiation...

    if( site.get_ballistic() ){
      for( unsigned int j=0; j<site.get_numNeigh(); j++ ){
        double inten = flux_to_send/site.get_numNeigh();
        site.addRadiationDiffOut( j, (float) inten ); 
      }
    }else{
      for( unsigned int m=0; m<number_of_directions; m++ ){
        double inten = flux_to_send/number_of_directions;
        site.addRadiationDiffOut( m, (float) inten ); 
      }
    }

    //wat je eigenlijk zou moeten doen is de neighDist op 0.5 zetten en solve_rate_eq aanroepen,maar
    //let er daarbij op dat je de intensityOut die er al is eerst op moet slaan en later terug moet geven!

  }else{

    //source doesn't ionise its own cell, so give radiation to neighbours


    //in case of ballistic transport, intensity has size of number of neighbours
    if( site.get_ballistic() ){
      for( unsigned int j=0; j<site.get_numNeigh(); j++ ){
        double inten = flux_to_send/site.get_numNeigh();        
        sites[ site.get_neighId(j) ].addRadiationDiffOut( site.get_outgoing(j), (float) inten ); 
      }
      //total_inten += (double) it->get_flux() * UNIT_T * UNIT_I;

    }else{
            
      //in case of direction conserving transport and combined transport, intensity has 
      //the size of the tesselation of the unit sphere
      numPixels = number_of_directions;
      for( unsigned int m=0; m<numPixels; m++ ){
        //photons to send
        double inten = flux_to_send/numPixels;

        unsigned int dir = (unsigned) site.get_outgoing(m);
        unsigned int neigh = site.get_neighId( dir );

        //check if the neighbour that gets the radiation is ballistic
        if(sites[ neigh ].get_ballistic() ){
          //if the neighbour is ballistic, find out which 
          //Delaunay line is associated with this direction
          bool found = 0;
          for( unsigned int j=0; !found && j<sites[ neigh ].get_numNeigh(); j++ ){
            if( sites[neigh].get_neighId(j) == site.get_site_id() ){
              found = 1;
              sites[neigh].addRadiationDiffOut( j, (float) inten );
            } 
          }
          if(!found){
            cerr << " error in source radiation" << endl;
            MPI::COMM_WORLD.Abort(-1);
          }
        }else{
    //if the neighbour is not ballistic, simply 
    //add the radiation in the same direction bin
          sites[ neigh ].addRadiationDiffOut( m, (float) inten ); 
        }

      }//for all pixels

    }//if ballistic

  }//if source ionised its own cell

}


void SimpleX::diffuse_transport( Site& site, double& N_out_total ) {

  //redistribute the photons to all neighbours

  //in case of ballistic transport, intensity has size of number of neighbours;
  //in case of direction conserving transport, intensity has 
  //the size of the tesselation of the unit sphere
  //numPixels = ( site.get_ballistic() ) ? site.get_numNeigh() : number_of_directions;


  if( site.get_ballistic() ){

    //if site is ballistic, divide intensity over all neighbours
    double inten = N_out_total/site.get_numNeigh();

    for(unsigned int j=0; j<site.get_numNeigh(); j++){

      //local id of the connection in neighbour's neighbour array
      unsigned int neighIdLoc = (unsigned) site.get_outgoing(j);
      //site id of the neighbour that will receive the photons
      unsigned int neigh = site.get_neighId( j );

      if( sites[neigh].get_ballistic() ){

        sites[ neigh ].addRadiationDiffIn( neighIdLoc, (float) inten ); 

      }else{

  //if the neighbour is not ballistic, find out in which direction 
  //bins the photons should go
        vector< unsigned int > dir_to_use;
  //this trick only works for neighbours on this proc,
  //otherwise outgoing array does not exist
        if( sites[ neigh ].get_process() == COMM_RANK ){
          for( unsigned int n=0; n<number_of_directions; n++ ){
            if( sites[ neigh ].get_outgoing(n) == neighIdLoc ){
              dir_to_use.push_back( n );
            }
          }//for all directions
        }//if neigh on this proc

  //if no direction is associated with this neighbour, find one
        if(dir_to_use.size() == 0){

    //find closest associated direction for this Delaunay line

    //directions are stored in header
          float **refVector;
    // Assign memory to the refVector
          refVector = new float*[number_of_directions];
          for( unsigned int n=0; n<number_of_directions; n++ ){
            refVector[n]=new float[3];
          }
    //assign the first orientation to refVector
          for( unsigned int n=0; n<number_of_directions; n++ ){
            refVector[n][0] = (float) orient[orientation_index][n][0];
            refVector[n][1] = (float) orient[orientation_index][n][1];
            refVector[n][2] = (float) orient[orientation_index][n][2];
          }

          float vectorNeigh[3];
          vectorNeigh[0] = site.get_x() - sites[ neigh ].get_x();
          vectorNeigh[1] = site.get_y() - sites[ neigh ].get_y();
          vectorNeigh[2] = site.get_z() - sites[ neigh ].get_z();

    // Find the neighVec that is closest (in angle) to the refVec
          float highestInprod=-FLT_MAX;
          int highestInprodNumber=-1;
          for(unsigned int n=0; n<number_of_directions; n++) {
            float tempInprod = inproduct(vectorNeigh, refVector[n], 3);
            if( tempInprod > highestInprod ){
              highestInprod = tempInprod;
              highestInprodNumber = n;
            }
          }
          dir_to_use.push_back( highestInprodNumber );

    // free memory of the refVector
          for(unsigned int n=0; n<number_of_directions; n++){
            delete [] refVector[n];
          }
          delete [] refVector;

        }//if no direction associated with this linw

        for(unsigned int n=0; n<dir_to_use.size() ; n++ ){
          sites[ neigh ].addRadiationDiffIn( dir_to_use[n], (float) inten/dir_to_use.size() );
        }//for all neighbours to use

        dir_to_use.clear();

      }//if neighbour is ballistic

    }//for all neighbours


  }else{

    //in case of direction conserving transport and combined transport, intensity has 
    //the size of the tesselation of the unit sphere
    numPixels = number_of_directions;
    for( unsigned int m=0; m<numPixels; m++ ){

      //photons to send
      double inten = N_out_total/numPixels;

      unsigned int dir = (unsigned) site.get_outgoing(m);
      unsigned int neigh = site.get_neighId( dir );

      //check if the neighbour that gets the radiation is ballistic
      if(sites[ neigh ].get_ballistic() ){
  //if the neighbour is ballistic, find out which 
  //Delaunay line is associated with this direction
        bool found = 0;
        for( unsigned int j=0; !found && j<sites[ neigh ].get_numNeigh(); j++ ){
          if( sites[neigh].get_neighId(j) == site.get_site_id() ){
            found = 1;
            sites[neigh].addRadiationDiffIn( j, (float) inten );
          } 
        }
        if(!found){
          cerr << " error in source radiation" << endl;
          MPI::COMM_WORLD.Abort(-1);
        }
      }else{
  //if the neighbour is not ballistic, simply 
  //add the radiation in the same direction bin
        sites[ neigh ].addRadiationDiffIn( m, (float) inten ); 
      }

    }//for all pixels

  }//if ballistic

}


void SimpleX::ballistic_transport( Site& site, double& N_out_total ) { 


  //redistribute the photons to the d most straightforward neighbours

  //in case of ballistic transport, intensity has size of number of neighbours;
  //in case of direction conserving transport, intensity has 
  //the size of the tesselation of the unit sphere
  numPixels = ( site.get_ballistic() ) ? site.get_numNeigh() : number_of_directions;

  //determine N_in_total
  double N_in_total = 0.0;
  for( unsigned int j=0; j < numPixels; j++) {
    if( site.get_intensityOut( j ) > 0.0 ) { 
      N_in_total += (double) site.get_intensityOut( j );
    } //if intensityOut
  } //for all neighbours


  //loop over all neighbours/directions
  for( unsigned int j=0; j < numPixels; j++) {
    //only include directions/neighbours that have intensity to send
    if( site.get_intensityOut( j ) > 0.0 ) { 

      //intensity to send out, get correct part from all added intensities
      double N_out = N_out_total * ( (double) site.get_intensityOut( j ) / N_in_total );

      if( site.get_ballistic() ){

  //in ballistic case, send intensity to d most straightforward Delaunay lines,
  //if they are not pointing backwards
  //in outgoing is the local id of this site in neighbour's neighbour array
        for( short int l=0; l<site.get_numStraight( j ); l++ ){
    //most straighforward neighbour of this direction of this site
          unsigned int straight = site.get_straight( j, l );
    //local id of the straightest connection in neighbour's neighbour array
          unsigned int neighIdLoc = (unsigned) site.get_outgoing(straight);
    //site id of the neighbour that will receive the photons
          unsigned int neigh = site.get_neighId( straight );
    //photons to send
          double inten = N_out/double(site.get_numStraight( j ));

    //if(!sites[ neigh ].get_border() ){
      //check if neighbour is ballistic
          if( sites[ neigh ].get_ballistic() ){

        //send the photons to the neighbour in the place in the neighbour vector
        //pointing to this site
            sites[ neigh ].addRadiationDiffIn(  neighIdLoc, (float) inten );

          }else{

        //if the neighbour is not ballistic, find out in which direction 
        //bins the photons should go
            vector< unsigned int > dir_to_use;
        //this trick only works for neighbours on this proc,
        //otherwise outgoing array does not exist
            if( sites[ neigh ].get_process() == COMM_RANK ){
              for( unsigned int n=0; n<number_of_directions; n++ ){
                if( sites[ neigh ].get_outgoing(n) == neighIdLoc ){
                  dir_to_use.push_back( n );
                }
              }//for all directions
            }//if neigh on this proc

        //if no direction is associated with this neighbour, find one
            if(dir_to_use.size() == 0){

    //find closest associated direction for this Delaunay line

    //directions are stored in header
              float **refVector;
    // Assign memory to the refVector
              refVector = new float*[number_of_directions];
              for( unsigned int n=0; n<number_of_directions; n++ ){
                refVector[n]=new float[3];
              }
    //assign the first orientation to refVector
              for( unsigned int n=0; n<number_of_directions; n++ ){
                refVector[n][0] = (float) orient[orientation_index][n][0];
                refVector[n][1] = (float) orient[orientation_index][n][1];
                refVector[n][2] = (float) orient[orientation_index][n][2];
              }

              float vectorNeigh[3];
              vectorNeigh[0] = site.get_x() - sites[ neigh ].get_x();
              vectorNeigh[1] = site.get_y() - sites[ neigh ].get_y();
              vectorNeigh[2] = site.get_z() - sites[ neigh ].get_z();

    // Find the neighVec that is closest (in angle) to the refVec
              float highestInprod=-FLT_MAX;
              int highestInprodNumber=-1;
              for(unsigned int n=0; n<number_of_directions; n++) {
                float tempInprod = inproduct(vectorNeigh, refVector[n], 3);
                if( tempInprod > highestInprod ){
                  highestInprod = tempInprod;
                  highestInprodNumber = n;
                }
              }
              dir_to_use.push_back( highestInprodNumber );

    // free memory of the refVector
              for(unsigned int n=0; n<number_of_directions; n++){
                delete [] refVector[n];
              }
              delete [] refVector;

            }//if no direction associated with this linw

            for(unsigned int n=0; n<dir_to_use.size() ; n++ ){
              sites[ neigh ].addRadiationDiffIn( dir_to_use[n], (float) inten/dir_to_use.size() );
            }//for all neighbours to use
            dir_to_use.clear();

          }//if neighbour is ballistic
      // }//if neighbour not in border
        }//for all straight

      }else{
  //DCT
        if(straight_from_tess){
    //the most straightforward neighbours are taken from the Delaunay edge with which 
    //this direction bin is associated

    //in direction conserving transport, send to the same direction bin
    //in outgoing is the neighbour that is associated with current direction bin
          for( short int l=0; l<site.get_numStraight( site.get_outgoing(j) ); l++ ){
      //neighbour associated with this direction bin
            unsigned int dir = (unsigned) site.get_outgoing(j);
      //straightest neighbour for this line
            unsigned int straight = site.get_straight( dir, l );
      //site id of the neighbour that will receive the photons
            unsigned int neigh = site.get_neighId( straight );
      //photons to send
            double inten = N_out/double(site.get_numStraight( site.get_outgoing(j) ));

      //check if neighbour is ballistic
            if(sites[ neigh ].get_ballistic() ){

        //if the neighbour is ballistic, find out which 
        //Delaunay line is associated with this direction
              bool found = 0;
              for( unsigned int n=0; !found && n<sites[ neigh ].get_numNeigh(); n++ ){
                if( sites[neigh].get_neighId(n) == site.get_site_id() ){
                  found = 1;
                  sites[neigh].addRadiationDiffIn( n, (float) inten );
                } 
              }
              if(!found){
                cerr << " Error in photon transport: neighbour not found! " << endl;
              }
            }else{

        //if not ballistic,
        //send the photons to the neighbour in the same direction bin
              sites[ neigh ].addRadiationDiffIn( j , (float) inten );
            }//if neighbour is ballistic
          }//for all straight

        }else{
    //calculate the straightest neigbours with respect to the current direction bin

    //arrays to hold the neighbour vector 
    //and the maximum inner product
          double max[3], vector1[3];
    //array to hold the indices of the most straight forward neighbours
          unsigned int index[3];
          double cosStraightAngle = cos( straightAngle );

    //initialise the max and index arrays
          for( short int q=0; q<dimension; q++ ){
            max[q] = -FLT_MAX; 
            index[q] = -1;
          }//for q

    //fill vector1 with jth direction bin vector 
          vector1[0] = orient[orientation_index][j][0];//(double) it->get_x() - sites[ neighId ].get_x();
          vector1[1] = orient[orientation_index][j][1];//(double) it->get_y() - sites[ neighId ].get_y();
          vector1[2] = orient[orientation_index][j][2];//(double) it->get_z() - sites[ neighId ].get_z();

    //loop over all neighbours except the jth neighbour to calculate 
    //the most straight forward ones with respect to neighbour j
          for( unsigned int k=0; k<site.get_numNeigh(); k++ ){

      //site id of this neighbour
            unsigned int neighId = site.get_neighId(k);

      //arrays to hold the neighbour vector of this neighbour
            double vector2[3];

      //fill vector2 with neighbour vector
            vector2[0] = (double) sites[ neighId ].get_x() - site.get_x();
            vector2[1] = (double) sites[ neighId ].get_y() - site.get_y();
            vector2[2] = (double) sites[ neighId ].get_z() - site.get_z();

      //calculate inner product of both vectors
            double inprod = inproduct(vector1, vector2, dimension);

      //store d largest inner products in max array
      //store the according neighbour indices in index
            if( inprod > max[0] ) {
              max[2] = max[1];
              max[1] = max[0]; 
              max[0] = inprod;
              index[2] = index[1];
              index[1] = index[0];
              index[0] = k;
            } else if ( inprod > max[1]) {
              max[2] = max[1];
              max[1] = inprod;
              index[2] = index[1];
              index[1] = k;
        // Third most straigthforward only needed in 3D case
            } else if ( dimension == 3 && inprod > max[2] ) { 
              max[2] = inprod;
              index[2] = k;
            }

          }//for k 


          short int num_straight = 1;
    // Second and third depending on angle
          for(int l=1; l < dimension; l++) {
      //add the neighbour if the maximum inner product
      //is larger than the cosine of the largest angle
            if( max[l] > cosStraightAngle ){
              num_straight++;
            }	
          }	      

    //photons to send
          double inten = N_out/double(num_straight);

          for(short int l =0;l<num_straight;l++){
      //neighbour id of straightest neighbour
            unsigned int neigh = site.get_neighId( index[l] );
      //check if neighbour is ballistic
            if(sites[ neigh ].get_ballistic() ){
        //if the neighbour is ballistic, find out which 
        //Delaunay line is associated with this direction
              bool found = 0;
              for( unsigned int n=0; !found && n<sites[ neigh ].get_numNeigh(); n++ ){
                if( sites[neigh].get_neighId(n) == site.get_site_id() ){
                  found = 1;
                  sites[neigh].addRadiationDiffIn( n, (float) inten );
                } 
              }
              if(!found){
                cerr << " Error in photon transport: neighbour not found! " << endl;
              }
            }else{

        //if not ballistic,
        //send the photons to the neighbour in the same direction bin
              sites[ neigh ].addRadiationDiffIn( j , (float) inten );
            }//if neighbour is ballistic

          }//for all straight

        }//if straight from tesselation

      }//if this site ballistic
    } //if intensityOut
  } //for all neighbours

  //loop over all neighbours/directions
  for( unsigned int j=0; j < numPixels; j++) {
    //intensity has been send away, so delete it
    site.set_intensityOut(j,0.0);
  }

} 



void SimpleX::radiation_transport( const unsigned int& run ){

  double t0 = MPI::Wtime();

  int dummy=0;
  double numRecombs = 0.0, allRecombs = 0.0; 
  double totalNumRecombs;

  //make sure no intensity is present at the start of the simulation
  if( run==0 ){
    for( SITE_ITERATOR it=sites.begin(); it!=sites.end(); it++ ){

      //in case of ballistic transport, intensity has size of number of neighbours;
      //in case of direction conserving transport, intensity has 
      //the size of the tesselation of the unit sphere
      numPixels = ( it->get_ballistic() ) ? it->get_numNeigh() : number_of_directions;

      for( unsigned int j=0; j<numPixels; j++ ) { 
        it->set_intensityOut( j, 0.0 );
        it->set_intensityIn( j, 0.0 );
      }
    }
  }//if run

  if(run==0 && fillChoice == AUTOMATIC){
    //set_source_in_centre();
  }  

  if( COMM_RANK == 0)
    cerr << endl; 

  for( unsigned int times=0; times<numSweeps; times++ ) { 

    if( COMM_RANK == 0 ){
      if( (int)floor( ( 100.0*times)/numSweeps ) > dummy ) { 
        dummy = (int) floor((100.0*times)/numSweeps); 
        cerr << "\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\t" << dummy+1 << "% completed..." << flush; 
      } 
    }

    //Do all sources
    //double total_inten = 0.0;
    for( SITE_ITERATOR it=sites.begin(); it!=sites.end(); it++ ){
      if( it->get_flux() > 0.0){
        source_transport( *it );
      }//if flux
    }//for all sites
    
    send_intensities();

    //Compute the recombinations
    double recombinations = 0.0;
    if( recombination && !photon_conservation ) { 
      //old way of doing recombinations
      double rec =  recombine();
      numRecombs += rec;
      MPI::COMM_WORLD.Allreduce(&numRecombs,&totalNumRecombs,1,MPI::DOUBLE,MPI::SUM);
      MPI::COMM_WORLD.Allreduce(&rec,&recombinations,1,MPI::DOUBLE,MPI::SUM);
      allRecombs += totalNumRecombs;
    } 


    //redistribute the photons
    for( SITE_ITERATOR it=sites.begin(); it!=sites.end(); it++ ){
      if( !it->get_border() && it->get_process() == COMM_RANK ) {

        //solve the rate equation to determine ionisations and recombinations
        double N_out = solve_rate_equation( *it );
        //if there is intensity to send, do so
        if(N_out > 0.0){
          ballistic_transport( *it, N_out );
        }
      }//if not in border   
    }//for all sites


    for( SITE_ITERATOR it=sites.begin(); it!=sites.end(); it++ ){

      //in case of ballistic transport, intensity has size of number of neighbours;
      //in case of direction conserving transport, intensity has 
      //the size of the tesselation of the unit sphere
      numPixels = ( it->get_ballistic() ) ? it->get_numNeigh() : number_of_directions;

      for( unsigned int j=0; j<numPixels; j++ ) { 
        //double inten = (double) it->get_intensityIn(j) + (double) it->get_intensityOut(j);
        it->set_intensityOut( j, it->get_intensityIn(j) );
        it->set_intensityIn( j, 0.0 );
      }

    }//for all sites

    //if direction conserving transport is possible, rotate solid angles or
    //do virtual vertex movement if set
    if( dirConsTransport || combinedTransport ){
      rotate_solid_angles();
    }//if direction conserving transport

  } //for all sweeps

  double t1 = MPI::Wtime();
  if( COMM_RANK == 0 )
    simpleXlog << endl << "  Calculating radiation transport took " << t1-t0 << " seconds" << endl << endl;

}


/********************************************************************************/
/*                             Output Routines                                  */
/********************************************************************************/

//generate output files
void SimpleX::generate_output( const unsigned int& run ){

  int output=numRuns;

  if(numOutputs > 0) 
    output = (int) numRuns/numOutputs;

  if( run==0 || (run+1)%output == 0 ){

    double t0 = MPI::Wtime();

    if( COMM_RANK == 0){
      cerr << endl << endl << " (" << COMM_RANK << ") Generating output" << endl;
    }

    float time = (run+1)*simTime*time_conversion/(numRuns*secondsPerMyr);

    char fileName[100][256];
    const char* outputBaseName1 = "SimpleX_data";

    //write output per processor to avoid disk waits
    for(unsigned int i=0; i<COMM_SIZE; i++){
      if(COMM_RANK == i){

#ifdef HDF5_PARALLEL
        sprintf(fileName[0], "%s_%f.hdf5", outputBaseName1, time);
#else
        sprintf(fileName[0], "%s_%f_%d.hdf5", outputBaseName1, time, COMM_RANK);
#endif

        write_hdf5_output(fileName[0], run);
      }
      MPI::COMM_WORLD.Barrier();
    }

    double t1 = MPI::Wtime();
    if( COMM_RANK == 0 ){
      simpleXlog << endl << "  Generating output took " << t1-t0 << " seconds" << endl;
    }

    if( COMM_RANK == 0 ){

 
      simpleXlog << endl << endl << " END OF RUN " << run+1 << endl << endl;
      simpleXlog << "  Output written to: " << fileName[0] << endl;

    }

  }//if output in this run

}


//Write hdf5 output. Difference for serial version is the proc dependence 
//of the offset
void SimpleX::write_hdf5_output(char *name, const unsigned int& run){


  //determine sites local to proc
  //using size_type is safer than using (unsigned) ints
  vector< unsigned long long int > local_sites;
  for( unsigned long long int i=0; i<sites.size(); i++ ){
    if( sites[i].get_process() == COMM_RANK ){
      local_sites.push_back( i );
    }
  }

  //simplices are local always, determine total number of simplices
  //across procs
  //size_type can't be send by mpi, so convert to double and then convert back

#ifdef HDF5_PARALLEL
  //get total number of simplices
  double numSimplLocal = (double) simplices.size();
  double temp;
  MPI::COMM_WORLD.Allreduce(&numSimplLocal, &temp, 1, MPI::DOUBLE, MPI::SUM);
  unsigned long long int numSimpl = (unsigned long long int) temp;

  //let every proc know the number of vertices and simplices on other procs
  vector< unsigned int > number_of_vertices( COMM_SIZE+1, 0 );
  vector< unsigned int > number_of_simplices( COMM_SIZE+1, 0 );
  vector< double > temp_number_of_vertices( COMM_SIZE+1, 0 );
  vector< double > temp_number_of_simplices( COMM_SIZE+1, 0 );
  vector< double > IDX( COMM_SIZE+1, 0 );
  vector< double > IDXX( COMM_SIZE+1, 0 );

  IDX[ COMM_RANK ] = (double) local_sites.size();
  IDXX[ COMM_RANK ] = (double) simplices.size();

  MPI::COMM_WORLD.Allreduce(&IDX[0], &temp_number_of_vertices[0], COMM_SIZE+1, MPI::DOUBLE, MPI::SUM);
  MPI::COMM_WORLD.Allreduce(&IDXX[0], &temp_number_of_simplices[0], COMM_SIZE+1, MPI::DOUBLE, MPI::SUM);

  IDX.clear();
  IDXX.clear();

  unsigned int sum = 0;
  unsigned int buffer = (unsigned int) temp_number_of_vertices[0];
  temp_number_of_vertices[0] = 0;
  for(unsigned int i=1;i<COMM_SIZE+1;++i){
    sum = (unsigned int) temp_number_of_vertices[i-1] + buffer;
    buffer = (unsigned int) temp_number_of_vertices[i];
    number_of_vertices[i] = sum;
  }

  unsigned int sum2 = 0;
  unsigned int buffer2 = (unsigned int) temp_number_of_simplices[0];
  temp_number_of_simplices[0] = 0;
  for(unsigned int i=1;i<COMM_SIZE+1;++i){
    sum2 = (unsigned int) temp_number_of_simplices[i-1] + buffer2;
    buffer2 = (unsigned int) temp_number_of_simplices[i];
    number_of_simplices[i] = sum2;
  }

  temp_number_of_vertices.clear();
  temp_number_of_simplices.clear();

#endif

  arr_1D<double> double_arr;
  arr_1D<unsigned int> int_arr;

  int dims[2];
  unsigned long long int offset[2];

  //open file  
  h5w file(name,'n');

  // write structure of the file
  file.make_group("/Header");
  file.make_group("/Vertices");
  file.make_group("/Simplices");


  //Vertex attributes
#ifdef HDF5_PARALLEL
  dims[0] = numSites;
#else
  dims[0] = sites.size();
#endif

  dims[1] = 3;                              // only write 3D data for now...
  file.make_dataset("/Vertices/coordinates","double",2,dims);
  file.write_attr("/Vertices/coordinates","var_name","coords");
  dims[1] = 1;
  file.make_dataset("/Vertices/vertex_index","unsigned int",1,dims);
  file.write_attr("/Vertices/vertex_index","var_name","vertex_id");
  file.make_dataset("/Vertices/border","unsigned int",1,dims);
  file.write_attr("/Vertices/border","var_name","border");
  file.make_dataset("/Vertices/H_neutral_fraction","double",1,dims);
  file.write_attr("/Vertices/H_neutral_fraction","var_name","HI");
  file.make_dataset("/Vertices/number_density","double",1,dims);
  file.write_attr("/Vertices/number_density","var_name","n");
  file.make_dataset("/Vertices/volume","double",1,dims);
  file.write_attr("/Vertices/volume","var_name","vol");
  file.make_dataset("/Vertices/luminosity","double",1,dims);
  file.write_attr("/Vertices/luminosity","var_name","lum");

  //Simplex attributes
#ifdef HDF5_PARALLEL
  dims[0] = numSimpl;
#else
  dims[0] = simplices.size();
#endif

  dims[1] = 4;                              
  file.make_dataset("/Simplices/indices","unsigned int",2,dims);
  file.write_attr("/Simplices/indices","var_name","indices");
  dims[1] = 1;
  file.make_dataset("/Simplices/volume","double",1,dims);
  file.write_attr("/Simplices/volume","var_name","vol");



  //=============== write header ======================================//

#ifdef HDF5_PARALLEL
  file.write_attr("/Header","number_of_sites", (unsigned int) numSites);
  file.write_attr("/Header","number_of_simplices", (unsigned int) numSimpl);
#else
  unsigned int numSitesLoc = local_sites.size();
  file.write_attr("/Header","local_number_of_sites", numSitesLoc );
  unsigned int numSimplLoc = simplices.size();
  file.write_attr("/Header","local_number_of_simplices", numSimplLoc );
  file.write_attr("/Header","number_of_sites", (unsigned int) numSitesLoc);
  file.write_attr("/Header","number_of_simplices", (unsigned int) numSimplLoc);
#endif

  file.write_attr("/Header","number_of_sweeps", numSweeps);
  file.write_attr("/Header","box_size_in_pc", sizeBox);
  //file.write_attr("/Header","redshift", redShift);
  double this_time = (run+1)*simTime/(numRuns*secondsPerMyr);
  file.write_attr("/Header","current_time_in_Myr", this_time);
  double total_time = simTime/secondsPerMyr;
  file.write_attr("/Header","total_simulation_time_in_Myr", total_time);

  if (recombination)
    file.write_attr("/Header","recombination_on","true");
  else
    file.write_attr("/Header","recombination_on","false");

  //======================= Write Vertices ===========================//

  // writing is done one variable at a time, through the arr_1D instance of size chunk_size
  // its not equal to numSites in order to preserve memory....

  unsigned int chunk_size_vert;

  //writing in chunks isn't implemented in parallel
#ifdef HDF5_PARALLEL
  chunk_size_vert = local_sites.size();
#else
  chunk_size_vert = chunk_size;
  if (chunk_size_vert > local_sites.size())
    chunk_size_vert = local_sites.size();
#endif  

  // do the writing !!!
  offset[1] = 0;
  for( unsigned long long int i=0; i<local_sites.size(); i+=chunk_size_vert ){

#ifdef HDF5_PARALLEL
    offset[0] = (int) i + number_of_vertices[COMM_RANK];
#else
    offset[0] = i;
#endif

    if (i+chunk_size_vert >= local_sites.size()) // make sure not to write outside of data range
      dims[0] = (int) local_sites.size()-i;
    else
      dims[0] = chunk_size_vert;

      // writing coordinates
    dims[1] = 3; 

    double_arr.reinit(2,dims);
    for( int j=0; j<dims[0]; j++ ){
      unsigned long long int index = local_sites[i+j];
      double_arr(j,0) = sites[index].get_x();
      double_arr(j,1) = sites[index].get_y();
      double_arr(j,2) = sites[index].get_z();

    }
    file.write_data("/Vertices/coordinates",offset, &double_arr);


    dims[1]=1;
      //write id
    int_arr.reinit(1,dims);
    for( int j=0; j<dims[0]; j++ )
      int_arr(j) = sites[ local_sites[i+j] ].get_vertex_id();
    file.write_data("/Vertices/vertex_index",offset, &int_arr);

      //write border
    int_arr.reinit(1,dims);
    for( int j=0; j<dims[0]; j++ )
      int_arr(j) = (unsigned int) sites[ local_sites[i+j] ].get_border();
    file.write_data("/Vertices/border",offset, &int_arr);

      // writing neutral fraction      
    double_arr.reinit(1,dims);
    for( int j=0; j<dims[0]; j++ ){
      double neutr_frac = (double) sites[ local_sites[i+j] ].get_n_HI()/( (double) sites[ local_sites[i+j] ].get_n_HI() + (double) sites[ local_sites[i+j] ].get_n_HII() );
      double_arr(j) = neutr_frac;
    }
    file.write_data("/Vertices/H_neutral_fraction",offset, &double_arr);

      // writing number density of gas

    double_arr.reinit(1,dims);
    for( int j=0; j<dims[0]; j++ )
      double_arr(j) = ( (double) sites[ local_sites[i+j] ].get_n_HI() + (double) sites[ local_sites[i+j] ].get_n_HII() ) * UNIT_D;
    file.write_data("/Vertices/number_density",offset, &double_arr);

      // writing cell volumes
    double_arr.reinit(1,dims);
    for( int j=0; j<dims[0]; j++ ) 
      double_arr(j) = (double) sites[ local_sites[i+j] ].get_volume() * UNIT_V;
    file.write_data("/Vertices/volume",offset, &double_arr);

      // writing Luminositites
    double_arr.reinit(1,dims);
    for( int j=0; j<dims[0]; j++ )
      double_arr(j) = (double) sites[ local_sites[i+j] ].get_flux() * UNIT_I;
    file.write_data("/Vertices/luminosity",offset, &double_arr);

  }//for all local sites


  local_sites.clear();

  //====================== Write Simplices =================================//

  unsigned int chunk_size_simpl = chunk_size;

  //writing in chunks isn't implemented in parallel
#ifdef HDF5_PARALLEL
  chunk_size_simpl = simplices.size();
#else
  chunk_size_simpl = chunk_size;
  if (chunk_size_simpl > simplices.size())
    chunk_size_simpl = simplices.size();
#endif


  // do the writing !!!
  offset[1] = 0;
  for (unsigned long long int i=0; i<simplices.size(); i+=chunk_size_simpl){

#ifdef HDF5_PARALLEL
    offset[0] = (int) i + number_of_simplices[COMM_RANK];
#else
    offset[0] = (int) i;
#endif

    if ( (i + chunk_size_simpl) >= simplices.size() ) // make sure not to write outside of data range
      dims[0] = (int) simplices.size() - i;
    else
      dims[0] = chunk_size_simpl;

      // writing indices
    dims[1] = 4; 

    int_arr.reinit(2,dims);
    for( int j=0; j<dims[0]; j++ ){
      int_arr(j,0) = simplices[i+j].get_id1();//sites[simplices[i+j].get_id1()].get_vertex_id();
      int_arr(j,1) = simplices[i+j].get_id2();//sites[simplices[i+j].get_id2()].get_vertex_id();
      int_arr(j,2) = simplices[i+j].get_id3();//sites[simplices[i+j].get_id3()].get_vertex_id();
      int_arr(j,3) = simplices[i+j].get_id4();//sites[simplices[i+j].get_id4()].get_vertex_id();
    }
    file.write_data("/Simplices/indices",offset, &int_arr);

    dims[1]=1;

      // writing simplex volumes
    double_arr.reinit(1,dims);
    for( int j=0; j<dims[0]; j++ ) 
      double_arr(j) = simplices[i+j].get_volume();

    file.write_data("/Simplices/volume",offset, &double_arr);

  }

  dims[0]=0;
  dims[1]=0;

  file.close();

#ifdef HDF5_PARALLEL
  number_of_vertices.clear();
  number_of_simplices.clear();
#endif

}

/********************************************************************************/
/*                             Generic Routines                                  */
/********************************************************************************/

//clear all temporary arrays
void SimpleX::clear_temporary(){

  //delete temporary
  for( SITE_ITERATOR it=sites.begin(); it!=sites.end(); it++ ){

    it->delete_straight();
    it->delete_neighId();
    it->delete_intensityIn();
    it->delete_intensityOut();
    it->delete_outgoing();

  }//for all sites

  simplices.clear();
  send_list.clear();

}

