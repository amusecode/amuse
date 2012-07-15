/*************************************************************************
description:
This file contains the implementation of the SimpleX class that does
the radiative transfer calculations

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

#include "SimpleX.h"
#include "Map16.h"
#include "Map21.h"
#include "Map32.h"
#include "Map42.h"
#include "Map64.h"
#include "Map84.h"

#include <algorithm>
#include <valarray>

#ifndef M_PI
#define M_PI           3.14159265358979323846
#endif

using namespace std;

// constructor
SimpleX::SimpleX(const string& output_path, const string& data_path){

  //rank of this processor and total number of processors
  COMM_RANK = MPI::COMM_WORLD.Get_rank();    
  COMM_SIZE = MPI::COMM_WORLD.Get_size();    

  //dimension of the simulation (only 3D is currently supported)
  dimension = 3;

  //conversion factors between numerical and physical units
  UNIT_T = 0.0;
  UNIT_L = 0.0;

  //total volume, total number of atoms and maximal resolution of/in simulation domain
  totalVolume = 0.0;
  totalAtoms = 0.0;
  localMaxRes = 0;
  time_conversion = 1.0;

  //euler angles
  euler_phi = 0.0;
  euler_theta = 0.0;
  euler_psi = 0.0;

  //number of pixels of HEAL_PIX sphere
  numPixels = 0;

  //orientation of the unit sphere tesselation
  orientation_index = 0;

  // Random number generator
  ran = gsl_rng_alloc (gsl_rng_taus);

  // Method for RT
  ballisticTransport = 0;
  dirConsTransport = 0;
  combinedTransport = 0;
  diffuseTransport = 0;
  
  // Value for effective temperature of source (should be read from HDF5)
  sourceTeff = 1e5;

  //output to log file
  if (output_path.size() < 1) 
    simpleXlog.open( "SimpleX.log" );
  else
    simpleXlog.open( (output_path + "/SimpleX.log").c_str() );
  
  dataPath = data_path;


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
  unsigned int total_numSweeps = (unsigned int) floor(simTime/UNIT_T);

  //set the number of run
  numRuns = numOutputs;
  
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
  read_vertex_list();

  //only the master proc creates the boundary around the domain
  if(COMM_RANK == 0){

    //create boundary
    create_boundary();
    
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
void SimpleX::read_parameters( char* initFileName ){

  //create the keyvalue variable
  string strFilename = initFileName;
  ConfigFile KeyValueFile(strFilename);

  //strings needed for reading in strings
  int freq_choice;

  //read in input files
  KeyValueFile.readInto(inputFileName, "inputFile", string("./Input.hdf5") );
  //random seed to be used
  KeyValueFile.readInto(randomSeed, "randomSeed");
  //Number of points in boundary
  KeyValueFile.readInto(borderSites, "borderPoints", (unsigned) 20000 );
  //Buffer around unity domain in which extra points are placed
  KeyValueFile.readInto(borderBox, "padding", (float) 0.1 );
  //hilbert order to determine domain decomposition
  KeyValueFile.readInto(hilbert_order, "hilbert_order", (unsigned int) 1 );
  //Buffer around subbox in which extra points are triangulated
  KeyValueFile.readInto(padding_subbox, "padding_subbox", (float) 0.15 );
  //Physical timestep per sweep
  KeyValueFile.readInto( UNIT_T, "time_step", (double) 0.05 );
  //Simulation time
  KeyValueFile.readInto(simTime, "simTime", (double) 500.0 );
  //number of frequencies
  KeyValueFile.readInto(numFreq, "number_of_frequencies", (short int) 1 );
  //frequency bin spacing
  KeyValueFile.readInto(freq_choice, "frequency_spacing", 0 );

  switch (freq_choice){

  case 0:
    freq_spacing = NO_WEIGHTS;
    break;
  case 1:
    freq_spacing = LOG_WEIGHTS;
    break;
  case 2:
    freq_spacing = ENERGY_WEIGHTS;
    break;
  case 3:
    freq_spacing = IONISATION_WEIGHTS;
    break;
  default: 
    freq_spacing = IONISATION_WEIGHTS;

  }

  //Temperature of ionized gas
  KeyValueFile.readInto(gasIonTemp, "gasIonTemp", (double) 1.0e4 );
  //Temperature of neutral gas
  KeyValueFile.readInto(gasNeutralTemp, "gasNeutralTemp", (double) 10.0 );
  //Units of source
  KeyValueFile.readInto(UNIT_I, "sourceUnits", (double) 1.0e48 );
  //Include recombination?
  KeyValueFile.readInto( recombination, "recombination", (bool) 1 );
  //Only ballistic transport, DCT or combined?
  KeyValueFile.readInto( RTmethod, "RTmethod", (short) 2 );
  
  switch(RTmethod){
  case 0:
    ballisticTransport = 1;
    break;
  case 1:
    dirConsTransport = 1;
    break;
  case 2:
    combinedTransport = 1;
    break;
  case 3:
    diffuseTransport = 1;
    break;
  default:
    if(COMM_RANK == 0)
      cerr << " (" << COMM_RANK << ") WARNING: No correct input for the RT method chosen. Defaulting to combined transport." << endl;
    combinedTransport = 1;
    break;
  }

  //Number of directions
  KeyValueFile.readInto( number_of_directions, "number_of_directions", (unsigned) 42 );
  //Number of orientations
  KeyValueFile.readInto( number_of_orientations, "number_of_orientations", (unsigned) 100 );
  //Number of outputs
  KeyValueFile.readInto(numOutputs, "outputs", (unsigned) 50 );
  //Output IFront for source in centre?
  KeyValueFile.readInto(give_IFront, "IFront", (bool) 0 );
  //Chunk size for hdf5 writing
  KeyValueFile.readInto( chunk_size, "chunkSize", (unsigned) 100000 );
  //Maximum number of messages to send in MPI routines
  KeyValueFile.readInto( max_msg_to_send, "max_msg_to_send", (unsigned) 100000 );
  //Number of reference pixels in the case HEAL_PIX is used
  KeyValueFile.readInto( num_ref_pix_HP, "number_of_reference_pixels_HealPix", 5);
  //Maximal angle between straightforward direction and Delaunay lines (second and third)
  KeyValueFile.readInto( straightAngle, "straight_angle", (float) 90.0 );
  straightAngle *= M_PI/180.0;
  //Calculate the straightest directions from the tessellation or not? 
  KeyValueFile.readInto(straight_from_tess, "straight_from_tess", (bool) 1 );
  //Vertex deletion criterium
  KeyValueFile.readInto(switchTau, "maxRatioDiff");
  //Use temporal photon conservation?
  KeyValueFile.readInto( photon_conservation, "temporal_photon_conservation", (bool) 1 );
  //Fraction of characteristic time scale at which subcycling is done
  KeyValueFile.readInto( subcycle_frac, "subcycle_fraction", (double) 0.05 );
  //include collisional ionisations?
  KeyValueFile.readInto( coll_ion, "coll_ion", (bool) 1 );
  
  //include heating and cooling?
  KeyValueFile.readInto( heat_cool, "heat_cool", (bool) 0 );
  if(COMM_RANK == 0)
    if(heat_cool&& (!blackBody)){
      cerr << "ERROR: Heating needs spectral information (e.g. bbSpectrum = 1)." << endl;
      MPI::COMM_WORLD.Abort(-1);
    }
  //include metal line cooling?
  KeyValueFile.readInto( metal_cooling, "metal_cooling", (bool) 0 );
  if(COMM_RANK == 0)
    if(metal_cooling&& (!heat_cool)){
      cerr << "ERROR: Metal cooling needs heat_cool as well." << endl;
      MPI::COMM_WORLD.Abort(-1);
    }
  

  //write all information to logfile
  if(COMM_RANK == 0){
    simpleXlog << endl <<" *****  This is SimpleX version 2.5  ***** " << endl << endl;
    simpleXlog << "  Number of processors: " << COMM_SIZE << endl;
    simpleXlog << "  Input file: " << inputFileName << endl;
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

    if(recombination){
      simpleXlog << "  Recombination is included" << endl;
    }

    if(ballisticTransport){
      simpleXlog << endl << "  Mode of transport: ballistic transport " << endl << endl;
    }
    if(dirConsTransport){
      simpleXlog << endl << "  Mode of transport: direction conserving transport " << endl << endl;
    }
    if(combinedTransport){
      simpleXlog << endl << "  Mode of transport: combined transport " << endl << endl;
    }
    if(diffuseTransport){
      simpleXlog << endl << "  Mode of transport: diffuse transport " << endl << endl;
    }

    if(dirConsTransport || combinedTransport){
      simpleXlog << "  Number of direction bins used        : " << number_of_directions << endl;
      simpleXlog << "  Number of orientations in header file: " << number_of_orientations << endl << endl;
    }


#ifdef HEALPIX
    simpleXlog << "  Healpix sphere with " << num_ref_pix_HP 
	       << " pixels is used for quick referencing in compute_solid_angles()" << endl;
#endif

    simpleXlog << "  Switch between ballistic and dirCons transport: " << switchTau << endl;
    simpleXlog << "  Maximal angle between forward direction and real Delaunay direction: " 
	       << straightAngle << " degrees. " << endl;

  }

  if(give_IFront){
    //output to IFront file
    IFront_output.open( "IFront.dat" );
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
}


/****  Read in vertices from hdf5 file  ****/
void SimpleX::read_vertex_list(){

  //open hdf5 file
  char fileName[200];
  sprintf(fileName, "%s", inputFileName.c_str());
  h5w file(fileName,'o');        // open file for reading

  //read in total number of sites
  unsigned int temp;
  file.read_attr("/Header","number_of_sites", &temp); 
  numSites = (unsigned long int) temp;  

  //read in conversion factors for dimensionless units 
  file.read_attr("/Header","UNIT_D",&UNIT_D);
  file.read_attr("/Header","UNIT_I",&UNIT_I);
  file.read_attr("/Header","UNIT_L",&UNIT_L);
  file.read_attr("/Header","UNIT_M",&UNIT_M);
  UNIT_V = pow( UNIT_L, 3.0 );


  //read in the size of the simulation box
  file.read_attr("/Header","box_size", &sizeBox);
  short int dummy;
  file.read_attr("/Header","sourceType", &dummy);
  blackBody = (bool)dummy;
  file.read_attr("/Header","sourceTemperature", &sourceTeff);

  //read if clumping factor is included
  // short int read_clumping;
  // file.read_attr("/Header","clumping", &read_clumping);
  // cerr << read_clumping << endl;

  //create vector containing the vertices, and the temporary lists containing fluxes and masses
  //the last two are necessary because they are needed in check_undersampled()
  vertices.resize(numSites);
  temp_n_HI_list.resize(numSites);
  temp_flux_list.resize(numSites);
  temp_n_HII_list.resize(numSites);
  temp_u_list.resize(numSites);
  temp_dudt_list.resize(numSites);
  temp_clumping_list.resize(numSites);
  temp_metallicity_list.resize(numSites);
  
  //structures needed for reading in the values from hdf5
  arr_1D<float> double_arr;
  arr_1D<unsigned int> int_arr;

  //arrays holding the dimensions of the data and the offset in case the data is read in in chunks
  unsigned long long int  dims[2], offset[2];

  unsigned int chunk_size_read = chunk_size;  
  if (chunk_size_read > numSites){
    chunk_size_read = numSites;
  }

  // do the reading !!!
  offset[1] = 0;
  for( unsigned long int i=0; i<numSites; i+=chunk_size_read ){

    offset[0] = i;
    if (i+chunk_size_read >= numSites) // make sure not to write outside of data range
      dims[0] = numSites - i;
    else
      dims[0] = chunk_size_read;

    // coordinates
    dims[1] = 3; 
    double_arr.reinit(2,dims);
    file.read_data("/Vertices/coordinates",offset, &double_arr);
    for(unsigned int j=0; j<dims[0]; j++){
      vertices[ j + i ].set_x( double_arr(j,0) ); // coordinates
      vertices[ j + i ].set_y( double_arr(j,1) ); 
      vertices[ j + i ].set_z( double_arr(j,2) );    
    }
    
    dims[1]=1;
    //number density
    double_arr.reinit(1,dims);
    file.read_data("Vertices/n_HI",offset, &double_arr);
    for(unsigned int j=0; j<dims[0]; j++ ){
      temp_n_HI_list[ j + i ] = double_arr(j);
    }

    //ionised fraction
    double_arr.reinit(1,dims);
    file.read_data("Vertices/n_HII",offset, &double_arr);
    for(unsigned int j=0; j<dims[0]; j++ ){
      temp_n_HII_list[ j + i ] = double_arr(j);
    }

    //flux
    double_arr.reinit(1,dims);
    file.read_data("Vertices/sourceNIon",offset, &double_arr);
    for(unsigned int j=0; j<dims[0]; j++ ){
      temp_flux_list[ j + i ] = double_arr(j);
    }

    //temperature
    double_arr.reinit(1,dims);
    file.read_data("Vertices/temperature",offset, &double_arr);
    for(unsigned int j=0; j<dims[0]; j++ ){
      temp_u_list[ j + i ] = T_to_u(double_arr(j),
       ( temp_n_HI_list[ j + i ] + temp_n_HII_list[ j + i ] )/
         (temp_n_HI_list[ j + i ] + 2*temp_n_HII_list[ j + i ]));
      temp_dudt_list[ j + i ] = 0.0;
    }

    //clumping
    double_arr.reinit(1,dims);
    file.read_data("Vertices/clumpingFactor",offset, &double_arr);
    for(unsigned int j=0; j<dims[0]; j++ ){
      temp_clumping_list[ j + i ] = double_arr(j);
    }
  }//for all sites to read

  //set the vertex id
  //perhaps in the future this should be read in as well?
  unsigned long int i=0;
  for( VERTEX_ITERATOR it=vertices.begin(); it!=vertices.end(); it++, i++ ){
    it->set_vertex_id( i );
  }

  origNumSites = vertices.size();
  vertex_id_max = vertices.size() + borderSites;

  file.close();

  if(COMM_RANK == 0){
    simpleXlog << "  Read " << numSites << " sites from file " << endl;
    simpleXlog << "  Size of the simulation domain: " << sizeBox << " pc" << endl;
    cerr << " (" <<COMM_RANK <<") Read " << numSites << " sites from file " << endl;
    cerr << " (" <<COMM_RANK <<") Size of the simulation domain: " << sizeBox << " pc" << endl;
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
  for( unsigned long int i=0; i<borderSites; i++) {
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
    tempVert.set_vertex_id( (unsigned long int ) origNumSites+i );
    tempVert.set_border(1);

    vertices.push_back( tempVert );

  }//for all border sites

  if( COMM_RANK == 0 )
    simpleXlog << "  Boundary points placed " << endl;

}



/****  Create octree of vertices  ****/
//size of octree is 4 * #subboxes, to ensure fast localisation of vertices
//in subbox and subbox boundaries
void SimpleX::create_vertex_tree(){

  //make sure tree is empty
  vert_tree.delete_octree();

  //number of subboxes in one dimension
  unsigned int base = 2;
  int power = (int) hilbert_order;
  unsigned int hilbert_resolution = pow( base, power );
  //tree size is 4*number of subboxes (in one dimension) to make search for boundary vertices most efficient
  unsigned int tree_size = 4*hilbert_resolution;

  //initilise the tree with this computed size
  double temp_borderBox = (double) borderBox;
  vert_tree.init_octree( tree_size, temp_borderBox );

  //put the vertices in the tree
  for( unsigned long int i=0; i<vertices.size(); i++ ){
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
  unsigned int numSites_per_proc = unsigned( ceil( (double) vertices.size()/(COMM_SIZE) ));

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
    unsigned long r = (unsigned long) i;

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
    unsigned long r = (unsigned long) i;

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
  unsigned long this_subbox = (unsigned long) start_number;

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
	      coord_CC[0] = (unsigned long) x_hilbert;
	      coord_CC[1] = (unsigned long) y_hilbert;
	      coord_CC[2] = (unsigned long) z_hilbert;

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
		    coord_vert[0] = (unsigned long) x_vert;
		    coord_vert[1] = (unsigned long) y_vert;
		    coord_vert[2] = (unsigned long) z_vert;

		    //calculate the hilbert number of the subbox the vertex is in
		    unsigned long r_vert = hilbert_c2i( dimension, hilbert_order, coord_vert );

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
		x_min = xCC - 1.001*radiusCC;
		correct = 0;
	      }
	      if( (xCC + radiusCC) > x_max ){
		x_max = xCC + 1.001*radiusCC;
		correct = 0;
	      }
	      if( (yCC - radiusCC) < y_min ){
		y_min = yCC - 1.001*radiusCC;
		correct = 0;
	      }
	      if( (yCC + radiusCC) > y_max ){
		y_max = yCC + 1.001*radiusCC;
		correct = 0;
	      }
	      if( (zCC - radiusCC) < z_min ){
		z_min = zCC - 1.001*radiusCC;
		correct = 0;
	      }
	      if( (zCC + radiusCC) > z_max ){
		z_max = zCC + 1.001*radiusCC;
		correct = 0;
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

      // if(!correct){
      //   cerr << " (" << COMM_RANK << ") Subbox " << i << " boundaries too small, retriangulating with subbox: (" << x_min << "," << x_max << ") (" 
      //     << y_min << "," << y_max << ") (" << z_min << "," << z_max << ")" << endl;
      // }

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
//since we count every site_index once. This is possible since duplicate
//sites outside the simulation domain point with their site_index to 
//the site that is in the simulation domain
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
    //cerr << " (" << COMM_RANK << ") Number of sites in triangulation is " << numSites << endl;
    simpleXlog << "  Final triangulation contains " << numSites << " sites " << endl;
  }

  //   if( COMM_RANK == 0 ){
  //     cerr << " (" << COMM_RANK << ") Number of simplices in triangulation is " << simplices.size() << endl;
  //   }

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


/****  Assign id's to sites  ****/
void SimpleX::assign_site_ids(){

  //set correct site id
  unsigned long int i=0;
  for( SITE_ITERATOR it=sites.begin(); it!=sites.end(); it++, i++ ){
    it->set_site_id(i);
  }

  //sort the sites
  sort( sites.begin(), sites.end(), compare_vertex_id_site );

  //vector to keep track of the indices
  vector<unsigned long int> indices( sites.size(), 0 );

  //store the new place of the sites in the indices array at the place of their 
  //former index
  for( unsigned long int i=0; i<sites.size(); i++ ){
    indices[ sites[i].get_site_id() ] = i;
  }   

  //loop over all simplices
  for( SIMPL_ITERATOR it=simplices.begin(); it!=simplices.end(); it++ ){

    //former place of the sites in sites array
    unsigned long int id1 = it->get_id1();
    unsigned long int id2 = it->get_id2();
    unsigned long int id3 = it->get_id3();
    unsigned long int id4 = it->get_id4();

    //set the id to the new id of the sites
    it->set_id1( (unsigned long int) indices[id1] );
    it->set_id2( (unsigned long int) indices[id2] );
    it->set_id3( (unsigned long int) indices[id3] );
    it->set_id4( (unsigned long int) indices[id4] );


  }

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

  //testing, simplices no longer needed here
  //simplices.clear();
  //vector< Simpl >().swap(simplices);

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
  compute_solid_angles(1);

  //send the properties of neighbours needed for ballistic transport to other procs
  send_neighbour_properties();

  //calculate the mean Delaunay length of every site
  calculate_line_lengths();

  //create the intensity arrays, in the case of diffuse transport just one value, in the case of 
  //ballistic transport the number of HEALPIX directions
  for( SITE_ITERATOR it=sites.begin();it!=sites.end();it++ ){

    //in case of ballistic transport, intensity has size of number of neighbours;
    //in case of direction conserving transport, intensity has 
    //the size of the tesselation of the unit sphere
    numPixels = ( it->get_ballistic() ) ? it->get_numNeigh() : number_of_directions;

    //create the intensity arrays
    it->create_intensityIn( numFreq, numPixels );
    it->create_intensityOut( numFreq, numPixels );

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
  vector< vector<unsigned int> > mapping_pixels(numPixels, vector<unsigned int>(numPixels_ref,0));
  //loop over all directions of the intensities
  for( unsigned int m=0; m<numPixels; m++ ){

    //storage for the maximum inner product and 
    //associated index in refVector_ref
    vector< float > max( numPixels_ref, -1.0 );
    vector< int > index( numPixels_ref, -1 );

    //fill max and index with inner products

    //loop over all pixels of the high res healpix sphere
    for(unsigned int p=0; p<numPixels_ref; p++ ){

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

    //sort the inner products and according indices with highest first
    quickSortPerm( max, index, 0, numPixels_ref-1 );

    //loop over all pixels of healpix sphere
    for(unsigned int p=0; p<numPixels_ref; p++ ){
      //make sure the index value is in correct interval
      if( (index[p] >=0) && (index[p] < (int) numPixels_ref) ){
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


void SimpleX::calculate_line_lengths(){

  double minDist = 666.;
  //vector<double> neigh_dist;
  for( SITE_ITERATOR it=sites.begin(); it!=sites.end(); it++ ){

    double totalLength=0.0;
    for( unsigned int j=0; j<it->get_numNeigh(); j++ ){    

      double length=0.0;

      double x1 = (double) it->get_x();
      double x2 = (double) sites[ it->get_neighId(j) ].get_x();
      length += pow( x1 - x2, 2);

      double y1 = (double) it->get_y();
      double y2 = (double) sites[ it->get_neighId(j) ].get_y();
      length += pow( y1 - y2, 2);

      double z1 = (double) it->get_z();
      double z2 = (double) sites[ it->get_neighId(j) ].get_z();
      length += pow( z1 - z2, 2);

      totalLength+=sqrt(length);

    }

    it->set_neigh_dist( totalLength/it->get_numNeigh() );

    if( it->get_neigh_dist() < minDist && !it->get_border() && it->get_process() == COMM_RANK){
      minDist = it->get_neigh_dist();
    }

  }//for all vertices

  localMaxRes = (int) ceil( 1.0/minDist );

  if( COMM_RANK == 0 ){
    simpleXlog << "  Delaunay line lengths computed " << endl;
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
    coord_CC[0] = (unsigned long) x_hilbert;
    coord_CC[1] = (unsigned long) y_hilbert;
    coord_CC[2] = (unsigned long) z_hilbert;

    //determine the hilbert number of the subbox the circumcentre is in
    unsigned long long r_CC = hilbert_c2i( dimension, hilbert_order, coord_CC );

    //if this simplex belongs to this proc, we can add it
    if( dom_dec[ r_CC ] == COMM_RANK ){

      if( !sites[ one ].get_border() && !sites[ two ].get_border() &&
       	  !sites[ three ].get_border() && !sites[ four ].get_border() ){
      //if( !sites[ one ].get_border() || !sites[ two ].get_border() ||
      //	  !sites[ three ].get_border() || !sites[ four ].get_border() ){

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
    vector< unsigned  long int >().swap( intens_ids );

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
      site_intensities.insert( site_intensities.end(), numFreq*numPixels, 0.0 );

      //loop over directions
      for( unsigned int j=0; j<numPixels; j++ ){
  //loop over frequencies
        for(short int f=0; f<numFreq; f++){

    //position of this direction in array
          unsigned int pos = f + j*numFreq + site_intensities.size() - numFreq*numPixels;

    //add all intensity to site_intensities in correct place
          site_intensities[pos] += it->get_intensityIn(f,j) + it->get_intensityOut(f,j);

    //now that they are stored, set them to zero
          it->set_intensityOut(f,j,0.0);
          it->set_intensityIn(f,j,0.0);

        }//for all freq
      }//for all pixels
    }//if
  }//for all sites

}

const vector<unsigned long int> SimpleX::get_ballistic_sites_to_store(){

  vector<unsigned long int> sites_to_store;

  //loop over all sites
  for( SITE_ITERATOR it=sites.begin(); it!=sites.end(); it++){
    //only consider vertices on this proc and inside simulation domain
    if( it->get_process() == COMM_RANK && !it->get_border() && it->get_ballistic() ){
      //loop over all neighbours to see if they have intensities
      bool use = 0;
      for(short int f=0; f<numFreq;f++){
	for( unsigned int j=0; !use && j<it->get_numNeigh(); j++ ){
	  if( it->get_intensityIn(f,j) > 0.0 || it->get_intensityOut(f,j) > 0.0 ){
	    use = 1;
	  }//if there is intensity in this site
	}//for all neighbours
      }//for all freq
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
void SimpleX::store_ballistic_intensities( const vector<unsigned long int>& sites_to_store ){

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
  for( unsigned long int i=0; i<sites_to_store.size(); i++ ){

    unsigned long int site_id = sites_to_store[i];

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
    site_intensities.insert( site_intensities.end(), numFreq*numPixels, 0.0 );

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
	    for(short int f=0; f<numFreq; f++){
	      unsigned int pos = f + m*numFreq + site_intensities.size() - numFreq*numPixels;
	      double inten = ( (double) sites[ site_id ].get_intensityIn(f,j) + (double) sites[ site_id ].get_intensityOut(f,j) )/count;
	      site_intensities[pos] += (float) inten;
	    }//for all freq
	  }
	}
	//now that they are stored, set them to zero
	for(short int f=0; f<numFreq;f++){
	  sites[ site_id ].set_intensityOut(f,j,0.0);
	  sites[ site_id ].set_intensityIn(f,j,0.0);
	}
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
	for(short int f=0;f<numFreq; f++){

	  unsigned long int pos = f + highestInprodNumber*numFreq + site_intensities.size() - numFreq*numPixels;

	  double inten = ( (double) sites[ site_id ].get_intensityIn(f,j) + (double) sites[ site_id ].get_intensityOut(f,j) );
	  site_intensities[pos] += (float) inten;

	  //now that they are stored, set them to zero
	  sites[ site_id ].set_intensityOut(f,j,0.0);
	  sites[ site_id ].set_intensityIn(f,j,0.0);

	}//for all freq
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
void SimpleX::store_ballistic_intensities( const vector<unsigned long int>& sites_to_store ){


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
  for( unsigned long int i=0; i<sites_to_store.size(); i++ ){

    //place in sites vector of this site
    unsigned long int site_id = sites_to_store[i];


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
    site_intensities.insert( site_intensities.end(), numFreq*numPixels, 0.0 );

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
	    //loop over frequencies
	    for(short int f=0; f<numFreq;f++){

	      unsigned long int pos = f + m*numFreq + site_intensities.size() - numFreq*numPixels;

	      double inten = ( (double) sites[ site_id ].get_intensityIn(f,j) + (double) sites[ site_id ].get_intensityOut(f,j) )/(double) count;
	      site_intensities[pos] += (float) inten;

	      //now that they are stored, set them to zero
	      sites[ site_id ].set_intensityOut(f,j,0.0);
	      sites[ site_id ].set_intensityIn(f,j,0.0);
	      
	    }//for all freq
	  }
	}

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

	//loop over frequencies
	for(short int f=0; f<numFreq; f++){

	  unsigned long int pos = f + highestInprodNumber*numFreq + site_intensities.size() - numFreq*numPixels;

	  double inten = ( (double) sites[ site_id ].get_intensityIn(f,j) + (double) sites[ site_id ].get_intensityOut(f,j) );
	  site_intensities[pos] += (float) inten;

	  //now that they are stored, set them to zero
	  sites[ site_id ].set_intensityOut(f,j,0.0);
	  sites[ site_id ].set_intensityIn(f,j,0.0);

	}//for all freq
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
    vector< unsigned long int > mapping_sites(this_chunk_size,vertex_id_max+1);
    bool stop = 0;
    while( it != sites.end() && !stop ){ 
      //only include direction conserving sites on this proc not in border
      if( it->get_process() == COMM_RANK && !it->get_border() && !it->get_ballistic() ){  
  //associate every vertex id with its place in the array 
        if( ( (long int) it->get_vertex_id() - (long int) start_vertex_id) >= 0 && 
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

      //get the closest associate with this neighbour from the mapping 
      //stored in the header file
            unsigned int pos = maps[orientation_index][orientation_index_old][j];

      //loop over frequencies
            for(short int f=0; f<numFreq;f++){

        //get the intensity that belongs to this site
              double inten = (double) site_intensities[ i*numPixels*numFreq + j*numFreq + f ];

        //in case there is already intensity in this direction
        //obsolete with the new header file
              inten += (double) sites[site_id].get_intensityOut(f, pos);

        //assign the intensity to the site
              sites[site_id].set_intensityOut( f, pos, (float) inten );
              sites[site_id].set_intensityIn( f, pos, 0.0 );

            }//for all freqs
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
    unsigned long int start_vertex_id = q*max_msg_to_send; 
    //make sure the final chunk has correct size
    if( start_vertex_id + this_chunk_size >= vertex_id_max ){
      this_chunk_size = vertex_id_max - start_vertex_id;
    }

    //mapping from vertex id to site id
    vector< unsigned long int > mapping_sites(this_chunk_size,vertex_id_max+1);
    bool stop = 0;
    while( it != sites.end() && !stop ){ 
      //only include direction conserving sites on this proc not in border                                                                                      
      if( it->get_process() == COMM_RANK && !it->get_border() && it->get_ballistic() ){  
	//associate every vertex id with its place in the array                                                                                                 
	if( ( (long int) it->get_vertex_id() - (long int) start_vertex_id ) >= 0 && 
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
    for( unsigned long int i=0; i<intens_ids.size(); i++ ){
      //only include vertex ids in current chunk
      if(intens_ids[i] >= start_vertex_id && intens_ids[i] < (start_vertex_id+this_chunk_size) ){
	//get site id from vertex id using pre-calculated mapping
	unsigned long int site_id = mapping_sites[ intens_ids[i] - start_vertex_id ];
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

		for(short int f=0; f<numFreq; f++){

		  //give this neighbour the intensity in site_intensities
		  double inten = site_intensities[ i*numPixels*numFreq + m*numFreq + f ];
		  inten += (double) sites[site_id].get_intensityOut(f,j);

		  sites[site_id].set_intensityOut(f,j, (float) inten);
		  sites[site_id].set_intensityIn(f,j,0.0);

		}//for all freq
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
    unsigned long int start_vertex_id = q*max_msg_to_send; 
    //make sure the final chunk has correct size
    if( start_vertex_id + this_chunk_size >= vertex_id_max ){
      this_chunk_size = vertex_id_max - start_vertex_id;
    }

    //mapping from vertex id to site id
    vector< unsigned long int > mapping_sites(this_chunk_size,vertex_id_max+1);
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
    for( unsigned long int i=0; i<intens_ids.size(); i++ ){
      //only include vertex ids in current chunk
      if(intens_ids[i] >= start_vertex_id && intens_ids[i] < (start_vertex_id+this_chunk_size) ){
	//get site id from vertex id using pre-calculated mapping
	unsigned long int site_id = mapping_sites[ intens_ids[i] - start_vertex_id ];
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

	    //loop over frequencies
	    for(short int f=0; f<numFreq; f++){

	      //give intensity of this pixel to Delaunay line with highest inner product
	      double inten = site_intensities[ i*numPixels*numFreq + m*numFreq + f ];

	      inten += (double) sites[site_id].get_intensityOut(f,j);

	      sites[site_id].set_intensityOut(f,j,(float) inten);
	      sites[site_id].set_intensityIn(f,j,0.0);

	    } //for all freqs
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
  vector< unsigned long int > sites_to_store;
  //vector to hold all sites in internal boundary that have
  //been changed
  vector< unsigned long int > sites_to_send;

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

      //calculate the mean optical depth over frequencies
      double tau = 0.0;
      for(short int f=0; f<numFreq; f++){
	tau += it->get_n_HI() * UNIT_D * it->get_neigh_dist() * UNIT_L * cross_H[f] * straight_correction_factor;
      }
      tau /= numFreq;

      //if optical depth is smaller than the switch set by user, switch from
      //ballistic to direction conserving
      //perhaps do this only if radiation has travelled through this site?
      if( tau < switchTau && is_ionised ){
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
	    unsigned long int neigh = it->get_neighId(j);
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
	if( it->get_ballistic() == 0 && !it->get_source() ){

	  it->set_ballistic( 1 );

	  //if this site has incarnation on other proc, 
	  //put it in the list to send
	  bool flag = 0;
	  //loop over all neighbours to check whether this site has at least
	  //one neighbour on another proc
	  for( unsigned int j=0; !flag && j<it->get_numNeigh(); j++ ) {
	    unsigned long int neigh = it->get_neighId(j);
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
	  site_intensities.insert( site_intensities.end(), numPixels*numFreq, 0.0 );
	  //loop over directions
	  for( unsigned int j=0; j<numPixels; j++ ){
	    //loop over frequencies
	    for(short int f=0; f<numFreq; f++){

	      //position of this direction in array
	      unsigned long int pos = f + j*numFreq + site_intensities.size() - numPixels*numFreq;

	      //add all intensity to site_intensities in correct place
	      site_intensities[pos] += it->get_intensityIn(f,j) + it->get_intensityOut(f,j);

	      //now that they are stored, set them to zero
	      it->set_intensityOut(f,j,0.0);
	      it->set_intensityIn(f,j,0.0);

	    }//for all freq
	  }//for all pixels

	  //delete the intensity arrays
	  it->delete_intensityOut(numFreq);
	  it->delete_intensityIn(numFreq);

	  //create new intensity arrays with correct size
	  it->create_intensityIn( numFreq, it->get_numNeigh() );
	  it->create_intensityOut( numFreq, it->get_numNeigh() );


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
	    unsigned long int neigh = it->get_neighId(j);
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
  for(unsigned long int i=0; i<sites_to_store.size(); i++ ){

    unsigned long int site_id = sites_to_store[i];

    sites[ site_id ].delete_intensityOut(numFreq);
    sites[ site_id ].delete_intensityIn(numFreq);

    sites[ site_id ].create_intensityOut( numFreq, number_of_directions );
    sites[ site_id ].create_intensityIn( numFreq, number_of_directions );

  }

  sites_to_store.clear();
  sites_to_send.clear();

  //return the intensities of the sites that were switched from 
  //direction conserving to ballistic
  return_ballistic_intensities();

}


/****  Store the relevant properties of the sites  ****/
// stored are:  vertex_id
//              site_id
//              numberDensity
//              ion_frac  
//              localIntensity
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

  //calculate the mean optical depth over frequencies
        double tau = 0.0;
        for(short int f=0; f<numFreq; f++){
          tau += it->get_n_HI() * UNIT_D * aver_neigh_dist * UNIT_L * cross_H[f] * straight_correction_factor;
        }
        tau /= numFreq;

  //if optical depth is smaller than the switch set by user, switch from
  //ballistic to direction conserving
        if( tau < switchTau && is_ionised ){
          it->set_ballistic( 0 );
        }else{
          it->set_ballistic( 1 );
        }
  //sources are always direction conserving
        if( it->get_source() ){
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

      if( it->get_source() ){
  //loop over frequency bins to get total flux
        double total_flux = 0.0;
        for(short int f=0; f<numFreq; f++){
          total_flux += (double) it->get_flux(f);
        }
        temp.set_flux( (float) total_flux );

      }else{
        temp.set_flux( 0.0 );
      }

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

  //cerr << " (" << COMM_RANK << ") number of vertices: " << vertices.size() << endl;

  if( COMM_RANK == 0 ){
    simpleXlog << "  Created new vertex list" << endl;
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

	    //do the sending and receiving
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

      level = ngrp - 1;// Keep the level loop from duplicate operations
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

  return;
}


/**** Send neighbour information  ****/
//send the properties of the neighbours of sites that belong to this proc
//but also exist on another proc to other procs
void SimpleX::send_neighbour_properties(){

  if( COMM_RANK == 0 ){
    simpleXlog << "  Sending neighbour properties" << endl;
  }

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
	  unsigned long int neigh = it->get_neighId(j);
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
	cerr << " (" << COMM_RANK << ") Error in send_neighbour_properties: found site " << site_it->get_vertex_id() << " but not neighbour " << sit->get_neighId() << endl;
	MPI::COMM_WORLD.Abort( -1 );
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
	//loop over frequencies
	for(short int f=0; f<numFreq;f++){
	  //if the site holds intensity, it needs to be send 
	  if( sites[index].get_intensityIn( f,j ) > 0.0 || sites[index].get_intensityOut( f,j ) > 0.0  ){

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
	    //frequency bin
	    temp.set_freq_bin(f);
	    //incoming intensity
	    temp.set_intensityIn( sites[index].get_intensityIn( f,j ) );
	    //ougoing intensity
	    temp.set_intensityOut( sites[index].get_intensityOut( f,j ) );

	    //store the information in vector
	    sites_to_send.push_back(temp);

	    //photons will be on other proc, so remove them here 
	    //to conserve photons
	    sites[index].set_intensityIn( f,j, 0.0 ); 
	    sites[index].set_intensityOut( f,j, 0.0 );

	    //keep track of the number of sites to send to
	    //which proc
	    nsend_local[ sites[index].get_process() ]++;

	  }//if site holds intensity
	}//for all freq
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
	sites[ it->get_id() ].addRadiationDiffIn( it->get_freq_bin(), it->get_neighId(), it->get_intensityIn() );
      }

      //outgoing intensities
      if( it->get_intensityOut() ){
	sites[ it->get_id() ].addRadiationDiffOut( it->get_freq_bin(), it->get_neighId(), it->get_intensityOut() );
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


/****  Send the updates list of vertices to master  ****/
//vertices no longer exist only on master proc as when the grid
//was first created, so every proc has to send its own list to
//master
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
void SimpleX::send_site_ballistics( const vector< unsigned long int >& sites_to_send ){

  //send the properties of the vertex on this proc to other proc
  Send_Neigh tempNeigh;
  vector< Send_Neigh > neigh_to_send;
  //vector to hold received neighbours
  vector<Send_Neigh> neigh_recv;

  //sets if one of the procs still has sites to send
  bool total_neigh_to_send = 1;

  //loop over all sites
  unsigned long int i=0;
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
      unsigned long int site_id = sites_to_send[i];

      //keep track of the procs the site has already been sent to
      vector< bool > proc(COMM_SIZE,0);
      //loop over all neighbours to check to which 
      //procs this site has to be send
      for( unsigned int j=0; j<sites[site_id].get_numNeigh(); j++ ){

	//place in sites array of this neighbour
	unsigned long int neigh = sites[site_id].get_neighId(j);
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
	    site_it->delete_intensityOut(numFreq);
	    site_it->delete_intensityIn(numFreq);

	    //create new intensity arrays with correct size
	    site_it->create_intensityIn( numFreq, site_it->get_numNeigh() );
	    site_it->create_intensityOut( numFreq, site_it->get_numNeigh() );

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
	  site_it->delete_intensityOut(numFreq);
	  site_it->delete_intensityIn(numFreq);

	  //create new intensity arrays with correct size
	  site_it->create_intensityIn( numFreq, number_of_directions );
	  site_it->create_intensityOut( numFreq, number_of_directions );

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

  return;
}

// Send physical quantities to sites that have moved to another processor
void SimpleX::send_site_physics( ){
  
  //sort the site_properties on vertex id
  sort( site_properties.begin(), site_properties.end(), compare_vertex_id_site_update );
  
  vector<Vertex>::iterator it=vertices.begin();
  vector<Site_Update>::iterator it2=site_properties.begin();
  
  //determine whether there are still sites to be sent
  bool isSitesToSend = 1;
  
  vector<Site_Update> site_properties_recv;



  while( isSitesToSend ){

    //vector to hold site updates to send
    vector<Site_Update> sitesToSend;

    //vector to hold the sites that should be send 
    //from this proc to other proc
    vector<unsigned long long int> nSendLocal( COMM_SIZE, 0 );

    //while the maximum number of sites to send is not yet reached, 
    //search for sites that need to be send
    //note that if a site should be send to more than one proc,
    //the chunk will be slightly bigger than the maximum size,
    //but never more than COMM_SIZE    
    while( sitesToSend.size() < max_msg_to_send && it2!=site_properties.end() ){

      if( it->get_vertex_id() < it2->get_vertex_id() ){

        it++;

      }else if( it->get_vertex_id() == it2->get_vertex_id() ){

        if( it->get_process() != COMM_RANK ){
          it2->set_process( it->get_process() );
          sitesToSend.push_back( *it2 );
          nSendLocal[ it->get_process() ]++;
        }

        it++;
        it2++;

      }else if( it->get_vertex_id() > it2->get_vertex_id() ){
        cerr << " (" << COMM_RANK << ") Warning: Vertex " << it2->get_vertex_id() << " has disappeared from the computational domain " << endl;
        //MPI::COMM_WORLD.Abort(-1);
      }

    }

    //sort the site updates on process
    sort( sitesToSend.begin(), sitesToSend.end(), compare_process_site_update );
        
    //define the offset of different procs
    vector< unsigned long long int > offset( COMM_SIZE, 0 );
    for( unsigned int p=1; p<COMM_SIZE; p++ ){
      offset[p] = offset[p - 1] + nSendLocal[p - 1];
    }

    //gather all the local numbers of sites to send to all procs
    vector< unsigned long long int > nSend( COMM_SIZE * COMM_SIZE, 0 );
    MPI::COMM_WORLD.Allgather(&nSendLocal[0], COMM_SIZE, MPI::UNSIGNED_LONG_LONG, &nSend[0], COMM_SIZE, MPI::UNSIGNED_LONG_LONG );

    //calculate number of sites to receive
    unsigned long long int nSitesToReceive = 0;
    for( unsigned int p=0; p<COMM_SIZE; p++){
      if( p != COMM_RANK ){
        nSitesToReceive += nSend[ p*COMM_SIZE + COMM_RANK ];
      }
    }

    //vector to hold received sites
    vector< Site_Update > sitePropertiesReceived( nSitesToReceive );

    int PTask,recvTask;
    int NTask = static_cast<int>( COMM_SIZE );
    int ThisTask = static_cast<int>( COMM_RANK );

    // calculate PTask,  the smallest integer that satisfies COMM_SIZE<=2^PTask
    for(PTask = 0; NTask > (1 << PTask); PTask++) ;

    //vector to hold buffer with number of sites
    //already send to that proc
    vector< unsigned long long int > nBuffer( COMM_SIZE, 0 );

      //loop over all tasks
    for( int ngrp = 1; ngrp < (1 << PTask); ngrp++ ){

      recvTask = ThisTask ^ ngrp;

      if(recvTask < NTask){

        //check if there's sites to receive or to send
        if(nSend[ ThisTask * NTask + recvTask ] > 0 || nSend[ recvTask * NTask + ThisTask ] > 0){

          //do the sending
          MPI::COMM_WORLD.Sendrecv(&sitesToSend[ offset[recvTask] ], nSendLocal[recvTask] * sizeof(Site_Update),
            MPI::BYTE, recvTask, 0, 
            &sitePropertiesReceived[ nBuffer[ThisTask] ], nSend[recvTask * NTask + ThisTask] * sizeof(Site_Update),
            MPI::BYTE, recvTask, 0);

        }//if there's sites to send or receive

      }//if receiving task is smaller than COMM_SIZE

      //update the buffer size
      for( int p = 0; p < NTask; p++){
        if( (p ^ ngrp) < NTask ){
          nBuffer[p] += nSend[(p ^ ngrp) * NTask + p];
        }
      }

    }//for ngrp

    nBuffer.clear();
    sitesToSend.clear();
    nSendLocal.clear();
    nSend.clear();
    offset.clear();

    site_properties_recv.insert( site_properties_recv.end(), sitePropertiesReceived.begin(), sitePropertiesReceived.end() );
        
    //clear memory
    sitePropertiesReceived.clear();

    //check if a proc still has sites to send left
    int nSitesToSendLocal = ( it2 == site_properties.end() ) ? 0 : 1;
    int nSitesToSend;

    MPI::COMM_WORLD.Allreduce(&nSitesToSendLocal, &nSitesToSend, 1, MPI::INT, MPI::SUM );
    if(nSitesToSend == 0){
      isSitesToSend = 0;
    }

  }//while there are sites to send
  
  //insert in vector
  site_properties.insert( site_properties.end(), site_properties_recv.begin(), site_properties_recv.end() );
  
  site_properties_recv.clear();
  
  //remove entries not on this proc
  vector<Site_Update> tmp = site_properties;
  site_properties.clear();
  for( vector<Site_Update>::iterator it=tmp.begin(); it!=tmp.end(); it++ ){
    if(it->get_process() == COMM_RANK){
      site_properties.push_back( *it );
    }
  }
  
  return;
}

//send the intensities from sites that changed process
void SimpleX::send_site_intensities(){
  
  //create permutations vector
  vector<unsigned long int> permutation( intens_ids.size(), 0 );
  for( unsigned int i=0; i<permutation.size(); i++ ){
    permutation[i] = i;
  }
  
  if( intens_ids.size() > 0 ){  
    //sort the intensity ids on vertex id
    quickSortPerm( intens_ids, permutation, 0, intens_ids.size()-1 );
  }
    
  //update the intensities for the new ordering
  vector<float> tmp = site_intensities;
    
  for( unsigned int i=0; i<permutation.size(); i++ ){
    unsigned int posOld = i*numFreq*numPixels;
    unsigned int posNew = permutation[i]*numFreq*numPixels;
    for( unsigned short int j=0; j<numFreq*numPixels; j++ ){
      site_intensities[ posOld + j ] = tmp[ posNew + j ];
    }
  } 
   
  //loop over all vertices and check if intensities need to be sent   
  vector<Vertex>::iterator it=vertices.begin();
  unsigned long long int intensityIndex = 0;
  
  //determine whether there are still sites to be sent
  bool isSitesToSend = 1;
  
  vector<Send_Intensity> intensities_recv;

  while( isSitesToSend ){

    //vector to hold site updates to send
    vector<Send_Intensity> intensitiesToSend;

    //vector to hold the sites that should be send 
    //from this proc to other proc
    vector<unsigned long long int> nSendLocal( COMM_SIZE, 0 );

    //while the maximum number of sites to send is not yet reached, 
    //search for sites that need to be send
    //note that if a site should be send to more than one proc,
    //the chunk will be slightly bigger than the maximum size,
    //but never more than COMM_SIZE    
    while( intensitiesToSend.size() < max_msg_to_send && intensityIndex < intens_ids.size() ){

      if( it->get_vertex_id() < intens_ids[intensityIndex] ){

        it++;

      }else if( it->get_vertex_id() == intens_ids[intensityIndex] ){
        //check whether vertex is on other proc. If so, it needs to be send
        if( it->get_process() != COMM_RANK ){
          //loop over directions
          for( unsigned int j=0; j<numPixels; j++ ){
            //loop over frequencies
            for(short int f=0; f<numFreq; f++){
              Send_Intensity tmpSend;
              tmpSend.set_process( it->get_process() );
              tmpSend.set_neighId( j );
              tmpSend.set_freq_bin( f );
              tmpSend.set_id( intens_ids[intensityIndex] );
              unsigned int pos = intensityIndex*numPixels*numFreq + j*numFreq + f;
              tmpSend.set_intensityIn( site_intensities[ pos ] );

              intensitiesToSend.push_back( tmpSend );
              nSendLocal[ it->get_process() ]++;
            }
          }
        }

        it++;
        intensityIndex++;

      }else if( it->get_vertex_id() > intens_ids[intensityIndex] ){
        cerr << " (" << COMM_RANK << ") Warning: Vertex " << intens_ids[intensityIndex] << " has disappeared from the computational domain " << endl;
        MPI::COMM_WORLD.Abort(-1);
      }

    }
           
    //sort the site updates on process
    sort( intensitiesToSend.begin(), intensitiesToSend.end(), compare_process_send_intensity );
        
    //define the offset of different procs
    vector< unsigned long long int > offset( COMM_SIZE, 0 );
    for( unsigned int p=1; p<COMM_SIZE; p++ ){
      offset[p] = offset[p - 1] + nSendLocal[p - 1];
    }

    //gather all the local numbers of sites to send to all procs
    vector< unsigned long long int > nSend( COMM_SIZE * COMM_SIZE, 0 );
    MPI::COMM_WORLD.Allgather(&nSendLocal[0], COMM_SIZE, MPI::UNSIGNED_LONG_LONG, &nSend[0], COMM_SIZE, MPI::UNSIGNED_LONG_LONG );

    //calculate number of sites to receive
    unsigned long long int nSitesToReceive = 0;
    for( unsigned int p=0; p<COMM_SIZE; p++){
      if( p != COMM_RANK ){
        nSitesToReceive += nSend[ p*COMM_SIZE + COMM_RANK ];
      }
    }

    //vector to hold received intensities
    vector< Send_Intensity > intensitiesReceived( nSitesToReceive );

    int PTask,recvTask;
    int NTask = static_cast<int>( COMM_SIZE );
    int ThisTask = static_cast<int>( COMM_RANK );

    // calculate PTask,  the smallest integer that satisfies COMM_SIZE<=2^PTask
    for(PTask = 0; NTask > (1 << PTask); PTask++) ;

    //vector to hold buffer with number of sites
    //already send to that proc
    vector< unsigned long long int > nBuffer( COMM_SIZE, 0 );

      //loop over all tasks
    for( int ngrp = 1; ngrp < (1 << PTask); ngrp++ ){

      recvTask = ThisTask ^ ngrp;

      if(recvTask < NTask){

        //check if there's sites to receive or to send
        if(nSend[ ThisTask * NTask + recvTask ] > 0 || nSend[ recvTask * NTask + ThisTask ] > 0){

          //do the sending
          MPI::COMM_WORLD.Sendrecv(&intensitiesToSend[ offset[recvTask] ], nSendLocal[recvTask] * sizeof(Send_Intensity),
            MPI::BYTE, recvTask, 0, 
            &intensitiesReceived[ nBuffer[ThisTask] ], nSend[recvTask * NTask + ThisTask] * sizeof(Send_Intensity),
            MPI::BYTE, recvTask, 0);

        }//if there's sites to send or receive

      }//if receiving task is smaller than COMM_SIZE

      //update the buffer size
      for( int p = 0; p < NTask; p++){
        if( (p ^ ngrp) < NTask ){
          nBuffer[p] += nSend[(p ^ ngrp) * NTask + p];
        }
      }

    }//for ngrp

    nBuffer.clear();
    intensitiesToSend.clear();
    nSendLocal.clear();
    nSend.clear();
    offset.clear();

    intensities_recv.insert( intensities_recv.end(), intensitiesReceived.begin(), intensitiesReceived.end() );
        
    //clear memory
    intensitiesReceived.clear();

    //check if a proc still has sites to send left
    int nSitesToSendLocal = ( intensityIndex == intens_ids.size() ) ? 0 : 1;
    int nSitesToSend;

    MPI::COMM_WORLD.Allreduce(&nSitesToSendLocal, &nSitesToSend, 1, MPI::INT, MPI::SUM );
    if(nSitesToSend == 0){
      isSitesToSend = 0;
    }

  }//while there are sites to send
  
  //sort the received intensities on vertex id
  sort( intensities_recv.begin(), intensities_recv.end(), compare_index_send_intensity );
  
  //insert in vectors
  vector<Send_Intensity>::iterator it2=intensities_recv.begin(); 
  while(it2!=intensities_recv.end() ){

    unsigned long long int index =  it2->get_id();
    intens_ids.push_back( index );

    //create space for the frequencies
    site_intensities.insert( site_intensities.end(), numFreq*numPixels, 0.0 );

    for(unsigned int j=0; j<numPixels*numFreq; j++){

      unsigned int pos = it2->get_freq_bin() + it2->get_neighId()*numFreq + site_intensities.size() - numFreq*numPixels;
          
      //check to be sure
      if( it2->get_id() != index ){
        cerr << " (" << COMM_RANK << ") Error in sendSiteIntensities(): Not all bins have been received " << endl;
        cerr << "     Direction bin: " << it2->get_neighId() << " Frequency bin: " << it2->get_freq_bin() << endl;
        cerr << "     Index: " << index << endl;
        cerr << "     Total number of intensities received: " << intensities_recv.size() << endl;
        MPI::COMM_WORLD.Abort(-1);
      }

      site_intensities[pos] =  it2->get_intensityIn();
      it2++;
    }
  }
  
  //remove entries not on this proc
  //probably not necessary :)
  return;
  
}


/****************************************************************************************/
/*                            Physics Functions                                         */
/****************************************************************************************/

void SimpleX::compute_physics( const unsigned int& run ){

  double t0 = MPI::Wtime();

  if( run == 0 ){

    initialise_physics_params();

    assign_read_properties();

    if(blackBody) {   // -> Blackbody spectrum
      black_body_source( sourceTeff );
    } else {             // -> Monochromatic spectrum
      //if recombination photons are included, use 2 bins
      if(rec_rad){
        
        cross_H.resize(2, cross_HI( nu0HI, nu0HI ) );
        
        //put the flux of the source in first bin instead of zeroth,
        //that is used for recombination radiation
        for( SITE_ITERATOR it=sites.begin(); it!=sites.end(); it++ ){
          if(it->get_source()) {
            float flux = it->get_flux( 0 );
            it->set_flux( 0, 0.0 ); 
            it->set_flux( 1, flux );
          }
        }

      }else{
        cross_H.resize(1, cross_HI( nu0HI, nu0HI ) );
      }
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
    vector< unsigned long int >().swap( intens_ids );

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
  UNIT_T = UNIT_T_MYR * secondsPerMyr;  
//  UNIT_T *= secondsPerMyr;  
  UNIT_L = sizeBox * parsecToCm; 
  UNIT_V = pow(UNIT_L, 3.0);

  //read in metal line cooling data 
  if(metal_cooling)
    read_metals();
  
  if( COMM_RANK == 0 )
    simpleXlog << "  Essential physical parameters calculated " << endl;

}

/****  Read metal line cooling data  ****/
void SimpleX::read_metals(){

  string datadir = dataPath + "/";

  // Read the curves
  cooling_curve hydrogen(datadir+"H.cool");
  cooling_curve helium(datadir+"He.cool");
  cooling_curve carbon(datadir+"C.cool");
  cooling_curve carbonH(datadir+"C.Hcool");
  cooling_curve nitrogen(datadir+"N.cool");
  cooling_curve oxigen(datadir+"O.cool");
  cooling_curve oxigenH(datadir+"O.Hcool");
  cooling_curve neon(datadir+"Ne.cool");
  cooling_curve silicon(datadir+"Si.cool");
  cooling_curve siliconH(datadir+"Si.Hcool");
  cooling_curve iron(datadir+"Fe.cool");
  cooling_curve ironH(datadir+"Fe.Hcool");

  // Place them in the cooling curve vector
  curves.push_back(hydrogen);
  curves.push_back(helium);
  curves.push_back(carbon);
  curves.push_back(carbonH);
  curves.push_back(nitrogen);
  curves.push_back(oxigen);
  curves.push_back(oxigenH);
  curves.push_back(neon);
  curves.push_back(silicon);
  curves.push_back(siliconH);
  curves.push_back(iron);
  curves.push_back(ironH);

  
  // Solar abundances (by number)
  // (/"H   ","He  ","C   ","N   ","O   ","Ne  ","Si  ","Fe  "/) 
  // 1., .085, 3.31e-4, 8.3e-5, 6.76e-4, 1.2e-4, 3.55e-5, 3.2e-5

  // Initialize the abundances relative to the total number density of atoms
  abundances.push_back(1.0);//H
  abundances.push_back(0.085);//He
  abundances.push_back(3.31e-4);//C
  abundances.push_back(3.31e-4);//C
  abundances.push_back(8.3e-5);//N
  abundances.push_back(6.76e-4);//O
  abundances.push_back(6.76e-4);//O
  abundances.push_back(1.2e-4);//Ne
  abundances.push_back(3.55e-5);//Si
  abundances.push_back(3.55e-5);//Si
  abundances.push_back(3.2e-5);//Fe
  abundances.push_back(3.2e-5);//Fe

  if( COMM_RANK == 0 ){
    simpleXlog << "  Metal line cooling curves initialized " << endl;
    //cerr << "  Metal line cooling curves initialized " << endl;
  }

  return;
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
  //set internal energy
        it->set_internalEnergy( temp_u_list[ it->get_vertex_id() ] );
  //set temperature
        double mu = compute_mu( *it );
        float T = static_cast<float>( u_to_T( temp_u_list[ it->get_vertex_id() ],  mu ) );
        it->set_temperature( T );

  //set dudt
        it->set_dinternalEnergydt( temp_dudt_list[ it->get_vertex_id() ] );

        //metallicity
        it->set_metallicity( temp_metallicity_list[it->get_vertex_id() ] );

  //set number of ionising photons
        if(temp_flux_list[ it->get_vertex_id() ] > 0.0){

    //site is source
          it->set_source(1);

    //create flux array
          it->create_flux(numFreq);
    //put total number of ionsing photons in first bin
          it->set_flux( 0, temp_flux_list[ it->get_vertex_id() ] );
        }else{
    //site is no source
          it->set_source(0);
        }

      }else{
        it->set_n_HI( 0.0 );
        it->set_n_HII( 0.0 );
        it->set_source( 0 );
        it->set_internalEnergy( 0.0 );
        it->set_temperature( 0.0 );
        it->set_metallicity( 0.0 );
      }
    }
  }

  //if clumping factor is included, add to sites
  if(temp_clumping_list.size() > 0){
    for( SITE_ITERATOR it=sites.begin(); it!=sites.end(); it++ ){
      if( it->get_process() == COMM_RANK && !it->get_border()){
        it->set_clumping( temp_clumping_list[ it->get_vertex_id() ] );
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
  temp_u_list.clear();
  vector< float >().swap( temp_u_list );
  temp_dudt_list.clear();
  vector< float >().swap( temp_dudt_list );
  temp_clumping_list.clear();
  vector< float >().swap( temp_clumping_list );
  temp_metallicity_list.clear();
  vector< float >().swap( temp_metallicity_list );

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
  unsigned long int i=0;
  for( SITE_ITERATOR it=sites.begin(); it!=sites.end(); it++ ){
    //only include sites on this proc that are not in border
    if( it->get_process() == COMM_RANK ){
      if( !it->get_border() ) {
        //check if the vertex ids match
        if(it->get_vertex_id() == site_properties[i].get_vertex_id() ){

	  //assign properties
	  it->set_n_HI( site_properties[ i ].get_n_HI() );
	  it->set_n_HII( site_properties[ i ].get_n_HII() ); 
	  it->set_ballistic( site_properties[ i ].get_ballistic() );
	  it->set_internalEnergy( site_properties[ i ].get_internalEnergy() );
    it->set_metallicity( site_properties[ i ].get_metallicity() );
	  //set temperature
    double mu = compute_mu( *it );
    float T = static_cast<float>( u_to_T( it->get_internalEnergy(),  mu ) );
    it->set_temperature( T );
  	
	  //if site is source, put flux in first bin
	  if( site_properties[ i ].get_flux() > 0.0 ){
	    it->set_source(1);
	    it->create_flux(numFreq);
	    it->set_flux( 0, site_properties[ i ].get_flux() );
 
	  }else{
	    it->set_source(0);
	  }

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


  if(blackBody){
    
    if(rec_rad){
      numFreq--;
    }
    
  //fill the source bins correctly
  //upper and lower bounds of integration in terms of the ionisation frequency for hydrogen
  //tweak the upper bound in case of evenly spaced bins
    const double upperBound = 1e2;
    const double lowerBound = 1.0;

  //calculate the frequency bounds
    vector<double> freq_bounds = calc_freq_bounds(lowerBound, upperBound);

  //Proper number of frequencies
  if(rec_rad){
    numFreq++;
    //add extra bin
    vector<double>::iterator it = freq_bounds.begin();
    freq_bounds.insert(it, 1, lowerBound);
  }
  
  //Planck curve divided by h * nu over entire domain
    double normBB = qromb(PlanckDivNu, lowerBound, upperBound, sourceTeff);



  //loop over all sites
    for( SITE_ITERATOR it=sites.begin(); it!=sites.end(); it++ ){
    //if site is source
      if( it->get_source() ){
        double total_flux = it->get_flux(0);
      //loop over frequencies
        for(short int f=0; f<numFreq; f++){

  // Integrated number of photons for this bin
          double bbOverNu = qromb(PlanckDivNu, freq_bounds[f], freq_bounds[f+1], sourceTeff);

  // Normalise to the source strength
          double flux = total_flux * bbOverNu/normBB;
          it->set_flux( f, (float) flux);

        }//for all freqs
      }//if source
    }//for all sites
  }else{
        //put the flux of the source in first bin instead of zeroth,
        //that is used for recombination radiation
    for( SITE_ITERATOR it=sites.begin(); it!=sites.end(); it++ ){
      if(it->get_source()) {
        float flux = it->get_flux( 0 );
        it->set_flux( 0, 0.0 ); 
        it->set_flux( 1, flux );
      }
    }
  }

  if( COMM_RANK == 0 ){
    simpleXlog << " Returned physics to proper sites" << endl;
  }
}


//calculate the bounds of the frequency bins
vector<double> SimpleX::calc_freq_bounds( const double& lowerBound, const double& upperBound ){

  //tolerance in calculation
  double TOL = 1e-8;
  //total of the quantity that determines the bin spacing
  double total = 0.0;

  //vector that holds the spacing of the frequency bins
  vector<double> freq_bounds(numFreq+1, 0.0);

  //in case of energy and ionisation weigths, calculate the total quantity within the range 
  switch(freq_spacing){

  case ENERGY_WEIGHTS:

    // Use the energy of the radiation field to choose bins such that they carry equal energy
    total = qromb(Planck, lowerBound, upperBound, sourceTeff);
    break;

  case IONISATION_WEIGHTS:

    // Use the enrgy times the cross section to choose bins
    total = crossIntHI( upperBound, sourceTeff );

    break;

  }

  // Cross section per bin
  double perBin = total/double(numFreq);

  /*  
  // Use minimisation routine to find the bounding frequencies for the bins
  */
  double deviation = FLT_MAX;

  // Store the lower and upper bound
  freq_bounds[0] = lowerBound;
  freq_bounds[numFreq] = upperBound;


  //determine the spacing of the rest of the bins    
  switch(freq_spacing){

  case ENERGY_WEIGHTS:

    for(short int f=1; f<numFreq; f++){
      double bound = zbrent(PlanckInt, upperBound, TOL, sourceTeff, f*perBin);
      freq_bounds[f] = bound;
    }


    // Check whether the last bin contains the correct number of ionisations per unit time 
    if(numFreq>1){
      deviation = ( qromb(Planck, freq_bounds[freq_bounds.size()-2], upperBound, sourceTeff) - 
		    perBin)/perBin;
    }else{
      deviation = 0.0;
    }

    break;

  case IONISATION_WEIGHTS:

    // Loop over bins
    for(short int f=1; f<numFreq; f++){
      double bound =  zbrent(crossIntHI, upperBound, TOL, sourceTeff, f*perBin);
      freq_bounds[f] = bound;
    }


    if(numFreq>1){
      // Check whether the last bin contains the correct number of ionisations per unit time 
      double total = qromb(fHI, freq_bounds[freq_bounds.size()-2], upperBound, sourceTeff);
      deviation = ( total - perBin)/perBin;
    } else {
      deviation = 0.0;
    }

    break;

  case LOG_WEIGHTS:

    // logarithmic bins
    double logLowerBound;
    logLowerBound = log10(lowerBound);
    double logUpperBound;
    logUpperBound = log10(upperBound);

    deviation = 0.0;
    break;

  case NO_WEIGHTS:

    //upperBound = 3;// This is not entirely correct but if spacing linearly to 100 everything happens in the first bin.

    if(numFreq>1){    
      for(short int f=1; f<numFreq; f++){
	freq_bounds[f] = f*(upperBound-lowerBound)/numFreq + lowerBound;
      }
    }

    deviation = 0.0;
    break;

  }

  if(deviation > TOL){
    cerr << "WARNING: Integration over last bin deviates more than TOL " << TOL << " from perBin by: " << deviation << endl;
    cerr << "This may signal inaccurate cross sections." << endl;
    MPI::COMM_WORLD.Abort(-1);
  }

  return freq_bounds;

}

//calculate effective cross section and spectrum of black body source
void SimpleX::black_body_source(const double& tempSource) {


  //Determine the frequency integrated cross section of hydrogen for a black bosy source of temperature tempSource
  //Units are scaled to close to 1

  //upper and lower bounds of integration in terms of the ionisation frequency for hydrogen
  //tweak the upper bound in case of evenly spaced bins
  const double upperBound = 1e2;
  const double lowerBound = 1.0;

  //if recombination radiation is included, zeroth bin is for radiation at Lyman limit
  //so is not included in calculation of bin spacing
  if(rec_rad){
    numFreq--;
  }
    
  //calculate the frequency bounds
  vector<double> freq_bounds = calc_freq_bounds(lowerBound, upperBound);

  //Proper number of frequencies
  if(rec_rad){
    numFreq++;
    //add extra bin
    vector<double>::iterator it = freq_bounds.begin();
    freq_bounds.insert(it, 1, lowerBound);
  }

 
  // cerr << " Frequency bounds: ";
  // for(int i=0; i<freq_bounds.size(); i++){
  //   cerr << freq_bounds[i] << " ";
  // }
  // cerr << endl;

  //***** Fill the source bins *****//

  //Planck curve divided by h * nu over entire domain
  double normBB = qromb(PlanckDivNu, lowerBound, upperBound, tempSource);

  //cerr << " Flux: ";

  //loop over all sites
  for( SITE_ITERATOR it=sites.begin(); it!=sites.end(); it++ ){
    //if site is source
    if( it->get_source() ){
      double total_flux = it->get_flux(0);
      //loop over frequencies
      for(short int f=0; f<numFreq; f++){

	// Integrated number of photons for this bin
	double bbOverNu = qromb(PlanckDivNu, freq_bounds[f], freq_bounds[f+1], tempSource);

	// Normalise to the source strength
	double flux = total_flux * bbOverNu/normBB;
	it->set_flux( f, (float) flux);

	//cerr << flux << " ";

      }//for all freqs
    }//if source
  }//for all sites

  //cerr << endl;

  //***** Determine cross sections ******//

  //give cross sections vector the correct size
  cross_H.resize(numFreq, 0.0);

  //store integral over cross section times source spectrum
  //vector<double> intHI(numFreq, 0.0);
  //loop over frequencies
  // for(short int f=0; f<numFreq; f++){
  //   //integral of hydrogen cross section times Planck curve
  //   intHI[f] = qromb(fHI, lowerBound, freqBounds[f+1], tempSource ) ;
  // }
  // // Take the differences of the entries computed above to obtain the integrals over the bins
  // for(short int f=numFreq-1; f>0; --f){
  //   intHI[f] -= intHI[f-1];
  // }
  // //loop over frequencies and normalise the cross sections
  // for(short int f=0; f<numFreq; f++){
  //   //PLanck curve divided by h*nu
  //   double bbOverNu = qromb(PlanckDivNu, freqBounds[f], freqBounds[f+1], tempSource);
  //   intHI[f] /= bbOverNu;
  //   //back to physical units
  //   cross_H[f] = intHI[f]*crossFactor;
  // }

  //loop over frequencies
  for(short int f=0; f<numFreq; f++){

    //integral of hydrogen cross section times Planck curve
    double intHI = qromb(fHI, freq_bounds[f], freq_bounds[f+1], tempSource ) ;

    //Planck curve divided by h*nu
    double bbOverNu = qromb(PlanckDivNu, freq_bounds[f], freq_bounds[f+1], tempSource);

    //cross section in physical units
    cross_H[f] = crossFactor * intHI/bbOverNu;

    //cerr << " Cross section: " << cross_H[f] << " crossFactor: " << crossFactor << " intHI: " << intHI << " bbOverNu: " << bbOverNu << endl;

  }

  //if recombination radiation is included, cross section is cross section at Lyman limit
  if(rec_rad){
    cross_H[0] = cross_HI(nu0HI,nu0HI);
  }
 
  //***** Determine excess energy of the photons ******//

  //give photon excess energy vector correct size
  photon_excess_energy.resize(numFreq, 0.0);

  //loop over frequencies
  for(short int f=0; f<numFreq; f++){

    //photoheating coefficient
    double photo_heat_coeff = qromb(enerfHI, freq_bounds[f], freq_bounds[f+1], tempSource );

    //photo-ionisation rate
    double Gamma_H = qromb(fHI, freq_bounds[f], freq_bounds[f+1], tempSource );

    //set the photon excess energy (convert units)
    photon_excess_energy[f] = (Gamma_H > 0.0) ? planckConstant * nu0HI * photo_heat_coeff/Gamma_H : 0.0;

  }

  //if recombination radiation is included, excess energy is k*T_gas
  //for now assume that all gas that emits photons is at 10^4 K
  if(rec_rad){
    photon_excess_energy[0] = k_B * 1.e4;
  }
  
  // cerr << " Heating: ";
  // for(int i=0; i<photon_excess_energy.size(); i++){
  //   cerr << photon_excess_energy[i] << " ";
  // }
  // cerr << endl;

  if( COMM_RANK == 0 ){
    simpleXlog << "  Calculated effective cross sections and photon excess energy for black body source" << endl;
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

  double tau_tot=0.0; 
  unsigned int count = 0; 
  for( SITE_ITERATOR it=sites.begin(); it!=sites.end(); it++ ){
    if( !it->get_border() && it->get_process() == COMM_RANK ){ 
      //calculate the mean optical depth over frequencies
      double tau=0.0;
      for(short int f=0; f<numFreq; f++){
	tau += it->get_n_HI() * UNIT_D * it->get_neigh_dist() * UNIT_L * cross_H[f] * straight_correction_factor;
      }
      tau /= numFreq;
      tau_tot += tau;
      count++;  
    }
  }  

  double totalTau;
  MPI::COMM_WORLD.Allreduce(&tau_tot,&totalTau,1,MPI::DOUBLE,MPI::SUM);
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

  int maxRes;

  MPI::COMM_WORLD.Reduce(&localMaxRes,&maxRes,1,MPI::INT,MPI::MIN,0);

  double mean_tau = output_optical_depth();

  if(COMM_RANK == 0){

    cerr << endl << endl
	 << "  Total volume: " << totalVolume << " (should be close to " << 1.0 << ")." << endl
	 << "  Max resolution: " << maxRes << "^3 effective (adaptive grid)." << endl
	 << endl;
    cerr << "  Effective Photo-Ionization Cross Section     : ";
    for(short int f=0; f<numFreq-1; f++){         
      cerr << cross_H[f] << " | ";
    }
    cerr << cross_H[numFreq-1] << " cm^2" << endl;

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
	       << "  Max resolution: " << maxRes << "^3 effective (adaptive grid)." << endl
	       << endl;         
    simpleXlog << "  Effective Photo-Ionization Cross Section     : ";
    for(short int f=0; f<numFreq-1; f++){         
      simpleXlog << cross_H[f] << " | ";
    }
    simpleXlog << cross_H[numFreq-1] << " cm^2" << endl;
    simpleXlog << endl  
	       << "  Average density SimpleX grid                 : " << totalAtoms/pow( UNIT_L, 3 ) << " cm^-3" << endl  
	       << "  Total number of atoms within SimpleX grid    : " << totalAtoms << endl  
	       << endl
	       << "  Average optical depth                        : " << mean_tau << endl
	       << endl;

  }


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

      //number of recombinations is this factor times recombination coefficient
      double num_rec_factor = (double) it->get_clumping() * N_HII * (double) it->get_n_HII() * UNIT_D * UNIT_T;

      //number of recombinations
      double num_rec = (rec_rad) ? num_rec_factor*recomb_coeff_HII_caseA( it->get_temperature() ) 
        : num_rec_factor*recomb_coeff_HII_caseB( it->get_temperature() );

      //make sure there are no more recombinations than ionised atoms
      if( num_rec > N_HII ) {
	      num_rec = N_HII;  
      } 

      //if recombination photons are included, their contribution to ionisations 
      //inside the cell needs to be taken into account
      //the rest is sent to all neighbours
      if( rec_rad ) { 

        //optical depth through the cell
        double tau = (double) it->get_n_HI() * UNIT_D * (double) it->get_neigh_dist() * UNIT_L * cross_H[0];
        //fraction of the photons that can escape
        double rec_escape = rec_rad_escape( tau );

        //total number of recombination photons produced in the cell
        double N_out_diffuse = num_rec_factor * ( recomb_coeff_HII_caseA( it->get_temperature() ) - recomb_coeff_HII_caseB( it->get_temperature() ) );

        //the fraction that is not escaping balances recombinations
        num_rec -= (1.0 - rec_escape) * N_out_diffuse;

        if(num_rec < 0.0){
          cerr << " (" << COMM_RANK << ") Negative number of recombinations! " << endl
            << " rec_escape: " << rec_escape << " tau: " << tau << endl;
        }

        //send the fraction that escapes to all neighbours
        vector<double> N_out(numFreq, 0.0);
        //put the recombination photons in zeroth bin
        N_out[0] = rec_escape * N_out_diffuse/UNIT_I;
        if(N_out_diffuse > 0.0){
          diffuse_transport( *it, N_out );
        }

      }


      totalRecombs += num_rec;

      //add the recombinations to the sites
      N_HI += num_rec;
      N_HII -= num_rec;
      double tmp = N_HI/(UNIT_D * (double) it->get_volume() * UNIT_V);
      it->set_n_HI( (float) tmp );
      tmp = N_HII/(UNIT_D * (double) it->get_volume() * UNIT_V);
      it->set_n_HII( (float) tmp );

    }
  } //for all vertices 


  return totalRecombs;

}

//rate of energy loss inside cell
double SimpleX::cooling_rate( Site& site ){

  //normalised cooling function
  double C = 0.0;
  //neutral hydrogen density
  double n_HI = site.get_n_HI() * UNIT_D;
  //ionised hydrogen density
  double n_HII = site.get_n_HII() * UNIT_D;
  //total number density
  double n_H = n_HI + n_HII;

  //electron density is the same in case of hydrogen only
  double n_e = n_HII;

  //recombination cooling
  if(rec_rad){
    C += recomb_cooling_coeff_HII_caseA( site.get_temperature() ) * n_e * n_HII * site.get_clumping();
  }else{
    C += recomb_cooling_coeff_HII_caseB( site.get_temperature() ) * n_e * n_HII * site.get_clumping();;
  }



  //collisional ionisation cooling
  C += coll_ion_cooling_coeff_HI( site.get_temperature() ) * n_e * n_HI;


  //collisional de-excitation cooling
  C += coll_excit_cooling_coeff_HI( site.get_temperature() ) * n_e * n_HI;
  

  //free-free cooling
  C += ff_cooling_coeff( site.get_temperature() ) * n_e * n_HII;


  //metal cooling (including helium for the moment)
  if(metal_cooling){
    
    double total = 0.0;
        
    //loop over metal line cooling curves
    for(unsigned int j=1; j<curves.size();j++){// j=1 because hydrogen is not taken into account here (but explicitly above)
      
      // abundance*cooling/rate
      if(withH[j]){
	      C += abundances[j]*curves[j].get_value(site.get_temperature()) * n_H * n_H * site.get_metallicity();
        total += abundances[j]*curves[j].get_value(site.get_temperature()) * n_H * n_H * site.get_metallicity();
      }else{
        C += abundances[j]*curves[j].get_value(site.get_temperature()) * n_H * n_e * site.get_metallicity();
        total += abundances[j]*curves[j].get_value(site.get_temperature()) * n_H * n_e * site.get_metallicity();
      }
    }
    // if(site.get_metallicity() > 0.0 ){
    //   
    //   cout << " nH: " << n_H << endl;
    //   cout << " nHI: " << n_HI << endl;
    //   cout << " nHII: " << n_HII << endl;
    //   cout << " nE: " << n_e << endl;
    //   
    //   cout << " Recombination cooling: " << recomb_cooling_coeff_HII_caseB( site.get_temperature() ) * n_e * n_HII * site.get_clumping() << endl;
    //   cout << " Collisional ionisation cooling: " << coll_ion_cooling_coeff_HI( site.get_temperature() ) * n_e * n_HI << endl;
    //   cout << " Collisional excitation cooling: " << coll_excit_cooling_coeff_HI( site.get_temperature() ) * n_e * n_HI << endl;
    //   cout << " Free-free cooling: " << ff_cooling_coeff( site.get_temperature() ) * n_e * n_HII << endl;
    //   
    //   cout << " Metal cooling: " << total << " metallicity: " << site.get_metallicity() << endl;
    //   cout << " Cooling: " << C * site.get_volume() * UNIT_V << endl;
    //   
    //   exit(-1);
    // }
  }
  
  //total cooling rate in cell
  double cooling = C * site.get_volume() * UNIT_V;


  return cooling;

}

//rate of energy gain inside cell
double SimpleX::heating_rate( const vector<double>& N_ion, const double& t_end ){

  //normalised heating function H = n_HI * Gamma * E/n_H^2
  //total heating in this cell is heating = n_H^2 * V * H,
  //which is given here

  //energy gain
  double heating = 0.0;

  //loop over frequencies
  for( short int f=0; f<numFreq; f++ ){
  
    //number of photons per second
    double N_ion_sec = N_ion[f]/t_end;

    //heating due to these photons
    heating += N_ion_sec * photon_excess_energy[f];
  }

  return heating;

}

double SimpleX::compute_mu(Site& site){
  
  double mu = ( site.get_n_HI() + site.get_n_HII() )/(site.get_n_HI() + 2*site.get_n_HII() );
  
  return mu;
  
}

//convert internal energy to temperature using ideal gas law
double SimpleX::u_to_T( const double& u, const double& mu ){

  //temperature
  double T = (mu * m_H * u)/(1.5 * k_B);

  return T;

}

//convert temperature to internal energy using ideal gas law
double SimpleX::T_to_u( const double& T, const double& mu ){

  //internal energy
  double u = (1.5 * k_B * T)/(mu * m_H);

  return u;

}

//update the temperature state of the gas
// function returns the heating/cooling time
double SimpleX::update_temperature( Site& site, const vector<double>& N_ion, const double& t_end ){

  // double uBefore = site.get_internalEnergy();
  // double Tbefore = site.get_temperature();

  //range of temperatures between which the rates are valid
  double T_min = 10.;
  double T_max = 1.e9;

  //time passed in heating cycling
  double t = 0.0;

  //heating/cooling time scale
  double t_heat = t_end;

  //mean molecular weight
  double mu = ( site.get_n_HI() + site.get_n_HII() )/(site.get_n_HI() + 2*site.get_n_HII() );

  double u_min = T_to_u( T_min, mu);
  double u_max = T_to_u( T_max, mu);

  //number of hydrogen atoms
  double N_H = ( site.get_n_HI() + site.get_n_HII() ) * UNIT_D * site.get_volume() * UNIT_V;

  //inclusion of adiabatic cooling term from hydro
  double duDtAdiabatic = site.get_dinternalEnergydt();
  
  //divide this by the internal energy to obtain a constant quantity
  //double u_0 = site.get_internalEnergy();

  // not a good idea (shock heating terms)
  // duDtAdiabatic = (u_0 > 0.0)? duDtAdiabatic/u_0 : 0.0;
  
  //loop until end time is reached 
  while( t < t_end ){

    //current internal energy of the cell per unit mass
    double u = site.get_internalEnergy();

    //total heating inside this cell per second
    double H = heating_rate( N_ion, t_end );

    //total cooling inside this cell per second
    double C = cooling_rate( site ); 
        
    //change in internal energy per second per unit mass
    double du = duDtAdiabatic + (H - C)/(N_H * m_H);
    // removed *u

    //determine heating time scale
    t_heat = ( fabs(du) > 0.0 ) ?  u/fabs(du) : 1.e5*UNIT_T; 

    //take time step fraction of the heating time
    double dt = subcycle_frac*t_heat;

    //check that we stay within the total time
    if( (t+dt) > t_end ){
      dt = t_end - t;
    }

    // if( site.get_vertex_id() == 85 ){
    //   
    //   cerr << " Density: " << ( site.get_n_HI() + site.get_n_HII() ) * UNIT_D
    //        << " Metallicity: " << site.get_metallicity() 
    //        << " Temperature: " << site.get_temperature() 
    //        << " Internal energy: " << u 
    //        << " Cooling rate per unit mass: " << C/(N_H * m_H) 
    //        << " du: " << du 
    //        << " dt: " << dt 
    //        << endl;
    // }
    
    //change in internal energy in total time 
    du *= dt;

    //make sure internalEnergy is in range where rates are valid
    u = u+du;

    if(u > u_max){
      u = u_max;
    }else if(u < u_min){
      u = u_min;
    }
    
    //set new internal energy
    site.set_internalEnergy( (float) u );
  	//set temperature
    float T = static_cast<float>( u_to_T( u,  mu ) );
    site.set_temperature( T );



    //add time step to total time
    t += dt; 

  }


  // double uAfter = site.get_internalEnergy();
  // double Tafter = site.get_temperature();

  // cout << " Temperature: " << Tbefore << " " << Tafter << endl;
  // cout << " Energy     : " << uBefore << " " << uAfter <<  " " << duDtAdiabatic << endl;
  // exit(-1);
    
  return t_heat;

}




/****************************************************************************************/
/*                       Radiative Transport Functions                                  */
/****************************************************************************************/



//Calculate the 
vector<double> SimpleX::solve_rate_equation( Site& site ){

  static int count = 0;
  
  //calculate number density and optical depth at beginning of time step
  double n_H = ( (double) site.get_n_HI() + (double) site.get_n_HII() ) * UNIT_D;

  vector<double> initial_tau(numFreq,0.0); 
  for(short int f=0; f<numFreq;f++){
    initial_tau[f] = (double) site.get_n_HI() * UNIT_D * (double) site.get_neigh_dist() * UNIT_L * cross_H[f] * straight_correction_factor;

    //dirty hack to avoid nans
    if(initial_tau[f] == 0.0 && n_H > 0.0){
      initial_tau[f] = 1.e-15;
    }
  }

  //total number of incoming and outgoing photons
  vector<double> N_in_total(numFreq, 0.0);
  vector<double> N_out_total(numFreq, 0.0);

  //number of neutral and ionised atoms in the cell
  double initial_N_HI = (double) site.get_n_HI() * UNIT_D * (double) site.get_volume() * UNIT_V;
  double initial_N_HII = (double) site.get_n_HII() * UNIT_D * (double) site.get_volume() * UNIT_V;

  //in case of ballistic transport, intensity has size of number of neighbours;
  //in case of direction conserving transport, intensity has 
  //the size of the tesselation of the unit sphere
  numPixels = ( site.get_ballistic() ) ? site.get_numNeigh() : number_of_directions;


  for( unsigned int j=0; j < numPixels; j++) {
    for(short int f=0; f<numFreq;f++){
      //incoming intensity in physical units
      N_in_total[f] += (double) site.get_intensityOut( f,j ) * UNIT_I;
    } //for all freqs
  } //for all neighbours

  //if rho is zero, no need to go into loop
  if(n_H == 0.0){
    for(short int f=0; f<numFreq;f++){
      N_out_total[f] = N_in_total[f];
    }
  } else {//this loop is only useful if there is gas in the cell
    
    // NOTE: if there are no collisional ionizations, this loop is only useful if there is radiation!

    //calculate the number of absorptions and ionisations if N_HI and N_HII were constant during RT time step
    //double initial_N_retained_total = (1.0 - exp(-initial_tau)) * N_in_total;
    //double initial_numIonised = ( initial_N_retained_total >= N_HI ) ? N_HI : initial_N_retained_total;

    //absorptions by hydrogen
    vector<double> initial_N_retained_H(numFreq, 0.0);
    //total number of photoionisations in this time step
    double initial_photo_ionisations = 0.0;
    for(short int f=0; f<numFreq;f++){
      initial_N_retained_H[f] = (1.0 - exp( -initial_tau[f] )) * N_in_total[f];
      initial_photo_ionisations += initial_N_retained_H[f];

      if(initial_photo_ionisations != initial_photo_ionisations ){
      	cerr << initial_photo_ionisations << " " << initial_N_retained_H[f] << " " << initial_tau[f] << " " << N_in_total[f] << endl;
      }

    }

    //total number of collisional ionisations in this time step
    double initial_coll_ionisations = 0.0;
    if(coll_ion){
      double n_e = (double) site.get_n_HII() * UNIT_D;
      //number of coll ionisations is coll ion coeff * n_e * number of atoms * time step
      initial_coll_ionisations = coll_ion_coeff_HI( site.get_temperature() ) * n_e * initial_N_HI * UNIT_T;
    }

    //the initial ionisations without accounting for the number of neutral atoms
    double initial_numIonised = initial_photo_ionisations + initial_coll_ionisations;

    //use temporal photon condervation or not?
    unsigned long int Nsteps = 1;
    if(photon_conservation){

      //calculate relevant time scales to see whether assumption of constant N_HI and N_HII is correct

      //ionisation time scale, if no ionisations it is taken to be much larger than time step
      double t_ion = (initial_numIonised > 1.e-20) ? UNIT_T * initial_N_HI/(initial_numIonised) : 1.e5 * UNIT_T; 

      //total #ions available for recombinations
      double tot_ions = initial_N_HII + initial_numIonised;

      //recombination time scale, if no recombinations it is taken to be much larger than time step
      double t_rec = (tot_ions > 0.0 ) ? ( site.get_volume() * pow( UNIT_L, 3.0 ) )/
	( tot_ions * site.get_clumping() * recomb_coeff_HII_caseB( site.get_temperature() ) ) : 1.e5 * UNIT_T;

      //'equality' time scale
      double t_eq = t_ion*t_rec/(t_ion + t_rec);

      if(t_eq != t_eq){
        cerr << endl << " (" << COMM_RANK << ") NAN detected! t_eq: " << t_eq << " t_ion: " << t_ion << " t_rec: " << t_rec << endl;
        cerr << "  initial_N_retained_H[0]: " << initial_N_retained_H[0]
          << "  initial_tau[0]: " << initial_tau[0] << endl;
        cerr << "  initial_N_HI " << initial_N_HI << " initial_N_HII " << initial_N_HII << endl;
        cerr << "  initial_neutral_fraction: " << initial_N_HI/(initial_N_HI+initial_N_HII)
          << " N_in_total[0]: " << N_in_total[0]
          << " n_H: " << n_H << endl;

        cerr << "  initial_N_retained_total[0] = (1.0 - exp(-"<<initial_tau[0]<<")) * "<< N_in_total[0] <<" " << endl;
        cerr << "  volume " << site.get_volume() * UNIT_V << " N_HII + initial_numIonised " 
          << + initial_numIonised << " fclump " 
          << site.get_clumping() << " alpha " << recomb_coeff_HII_caseB( site.get_temperature() ) << endl;
        cerr << "  coll_ion " <<  initial_coll_ionisations << endl; 
        cerr << "  tot_ions " <<  tot_ions << endl;
        cerr << " UNIT_T " << UNIT_T << endl;
        exit(1);

      }
       
      //subcycling is necessary when the RT time step is larger than 
      //some factor times ionisation time scale or recombination time scale 
      if( UNIT_T > subcycle_frac*t_eq ){
      	Nsteps = (unsigned long int) ceil( UNIT_T/(subcycle_frac*t_eq));
      }

      //keep track of the total number of photo-ionisations
      vector<double> number_of_photo_ionisations(numFreq, 0.0);
      double number_of_recombinations = 0.0;

      //total number of diffuse recombination photons during this time step
      double N_rec_phot = 0.0;

      //number of neutral and ionised atoms during subcycle step
      double N_HI_step = initial_N_HI;
      double N_HII_step = initial_N_HII;

      //time step during subcycle step
      double dt_ss = UNIT_T/(double) Nsteps;

      //number of ionisation in subcycle step
      double numIonised_step = initial_numIonised;

      //effective time step in calculations of ionised atoms
      double dt_eff = dt_ss/UNIT_T;
      //one over the initial neuatral atoms
      double one_over_initial_N_HI = (initial_N_HI > 0.0 ) ? 1.0/initial_N_HI : 1.0;
      //one over the physical volume
      double one_over_volume = 1.0/( site.get_volume() * UNIT_V );
      //convert number of atoms to number density in code units
      double num_to_dens = 1.0/( UNIT_D * site.get_volume() * UNIT_V ); 

      //optical depth in time step
      vector<double> tau_step(numFreq,0.0);
      //number of ionisations in time step
      vector<double> photo_ionisations_step(numFreq,0.0);

      //-------------------------  Loop over subcycle steps  ---------------------------------//
      double time_passed = 0.0;
      for( unsigned int i=0; i<Nsteps; i++ ){

	//make sure we're not too long in subcycling
	if( (time_passed+dt_ss) > UNIT_T ){
	  dt_ss = UNIT_T - time_passed;
	  dt_eff = dt_ss/UNIT_T;
	}

	time_passed += dt_ss;
	
	double total_photo_ionisations_step = 0.0;
	for(short int f=0; f<numFreq; f++){

	  //optical depth in this step
	  tau_step[f] = initial_tau[f]*N_HI_step*one_over_initial_N_HI;
	  //number of ionisations in this step
	  photo_ionisations_step[f] = N_in_total[f] * ( 1.0 - exp(-tau_step[f]) ) * dt_eff;//dt_ss/UNIT_T;

	  //keep track of all photo-ionisations
	  number_of_photo_ionisations[f] += photo_ionisations_step[f];

	  //total number of photoionisations in this step
	  total_photo_ionisations_step += photo_ionisations_step[f];

	  if( tau_step[f] != tau_step[f] ){
	    cerr << tau_step[f] << " " << initial_tau[f] << " " << N_HI_step << " " << i << endl;
	    cerr << (double) site.get_n_HI() << " " << UNIT_D << " " << (double) site.get_neigh_dist() << " " 
		 << UNIT_L << " " << " " << straight_correction_factor << endl;
	    MPI::COMM_WORLD.Abort(-1);
	  }
	}
	
	//number of collisional ionisations in this step
	double coll_ionisations_step = 0.0;
	if(coll_ion){
	  //number of coll ionisations is coll ion coeff * n_e * number of atoms * time step
	  coll_ionisations_step = coll_ion_coeff_HI( site.get_temperature() ) * N_HII_step * one_over_volume * N_HI_step * dt_ss;
	}

	numIonised_step = total_photo_ionisations_step + coll_ionisations_step;

	if( numIonised_step > N_HI_step){
	  cerr << " Warning, subcycle step too large: numIonised: " << numIonised_step << " N_HI_step: " << N_HI_step << endl;
	  numIonised_step = N_HI_step;
	}

	//adjust number of neutrals and ioniseds accordingly
	N_HI_step  -= numIonised_step;
	N_HII_step += numIonised_step;

	//number of recombinations is this factor times recombination coefficient
	double num_rec_factor = N_HII_step * N_HII_step * one_over_volume * site.get_clumping() * dt_ss;

	//if flag for recombinations is not set, set dx to zero
	if(!recombination){
	  num_rec_factor = 0;
	}

  //number of recombinations
  double num_rec = (rec_rad) ? num_rec_factor*recomb_coeff_HII_caseA( site.get_temperature() ) 
    : num_rec_factor*recomb_coeff_HII_caseB( site.get_temperature() );


  //number of diffuse recombination photons in this subcycle step
  double N_rec_phot_step = 0.0;
  //number of ionisations by diffuse recombination photons in this subcycle step
  double N_rec_ion_step = 0.0;
  //calculate contribution of recombination radiation to ionisation inside cell
  //calculate the number of photons that is transported to neighbours

  if( rec_rad ) { 

    //fraction of the photons that can escape
    double rec_escape = rec_rad_escape( tau_step[0]/straight_correction_factor );

    //total number of recombination photons produced in the cell
    N_rec_phot_step = num_rec_factor * ( recomb_coeff_HII_caseA( site.get_temperature() ) - recomb_coeff_HII_caseB( site.get_temperature() ) );
    //N_rec_phot_step = num_rec_factor * ( recomb_coeff_HII_caseA( 10000. ) - recomb_coeff_HII_caseB( 10000. ) );


    //the fraction of diffuse photons that is not escaping 
    N_rec_ion_step = (1.0 - rec_escape) * N_rec_phot_step;
    //these photons balance recombinations
    num_rec -= N_rec_ion_step;
    //they also add to the heating, so add to photo-ionisations in this step
    photo_ionisations_step[0] += N_rec_ion_step;

    //number of outgoing diffuse photons
    N_rec_phot_step *= rec_escape;

    //keep track of total number of recombination photons
    N_rec_phot += N_rec_phot_step;

    if(N_rec_phot_step < 0.0){
      cerr << " (" << COMM_RANK << ") Negative number of diffuse photons!" << endl;
      cerr << num_rec_factor << " " << rec_escape << " " << site.get_temperature() << " " << ( recomb_coeff_HII_caseA( site.get_temperature() ) - recomb_coeff_HII_caseB( site.get_temperature() ) ) << endl;
      MPI::COMM_WORLD.Abort(-1);
    }

  }

	if( num_rec > N_HII_step){
	  num_rec = N_HII_step;
	  cerr << " Warning, subcycle step too large: dx: " << num_rec << " N_HII_step: " << N_HII_step << endl;
	}

	number_of_recombinations += num_rec;

	//check to be sure
	if( numIonised_step > N_HI_step || num_rec > N_HII_step || num_rec < 0.0 ){

	  cerr << " (" << COMM_RANK << ") Warning: time step in subcycle step " << i << " too big! " << endl
		 << " numIonised_step: " << numIonised_step << " N_HI_step: " << N_HI_step << endl
		 << " num_rec: " << num_rec << " N_HII_step: " << N_HII_step << endl
		 << " total_photo_ionisations_step: " << total_photo_ionisations_step << endl;

	    cerr << " photo_ionisations_step: ";
	    for(short int f=0; f<numFreq; f++){
	      cerr << " " << photo_ionisations_step[f];
	    }
	    cerr << endl;
	    cerr << " N_in_total: ";
	    for(short int f=0; f<numFreq; f++){
	      cerr << " " << N_in_total[f];
	    }
	    cerr << endl;

	    //MPI::COMM_WORLD.Abort(-1);
	    N_HI_step = 0.0;
	    N_HII_step = n_H/one_over_volume;

	}

	N_HI_step += num_rec;
	N_HII_step -= num_rec;

	//update the site for temperature calculation
	double tmp = N_HI_step * num_to_dens;
	site.set_n_HI( (float) tmp );
	tmp = N_HII_step * num_to_dens;
	site.set_n_HII( (float) tmp );

	//update temperature	
	double t_heat = 0.0;
	if(heat_cool){
	  
	  t_heat = update_temperature( site, photo_ionisations_step, dt_ss);

	}else{
	  //if gas gets ionised set temperature to gasIonTemp
	  if( (double) site.get_n_HII()/n_H > 0.1){
	    site.set_temperature( gasIonTemp );
	  } else {
	    site.set_temperature( gasNeutralTemp );
	  }
	}

	//if photoionisation equilibrium is reached, no need to recalculate these
	if( fabs(numIonised_step - num_rec)/num_rec < 1e-5 && (int) i < (int) (Nsteps-2) ){

	  //number of steps left to take
	  unsigned int Nleft = Nsteps - i - 1;

    //remove extra photo-ionisations by diffuse photons to avoid counting them double
    //in equilibrium
    if(rec_rad){
      photo_ionisations_step[0] -= N_rec_ion_step ;
    }

	  //in case of heating and cooling subcycle on heating time scale
	  if(heat_cool){

	    //time counter
	    double t2 = 0.0;
	    double t_left = Nleft * dt_ss;//t_end - t_ss;

	    //while( t2 < t_end2 ){
	    while( t2 < t_left ){
	      //time scale is heating time scale
	      double dt2 = subcycle_frac*t_heat;
	      //make sure time is not too big
	      if( (t2+dt2) > t_left ){
		dt2 = t_left - t2;
	      }

	      //set equilibrium fractions at this temperature
	      //assume photoionisation rate is constant

	      //total number of atoms
	      double N_H = n_H * site.get_volume() * UNIT_V;
	      //photo-ionisation rate
	      double Gamma_H = total_photo_ionisations_step/(dt_ss * N_HI_step);
	      //collisional ionisation rate
	      double Gamma_c = 0.0;
	      if(coll_ion){
		Gamma_c = coll_ion_coeff_HI( site.get_temperature() ) * (double) site.get_n_HII() * UNIT_D;
	      }
	      //total ionisation rate
	      double Gamma = Gamma_H + Gamma_c;
	      //electron density times recombination coefficient
	      double n_e_times_alpha=0.0;
  
        if(rec_rad){
        //for now assume that in this case the optical depth is low, so recombination photons almost all
        //escape. This makes it somewhat safe to use alpha_A in calculation of the equilibrium
        //Error in number of diffuse recombination photons < 0.1%
        //If case B would be used: Error > 50%
        //this error is comparable to error without diffuse photons!
          n_e_times_alpha = site.get_clumping()*recomb_coeff_HII_caseA( site.get_temperature() ) * site.get_n_HII() * UNIT_D;
        }else{
          n_e_times_alpha = site.get_clumping()*recomb_coeff_HII_caseB( site.get_temperature() ) * site.get_n_HII() * UNIT_D;
        }
                
	      //equilibrium number of neutral and ionised atoms
	      N_HI_step = n_e_times_alpha * N_H/( n_e_times_alpha + Gamma);
	      N_HII_step = N_H - N_HI_step;

        double N_rec_ion_eq = 0.0;
            //number of recombination photons
        if( rec_rad ) { 

        //optical depth through the cell
          double tau = N_HI_step * one_over_volume * (double) site.get_neigh_dist() * UNIT_L * cross_H[0] * straight_correction_factor;
        //fraction of the photons that can escape
          double rec_escape = rec_rad_escape( tau );

        //total number of recombination photons produced in the cell, error < 0.1% 
          double N_rec_phot_eq = num_rec_factor * 
            ( recomb_coeff_HII_caseA( site.get_temperature() ) - recomb_coeff_HII_caseB( site.get_temperature() ) );

        //keep track of total number of recombination photons
          N_rec_phot += rec_escape * N_rec_phot_eq*dt2/dt_ss;

        //number of ionisations by recombination photons
          N_rec_ion_eq = (1.0 - rec_escape)*N_rec_phot_eq*dt2/dt_ss;
        }
        
	      //update the site for temperature calculation
	      double tmp = N_HI_step * num_to_dens;
	      site.set_n_HI( (float) tmp );
	      tmp = N_HII_step * num_to_dens;
	      site.set_n_HII( (float) tmp );

	      //keep track of photoionisations
	      //photo_ionisations_step needs to stay constant, so make new vector
	      vector<double> photo_ionisations_heating_step(numFreq,0.0);
	      for(short int f=0; f<numFreq; f++){
		photo_ionisations_heating_step[f] = photo_ionisations_step[f]*dt2/dt_ss;
		//total number of photo-ionisations
		number_of_photo_ionisations[f] += photo_ionisations_heating_step[f];
	      }

        //add recombination photons that add to the heating
        if(rec_rad){
          photo_ionisations_heating_step[0] += N_rec_ion_eq;
        }
                
	      //calculate new temperature	    
	      t_heat = update_temperature( site, photo_ionisations_step, dt2);

	      //this is no longer needed, so clear it
	      photo_ionisations_heating_step.clear();

	      t2 += dt2;
	    }
	  }else{

	    //keep track of photoionisations
	    for(short int f=0; f<numFreq; f++){
	      number_of_photo_ionisations[f] += photo_ionisations_step[f]*Nleft;//*t_left/dt_ss;//*Nleft;
	    }

      //if diffuse recombination radiation is included, keep track of photons that need to be sent
      if( rec_rad ) { 
        //double total_recombinations = num_rec*Nleft;//*t_left/dt_ss;//Nleft;     
        N_rec_phot += N_rec_phot_step*Nleft;

      }//if diffuse
	  }

	  i = Nsteps-1;

	}//if photoionsation equilibrium
      }//for all steps

      count++;

      //calculate the number of outgoing photons after the subcycling
      for(short int f=0; f<numFreq; f++){

  //make sure N_out is bigger than zero, when subcycle time step is small this might 
  //not be the case due to numerical errors
        N_out_total[f] = ( (N_in_total[f] - number_of_photo_ionisations[f]) > 0.0 ) ? 
          (N_in_total[f] - number_of_photo_ionisations[f]) : 0.0;

  //check for NANs
        if(N_out_total[f] != N_out_total[f] ){
          cerr << endl << " (" << COMM_RANK << ") NAN detected! N_out_total: " << N_out_total[f] << " " 
            << N_in_total[f] << " " << number_of_photo_ionisations[f] << endl;
          MPI::COMM_WORLD.Abort(-1);
        }
      }

      //send the recombination photons
      if(N_rec_phot > 0.0){
        vector<double> N_out(numFreq,0.0);
        N_out[0] = N_rec_phot/(UNIT_I);
        diffuse_transport( site, N_out );
      }

    }else{ //don't use temporal photon conservation scheme

      //make sure number of ionisations is less than number of neutral atoms. Only necessary if subcycle step is too large
      if( initial_numIonised > initial_N_HI){
        initial_numIonised = initial_N_HI;
      }

      //adjust number of neutrals and ioniseds accordingly
      initial_N_HI  -= initial_numIonised;
      initial_N_HII += initial_numIonised;

      double tmp = initial_N_HI/( UNIT_D * site.get_volume() * UNIT_V );
      site.set_n_HI( (float) tmp );
      tmp = initial_N_HII/( UNIT_D * site.get_volume() * UNIT_V );
      site.set_n_HII( (float) tmp );

      //total ionisations
      double tot_ion = initial_photo_ionisations + initial_coll_ionisations;
      //fraction of ionisations caused by photo-ionisations
      double phot_ion_frac = (tot_ion > 0.0) ? initial_photo_ionisations /tot_ion : 0.0;
      double one_over_init_phot_ion = (initial_photo_ionisations > 0.0) ? 1.0/initial_photo_ionisations : 0.0;
      for(short int f=0; f<numFreq; f++){
        N_out_total[f] = N_in_total[f] - (initial_numIonised * phot_ion_frac * initial_N_retained_H[f] * one_over_init_phot_ion);

        if(N_out_total[f] != N_out_total[f] ){
          cerr << endl << " (" << COMM_RANK << ") NAN detected! N_out_total: " << N_out_total[f] << " " << N_in_total[f] << " " << initial_N_retained_H[f] << endl
	       << initial_numIonised << " " << phot_ion_frac << " " << one_over_init_phot_ion << endl;
          MPI::COMM_WORLD.Abort(-1);
        }
      }

      //set the correct temperature of the gas
      double t_heat = 0.0;
      if(heat_cool){
        t_heat = update_temperature( site, initial_N_retained_H, UNIT_T);
	//if(site.get_temperature() > 100)
	//cerr << " New temperature is: " << site.get_temperature() << " K" << endl;
      }else{
	//if gas gets ionised set temperature to gasIonTemp
        if( (double) site.get_n_HII()/n_H > 0.1){
          site.set_temperature( gasIonTemp );
	} else {
	  site.set_temperature( gasNeutralTemp );
	}
      }
    }
  }

  //convert back to numerical units
  for(short int f=0; f<numFreq; f++){
    N_in_total[f] /= UNIT_I;
    N_out_total[f] /= UNIT_I;
  }


  return N_out_total;

}


void SimpleX::source_transport( Site& site ){

  //Before transporting, source should ionise its own cell!

  //in case of ballistic transport, intensity has size of number of neighbours
  if( site.get_ballistic() ){
    for( unsigned int j=0; j<site.get_numNeigh(); j++ ){
      //loop over frequencies
      for(short int f=0; f<numFreq; f++){
        double inten = (double) site.get_flux(f) * UNIT_T/site.get_numNeigh();
        sites[ site.get_neighId(j) ].addRadiationDiffOut( f, site.get_outgoing(j), (float) inten ); 
      }
    }
    //total_inten += (double) it->get_flux() * UNIT_T * UNIT_I;

  }else{
    //in case of direction conserving transport and combined transport, intensity has 
    //the size of the tesselation of the unit sphere
    numPixels = number_of_directions;
    for( unsigned int m=0; m<numPixels; m++ ){
      //unsigned int m = 0;
      //double inten = (double) it->get_flux()* UNIT_T;

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
	    //loop over frequencies
            for(short int f=0; f<numFreq; f++){
              double inten = (double) site.get_flux(f) * UNIT_T/numPixels;
              sites[neigh].addRadiationDiffOut( f, j, (float) inten );
            }
          } 
        }
        if(!found){
          cerr << " error in source radiation" << endl;
          MPI::COMM_WORLD.Abort(-1);
        }
      }else{
	//if the neighbour is not ballistic, simply 
	//add the radiation in the same direction bin
        for(short int f=0; f<numFreq; f++){
          double inten = (double) site.get_flux(f) * UNIT_T/numPixels;
          sites[ neigh ].addRadiationDiffOut( f, m, (float) inten ); 
        }
      }

    }//for all pixels

  }//if ballistic
}


void SimpleX::diffuse_transport( Site& site, vector<double>& N_out_total ){

  //in case of ballistic transport, intensity has size of number of neighbours
  if( site.get_ballistic() ){

    //loop over all neighbours
    for( unsigned int j=0; j<site.get_numNeigh(); j++ ){

      //neighbour id
      unsigned int neigh = site.get_neighId( j );
      if( sites[ neigh ].get_ballistic() ){

	//loop over frequencies
        for(short int f=0; f<numFreq; f++){
          double inten = N_out_total[f]/site.get_numNeigh();
          sites[ neigh ].addRadiationDiffIn( f, site.get_outgoing(j), (float) inten ); 
        }

      }else{

	//if the neighbour is not ballistic, find out in which direction 
	//bins the photons should go
        vector< unsigned int > dir_to_use;
	//this trick only works for neighbours on this proc,
	//otherwise outgoing array does not exist
        if( sites[ neigh ].get_process() == COMM_RANK ){
          for( unsigned int n=0; n<number_of_directions; n++ ){
            if( sites[ neigh ].get_outgoing(n) == site.get_outgoing(j) ){
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
          for(short int f=0; f<numFreq;f++){
	    //intensity to send out
            double inten = N_out_total[f]/site.get_numNeigh();
            sites[ neigh ].addRadiationDiffIn( f,  dir_to_use[n], (float) inten/dir_to_use.size() );
          }
        }//for all neighbours to use
        dir_to_use.clear();

      }

    }

  }else{
    //in case of direction conserving transport and combined transport, intensity has 
    //the size of the tesselation of the unit sphere
    numPixels = number_of_directions;
    for( unsigned int m=0; m<numPixels; m++ ){
      //unsigned int m = 0;
      //double inten = (double) it->get_flux()* UNIT_T;

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
	    //loop over frequencies
            for(short int f=0; f<numFreq; f++){
	      //double inten = (double) site.get_flux(f) * UNIT_T/numPixels;
              double inten = (double) N_out_total[f]/numPixels;
              sites[neigh].addRadiationDiffIn( f, j, (float) inten );
            }
          } 
        }
        if(!found){
          cerr << " error in source radiation" << endl;
          MPI::COMM_WORLD.Abort(-1);
        }
      }else{
	//if the neighbour is not ballistic, simply 
	//add the radiation in the same direction bin
        for(short int f=0; f<numFreq; f++){
	  //double inten = (double) site.get_flux(f) * UNIT_T/numPixels;
          double inten = (double) N_out_total[f]/numPixels;
          sites[ neigh ].addRadiationDiffIn( f, m, (float) inten ); 
        }
      }

    }//for all pixels
  }//if ballistic

}


//redistribute the photons to the d most straightforward neighbours
void SimpleX::non_diffuse_transport( Site& site, vector<double>& N_out_total ) { 

  //in case of ballistic transport, intensity has size of number of neighbours;
  //in case of direction conserving transport, intensity has 
  //the size of the tesselation of the unit sphere
  numPixels = ( site.get_ballistic() )? site.get_numNeigh() : number_of_directions;

  //determine N_in_total
  vector<double> N_in_total(numFreq, 0.0);
  vector<bool> inten_to_send(numPixels, 0);
  for( unsigned int j=0; j < numPixels; j++) {
    for(short int f=0; f<numFreq;f++){
      if( site.get_intensityOut( f,j ) > 0.0 ) { 
        N_in_total[f] += (double) site.get_intensityOut( f,j );
        inten_to_send[j] = 1;
      }
    } //if intensityOut
  } //for all neighbours

  //loop over all neighbours/directions
  for( unsigned int j=0; j < numPixels; j++) {
    //only include directions/neighbours that have intensity to send
    if( inten_to_send[j] ) { 

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

	  //if(!sites[ neigh ].get_border() ){
	  //check if neighbour is ballistic
          if( sites[ neigh ].get_ballistic() ){

	    //send the photons to the neighbour in the place in the neighbour vector
	    //pointing to this site
            for(short int f=0; f<numFreq;f++){
	      //intensity to send out, get correct part from all added intensities
              double inten = (N_in_total[f] > 0.0) ? N_out_total[f] * ( (double) site.get_intensityOut( f,j ) / N_in_total[f] ) : 0.0;
              inten /= double(site.get_numStraight( j ));
              sites[ neigh ].addRadiationDiffIn( f, neighIdLoc, (float) inten );
            }

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
              for(short int f=0; f<numFreq;f++){
		//intensity to send out, get correct part from all added intensities
                double inten = (N_in_total[f] > 0.0) ? N_out_total[f] * ( (double) site.get_intensityOut( f,j ) / N_in_total[f] ) : 0.0;
                inten /= double(site.get_numStraight( j ));
                sites[ neigh ].addRadiationDiffIn( f,  dir_to_use[n], (float) inten/dir_to_use.size() );
              }
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

	    //check if neighbour is ballistic
            if(sites[ neigh ].get_ballistic() ){

	      //if the neighbour is ballistic, find out which 
	      //Delaunay line is associated with this direction
              bool found = 0;
              for( unsigned int n=0; !found && n<sites[ neigh ].get_numNeigh(); n++ ){
                if( sites[neigh].get_neighId(n) == site.get_site_id() ){
                  found = 1;
                  for(short int f=0; f<numFreq;f++){
                    double inten = (N_in_total[f] > 0.0) ? N_out_total[f] * ( (double) site.get_intensityOut( f,j ) / N_in_total[f] ) : 0.0;
                    inten /= double(site.get_numStraight( site.get_outgoing(j) ));
                    sites[neigh].addRadiationDiffIn( f, n, (float) inten );
                  }//for all freq
                } 
              }
              if(!found){
                cerr << " Error in photon transport: neighbour not found! " << endl;
              }
            }else{

	      //if not ballistic,
	      //send the photons to the neighbour in the same direction bin
              for(short int f=0; f<numFreq;f++){
                double inten = (N_in_total[f] > 0.0) ? N_out_total[f] * ( (double) site.get_intensityOut( f,j ) / N_in_total[f] ) : 0.0;
                inten /= double(site.get_numStraight( site.get_outgoing(j) ));
                sites[neigh].addRadiationDiffIn( f, j, (float) inten );
              }//for all freq

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
                  for(short int f=0; f<numFreq;f++){
                    double inten = (N_in_total[f] > 0.0) ? N_out_total[f] * ( (double) site.get_intensityOut( f,j ) / N_in_total[f] ) : 0.0;
                    inten /= double(num_straight);
                    sites[neigh].addRadiationDiffIn( f, n, (float) inten );
                  }
                } 
              }
              if(!found){
                cerr << " Error in photon transport: neighbour not found! " << endl;
              }
            }else{

	      //if not ballistic,
	      //send the photons to the neighbour in the same direction bin
              for(short int f=0; f<numFreq;f++){
                double inten = (N_in_total[f] > 0.0) ? N_out_total[f] * ( (double) site.get_intensityOut( f,j ) / N_in_total[f] ) : 0.0;
                inten /= double(num_straight);
                sites[neigh].addRadiationDiffIn( f, j, (float) inten );
              }
            }//if neighbour is ballistic
          }//for all straight
        }//if straight from tesselation
      }//if this site ballistic
    } //if intensityOut
  } //for all neighbours

  //loop over all neighbours/directions
  for( unsigned int j=0; j < numPixels; j++) {
    for(short int f=0; f<numFreq;f++){
      //intensity has been send away, so delete it
      site.set_intensityOut(f,j,0.0);
    }
  }

  numPixels = ( site.get_ballistic() ) ? site.get_numNeigh() : number_of_directions;
  for( unsigned int j=0; j<numPixels; j++ ) { 
    for(short int f=0; f<numFreq; f++){
      if( site.get_intensityIn( f, j ) != site.get_intensityIn( f, j )){
        cerr << " NAN detected in intensityIn " << site.get_vertex_id() << " " << f << endl;
        MPI::COMM_WORLD.Abort(-1);
      }
    }
  }


} 



void SimpleX::radiation_transport( const unsigned int& run ){

  double t0 = MPI::Wtime();
  //int dummy=0;
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
        for(short int f=0; f<numFreq; f++){
          it->set_intensityOut( f, j, 0.0 );
          it->set_intensityIn( f, j, 0.0 );
        }
      }
    }
  }//if run

  // if( COMM_RANK == 0)
  //   cerr << endl; 

  for( unsigned int times=0; times<numSweeps; times++ ) { 

    // if( COMM_RANK == 0 ){
    //   if( (int)floor( ( 100.0*times)/numSweeps ) > dummy ) { 
    //     dummy = (int) floor((100.0*times)/numSweeps); 
    //     cerr << "\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\t" << dummy+1 << "% completed..." << flush; 
    //   } 
    // }

    //Do all sources
    //double total_inten = 0.0;
    for( SITE_ITERATOR it=sites.begin(); it!=sites.end(); it++ ){
      if( it->get_source() ){
        source_transport( *it );
      }//if flux
    }//for all sites

    //double source_flux;
    //MPI::COMM_WORLD.Allreduce(&total_inten,&source_flux,1,MPI::DOUBLE,MPI::SUM);

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

    // if( COMM_RANK == 0 ){
    //   cerr << " (" << COMM_RANK << ") Total number of photons sent: " << source_flux << " total recombinations: " << recombinations << endl;
    // }

    //redistribute the photons
    for( SITE_ITERATOR it=sites.begin(); it!=sites.end(); it++ ){
      if( !it->get_border() && it->get_process() == COMM_RANK ) {

        //solve the rate equation to determine ionisations and recombinations
        vector<double> N_out = solve_rate_equation( *it );

        if(!diffuseTransport){
          //transport photons with ballistic transport or DCT
          non_diffuse_transport( *it, N_out );
        } else {
          //transport photons diffusely
          diffuse_transport( *it, N_out );
        }

      }//if not in border   
    }//for all sites

    double total_diffuse_proc=0.0;
    for( SITE_ITERATOR it=sites.begin(); it!=sites.end(); it++ ){

      //in case of ballistic transport, intensity has size of number of neighbours;
      //in case of direction conserving transport, intensity has 
      //the size of the tesselation of the unit sphere
      numPixels = ( it->get_ballistic() ) ? it->get_numNeigh() : number_of_directions;

      for( unsigned int j=0; j<numPixels; j++ ) { 
        for(short int f=0; f<numFreq; f++){
	        //double inten = (double) it->get_intensityIn(j) + (double) it->get_intensityOut(j);
	        //float inten = it->get_intensityIn(f,j);
          it->set_intensityOut( f, j, it->get_intensityIn(f,j) );
          it->set_intensityIn( f, j, 0.0 );

        }
        //total_diffuse_proc += (double) it->get_intensityOut(0,j);
      }//for all pixels
    }//for all sites

    // double total_diffuse;
    // MPI::COMM_WORLD.Allreduce(&total_diffuse_proc,&total_diffuse,1,MPI::DOUBLE,MPI::SUM);
    // 
    // if(COMM_RANK == 0){
    //   cerr << " Number of diffuse photons: " << total_diffuse*UNIT_I << endl;
    // }

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

    double IFront_pos = 0.0;
    if(give_IFront){
      IFront_pos = calc_IFront( run );
    }

    float time = (run+1)*simTime*time_conversion/(numRuns*secondsPerMyr);

    char fileName[100][256];
    const char* outputBaseName1 = "SimpleX_data";

    //write output per processor to avoid disk waits
    for(unsigned int i=0; i<COMM_SIZE; i++){
      if(COMM_RANK == i){


        sprintf(fileName[0], "%s_%f_%d.hdf5", outputBaseName1, time, COMM_RANK);

        write_hdf5_output(fileName[0], run);
      }
      MPI::COMM_WORLD.Barrier();
    }

    double t1 = MPI::Wtime();
    if( COMM_RANK == 0 ){
      simpleXlog << endl << "  Generating output took " << t1-t0 << " seconds" << endl;
    }

    if( COMM_RANK == 0 ){

      if(give_IFront){
        cerr << "  Position of the I-Front            : " << IFront_pos << " cm" << endl;
        cerr << endl; 
      }


      simpleXlog << endl << endl << " END OF RUN " << run+1 << endl << endl;
      if(give_IFront){
        simpleXlog << "  Position of the I-Front            : " << IFront_pos << " cm"<< endl;
      }
      simpleXlog << "  Output written to: " << fileName[0] << endl;

    }
  }//if output in this run

  return;
}


//Write hdf5 output. Difference for serial version is the proc dependence 
//of the offset
void SimpleX::write_hdf5_output(char *name, const unsigned int& run){

  //determine sites local to proc
  vector< unsigned long int > local_sites;
  //vector< unsigned long int > indices(vertex_id_max,vertex_id_max);
  for( unsigned long int i=0; i<sites.size(); i++ ){
    if( sites[i].get_process() == COMM_RANK ){
      local_sites.push_back( i );
      //indices[sites[i].get_vertex_id()]=i;
    }
  }


  //simplices are local always, determine total number of simplices
  //across procs

  arr_1D<double> double_arr;
  arr_1D<unsigned int> int_arr;

  unsigned long long int dims[2];
  unsigned long long int offset[2];

  //open file  
  h5w file(name,'n');

  // write structure of the file
  file.make_group("/Header");
  file.make_group("/Vertices");
  file.make_group("/Simplices");


  //Vertex attributes
  dims[0] = local_sites.size();

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
  file.make_dataset("/Vertices/temperature","double",1,dims);
  file.write_attr("/Vertices/temperature","var_name","temperature");

  //Simplex attributes
  dims[0] = simplices.size();

  dims[1] = 4;                              
  file.make_dataset("/Simplices/indices","unsigned int",2,dims);
  file.write_attr("/Simplices/indices","var_name","indices");
  dims[1] = 1;
  file.make_dataset("/Simplices/volume","double",1,dims);
  file.write_attr("/Simplices/volume","var_name","vol");



  //=============== write header ======================================//


  file.write_attr("/Header","local_number_of_sites",(unsigned int) local_sites.size() );
  file.write_attr("/Header","local_number_of_simplices", (unsigned int) simplices.size() );
  //file.write_attr("/Header","number_of_sites", (unsigned int) sites.size());
  file.write_attr("/Header","number_of_sites", (unsigned int)numSites);
  file.write_attr("/Header","number_of_simplices", (unsigned int) simplices.size());


  file.write_attr("/Header","number_of_sweeps", numSweeps);
  file.write_attr("/Header","box_size_in_pc", sizeBox);
  //file.write_attr("/Header","redshift", redShift);
  double this_time = (run+1)*simTime/(numRuns*secondsPerMyr);
  file.write_attr("/Header","current_time_in_Myr", this_time);
  double total_time = simTime/secondsPerMyr;
  file.write_attr("/Header","total_simulation_time_in_Myr", total_time);


  if (recombination){
    file.write_attr("/Header","recombination_on","true");
  }else{
    file.write_attr("/Header","recombination_on","false");
  }


  //======================= Write Vertices ===========================//

  // writing is done one variable at a time, through the arr_1D instance of size chunk_size
  // its not equal to numSites in order to preserve memory....

  unsigned int chunk_size_vert;

  //writing in chunks isn't implemented in parallel
  chunk_size_vert = chunk_size;
  if (chunk_size_vert > local_sites.size())
    chunk_size_vert = local_sites.size();

  // do the writing !!!
  offset[1] = 0;
  for( unsigned long int i=0; i<local_sites.size(); i+=chunk_size_vert ){


    offset[0] = i;

    if (i+chunk_size_vert >= local_sites.size()) // make sure not to write outside of data range
      dims[0] = (int) local_sites.size()-i;
    else
      dims[0] = chunk_size_vert;

    // writing coordinates
    dims[1] = 3; 

    double_arr.reinit(2,dims);
    for(unsigned int j=0; j<dims[0]; j++ ){
      unsigned long int index = local_sites[i+j];
      double_arr(j,0) = sites[index].get_x();
      double_arr(j,1) = sites[index].get_y();
      double_arr(j,2) = sites[index].get_z();

    }
    file.write_data("/Vertices/coordinates",offset, &double_arr);


    dims[1]=1;
    //write id
    int_arr.reinit(1,dims);
    for(unsigned int j=0; j<dims[0]; j++ )
      int_arr(j) = sites[ local_sites[i+j] ].get_vertex_id();
    file.write_data("/Vertices/vertex_index",offset, &int_arr);

    //write border
    int_arr.reinit(1,dims);
    for(unsigned int j=0; j<dims[0]; j++ )
      int_arr(j) = (unsigned int) sites[ local_sites[i+j] ].get_border();
    file.write_data("/Vertices/border",offset, &int_arr);

    // writing neutral fraction      
    double_arr.reinit(1,dims);
    for(unsigned int j=0; j<dims[0]; j++ ){
      double neutr_frac = (double) sites[ local_sites[i+j] ].get_n_HI()/( (double) sites[ local_sites[i+j] ].get_n_HI() + (double) sites[ local_sites[i+j] ].get_n_HII() );
      double_arr(j) = neutr_frac;
    }
    file.write_data("/Vertices/H_neutral_fraction",offset, &double_arr);

    // writing number density of gas

    double_arr.reinit(1,dims);
    for(unsigned int j=0; j<dims[0]; j++ )
      double_arr(j) = ( (double) sites[ local_sites[i+j] ].get_n_HI() + (double) sites[ local_sites[i+j] ].get_n_HII() ) * UNIT_D;
    file.write_data("/Vertices/number_density",offset, &double_arr);

    // writing cell volumes
    double_arr.reinit(1,dims);
    for(unsigned int j=0; j<dims[0]; j++ ) 
      double_arr(j) = (double) sites[ local_sites[i+j] ].get_volume() * UNIT_V;
    file.write_data("/Vertices/volume",offset, &double_arr);

    // writing Luminositites
    double_arr.reinit(1,dims);
    for(unsigned int j=0; j<dims[0]; j++ ){
      double flux = 0.0;
      if(sites[ local_sites[i+j] ].get_source()){
        for(short int f=0; f<numFreq; f++){
          flux += (double) sites[ local_sites[i+j] ].get_flux(f) * UNIT_I;
        }
      }
      double_arr(j) = flux;
    }
    file.write_data("/Vertices/luminosity",offset, &double_arr);

    // writing temperature
    double_arr.reinit(1,dims);
    for(unsigned int j=0; j<dims[0]; j++ )
      double_arr(j) = (double) sites[ local_sites[i+j] ].get_temperature();
    file.write_data("/Vertices/temperature",offset, &double_arr);

  }//for all local sites


  local_sites.clear();

  //====================== Write Simplices =================================//

  unsigned int chunk_size_simpl = chunk_size;

  //writing in chunks isn't implemented in parallel

  chunk_size_simpl = chunk_size;
  if (chunk_size_simpl > simplices.size())
    chunk_size_simpl = simplices.size();


  // do the writing !!!
  offset[1] = 0;
  for (unsigned long int i=0; i<simplices.size(); i+=chunk_size_simpl){

    offset[0] = (int) i;

    if ( (i + chunk_size_simpl) >= simplices.size() ) // make sure not to write outside of data range
      dims[0] = (int) simplices.size() - i;
    else
      dims[0] = chunk_size_simpl;

    // writing indices
    dims[1] = 4; 

    int_arr.reinit(2,dims);
    for(unsigned int j=0; j<dims[0]; j++ ){
      int_arr(j,0) = simplices[i+j].get_id1();
      int_arr(j,1) = simplices[i+j].get_id2();
      int_arr(j,2) = simplices[i+j].get_id3();
      int_arr(j,3) = simplices[i+j].get_id4();
    }
    file.write_data("/Simplices/indices",offset, &int_arr);

    dims[1]=1;

    // writing simplex volumes
    double_arr.reinit(1,dims);
    for(unsigned int j=0; j<dims[0]; j++ ) 
      double_arr(j) = simplices[i+j].get_volume();

    file.write_data("/Simplices/volume",offset, &double_arr);

  }

  dims[0]=0;
  dims[1]=0;

  file.close();

}

//write position of the IFront to file
double SimpleX::calc_IFront( const unsigned int& run ){

  int teller=0; 
  int tellerTotal;
  double averDist=0.0; 
  double averDistTotal;
  double source_x = 0.5;
  double source_y = 0.5;
  double source_z = 0.5;

  //first update the local copies of the objects
  //update_sites();

  for( SITE_ITERATOR it=sites.begin(); it!=sites.end(); it++ ){

    if( !it->get_border() && it->get_process() == COMM_RANK ){

      double frac1 = (double) it->get_n_HII()/( (double) it->get_n_HI() + (double) it->get_n_HII() );

      //if( frac1 >= 0.49 && frac1 <= 0.51 ){  
      if( frac1 >= 0.45 && frac1 <= 0.55 ){  

        averDist += sqrt( pow( (double) it->get_x() - (double) source_x, 2 ) + 
			  pow( (double) it->get_y() - (double) source_y, 2 ) + 
			  pow( (double) it->get_z() - (double) source_z, 2 ) );
        teller++;  

      } //if 
    }  //if
  }  //for all vertices

  MPI::COMM_WORLD.Allreduce(&averDist,&averDistTotal,1,MPI::DOUBLE,MPI::SUM);
  MPI::COMM_WORLD.Allreduce(&teller,&tellerTotal,1,MPI::INT,MPI::SUM);

  if( COMM_RANK == 0 ){
    float time = (run+1)*simTime/(numRuns*secondsPerMyr);
    if( tellerTotal ){ 
      IFront_output << time << " " << ( averDistTotal/tellerTotal ) * UNIT_L << endl; 
    } else { 
      IFront_output << time << " " << 0.0 << endl; 
    }
  }

  if(tellerTotal) { 
    return (averDistTotal/tellerTotal) * UNIT_L; 
  } else { 
    return -1.0;
  }

  return -1.0;

}

double SimpleX::calc_escape_fraction( const unsigned int& run, const unsigned int& times ){


  double numPhotons = 0.0;
  double flux=0.0;
  for( SITE_ITERATOR it=sites.begin(); it!=sites.end(); it++ ){

    if( it->get_source() ){
      for(short int f=0; f<numFreq; f++){
        flux += it->get_flux(f)*UNIT_I;
      }
    }

    //in case of ballistic transport, intensity has size of number of neighbours;
    //in case of direction conserving transport, intensity has 
    //the size of the tesselation of the unit sphere
    numPixels = ( it->get_ballistic() ) ? it->get_numNeigh() : number_of_directions;

    //if( it->get_process() == COMM_RANK && it->get_border() ){
    if(it->get_border()){
      for( unsigned int j=0; j<numPixels; j++ ){
        for(short int f=0; f<numFreq; f++){

	  // if( ((double) it->get_intensityOut(j) * UNIT_I) != ((double) it->get_intensityOut(j) * UNIT_I)){
	  //   cerr << endl << " (" << COMM_RANK << ") intensityOut NAN!" << (double) it->get_intensityOut(j) << endl;
	  // }
	  // if( ((double) it->get_intensityIn(j) * UNIT_I) != ((double) it->get_intensityIn(j) * UNIT_I)){
	  //   cerr << endl << " (" << COMM_RANK << ") intensityIn NAN! " << (double) it->get_intensityIn(j) << endl;
	  // }
          numPhotons += (double) it->get_intensityOut(f,j) * UNIT_I;
          numPhotons += (double) it->get_intensityIn(f,j) * UNIT_I;
          it->set_intensityOut(f,j, 0.0);
          it->set_intensityIn(f,j, 0.0);
        }
      }
    }
  }//for all sites

  //MPI::COMM_WORLD.Barrier();

  if(numPhotons != numPhotons){
    cerr << " (" << COMM_RANK << ") numPhotons NAN!" << endl;
  }

  double total_numPhotons=0.0;
  //MPI::COMM_WORLD.Allreduce(&numPhotons,&total_numPhotons,1,MPI::DOUBLE,MPI::SUM);
  MPI::COMM_WORLD.Reduce(&numPhotons,&total_numPhotons,1,MPI::DOUBLE,MPI::SUM,0);
  double total_flux;
  //MPI::COMM_WORLD.Allreduce(&flux,&total_flux,1,MPI::DOUBLE,MPI::SUM);
  MPI::COMM_WORLD.Reduce(&flux,&total_flux,1,MPI::DOUBLE,MPI::SUM,0);

  double f_esc = total_numPhotons/(UNIT_T);
  if(total_flux > 0.0){
    f_esc /= total_flux;
  }else{
    f_esc = -1.0;    
  }




  // if( COMM_RANK == 0 ){
  //   cerr << " (" << COMM_RANK << ") escape fraction: " << f_esc << " total photons: " << total_numPhotons << " total flux: " << total_flux << endl;
  // }

  if( COMM_RANK == 0 ){

    double time = UNIT_T * times;
    time += (run)*simTime/numRuns;
  }

  return f_esc;

}



/********************************************************************************/
/*                             Generic Routines                                  */
/********************************************************************************/

//clear all temporary arrays
void SimpleX::clear_temporary(){

  //delete temporary
  for( SITE_ITERATOR it=sites.begin(); it!=sites.end(); it++ ){

    it->delete_flux();
    it->delete_straight();
    it->delete_neighId();
    it->delete_intensityIn(numFreq);
    it->delete_intensityOut(numFreq);
    it->delete_outgoing();

  }//for all sites

  simplices.clear();
  send_list.clear();

}

double SimpleX::count_photons(){

  //send_intensities();

  double total_numPhotons;
  double numPhotons = 0.0;
  for( SITE_ITERATOR it=sites.begin(); it!=sites.end(); it++ ){

    //in case of ballistic transport, intensity has size of number of neighbours;
    //in case of direction conserving transport, intensity has 
    //the size of the tesselation of the unit sphere
    numPixels = ( it->get_ballistic() ) ? it->get_numNeigh() : number_of_directions;

    if( it->get_process() == COMM_RANK ){
      for( unsigned int j=0; j<numPixels; j++ ){
        for(short int f=0; f<numFreq; f++){
          numPhotons += (double) it->get_intensityOut(f,j) * UNIT_I;
          numPhotons += (double) it->get_intensityIn(f,j) * UNIT_I;
        }
      }

    }else if( it->get_border() ){
      for( unsigned int j=0; j<numPixels; j++ ){
        for(short int f=0; f<numFreq; f++){
          numPhotons += (double) it->get_intensityOut(f,j) * UNIT_I;
          numPhotons += (double) it->get_intensityIn(f,j) * UNIT_I;
        }
      }
    }
  }

  MPI::COMM_WORLD.Reduce(&numPhotons,&total_numPhotons,1,MPI::DOUBLE,MPI::SUM,0);

  //   if( COMM_RANK == 0 ){
  //     cerr << " (" << COMM_RANK << ") Total number of photons: " << total_numPhotons << endl;
  //   }

  double total_numPhotons_in_domain = 0.0;
  numPhotons = 0.0;
  for( SITE_ITERATOR it=sites.begin(); it!=sites.end(); it++ ){
    if( it->get_process() == COMM_RANK && !it->get_border() ){

      //in case of ballistic transport, intensity has size of number of neighbours;
      //in case of direction conserving transport, intensity has 
      //the size of the tesselation of the unit sphere
      numPixels = ( it->get_ballistic() ) ? it->get_numNeigh() : number_of_directions;

      for( unsigned int j=0; j<numPixels; j++ ){
        for(short int f=0; f<numFreq; f++){
          numPhotons += (double) it->get_intensityOut(f,j) * UNIT_I;
          numPhotons += (double) it->get_intensityIn(f,j) * UNIT_I;
        }
      }
    }
  }

  MPI::COMM_WORLD.Reduce(&numPhotons,&total_numPhotons_in_domain,1,MPI::DOUBLE,MPI::SUM,0);

  //   if( COMM_RANK == 0 )
  //     cerr << " (" << COMM_RANK << ") Total number of photons in simulation domain: " << total_numPhotons_in_domain << endl;

  return total_numPhotons_in_domain;

}
