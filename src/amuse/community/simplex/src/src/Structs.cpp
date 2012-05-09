/*************************************************************************
file:         Sites.cpp
author:       Jan-Pieter Paardekooper
mail:         jppaarde@strw.leidenuniv.nl
version:      0.1
last change:  03.04.2008
---------------------------------------------------------------------
description:
This file contains the class that has all the information needed for the 
radiative transfer calculations at every site.
**************************************************************************/
/*
 * Date: Name
 * Put additional comments here
 *
 */

/***** To Do *******
 *
 *
 *******************/

#include "Structs.h"

using namespace std;


Vertex::Vertex(){

  border = 0;
  process = 0;
  vertex_id = (unsigned long long) 0;
  x = 0.0;
  y = 0.0;
  z = 0.0;

}

Vertex& Vertex::operator=(const Site& p2){

  x         = p2.x;
  y         = p2.y;
  z         = p2.z;
  process   = p2.process;
  vertex_id = p2.vertex_id;
  border    = p2.border;  

  return *this;

}

MPI::Datatype Vertex::MPI_Type;
void Vertex::construct_datatype(){

  const unsigned int entries = 6;
  Vertex dummy[1];

  int block[entries] = {1,1,1,1,1,1};

  MPI::Aint displ[entries];
  displ[0] = MPI::Get_address( &dummy[0].x );
  displ[1] = MPI::Get_address( &dummy[0].y );
  displ[2] = MPI::Get_address( &dummy[0].z );
  displ[3] = MPI::Get_address( &dummy[0].process );
  displ[4] = MPI::Get_address( &dummy[0].vertex_id );
  displ[5] = MPI::Get_address( &dummy[0].border );

  MPI::Aint base = displ[0];
  for(unsigned int i=0; i<entries; i++ ){
    displ[i] -= base;
  }

  MPI::Datatype types[entries] = { MPI::FLOAT, MPI::FLOAT, MPI::FLOAT, 
				   MPI::UNSIGNED, MPI::UNSIGNED, MPI::BOOL };

  MPI_Type = MPI_Type.Create_struct( entries, block, displ, types );
  MPI_Type.Commit();

}

//! constructor for cooling_curve object
cooling_curve::~cooling_curve(){
  
  Cooling.clear();

  return;
}

cooling_curve::cooling_curve(string fileName){

  // Open file
  ifstream ifs(fileName.c_str(),ifstream::in);
  if(ifs.is_open()){

    // Read first line (number of elements)
    ifs >> elements;

    // Temporary variables
    double T,C;

    // Fill vectors
    while(ifs.good()){
      ifs >> T >> C;
      Temperature.push_back(T);
      Cooling.push_back(C);
    }
    ifs.close();

    // Last entry is repetition
    Temperature.pop_back();
    Cooling.pop_back();

  } else {
    
    cerr << "Error reading file " << fileName << endl;
  }

  // Set the min/max and delta values
  max = Temperature.back();
  min = Temperature.front();
  delta = max - min;

  // Clear the temperature vector
  Temperature.clear();

  return;

}

double cooling_curve::get_value(double T){

  double logT = log10(T);

  if(logT>max)
    return pow(10.0, Cooling.back());// return the last entry

  if(logT<min)
    return pow(10.0, Cooling.front());// return the first entry

  double pos = ((logT-min)/delta)*double(elements-1);
  int posInt = floor(pos);
  double diff = pos - double(posInt);
  double value = Cooling[posInt]*(1-diff)+diff*Cooling[posInt+1];
  
  return pow(10.0,value);

}



//! constructor for Simpl array
Simpl::Simpl(){

  id1=0;
  id2=0;
  id3=0;
  id4=0;
  volume=0.0;

} 

//! Constructor for Sites array
Site::Site(){ 

    numNeigh = 0;
    border = 0;
    volume = 0.0;
    neigh_dist = 0.0;   
    n_HI = 0.0;
    n_HII = 0.0;
    ballistic = 1;
    internalEnergy = 0.0;
    dinternalEnergydt = 0.0;
    clumping = 1.0;
    metallicity=0.0;
    source = 0;

    flux = NULL; 
    neighId = NULL;
    intensityIn = NULL;
    intensityOut = NULL;
    outgoing = NULL;  

}

Site::~Site(){


}

//! Constructor for Site_Update
Site_Update::Site_Update(){

    vertex_id = 0;   
    site_id = 0;   
    n_HI = 0.0;   
    n_HII = 0.0; 
    internalEnergy = 0.0;
    dinternalEnergydt = 0.0;
    //flux = 0.0;            
    ballistic = 0;
  
}


Site& Site::operator=(const Vertex& p2){

  x         = p2.x;
  y         = p2.y;
  z         = p2.z;
  process   = p2.process;
  vertex_id = p2.vertex_id;
  border    = p2.border;

  return *this;

}

Site& Site::operator=(const Site_Update& p2){

  vertex_id      = p2.vertex_id;
  site_id        = p2.site_id;
  process        = p2.process;
  n_HI           = p2.n_HI;
  n_HII          = p2.n_HII;
  internalEnergy = p2.internalEnergy;
  dinternalEnergydt = p2.dinternalEnergydt;
  //flux         = p2.flux;
  ballistic      = (bool) p2.ballistic;

  return *this;

}

void Site::addRadiationDiffIn(const short int& f, const unsigned int& id, const double& intensity) {

  double intensity_old = (double) intensityIn[f][id];
  double intensity_new = intensity_old + intensity;
  intensityIn[f][id] = (float) intensity_new;

}

void Site::addRadiationDiffOut(const short int& f, const unsigned int& id, const double& intensity) {

  double intensity_old = (double) intensityOut[f][id];
  double intensity_new = intensity_old + intensity;
  intensityOut[f][id] = (float) intensity_new;

}

Send_Intensity::Send_Intensity(){

    process = 0;
    neighId = 0;
    id = 0;
    intensityIn = 0.0;
    intensityOut = 0.0;

}

Site_Update& Site_Update::operator=(const Site& p2){

  process        = p2.process;
  vertex_id      = p2.vertex_id;
  site_id        = p2.site_id;
  n_HI           = p2.n_HI;
  n_HII          = p2.n_HII;
  internalEnergy = p2.internalEnergy;
  dinternalEnergydt = p2.dinternalEnergydt;
  //flux           = p2.flux;
  ballistic      = (unsigned int) p2.ballistic;

  return *this;

}

//!compare Vertex class by vertex id for sort routine
bool compare_vertex_id_vertex( const Vertex& x, const Vertex& y){
  return x.get_vertex_id() < y.get_vertex_id();
}

bool compare_vertex_id_site(const Site& x, const Site& y){
  return x.get_vertex_id() < y.get_vertex_id();
}

//overload < operator to use sort function
bool compare_process_send_site(const Send_Site& x, const Send_Site& y){
  return x.get_process() < y.get_process();
}

//compare processes for sort
bool compare_process_send_neigh(const Send_Neigh& x, const Send_Neigh& y){
  return x.get_process() < y.get_process();
}
//compare vertex ids for sort
bool compare_vertex_id_send_neigh(const Send_Neigh& x, const Send_Neigh& y){
  return x.get_vertex_id() < y.get_vertex_id();
}


//overload < operator to use sort function
bool compare_process_send_intensity(const Send_Intensity& x, const Send_Intensity& y){
  return x.get_process() < y.get_process();
}
//overload < operator to use sort function
bool compare_index_send_intensity(const Send_Intensity& x, const Send_Intensity& y){
  return x.get_id() < y.get_id();
}

//overload < operator to use sort function
bool compare_process_site_remove(const Site_Remove& x, const Site_Remove& y){
  return x.get_process() < y.get_process();
}


bool compare_vertex_id_site_update( const Site_Update& x, const Site_Update& y){
  return x.get_vertex_id() < y.get_vertex_id();
}

//overload < operator to use sort function
bool compare_process_site_update( const Site_Update& x, const Site_Update& y){
  return x.get_process() < y.get_process();
}
