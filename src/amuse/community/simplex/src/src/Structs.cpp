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
    flux = 0.0; 
    ballistic = 1;


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
    flux = 0.0;            
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
  n_HI           = p2.n_HI;
  n_HII          = p2.n_HII;
  flux           = p2.flux;
  ballistic      = (bool) p2.ballistic;

  return *this;

}

void Site::addRadiationDiffIn(const unsigned int& id, const double& intensity) {

  double intensity_old = (double) intensityIn[id];
  double intensity_new = intensity_old + intensity;
  intensityIn[id] = (float) intensity_new;

}

void Site::addRadiationDiffOut(const unsigned int& id, const double& intensity) {

  double intensity_old = (double) intensityOut[id];
  double intensity_new = intensity_old + intensity;
  intensityOut[id] = (float) intensity_new;

}

Send_Intensity::Send_Intensity(){

    process = 0;
    neighId = 0;
    id = 0;
    intensityIn = 0.0;
    intensityOut = 0.0;

}

Site_Update& Site_Update::operator=(const Site& p2){

  vertex_id      = p2.vertex_id;
  site_id        = p2.site_id;
  n_HI           = p2.n_HI;
  n_HII          = p2.n_HII;
  flux           = p2.flux;
  ballistic      = (unsigned int) p2.ballistic;

  return *this;

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
bool compare_process_site_remove(const Site_Remove& x, const Site_Remove& y){
  return x.get_process() < y.get_process();
}


bool compare_vertex_id_site_update( const Site_Update& x, const Site_Update& y){
  return x.get_vertex_id() < y.get_vertex_id();
}
