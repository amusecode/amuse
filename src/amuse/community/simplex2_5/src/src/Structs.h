/*************************************************************************
file:         Structs.h
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

 */

/***** TO DO *****
 *
 *
 *****************/

#ifndef STRUCTS_H
#define STRUCTS_H

#include "Common.h"
#include "mpi.h"

using namespace std;

//forward declaration
class Vertex;
class Site;
class Site_Update;
class Simpl;
class SimpleX;

typedef vector< Vertex >::size_type vert_type;
typedef vector< Site >::size_type site_type;
typedef vector< Simpl >::size_type simpl_type;

//! Class that holds the information of a vertex
class Vertex{

 public:

  friend class Site;
  friend class Site_Update;

  Vertex();

  //these two should be statc
  static MPI::Datatype MPI_Type;
  static void construct_datatype();

  //! set the process the vertex lives on
  void set_process(const unsigned int& proc ){ process = proc; }
  //! get the process the vertex lives on
  const unsigned int& get_process() const{ return process; }

  //! set global id of the vertex
  void set_vertex_id( const unsigned long long int& index ){ vertex_id = index; }
  //! get global id of the vertex
  const unsigned long long int& get_vertex_id() const{ return vertex_id;  }

  //! set if vertex is in border or not
  void set_border( const unsigned int& brdr ){ border = brdr; }
  //! check whether site is in the border of the domain
  const unsigned int& get_border() const{ return border; }

  //! set x coordinate of the vertex
  void set_x( const float& x1 ){ x = x1; }
  //! get x coordinate
  const float& get_x() const{ return x; }

  //! set y coordinate of the vertex
  void set_y( const float& x2 ){ y = x2; }
  //! get x coordinate of the vertex
  const float& get_y() const{ return y; }

  //! set z coordinate of the vertex
  void set_z( const float& x3 ){ z = x3; }
  //! get z coordinate of the vertex
  const float& get_z() const{ return z; }


 private:

  unsigned int process;   //!< Process the site lives on
  unsigned long long int vertex_id;    //!< Global unique index of the vertex
  unsigned int border;            //!< Vertex is in border or not
  float x;                    //!< x-coordinate of the vertex
  float y;                    //!< y-coordinate of the vertex
  float z;                    //!< z-coordinate of the vertex

 public:
  Vertex& operator=(const Site& p2); 

};


//! Class that holds the information of a simplex

//! A simplex consists of 4 vertices in 3d
class Simpl{

 public:

  Simpl();

  //! set id of first vertex
  void set_id1( const unsigned long long int& id ){ id1 = id; }
  //! get id of first vertex
  const unsigned long long int& get_id1() const{ return id1; }

  //! set id of second vertex
  void set_id2( const unsigned long long int& id ){ id2 = id; }
  //! get id of second vertex
  const unsigned long long int& get_id2() const{ return id2; }

  //! set id of third vertex
  void set_id3( const unsigned long long int& id ){ id3 = id; }
  //! get id of third vertex
  const unsigned long long int& get_id3() const{ return id3; }

  //! set id of fourth vertex
  void set_id4( const unsigned long long int& id ){ id4 = id; }
  //! get id of fourth vertex
  const unsigned long long int& get_id4() const{ return id4; }

  //! set volume of simplex
  void set_volume( const float& vol ){ volume = vol; }
  //! get volume of simplex
  const float& get_volume() const{ return volume; }

 private:

  unsigned long long int id1;     //!< Index of first vertex
  unsigned long long int id2;     //!< Index of second vertex
  unsigned long long int id3;     //!< Index of third vertex
  unsigned long long int id4;     //!< Index of fourth vertex
  float volume;     //!< Volume of the simplex

};

//! Class that holds information needed to do radiation transport between vertices
class Site : public Vertex {

 public:

  friend class Vertex;
  friend class Site_Update;

  //!< Constructor
  Site();
  //!< Destructor
  ~Site();

  //! set id of the site

  //! This is the position of the site in the local sites array
  void set_site_id( const unsigned long long int& index ){ site_id = index; }
  //! get local id of the vertex
  const unsigned long long int& get_site_id() const{ return site_id;  }

  //! set number of neighbours of vertex
  void set_numNeigh( const unsigned char& neigh ){ numNeigh = neigh; }
  //! get number of neighbours of vertex
  const unsigned char& get_numNeigh() const{ return numNeigh; }

  //! set volume of the Voronoi cell around vertex
  void set_volume( const float& vol ){ volume = vol; }
  //! get volume of the Voronoi cell around vertex
  const float& get_volume() const{ return volume; }

  //! assign optical depth to site
  void set_neigh_dist(const float& dist){ neigh_dist = dist; }
  //! get theoptical depth to the site
  const float& get_neigh_dist() const{ return neigh_dist; }

  /* //! assign the total number density to site */
  /* void set_number_density(const float& nH){ numberDensity = nH; } */
  /* //! get the total number of atoms the site */
  /* const float& get_number_density() const{ return numberDensity; } */

  /* //! assign the ionised fraction to site */
  /* void set_ionised_fraction(const float& ionFrac){ ion_frac = ionFrac; } */
  /* //! get the ionised fraction in the site */
  /* const float& get_ionised_fraction() const{ return ion_frac; } */

  //! assign the number density of neutral hydrogen to site
  void set_n_HI(const float& nH){ n_HI = nH; }
  //! get the number density of neutral hydrogen of site
  const float& get_n_HI() const{ return n_HI; }

  //! assign the number density of ionised hydrogen to site
  void set_n_HII(const float& nH){ n_HII = nH; }
  //! get the number density of ionised hydrogen of site
  const float& get_n_HII() const{ return n_HII; }

  //! assign internalEnergy to site
  void set_internalEnergy(const float& _u){ internalEnergy = _u; }
  //! get the u in the site
  const float& get_internalEnergy() const{ return internalEnergy; }

  //! assign dudt to site
  void set_dinternalEnergydt(const float& _du){ dinternalEnergydt = _du; }
  //! get the dudt in the site
  const float& get_dinternalEnergydt() const{ return dinternalEnergydt; }

  //! assign internalEnergy from temperature
  void set_temperature(const float& _T){ internalEnergy =  (1.5 * k_B * _T)/((n_HI+n_HII)/(n_HI+2*n_HII) * m_H); }
  //! get the temperature
  float get_temperature() { return (((n_HI+n_HII)/(n_HI+2*n_HII)) * m_H * internalEnergy)/(1.5 * k_B); }

  //! assign clumping factor to site
  void set_clumping(const float& _C){ clumping = _C; }
  //! get the clumping factor of the site
  const float& get_clumping() const{ return clumping; }

  //! assign ballistic 
  void set_ballistic(const bool& ballistic_){ ballistic = ballistic_; }
  //! get ballistic
  const bool& get_ballistic() const{ return ballistic; }

  //! assign source
  void set_source(const bool& s_){ source = s_; }
  //! get ballistic
  const bool& get_source() const{ return source; }

  //!delete flux array
  void delete_flux(){
    if(flux){
      delete [] flux;
      flux = NULL;
    }
  }
  //!create flux array
  void create_flux(const short int& numFreq){ flux = new float[numFreq]; }
  //! assign flux to site
  void set_flux(const short int& f, const float& phi){ flux[f] = phi; }
  //! get the flux in the site
  const float& get_flux(const short int& f) const{ return flux[f]; }

  //! delete the neighId array
  void delete_neighId(){
    if(neighId){
      delete [] neighId;
      neighId = NULL;
    }
  }
  //! create the neighId array
  void create_neighId(){ neighId = new unsigned int[numNeigh]; }
  //! store the id's of neighbours of site
  void set_neighId(const unsigned int& j, const unsigned int& id){ neighId[j]=id; }
  //! get the id of neighbour of site
  const unsigned int& get_neighId(const unsigned int& i) const{ return neighId[i]; }

  //! create straight array
  void create_straight(){
    straight.resize( numNeigh );
  } 
  //! delete straight array
  void delete_straight(){
    straight.clear();
  }
  //! assign the id's of straightest neighbours to straight array
  void add_straight( const unsigned int& _neigh, const unsigned int& _neighNr ){ straight[ _neigh ].push_back( _neighNr);}
  //! remove id's of straightest neighbours to straight array
  void remove_straight( const unsigned int& _neigh, const unsigned int& _neighNr ){ straight[ _neigh ].erase( straight[ _neigh ].begin() + _neighNr );}
  //! get the id's of the most straight neighbours of the site
  const unsigned int& get_straight(const unsigned int& j, const short int& i) const{ return straight[j][i]; }
  //! get the length of the straight vector
  short int get_numStraight(const unsigned char& i) { return (short) straight[i].size(); }

  //! delete the intensityIn array 
  void delete_intensityIn(const short int numFreq){
    if(intensityIn){
      for(short int f=0; f<numFreq;f++){
        delete [] intensityIn[f];
      }
      delete [] intensityIn;
      intensityIn = NULL;
    }
  }
  //! create the Intensity array
  void create_intensityIn( const short int& numFreq, const unsigned int& directions ){ 
    intensityIn = new float*[numFreq];
    for(short int f=0;f<numFreq;f++){
      intensityIn[f] = new float[directions];
      for( unsigned int j=0; j<directions; j++ ){
        intensityIn[f][j] = 0.0;
      }
    }
  }  
  //! store
  void set_intensityIn(const short int& f, const unsigned int& j, const float& in){ intensityIn[f][j]=in; }
  //! get the id of neighbour of site
  const float& get_intensityIn(const short int& f, const unsigned int& i) const{ return intensityIn[f][i]; }

  //! delete the intensityIn array 
  void delete_intensityOut(const short int& numFreq){
    if(intensityOut){
      for(short int f=0; f<numFreq; f++){
        delete [] intensityOut[f];
      }
      delete [] intensityOut;
      intensityOut = NULL;
    }
  }
  //! create the IntensityOut array
  void create_intensityOut( const short int& numFreq, const unsigned int& directions ){ 
    intensityOut = new float*[numFreq];
    for(short int f=0; f<numFreq; f++){
      intensityOut[f] = new float[directions]; 
      for( unsigned int j=0; j<directions; j++ ){
        intensityOut[f][j] = 0.0;
      }
    }
  }  
  //! store
  void set_intensityOut(const short int& f, const unsigned int& j, const float& out){ intensityOut[f][j]=out; }
  //! get the id of neighbour of site
  const float& get_intensityOut(const short int& f, const unsigned int& i) const{ return intensityOut[f][i]; }

  //! delete the outgoing vector 
  void delete_outgoing(){
    if(outgoing){
      delete [] outgoing;
      outgoing = NULL;
    }
  }
  //! create the outgoing vector
  void create_outgoing( const unsigned int& directions ){ outgoing = new unsigned char[directions]; }  
  //! set outgoing directions for solid angles
  void set_outgoing( const unsigned int& entry, const unsigned char& value ){ outgoing[entry] = value; }
  //! get outgoing directions for solid angles
  unsigned char& get_outgoing( const unsigned int& entry ){ return outgoing[entry]; }

  //! add the intensity from neighbour in ingoing intensity array
  void addRadiationDiffIn(const short int& f, const unsigned int& id, const double& intensity);
  //! add the intensity from neighbour in outgoing intensity array
  void addRadiationDiffOut(const short int& f, const unsigned int& id, const double& intensity);


 private:

  unsigned long long int site_id;   //!< Position of the site in the local sites vector
  unsigned char numNeigh;           //!< Number of neighbours of the vertex
  float volume;                     //!< Volume of the Voronoi cell around vertex
  float neigh_dist;                 //!< Mean neighbour distance
  float n_HI;                       //!< Number density of neutral hydrogen
  float n_HII;                      //!< Number density of ionised hydrogen
  float internalEnergy;             //!< Internal energy of the gas
  float dinternalEnergydt;          //!< Rate of change of internal energy (adiabatic)
  float clumping;                   //!< Clumping factor of the gas
  bool ballistic;                   //!< Ballistic transport or direction conserving?
  bool source;                      //!< Source or not?

  float* flux;                      //!< Flux of the site in case its a source

  unsigned int* neighId;                     //!< vector holding id's of the neighbours of the vertex
  vector< vector< unsigned int > > straight; //!< Straightest neighbours
  float** intensityIn;                       //!< Array of incoming intensities
  float** intensityOut;                      //!< Array of outgoing intensities
  unsigned char* outgoing;                   //!< Array of numbers associating the angular direction (48 in total) with a vertex

 public:
  Site& operator=(const Vertex& p2); 
  Site& operator=(const Site_Update& p2); 

};

class cooling_curve{

  friend class SimpleX;

 public:
  cooling_curve(string fileName);
  ~cooling_curve();

  double get_value(double T);

 private:

  int elements;
  vector<double> Temperature;
  vector<double> Cooling;
  double min;
  double max;
  double delta;

};

class Send_Intensity{

 public:

  friend class Site;

  Send_Intensity();
  ~Send_Intensity(){};

  void set_process(const unsigned int& proc){ process=proc;}
  const unsigned int& get_process() const{ return process; }

  void set_neighId(const unsigned long long int& iden){ neighId=iden; }
  //! get the id of neighbour of site
  const unsigned long long int& get_neighId() const{ return neighId; }

  void set_id(const unsigned long long int& index){ id=index; }
  //! get the id of neighbour of site
  const unsigned long long int& get_id() const{ return id; }

  void set_freq_bin(const short int& index){ freq_bin=index; }
  //! get the id of neighbour of site
  const short int& get_freq_bin() const{ return freq_bin; }

  void set_intensityIn(const float& inten){ intensityIn=inten; }
  //! get the id of neighbour of site
  const float& get_intensityIn() const{ return intensityIn; }

  void set_intensityOut(const float& inten){ intensityOut=inten; }
  //! get the id of neighbour of site
  const float& get_intensityOut() const{ return intensityOut; }

 private:

  unsigned int process;
  unsigned long long int neighId;
  unsigned long long int id;
  short int freq_bin;
  float intensityIn;
  float intensityOut;

};


class Send_Site{

 public:

  Send_Site(){};
  ~Send_Site(){};

  void set_process(const unsigned int& proc){ process=proc;}
  const unsigned int& get_process() const{ return process; }

  void set_vertex_id(const unsigned long long int& index){ vertex_id=index; }
  const unsigned long long int& get_vertex_id() const{ return vertex_id; }

  void set_site_id(const unsigned long long int& index){ site_id=index; }
  const unsigned long long int& get_site_id() const{ return site_id; }

  //! assign ballistic 
  void set_ballistic(const unsigned int& ballistic_){ ballistic = ballistic_; }
  //! get the flux in the site
  const unsigned int& get_ballistic() const{ return ballistic; }

 private:

  unsigned int process;
  unsigned long long int vertex_id;
  unsigned long long int site_id;
  unsigned int ballistic;  

};

class Send_Neigh{

 public:

  friend class Site;

  Send_Neigh(){};
  ~Send_Neigh(){};

  void set_ballistic(const unsigned int& proc){ ballistic=proc;}
  const unsigned int& get_ballistic() const{ return ballistic; }

  void set_process(const unsigned int& proc){ process=proc;}
  const unsigned int& get_process() const{ return process; }

  void set_vertex_id(const unsigned long long int& index){ vertex_id=index; }
  const unsigned long long int& get_vertex_id() const{ return vertex_id; }

  void set_neighId(const unsigned long long int& index){ neighId=index; }
  const unsigned long long int& get_neighId() const{ return neighId; }

  void set_neighIdLoc(const unsigned int& index){ neighIdLoc=index; }
  const unsigned int& get_neighIdLoc() const{ return neighIdLoc; }

 private:

  unsigned int ballistic;
  unsigned int process;
  unsigned long long int vertex_id;
  unsigned long long int neighId;
  unsigned int neighIdLoc;


};

//! Class that holds information of updated sites
class Site_Update{

 public:

  friend class Vertex;
  friend class Site;

  //!< Constructor
  Site_Update();  

  //! set process
  void set_process(const unsigned int& proc){ process=proc;}
  //! get process
  const unsigned int& get_process() const{ return process; }


  //! set id of the site

  //! This is the position of the site in the local sites array
  void set_site_id( const unsigned long long int& index ){ site_id = index; }
  //! get local id of the vertex
  const unsigned long long int& get_site_id() const{ return site_id;  }

  //! set id of the vertex

  //! This is the global id of the vertex
  void set_vertex_id( const unsigned long long int& index ){ vertex_id = index; }
  //! get local id of the vertex
  const unsigned long long int& get_vertex_id() const{ return vertex_id;  }

  //! assign the number density of neutral hydrogen to site
  void set_n_HI(const float& nH){ n_HI = nH; }
  //! get the number density of neutral hydrogen of site
  const float& get_n_HI() const{ return n_HI; }

  //! assign the number density of ionised hydrogen to site
  void set_n_HII(const float& nH){ n_HII = nH; }
  //! get the number density of ionised hydrogen of site
  const float& get_n_HII() const{ return n_HII; }

  //! set the 
  void set_internalEnergy(const float& _u){ internalEnergy = _u; }
  //! get the 
  const float& get_internalEnergy() const{ return internalEnergy; }

  //! assign dudt to site
  void set_dinternalEnergydt(const float& _du){ dinternalEnergydt = _du; }
  //! get the dudt in the site
  const float& get_dinternalEnergydt() const{ return dinternalEnergydt; }

  //! assign flux to site
  void set_flux(const float& phi){ flux = phi; }
  //! get the flux in the site
  const float& get_flux() const{ return flux; }

  //! assign ballistic 
  void set_ballistic(const unsigned int& ballistic_){ ballistic = ballistic_; }
  //! get the flux in the site
  const unsigned int& get_ballistic() const{ return ballistic; }


 private:

  unsigned int process;
  unsigned long long int vertex_id; //!< Position of the site in the local sites vector
  unsigned long long int site_id;   
  float n_HI;
  float n_HII;
  float internalEnergy;             //!< internalEnergy in site
  float dinternalEnergydt;          //!< dudt in site
  float flux;                       //!< Flux of the site in case it's a source
  unsigned int ballistic;           //!< Ballistic transport or direction conserving?

 public:
  Site_Update& operator=(const Site& p2); 

};


class Site_Remove{

 public:

  void set_process(const unsigned int& proc){ process=proc;}
  const unsigned int& get_process() const{ return process; }

  //! This is the position of the site in the local sites array
  void set_site_id( const unsigned long long int& index ){ site_id = index; }
  //! get local id of the vertex
  const unsigned long long int& get_site_id() const{ return site_id;  }

  //! set id of the vertex

  //! This is the global id of the vertex
  void set_vertex_id( const unsigned long long int& index ){ vertex_id = index; }
  //! get local id of the vertex
  const unsigned long long int& get_vertex_id() const{ return vertex_id;  }

  //! set the update property
  void set_update_property(const float& prop){ update_property = prop; }
  //! get the update property
  const float& get_update_property() const{ return update_property; }


 private:

  unsigned int process;
  unsigned long long int vertex_id;
  unsigned long long int site_id;
  float update_property;

};


//!compare Vertex class by vertex id for sort routine
bool compare_vertex_id_vertex( const Vertex&, const Vertex& );

//!compare Site class by vertex id for sort routine
bool compare_vertex_id_site(const Site& x, const Site& y);

//!compare Send_Site class by process for sort routine
bool compare_process_send_site(const Send_Site& x, const Send_Site& y);

//!compare Send_Neigh class by vertex id for sort routine
bool compare_vertex_id_send_neigh(const Send_Neigh& x, const Send_Neigh& y);
//!compare Send_Neigh class by process for sort routine
bool compare_process_send_neigh(const Send_Neigh& x, const Send_Neigh& y);

//!compare Send_Intensity class by process for sort routine
bool compare_process_send_intensity(const Send_Intensity& x, const Send_Intensity& y);

//!compare Site_Remove class by process for sort routine
bool compare_process_site_remove(const Site_Remove& x, const Site_Remove& y);

//!compare Site_Update class by vertex id for sort routine
bool compare_vertex_id_site_update( const Site_Update& x, const Site_Update& y);

//!compare Site_Update class by process for sort routine
bool compare_process_site_update( const Site_Update& x, const Site_Update& y);


#endif
