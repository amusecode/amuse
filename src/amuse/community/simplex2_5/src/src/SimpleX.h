/*************************************************************************

description:
This file contains the class that can do the actual radiative transport, 
including routines for grid calculation, physics and the actual transport

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

#ifndef SIMPLEX_H
#define SIMPLEX_H

#include "rates.h"
#include "Common.h"
#include "Structs.h"
#include "mpi.h"
#include <gsl/gsl_rng.h>  //random number generator
#include "configfile.h"   //keyvalue input file

#ifdef HEAL_PIX
  #include "healpix_base.h" //healpix header
#endif


#include "tree_structures.h" //octree
#include "hilbert.h"         //hilbert curve

#include "h5w_serial.h"      //hdf5 header


#include <algorithm>

#if defined(__cplusplus)
extern "C"
{
#endif
#include <stdio.h>
#include <stdlib.h>
#include <qhull.h>
#include <mem.h>
#include <qset.h>
#include <geom.h>
#include <merge.h>
#include <poly.h>
#include <io.h>
#include <stat.h>
#if defined(__cplusplus)
}
#endif

using namespace std;

#define VERTEX_ITERATOR vector< Vertex >::iterator
#define SITE_ITERATOR vector< Site >::iterator
#define SIMPL_ITERATOR vector< Simpl >::iterator

//! Class in which all functions that work on the SimpleX grid are stored

//! This is the main class of the simulation!
class SimpleX{

  public:
    //! constructor
    SimpleX(const string& output_path = ".");
    //! destructor
    ~SimpleX();

  //================ Grid Initilisation Functions ==================================//

    //! Create triangulation

    //! Initialize the simulation by creating/reading the vertex list and compute the triangulation
    //! Result is a list of Sites and a list of Simplices
    void init_triangulation(char* inputName);


    //TEST
    void reinit_triangulation();
    //END TEST


    //! Read in parameter file
    void read_parameters( char* inputName );

    //! Create orientations and mappings from correct header file
    void set_direction_bins();

    //! Create a simple homogeneous point distribution 
    
    //! Points are placed using the gsl randoom number generator
    void poisson_square();

    //! Read in hdf5 file with site information
    void read_vertex_list();
  
    //! Create boundary around computational domain

    //! Width of the boundary is variable borderBox. 
    //! Number of points in the boundary is variable borderSites.
    void create_boundary();


    //!Create octree containing all vertices

    //!Size of the tree depends on the number of subboxes used, which
    //!is determined by the hilbert order m
    void create_vertex_tree();

    //! Decompose the domain

    //! Count the number of sites contained in one subbox and add those until 
    //! number of points is approximately equal per processor
    void decompose_domain();
    
    //! Assign the correct process to vertices
    void assign_process();

    //! Check if site is in unit domain
    bool inDomain(const float& x, const float& y, const float& z, const unsigned int& rank);

    //! Compute the triangulation

    //! Call to QHull to perform the triangulation of vertices in one subbox
    //! The result is a vector of simplices
    //! Also, a check is done whether the subbox and domain boundaries are
    //! sufficiently large. If not, subbox boundaries are extended. If the
    //! boundary around the unity domain is too small, the code exits,
    //! to avoid memory issues when this number would be increased
    void compute_triangulation();

    //! Create the sites array on which radiative transfer is performed

    //! From the list of simplices the vertices relevant for this proc are selected and
    //! put in the sites vector
    void create_sites();

    //! Give each site an id corresponding to place in local sites vector
    void assign_site_ids();

    //! Determine which sites use ballistic transport and which sites 
    //! should use direction conserving transport at the start of the
    //! simulation. 
    void initiate_ballistic_sites();

    //! Shuffle list of sites

    //! This is necessary since the sites list is created from the simplices returned by QHull,
    //! so they are ordered.
    void shuffle_sites();


  //================ Site Properties Functions ==================================//

    //! Compute properties of the sites from the triangulation

    //! Important properties like the neighbours, the volumes and 
    //! the most straight paths are computed
    void compute_site_properties();

    //! Compute the neighbours of vertex

    //! Generate the neighbour vectors for each vertex from the simplex vector
    //! Care should be taken so that no vertex is put twice in the neighbour list. 
    //! Also get rid of vertices in the sub box boundary that are not connected to a site
    //! in the computational domain
    void compute_neighbours();


    //! Compute the volume of the Voronoi cell

    //! Use the volumes of the simplices that the vertex is part of to compute
    //! the volume of the Voronoi cell around the vertex
    void compute_volumes();


    //! Calculate the most straight paths for every incoming direction

    //! Needed for fast radiative transport
    void calculate_straight();

    //! Identify the place of every ballistic site in neighbour array
    //! of neighbouring site
    void match_neighbours();
 
    //! Compute the solid angles into which the intensity is distributed
    void compute_solid_angles( const bool& rotate );


    //! Check for cyclic connections in the grid and remove them if possible
    void check_cyclic_connections();


    //! Calculate mean distance to neighbours
    void calculate_line_lengths();


    //! Create list without border simplicesc
    void remove_border_simplices();


    //================ Grid Update Routines ============================//

    //! Rotate the solid angles of the unit sphere 
    //! to avoid preferential directions
    void rotate_solid_angles();

    //! calculate normal distributed random numbers using the Box-Muller transformation
    void randGauss( float& y1, float& y2 );

    //! Store the intensities for use in later run, directions are kept.
    void store_intensities();

    //! Get the ballistic sites that have intensity to be stored
    const vector< unsigned long int > get_ballistic_sites_to_store();

    //! Store ballistic intensities in unit sphere tesselation
    void store_ballistic_intensities( const vector< unsigned long int >& sites_to_store );

    //! Return the intensities from the previous run
    void return_intensities();

    //! Return the intensities to ballistic sites
    void return_ballistic_intensities();

    //! Check which sites are ballistic and which not
    //! according to user specified switch

    //! Only used is case of combined transport
    void check_ballistic_sites();

    //! Create new sites array without the removed sites
    void store_site_properties();

    //! Create updated vertex list from site_properties vector
    void create_new_vertex_list();


    //================ Send Routines ===================================//
   

    //! Send vertices to processors
    void send_vertices();

    //! Send the domain decomposition to all procs
    void send_dom_dec();

    //! Fill the list of indices of sites that contain information to be send 
    //! to other procs
    void fill_send_list();

    //! Give every vertex its correct process (needed for vertices in boundary between procs)
    void send_site_properties();

    //! Give every ballistic site its correct neighbour information 
    
    //! This is only needed for ballistic sites in boundary between procs, and is necessary because
    //! the order of the neighbours in the neighbour vector is not the same among procs.
    void send_neighbour_properties();

    //! Send teh updated ballistic sites to all procs
    void send_site_ballistics( const vector< unsigned long int >& sites_to_send );

     //! Send the intensities among procs
    void send_intensities();

    //! Send site physics to all procs

    //! In case sites change from proc send physics
    void send_site_physics();

    //! In case sites change from proc send intensities
    void send_site_intensities();

    //! Send stored intensities among procs
    void send_stored_intensities();

    //! Send the updated vertex list to all procs
    void send_new_vertex_list();


    //================ Physics Functions ==================================//


    //! Compute all the physics needed to do radiative transport
    void compute_physics( const unsigned int& run );

    //! Initialize physical parameters
    void initialise_physics_params();

    //! Read metal line cooling data and initialize vectors
    void read_metals();

    //! set homogeneous number density. 

    //! Use only in combination with Poisson_Square()
    void set_homogeneous_number_density( const float& nH );


    //! Give the masses and fluxes that were read in to the sites

    //! Use only in combination with read_vertex_list()
    void assign_read_properties();

    //! Give the updates sites their correct physical quantities from stored values
    void return_physics();

    //! Calculate the bounds of the frequency bins
    vector<double> calc_freq_bounds( const double& lowerBound, const double& upperBound );

    //! calculate effective cross section and spectrum of black body source
    void black_body_source( const double& tempSource );

    //! Compute the total number of atoms on this proc, used for output
    void output_number_of_atoms();

    //! Return mean optical depth of total grid
    double output_optical_depth();

    //! output relevant physical parameters to screen 
    void parameter_output( const unsigned int& run );
    
    //! Add more neighbours to vertices that have a flux, by adding the neighbours
    //! of neighbours.
    void make_source_isotropic();

    //! calculate recombination
    double recombine();

    //! Rate of energy loss in cell
    double cooling_rate( Site& site );

    //! Rate of energy gain in cell
    double heating_rate( const vector<double>& N_ion, const double& t_end );

    //!Convert internal energy to temperature
    double u_to_T( const double& u, const double& mu );

    //!Convert temperature to internal energy
    double T_to_u( const double& T, const double& mu );

    //! update the temperature by including heating and cooling processes
    double update_temperature( Site& site, const vector<double>& N_ion, const double& t_end );

 //================= Radiative Transfer Functions ======================//
    
    //! Solve the rate equation to determine ionisations and recombinations
    vector<double> solve_rate_equation( Site& site );

    //! Send source photons
    void source_transport( Site& site );

    //! Transport photons diffusely
    void diffuse_transport( Site& site, vector<double>& N_out_total );

    //! Redistribute intensity: actual radiative transport
    void non_diffuse_transport( Site& site, vector<double>& N_out_total );

    //! Call to do radiative transfer
    void radiation_transport( const unsigned int& run );

  //================== Output Functions ================================//


    //! Call the output routines
    void generate_output( const unsigned int& run );    

    //! Write hdf5 output file
    void write_hdf5_output(char* name, const unsigned int& run);
    
    //! calculate position of the I-Front for source in centre
    double calc_IFront( const unsigned int& run );

    //! calculate escape fraction
    double calc_escape_fraction( const unsigned int& run, const unsigned int& times );


    //================ Generic Functions ==============================//
   
    //! clear all temporary structures for the next run
    void clear_temporary();
    
    //! Function to give diagnostic output about total number of photons in simulation
    double count_photons();

    //! return number of runs
    const unsigned int& get_numRuns() const{ return numRuns; }
  
    //! Return number of outputs
    const unsigned int& get_numOutputs() const{ return numOutputs; }

    // ================ Public lists and variables ===================//

    //! List of vertices that will be used to perform triangulation
    vector< Vertex > vertices;

    //! List of sites for doing radiative transfer
    vector< Site > sites;

    //! Vector holding the domain decomposition
    vector< unsigned int > dom_dec;

    //! List of sites from the previous run, whose properties are still needed
    vector< Site_Update > site_properties;

    //! List of entries in site_intensities array of every site
    vector<unsigned long int> intens_ids;
    //! List of site intensities from previous run, that are still needed
    vector<float> site_intensities;

    //! List of indices of sites needed for transport along procs
    vector< unsigned long int > send_list;

    //! List of site indices that might be removed
    vector< unsigned long int > sites_marked_for_removal;
    //! List of sites and properties to be updated
    vector< Site_Remove > total_sites_marked_for_removal;

    //! List of simplices that results from the QHull call
    vector< Simpl > simplices;
    //! Temporary list of neutral number densities that are read in
    vector< float > temp_n_HI_list;
    //! Temporary list of ionised number densities that are read in
    vector< float > temp_n_HII_list;
    //! Temporary list of fluxes that are read in
    vector< float > temp_flux_list;
    //! Temporary list of temperatures that are read in
    vector< float > temp_T_list;
    //! Temporary list of clumping factors that are read in
    vector< float > temp_clumping_list;
    //! Cross section of hydrogen
    vector<double> cross_H;
    //! The photon excess energy used to heat the gas
    vector<double> photon_excess_energy;  
    //! The elemental abundances of metals (and Helium)
    vector<double> abundances;
    //! Vector of cooling curves
    vector< cooling_curve > curves;

    //! Variable used by random number generator
    gsl_rng * ran;
    //! Output to logfile
    ofstream simpleXlog; 
    
  protected:

    unsigned long int numSites;            //!< Total number of sites in the simulation
    unsigned long int origNumSites;        //!< Original number of sites
    unsigned int borderSites;      //!< Number of sites in the boundary around computational domain
    unsigned int hilbert_order;    //!< Hilbert order that determines the number of subboxes
    unsigned int numSweeps;        //!< Number of grid sweeps during simulation
    unsigned int numRuns;          //!< Number of runs of the simulation
    unsigned int numOutputs;       //!< Number of outputs of the simulation
    unsigned int numPixels;        //!< Number of pixels in HEALPix calculation
    unsigned int numStraight;      //!< Number of straight directions
    int randomSeed;                //!< Seed for random number generator
    short dimension;               //!< Dimension of the simulation (default = 3)
    string inputFileName;          //!< File which holds densities and fluxes to be read in
    bool blackBody;                //!< Use blackbody for source spectrum?
    bool recombination;            //!< Include recombination or not
    bool coll_ion;                 //!< Include collisional ionisations?
    bool heat_cool;                //!< Include heating and cooling?
    bool metal_cooling;            //!< Include metal line cooling?
    short int numFreq;             //!< Number of frequencies
    short int freq_spacing;        //!< Spacing of the frequency bins
    short RTmethod;                //!< Method for RT
    bool diffuseTransport;       //!< Set this to 1 to do only ballistic transport
    bool ballisticTransport;       //!< Set this to 1 to do only ballistic transport
    bool dirConsTransport;         //!< Set this to 1 to do only direction conserving transport
    bool combinedTransport;        //!< Set this to 1 to do combined transport
    bool photon_conservation;      //!< Use temporal photon conservation scheme or not?
    double subcycle_frac;          //!< Fraction of the characteristic time scale at which subcycling is done
    double sizeBox;                //!< Size of the box in parsec
    float borderBox;               //!< Size of boundary around computational domain
    float padding_subbox;          //!< Subbox boundary in which extra points are triangulated
    double simTime;                //!< Total simulation time in Myr
    double time_conversion;        //!< Time conversion for output
    double sourceTeff;             //!< Effective temperature of the sources in K
    double gasIonTemp;             //!< Temperature of the ionised gas in K
    double switchTau;              //!< Optical depth below which ballistic transport is no longer trusted
    double totalVolume;            //!< Total volume of all Voronoi cells added up
    double totalAtoms;             //!< Total number of atoms in the simulation domain
    int localMaxRes;               //!< Maximum resolution on the local processor
    unsigned int COMM_RANK;        //!< Rank of this proc
    unsigned int COMM_SIZE;        //!< Total number of procs
    unsigned int chunk_size;       //!< Maximum size of chunks written in one piece to hdf5 file
    unsigned int max_msg_to_send;  //!< Maximum number of messages to send in MPI call
    unsigned long int vertex_id_max;    //!< Maximum vertex index
    short int orientation_index;      //!< Orientation of the unit sphere tesselation
    short int orientation_index_old;  //!< Orientation of the unit sphere tesselation in previous run
    float euler_phi;               //!< Euler angle phi for random rotation of coordinate system
    float euler_theta;             //!< Euler angle theta for random rotation of coordinate system
    float euler_psi;               //!< Euler angle psi for random rotation of coordinate system
    bool heal_pix;                 //!< Use heal_pix for associating directions and neighbours or not?
    int num_ref_pix_HP;            //!< Number of pixels to use for associating neighbours in compute_solid_angles()
    float straightAngle;           //!< Maximum angle allowed between straight direction and Delaunay direction.
    bool straight_from_tess;       //!< In DCT, calculate straightest neighbours wrt Delaunay edge or direction bin?
    double UNIT_L;                 //!< Unit length in box
    double UNIT_M;                 //!< Unit mass in box, converts number density into number of atoms
    double UNIT_I;                 //!< Unit number of ionising photons in box
    double UNIT_T;                 //!< Time of 1 sweep in seconds
    double UNIT_T_MYR;             //!< Time of 1 sweep in Myr
    double UNIT_D;                 //!< Unit number density
    double UNIT_V;                 //!< Unit volume
    float straight_correction_factor; //!< Correction factor to correct for the fact that no straight lines exist on the grid
    vertex_tree vert_tree;         //!< Octree that will contain the vertices
    bool give_IFront;                   //!< Calculate IFront position for a source in centre?
    ofstream IFront_output;        //!< Output of IFront data to file
    float*** orient;               //!< Array to hold the direction bins
    unsigned int*** maps;          //!< Array to hold the mappings between direction bins
    unsigned int number_of_directions; //!< Number of discretizations of the unit sphere for DCT
    unsigned int number_of_orientations;//!< Number of random rotations of the DCT discretization
};

#endif
