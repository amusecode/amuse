#include "Main.h"
#include "configfile.h"

#include "amuse_interface.h"

#ifndef M_PI
#define M_PI           3.14159265358979323846
#endif

using namespace std;

//initialise parameters that are usually read in from file
void AMUSE_SimpleX::read_parameters(){

//----------- Change these when necessary  ----------------//

  //random seed
  randomSeed = 25061977;

  //hilbert order, determines number of subboxes. 
  //should be > 0 for parallel runs!
  hilbert_order = 1;

  //time step
  UNIT_T_MYR = 0.05;
  //size of simulation domain in pc
  sizeBox = 13200;

  //use blackbody for source? Otherwise monochromatic
  blackBody = 0;
  //effective temperature of source
  sourceTeff = 1.e5;
  
  //units of source
  UNIT_I = 1.e48;
  
  //number of frequencies
  numFreq = 1;
  
  //include collisional ionisations?
  coll_ion = 0;

  //include heating and cooling?
  heat_cool = 0;

  //include metal line cooling?
  metal_cooling = 0;
  
  //include recombination radiation?
  rec_rad = 0;
  
//----------- Don't change these unless you know what you're doing -------------------------//

 //dimension should always be 3
  dimension = 3;
  
  //number of points in border, set this as you like
  borderSites = 25000;
  //region in which border points are placed
  borderBox = 0.1;
  //size of border around subbox
  padding_subbox = 0.25;
  
  //temperature of the ionised gas
  gasIonTemp = 1.e4;
  
  //temperature of neutral gas
  gasNeutralTemp = 100.0;

  //include recombinations?
  recombination = 1;

  subcycle_frac=0.05;

  //type of transport
  ballisticTransport = 0;
  dirConsTransport = 0;
  combinedTransport = 1;  

  //Chunk size for hdf5 writing
  chunk_size = 100000;
  //Maximum number of messages to send in MPI routines
  max_msg_to_send = 100000;

  //number of reference pixels in the case HEAL_PIX is used
  num_ref_pix_HP = 4;

  // Maximal angle in degrees between straightforward direction and Delaunay lines (second and third)
  straightAngle = 90;
  //convert to radians
  straightAngle *= M_PI/180.0;

  //Switch between dct and ballistic transport
  switchTau = 0.5;
  
  //Use temporal photon conservation?
  photon_conservation = 1;

  //Number of directions
  number_of_directions = 42;
  //Number of orientations, fixed in header files
  number_of_orientations = 100;
  
  //calculate straight from tesselation instead of direction bins
  //is faster in postprocessing, but less precise
  straight_from_tess = 1;

  freq_spacing = ENERGY_WEIGHTS;
//----------------------------------------------------------------------//


}

//add one vertex to the list of vertices that needs to be triangulated
int AMUSE_SimpleX::add_vertex(long *id, double x,double y,double z,double rho,
                                           double flux,double xion,double uInt, double metallicity){

  Vertex tempVert;
  if( COMM_RANK ==0){
    tempVert.set_x(  x );
    tempVert.set_y(  y );
    tempVert.set_z(  z );
    unsigned long long int tempid=*id;
    tempVert.set_vertex_id( tempid);
    tempVert.set_border( 0 );
    tempVert.set_process( COMM_RANK );

    vertices.push_back( tempVert );
    numSites=vertices.size();
  }
  
  double n_HI = rho*(1.0-xion);
  double n_HII = rho*xion;
  temp_n_HI_list.push_back( (float) n_HI);
  temp_n_HII_list.push_back( (float) n_HII );
  temp_flux_list.push_back(flux);
  temp_u_list.push_back(uInt);
  temp_dudt_list.push_back(0.0);
  temp_metallicity_list.push_back(metallicity);

  return 0;

}

//add a site to the sites list
int AMUSE_SimpleX::add_site(long *id, double x,double y,double z,double rho,
                                           double flux,double xion, double uInt, double metallicity){

  Site tempSite;

  double n_HI = rho*(1.0-xion);
  double n_HII = rho*xion;

  if(inDomain(x,y,z,COMM_RANK)){
    tempSite.set_x( x );
    tempSite.set_y( y );
    tempSite.set_z( z );
    tempSite.set_vertex_id( (unsigned long long int) *id );
    tempSite.set_border( 0 );
    tempSite.set_process( COMM_RANK );
    tempSite.set_n_HI( (float) n_HI );
    tempSite.set_n_HII( (float) n_HII );
    if(flux > 0.0){
      tempSite.set_source(1);
      tempSite.create_flux(numFreq);
      tempSite.set_flux( 0, flux );
    }else{
      tempSite.set_source(0);
    }
    tempSite.set_neigh_dist(1.);
    tempSite.set_internalEnergy(uInt);
    tempSite.set_dinternalEnergydt(0.0);
    tempSite.set_metallicity(metallicity);
    
    sites.push_back( tempSite );
    numSites++;
    
    // cerr << " Add site " << tempSite.get_vertex_id() << " x: " << tempSite.get_x() << " y: " << tempSite.get_y() << " z: " << tempSite.get_z() 
    //      << " n_HI: " << tempSite.get_n_HI() << " n_HII: " << tempSite.get_n_HII() << endl;
        
    return 0;
  } else
  {
    numSites++;
    return -1;
  }
}


//remove a site from the list
int AMUSE_SimpleX::remove_site(long id) {
   SITE_ITERATOR p;
   Site tmp;
   
   tmp.set_vertex_id((unsigned long long) id);

   p=lower_bound(sites.begin(),sites.end(), tmp, compare_vertex_id_site);
   if(p->get_vertex_id() == (unsigned long long int) id){
     if (p->get_process() == COMM_RANK){
       p->set_neigh_dist(-1.);
       return 1;
     }
   }  
   return 0;
}

//convert the physical units to code units
void AMUSE_SimpleX::convert_units(){

  //convert parsec to cm
  UNIT_L = sizeBox * 3.08568025e18; 
  UNIT_V = pow(UNIT_L, 3);
  
  
  //calculate mean number density
  //double n_H = 0.0;
  //for(unsigned int i=0; i<temp_n_HI_list.size(); i++){
    //n_H += temp_n_HI_list[i] + temp_n_HII_list[i];
  //}
  //n_H/=temp_n_HI_list.size();
  
  //set UNIT_D and UNIT_M
  UNIT_D = 1.0;//n_H;
  UNIT_M = UNIT_D * UNIT_V;

  //set the neutral and ionised hydrogen fractions
  //for(unsigned int i=0; i<temp_n_HI_list.size(); i++){
    //temp_n_HI_list[i] /= UNIT_D;
    //temp_n_HII_list[i] /= UNIT_D;  
  //}
  
}

//initialise the code:
// create boundary around domain
// do domain decomposition
// compute triangulation
// create the sites for RT
// compute the physical properties of the sites
int AMUSE_SimpleX::initialize(double current_time) {
  
  // if(COMM_RANK == 0){
  //   cerr << "AMUSE_SimpleX: initialising...";
  // }

  if(rec_rad){
    numFreq++;
  }
  
  // initialise random number generator with random seed
  gsl_rng_set ( ran , randomSeed );

  //set the correct directions headers
  set_direction_bins();
  
  //first set the relevant units
  convert_units();
    
  //only the master proc creates the point distribution (in case of automatic filling) and
  //the boundary around the domain (always)
  
  if(COMM_RANK == 0){
    numSites = vertices.size();  
    origNumSites = vertices.size();
    vertex_id_max = vertices.size() + borderSites;
  }
  
  MPI::COMM_WORLD.Bcast( &vertex_id_max, 1, MPI::UNSIGNED_LONG_LONG, 0 );
  MPI::COMM_WORLD.Bcast( &origNumSites, 1, MPI::UNSIGNED_LONG_LONG, 0 );
  MPI::COMM_WORLD.Bcast( &numSites, 1, MPI::UNSIGNED_LONG_LONG, 0 );
    
  total_time=current_time;
  
  if(COMM_RANK == 0){

    //set boundary around unity domain
    create_boundary();
    
    //create octree on master proc
    create_vertex_tree();

    //decompose the domain
    decompose_domain();
    
    //assign process to vertices
    assign_process();

  }

  //send the domain decomposition to other procs
  send_dom_dec();
  //send the vertex positions round to all procs, so every proc (temporary!) has its own copy
  send_vertices();

  //now that all procs have a list of vertices, create octree
  if( COMM_RANK != 0 ){
    create_vertex_tree();
  }

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
  MPI::COMM_WORLD.Barrier();
  // from main:
  compute_site_properties();
  compute_physics( 0 );
  remove_border_simplices();
  syncflag=0;
 
  // if(COMM_RANK == 0){
  //   cerr << " Done" << endl;
  // }
 
  return 0;
  
}

//set up a simulation
int AMUSE_SimpleX::setup_simplex(){

 read_parameters();
 vertices.clear();

 return 0;
}

//evolve the radiative transfer over time t_target
int AMUSE_SimpleX::evolve(double t_target, int sync) {
  

    
  double dt= t_target*secondsPerMyr - total_time;
  
  //cout << "time: " << t_target << " " << dt << "\n"; 
    
  if(syncflag==1){
    reinitialize();
  }
  
  // if(COMM_RANK == 0){
  //   cerr << "AMUSE_SimpleX: performing radiation transport...";
  // }
    
  if(COMM_RANK == 0){
    //cerr << "start sweeping" << endl;
    simpleXlog << endl << "  start sweeping till " << t_target << endl;  
  }

  numSweeps=1;
  while(total_time < t_target*secondsPerMyr - UNIT_T/2)
  {
  radiation_transport(1);
  total_time+=UNIT_T;
  }

  if(COMM_RANK == 0){
    simpleXlog << "   did " << numSweeps << ") sweeps"<< endl;
  }
  
  if(sync == 1) {
    syncflag=1;     
  }
  
  // if(COMM_RANK == 0){
  //   cerr << " Done" << endl;
  // }
  
  return 0;

}


//store everything in the simulation and recompute triangulation and physics
int AMUSE_SimpleX::reinitialize(){

  // if(COMM_RANK == 0){
  //   cerr << "AMUSE_SimpleX: recomputing triangulation...";
  // }
  
  // for( SITE_ITERATOR it=sites.begin(); it!=sites.end(); it++ ){
  //   if(it->get_process() == COMM_RANK && !it->get_border() ){
  //     cerr << " Site: " << it->get_vertex_id() << " nHI: " << it->get_n_HI() << " nHII: " << it->get_n_HII() << endl; 
  //   }
  // }  
  // 
   //make sure that the vectors that will be filled are empty 
    site_intensities.clear();
    intens_ids.clear();
    
    //make sure that photons have left all ghost vertices 
    //before calculating local intensities
    send_intensities();

    //store the intensities in big array
    store_intensities();

    //in this case, also store the ballistic intensities
    vector< unsigned long int > sites_to_store = get_ballistic_sites_to_store();

    store_ballistic_intensities( sites_to_store );
    sites_to_store.clear();

    //remember the orientation index with which the intensities were stored
    orientation_index_old = orientation_index;

    //store the relevant properties of the sites to be used 
    //in the coming run
    store_site_properties();

    //create a list of vertices to be triangulated
    create_new_vertex_list();

    //clear temporary structures in the sites 
    //and completely clear the sites vector
    clear_temporary();
    sites.clear();
    vector< Site >().swap(sites);

    simplices.clear();
    vector< Simpl >().swap(simplices);

    //send the list to master proc
    send_new_vertex_list();
 
    // if( COMM_RANK == 0 ){
    //   cerr << " (" << COMM_RANK << ") Computing triangulation" << endl;
    // }
     
    if(COMM_RANK == 0){

      //set boundary around unity domain
      create_boundary();
      
      // Sort the vertices
      sort( vertices.begin(), vertices.end(), compare_vertex_id_vertex );

      //create octree
      create_vertex_tree();

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

    //now that all procs have a list of vertices, create octree
    if( COMM_RANK != 0 ){
      // Sort the vertices
      sort( vertices.begin(), vertices.end(), compare_vertex_id_vertex );

      // Create the tree
      create_vertex_tree();
    }

    //check if vertices have moved to different process and send the
    //info accordingly
    if(COMM_SIZE > 1){
      send_site_physics();  
      send_site_intensities();
    }
            
    //compute the triangulation
    compute_triangulation();

    //create the sites vector from the vertex list
    create_sites();

    //assign the correct site ids to the sites
    assign_site_ids();  

    //return the physical properties to the sites
    return_physics();
            
    compute_site_properties();
    
    compute_physics( 1 );
    
    remove_border_simplices();
      
    // if(COMM_RANK == 0){
    //   cerr << " Done" << endl;
    // }
    
    // for( SITE_ITERATOR it=sites.begin(); it!=sites.end(); it++ ){
    //   if(it->get_process() == COMM_RANK && !it->get_border() ){
    //     cerr << " Site: " << it->get_vertex_id() << " nHI: " << it->get_n_HI() << " nHII: " << it->get_n_HII() << endl; 
    //   }
    // }
    
  syncflag=0;
  return 0;
}

//return properties of specified site
int AMUSE_SimpleX::get_site(int id, double *x,double *y,double *z,double *rho,
                                              double *flux,double *xion, double *uInt, double *metallicity){
   SITE_ITERATOR p;
   Site tmp;
   
   tmp.set_vertex_id((unsigned long long) id);
   p=lower_bound(sites.begin(),sites.end(), tmp, compare_vertex_id_site);
   if(p->get_vertex_id() == (unsigned long long int)id){
     if (p->get_process() == COMM_RANK){
       *x=p->get_x();
       *y=p->get_y();
       *z=p->get_z();
       double nH = (double) p->get_n_HI() + (double) p->get_n_HII();
       *rho = nH;
       double totalFlux = 0.0;
       if( p->get_source() ){
         for(short int f=0; f<numFreq;f++){
           totalFlux += p->get_flux(f);
         }
       }
       *flux = totalFlux;
       *xion = (double) p->get_n_HII()/nH;
       *uInt = p->get_internalEnergy();
       *metallicity = p->get_metallicity();
       
       return 1;
     }
   }
   return 0;
}

int AMUSE_SimpleX::get_position(int id, double *x,double *y,double *z){
    SITE_ITERATOR p;
    Site tmp;
    
    tmp.set_vertex_id((unsigned long long) id);
    p = lower_bound(sites.begin(), sites.end(), tmp, compare_vertex_id_site);
    if(p->get_vertex_id() == (unsigned long long int)id){
        if (p->get_process() == COMM_RANK){
            *x=p->get_x();
            *y=p->get_y();
            *z=p->get_z();
            return 1;
        }
    }
    return 0;
}

int AMUSE_SimpleX::get_density(int id, double *rho){
    SITE_ITERATOR p;
    Site tmp;
    
    tmp.set_vertex_id((unsigned long long) id);
    p=lower_bound(sites.begin(), sites.end(), tmp, compare_vertex_id_site);
    if(p->get_vertex_id() == (unsigned long long int)id){
        if (p->get_process() == COMM_RANK){
            *rho = (double) p->get_n_HI() + (double) p->get_n_HII();
            return 1;
        }
    }
    return 0;
}

int AMUSE_SimpleX::get_flux(int id, double *flux){
    SITE_ITERATOR p;
    Site tmp;
    
    tmp.set_vertex_id((unsigned long long) id);
    p=lower_bound(sites.begin(), sites.end(), tmp, compare_vertex_id_site);
    if(p->get_vertex_id() == (unsigned long long int)id){
        if (p->get_process() == COMM_RANK){
          double totalFlux = 0.0;
          if( p->get_source() ){
            for(short int f=0; f<numFreq;f++){
              totalFlux += p->get_flux(f);
            }
          }
          *flux = totalFlux;
            return 1;
        }
    }
    return 0;
}

int AMUSE_SimpleX::get_mean_intensity(int id, double *mean_intensity){
  
  SITE_ITERATOR p;
  Site tmp;
    
  tmp.set_vertex_id((unsigned long long) id);
  p=lower_bound(sites.begin(), sites.end(), tmp, compare_vertex_id_site);
  if(p->get_vertex_id() == (unsigned long long int)id){
      if (p->get_process() == COMM_RANK){
        
        double meanIntensity = 0.0;
        //in case of ballistic transport, intensity has size of number of neighbours;
        //in case of direction conserving transport, intensity has 
        //the size of the tesselation of the unit sphere
        numPixels = ( p->get_ballistic() ) ? p->get_numNeigh() : number_of_directions;
        
        for(unsigned int i=0; i<numPixels; i++){
          short int start = (rec_rad) ? 1 : 0;
          for(short int f=start;f<numFreq;f++){
            meanIntensity += p->get_intensityOut(f,i) + p->get_intensityIn(f,i);
          }
        }
        *mean_intensity = meanIntensity;
        return 1;
      }
  }
  return 0;
   
}

int AMUSE_SimpleX::get_diffuse_intensity(int id, double *diffuse_intensity){
  
  SITE_ITERATOR p;
  Site tmp;
    
  tmp.set_vertex_id((unsigned long long) id);
  p=lower_bound(sites.begin(), sites.end(), tmp, compare_vertex_id_site);
  if(p->get_vertex_id() == (unsigned long long int)id){
      if (p->get_process() == COMM_RANK){
        
        double diffuseIntensity = 0.0;
        numPixels = ( p->get_ballistic() ) ? p->get_numNeigh() : number_of_directions;
        
        for(unsigned int i=0; i<numPixels; i++){
          short int f = 0;
          diffuseIntensity += p->get_intensityOut(f,i) + p->get_intensityIn(f,i);
        }
        *diffuse_intensity = diffuseIntensity;
        return 1;
      }
  }
  return 0;
   
}

int AMUSE_SimpleX::get_ionisation(int id, double *xion){
    SITE_ITERATOR p;
    Site tmp;
    
    tmp.set_vertex_id((unsigned long long) id);
    p=lower_bound(sites.begin(), sites.end(), tmp, compare_vertex_id_site);
    if(p->get_vertex_id() == (unsigned long long int)id){
        if (p->get_process() == COMM_RANK){
            *xion = (double) p->get_n_HII() / ((double) p->get_n_HI() + (double) p->get_n_HII());
            return 1;
        }
    }
    return 0;
}

int AMUSE_SimpleX::get_metallicity(int id, double *metallicity){
    SITE_ITERATOR p;
    Site tmp;
    
    tmp.set_vertex_id((unsigned long long) id);
    p=lower_bound(sites.begin(), sites.end(), tmp, compare_vertex_id_site);
    if(p->get_vertex_id() == (unsigned long long int)id){
        if (p->get_process() == COMM_RANK){
            *metallicity = (double) p->get_metallicity();
            return 1;
        }
    }
    return 0;
}

int AMUSE_SimpleX::get_internalEnergy(int id, double *uInt){

  SITE_ITERATOR p;
  Site tmp;

  tmp.set_vertex_id((unsigned long long) id);
  p=lower_bound(sites.begin(), sites.end(), tmp, compare_vertex_id_site);
  if(p->get_vertex_id() == (unsigned long long int)id){
    if (p->get_process() == COMM_RANK){
      *uInt = p->get_internalEnergy();
      return 1;

    }
  }
  return 0;

  
}

int AMUSE_SimpleX::get_dinternalEnergydt(int id, double *uInt){

  SITE_ITERATOR p;
  Site tmp;

  tmp.set_vertex_id((unsigned long long) id);
  p=lower_bound(sites.begin(), sites.end(), tmp, compare_vertex_id_site);
  if(p->get_vertex_id() == (unsigned long long int)id){
    if (p->get_process() == COMM_RANK){
      *uInt = p->get_dinternalEnergydt();
      return 1;

    }
  }
  return 0;

  
}

//set properties of specified site
int AMUSE_SimpleX::set_site(int id, double x, double y, double z, double rho,
                                              double flux, double xion, double uInt, double metallicity){
    SITE_ITERATOR p;
    Site tmp;
    
    tmp.set_vertex_id((unsigned long long) id);
    p = lower_bound(sites.begin(), sites.end(), tmp, compare_vertex_id_site);
    if(p->get_vertex_id() == (unsigned long long int)id){
        p->set_x(x);
        p->set_y(y);
        p->set_z(z);
        //set number of ionising photons
        if(flux > 0.0){
          if(p->get_source()==0) p->create_flux(numFreq);
          p->set_source(1);
          p->set_flux( 0, flux );
        }else{
          if(p->get_source()==1) p->delete_flux();
          p->set_source(0);
        }
        p->set_n_HI((1-xion)*rho);
        p->set_n_HII(xion*rho);
	      p->set_internalEnergy( uInt );
        p->set_metallicity( metallicity );
        
        return 1;
    }
    return 0;
}

int AMUSE_SimpleX::set_position(int id, double x, double y, double z){
    SITE_ITERATOR p;
    Site tmp;
    
    tmp.set_vertex_id((unsigned long long) id);
    p = lower_bound(sites.begin(), sites.end(), tmp, compare_vertex_id_site);
    if(p->get_vertex_id() == (unsigned long long int)id){
        p->set_x(x);
        p->set_y(y);
        p->set_z(z);
        return 1;
    }
    return 0;
}

int AMUSE_SimpleX::set_density(int id, double rho){
    SITE_ITERATOR p;
    Site tmp;
    
    tmp.set_vertex_id((unsigned long long) id);
    p=lower_bound(sites.begin(), sites.end(), tmp, compare_vertex_id_site);
    if(p->get_vertex_id() == (unsigned long long int)id){
        double old_rho = (double) p->get_n_HI() + (double) p->get_n_HII();
        p->set_n_HI((double) p->get_n_HI() * rho/old_rho);
        p->set_n_HII((double) p->get_n_HII() * rho/old_rho);
        return 1;
    }
    return 0;
}

int AMUSE_SimpleX::set_flux(int id, double flux){
    SITE_ITERATOR p;
    Site tmp;
    
    tmp.set_vertex_id((unsigned long long) id);
    p=lower_bound(sites.begin(), sites.end(), tmp, compare_vertex_id_site);
    if(p->get_vertex_id() == (unsigned long long int)id){
        if(flux > 0.0){
          if(p->get_source()==0) p->create_flux(numFreq);
          p->set_source(1);
          p->set_flux( 0, flux );
        }else{
          if(p->get_source()==1) p->delete_flux();
          p->set_source(0);
        }
        return 1;
    }
    return 0;
}

int AMUSE_SimpleX::set_ionisation(int id, double xion){
    SITE_ITERATOR p;
    Site tmp;
    
    tmp.set_vertex_id((unsigned long long) id);
    p=lower_bound(sites.begin(), sites.end(), tmp, compare_vertex_id_site);
    if(p->get_vertex_id() == (unsigned long long int)id){
        double rho = (double) p->get_n_HI() + (double) p->get_n_HII();
        p->set_n_HI((1-xion)*rho);
        p->set_n_HII(xion*rho);
        return 1;
    }
    return 0;
}

int AMUSE_SimpleX::set_metallicity(int id, double metallicity){
    SITE_ITERATOR p;
    Site tmp;
    
    tmp.set_vertex_id((unsigned long long) id);
    p=lower_bound(sites.begin(), sites.end(), tmp, compare_vertex_id_site);
    if(p->get_vertex_id() == (unsigned long long int)id){
      p->set_metallicity( metallicity );
      return 1;
    }
    return 0;
}

int AMUSE_SimpleX::set_internalEnergy(int id, double uInt){
    SITE_ITERATOR p;
    Site tmp;
    
    tmp.set_vertex_id((unsigned long long) id);
    p=lower_bound(sites.begin(), sites.end(), tmp, compare_vertex_id_site);
    if(p->get_vertex_id() == (unsigned long long int)id){
      p->set_internalEnergy( uInt );
      return 1;
    }
    return 0;
}

int AMUSE_SimpleX::set_dinternalEnergydt(int id, double uInt){
    SITE_ITERATOR p;
    Site tmp;
    
    tmp.set_vertex_id((unsigned long long) id);
    p=lower_bound(sites.begin(), sites.end(), tmp, compare_vertex_id_site);
    if(p->get_vertex_id() == (unsigned long long int)id){
      p->set_dinternalEnergydt( uInt );
      return 1;
    }
    return 0;
}


//relevant global parameters for run
AMUSE_SimpleX *SimpleXGrid;
int COMM_RANK;
int lastid = 0;
string global_output_path = ".";
string global_data_path = ".";

int initialize_code() {
  COMM_RANK = MPI::COMM_WORLD.Get_rank();
  SimpleXGrid=new AMUSE_SimpleX(global_output_path, global_data_path);
  return (*SimpleXGrid).setup_simplex();  
}

int set_output_directory(char *output_path){
    global_output_path = output_path;
    return 0;
}

int set_data_directory( char *data_path ){
  global_data_path = data_path;
  return 0;
}

char buf[1024];
int get_data_directory(char **data_path ){
  strncpy(buf,global_data_path.c_str(),1024);
  *data_path=buf;
  return 0;
}

int commit_particles() {
    return (*SimpleXGrid).initialize(0.0);
}

int recommit_particles() {
    return (*SimpleXGrid).reinitialize();
}

int new_particle(int *id, double x,double y,double z,double rho,
                                        double flux,double xion, double uInt, double metallicity ){
                                          
                                                                                  
    long tmp_id;
    double bs;
    
    (*SimpleXGrid).get_sizeBox(&bs);
    if(bs==0) return -2;
    x=(x/bs)+0.5;y=y/bs+0.5;z=z/bs+0.5;
    if (x<0 || x>1 || 
        y<0 || y>1 ||
        z<0 || z>1 ) 
        {
            return -3;
        }
    if((*SimpleXGrid).get_syncflag() == 0)
        return -1;
    *id = lastid++;
    tmp_id = *id;
    if((*SimpleXGrid).get_syncflag() == -1) {
        return (*SimpleXGrid).add_vertex(&tmp_id, x, y, z, rho, flux, xion, uInt, metallicity);
    }
    if((*SimpleXGrid).get_syncflag() == 1) {
        return (*SimpleXGrid).add_site(&tmp_id, x, y, z, rho, flux, xion, uInt, metallicity);
    }
    return -1;
}

int delete_particle(int id) {
 int returnvalue,returnvalues;
 if((*SimpleXGrid).get_syncflag() != 1) return -1;
 returnvalue=(*SimpleXGrid).remove_site(id);
 MPI::COMM_WORLD.Reduce(&returnvalue,&returnvalues,1,MPI::INT,MPI::SUM,0);
 MPI::COMM_WORLD.Barrier();
 return returnvalues-1;
}
  
//int evolve(double t_target,int sync) {
 //return (*SimpleXGrid).evolve(t_target, sync);
//}
int evolve_model(double t_target) {
 return (*SimpleXGrid).evolve(t_target, 1);
}

int get_time(double *t){
    return (*SimpleXGrid).get_time(t);
}

int set_time(double t){
    return (*SimpleXGrid).set_time(t);
}

int get_state(int id, double *x, double *y, double *z, double *rho,
                                           double *flux, double *xion, double *uInt, double *metallicity){
    double fx=0.0, fy=0.0, fz=0.0, frho=0.0, fflux=0.0, fxion=0.0, fuInt=0.0, fmetallicity=0.0;
    double send[8], recv[8];
    int ret, totalret=0;
    double bs;
    
    (*SimpleXGrid).get_sizeBox(&bs);
    if(bs==0) return -2;  
    
    ret=(*SimpleXGrid).get_site(id, &fx, &fy, &fz, &frho, &fflux, &fxion, &fuInt, &fmetallicity);
    MPI::COMM_WORLD.Reduce(&ret,&totalret,1,MPI::INT,MPI::SUM,0); 
    send[0]=fx;send[1]=fy;send[2]=fz;send[3]=frho;send[4]=fflux;send[5]=fxion;send[6]=fuInt;send[7]=fmetallicity;
    MPI::COMM_WORLD.Reduce(&send[0],&recv[0],8,MPI::DOUBLE,MPI::SUM,0); 
    MPI::COMM_WORLD.Barrier();
    fx=recv[0];fy=recv[1];fz=recv[2];
    frho=recv[3];
    fflux=recv[4];
    fxion=recv[5];
    fuInt=recv[6];
    fmetallicity=recv[7];
    *x=(fx-0.5)*bs;*y=(fy-0.5)*bs;*z=(fz-0.5)*bs;
    *rho=frho;
    *flux=fflux;
    *xion=fxion;
    *uInt=fuInt;
    *metallicity=fmetallicity;
    return totalret-1;
}

int get_position(int id, double *x, double *y, double *z){
    double fx=0.0, fy=0.0, fz=0.0;
    double send[3], recv[3];
    int ret, totalret;
    double bs;
    
    (*SimpleXGrid).get_sizeBox(&bs);
    if(bs==0) return -2;     
    
    ret = (*SimpleXGrid).get_position(id, &fx, &fy, &fz);
    MPI::COMM_WORLD.Reduce(&ret, &totalret, 1, MPI::INT, MPI::SUM, 0);
    send[0] = fx;
    send[1] = fy;
    send[2] = fz;
    MPI::COMM_WORLD.Reduce(&send[0], &recv[0], 3, MPI::DOUBLE, MPI::SUM, 0);
    MPI::COMM_WORLD.Barrier();
    *x = (recv[0]-0.5)*bs;
    *y = (recv[1]-0.5)*bs;
    *z = (recv[2]-0.5)*bs;
    return totalret-1;
}

int get_density(int id, double *rho){
    double frho=0.0;
    double send, recv;
    int ret, totalret;
    
    ret = (*SimpleXGrid).get_density(id, &frho);
    MPI::COMM_WORLD.Reduce(&ret, &totalret, 1, MPI::INT, MPI::SUM, 0);
    send = frho;
    MPI::COMM_WORLD.Reduce(&send, &recv, 1, MPI::DOUBLE, MPI::SUM, 0);
    MPI::COMM_WORLD.Barrier();
    *rho = recv;
    return totalret-1;
}

int get_flux(int id, double *flux){
    double fflux=0.0;
    double send, recv;
    int ret, totalret;
    
    ret = (*SimpleXGrid).get_flux(id, &fflux);
    MPI::COMM_WORLD.Reduce(&ret, &totalret, 1, MPI::INT, MPI::SUM, 0);
    send = fflux;
    MPI::COMM_WORLD.Reduce(&send, &recv, 1, MPI::DOUBLE, MPI::SUM, 0);
    MPI::COMM_WORLD.Barrier();
    *flux = recv;
    return totalret-1;
}

int get_mean_intensity(int id, double *mean_intensity){
    double fmean_intensity = 0.0;
    double send, recv;
    int ret, totalret;
    
    ret = (*SimpleXGrid).get_mean_intensity(id, &fmean_intensity);
    MPI::COMM_WORLD.Reduce(&ret, &totalret, 1, MPI::INT, MPI::SUM, 0);
    send = fmean_intensity;
    MPI::COMM_WORLD.Reduce(&send, &recv, 1, MPI::DOUBLE, MPI::SUM, 0);
    MPI::COMM_WORLD.Barrier();
    *mean_intensity = recv;
    return totalret-1;
}

int get_diffuse_intensity(int id, double *diffuse_intensity){
    double fdiffuse_intensity = 0.0;
    double send, recv;
    int ret, totalret;
    
    ret = (*SimpleXGrid).get_diffuse_intensity(id, &fdiffuse_intensity);
    MPI::COMM_WORLD.Reduce(&ret, &totalret, 1, MPI::INT, MPI::SUM, 0);
    send = fdiffuse_intensity;
    MPI::COMM_WORLD.Reduce(&send, &recv, 1, MPI::DOUBLE, MPI::SUM, 0);
    MPI::COMM_WORLD.Barrier();
    *diffuse_intensity = recv;
    return totalret-1;
}

int get_ionisation(int id, double *xion){
    double fxion=0.0;
    double send, recv;
    int ret, totalret;
    
    ret = (*SimpleXGrid).get_ionisation(id, &fxion);
    MPI::COMM_WORLD.Reduce(&ret, &totalret, 1, MPI::INT, MPI::SUM, 0);
    send = fxion;
    MPI::COMM_WORLD.Reduce(&send, &recv, 1, MPI::DOUBLE, MPI::SUM, 0);
    MPI::COMM_WORLD.Barrier();
    *xion = recv;
    return totalret-1;
}

int get_metallicity(int id, double *metallicity){
    double fmetallicity=0.0;
    double send, recv;
    int ret, totalret;
    
    ret = (*SimpleXGrid).get_metallicity(id, &fmetallicity);
    MPI::COMM_WORLD.Reduce(&ret, &totalret, 1, MPI::INT, MPI::SUM, 0);
    send = fmetallicity;
    MPI::COMM_WORLD.Reduce(&send, &recv, 1, MPI::DOUBLE, MPI::SUM, 0);
    MPI::COMM_WORLD.Barrier();
    *metallicity = recv;
    return totalret-1;
}

int get_internal_energy(int id, double *uInt){
  double fuInt=0.0;
  double send, recv;
  int ret, totalret;
  
  ret = (*SimpleXGrid).get_internalEnergy(id, &fuInt);
  MPI::COMM_WORLD.Reduce(&ret, &totalret, 1, MPI::INT, MPI::SUM, 0);
  send = fuInt;
  MPI::COMM_WORLD.Reduce(&send, &recv, 1, MPI::DOUBLE, MPI::SUM, 0);
  MPI::COMM_WORLD.Barrier();
  *uInt= recv;
  return totalret-1;
}

int get_dinternal_energy_dt(int id, double *dudt){
  double fdudt=0.0;
  double send, recv;
  int ret, totalret;
  
  ret = (*SimpleXGrid).get_dinternalEnergydt(id, &fdudt);
  MPI::COMM_WORLD.Reduce(&ret, &totalret, 1, MPI::INT, MPI::SUM, 0);
  send = fdudt;
  MPI::COMM_WORLD.Reduce(&send, &recv, 1, MPI::DOUBLE, MPI::SUM, 0);
  MPI::COMM_WORLD.Barrier();
  *dudt= recv;
  return totalret-1;
}


int set_state(int id, double x, double y, double z, double rho,
                                           double flux, double xion, double uInt, double metallicity){
    int ret,totalret;
    double bs;
    
    (*SimpleXGrid).get_sizeBox(&bs);
    if(bs==0) return -2;
    x=(x/bs)+0.5;y=y/bs+0.5;z=z/bs+0.5;
    if (x<0 || x>1 || 
        y<0 || y>1 ||
        z<0 || z>1 ) return -3;
          
    ret = (*SimpleXGrid).set_site(id, x, y, z, rho, flux, xion, uInt, metallicity);
    MPI::COMM_WORLD.Reduce(&ret, &totalret, 1, MPI::INT, MPI::SUM, 0);
    MPI::COMM_WORLD.Barrier();
    return totalret-1;
}

int set_position(int id, double x, double y, double z){
    int ret, totalret;
    double bs;
    
    (*SimpleXGrid).get_sizeBox(&bs);    
    if(bs==0) return -2;
    x=(x/bs)+0.5;y=y/bs+0.5;z=z/bs+0.5;
    if (x<0 || x>1 || 
        y<0 || y>1 ||
        z<0 || z>1 ) return -3;
        
    ret = (*SimpleXGrid).set_position(id, x, y, z);
    MPI::COMM_WORLD.Reduce(&ret, &totalret, 1, MPI::INT, MPI::SUM, 0);
    MPI::COMM_WORLD.Barrier();
    return totalret-1;
}

int set_density(int id, double rho){
    int ret, totalret;
    
    ret = (*SimpleXGrid).set_density(id, rho);
    MPI::COMM_WORLD.Reduce(&ret, &totalret, 1, MPI::INT, MPI::SUM, 0);
    MPI::COMM_WORLD.Barrier();
    return totalret-1;
}

int set_flux(int id, double flux){
    int ret, totalret;
    
    ret = (*SimpleXGrid).set_flux(id, flux);
    MPI::COMM_WORLD.Reduce(&ret, &totalret, 1, MPI::INT, MPI::SUM, 0);
    MPI::COMM_WORLD.Barrier();
    return totalret-1;
}

int set_ionisation(int id, double xion){
    int ret, totalret;
    
    ret = (*SimpleXGrid).set_ionisation(id, xion);
    MPI::COMM_WORLD.Reduce(&ret, &totalret, 1, MPI::INT, MPI::SUM, 0);
    MPI::COMM_WORLD.Barrier();
    return totalret-1;
}

int set_metallicity(int id, double metallicity){
    int ret, totalret;
    
    ret = (*SimpleXGrid).set_metallicity(id, metallicity);
    MPI::COMM_WORLD.Reduce(&ret, &totalret, 1, MPI::INT, MPI::SUM, 0);
    MPI::COMM_WORLD.Barrier();
    return totalret-1;
}

int set_internal_energy(int id, double uInt){
    int ret, totalret;
    
    ret = (*SimpleXGrid).set_internalEnergy(id, uInt);
    MPI::COMM_WORLD.Reduce(&ret, &totalret, 1, MPI::INT, MPI::SUM, 0);
    MPI::COMM_WORLD.Barrier();
    return totalret-1;
}

int set_dinternal_energy_dt(int id, double dut){
    int ret, totalret;
    
    ret = (*SimpleXGrid).set_dinternalEnergydt(id, dut);
    MPI::COMM_WORLD.Reduce(&ret, &totalret, 1, MPI::INT, MPI::SUM, 0);
    MPI::COMM_WORLD.Barrier();
    return totalret-1;
}


int cleanup_code(void){
 (*SimpleXGrid).clear_temporary();
 return 0;
}

int set_box_size(double bs){
  return (*SimpleXGrid).set_sizeBox(bs);
}

int get_box_size(double *bs){
  return (*SimpleXGrid).get_sizeBox(bs);
}

int set_timestep(double ts){
  return (*SimpleXGrid).set_UNIT_T(ts);
}

int get_timestep(double *ts){
  return (*SimpleXGrid).get_UNIT_T(ts);
}

int set_hilbert_order(int ho){
  return (*SimpleXGrid).set_hilbert_order(ho);
}

int get_hilbert_order(int *ho){
  return (*SimpleXGrid).get_hilbert_order(ho);
}


int set_number_frequency_bins(int ts){
  return (*SimpleXGrid).set_numFreq(ts);
}

int get_number_frequency_bins(int *ts){
  return (*SimpleXGrid).get_numFreq(ts);
}


int set_number_of_border_sites(int ts){
  return SimpleXGrid->set_numBorderSites(ts);
}

int get_number_of_border_sites(int *ts){
  return SimpleXGrid->get_numBorderSites(ts);
}


int set_thermal_evolution(int ts){
  return (*SimpleXGrid).set_heat_cool(ts);
}

int get_thermal_evolution(int *ts){
  return (*SimpleXGrid).get_heat_cool(ts);
}

int set_metal_cooling(int ts){
  return (*SimpleXGrid).set_metal_cooling(ts);
}

int get_metal_cooling(int *ts){
  return (*SimpleXGrid).get_metal_cooling(ts);
}

int set_recombination_radiation(int ts){
  return (*SimpleXGrid).set_recombination_radiation(ts);
}

int get_recombination_radiation(int *ts){
  return (*SimpleXGrid).get_recombination_radiation(ts);
}


int set_source_Teff(double ts){
  return (*SimpleXGrid).set_sourceTeff(ts);
}

int get_source_Teff(double *ts){
  return (*SimpleXGrid).get_sourceTeff(ts);
}

int set_collisional_ionization(int ts){
  return (*SimpleXGrid).set_coll_ion(ts);
}

int get_collisional_ionization(int *ts){
  return (*SimpleXGrid).get_coll_ion(ts);
}

int set_blackbody_spectrum(int ts){
  return (*SimpleXGrid).set_blackBody(ts);
}

int get_blackbody_spectrum(int *ts){
  return (*SimpleXGrid).get_blackBody(ts);
}




int commit_parameters(){
    return 0;
}

int recommit_parameters(){
    return 0;
}
