#include "Main.h"
#include "configfile.h"

#include "amuse_interface.h"

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
  UNIT_T = 0.05;
  //size of simulation domain in pc
  sizeBox = 13200;

  //use blackbody for source? Otherwise monochromatic
  blackBody = 0;
  //effective temperature of source
  sourceTeff = 1.e5;
  
  //units of source
  UNIT_I = 1.e48;
  

//----------- Don't change these unless you know what you're doing -------------------------//

  //fill choice should always be read
  fillChoice = READ;
  //dimension should always be 3
  dimension = 3;
  //periodic boundaries don't work in parallel
  periodic = 0;
  
  //number of points in border, st this as you like
  borderSites = 25000;
  //region in which border points are places
  borderBox = 0.1;
  //size of border around subbox
  padding_subbox = 0.25;
  
  //temperature of the ionised gas
  gasIonTemp = 1.e4;
  
  //source ionises its own cell or not
  source_inside_cell = 0;
  
  //include recombinations?
  recombination = 1;

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

  // Do cyclic check?
  cyclic_check = 0;
  
  //Use dynamical updates?
  updates = 0;
  //Number of updates
  numGridDyn = 100;
  //Switch between dct and ballistic transport
  switchTau = 0.5;
  //minimum resolution 
  minResolution = 32;
  //number of bins in update routine
  nbins = 100;  
  
  //Use temporal photon conservation?
  photon_conservation = 1; 

  //Number of directions
  number_of_directions = 42;
  //Number of orientations, fixed in header files
  number_of_orientations = 100;
  
  //calculate straight from tesselation instead of direction bins
  //is faster in postprocessing, but less precise
  straight_from_tess = 1;

//----------------------------------------------------------------------//


}

//add one vertex to the list of vertices that needs to be triangulated
int AMUSE_SimpleX::add_vertex(long *id, double x,double y,double z,double rho,
                                           double flux,double xion){

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


  return 0;

}

//add a site to the sites list
int AMUSE_SimpleX::add_site(long *id, double x,double y,double z,double rho,
                                           double flux,double xion){

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
    tempSite.set_n_HI( (float) n_HII );
    tempSite.set_flux( flux );
    tempSite.set_neigh_dist(1.);
    sites.push_back( tempSite );
    numSites++;
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
    
  double dt= t_target*secondsPerMyr -total_time;
    
  if(syncflag==1){
    reinitialize();
  }
    
  simpleXlog << endl << "  start sweeping" << endl;  
    
  if(dt > 0){
    numSweeps=dt/UNIT_T+1;
    printf("proc %d working\n",COMM_RANK);
    radiation_transport(1);
    total_time+=numSweeps*UNIT_T;
  }
  
  simpleXlog << endl << "  done sweeping" << endl;
  
  if(sync == 1) {
    syncflag=1;     
  }
  
  return 0;

}

//store everything in the simulation and recompute triangulation and physics
int AMUSE_SimpleX::reinitialize(){

   //make sure that the vectors that will be filled are empty 
    site_intensities.clear();
    intens_ids.clear();
    
    //make sure that photons have left all ghost vertices 
    //before calculating local intensities
    send_intensities();

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
      create_vertex_tree();
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
    
  syncflag=0;
  return 0;
}

//return properties of specified site
int AMUSE_SimpleX::get_site(int id, double *x,double *y,double *z,double *rho,
                                              double *flux,double *xion){
   SITE_ITERATOR p;
   Site tmp;
   
   tmp.set_vertex_id((unsigned long long) id);
   p=lower_bound(sites.begin(),sites.end(), tmp, compare_vertex_id_site);
   if(p->get_vertex_id() == (unsigned long long int)id){
     if (p->get_process() == COMM_RANK){
       *x=p->get_x();
       *y=p->get_y();
       *z=p->get_z();
       double tmp = (double) p->get_n_HI() + (double) p->get_n_HII();
       *rho = tmp;
       *flux = p->get_flux();
       *xion = (double) p->get_n_HII()/tmp;
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
            *flux = p->get_flux();
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

//set properties of specified site
int AMUSE_SimpleX::set_site(int id, double x, double y, double z, double rho,
                                              double flux, double xion){
    SITE_ITERATOR p;
    Site tmp;
    
    tmp.set_vertex_id((unsigned long long) id);
    p = lower_bound(sites.begin(), sites.end(), tmp, compare_vertex_id_site);
    if(p->get_vertex_id() == (unsigned long long int)id){
        p->set_x(x);
        p->set_y(y);
        p->set_z(z);
        p->set_flux(flux);
        p->set_n_HI((1-xion)*rho);
        p->set_n_HII(xion*rho);
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
        p->set_flux(flux);
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

//relevant global parameters for run
AMUSE_SimpleX *SimpleXGrid;
int COMM_RANK;
int lastid = 0;
string global_output_path = ".";

int initialize_code() {
  COMM_RANK = MPI::COMM_WORLD.Get_rank();
  SimpleXGrid=new AMUSE_SimpleX(global_output_path);
  return (*SimpleXGrid).setup_simplex();  
}

int set_output_directory(char *output_path){
    global_output_path = output_path;
    return 0;
}

int commit_particles() {
    return (*SimpleXGrid).initialize(0.0);
}

int recommit_particles() {
    return (*SimpleXGrid).reinitialize();
}

int new_particle(int *id, double x,double y,double z,double rho,
                                        double flux,double xion){
    long tmp_id;
    double bs;
    
    (*SimpleXGrid).get_sizeBox(&bs);
    if(bs==0) return -2;
    x/=bs;y/=bs;z/=bs;
    if (x<0 || x>1 || 
        y<0 || y>1 ||
        z<0 || z>1 ) return -3;
    
    if((*SimpleXGrid).get_syncflag() == 0)
        return -1;
    *id = lastid++;
    tmp_id = *id;
    if((*SimpleXGrid).get_syncflag() == -1)
        return (*SimpleXGrid).add_vertex(&tmp_id, x, y, z, rho, flux, xion);
    if((*SimpleXGrid).get_syncflag() == 1)
        return (*SimpleXGrid).add_site(&tmp_id, x, y, z, rho, flux, xion);
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
int evolve(double t_target) {
 return (*SimpleXGrid).evolve(t_target, 1);
}

int get_state(int id, double *x, double *y, double *z, double *rho,
                                           double *flux, double *xion){
    double fx=0.0, fy=0.0, fz=0.0, frho=0.0, fflux=0.0, fxion=0.0;
    double send[6], recv[6];
    int ret, totalret=0;
    double bs;
    
    (*SimpleXGrid).get_sizeBox(&bs);
    if(bs==0) return -2;  
    
    ret=(*SimpleXGrid).get_site(id, &fx, &fy, &fz, &frho, &fflux, &fxion);
    MPI::COMM_WORLD.Reduce(&ret,&totalret,1,MPI::INT,MPI::SUM,0); 
    send[0]=fx;send[1]=fy;send[2]=fz;send[3]=frho;send[4]=fflux;send[5]=fxion;
    MPI::COMM_WORLD.Reduce(&send[0],&recv[0],6,MPI::DOUBLE,MPI::SUM,0); 
    MPI::COMM_WORLD.Barrier();
    fx=recv[0];fy=recv[1];fz=recv[2];frho=recv[3];fflux=recv[4];fxion=recv[5];
    *x=fx*bs;*y=fy*bs;*z=fz*bs;*rho=frho;*flux=fflux;*xion=fxion;
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
    *x = recv[0]*bs;
    *y = recv[1]*bs;
    *z = recv[2]*bs;
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

int set_state(int id, double x, double y, double z, double rho,
                                           double flux, double xion){
    int ret,totalret;
    double bs;
    
    (*SimpleXGrid).get_sizeBox(&bs);
    if(bs==0) return -2;
    x/=bs;y/=bs;z/=bs;
    if (x<0 || x>1 || 
        y<0 || y>1 ||
        z<0 || z>1 ) return -3;
          
    ret = (*SimpleXGrid).set_site(id, x, y, z, rho, flux, xion);
    MPI::COMM_WORLD.Reduce(&ret, &totalret, 1, MPI::INT, MPI::SUM, 0);
    MPI::COMM_WORLD.Barrier();
    return totalret-1;
}

int set_position(int id, double x, double y, double z){
    int ret, totalret;
    double bs;
    
    (*SimpleXGrid).get_sizeBox(&bs);    
    if(bs==0) return -2;
    x/=bs;y/=bs;z/=bs;
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

int cleanup_code(void){
 (*SimpleXGrid).clear_temporary();
 return 0;
}

int set_box_size_parameter(double bs){
  return (*SimpleXGrid).set_sizeBox(bs);
}

int get_box_size_parameter(double *bs){
  return (*SimpleXGrid).get_sizeBox(bs);
}

int set_timestep_parameter(double ts){
  return (*SimpleXGrid).set_UNIT_T(ts);
}

int get_timestep_parameter(double *ts){
  return (*SimpleXGrid).get_UNIT_T(ts);
}

int set_hilbert_order_parameter(int ho){
  return (*SimpleXGrid).set_hilbert_order(ho);
}

int get_hilbert_order_parameter(int *ho){
  return (*SimpleXGrid).get_hilbert_order(ho);
}


int commit_parameters(){
    return 0;
}

int recommit_parameters(){
    return 0;
}
