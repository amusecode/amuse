#include "octree.h"

#define NMAXSAMPLE 20000



void octree::mpiInit(int argc,char *argv[], int &procId, int &nProcs)
{
    int  namelen;   
    char processor_name[MPI_MAX_PROCESSOR_NAME];

    int bla;
    MPI_Initialized(&bla);    
    if(bla)
    {
	    cout << "Mpi al geinit\n";
    }
    else
    {
	    MPI_Init(&argc,&argv);
    }


    MPI_Comm_size(MPI_COMM_WORLD, &nProcs);
    MPI_Comm_rank(MPI_COMM_WORLD, &procId);
    MPI_Get_processor_name(processor_name,&namelen);

    #ifdef PRINT_MPI_DEBUG
      fprintf(stderr, "Proc id: %d @ %s , total processes: %d (mpiInit) \n", procId, processor_name, nProcs);
    #endif


//TODO write delete functions in destructor
    
    //Allocate memory for the used buffers
    cur_xlow  = new double4[nProcs];
    cur_xhigh = new double4[nProcs];  
    
    cur_xlow_xhigh = new double4[nProcs*2];
    
    xlow      = new double4[nProcs];
    xhigh     = new double4[nProcs];
    
    let_xlow      = new double4[nProcs];
    let_xhigh     = new double4[nProcs];    
}



//Utility functions
void octree::mpiSync(){
  MPI_Barrier(MPI_COMM_WORLD);
}

int octree::mpiGetRank(){ 
  return procId;
}

int octree::mpiGetNProcs(){
  return nProcs;
}

void octree::AllSum(double &value)
{
  double tmp = -1;
  MPI_Allreduce(&value,&tmp,1, MPI_DOUBLE, MPI_SUM,MPI_COMM_WORLD);
  value = tmp;
}

int octree::SumOnRootRank(int &value)
{
  int temp;
  MPI_Reduce(&value,&temp,1, MPI_INT, MPI_SUM,0, MPI_COMM_WORLD);
  return temp;
}

//TODO make this workable in one communication by
//storing data one buffer
void octree::mpiRadiusFind(real4 &rmin, real4 &rmax)
{
  int n = nProcs;
  
  //Store the real4 in double4 for communication
  double4 rmin_d4 = (double4){rmin.x, rmin.y, rmin.z, rmin.w};
  double4 rmax_d4 = (double4){rmax.x, rmax.y, rmax.z, rmax.w};
  
  double4 minmax[2];
  minmax[0] = (double4){rmin.x, rmin.y, rmin.z, rmin.w};
  minmax[1] = (double4){rmax.x, rmax.y, rmax.z, rmax.w};    
  

//   cur_xlow_xhigh = new double4[nProcs*2];
//   double4 minmax[2];
  minmax[0] = (double4){rmin.x, rmin.y, rmin.z, rmin.w};
  minmax[1] = (double4){rmax.x, rmax.y, rmax.z, rmax.w};  
  MPI_Allgather(&minmax, sizeof(double4)*2, MPI_BYTE,  cur_xlow_xhigh,
                sizeof(double4)*2, MPI_BYTE, MPI_COMM_WORLD);  
                
                
                
  rmin.x = cur_xlow[0].x = cur_xlow_xhigh[0].x;
  rmin.y = cur_xlow[0].y = cur_xlow_xhigh[0].y;
  rmin.z = cur_xlow[0].z = cur_xlow_xhigh[0].z;
  cur_xlow[0].w = cur_xlow_xhigh[0].w;
  for(int i=1; i < n; i++)
  {
    rmin.x = fmin(rmin.x, cur_xlow_xhigh[i*2].x);
    rmin.y = fmin(rmin.y, cur_xlow_xhigh[i*2].y);
    rmin.z = fmin(rmin.z, cur_xlow_xhigh[i*2].z);
    
    cur_xlow[i].x = cur_xlow_xhigh[i*2].x;
    cur_xlow[i].y = cur_xlow_xhigh[i*2].y;
    cur_xlow[i].z = cur_xlow_xhigh[i*2].z;
    cur_xlow[i].w = cur_xlow_xhigh[i*2].w;
  }                
          
  rmax.x = cur_xhigh[0].x = cur_xlow_xhigh[1].x;
  rmax.y = cur_xhigh[0].y = cur_xlow_xhigh[1].y;
  rmax.z = cur_xhigh[0].z = cur_xlow_xhigh[1].z;
  cur_xhigh[0].w = cur_xlow_xhigh[1].w;
  for(int i=1; i < n; i++)
  {
    rmax.x = fmax(rmax.x, cur_xlow_xhigh[i*2+1].x);
    rmax.y = fmax(rmax.y, cur_xlow_xhigh[i*2+1].y);
    rmax.z = fmax(rmax.z, cur_xlow_xhigh[i*2+1].z);
    
    cur_xhigh[i].x = cur_xlow_xhigh[i*2+1].x;
    cur_xhigh[i].y = cur_xlow_xhigh[i*2+1].y;
    cur_xhigh[i].z = cur_xlow_xhigh[i*2+1].z;
    cur_xhigh[i].w = cur_xlow_xhigh[i*2+1].w;
  }                

      
  
/*
  MPI_Allgather(&rmin_d4, sizeof(double4), MPI_BYTE,  cur_xlow,
                sizeof(double4), MPI_BYTE, MPI_COMM_WORLD);

  rmin.x = cur_xlow[0].x;
  rmin.y = cur_xlow[0].y;
  rmin.z = cur_xlow[0].z;
  for(int i=1; i < n; i++)
  {
    rmin.x = fmin(rmin.x, cur_xlow[i].x);
    rmin.y = fmin(rmin.y, cur_xlow[i].y);
    rmin.z = fmin(rmin.z, cur_xlow[i].z);
  }

  MPI_Allgather(&rmax_d4, sizeof(double4), MPI_BYTE,  cur_xhigh,
                sizeof(double4), MPI_BYTE, MPI_COMM_WORLD);

  rmax.x = cur_xhigh[0].x;
  rmax.y = cur_xhigh[0].y;
  rmax.z = cur_xhigh[0].z;  
  for(int i=1; i < n; i++)
  {
    rmax.x = fmax(rmax.x, cur_xhigh[i].x);
    rmax.y = fmax(rmax.y, cur_xhigh[i].y);
    rmax.z = fmax(rmax.z, cur_xhigh[i].z);
  }
  
    fprintf(stderr, "Proc: %d recv: %f %f %f %f \t %f %f %f %f \n", mpiGetRank(),
            rmin.x, rmin.y, rmin.z, rmin.w,
            rmax.x, rmax.y, rmax.z, rmax.w);
            
/*
  minmax[0] = (double4){rmin.x, rmin.y, rmin.z, rmin.w};
  minmax[1] = (double4){rmax.x, rmax.y, rmax.z, rmax.w}; */ 
/*  MPI_Allgather(&minmax, sizeof(double4)*2, MPI_BYTE,  cur_xlow_xhigh,
                sizeof(double4)*2, MPI_BYTE, MPI_COMM_WORLD);  
          
  for(int i=0; i < nProcs; i++)
  {
    fprintf(stderr, "Proc2: %d recv: %f %f %f %f \t %f %f %f %f \n", mpiGetRank(),
            cur_xlow_xhigh[i*2].x, cur_xlow_xhigh[i*2].y, cur_xlow_xhigh[i*2].z, cur_xlow_xhigh[i*2].w,
            cur_xlow_xhigh[i*2+1].x, cur_xlow_xhigh[i*2+1].y, cur_xlow_xhigh[i*2+1].z, cur_xlow_xhigh[i*2+1].w);
    
  }
  
  rmin.x = cur_xlow[0].x = cur_xlow_xhigh[0].x;
  rmin.y = cur_xlow[0].y = cur_xlow_xhigh[0].y;
  rmin.z = cur_xlow[0].z = cur_xlow_xhigh[0].z;
  cur_xlow[0].w = cur_xlow_xhigh[0].w;
  for(int i=1; i < n; i++)
  {
    rmin.x = fmin(rmin.x, cur_xlow_xhigh[i*2].x);
    rmin.y = fmin(rmin.y, cur_xlow_xhigh[i*2].y);
    rmin.z = fmin(rmin.z, cur_xlow_xhigh[i*2].z);
    
    cur_xlow[i].x = cur_xlow_xhigh[i*2].x;
    cur_xlow[i].y = cur_xlow_xhigh[i*2].y;
    cur_xlow[i].z = cur_xlow_xhigh[i*2].z;
    cur_xlow[i].w = cur_xlow_xhigh[i*2].w;
  }                
          
  rmax.x = cur_xhigh[0].x = cur_xlow_xhigh[1].x;
  rmax.y = cur_xhigh[0].y = cur_xlow_xhigh[1].y;
  rmax.z = cur_xhigh[0].z = cur_xlow_xhigh[1].z;
  cur_xhigh[0].w = cur_xlow_xhigh[1].w;
  for(int i=1; i < n; i++)
  {
    rmax.x = fmax(rmax.x, cur_xlow_xhigh[i*2+1].x);
    rmax.y = fmax(rmax.y, cur_xlow_xhigh[i*2+1].y);
    rmax.z = fmax(rmax.z, cur_xlow_xhigh[i*2+1].z);
    
    cur_xhigh[i].x = cur_xlow_xhigh[i*2+1].x;
    cur_xhigh[i].y = cur_xlow_xhigh[i*2+1].y;
    cur_xhigh[i].z = cur_xlow_xhigh[i*2+1].z;
    cur_xhigh[i].w = cur_xlow_xhigh[i*2+1].w;  
  }
    
    fprintf(stderr, "Proc3: %d recv: %f %f %f %f \t %f %f %f %f \n", mpiGetRank(),
            rmin.x, rmin.y, rmin.z, rmin.w,
            rmax.x, rmax.y, rmax.z, rmax.w);    

exit(0);
                */
                /*
                Finished!!! Took in total: 19.1085 sec
Exchange_copy_h2d took: 18.682943        millisecond
iter=15 : time= 1  Etot= -0.2455571191  Ekin= 0.242635   Epot= -0.488192 : de= -6.47833e-05 ( 6.47833e-05 ) d(de)= 2.26062e-16 ( 6.8184e-06 ) t_sim=  17.7919 sec
 Finished: 1  > 1 loop alone took: 17.792
Finished!!! Took in total: 19.1085 sec
*/
//   exit(0);                
  
  
  
}

/*
  double4 tempMin;
  tempMin.x = bMin.x; tempMin.y = bMin.y; tempMin.z = bMin.z; tempMin.w = bMin.w;

  double4 tempMax;
  tempMax.x = bMax.x; tempMax.y = bMax.y; tempMax.z = bMax.z; tempMax.w = bMax.w;

  MPI_Allgather(&tempMin, sizeof(double4), MPI_BYTE,  let_xlow, sizeof(double4), MPI_BYTE, MPI_COMM_WORLD);
  MPI_Allgather(&tempMax, sizeof(double4), MPI_BYTE,  let_xhigh, sizeof(double4), MPI_BYTE, MPI_COMM_WORLD);
*/

//Main functions

//Copied from Makino code
void octree::createORB()
{      
  int n0, n1;
  n0 = (int)pow(nProcs+0.1,0.33333333333333333333);
  while(nProcs % n0)
    n0--;

  nx = n0;
  n1 = nProcs/nx;
  n0 = (int)sqrt(n1+0.1);
  while(n1 % n0)
    n0++;

  ny = n0; nz = n1/n0;
  int ntmp;
  if (nz > ny){
      ntmp = nz; nz = ny; ny = ntmp;
  }
  if (ny > nx){
      ntmp = nx; nx = ny; ny = ntmp;
  }
  if (nz > ny){
      ntmp = nz; nz = ny; ny = ntmp;
  }
  if (nx*ny*nz != nProcs){
      cerr << "create_division: Intenal Error " << nProcs << " " << nx
            << " " << ny << " " << nz <<endl;
  }
  
  #ifdef PRINT_MPI_DEBUG
    if(procId == 0) cout << "Division: nx: " << nx << " ny: " << ny << " nz: " << nz << endl;
  #endif
}


void octree::determine_sample_freq(int numberOfParticles)
{
    //Sum the number of particles on all processes
    int tmp;
    MPI_Allreduce(&numberOfParticles,&tmp,1, MPI_INT, MPI_SUM,MPI_COMM_WORLD);

    #ifdef PRINT_MPI_DEBUG
    if(procId == 0)
      cout << "Total number of particles: " << tmp << endl;
    #endif
   
    int maxsample = (int)(NMAXSAMPLE*0.8); // 0.8 is safety factor    
    sampleFreq = (tmp+maxsample-1)/maxsample;
}


void octree::collect_sample_particles(real4 *bodies,
                                      int nbody,
                                      int sample_freq,
                                      vector<real4> &sampleArray,                              
                                      int &nsample,
                                      double &rmax)
{
    //Select the sample particles    
    int ii, i;     
    for(i = ii= 0;ii<nbody; i++,ii+=sample_freq)
    {
        sampleArray.push_back(bodies[ii]);
    }
    nsample = i;

    //Now gather the particles at process 0
    //NOTE: Im using my own implementation instead of makino's which
    //can grow out of the memory array (I think...)

    //Sum the total amount of sample particles
    unsigned int totalNumberOfSamples;
    MPI_Reduce(&nsample,&totalNumberOfSamples,1, MPI_INT, MPI_SUM,0, MPI_COMM_WORLD);
    
    //Increase the size of the result buffer if needed
    if(procId == 0)
    {
      if(totalNumberOfSamples > sampleArray.size())
      {
        sampleArray.resize(totalNumberOfSamples);
      }
    }

    //Send an receive the sample data
    MPI_Status status;
    if(procId != 0)
    {
        //First send the number of samples, then the actual sample data
        MPI_Send( &nsample, 1, MPI_INT, 0, procId*2 , MPI_COMM_WORLD);
        MPI_Send( (real*)&sampleArray[0], nsample*sizeof(real)*4, MPI_BYTE, 0, procId*2+1, MPI_COMM_WORLD);
    }
    else
    {
        for(int i=1; i<nProcs; i++)
        {
            int nreceive;    
            MPI_Recv( &nreceive, 1, MPI_INT, i,i*2, MPI_COMM_WORLD,&status);            
            MPI_Recv((real*)(&sampleArray[nsample]), sizeof(real)*4*nreceive, MPI_BYTE, i,i*2+1, MPI_COMM_WORLD,&status);

            nsample+=nreceive;
        }
    }
    
    //TODO, NOTE this can be extracted after the update position / compute properties on the GPU
    //Find the maximum particle position 
    double tmp = 0;
    for(i = 0;i<nbody; i++)
    {
        real4 r = bodies[i];
        //check x,y and z       
        if(fabs(r.x)>tmp)  tmp=fabs(r.x);
        if(fabs(r.y)>tmp)  tmp=fabs(r.y);
        if(fabs(r.z)>tmp)  tmp=fabs(r.z);
    }
    
    //Find the global maximum
    MPI_Allreduce(&tmp, &rmax,1, MPI_DOUBLE, MPI_MAX,MPI_COMM_WORLD);
}

void octree::createDistribution(real4 *bodies, int n_bodies)
{
  determine_sample_freq(n_bodies);

  vector<real4> sampleArray;
  sampleArray.reserve(NMAXSAMPLE);

  int     nsample;  //Number of samples for this process
  double  rmax;     //TODO: wat doet dit precies...

  //Get the sample particles from the other processes
  collect_sample_particles(bodies, n_bodies, sampleFreq, sampleArray, nsample, rmax);
 
  //Now that we have a sample from all proces we setup the space division
  //Processor 0 determines the division
  if(procId == 0)
    determine_division(nsample, sampleArray,nx, ny, nz, rmax,xlow, xhigh);
  
  if(procId == 0)
  {
//     xhigh[0].z = 0;
//     xlow[1].x = 0;
//     
//     printf("TEST low: %f \t %f \t %f \t high: %f \t %f \t %f\n", xlow[0].x, xlow[0].y, xlow[0].z, xhigh[0].x, xhigh[0].y, xhigh[0].z);
//     printf("TEST 2 low: %f \t %f \t %f \t high: %f \t %f \t %f\n", xlow[1].x, xlow[1].y, xlow[1].z, xhigh[1].x, xhigh[1].y, xhigh[1].z);
  }
  //Now broadcast the results to all other processes
  MPI_Bcast(xlow,  sizeof(double4)*nProcs,MPI_BYTE,0,MPI_COMM_WORLD);
  MPI_Bcast(xhigh, sizeof(double4)*nProcs,MPI_BYTE,0,MPI_COMM_WORLD);


  return;
}

//Updates the box dimensions
//TODO how is this different from createDistribution?
void octree::updateDistribution(real4 *bodies, int n_bodies)
{  
  determine_sample_freq(n_bodies);
  
  vector<real4> sampleArray;
  sampleArray.reserve(NMAXSAMPLE);  
  
  int     nsample;                                 //Number of samples for this process
  double  rmax;                                    //TODO: wat doet dit precies...
 
  collect_sample_particles(bodies, n_bodies, sampleFreq, sampleArray, nsample, rmax);
 
  //Now that we have a sample from all proces we setup the space division

  //Processor 0 determines the division
  if(procId == 0)
    determine_division(nsample, sampleArray,nx, ny, nz, rmax,xlow, xhigh);

  //Now broadcast the results to all other processes
  MPI_Bcast(xlow,  sizeof(double4)*nProcs,MPI_BYTE,0,MPI_COMM_WORLD);
  MPI_Bcast(xhigh, sizeof(double4)*nProcs,MPI_BYTE,0,MPI_COMM_WORLD);


  return;
}



//Sort function based on Makinos function
//Sorts (a part) of the coordinate array
//containing the sample particles
//Either sorts the x,y or z direction
//lo is the lower bound of the to sorted part
//up is the upper bound of the to sorted part
//cid is the index/axes to sort
//cid=0=x, cid=1=y and cid=2=z
void octree::sortCoordinates(real4 *r, int lo, int up, int cid )
{
    int i, j;
    real4 tempr;
    while ( up>lo ) {
        i = lo;
        j = up;
        tempr = r[lo];
        /*** Split file in two ***/
        while ( i<j )
        {
            if(cid==0)
              for ( ; r[j].x > tempr.x; j-- );
            else if(cid==1)
              for ( ; r[j].y > tempr.y; j-- );
            else
              for ( ; r[j].z > tempr.z; j-- );

            if(cid==0)
              for ( r[i]=r[j]; i<j && r[i].x <= tempr.x; i++ );
            else if(cid==1)
              for ( r[i]=r[j]; i<j && r[i].y <= tempr.y; i++ );
            else
              for ( r[i]=r[j]; i<j && r[i].z <= tempr.z; i++ );
            
            r[j] = r[i];
        }
        r[i] = tempr;
        /*** Sort recursively, the smallest first ***/
        if ( i-lo < up-i ) 
        { 
          sortCoordinates(r,lo,i-1,cid);  
          lo = i+1; 
        }
        else
        { 
          sortCoordinates(r,i+1,up,cid);
          up = i-1; 
        }
    }
}


//Calculates the dimension of the box
//np number of sample particles
//pos the sample particle positions
//cid the coordinate index, 0=x, 1=y, 2=z
//istart/iend the start and end position of the sorted array
//rmax the maximum coordinate
//xlow/xhigh the box coordinates
void octree::calculate_boxdim(int np, real4 pos[], int cid, int istart, int iend,
                      double rmax, double & xlow, double & xhigh)
{
    if(istart == 0) 
    {
        xlow = -rmax;
    }
    else
    {
      if(cid==0)
        xlow = (pos[istart].x + pos[istart-1].x)/2;
      else if(cid==1)
        xlow = (pos[istart].y + pos[istart-1].y)/2;
      else
        xlow = (pos[istart].z + pos[istart-1].z)/2;
    }

    if(iend == np-1)
    {
        xhigh = rmax;
    }
    else
    {
      if(cid==0)
        xhigh = (pos[iend].x + pos[iend+1].x)/2;
      else if(cid==1)
        xhigh = (pos[iend].y + pos[iend+1].y)/2;
      else
        xhigh = (pos[iend].z + pos[iend+1].z)/2;
    }
}


void octree::determine_division(int np,         // number of particles
                        vector<real4> &pos,     // positions of particles
                        int nx,
                        int ny,
                        int nz,
                        double rmax,
                        double4 xlow[],         // left-bottom coordinate of divisions
                        double4 xhigh[])        // size of divisions
{
    int numberOfProcs = nProcs;
    int *istart  = new int[numberOfProcs+1];
    int *iend    = new int[numberOfProcs+1];
    int n = nx*ny*nz;

    sortCoordinates(&pos[0], 0, np-1, 0);
    
    //Split the array in more or less equal parts
    //TODO Waarom is dit precies
    for(int i = 0;i<n;i++)
    {
        istart[i] = (i*np)/n;
         //TODO: het was, maar schrijft voor array...  if(i>=0) 
        if(i > 0 )
         iend[i-1]=istart[i]-1;
    }
    iend[n-1] = np-1; 
 
    //Split the x-axis
    for(int ix = 0;ix<nx;ix++)
    {
        double x0, x1;
        int ix0 = ix*ny*nz;
        int ix1 = (ix+1)*ny*nz;       
        calculate_boxdim(np, &pos[0], 0,istart[ix0],iend[ix1-1],rmax,x0,x1);                
        for(int i=ix0; i<ix1; i++)
        {
            xlow[i].x   = x0;
            xhigh[i].x  = x1;
        }
    }

    //For each x split the various y parts
    for(int ix = 0;ix<nx;ix++)
    {
        int ix0 = ix*ny*nz;
        int ix1 = (ix+1)*ny*nz;
        int npy = iend[ix1-1] - istart[ix0] + 1;
        sortCoordinates(&pos[0], istart[ix0],iend[ix1-1], 1);
        for(int iy = 0;iy<ny;iy++){
            double y0, y1;
            int iy0 = ix0+iy*nz;
            int iy1 = ix0+(iy+1)*nz;
            calculate_boxdim(npy, &pos[istart[ix0]], 1,istart[iy0]-istart[ix0],
                             iend[iy1-1]-istart[ix0], rmax, y0,y1);
            for(int i=iy0; i<iy1; i++)
            {
                xlow[i].y  = y0;
                xhigh[i].y = y1;
            }
        }
    }
    
    //For each x and for each y split the z axis
    for(int ix = 0;ix<nx;ix++){
        int ix0 = ix*ny*nz;
        for(int iy = 0;iy<ny;iy++){
            int iy0 = ix0+iy*nz;
            int iy1 = ix0+(iy+1)*nz;
            int npz = iend[iy1-1] - istart[iy0] + 1;
            sortCoordinates(&pos[0], istart[iy0],iend[iy1-1], 2);          
            for(int iz = 0;iz<nz;iz++){
                double z0, z1;
                int iz0 = iy0+iz;
                calculate_boxdim(npz, &pos[istart[iy0]], 2,istart[iz0]-istart[iy0],
                                 iend[iz0]-istart[iy0], rmax, z0,z1);
                xlow[iz0].z   = z0;
                xhigh[iz0].z  = z1;
            }
        }
    }


//TODO: Why do we get errors if we use this?
  delete[] istart;
  delete[] iend;

  return;
}


template<class T>
int octree::MP_exchange_particle_with_overflow_check(int ibox,
                                                    T *source_buffer,
                                                    vector<T> &recv_buffer,
                                                    int firstloc,
                                                    int nparticles,
                                                    int isource,                                                   
                                                    int &nsend,
                                                    unsigned int &recvCount)
{
    int iret = 0;
    MPI_Status status;
    int local_proc_id = procId;
    nsend = nparticles;
    
    //first send&get the number of particles to send&get
    unsigned int nreceive;
    
    //Send and get the number of particles that are exchanged
    MPI_Sendrecv(&nsend,1,MPI_INT,ibox,local_proc_id*10,
                 &nreceive,1,MPI_INT,isource,isource*10,MPI_COMM_WORLD,
                 &status);                

    int ss         = sizeof(T);
    int sendoffset = nparticles-nsend;
    
    //Resize the receive buffer    
    if((nreceive + recvCount) > recv_buffer.size())      
    {      
      recv_buffer.resize(nreceive + recvCount);      
    }

    
    //Send the actual particles
    MPI_Sendrecv(&source_buffer[firstloc+sendoffset],ss*nsend,MPI_BYTE,ibox,local_proc_id*10+1,
                 &recv_buffer[recvCount],ss*nreceive,MPI_BYTE,isource,isource*10+1,
                 MPI_COMM_WORLD,&status);  
                 
    recvCount += nreceive;

    int giret;
    MPI_Allreduce(&iret, &giret,1, MPI_INT, MPI_MAX,MPI_COMM_WORLD);
    return giret;
} 


//Checks if the position falls within the specified box
inline int isinbox(real4 pos, double4 xlow, double4 xhigh)
{  
    if((pos.x < xlow.x)||(pos.x > xhigh.x))          
      return 0;
    if((pos.y < xlow.y)||(pos.y > xhigh.y))          
      return 0;
    if((pos.z < xlow.z)||(pos.z > xhigh.z))          
      return 0;
    
    return 1;
}


//TODO: UItleg
//waarschijnlijk wordt voor elke particle bepaald naar welke processor hij gestuurd moet worden
int octree::exchange_particles_with_overflow_check(tree_structure &tree)
{
  int myid      = procId;
  int nproc     = nProcs;
  int iloc      = 0;
  int totalsent = 0;
  int nbody     = tree.n;
  
    
  real4  *bodiesPositions = &tree.bodies_pos[0];
  real4  *velocities      = &tree.bodies_vel[0];  
  real4  *bodiesAcc0      = &tree.bodies_acc0[0];
  real4  *bodiesAcc1      = &tree.bodies_acc1[0];
  float2 *bodiesTime      = &tree.bodies_time[0];
  int    *bodiesIds       = &tree.bodies_ids[0];

  real4  tmpp; 
  float2 tmpp2;
  int    tmpp3;
  int *firstloc   = new int[nProcs+1];
  int *nparticles = new int[nProcs+1];

  // Loop over particles and determine which particle needs to go where
  // reorder the bodies in such a way that bodies that have to be send
  // away are stored after each other in the array
  double t1 = get_time();

  //Array reserve some memory at forehand , 1%
  vector<bodyStruct> array2Send;
  //vector<bodyStruct> array2Send(((int)(tree.n * 0.01)));

  for(int ib=0;ib<nproc;ib++)
  {
    int ibox       = (ib+myid)%nproc;
    firstloc[ibox] = iloc;      //Index of the first particle send to proc: ibox

    for(int i=iloc; i<nbody;i++)
    {
      //      if(myid == 0){PRC(i); PRC(pb[i].get_pos());}
      if(isinbox(bodiesPositions[i], xlow[ibox], xhigh[ibox]))
      {
        //Position
        tmpp                  = bodiesPositions[iloc];
        bodiesPositions[iloc] = bodiesPositions[i];
        bodiesPositions[i]    = tmpp;
        //Velocity
        tmpp             = velocities[iloc];
        velocities[iloc] = velocities[i];
        velocities[i]    = tmpp;
        //Acc0
        tmpp             = bodiesAcc0[iloc];
        bodiesAcc0[iloc] = bodiesAcc0[i];
        bodiesAcc0[i]    = tmpp;
        //Acc1
        tmpp             = bodiesAcc1[iloc];
        bodiesAcc1[iloc] = bodiesAcc1[i];
        bodiesAcc1[i]    = tmpp;    
        //Acc1
        tmpp2            = bodiesTime[iloc];
        bodiesTime[iloc] = bodiesTime[i];
        bodiesTime[i]    = tmpp2;              
        //IDs
        tmpp3            = bodiesIds[iloc];
        bodiesIds[iloc]  = bodiesIds[i];
        bodiesIds[i]     = tmpp3;   

        //Put the particle in the array of to send particles
        if(ibox != myid)
        {
          bodyStruct body;
          body.pos  = bodiesPositions[iloc];
          body.vel  = velocities[iloc];
          body.acc0 = bodiesAcc0[iloc];
          body.acc1 = bodiesAcc1[iloc];
          body.time = bodiesTime[iloc];
          body.id   = bodiesIds[iloc];
        
          array2Send.push_back(body);
        }

        iloc++;
      }// end if
    }//for i=iloc
    nparticles[ibox] = iloc-firstloc[ibox];//Number of particles that has to be send to proc: ibox
  } // for(int ib=0;ib<nproc;ib++)
  
  printf("Benodigde zoek tijd: %lg ,proc: %d gevonden in eigen box: %d n: %d  en naar andere: %ld \n", 
         get_time()-t1, myid, nparticles[myid], tree.n, array2Send.size());
  
/*  if(procId == 0)
  {
    for(int i=0; i < 4; i++)
    {
      fprintf(stderr, "%d ->  %d = %d en in bron: %d = %d \n",  
              procId, i, array2Send[i].id, i, bodiesIds[i]);
    }
  }
  mpiSync();
  if(procId == 1)
  {
    for(int i=0; i < 4; i++)
    {
      fprintf(stderr, "%d ->  %d = %d en in bron: %d = %d \n",  
              procId, i, array2Send[i].id, i, bodiesIds[i]);
    }
  }  */

  
  t1 = get_time();

  totalsent = nbody - nparticles[myid];

  int tmp;
  MPI_Reduce(&totalsent,&tmp,1, MPI_INT, MPI_SUM,0,MPI_COMM_WORLD);
  if(procId == 0)
  {
    totalsent = tmp;
    cout << "Exchanged particles = " << totalsent << endl;
  }
  
  if(iloc < nbody)
  {
      cerr << procId <<" exchange_particle error: particle in no box...iloc: " << iloc 
                     << " and nbody: " << nbody << "\n";
  }
    
  vector<bodyStruct> recv_buffer3(nbody- nparticles[myid]);
  
//   cout << myid << "ik ontvang: " << nbody- nparticles[myid] << "\t size"  << recv_buffer3.size() << endl;
  
  int tempidFirst, tempRecvCount; 
  unsigned int recvCount = 0;
  
  //Exchange the data with the other processors
  int ibend = -1;
  int nsend;
  int isource = 0;
  for(int ib=nproc-1;ib>0;ib--)
  {
    int ibox = (ib+myid)%nproc; //index to send...
      
    if (ib == nproc-1)
    {
      isource= (myid+1)%nproc;
    }
    else
    {
      isource = (isource+1)%nproc;
      if (isource == myid)isource = (isource+1)%nproc;
    }

    if(MP_exchange_particle_with_overflow_check<bodyStruct>(ibox, &array2Send[0],
                                                    recv_buffer3, firstloc[ibox] - nparticles[myid],
                                                    nparticles[ibox], isource, 
                                                    nsend, recvCount))
    {
      ibend = ibox; //Here we get if exchange failed
      ib = 0;
    }//end if mp exchang
  }//end for all boxes
  
  printf("Benodigde inter-process communicatie tijd: %lg ,proc: %d\n", 
         get_time()-t1, myid);  
  t1 = get_time();       
  double t2= t1;
    
  //    ... should do something different for nsend...                                           
  int idfirst;
  if(ibend >= 0)
  {
      idfirst = firstloc[ibend]+nparticles[ibend]-nsend;
  }
  else
  {
      idfirst = nparticles[myid];
  }

  //Have to resize the bodies vector to keep the numbering correct
  //cerr << "idfirst: " << idfirst << " nbmax " << nbmax << " iloc: " << iloc << endl;      
  tree.setN(idfirst+recvCount);    
  tree.bodies_pos.cresize (idfirst+recvCount + 1, false);  
  tree.bodies_acc0.cresize(idfirst+recvCount,     false);
  tree.bodies_acc1.cresize(idfirst+recvCount,     false);
  tree.bodies_vel.cresize (idfirst+recvCount,     false);
  tree.bodies_time.cresize(idfirst+recvCount,     false);
  tree.bodies_ids.cresize (idfirst+recvCount + 1, false);
  tree.bodies_Ppos.cresize(idfirst+recvCount + 1, false);  
  
  tree.generalBuffer1.cresize(3*(idfirst+recvCount)*4, false);
  
  printf("Benodigde gpu malloc tijd stap 1: %lg \t Size: %d \tRank: %d \t Size: %d \n", 
         get_time()-t1, idfirst+recvCount, mpiGetRank(), tree.bodies_Ppos.get_size()); 
  t1 = get_time();
//   cerr << myid << "ik ontvang: " << nbody- nparticles[myid] << "\t size"  << recv_buffer3.size() << "\t" << recvCount << endl;
    
  tempidFirst = idfirst; tempRecvCount = recvCount;

  //Copy data from struct into the main arrays
  for(unsigned int P=0; P < recvCount; P++)
  {
    //  cerr << P << "\t" << recv_buffer3[P].pos.x << "\t" << recv_buffer3[P].vel.x << "\t" << recv_buffer3[P].acc1.x;
    // cerr << "\t" << recv_buffer3[P].id << "\t" << recv_buffer3[P].acc0.z << endl;       
    tree.bodies_pos[idfirst+P]  = recv_buffer3[P].pos;        tree.bodies_vel[idfirst+P]      = recv_buffer3[P].vel;
    tree.bodies_acc0[idfirst+P] = recv_buffer3[P].acc0;       tree.bodies_acc1[idfirst+P]     = recv_buffer3[P].acc1;
    tree.bodies_time[idfirst+P] = recv_buffer3[P].time;       tree.bodies_ids[idfirst+P]      = recv_buffer3[P].id;

//     tree.bodies_Ppos[idfirst+P]  = recv_buffer3[P].pos; 
  }

  printf("Benodigde DATA in struct copy time: %lg \n", get_time()-t1); t1 = get_time();

 
  if(ibend == -1){
    
  }else{
      //Something went wrong
    cerr << "ERROR in exchange_particles_with_overflow_check! \n"; exit(0);
  }

  
  //Resize the arrays of the tree    
  reallocateParticleMemory(tree);   
     
  printf("Benodigde gpu malloc tijd stap 2: %lg \n", get_time()-t1);
  printf("Totale GPU interactie tijd: %lg \n", get_time()-t2);

  int retValue = 0;


  delete[] firstloc;
  delete[] nparticles;

  return retValue;
}


//Local essential tree functions


//TODO change this into one mpi all gather
void octree::getAllLETBoxes(real4 bMin, real4 bMax)
{
  //Gathers the box sizes of all processors and sends them to all processors
  //we use this since we predict particles after we set the domain sizes, so 
  //predict can move them over the domain size

  double4 tempMin;
  tempMin.x = bMin.x; tempMin.y = bMin.y; tempMin.z = bMin.z; tempMin.w = bMin.w;

  double4 tempMax;
  tempMax.x = bMax.x; tempMax.y = bMax.y; tempMax.z = bMax.z; tempMax.w = bMax.w;

  MPI_Allgather(&tempMin, sizeof(double4), MPI_BYTE,  cur_xlow, sizeof(double4), MPI_BYTE, MPI_COMM_WORLD);
  MPI_Allgather(&tempMax, sizeof(double4), MPI_BYTE,  cur_xhigh, sizeof(double4), MPI_BYTE, MPI_COMM_WORLD);
}


void octree::essential_tree_exchange(vector<real4> &treeStructure, tree_structure &tree, tree_structure &remote)
{
  int myid = procId;
  int nproc = nProcs;
  int isource = 0;

  real4  *bodies              = &tree.bodies_Ppos[0];
  real4  *velocities          = &tree.bodies_Pvel[0];
  real4  *multipole           = &tree.multipole[0];
  real4  *nodeSizeInfo        = &tree.boxSizeInfo[0];
  real4  *nodeCenterInfo      = &tree.boxCenterInfo[0];
  
  vector<real4> recv_particles;
  vector<real4> recv_multipoleData;
  vector<real4> recv_nodeSizeData;
  vector<real4> recv_nodeCenterData;  
  
  real4 **treeBuffers;
  
  //   treeBuffers  = new real4*[mpiGetNProcs()-1]; //creates a new array of pointers to int objects
  //Own tree add test
  treeBuffers  = new real4*[mpiGetNProcs()]; //creates a new array of pointers to int objects
  int recvTree = 0;
 
  //For each process
  for(int ib=nproc-1;ib>0;ib--)
  {

    int ibox = (ib+myid)%nproc; //index to send...
    if (ib == nproc-1){
        isource= (myid+1)%nproc;
    }else{
        isource = (isource+1)%nproc;
        if (isource == myid)isource = (isource+1)%nproc;
    }

  /*
    cerr << "\nibox: " << ibox << endl;
    cerr << "Other proc has box: low: " << let_xlow[ibox].x << "\t" <<  let_xlow[ibox].y  << "\t" <<  let_xlow[ibox].z 
                                    << "\thigh: "  << let_xhigh[ibox].x << "\t" <<  let_xhigh[ibox].y  << "\t" <<  let_xhigh[ibox].z << endl;*/

    double4 boxCenter = {     0.5*(cur_xlow[ibox].x  + cur_xhigh[ibox].x),
                              0.5*(cur_xlow[ibox].y  + cur_xhigh[ibox].y),
                              0.5*(cur_xlow[ibox].z  + cur_xhigh[ibox].z), 0};
    double4 boxSize   = {fabs(0.5*(cur_xhigh[ibox].x - cur_xlow[ibox].x)),
                         fabs(0.5*(cur_xhigh[ibox].y - cur_xlow[ibox].y)),
                         fabs(0.5*(cur_xhigh[ibox].z - cur_xlow[ibox].z)), 0};  
                          
                          
//   printf("Other proc center and size: [%f %f %f] \t [%f %f %f] \n", boxCenter.x, boxCenter.y, boxCenter.z,
//          boxSize.x, boxSize.y, boxSize.z);

    uint2 node_begend;
    int level_start = 2;
    node_begend.x = tree.level_list[level_start].x;
    node_begend.y = tree.level_list[level_start].y;
    
    int particleCount, nodeCount;
    
//     double t1 = get_time();    
    create_local_essential_tree_count(bodies, multipole, nodeSizeInfo, nodeCenterInfo,
                                boxCenter, boxSize, cur_xlow[ibox].w, node_begend.x, node_begend.y,
                                particleCount, nodeCount);    

    //Buffer that will contain all the data:
    //|real4| 2*particleCount*real4| nodes*real4 | nodes*real4 | nodes*3*real4 |
    //1 + 2*particleCount + nodeCount + nodeCount + 3*nodeCount
    
    //0-1 )                               Info about #particles, #nodes, start and end of tree-walk
    //1- Npart)                           The particle positions
    //1+Npart-Npart )                     The particle velocities
    //1+2*Npart-Nnode )                   The nodeSizeData
    //1+*2Npart+Nnode - Npart+2*Nnode )   The nodeCenterData
    //1+2*Npart+2*Nnode - Npart+5*Nnode ) The multipole data, is 3x number of nodes (mono and quadrupole data)
    int bufferSize = 1 + 2*particleCount + 5*nodeCount;
    real4 *letDataBuffer = new real4[bufferSize];
         
    create_local_essential_tree_fill(bodies, velocities, multipole, nodeSizeInfo, nodeCenterInfo,
                                boxCenter, boxSize, cur_xlow[ibox].w, node_begend.x, node_begend.y,
                                particleCount, nodeCount, letDataBuffer);        
/*    
    printf("Speciaal: %lg\n", get_time()-t1);   
                                t1 = get_time();    
    create_local_essential_tree(bodies, multipole, nodeSizeInfo, nodeCenterInfo,
                                boxCenter, boxSize, node_begend.x, node_begend.y,
                                particles, multipoleData, nodeSizeData, nodeCenterData);
    printf("Gewoon: %lg\n", get_time()-t1); */                                
                                
    //Set the tree properties, before we exchange the data
    letDataBuffer[0].x = particleCount;         //Number of particles in the LET
    letDataBuffer[0].y = nodeCount;             //Number of nodes     in the LET
    letDataBuffer[0].z = node_begend.x;         //First node on the level that indicates the start of the tree walk
    letDataBuffer[0].w = node_begend.y;         //last node on the level that indicates the start of the tree walk

    //Exchange the data of the tree structures in between the nodes
    treeBuffers[recvTree] = MP_exchange_bhlist(ibox, isource, 
                                               bufferSize, letDataBuffer);
    recvTree++; //Increase the index by 1                                                

    delete[] letDataBuffer;       
  }//end for each process
  
  //We now have nProcs-1 partial local essential trees, we combine them into one tree
  //TODO As extra test we could also combine them with the main tree and see 
  //what the exectuion time would be. If we would do 1 tree-traverse instead of 2.
  //But lets first combine the results of one tree
  
  //Add the processors own tree to the LET tree

  int level_start   = 2;
  int particleCount = tree.n;
  int nodeCount = tree.n_nodes;
  int bufferSizeLocal = 1 + 2*particleCount + 5*nodeCount;
  
  //Only alloc if we have NPROC set to the total number of processes
  int PROCS = mpiGetNProcs();
 PROCS -= 1; //DO -1 when NOT including the own tree
  
  if(PROCS == mpiGetNProcs())
  {
    treeBuffers[mpiGetNProcs()-1] = new real4[bufferSizeLocal];

    int idx = 1;
    memcpy(&treeBuffers[mpiGetNProcs()-1][idx], &bodies[0], sizeof(real4)*particleCount);
    idx += particleCount;
    memcpy(&treeBuffers[mpiGetNProcs()-1][idx], &velocities[0], sizeof(real4)*particleCount);
    idx += particleCount;
    memcpy(&treeBuffers[mpiGetNProcs()-1][idx], &nodeSizeInfo[0], sizeof(real4)*nodeCount);
    idx += nodeCount;
    memcpy(&treeBuffers[mpiGetNProcs()-1][idx], &nodeCenterInfo[0], sizeof(real4)*nodeCount); 
    idx += nodeCount;
    memcpy(&treeBuffers[mpiGetNProcs()-1][idx], &multipole[0], sizeof(real4)*nodeCount*3);   
    
    treeBuffers[mpiGetNProcs()-1][0].x = particleCount;
    treeBuffers[mpiGetNProcs()-1][0].y = nodeCount;
    treeBuffers[mpiGetNProcs()-1][0].z = tree.level_list[level_start].x;
    treeBuffers[mpiGetNProcs()-1][0].w = tree.level_list[level_start].y;  
  }
  
  //Arrays to store and compute the offsets
  int *particleSumOffsets  = new int[mpiGetNProcs()+1];
  int *nodeSumOffsets      = new int[mpiGetNProcs()+1];
  int *startNodeSumOffsets = new int[mpiGetNProcs()+1];    
  uint2 *nodesBegEnd       = new uint2[mpiGetNProcs()+1];
  
  //Offsets start at 0 and then are increased by the number of nodes of each LET tree
  particleSumOffsets[0]           = 0;
  nodeSumOffsets[0]               = 0;
  startNodeSumOffsets[0]          = 0;  
//   nodesBegEnd[mpiGetNProcs()-1].x = nodesBegEnd[mpiGetNProcs()-1].y = 0; //Make valgrind happy
  nodesBegEnd[mpiGetNProcs()].x = nodesBegEnd[mpiGetNProcs()].y = 0; //Make valgrind happy
  int totalTopNodes               = 0;
  

  //Calculate the offsets
//   for(int i=0; i < mpiGetNProcs() -1; i++)
  for(int i=0; i < PROCS ; i++) //adding own tree test
  {
    int particles = (int)treeBuffers[i][0].x;
    int nodes     = (int)treeBuffers[i][0].y;
   
    nodesBegEnd[i].x = (int)treeBuffers[i][0].z;
    nodesBegEnd[i].y = (int)treeBuffers[i][0].w;
    
    totalTopNodes += nodesBegEnd[i].y-nodesBegEnd[i].x;
    
    particleSumOffsets[i+1]     = particleSumOffsets[i]  + particles;
    nodeSumOffsets[i+1]         = nodeSumOffsets[i]      + nodes - nodesBegEnd[i].y;    //Without the top-nodes
    startNodeSumOffsets[i+1]    = startNodeSumOffsets[i] + nodesBegEnd[i].y-nodesBegEnd[i].x;
  }

  //Compute total particles and total nodes, totalNodes is WITHOUT topNodes
//   int totalParticles    = particleSumOffsets[mpiGetNProcs()-1];
//   int totalNodes        = nodeSumOffsets[mpiGetNProcs()-1];
  int totalParticles    = particleSumOffsets[PROCS];
  int totalNodes        = nodeSumOffsets[PROCS];
  
  //To bind parts of the memory to different textures, the memory start address
  //has to be aligned with XXX bytes, so totalParticles*sizeof(real4) has to be
  //increased by an offset, so that the node data starts at a XXX byte boundary
  //same with the node information. XXX is a architecture specific value
  
  const int texBoundary = 512; //Fermi
  
  int particleTextOffset = 0;
  //Compute the number of bytes  
  particleTextOffset = totalParticles*sizeof(real4); 
  //Compute number of 256 byte blocks  
  particleTextOffset = (particleTextOffset / texBoundary) + (((particleTextOffset % texBoundary) > 0) ? 1 : 0); 
  //Compute the number of bytes padded / offset 
  particleTextOffset = (particleTextOffset * texBoundary) - totalParticles*sizeof(real4); 
  //Back to the actual number of elements
  particleTextOffset = particleTextOffset / sizeof(real4); 
  
  //Same steps for the node data
  int nodeTextOffset = 0;
  nodeTextOffset = (totalNodes+totalTopNodes)*sizeof(real4); 
  nodeTextOffset = (nodeTextOffset / texBoundary) + (((nodeTextOffset % texBoundary) > 0) ? 1 : 0); 
  nodeTextOffset = (nodeTextOffset * texBoundary) - (totalNodes+totalTopNodes)*sizeof(real4); 
  nodeTextOffset = nodeTextOffset / sizeof(real4); 

//   particleTextOffset = 0;
//   nodeTextOffset = 0;

  //Compute the total size of the buffer
  //int bufferSize      = totalParticles + 5*(totalNodes+totalTopNodes);
  int bufferSize        = 2*(totalParticles+particleTextOffset) + 5*(totalNodes+totalTopNodes + nodeTextOffset);
  
  
  //Allocate memory on host and device to store the merged tree-structure
  remote.fullRemoteTest.cresize(bufferSize, false);  //Change the size but ONLY if we need more memory
    
  real4 *combinedRemoteTree = &remote.fullRemoteTest[0];

  fprintf(stderr, "Total (%d): %d \t %d topNodes: %d BufferSize: %d\n", mpiGetRank(), totalParticles, totalNodes, totalTopNodes, bufferSize);
  
  //Copy all the pieces of the different trees at the correct memory offsets
//   for(int i=0; i < mpiGetNProcs() -1; i++)
  for(int i=0; i < PROCS; i++) //Adding own tree test
  {
    //Get the properties of the LET
    int remoteP = (int) treeBuffers[i][0].x;    //Number of particles
    int remoteN = (int) treeBuffers[i][0].y;    //Number of nodes
    int remoteB = (int) treeBuffers[i][0].z;    //Begin id of top nodes
    int remoteE = (int) treeBuffers[i][0].w;    //End   id of top nodes
    int remoteNstart = remoteE-remoteB;

    //Particles
    memcpy(&combinedRemoteTree[particleSumOffsets[i]],   &treeBuffers[i][1], sizeof(real4)*remoteP);
    
    //Velocities
    memcpy(&combinedRemoteTree[(particleTextOffset + totalParticles) + particleSumOffsets[i]],   
           &treeBuffers[i][1+remoteP], sizeof(real4)*remoteP);
 
    //The start nodes, nodeSizeInfo
    memcpy(&combinedRemoteTree[2*(particleTextOffset + totalParticles) + startNodeSumOffsets[i]],  
           &treeBuffers[i][1+2*remoteP+remoteB], //From the start node onwards
           sizeof(real4)*remoteNstart);
    
    //Non start nodes, nodeSizeInfo       
    memcpy(&combinedRemoteTree[2*(particleTextOffset + totalParticles) +  totalTopNodes + nodeSumOffsets[i]],  
           &treeBuffers[i][1+2*remoteP+remoteE], //From the last start node onwards
           sizeof(real4)*(remoteN-remoteE));    
   
    //The start nodes, nodeCenterInfo
    memcpy(&combinedRemoteTree[2*(particleTextOffset + totalParticles) + startNodeSumOffsets[i]
                                + (totalNodes + totalTopNodes + nodeTextOffset)],  
           &treeBuffers[i][1+2*remoteP+remoteB + remoteN], //From the start node onwards
           sizeof(real4)*remoteNstart);
    
    //Non start nodes, nodeCenterInfo       
    memcpy(&combinedRemoteTree[2*(particleTextOffset + totalParticles) +  totalTopNodes 
           + nodeSumOffsets[i] + (totalNodes + totalTopNodes + nodeTextOffset)],  
           &treeBuffers[i][1+2*remoteP+remoteE + remoteN], //From the last start node onwards
           sizeof(real4)*(remoteN-remoteE));   
        
    //The start nodes, multipole       
    memcpy(&combinedRemoteTree[2*(particleTextOffset + totalParticles) + 3*startNodeSumOffsets[i] +
           2*(totalNodes+totalTopNodes + nodeTextOffset)],  
           &treeBuffers[i][1+2*remoteP+2*remoteN + 3*remoteB], //From the start node onwards
           sizeof(real4)*remoteNstart*3);  
                     
    //Non start nodes, multipole       
    memcpy(&combinedRemoteTree[2*(particleTextOffset + totalParticles) +  3*totalTopNodes + 
           3*nodeSumOffsets[i] + 2*(totalNodes+totalTopNodes+nodeTextOffset)],  
           &treeBuffers[i][1+2*remoteP+remoteE*3 + 2*remoteN], //From the last start node onwards
           sizeof(real4)*(remoteN-remoteE)*3);               
    /*
    |real4| 2*particleCount*real4| nodes*real4 | nodes*real4 | nodes*3*real4 |
    1 + 2*particleCount + nodeCount + nodeCount + 3*nodeCount
    
    Info about #particles, #nodes, start and end of tree-walk
    The particle positions
    velocities
    The nodeSizeData
    The nodeCenterData
    The multipole data, is 3x number of nodes (mono and quadrupole data)  
    
    Now that the data is copied, modify the offsets of the tree so that everything works
    with the new correct locations and references. This takesp lace in two steps:
    First  the top nodes since
    Second the normal nodes
    Has to be done in two steps since they are not continous in memory if NPROCS > 2
      */  

    //Modify the top nodes
    int modStart = 2*(particleTextOffset + totalParticles) + startNodeSumOffsets[i];
    int modEnd   = modStart       + remoteNstart;

    for(int j=modStart; j < modEnd; j++)
    {
      real4 nodeCenter = combinedRemoteTree[j+totalTopNodes+totalNodes+nodeTextOffset];
      real4 nodeSize   = combinedRemoteTree[j];
      bool leaf        = nodeCenter.w <= 0;

      int childinfo = host_float_as_int(nodeSize.w);          
      int child, nchild;

      if(!leaf)
      {
        //Node
        child    =    childinfo & 0x0FFFFFFF;                  //Index to the first child of the node
        nchild   = (((childinfo & 0xF0000000) >> 28)) ;        //The number of children this node has          
        
        child = child - nodesBegEnd[i].y + totalTopNodes + nodeSumOffsets[i]; //Calculate the new start (non-leaf)        
        child = child | (nchild << 28);                                       //Merging back in one int
        
        if(nchild == 0) child = 0;                             //To prevent incorrect negative values
      }else{ //Leaf        
        child   =   childinfo & BODYMASK;                      //the first body in the leaf
        nchild  = (((childinfo & INVBMASK) >> LEAFBIT)+1);     //number of bodies in the leaf masked with the flag    
        
        child   =  child + particleSumOffsets[i];               //Increasing offset
        child   = child | ((nchild-1) << LEAFBIT);              //Merging back to one int
        
        if(nchild == 0) child = 0;                              //To prevent incorrect negative values
      }//end !leaf      
      combinedRemoteTree[j].w =  host_int_as_float(child);      //store the modified offset
    } 
    
    //Now the non-top nodes for this process
    modStart =  totalTopNodes + nodeSumOffsets[i] + 2*(particleTextOffset + totalParticles);
    modEnd   =  modStart      + remoteN-remoteE;
    for(int j=modStart; j < modEnd; j++)
    {
      real4 nodeCenter = combinedRemoteTree[j+totalTopNodes+totalNodes+nodeTextOffset];
      real4 nodeSize   = combinedRemoteTree[j];
      bool leaf        = nodeCenter.w <= 0;
          
      int childinfo = host_float_as_int(nodeSize.w);          
      int child, nchild;

      if(!leaf) {  //Node       
        child    =    childinfo & 0x0FFFFFFF;                   //Index to the first child of the node
        nchild   = (((childinfo & 0xF0000000) >> 28)) ;         //The number of children this node has              

        //Calculate the new start (non-leaf)
        child = child - nodesBegEnd[i].y + totalTopNodes + nodeSumOffsets[i];  ;      
        
        //Combine and store
        child = child | (nchild << 28);  
        
        //if(nchild == 0xFFFFFFFF) child = 0;
        if(nchild == 0) child = 0;                              //To prevent incorrect negative values
      }else{ //Leaf
        child   =   childinfo & BODYMASK;                       //the first body in the leaf
        nchild  = (((childinfo & INVBMASK) >> LEAFBIT)+1);      //number of bodies in the leaf masked with the flag    
        
        child =  child + particleSumOffsets[i];                 //Modify the particle offsets
        child = child | ((nchild-1) << LEAFBIT);                //Merging the data back into one int
        
        //if(childinfo == 0xFFFFFFFF) child = 0;
        if(nchild == 0) child = 0;                              //To prevent incorrect negative values
      }//end !leaf
      combinedRemoteTree[j].w =  host_int_as_float(child);      //Store the modified value
    }   
    
     delete[] treeBuffers[i];    //Free the memory of this part of the LET
  }
  
  /*  
  The final tree structure looks as follows:
  particlesT1, partcilesT2,...mparticlesTn |,
  topNodeSizeT1, topNodeSizeT2,..., topNodeSizeT2 | nodeSizeT1, nodeSizeT2, ...nodeSizeT3 |,
  topNodeCentT1, topNodeCentT2,..., topNodeCentT2 | nodeCentT1, nodeCentT2, ...nodeCentT3 |, 
  topNodeMultT1, topNodeMultT2,..., topNodeMultT2 | nodeMultT1, nodeMultT2, ...nodeMultT3 
  
  NOTE that the Multipole data consists of 3 float4 values per node
    
  */
  
  //Combine the offsets into one int32 so we can extract it in the other
  //part of the program
  particleTextOffset = (particleTextOffset  << 16) | nodeTextOffset;
  
  //Store the tree properties (number of particles, number of nodes, start and end topnode)
  remote.remoteTreeStruct.x = totalParticles;
  remote.remoteTreeStruct.y = totalNodes+totalTopNodes;
  remote.remoteTreeStruct.z = particleTextOffset;
  remote.remoteTreeStruct.w = totalTopNodes;
    
  cerr << "Aantal local bodies: " << tree.n << " aantal LET bodies: " << totalParticles << endl;
  
  delete[] particleSumOffsets;
  delete[] nodeSumOffsets;
  delete[] startNodeSumOffsets;
  delete[] nodesBegEnd;
  delete[] treeBuffers;
}


//Improved Barnes Hut criterium
#ifdef INDSOFT
bool split_node_grav_impbh(float4 nodeCOM, double4 boxCenter, double4 boxSize,
                           float group_eps, float node_eps)                                      
#else
bool split_node_grav_impbh(float4 nodeCOM, double4 boxCenter, double4 boxSize)
#endif
{
  //Compute the distance between the group and the cell
  float3 dr = {fabsf(boxCenter.x - nodeCOM.x) - (boxSize.x),
               fabsf(boxCenter.y - nodeCOM.y) - (boxSize.y),
               fabsf(boxCenter.z - nodeCOM.z) - (boxSize.z)};

  dr.x += fabsf(dr.x); dr.x *= 0.5f;
  dr.y += fabsf(dr.y); dr.y *= 0.5f;
  dr.z += fabsf(dr.z); dr.z *= 0.5f;

  //Distance squared, no need to do sqrt since opening criteria has been squared
  float ds2    = dr.x*dr.x + dr.y*dr.y + dr.z*dr.z;

  #ifdef INDSOFT
    //Naar idee van Inti nu minder overbodige openingen                                                                                                       
    if(ds2 <= ((group_eps + node_eps ) * (group_eps + node_eps) )) return true;              
  #endif

  return (ds2 <= fabs(nodeCOM.w));
}

//Minimal Distance version

//Minimum distance opening criteria
#ifdef INDSOFT
  bool split_node(real4 nodeCenter, real4 nodeSize, double4 boxCenter, double4 boxSize,
                                     float group_eps, float node_eps)
#else
  bool split_node(real4 nodeCenter, real4 nodeSize, double4 boxCenter, double4 boxSize)
#endif
{
  //Compute the distance between the group and the cell
  float3 dr = {fabsf(boxCenter.x - nodeCenter.x) - (boxSize.x + nodeSize.x),
               fabsf(boxCenter.y - nodeCenter.y) - (boxSize.y + nodeSize.y),
               fabsf(boxCenter.z - nodeCenter.z) - (boxSize.z + nodeSize.z)};

  dr.x += fabsf(dr.x); dr.x *= 0.5f;
  dr.y += fabsf(dr.y); dr.y *= 0.5f;
  dr.z += fabsf(dr.z); dr.z *= 0.5f;

  
  //Do the boxes overlap? TODO check if this is correct method
//   if(dr.x == 0) return true;
//   if(dr.y == 0) return true;
//   if(dr.z == 0) return true;
  

  float ds2    = dr.x*dr.x + dr.y*dr.y + dr.z*dr.z;
  
  #ifdef INDSOFT
    //Naar idee van Inti nu minder overbodige openingen                                                                                                       
    if(ds2 <= ((group_eps + node_eps ) * (group_eps + node_eps) )) return true;              
  #endif  
  
//   if(dr.x == 0 || dr.y == 0 || dr.z == 0)
//   {
//     cerr << "Test overlap! \n";
//     fprintf(stderr, "dr: [%f, %f, %f] [ds: %f crit: %f open: %f ]\n", dr.x,dr.y,dr.z, ds2, nodeCenter.w, ds2 < fabs(nodeCenter.w));
//     fprintf(stderr, "Node: [%f, %f, %f] en [%f, %f, %f] \n", nodeCenter.x,nodeCenter.y,nodeCenter.z,
//                                                           nodeSize.x, nodeSize.y, nodeSize.z);
//     fprintf(stderr, "Group: [%f, %f, %f] en [%f, %f, %f] \n", boxCenter.x,boxCenter.y,boxCenter.z,
//                                                           boxSize.x, boxSize.y, boxSize.z);    
//     
//     
//     
//     exit(0);
//   }  

  //Cell opening magic
  return (ds2 < fabs(nodeCenter.w));
 }

void octree::create_local_essential_tree(real4* bodies, real4* multipole, real4* nodeSizeInfo, real4* nodeCenterInfo,   
                                         double4 boxCenter, double4 boxSize, float group_eps, int start, int end, 
                                         vector<real4> &particles, vector<real4> &multipoleData,
                                         vector<real4> &nodeSizeData, vector<real4> &nodeCenterData)
{
    //Walk the tree as is done on device, level by level
    vector<int> curLevel;
    vector<int> nextLevel;
    
 
    double t1 = get_time();

    double massSum = 0;

    int nodeCount       = 0;
    
    //Add the initial nodes to the curLevel list
    for(int i=start; i < end; i++)
    {
      curLevel.push_back(i);
    }
    
    //Add the nodes before the start and end to the node list
    for(int i=0; i < start; i++)
    {
      nodeSizeData.push_back(nodeSizeInfo[i]);
      nodeCenterData.push_back(nodeCenterInfo[i]);   
      
      multipoleData.push_back(multipole[i*3 + 0]);
      multipoleData.push_back(multipole[i*3 + 1]);
      multipoleData.push_back(multipole[i*3 + 2]);
      nodeCount++;
    }
      
    //Start the tree-walk
    printf("Start: %d end: %d \n", start, end);
    cout << "Starting walk on: " << curLevel.size() << " items! \n"; 
    
    
    int childNodeOffset         = end;
    int childParticleOffset     = 0;
    
    while(curLevel.size() > 0)
    {
      for(unsigned int i=0; i < curLevel.size(); i++)
      {
        //Read node data
        int node         = curLevel[i];
        real4 nodeCenter = nodeCenterInfo[node];
        real4 nodeSize   = nodeSizeInfo[node];
        bool leaf        = nodeCenter.w <= 0;
        
        union{float f; int i;} u; //__float_as_int
        u.f           = nodeSize.w;      
        int childinfo = u.i;
        
        int child, nchild;
        if(!leaf)
        {
          //Node
          child    =    childinfo & 0x0FFFFFFF;                         //Index to the first child of the node
          nchild   = (((childinfo & 0xF0000000) >> 28)) ;         //The number of children this node has              
        }
        else
        {
          //Leaf
          child   =   childinfo & BODYMASK;                                     //thre first body in the leaf
          nchild  = (((childinfo & INVBMASK) >> LEAFBIT)+1);     //number of bodies in the leaf masked with the flag        
        }
        
        #ifdef INDSOFT
          //Very inefficient this but for testing I have to live with it...                                                                                   
          float node_eps_val = multipole[node*3 + 1].w; 
        #endif         
        
        #ifdef IMPBH
          //Improved barnes hut version
          float4 nodeCOM     = multipole[node*3 + 0];
          nodeCOM.w = nodeCenter.w;
          
          #ifdef INDSOFT
             bool split = split_node_grav_impbh(nodeCOM, boxCenter, boxSize, group_eps, node_eps_val);
          #else
             bool split = split_node_grav_impbh(nodeCOM, boxCenter, boxSize);  
          #endif
          
        #else
          //Minimal distance version
          
          #ifdef INDSOFT
            bool split = split_node(nodeCenter, nodeSize, boxCenter, boxSize, group_eps, node_eps_val);  //Check if node should be split
          #else
            bool split = split_node(nodeCenter, nodeSize, boxCenter, boxSize);
          #endif          
          
        #endif 
  //       printf("Node %d is a leaf: %d  en childinfo: %d  \t\t-> %d en %d \t split: %d\n", node, leaf, childinfo, child, nchild, split);
//          split = false;
        uint temp =0;
        //if split & node add children to next lvl stack
        if(split && !leaf)
        {       
          for(int i=child; i < child+nchild; i++)
          {
            nextLevel.push_back(i);            
          }
          
          temp = childNodeOffset | (nchild << 28);  
          //Update reference to children
          childNodeOffset += nchild;
        }
        
        //if split & leaf add particles to particle list
        if(split && leaf)
        { 
          for(int i=child; i < child+nchild; i++)
          {
            particles.push_back(bodies[i]);
            massSum += bodies[i].w;
          }
          
          temp = childParticleOffset | ((nchild-1) << LEAFBIT);
          childParticleOffset += nchild;          
        }
        
        
        //Add the node data to the appropriate arrays
        //and modify the node reference
        //start ofset for its children, should be nodeCount at start of this level +numberofnodes on this level
        //plus a counter that counts the number of childs of the nodes we have tested

        //New childoffset:
        union{int i; float f;} itof; //__int_as_float
        itof.i           = temp;            
        float tempConv = itof.f;
        
        //Add node properties and update references
        real4 nodeSizeInfoTemp  = nodeSizeInfo[node];
        nodeSizeInfoTemp.w      = tempConv;             //Replace child reference
        nodeSizeData.push_back(nodeSizeInfoTemp);
      
        
        multipoleData.push_back(multipole[node*3 + 0]);
        multipoleData.push_back(multipole[node*3 + 1]);
        multipoleData.push_back(multipole[node*3 + 2]);        
        
        if(!split)
        {          
          massSum += multipole[node*3 + 0].w;
/*          
          if(mpiGetRank() == 0)
            cout << "Mass van " << node << " is: " << multipole[node*3 + 0].w << endl;
          */
        }


      } //end for curLevel.size
      
      
      //Put next level stack into current level and continue
      curLevel.clear();
      
//       cout << "Next level: " << nextLevel.size() << endl;
      curLevel.assign(nextLevel.begin(), nextLevel.end());
      nextLevel.clear();
      
      
    }//end while
  
  cout << "New tree structure: bodies: " << particles.size() << "\tnodes: " << nodeSizeData.size() << "\t took: " << get_time() -t1 << endl;
  cout << "Mass sum: " << massSum << endl;
  cout << "Mass sumtest: " << multipole[0*0 + 0].w << endl;
  
}


void octree::create_local_essential_tree_fill(real4* bodies, real4* velocities, real4* multipole, real4* nodeSizeInfo, real4* nodeCenterInfo,   
                                         double4 boxCenter, double4 boxSize, float group_eps, int start, int end, 
                                         int particleCount, int nodeCount, real4 *dataBuffer)
{
    //Walk the tree as is done on device, level by level
    vector<int> curLevel;
    vector<int> nextLevel;
  
    
//     double t1 = get_time();

    double massSum = 0;

    int particleOffset   = 1;
    int velParticleOffset  = particleOffset      + particleCount; 
    int nodeSizeOffset     = velParticleOffset   + particleCount;
    int nodeCenterOffset   = nodeSizeOffset      + nodeCount;
    int multiPoleOffset    = nodeCenterOffset    + nodeCount;
    
    //|real4| 2*particleCount*real4| nodes*real4 | nodes*real4 | nodes*3*real4 |    
    //Info about #particles, #nodes, start and end of tree-walk
    //The particle positions and velocities
    //The nodeSizeData
    //The nodeCenterData
    //The multipole data    
    
    //Add the initial nodes to the curLevel list
    for(int i=start; i < end; i++)
    {
      curLevel.push_back(i);
    }
    
    //Add the nodes before the start and end to the node list
    for(int i=0; i < start; i++)
    {      
      dataBuffer[nodeSizeOffset++]   = nodeSizeInfo[i];
      dataBuffer[nodeCenterOffset++] = nodeCenterInfo[i];
      
//       //TEST CODE!!! REMOVE!!
//       multipole[i*3 + 0].y = i;            
//       printf("TEST: %d\t%f w: %f\n", i, multipole[i*3 + 0].x, multipole[i*3 + 0].w);
      
      dataBuffer[multiPoleOffset++]  = multipole[i*3 + 0];
      dataBuffer[multiPoleOffset++]  = multipole[i*3 + 1];
      dataBuffer[multiPoleOffset++]  = multipole[i*3 + 2];        
    }
    
    //Start the tree-walk
 /*   printf("Start: %d end: %d \n", start, end);
    cout << "Starting walk on: " << curLevel.size() << " items! \n"; 
    fprintf(stderr,"Group: %f %f %f by %f %f %f \n", boxCenter.x, boxCenter.y, boxCenter.z, boxSize.x, boxSize.y, boxSize.z);
   */ 
    //Variables to rewrite the tree-structure indices
    int childNodeOffset         = end;
    int childParticleOffset     = 0;
    

    while(curLevel.size() > 0)
    {
      for(unsigned int i=0; i < curLevel.size(); i++)
      {
        //Read node data
        int node         = curLevel[i];
        real4 nodeCenter = nodeCenterInfo[node];
        real4 nodeSize   = nodeSizeInfo[node];
        bool leaf        = nodeCenter.w <= 0;
        
        union{float f; int i;} u; //__float_as_int
        u.f           = nodeSize.w;      
        int childinfo = u.i;
        
       
        
        int child, nchild;
        if(!leaf)
        {
          //Node
          child    =    childinfo & 0x0FFFFFFF;                         //Index to the first child of the node
          nchild   = (((childinfo & 0xF0000000) >> 28)) ;         //The number of children this node has              
        }
        else
        {
          //Leaf
          child   =   childinfo & BODYMASK;                                     //thre first body in the leaf
          nchild  = (((childinfo & INVBMASK) >> LEAFBIT)+1);     //number of bodies in the leaf masked with the flag        
        }
  
  //TODO the opening check for Individual softening
        
        #ifdef INDSOFT
          //Very inefficient this but for testing I have to live with it...                                                                                   
          float node_eps_val = multipole[node*3 + 1].w; 
        #endif         
        
        #ifdef IMPBH
          //Improved barnes hut version
          float4 nodeCOM     = multipole[node*3 + 0];
          nodeCOM.w = nodeCenter.w;
          
          #ifdef INDSOFT
             bool   split   = split_node_grav_impbh(nodeCOM, boxCenter, boxSize, group_eps, node_eps_val);
          #else
             bool split = split_node_grav_impbh(nodeCOM, boxCenter, boxSize);  
          #endif
          
        #else
          //Minimal distance version
          
          #ifdef INDSOFT
            bool split = split_node(nodeCenter, nodeSize, boxCenter, boxSize, group_eps, node_eps_val);  //Check if node should be split
          #else
            bool split = split_node(nodeCenter, nodeSize, boxCenter, boxSize);
          #endif          
          
        #endif         
        
//         bool split = split_node(nodeCenter, nodeSize, boxCenter, boxSize);
       
//         bool split = true;
//         uint temp =childinfo;
//         uint temp = 0xFFFFFFFF;
        uint temp = 0;
        //if split & node add children to next lvl stack
        if(split && !leaf)
        {       
          for(int i=child; i < child+nchild; i++)
          {
            nextLevel.push_back(i);            
          }
          
          temp = childNodeOffset | (nchild << 28);  
          //Update reference to children
          childNodeOffset += nchild;
        }
        
        //if split & leaf add particles to particle list
        if(split && leaf)
        { 
          for(int i=child; i < child+nchild; i++)
          {
             dataBuffer[particleOffset++] = bodies[i];
             dataBuffer[velParticleOffset++] = velocities[i];
             massSum += bodies[i].w;
          }
          
          temp = childParticleOffset | ((nchild-1) << LEAFBIT);
          childParticleOffset += nchild;          
        }
        
        
        //Add the node data to the appropriate arrays
        //and modify the node reference
        //start ofset for its children, should be nodeCount at start of this level +numberofnodes on this level
        //plus a counter that counts the number of childs of the nodes we have tested

        //New childoffset:
        union{int i; float f;} itof; //__int_as_float
        itof.i           = temp;            
        float tempConv = itof.f;
        
        //Add node properties and update references
        real4 nodeSizeInfoTemp  = nodeSizeInfo[node];
        nodeSizeInfoTemp.w      = tempConv;             //Replace child reference
        
//         //TEST CODE!!! REMOVE!!
//         multipole[node*3 + 0].y = node;      
//         if(node < 20)
//         printf("TEST (%d): %d cent: %f %f %f size: %f %f %f w: %f\n", mpiGetRank(), node, 
//                nodeCenterInfo[node].x, nodeCenterInfo[node].y, nodeCenterInfo[node].z,
//                nodeSizeInfoTemp.x, nodeSizeInfoTemp.y, nodeSizeInfoTemp.z,
//                multipole[node*3 + 0].w);        
//         
        dataBuffer[nodeSizeOffset++]   = nodeSizeInfoTemp;
        dataBuffer[nodeCenterOffset++] = nodeCenterInfo[node];
        dataBuffer[multiPoleOffset++]  = multipole[node*3 + 0];
        dataBuffer[multiPoleOffset++]  = multipole[node*3 + 1];
        dataBuffer[multiPoleOffset++]  = multipole[node*3 + 2];       
        
   /*   dataBuffer[nodeSizeOffset-1].x = mpiGetRank()*10+1;
      dataBuffer[nodeCenterOffset-1].x = mpiGetRank()*10+2;
      
      dataBuffer[multiPoleOffset-3].x = mpiGetRank()*10+3;
      dataBuffer[multiPoleOffset-2].x = mpiGetRank()*10+4;
      dataBuffer[multiPoleOffset-1].x = mpiGetRank()*10+5;           
  */
        
        if(!split)
        {          
          massSum += multipole[node*3 + 0].w;
        }
      } //end for curLevel.size
      
      
      //Put next level stack into current level and continue
      curLevel.clear();
      
//       cout << "Next level: " << nextLevel.size() << endl;
      curLevel.assign(nextLevel.begin(), nextLevel.end());
      nextLevel.clear();
      
    }//end while
 
//   cout << "New offsets: "  << particleOffset << " \t" << nodeSizeOffset << " \t" << nodeCenterOffset << endl;
   cout << "Mass sum: " << massSum  << endl;
//   cout << "Mass sumtest: " << multipole[0*0 + 0].w << endl;
  
}




void octree::create_local_essential_tree_count(real4* bodies, real4* multipole, real4* nodeSizeInfo, real4* nodeCenterInfo,   
                                         double4 boxCenter, double4 boxSize, float group_eps, int start, int end, 
                                         int &particles, int &nodes)
{
    //Walk the tree as is done on device, level by level
    vector<int> curLevel;
    vector<int> nextLevel;

    int particleCount   = 0;
    int nodeCount       = 0;
    
    //Add the initial nodes to the curLevel list
    for(int i=start; i < end; i++)
    {
      curLevel.push_back(i);
    }
    
    //Add the nodes before the start and end to the node list
    for(int i=0; i < start; i++)
    {
      nodeCount++;
    }
      
    //Start the tree-walk       
    while(curLevel.size() > 0)
    {
      for(unsigned int i=0; i < curLevel.size(); i++)
      {
        //Read node data
        int node         = curLevel[i];
        real4 nodeCenter = nodeCenterInfo[node];
        real4 nodeSize   = nodeSizeInfo[node];
        bool leaf        = nodeCenter.w <= 0;
        
        union{float f; int i;} u; //__float_as_int
        u.f           = nodeSize.w;      
        int childinfo = u.i;
        
        int child, nchild;
        if(!leaf)
        {
          //Node
          child    =    childinfo & 0x0FFFFFFF;                         //Index to the first child of the node
          nchild   = (((childinfo & 0xF0000000) >> 28)) ;         //The number of children this node has              
        }
        else
        {
          //Leaf
          child   =   childinfo & BODYMASK;                                     //thre first body in the leaf
          nchild  = (((childinfo & INVBMASK) >> LEAFBIT)+1);     //number of bodies in the leaf masked with the flag        
        }
        
        #ifdef INDSOFT
          //Very inefficient this but for testing I have to live with it...                                                                                   
          float node_eps_val = multipole[node*3 + 1].w; 
        #endif         
        
        #ifdef IMPBH
          //Improved barnes hut version
          float4 nodeCOM     = multipole[node*3 + 0];
          nodeCOM.w = nodeCenter.w;
          
          #ifdef INDSOFT
             bool   split   = split_node_grav_impbh(nodeCOM, boxCenter, boxSize, group_eps, node_eps_val);
          #else
             bool split = split_node_grav_impbh(nodeCOM, boxCenter, boxSize);  
          #endif
          
        #else
          //Minimal distance version
          
          #ifdef INDSOFT
            bool split = split_node(nodeCenter, nodeSize, boxCenter, boxSize, group_eps, node_eps_val);  //Check if node should be split
          #else
            bool split = split_node(nodeCenter, nodeSize, boxCenter, boxSize);
          #endif          
          
        #endif 
        
//         bool split = split_node(nodeCenter, nodeSize, boxCenter, boxSize);

//         split = true;
        //if split & node add children to next lvl stack
        if(split && !leaf)
        {       
          for(int i=child; i < child+nchild; i++)
          {
            nextLevel.push_back(i);            
          }
        }
        
        //if split & leaf add particles to particle list
        if(split && leaf)
        { 
          for(int i=child; i < child+nchild; i++)
          {
            particleCount++;
          }
        }

        //Increase the nodeCount, since this node will be part of the tree-structure
        nodeCount++;        
      } //end for curLevel.size
      
      
      //Put next level stack into current level and continue
      curLevel.clear();
      
//       cout << "Next level: " << nextLevel.size() << endl;
      curLevel.assign(nextLevel.begin(), nextLevel.end());
      nextLevel.clear();
    }//end while
    
    particles = particleCount;
    nodes = nodeCount;  
    
    fprintf(stderr, "Count found: %d particles and %d nodes. Boxsize: (%f %f %f ) BoxCenter: (%f %f %f)\n", 
                    particles, nodes, boxSize.x ,boxSize.y, boxSize.z, boxCenter.x, boxCenter.y, boxCenter.z );   
}


real4* octree::MP_exchange_bhlist(int ibox, int isource, 
                                int bufferSize, real4 *letDataBuffer)
{
    MPI_Status status;
    int nrecvlist;
    int nlist = bufferSize;
    
    //first send&get the number of particles to send&get
    MPI_Sendrecv(&nlist,1,MPI_INT,ibox,procId*10, &nrecvlist,
                 1,MPI_INT,isource,isource*10,MPI_COMM_WORLD, &status);

    //Resize the buffer so it has the correct size and then exchange the tree 
    real4 *recvDataBuffer = new real4[nrecvlist];
    
    //Particles
    MPI_Sendrecv(&letDataBuffer[0], nlist*sizeof(real4), MPI_BYTE, ibox, procId*10+1,
                 &recvDataBuffer[0], nrecvlist*sizeof(real4), MPI_BYTE, isource, isource*10+1,
                 MPI_COMM_WORLD, &status);        
    
    return recvDataBuffer;             
} 



void octree::ICSend(int destination, real4 *bodyPositions, real4 *bodyVelocities,  int *bodiesIDs, int toSend)
{
    //First send the number of particles, then the actual sample data
    MPI_Send(&toSend, 1, MPI_INT, destination, destination*2 , MPI_COMM_WORLD);
    
    //Send the positions, velocities and ids
    MPI_Send( bodyPositions,  toSend*sizeof(real)*4, MPI_BYTE, destination, destination*2+1, MPI_COMM_WORLD);
    MPI_Send( bodyVelocities, toSend*sizeof(real)*4, MPI_BYTE, destination, destination*2+2, MPI_COMM_WORLD);
    MPI_Send( bodiesIDs,      toSend*sizeof(int),    MPI_BYTE, destination, destination*2+3, MPI_COMM_WORLD);
    
/*    MPI_Send( (real*)&bodyPositions[0],  toSend*sizeof(real)*4, MPI_BYTE, destination, destination*2+1, MPI_COMM_WORLD);
    MPI_Send( (real*)&bodyVelocities[0], toSend*sizeof(real)*4, MPI_BYTE, destination, destination*2+2, MPI_COMM_WORLD);
    MPI_Send( (int *)&bodiesIDs[0],      toSend*sizeof(int),    MPI_BYTE, destination, destination*2+3, MPI_COMM_WORLD);*/
}

void octree::ICRecv(int recvFrom, vector<real4> &bodyPositions, vector<real4> &bodyVelocities,  vector<int> &bodiesIDs)
{
   MPI_Status status;
   int nreceive;    
   int procId = mpiGetRank();
    
   //First send the number of particles, then the actual sample data
   MPI_Recv(&nreceive, 1, MPI_INT, recvFrom, procId*2, MPI_COMM_WORLD,&status);       
   
   bodyPositions.resize(nreceive);
   bodyVelocities.resize(nreceive);
   bodiesIDs.resize(nreceive);
    
   //Recv the positions, velocities and ids
   MPI_Recv( (real*)&bodyPositions[0],  nreceive*sizeof(real)*4, MPI_BYTE, recvFrom, procId*2+1, MPI_COMM_WORLD,&status);
   MPI_Recv( (real*)&bodyVelocities[0], nreceive*sizeof(real)*4, MPI_BYTE, recvFrom, procId*2+2, MPI_COMM_WORLD,&status);
   MPI_Recv( (int *)&bodiesIDs[0],      nreceive*sizeof(int),    MPI_BYTE, recvFrom, procId*2+3, MPI_COMM_WORLD,&status);
}
