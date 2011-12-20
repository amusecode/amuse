#include "octree.h"

void octree::mpiInit(int argc,char *argv[], int &procId, int &nProcs)
{
    int  namelen;   
    char processor_name[MPI_MAX_PROCESSOR_NAME];
    
    int mpiInitialized = 0;
    MPI_Initialized(&mpiInitialized);
    
    if(!mpiInitialized)
      MPI_Init(&argc,&argv);

    MPI_Comm_size(MPI_COMM_WORLD, &nProcs);
    MPI_Comm_rank(MPI_COMM_WORLD, &procId);
    MPI_Get_processor_name(processor_name,&namelen);

    #ifdef PRINT_MPI_DEBUG
      fprintf(stderr, "Proc id: %d @ %s , total processes: %d (mpiInit) \n", procId, processor_name, nProcs);
    #endif


    //Allocate memory for the used buffers
    domainRLow  = new double4[nProcs];
    domainRHigh = new double4[nProcs];

    domHistoryLow   = new int4[nProcs];
    domHistoryHigh  = new int4[nProcs];
    
    //Fill domainRX with constants so we can check if its initialized before
    for(int i=0; i < nProcs; i++)
    {
      domainRLow[i] = domainRHigh[i] = (double4){1e10, 1e10, 1e10, 1e10};
      
      domHistoryLow[i] = domHistoryHigh[i] = (int4){0,0,0,0};
    }
    
    currentRLow  = new double4[nProcs];
    currentRHigh = new double4[nProcs];
    
    xlowPrev  = new double4[nProcs];
    xhighPrev = new double4[nProcs];
    


    nSampleAndSizeValues    = new int2[nProcs];  
    curSysState             = new sampleRadInfo[nProcs];
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
//end utility



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
    
    nTotalFreq = tmp;
    
    

    #ifdef PRINT_MPI_DEBUG
    if(procId == 0)
      cout << "Total number of particles: " << tmp << endl;
    #endif
   
    int maxsample = (int)(NMAXSAMPLE*0.8); // 0.8 is safety factor    
    sampleFreq = (tmp+maxsample-1)/maxsample;
    
    fprintf(stderr,"Sample FREQ: %d \n", sampleFreq);
    
    prevSampFreq = sampleFreq;
    
}

//Uses one communication by storing data in one buffer
//nsample can be set to zero if this call is only used
//to get updated domain information
//This one is used when only updating the currentbox sizes
void octree::sendCurrentRadiusInfo(real4 &rmin, real4 &rmax)
{
  sampleRadInfo curProcState;
  
  int nsample               = 0; //Place holder to just use same datastructure
  curProcState.nsample      = nsample;
  curProcState.rmin         = (double4){rmin.x, rmin.y, rmin.z, rmin.w};
  curProcState.rmax         = (double4){rmax.x, rmax.y, rmax.z, rmax.w};
  
  //Get the number of sample particles and the domain size information
  MPI_Allgather(&curProcState, sizeof(sampleRadInfo), MPI_BYTE,  curSysState,
                sizeof(sampleRadInfo), MPI_BYTE, MPI_COMM_WORLD);  
  mpiSync();

  rmin.x                 = currentRLow[0].x = curSysState[0].rmin.x;
  rmin.y                 = currentRLow[0].y = curSysState[0].rmin.y;
  rmin.z                 = currentRLow[0].z = curSysState[0].rmin.z;
                           currentRLow[0].w = curSysState[0].rmin.w;
                           
  rmax.x                 = currentRHigh[0].x = curSysState[0].rmax.x;
  rmax.y                 = currentRHigh[0].y = curSysState[0].rmax.y;
  rmax.z                 = currentRHigh[0].z = curSysState[0].rmax.z;  
                           currentRHigh[0].w = curSysState[0].rmax.w;

  
  for(int i=1; i < nProcs; i++)
  {
    rmin.x = fmin(rmin.x, curSysState[i].rmin.x);
    rmin.y = fmin(rmin.y, curSysState[i].rmin.y);
    rmin.z = fmin(rmin.z, curSysState[i].rmin.z);
    
    rmax.x = fmax(rmax.x, curSysState[i].rmax.x);
    rmax.y = fmax(rmax.y, curSysState[i].rmax.y);
    rmax.z = fmax(rmax.z, curSysState[i].rmax.z);    
    
    currentRLow[i].x = curSysState[i].rmin.x;
    currentRLow[i].y = curSysState[i].rmin.y;
    currentRLow[i].z = curSysState[i].rmin.z;
    currentRLow[i].w = curSysState[i].rmin.w;
    
    currentRHigh[i].x = curSysState[i].rmax.x;
    currentRHigh[i].y = curSysState[i].rmax.y;
    currentRHigh[i].z = curSysState[i].rmax.z;
    currentRHigh[i].w = curSysState[i].rmax.w;    
  }                
}

//Uses one communication by storing data in one buffer
//nsample can be set to zero if this call is only used
//to get updated domain information
void octree::sendSampleAndRadiusInfo(int nsample, real4 &rmin, real4 &rmax)
{
  sampleRadInfo curProcState;
  
  curProcState.nsample      = nsample;
  curProcState.rmin         = (double4){rmin.x, rmin.y, rmin.z, rmin.w};
  curProcState.rmax         = (double4){rmax.x, rmax.y, rmax.z, rmax.w};
  
  globalRmax            = 0;
  totalNumberOfSamples  = 0;

  //Get the number of sample particles and the domain size information
  MPI_Allgather(&curProcState, sizeof(sampleRadInfo), MPI_BYTE,  curSysState,
                sizeof(sampleRadInfo), MPI_BYTE, MPI_COMM_WORLD);  

  rmin.x                 =  curSysState[0].rmin.x;
  rmin.y                 =  curSysState[0].rmin.y;
  rmin.z                 =  curSysState[0].rmin.z;
                           
  rmax.x                 =  curSysState[0].rmax.x;
  rmax.y                 =  curSysState[0].rmax.y;
  rmax.z                 =  curSysState[0].rmax.z;  

  totalNumberOfSamples   = curSysState[0].nsample;
  
  
  for(int i=1; i < nProcs; i++)
  {
    rmin.x = fmin(rmin.x, curSysState[i].rmin.x);
    rmin.y = fmin(rmin.y, curSysState[i].rmin.y);
    rmin.z = fmin(rmin.z, curSysState[i].rmin.z);
    
    rmax.x = fmax(rmax.x, curSysState[i].rmax.x);
    rmax.y = fmax(rmax.y, curSysState[i].rmax.y);
    rmax.z = fmax(rmax.z, curSysState[i].rmax.z);    
    
   
    totalNumberOfSamples   += curSysState[i].nsample;
  }                

  if(procId == 0)
  {
    if(fabs(rmin.x)>globalRmax)  globalRmax=fabs(rmin.x);
    if(fabs(rmin.y)>globalRmax)  globalRmax=fabs(rmin.y);
    if(fabs(rmin.z)>globalRmax)  globalRmax=fabs(rmin.z);
    if(fabs(rmax.x)>globalRmax)  globalRmax=fabs(rmax.x);
    if(fabs(rmax.y)>globalRmax)  globalRmax=fabs(rmax.y);
    if(fabs(rmax.z)>globalRmax)  globalRmax=fabs(rmax.z);
    
    if(totalNumberOfSamples > sampleArray.size())
    {
      sampleArray.resize(totalNumberOfSamples);
    }
  } 
}

void octree::gpu_collect_sample_particles(int nSample, real4 *sampleParticles)
{
  int *nReceiveCnts  = new int[nProcs];
  int *nReceiveDpls  = new int[nProcs];
  nReceiveCnts[0] = nSample*sizeof(real4);
  nReceiveDpls[0] = 0;
  
  if(procId == 0)
  {
    for(int i=1; i < nProcs; i++)
    {
      nReceiveCnts[i] = curSysState[i].nsample*sizeof(real4);
      nReceiveDpls[i] = nReceiveDpls[i-1] + nReceiveCnts[i-1];
    } 
  }

  //Collect sample particles
  MPI_Gatherv(&sampleParticles[0], nSample*sizeof(real4), MPI_BYTE,
              &sampleArray[0], nReceiveCnts, nReceiveDpls, MPI_BYTE,
              0, MPI_COMM_WORLD);
              
  delete[] nReceiveCnts;
  delete[] nReceiveDpls;             
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
    
    //Instead of using mpi-reduce but broadcast or something we can receive the 
    //individual values,so we dont haveto send them in the part below, saves 
    //communication time!!!  This function is only used once so no problem
    int *nSampleValues = new int[nProcs];
    int *nReceiveCnts  = new int[nProcs];
    int *nReceiveDpls  = new int[nProcs];
    MPI_Gather(&nsample, 1, MPI_INT, nSampleValues, 1, MPI_INT, 0, MPI_COMM_WORLD);
    
    //Increase the size of the result buffer if needed
    if(procId == 0)
    {
      //Sum the total amount of sample particles
      unsigned int totalNumberOfSamples = 0;
      
      for(int i=0; i < nProcs; i++)
      {
        totalNumberOfSamples += nSampleValues[i];
      }
      
      if(totalNumberOfSamples > sampleArray.size())
      {
        sampleArray.resize(totalNumberOfSamples);
      }
    }

    //Compute buffer and displacements for MPI_Gatherv
    nReceiveCnts[0] = nsample*sizeof(real4);
    nReceiveDpls[0] = 0;
    
    if(procId == 0)
    {
      for(int i=1; i < nProcs; i++)
      {
        nReceiveCnts[i] = nSampleValues[i]*sizeof(real4);
        nReceiveDpls[i] = nReceiveDpls[i-1] + nReceiveCnts[i-1];
      } 
    }

    //Collect sample particles, note the MPI_IN_PLACE to prevent MPI errors
    MPI_Gatherv((procId ? &sampleArray[0] : MPI_IN_PLACE), nsample*sizeof(real4), MPI_BYTE,
                &sampleArray[0], nReceiveCnts, nReceiveDpls, MPI_BYTE,
                0, MPI_COMM_WORLD);
 
    nsample = (nReceiveCnts[mpiGetNProcs()-1] +  nReceiveDpls[mpiGetNProcs()-1]) / sizeof(real4);
    
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
    
    delete[] nSampleValues;
    delete[] nReceiveCnts;
    delete[] nReceiveDpls;
}

void octree::createDistribution(real4 *bodies, int n_bodies)
{
  determine_sample_freq(n_bodies);

  vector<real4> sampleArray;
  sampleArray.reserve(NMAXSAMPLE);

  int     nsample;  //Number of samples for this process
  double  rmax;     //maximum coordinate used to create a box

  //Get the sample particles from the other processes
  collect_sample_particles(bodies, n_bodies, sampleFreq, sampleArray, nsample, rmax);
 
  //Now that we have a sample from all proces we setup the space division
  //Processor 0 determines the division
  if(procId == 0)
    determine_division(nsample, sampleArray,nx, ny, nz, rmax,domainRLow, domainRHigh);
  
  //Now broadcast the results to all other processes
  MPI_Bcast(domainRLow,  sizeof(double4)*nProcs,MPI_BYTE,0,MPI_COMM_WORLD);
  MPI_Bcast(domainRHigh, sizeof(double4)*nProcs,MPI_BYTE,0,MPI_COMM_WORLD);

  return;
}


/*
  Only update the box-sizes, box-boundaries
  of the different processes, do not do 
  anything related to sample particles
*/
void octree::gpu_updateDomainOnly()
{
  real4 r_min, r_max;
  //Get the current system/particle boundaries
  this->getBoundaries(localTree, r_min, r_max);
            
  int nSamples = 0;
  this->sendSampleAndRadiusInfo(nSamples, r_min, r_max);
  rMinGlobal = r_min;
  rMaxGlobal = r_max;  
}

/*
Get the domain boundaries 
Get the sample particles
Send them to process 0
Process 0 computes the new domain decomposition and broadcasts this
*/
void octree::gpu_updateDomainDistribution(double timeLocal)
{  
  my_dev::dev_stream aSyncStream;
  
  real4 r_min, r_max;                          
//    double t1 = get_time();
  //Get the current system/particle boundaries
  this->getBoundaries(localTree, r_min, r_max);
  
  int finalNRate; 
  
  //Average the previous and current execution time to make everything smoother
  //results in much better load-balance
  prevDurStep = (prevDurStep <= 0) ? timeLocal : prevDurStep;
  timeLocal   = (timeLocal + prevDurStep) / 2;

  double nrate = 0;
  #if 1
    //Only base load balancing on the computation time 
    double timeSum   = 0.0;

    //Sum the execution times
    MPI_Allreduce( &timeLocal, &timeSum, 1,MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);  
    
    nrate = timeLocal / timeSum; 
    
    if(1)       //Don't fluctuate particles too much
    {
      #define SAMPLING_LOWER_LIMIT_FACTOR  (1.9)

      double nrate2 = (double)localTree.n / (double) nTotalFreq;       
      nrate2       /= SAMPLING_LOWER_LIMIT_FACTOR;      
      
      if(nrate < nrate2)
      {
        nrate = nrate2;
      }
      
      double nrate2_sum = 0.0;
      MPI_Allreduce( &nrate, &nrate2_sum, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
      nrate /= nrate2_sum;        
    }
  #else
    //Equal number of particles
    nrate = (double)localTree.n / (double)nTotalFreq;    
  #endif
   
  int    nsamp    = (nTotalFreq *0.001) + 1;  //Total number of sample particles, global
  int nsamp_local = (nsamp*nrate) + 1;
  int nSamples    = nsamp_local;
  
  finalNRate      = localTree.n / nsamp_local;
  
  fprintf(stderr, "NSAMP [%d]: sample: %d nrate: %f finalrate: %d localTree.n: %d  \
                   previous: %d timeLocal: %f prevTimeLocal: %f \n",
                procId, nsamp_local, nrate, finalNRate, localTree.n, prevSampFreq,
                timeLocal, prevDurStep);  
                
  prevDurStep  = timeLocal; 
  prevSampFreq = finalNRate;
                

  my_dev::dev_mem<real4>  samplePositions(devContext);
  samplePositions.cmalloc_copy(localTree.generalBuffer1.get_pinned(), 
                          localTree.generalBuffer1.get_flags(), 
                          localTree.generalBuffer1.get_devMem(),
                          &localTree.generalBuffer1[0], 0,
                          nSamples, getAllignmentOffset(0));          
         

//   double t2 = get_time();
//   fprintf(stderr, "TIME1 (boundaries) %g \n", t2 - t1);  

  //Get the sample particles from the device and only copy 
  //the number of particles that is used
  //Action overlaps with the communication of domain boundary
  extractSampleParticles.set_arg<int>(0,     &localTree.n);
  extractSampleParticles.set_arg<int>(1,     &finalNRate);
  extractSampleParticles.set_arg<cl_mem>(2,  localTree.bodies_Ppos.p());
  extractSampleParticles.set_arg<cl_mem>(3,  samplePositions.p());
  extractSampleParticles.setWork(nSamples, 256);
  extractSampleParticles.execute(aSyncStream.s());            
  samplePositions.d2h(nSamples, false, aSyncStream.s());
  

  //Get number of sample particles per process and domain size information
  this->sendSampleAndRadiusInfo(nSamples, r_min, r_max);
  rMinGlobal = r_min;
  rMaxGlobal = r_max;
   
//    double t3 = get_time();
//   fprintf(stderr, "TIME2 (get and send sample info) %g \t %g \n", t3 - t2, t3-t1);            
  aSyncStream.sync();
  gpu_collect_sample_particles(nSamples, &samplePositions[0]);

  //double t4 = get_time();
  //fprintf(stderr, "TIME3 (get and send sample particles) %g \t %g \n", t4 - t3, t4-t1);            

  //Processor 0 determines the division
  if(procId == 0)
    determine_division(totalNumberOfSamples, sampleArray,nx, ny, nz, globalRmax, domainRLow, domainRHigh);
  
//   double t5 = get_time();
//   fprintf(stderr, "TIME4 (determ div ) %g \t %g \n", t5 - t4, t5-t1);

  //Now broadcast the results to all other processes
  MPI_Bcast(domainRLow,  sizeof(double4)*nProcs,MPI_BYTE,0,MPI_COMM_WORLD);
  MPI_Bcast(domainRHigh, sizeof(double4)*nProcs,MPI_BYTE,0,MPI_COMM_WORLD);
  
//   double t5 = get_time();
//   fprintf(stderr, "TIME4 (determ div and bcast) %g \t %g \n", t5 - t4, t5-t1);
//   fprintf(stderr, "TIME4 (Total sample part)  %g \n", t5-t1);       

//   if(this->nProcs > 1)
//   {
//     if(this->procId == 0)
//       for(int i = 0;i< this->nProcs;i++)     
//       {      
//         cerr << "Domain: " << i << " " << this->domainRLow[i].x << " " << this->domainRLow[i].y << " " << this->domainRLow[i].z << " " 
//                                        << this->domainRHigh[i].x << " " << this->domainRHigh[i].y << " " << this->domainRHigh[i].z <<endl;
//       }
//   }


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

bool sortByX (real4 i,real4 j) { return (i.x<j.x); }
bool sortByY (real4 i,real4 j) { return (i.y<j.y); }
bool sortByZ (real4 i,real4 j) { return (i.z<j.z); }


void octree::sortCoordinates2(real4 *r, int lo, int up, int cid )
{
  up += 1;
  if(cid == 0)
    std::sort(&r[lo], &r[up], sortByX);
  else if(cid == 1)
    std::sort(&r[lo], &r[up], sortByY);  
  else
    std::sort(&r[lo], &r[up], sortByZ);  
  
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

inline double computeNewCoordinate(double old, double newc, int prevChange)
{

  //return newc;

  int curChange = (fabs(old) > fabs(newc));

  double factor1 = 1, factor2 = 2, factor3 = 1;
  
  if(prevChange != curChange)
  {
    //Different direction take half a step                
    factor1 = 1; factor2 = 2;
  }
  else
  {
    //Same direction, take some bigger step (3/4th)
    factor1 = 3; factor2 = 4;

   //Same direction, take full step
//    factor3 = 0; factor1 = 1; factor2 = 1;
  }

 double temp = (factor3*old + factor1*newc) / factor2;
 
 return temp;
  
 //Default
 //return newc;
//avg
//  return (old+newc)/2;
//3/4:
//  return (old + 3*newc) / 4;

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
    
//     fprintf(stderr, "TIME4 TEST:  %d  %d \n", np-1, pos.size());
// 
//     double t1 = get_time();
    sortCoordinates(&pos[0], 0, np-1, 0);
//     sortCoordinates2(&pos[0], 0, np-1, 0);

//     double t2 = get_time();
//     fprintf(stderr, "TIME4 TEST: %g  %d \n", t2-t1, np-1);

    //Split the array in more or less equal parts
    for(int i = 0;i<n;i++)
    {
        istart[i] = (i*np)/n;
         //NOTE: It was i>= 0, changed this otherwise it writes before begin array...  
        if(i > 0 )
         iend[i-1]=istart[i]-1;
    }
    iend[n-1] = np-1; 
    
//     borderCnt++;
 
    //Split the x-axis
    for(int ix = 0;ix<nx;ix++)
    {
        double x0, x1;
        int ix0 = ix*ny*nz;
        int ix1 = (ix+1)*ny*nz;       
        calculate_boxdim(np, &pos[0], 0,istart[ix0],iend[ix1-1],rmax,x0,x1);                
        for(int i=ix0; i<ix1; i++)
        {
          //Check the domain borders and set a constant
          if(istart[ix0] == 0)
          {
            xlowPrev[i].x = 10e10;
          }
          if(iend[ix1-1] == np-1)
          {
            xhighPrev[i].x = 10e10;
          }

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
              //Check domain borders and set a constant
              if(istart[iy0]-istart[ix0] == 0)
              {
                xlowPrev[i].y = 10e10;
              }
              if( iend[iy1-1]-istart[ix0] == npy-1)
              {
                xhighPrev[i].y = 10e10;
              }              

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
                                 
                //Check the domain borders
                if(istart[iz0]-istart[iy0] == 0)
                {
                  xlowPrev[iz0].z = 10e10;
                }
                if(iend[iz0]-istart[iy0] == npz-1)
                {
                  xhighPrev[iz0].z = 10e10;
                }                                          
                                 
                xlow[iz0].z   = z0;
                xhigh[iz0].z  = z1;
            }
        }
    }


  //Do the magic to get a better load balance by changing the decompositon
  //slightly
  static bool isFirstStep  = true;
  
  if(!isFirstStep)
  {
    double temp;
    for(int i=0; i < nProcs; i++)
    {
//       fprintf(stderr,"DOMAIN LOW  %d || CUR: %f %f %f \tPREV: %f %f %f\n", i,
//               xlow[i].x, xlow[i].y, xlow[i].z,
//               xlowPrev[i].x, xlowPrev[i].y, xlowPrev[i].z);
//       fprintf(stderr,"DOMAIN HIGH %d || CUR: %f %f %f \tPREV: %f %f %f\n", i,
//               xhigh[i].x, xhigh[i].y, xhigh[i].z,
//               xhighPrev[i].x, xhighPrev[i].y, xhighPrev[i].z);
//               
       //Do magic !
       if(xlowPrev[i].x != 10e10)
       {          
         temp               = computeNewCoordinate(xlowPrev[i].x, xlow[i].x, domHistoryLow[i].x);
         domHistoryLow[i].x = (abs(xlowPrev[i].x) > fabs(xlow[i].x));
         xlow[i].x          = temp;
       }
       if(xhighPrev[i].x != 10e10)
       {         
         temp                 = computeNewCoordinate(xhighPrev[i].x, xhigh[i].x, domHistoryHigh[i].x);
         domHistoryHigh[i].x  = (abs(xhighPrev[i].x) > fabs(xhigh[i].x));
         xhigh[i].x           = temp;
       }
                         
       if(xlowPrev[i].y != 10e10)
       {         
         temp               = computeNewCoordinate(xlowPrev[i].y, xlow[i].y, domHistoryLow[i].y);
         domHistoryLow[i].y = (abs(xlowPrev[i].y) > fabs(xlow[i].y));
         xlow[i].y          = temp;         
       }
       if(xhighPrev[i].y != 10e10)
       {         
         temp                 = computeNewCoordinate(xhighPrev[i].y, xhigh[i].y, domHistoryHigh[i].y);
         domHistoryHigh[i].y  = (abs(xhighPrev[i].y) > fabs(xhigh[i].y));
         xhigh[i].y           = temp;
       }
       
       if(xlowPrev[i].z != 10e10)
       {         
         temp               = computeNewCoordinate(xlowPrev[i].z, xlow[i].z, domHistoryLow[i].z);
         domHistoryLow[i].z = (abs(xlowPrev[i].z) > fabs(xlow[i].z));
         xlow[i].z          = temp;         
       }
       if(xhighPrev[i].z != 10e10)
       {         
         temp                 = computeNewCoordinate(xhighPrev[i].z, xhigh[i].z, domHistoryHigh[i].z);
         domHistoryHigh[i].z  = (abs(xhighPrev[i].z) > fabs(xhigh[i].z));
         xhigh[i].z           = temp;
       }
      
    }
  }


  //Copy the current decomposition to the previous for next step
  for(int i=0; i < nProcs; i++)
  {
    xlowPrev[i]   = xlow[i];
    xhighPrev[i]  = xhigh[i];
  }
    
  isFirstStep = false;

  //Free up memory
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

//     int iret = 0;
//     int giret;
//     MPI_Allreduce(&iret, &giret,1, MPI_INT, MPI_MAX,MPI_COMM_WORLD);
//     return giret;
    return 0;
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


//Send particles to the appropriate processors
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
  real4  *predictedBodiesPositions = &tree.bodies_Ppos[0];
  real4  *predictedVelocities      = &tree.bodies_Pvel[0];    

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
      if(isinbox(predictedBodiesPositions[i], domainRLow[ibox], domainRHigh[ibox]))
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
        //Predicted position
        tmpp                           = predictedBodiesPositions[iloc];
        predictedBodiesPositions[iloc] = predictedBodiesPositions[i];
        predictedBodiesPositions[i]    = tmpp;
        //Predicted velocity
        tmpp                  = predictedVelocities[iloc];
        predictedVelocities[iloc] = predictedVelocities[i];
        predictedVelocities[i]    = tmpp;                
        //Time-step
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
          body.Ppos  = predictedBodiesPositions[iloc];          
          body.Pvel  = predictedVelocities[iloc];        
          array2Send.push_back(body);
        }

        iloc++;
      }// end if
    }//for i=iloc
    nparticles[ibox] = iloc-firstloc[ibox];//Number of particles that has to be send to proc: ibox
  } // for(int ib=0;ib<nproc;ib++)
  
  printf("Requires search time: %lg ,proc: %d found in our own box: %d n: %d  send to others: %ld \n", 
         get_time()-t1, myid, nparticles[myid], tree.n, array2Send.size());
  
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
  
  printf("Required inter-process communication time: %lg ,proc: %d\n", get_time()-t1, myid);  
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
  tree.setN(idfirst+recvCount);    
  tree.bodies_pos.cresize (idfirst+recvCount + 1, false);  
  tree.bodies_acc0.cresize(idfirst+recvCount,     false);
  tree.bodies_acc1.cresize(idfirst+recvCount,     false);
  tree.bodies_vel.cresize (idfirst+recvCount,     false);
  tree.bodies_time.cresize(idfirst+recvCount,     false);
  tree.bodies_ids.cresize (idfirst+recvCount + 1, false);
  tree.bodies_Ppos.cresize(idfirst+recvCount + 1, false);  
  tree.bodies_Pvel.cresize(idfirst+recvCount + 1, false);  
  
  //This one has to be at least the same size as the number of particles inorder to
  //have enough space to store the other buffers 
  tree.generalBuffer1.cresize(3*(idfirst+recvCount)*4, false);
  
  printf("Benodigde gpu malloc tijd stap 1: %lg \t Size: %d \tRank: %d \t Size: %d \n", 
         get_time()-t1, idfirst+recvCount, mpiGetRank(), tree.bodies_Ppos.get_size()); 
  t1 = get_time();

  tempidFirst = idfirst; tempRecvCount = recvCount;

  //Copy data from struct into the main arrays
  for(unsigned int P=0; P < recvCount; P++)
  {
    tree.bodies_pos[idfirst+P]  = recv_buffer3[P].pos;        tree.bodies_vel[idfirst+P]      = recv_buffer3[P].vel;
    tree.bodies_acc0[idfirst+P] = recv_buffer3[P].acc0;       tree.bodies_acc1[idfirst+P]     = recv_buffer3[P].acc1;
    tree.bodies_time[idfirst+P] = recv_buffer3[P].time;       tree.bodies_ids[idfirst+P]      = recv_buffer3[P].id;
    tree.bodies_Ppos[idfirst+P] = recv_buffer3[P].Ppos;       tree.bodies_Pvel[idfirst+P]     = recv_buffer3[P].Pvel;
  }

  printf("Required DATA in struct copy time: %lg \n", get_time()-t1); t1 = get_time();

 
  if(ibend == -1){
    
  }else{
      //Something went wrong
    cerr << "ERROR in exchange_particles_with_overflow_check! \n"; exit(0);
  }

  
  //Resize the arrays of the tree    
  reallocateParticleMemory(tree);   
     
  printf("Required gpu malloc time step 2: %lg \n", get_time()-t1);
  printf("Total GPU interaction time: %lg \n", get_time()-t2);

  int retValue = 0;


  delete[] firstloc;
  delete[] nparticles;

  return retValue;
}

//Function that uses the GPU to get a set of particles that have to be 
//send to other processes
void octree::gpuRedistributeParticles()
{
  //Memory buffers to hold the extracted particle information
  my_dev::dev_mem<uint>  validList(devContext);
  my_dev::dev_mem<uint>  compactList(devContext);

  compactList.cmalloc_copy(localTree.generalBuffer1.get_pinned(), 
                          localTree.generalBuffer1.get_flags(), 
                          localTree.generalBuffer1.get_devMem(),
                          &localTree.generalBuffer1[0], 0,
                          localTree.n, getAllignmentOffset(0));
                                    
  validList.cmalloc_copy(localTree.generalBuffer1.get_pinned(), 
                            localTree.generalBuffer1.get_flags(), 
                            localTree.generalBuffer1.get_devMem(),
                            &localTree.generalBuffer1[localTree.n], localTree.n,
                            localTree.n, getAllignmentOffset(localTree.n));
                            
  double4 thisXlow  = domainRLow [this->procId];
  double4 thisXhigh = domainRHigh[this->procId];
                                                                          
  domainCheck.set_arg<int>(0,     &localTree.n);
  domainCheck.set_arg<double4>(1, &thisXlow);
  domainCheck.set_arg<double4>(2, &thisXhigh);
  domainCheck.set_arg<cl_mem>(3,  localTree.bodies_Ppos.p());          
  domainCheck.set_arg<cl_mem>(4,  validList.p());
  domainCheck.setWork(localTree.n, 128);
  domainCheck.execute();            
  
  //Create a list of valid and invalid particles
  int validCount;
  gpuSplit(devContext, validList, compactList, localTree.n, &validCount);                 
       

  //Check if the memory size, of the generalBuffer is large enough to store the exported particles
  int tempSize = localTree.generalBuffer1.get_size() - localTree.n;
  int needSize = 1.01*(validCount*(sizeof(bodyStruct)/sizeof(int)));

  if(tempSize < needSize)
  {
    int itemsNeeded = needSize + localTree.n + 4092; //Slightly larger as before for offset space
    
    //Copy the compact list to the host we need this list intact
    compactList.d2h();
    int *tempBuf = new int[localTree.n];
    memcpy(tempBuf, &compactList[0], localTree.n*sizeof(int));

    //Resize the general buffer
    localTree.generalBuffer1.cresize(itemsNeeded, false);
    //Reset memory
    compactList.cmalloc_copy(localTree.generalBuffer1.get_pinned(),
                          localTree.generalBuffer1.get_flags(),
                          localTree.generalBuffer1.get_devMem(),
                          &localTree.generalBuffer1[0], 0,
                          localTree.n, getAllignmentOffset(0));

    //Restore the compactList
    memcpy(&compactList[0], tempBuf, localTree.n*sizeof(int));
    compactList.h2d();

    delete[] tempBuf;
  }
  
  my_dev::dev_mem<bodyStruct>  bodyBuffer(devContext);
  
  bodyBuffer.cmalloc_copy(localTree.generalBuffer1.get_pinned(), 
                  localTree.generalBuffer1.get_flags(), 
                  localTree.generalBuffer1.get_devMem(),
                  &localTree.generalBuffer1[localTree.n], localTree.n, 
                  validCount, getAllignmentOffset(localTree.n));                                 

  extractOutOfDomainBody.set_arg<int>(0,    &validCount);
  extractOutOfDomainBody.set_arg<cl_mem>(1, compactList.p());          
  extractOutOfDomainBody.set_arg<cl_mem>(2, localTree.bodies_Ppos.p());          
  extractOutOfDomainBody.set_arg<cl_mem>(3, localTree.bodies_Pvel.p());     
  extractOutOfDomainBody.set_arg<cl_mem>(4, localTree.bodies_pos.p());     
  extractOutOfDomainBody.set_arg<cl_mem>(5, localTree.bodies_vel.p());     
  extractOutOfDomainBody.set_arg<cl_mem>(6, localTree.bodies_acc0.p());
  extractOutOfDomainBody.set_arg<cl_mem>(7, localTree.bodies_acc1.p());
  extractOutOfDomainBody.set_arg<cl_mem>(8, localTree.bodies_time.p());
  extractOutOfDomainBody.set_arg<cl_mem>(9, localTree.bodies_ids.p());
  extractOutOfDomainBody.set_arg<cl_mem>(10, bodyBuffer.p());
  extractOutOfDomainBody.setWork(validCount, 128);
  extractOutOfDomainBody.execute();
  
  bodyBuffer.d2h(validCount);      
  
  //Now we have to move particles from the back of the array to the invalid spots
  //this can be done in parallel with exchange operation to hide some time

  //One integer for counting, true-> initialize to zero so counting starts at 0
  my_dev::dev_mem<uint>  atomicBuff(devContext, 1, true);

  //Internal particle movement
  internalMove.set_arg<int>(0,    &validCount);
  internalMove.set_arg<int>(1,    &localTree.n);
  internalMove.set_arg<double4>(2,    &thisXlow);
  internalMove.set_arg<double4>(3,    &thisXhigh);
  internalMove.set_arg<cl_mem>(4, compactList.p());
  internalMove.set_arg<cl_mem>(5, atomicBuff.p());               
  internalMove.set_arg<cl_mem>(6, localTree.bodies_Ppos.p());          
  internalMove.set_arg<cl_mem>(7, localTree.bodies_Pvel.p());     
  internalMove.set_arg<cl_mem>(8, localTree.bodies_pos.p());     
  internalMove.set_arg<cl_mem>(9, localTree.bodies_vel.p());     
  internalMove.set_arg<cl_mem>(10, localTree.bodies_acc0.p());
  internalMove.set_arg<cl_mem>(11, localTree.bodies_acc1.p());
  internalMove.set_arg<cl_mem>(12, localTree.bodies_time.p());
  internalMove.set_arg<cl_mem>(13, localTree.bodies_ids.p());
  internalMove.setWork(validCount, 128);
  internalMove.execute(execStream->s());    
  
  this->gpu_exchange_particles_with_overflow_check(localTree, &bodyBuffer[0], compactList, validCount);  

} //End gpuRedistributeParticles



//Exchange particles with other processes
int octree::gpu_exchange_particles_with_overflow_check(tree_structure &tree, 
                                                       bodyStruct *particlesToSend,
                                                       my_dev::dev_mem<uint> &extractList,
                                                       int nToSend)
{
  int myid      = procId;
  int nproc     = nProcs;
  int iloc      = 0;
  int nbody     = nToSend;


  bodyStruct  tmpp; 

  int *firstloc   = new int[nProcs+1];
  int *nparticles = new int[nProcs+1];

  // Loop over particles and determine which particle needs to go where
  // reorder the bodies in such a way that bodies that have to be send
  // away are stored after each other in the array
  double t1 = get_time();

  //Array reserve some memory at forehand , 1%
  vector<bodyStruct> array2Send;
  array2Send.reserve((int)(nToSend*1.5));
  
  for(int ib=0;ib<nproc;ib++)
  {
    int ibox       = (ib+myid)%nproc;
    firstloc[ibox] = iloc;      //Index of the first particle send to proc: ibox

    for(int i=iloc; i<nbody;i++)
    {
      if(isinbox(particlesToSend[i].Ppos, domainRLow[ibox], domainRHigh[ibox]))
      {
        //Reorder the particle information
        tmpp                  = particlesToSend[iloc];
        particlesToSend[iloc] = particlesToSend[i];
        particlesToSend[i]    = tmpp;

        //Put the particle in the array of to send particles
        array2Send.push_back(particlesToSend[iloc]);
        
        iloc++;
      }// end if
    }//for i=iloc
    nparticles[ibox] = iloc-firstloc[ibox];//Number of particles that has to be send to proc: ibox
  } // for(int ib=0;ib<nproc;ib++)

  
//   printf("Required search time: %lg ,proc: %d found in our own box: %d n: %d  to others: %ld \n", 
//          get_time()-t1, myid, nparticles[myid], tree.n, array2Send.size());

  
  if(iloc < nbody)
  {
      cerr << procId <<" exchange_particle error: particle in no box...iloc: " << iloc 
                     << " and nbody: " << nbody << "\n";
      exit(0);
  }
             

  /*totalsent = nbody - nparticles[myid];

  int tmp;
  MPI_Reduce(&totalsent,&tmp,1, MPI_INT, MPI_SUM,0,MPI_COMM_WORLD);
  if(procId == 0)
  {
    totalsent = tmp;
    cout << "Exchanged particles = " << totalsent << endl;
  }*/

  t1 = get_time();  
  
  //Allocate two times the amount of memory of that which we send
  vector<bodyStruct> recv_buffer3(nbody*2);
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
  
 
  if(ibend == -1){
    
  }else{
      //Something went wrong
    cerr << "ERROR in exchange_particles_with_overflow_check! \n"; exit(0);
  }
  
  
  printf("Required inter-process communication time: %lg ,proc: %d\n", 
         get_time()-t1, myid);  

  //Compute the new number of particles:
  int newN = tree.n + recvCount - nToSend;
  
  execStream->sync();   //make certain that the particle movement on the device
                        //is complete before we resize
                        
  //Have to resize the bodies vector to keep the numbering correct 
  //but do not reduce the size since we need to preserve the particles
  //in the oversized memory
  tree.bodies_pos.cresize (newN + 1, false);  
  tree.bodies_acc0.cresize(newN,     false);
  tree.bodies_acc1.cresize(newN,     false);
  tree.bodies_vel.cresize (newN,     false);
  tree.bodies_time.cresize(newN,     false);
  tree.bodies_ids.cresize (newN + 1, false);
  tree.bodies_Ppos.cresize(newN + 1, false);  
  tree.bodies_Pvel.cresize(newN + 1, false);  
  
  //This one has to be at least the same size as the number of particles inorder to
  //have enough space to store the other buffers 
  //Can only be resized after we are done since we still have
  //parts of memory pointing to that buffer (extractList)
  //Note that we allocate some extra memory to make everything texture/memory alligned
  tree.generalBuffer1.cresize(3*(newN)*4 + 4096, false);  

  //Now we have to copy the data in batches incase the generalBuffer1 is not large enough
  //Amount we can store:
  int spaceInIntSize    = 3*(newN)*4; 
  int newParticleSpace  = spaceInIntSize / (sizeof(bodyStruct) / sizeof(int));        
  int stepSize = newParticleSpace;

  my_dev::dev_mem<bodyStruct>  bodyBuffer(devContext);
  bodyBuffer.cmalloc_copy(localTree.generalBuffer1.get_pinned(), 
                  localTree.generalBuffer1.get_flags(), 
                  localTree.generalBuffer1.get_devMem(),
                  &localTree.generalBuffer1[0], 0, 
                  stepSize, getAllignmentOffset(0));           
 
  fprintf(stderr, "Exchange, received particles: (%d): %d \tnewN: %d\tItems that can be insert in one step: %d\n", 
                   procId, recvCount, newN, stepSize);     

  int insertOffset = 0;                   
  for(unsigned int i=0; i < recvCount; i+= stepSize)
  {
    int items = min(stepSize, (int)(recvCount-i));
      
    if(items > 0)
    {
      //Copy the data from the MPI receive buffers into the GPU-send buffer
      memcpy(&bodyBuffer[0], &recv_buffer3[insertOffset], sizeof(bodyStruct)*items);

      bodyBuffer.h2d(items);
  
//       int threads = max(nToSend, (int)recvCount);
            
      //Start the kernel that puts everything in place
      insertNewParticles.set_arg<int>(0,    &nToSend);
      insertNewParticles.set_arg<int>(1,    &items);
      insertNewParticles.set_arg<int>(2,    &tree.n);
      insertNewParticles.set_arg<int>(3,    &insertOffset);
      insertNewParticles.set_arg<cl_mem>(4, localTree.bodies_Ppos.p());          
      insertNewParticles.set_arg<cl_mem>(5, localTree.bodies_Pvel.p());     
      insertNewParticles.set_arg<cl_mem>(6, localTree.bodies_pos.p());     
      insertNewParticles.set_arg<cl_mem>(7, localTree.bodies_vel.p());     
      insertNewParticles.set_arg<cl_mem>(8, localTree.bodies_acc0.p());
      insertNewParticles.set_arg<cl_mem>(9, localTree.bodies_acc1.p());
      insertNewParticles.set_arg<cl_mem>(10, localTree.bodies_time.p());
      insertNewParticles.set_arg<cl_mem>(11, localTree.bodies_ids.p());
      insertNewParticles.set_arg<cl_mem>(12, bodyBuffer.p());
      insertNewParticles.setWork(items, 128);
      insertNewParticles.execute(); 
    }
    
    insertOffset += items;    
  }

//   printf("Benodigde gpu malloc tijd stap 1: %lg \t Size: %d \tRank: %d \t Size: %d \n", 
//          get_time()-t1, newN, mpiGetRank(), tree.bodies_Ppos.get_size()); 
//   t1 = get_time();
  
   
  tree.setN(newN);   
 
  //Resize the arrays of the tree    
  reallocateParticleMemory(tree);   
     
//   printf("Benodigde gpu malloc tijd stap 2: %lg \n", get_time()-t1);
//   printf("Totale GPU interactie tijd: %lg \n", get_time()-t2);

  int retValue = 0;

  delete[] firstloc;
  delete[] nparticles;

  return retValue;
}

//Local essential tree functions


void octree::essential_tree_exchange(vector<real4> &treeStructure, tree_structure &tree, tree_structure &remote)
{
  int myid    = procId;
  int nproc   = nProcs;
  int isource = 0;
 
  bool mergeOwntree = false;          //Default do not include our own tree-structre, thats mainly used for testing 
  int step          = nProcs - 1;     //Default merge all remote trees into one structure
//   step              = 1;
  int level_start   = 2;              //Default depth of where to start the tree-walk, default 2
  int procTrees     = 0;              //Number of trees that we've received and processed

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

  //creates a new array of pointers to int objects, with space for the local tree
  treeBuffers  = new real4*[mpiGetNProcs()]; 
  
  //Timers for the LET Exchange
  static double totalLETExTime    = 0;
//   double thisPartLETExTime = 0;
  thisPartLETExTime = 0;
  double tStart = 0;

//   for(int z=nproc-1; z > 0; z-=step)
  for(int z=nproc-1; z > 0; )
  {
    tStart = get_time();
    
    step = min(step, z);
    
    int recvTree = 0;    
    //For each process
    for(int ib = z; recvTree < step; recvTree++, ib--)
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

      double4 boxCenter = {     0.5*(currentRLow[ibox].x  + currentRHigh[ibox].x),
                                0.5*(currentRLow[ibox].y  + currentRHigh[ibox].y),
                                0.5*(currentRLow[ibox].z  + currentRHigh[ibox].z), 0};
      double4 boxSize   = {fabs(0.5*(currentRHigh[ibox].x - currentRLow[ibox].x)),
                           fabs(0.5*(currentRHigh[ibox].y - currentRLow[ibox].y)),
                           fabs(0.5*(currentRHigh[ibox].z - currentRLow[ibox].z)), 0};  
                              
                              
//       printf("Other proc min and max: [%f %f %f] \t [%f %f %f] \n", currentRLow[ibox].x, currentRLow[ibox].y, currentRLow[ibox].z,
//               currentRHigh[ibox].x, currentRHigh[ibox].y, currentRHigh[ibox].z);
    
    //   printf("Other proc center and size: [%f %f %f] \t [%f %f %f] \n", boxCenter.x, boxCenter.y, boxCenter.z,
    //          boxSize.x, boxSize.y, boxSize.z);

      uint2 node_begend;
      node_begend.x   = tree.level_list[level_start].x;
      node_begend.y   = tree.level_list[level_start].y;
      
      int particleCount, nodeCount;
      
  //     double t1 = get_time();    
      create_local_essential_tree_count(bodies, multipole, nodeSizeInfo, nodeCenterInfo,
                                  boxCenter, boxSize, currentRLow[ibox].w, node_begend.x, node_begend.y,
                                  particleCount, nodeCount);    
  //     printf("LET count: %lg\n", get_time()-t1);   
      //Buffer that will contain all the data:
      //|real4| 2*particleCount*real4| nodes*real4 | nodes*real4 | nodes*3*real4 |
      //1 + 2*particleCount + nodeCount + nodeCount + 3*nodeCount
      

      //Increase the number of particles and the number of nodes by the texture-offset such that these are correctly
      //aligned in memory
      particleCount += getTextureAllignmentOffset(particleCount, sizeof(real4));
      nodeCount     += getTextureAllignmentOffset(nodeCount    , sizeof(real4));
      
      //0-1 )                               Info about #particles, #nodes, start and end of tree-walk
      //1- Npart)                           The particle positions
      //1+Npart-Npart )                     The particle velocities
      //1+2*Npart-Nnode )                   The nodeSizeData
      //1+*2Npart+Nnode - Npart+2*Nnode )   The nodeCenterData
      //1+2*Npart+2*Nnode - Npart+5*Nnode ) The multipole data, is 3x number of nodes (mono and quadrupole data)
      int bufferSize = 1 + 2*particleCount + 5*nodeCount;
      real4 *letDataBuffer = new real4[bufferSize];
          
      create_local_essential_tree_fill(bodies, velocities, multipole, nodeSizeInfo, nodeCenterInfo,
                                  boxCenter, boxSize, currentRLow[ibox].w, node_begend.x, node_begend.y,
                                  particleCount, nodeCount, letDataBuffer);        
                                  
  /*    
      printf("LET count&fill: %lg\n", get_time()-t1);   
      t1 = get_time();   
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

      //Exchange the data of the tree structures  between the processes
      treeBuffers[recvTree] = MP_exchange_bhlist(ibox, isource, bufferSize, letDataBuffer);

      delete[] letDataBuffer;

      //This determines if we interrupt the exchange by starting a gravity kernel on the GPU
      if(execStream->isFinished())
      {
        fprintf(stderr,"GRAVFINISHED %d recvTree: %d\n", procId, recvTree);
        recvTree++;
        break;
      }
    }//end for each process
    

  
    z-=recvTree;
  
    //Now we have to merge the seperate tree-structures into one process

//     double t1 = get_time();

    int PROCS = recvTree;
    
    procTrees += recvTree;

    if(mergeOwntree) 
    {      
      //Add the processors own tree to the LET tree
      int particleCount   = tree.n;
      int nodeCount       = tree.n_nodes;
      
      int realParticleCount = tree.n;
      int realNodeCount     = tree.n_nodes;
      
      particleCount += getTextureAllignmentOffset(particleCount, sizeof(real4));
      nodeCount     += getTextureAllignmentOffset(nodeCount    , sizeof(real4));
    
      int bufferSizeLocal = 1 + 2*particleCount + 5*nodeCount;
      
      treeBuffers[PROCS]  = new real4[bufferSizeLocal];

      //Note that we use the real*Counts otherwise we read out of the array boundaries!!
      int idx = 1;
      memcpy(&treeBuffers[PROCS][idx], &bodies[0],         sizeof(real4)*realParticleCount);
      idx += particleCount;
      memcpy(&treeBuffers[PROCS][idx], &velocities[0],     sizeof(real4)*realParticleCount);
      idx += particleCount;
      memcpy(&treeBuffers[PROCS][idx], &nodeSizeInfo[0],   sizeof(real4)*realNodeCount);
      idx += nodeCount;
      memcpy(&treeBuffers[PROCS][idx], &nodeCenterInfo[0], sizeof(real4)*realNodeCount); 
      idx += nodeCount;
      memcpy(&treeBuffers[PROCS][idx], &multipole[0],      sizeof(real4)*realNodeCount*3);   
      
      treeBuffers[PROCS][0].x = particleCount;
      treeBuffers[PROCS][0].y = nodeCount;
      treeBuffers[PROCS][0].z = tree.level_list[level_start].x;
      treeBuffers[PROCS][0].w = tree.level_list[level_start].y;  
      PROCS                   = PROCS + 1; //Signal that we added one more tree-structure
      mergeOwntree            = false;     //Set it to false incase we do not merge all trees at once, we only inlcude our own once
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
    nodesBegEnd[mpiGetNProcs()].x   = nodesBegEnd[mpiGetNProcs()].y = 0; //Make valgrind happy
    int totalTopNodes               = 0;
    

    //Calculate the offsets
    for(int i=0; i < PROCS ; i++)
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
    int totalParticles    = particleSumOffsets[PROCS];
    int totalNodes        = nodeSumOffsets[PROCS];
    
    //To bind parts of the memory to different textures, the memory start address
    //has to be aligned with XXX bytes, so nodeInformation*sizeof(real4) has to be
    //increased by an offset, so that the node data starts at a XXX byte boundary
    //this is already done on the sending process, but since we modify the structure
    //it has to be done again
    int nodeTextOffset = getTextureAllignmentOffset(totalNodes+totalTopNodes, sizeof(real4));
    
    //Compute the total size of the buffer
    int bufferSize     = 2*(totalParticles) + 5*(totalNodes+totalTopNodes + nodeTextOffset);
    
    thisPartLETExTime += get_time() - tStart;
    //Allocate memory on host and device to store the merged tree-structure
    if(bufferSize > remote.fullRemoteTree.get_size())
    {
      //Can only resize if we are sure the LET is not running
      if(letRunning)
      {
        execStream->sync();     //Wait till the LET run is finished        
      }      
      remote.fullRemoteTree.cresize(bufferSize, false);  //Change the size but ONLY if we need more memory  
    }
    tStart = get_time();

    real4 *combinedRemoteTree = &remote.fullRemoteTree[0];
    
    //Copy all the pieces of the different trees at the correct memory offsets
    for(int i=0; i < PROCS; i++)
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
      memcpy(&combinedRemoteTree[(totalParticles) + particleSumOffsets[i]],   
            &treeBuffers[i][1+remoteP], sizeof(real4)*remoteP);
  
      //The start nodes, nodeSizeInfo
      memcpy(&combinedRemoteTree[2*(totalParticles) + startNodeSumOffsets[i]],  
            &treeBuffers[i][1+2*remoteP+remoteB], //From the start node onwards
            sizeof(real4)*remoteNstart);
      
      //Non start nodes, nodeSizeInfo       
      memcpy(&combinedRemoteTree[2*(totalParticles) +  totalTopNodes + nodeSumOffsets[i]],  
            &treeBuffers[i][1+2*remoteP+remoteE], //From the last start node onwards
            sizeof(real4)*(remoteN-remoteE));    
    
      //The start nodes, nodeCenterInfo
      memcpy(&combinedRemoteTree[2*(totalParticles) + startNodeSumOffsets[i]
                                  + (totalNodes + totalTopNodes + nodeTextOffset)],  
            &treeBuffers[i][1+2*remoteP+remoteB + remoteN], //From the start node onwards
            sizeof(real4)*remoteNstart);
      
      //Non start nodes, nodeCenterInfo       
      memcpy(&combinedRemoteTree[2*(totalParticles) +  totalTopNodes 
            + nodeSumOffsets[i] + (totalNodes + totalTopNodes + nodeTextOffset)],  
            &treeBuffers[i][1+2*remoteP+remoteE + remoteN], //From the last start node onwards
            sizeof(real4)*(remoteN-remoteE));   
          
      //The start nodes, multipole       
      memcpy(&combinedRemoteTree[2*(totalParticles) + 3*startNodeSumOffsets[i] +
            2*(totalNodes+totalTopNodes + nodeTextOffset)],  
            &treeBuffers[i][1+2*remoteP+2*remoteN + 3*remoteB], //From the start node onwards
            sizeof(real4)*remoteNstart*3);  
                      
      //Non start nodes, multipole       
      memcpy(&combinedRemoteTree[2*(totalParticles) +  3*totalTopNodes + 
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
      with the new correct locations and references. This takes place in two steps:
      First  the top nodes 
      Second the normal nodes
      Has to be done in two steps since they are not continous in memory if NPROCS > 2
      */  

      //Modify the top nodes
      int modStart = 2*(totalParticles) + startNodeSumOffsets[i];
      int modEnd   = modStart           + remoteNstart;

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
        }//end !leaf      
        combinedRemoteTree[j].w =  host_int_as_float(child);      //store the modified offset
      } 
      
      //Now the non-top nodes for this process
      modStart =  totalTopNodes + nodeSumOffsets[i] + 2*(totalParticles);
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
   
          if(nchild == 0) child = 0;                              //To prevent incorrect negative values
        }else{ //Leaf
          child   =   childinfo & BODYMASK;                       //the first body in the leaf
          nchild  = (((childinfo & INVBMASK) >> LEAFBIT)+1);      //number of bodies in the leaf masked with the flag    
          
          child =  child + particleSumOffsets[i];                 //Modify the particle offsets
          child = child | ((nchild-1) << LEAFBIT);                //Merging the data back into one int
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
    
    //Store the tree properties (number of particles, number of nodes, start and end topnode)
    remote.remoteTreeStruct.x = totalParticles;
    remote.remoteTreeStruct.y = totalNodes+totalTopNodes;
    remote.remoteTreeStruct.z = nodeTextOffset;  
    totalTopNodes             = (0 << 16) | (totalTopNodes);  //If its a merged tree we start at 0  
    remote.remoteTreeStruct.w = totalTopNodes;
    
//     fprintf(stderr,"Modifying the LET took: %g \n", get_time()-t1);
    fprintf(stderr,"Number of local bodies: %d number LET bodies: %d , number LET nodes: %d top nodes: %d Processed trees: %d (%d) \n",
                    tree.n, totalParticles, totalNodes, totalTopNodes, PROCS, procTrees);

    delete[] particleSumOffsets;
    delete[] nodeSumOffsets;
    delete[] startNodeSumOffsets;
    delete[] nodesBegEnd;
    
    
    thisPartLETExTime += get_time() - tStart;

    
    //Check if we need to summarize which particles are active, 
    //only done during the last approximate_gravity_let call
    bool doActivePart = (procTrees == mpiGetNProcs() -1);

    approximate_gravity_let(this->localTree, this->remoteTree, bufferSize, doActivePart);    
    
 
    
  } //end z
  delete[] treeBuffers;  
  
  totalLETExTime += thisPartLETExTime;
  
  fprintf(stderr,"LETEX [%d] curStep: %g\t   Total: %g \n", procId, thisPartLETExTime, totalLETExTime);
  
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
  float3 dr = {fabs(boxCenter.x - nodeCOM.x) - (boxSize.x),
               fabs(boxCenter.y - nodeCOM.y) - (boxSize.y),
               fabs(boxCenter.z - nodeCOM.z) - (boxSize.z)};

  dr.x += fabs(dr.x); dr.x *= 0.5f;
  dr.y += fabs(dr.y); dr.y *= 0.5f;
  dr.z += fabs(dr.z); dr.z *= 0.5f;

  //Distance squared, no need to do sqrt since opening criteria has been squared
  float ds2    = dr.x*dr.x + dr.y*dr.y + dr.z*dr.z;

  #ifdef INDSOFT    
    if(ds2      <= ((group_eps + node_eps ) * (group_eps + node_eps) ))           return true;
    //Limited precision can result in round of errors. Use this as extra safe guard
    if(fabs(ds2 -  ((group_eps + node_eps ) * (group_eps + node_eps) )) < 10e-04) return true;
  #endif

   if (ds2     <= fabs(nodeCOM.w))           return true;
   if (fabs(ds2 - fabs(nodeCOM.w)) < 10e-04) return true; //Limited precision can result in round of errors. Use this as extra safe guard
   
   return false;
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
  float3 dr = {fabs(boxCenter.x - nodeCenter.x) - (boxSize.x + nodeSize.x),
               fabs(boxCenter.y - nodeCenter.y) - (boxSize.y + nodeSize.y),
               fabs(boxCenter.z - nodeCenter.z) - (boxSize.z + nodeSize.z)};

  dr.x += fabs(dr.x); dr.x *= 0.5f;
  dr.y += fabs(dr.y); dr.y *= 0.5f;
  dr.z += fabs(dr.z); dr.z *= 0.5f;

  float ds2    = dr.x*dr.x + dr.y*dr.y + dr.z*dr.z;
  
  #ifdef INDSOFT
    if(ds2 <=      ((group_eps + node_eps ) * (group_eps + node_eps) ))           return true;
    if(fabs(ds2 -  ((group_eps + node_eps ) * (group_eps + node_eps) )) < 10e-04) return true;    
  #endif  

  if (ds2     <= fabs(nodeCenter.w))           return true;
  if (fabs(ds2 - fabs(nodeCenter.w)) < 10e-04) return true; //Limited precision can result in round of errors. Use this as extra safe guard

  return false;
     
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

    int particleOffset     = 1;
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
      
      dataBuffer[multiPoleOffset++]  = multipole[i*3 + 0];
      dataBuffer[multiPoleOffset++]  = multipole[i*3 + 1];
      dataBuffer[multiPoleOffset++]  = multipole[i*3 + 2];        
    }
    
    //Start the tree-walk
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
          child    =    childinfo & 0x0FFFFFFF;                   //Index to the first child of the node
          nchild   = (((childinfo & 0xF0000000) >> 28)) ;         //The number of children this node has              
        }
        else
        {
          //Leaf
          child   =    childinfo & BODYMASK;                     //the first body in the leaf
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
        

        uint temp = 0;  //A node that is not split and is not a leaf will get childinfo 0
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
              
        
        
        //Add the node data to the appropriate arrays and modify the node reference
        //start ofset for its children, should be nodeCount at start of this level +numberofnodes on this level
        //plus a counter that counts the number of childs of the nodes we have tested

        //New childoffset:
        union{int i; float f;} itof; //__int_as_float
        itof.i         = temp;            
        float tempConv = itof.f;
        
        //Add node properties and update references
        real4 nodeSizeInfoTemp  = nodeSizeInfo[node];
        nodeSizeInfoTemp.w      = tempConv;             //Replace child reference

        dataBuffer[nodeSizeOffset++]   = nodeSizeInfoTemp;
        dataBuffer[nodeCenterOffset++] = nodeCenterInfo[node];
        dataBuffer[multiPoleOffset++]  = multipole[node*3 + 0];
        dataBuffer[multiPoleOffset++]  = multipole[node*3 + 1];
        dataBuffer[multiPoleOffset++]  = multipole[node*3 + 2];       
        
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
//    cout << "Mass sum: " << massSum  << endl;
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
    
/*    fprintf(stderr, "Count found: %d particles and %d nodes. Boxsize: (%f %f %f ) BoxCenter: (%f %f %f)\n", 
                    particles, nodes, boxSize.x ,boxSize.y, boxSize.z, boxCenter.x, boxCenter.y, boxCenter.z );  */ 
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
