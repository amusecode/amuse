#include <iostream>
#include <stdlib.h>
#include <vector>
#include <fstream>

using namespace std;


#include "octree.h"

void read_dumbp_file_parallel(vector<real4> &bodyPositions, vector<real4> &bodyVelocities,  vector<int> &bodiesIDs,  float eps2,
                     string fileName, int rank, int procs, int &NTotal2, int &NFirst, int &NSecond, int &NThird, octree *tree)  
{
  //Process 0 does the file reading and sends the data
  //to the other processes
  
  //Now we have different types of files, try to determine which one is used
  /*****
  If individual softening is on there is only one option:
  Header is formatted as follows: 
  N     #       #       #
  so read the first number and compute how particles should be distributed
  
  If individual softening is NOT enabled, i can be anything, but for ease I assume standard dumbp files:
  no Header
  ID mass x y z vx vy vz
  now the next step is risky, we assume mass adds up to 1, so number of particles will be : 1 / mass
  use this as initial particle distribution
  
  */
  
  
  char fullFileName[256];
  sprintf(fullFileName, "%s", fileName.c_str());

  cout << "Trying to read file: " << fullFileName << endl;

  ifstream inputFile(fullFileName, ios::in);

  if(!inputFile.is_open())
  {
    cout << "Can't open input file \n";
    exit(0);
  }
  
  int NTotal;
  int idummy;
  real4 positions;
  real4 velocity;

  #ifndef INDSOFT
     inputFile >> idummy >> positions.w;
     inputFile.seekg(0, ios::beg); //Reset file pointer
     NTotal = 1 / positions.w;
  #else
     //Read the Ntotal from the file header
     inputFile >> NTotal >> NFirst >> NSecond >> NThird;
  #endif
  
  
  
  //Rough divide
  uint perProc = NTotal / procs;
  bodyPositions.reserve(perProc+10);
  bodyVelocities.reserve(perProc+10);
  bodiesIDs.reserve(perProc+10);
  perProc -= 1;

  //Start reading
  int particleCount = 0;
  int procCntr = 1;
  while(!inputFile.eof()) {
    
    inputFile >> idummy
              >> positions.w >> positions.x >> positions.y >> positions.z
              >> velocity.x >> velocity.y >> velocity.z;    
    
    #ifndef INDSOFT
      velocity.w = sqrt(eps2);
    #else
      inputFile >> velocity.w; //Read the softening from the input file
    #endif
    
    bodyPositions.push_back(positions);
    bodyVelocities.push_back(velocity);
    
    #ifndef INDSOFT    
      idummy = particleCount;
    #endif
    
    bodiesIDs.push_back(idummy);  
    
    particleCount++;
  
  
    if(bodyPositions.size() > perProc && procCntr != procs)
    { 
      tree->ICSend(procCntr,  &bodyPositions[0], &bodyVelocities[0],  &bodiesIDs[0], bodyPositions.size());
      procCntr++;
      
      bodyPositions.clear();
      bodyVelocities.clear();
      bodiesIDs.clear();
    }
  }//end while
  
  inputFile.close();
  
  //Clear the last one since its double
  bodyPositions.resize(bodyPositions.size()-1);  
  NTotal2 = particleCount-1;
  
  cerr << "NTotal: " << NTotal << "\tper proc: " << perProc << "\tFor ourself:" << bodiesIDs.size() << endl;
}


int read_snapshot(vector<real4> &bodyPositions, vector<real4> &bodyVelocities, vector<int> &bodiesIDs, float eps2, string fileName) {
  int N;

  double massTemp = 0;

  bodyPositions.clear();


  char fullFileName[256];
  sprintf(fullFileName, "%s", fileName.c_str());

  cerr << "Trying to read file: " << fullFileName << endl;

  ifstream inputFile(fullFileName, ios::in);

  if(!inputFile.is_open())
  {
    cout << "Can't open input file \n";
    exit(0);
  }


  // reading header
  inputFile >> N;
  bodyPositions.resize(N);
  bodyVelocities.resize(N);
  bodiesIDs.resize(N);

  int ndummy;
  float tnow;

  inputFile >> ndummy;
  inputFile >> tnow;

  // reading data
  for(int i = 0; i < N; i++) {
    inputFile >> bodyPositions[i].w;

    massTemp += bodyPositions[i].w;

  }
  for(int i = 0;i < N; i++) {
    inputFile >> bodyPositions[i].x
	>> bodyPositions[i].y
	>> bodyPositions[i].z;
  }
  for(int i = 0;i < N; i++) {
    inputFile >> bodyVelocities[i].x
	>> bodyVelocities[i].y
	>> bodyVelocities[i].z;
    bodyVelocities[i].w = sqrt(eps2);
    bodiesIDs[i] = i;
  }

  cerr << "Number of particles read: " << N << "\tMass: " << massTemp << endl;

  return 0;
}




void read_dumbp_file(vector<real4> &bodyPositions, vector<real4> &bodyVelocities, vector<int> &bodiesIDs, float eps2, string fileName, int offset) {
  bodyPositions.clear();


  char fullFileName[256];
  sprintf(fullFileName, "%s", fileName.c_str());

  cout << "Trying to read file: " << fullFileName << endl;

  ifstream inputFile(fullFileName, ios::in);

  if(!inputFile.is_open())
  {
    cout << "Can't open input file \n";
    exit(0);
  }

  int idummy;
  real4 positions;
  real4 velocity;

  int cntr = 0;
  while(!inputFile.eof()) {
    inputFile >> idummy
        >> positions.w >> positions.x >> positions.y >> positions.z
        >> velocity.x >> velocity.y >> velocity.z;

    velocity.w = sqrt(eps2);

    bodyPositions.push_back(positions);
    bodyVelocities.push_back(velocity);
//    bodiesIDs.push_back(idummy);
    bodiesIDs.push_back(cntr++ + offset);
    
    
    assert(cntr < 1024*1024*50);
  }

  inputFile.close();
  int n = bodyPositions.size();
  bodyPositions.resize(n-1);

  fprintf(stdout, "read %d bodies from dump file \n", n-1);
};


void read_dumbp_bd_file(vector<real4> &bodyPositions, vector<real4> &bodyVelocities, vector<int> &bodiesIDs, float eps2, string fileName) {
  bodyPositions.clear();


  char fullFileName[256];
  sprintf(fullFileName, "%s", fileName.c_str());

  cout << "Trying to read file: " << fullFileName << endl;

  ifstream inputFile(fullFileName, ios::in);

  if(!inputFile.is_open())
  {
    cout << "Can't open input file \n";
    exit(0);
  }

  int idummy;
  real4 positions;
  real4 velocity;

  inputFile >> idummy >> idummy >> idummy >> idummy;

  int cntr = 0;
  while(!inputFile.eof()) {
    assert(cntr < 1024*1024*50);
    
    inputFile >> idummy
        >> positions.w >> positions.x >> positions.y >> positions.z
        >> velocity.x >> velocity.y >> velocity.z >> velocity.w;

    //cout << cntr << "\t" <<  positions.w << "\t" << velocity.z << endl;

    velocity.w = sqrt(eps2);
//     velocity.w = 0.0005;
    bodyPositions.push_back(positions);
    bodyVelocities.push_back(velocity);
//    bodiesIDs.push_back(idummy);
    bodiesIDs.push_back(cntr++);
  }

  inputFile.close();
  int n = bodyPositions.size();
  bodyPositions.resize(n-1);

  fprintf(stdout, "read %d bodies from dump file \n", n-1);
};

void read_bd_dumbp_file(vector<real4> &bodyPositions, vector<real4> &bodyVelocities,
                     vector<int> &bodiesIDs, float eps2, string fileName,
                     int &NTotal, int &NStars, int &NGas, int &NBH) {
  bodyPositions.clear();


  char fullFileName[256];
  sprintf(fullFileName, "%s", fileName.c_str());

  cout << "Trying to read file: " << fullFileName << endl;

  ifstream inputFile(fullFileName, ios::in);

  if(!inputFile.is_open())
  {
    cout << "Can't open input file \n";
    exit(0);
  }


  int idummy;
  inputFile >> NTotal >> NGas >> NStars >> NBH;
//   inputFile >> NTotal >> NGas >> NStars;

  real4 positions;
  real4 velocity;

   int cntr = 0;
  while(!inputFile.eof()) {
    inputFile >> idummy
        >> positions.w >> positions.x >> positions.y >> positions.z
        >> velocity.x >> velocity.y >> velocity.z >> velocity.w;

    bodyPositions.push_back(positions);
    bodyVelocities.push_back(velocity);
    bodiesIDs.push_back(idummy);

    cntr++;
    assert(cntr++ < 1024*1024*50);
  }

  inputFile.close();
  int n = bodyPositions.size();
  bodyPositions.resize(n-1);

  //Softening of one of the last stars

  fprintf(stdout, "read %d bodies from dump file \n", n-1);
};


//Function that mangles the filename into the correct format 
//incase we use more than 1 process
void read_dumbp_file(vector<real4> &bodyPositions, vector<real4> &bodyVelocities,  vector<int> &bodiesIDs,  float eps2,
                     string fileName, int rank, int procs, int &NTotal, int &NFirst, int &NSecond, int &NThird)  
{

  
  char fullFileName[256];
  if(procs > 0)
    sprintf(fullFileName, "%sMP%d-%d", fileName.c_str(), procs, rank);
//  else
//    sprintf(fullFileName, "%s", fileName.c_str());
  
  #ifdef INDSOFT
    //read_bd_dumbp_file(bodyPositions, bodyVelocities, bodiesIDs, eps2, 
    //                   fileName.c_str(), NTotal, NStars, NGas, NBH);
    read_bd_dumbp_file(bodyPositions, bodyVelocities, bodiesIDs, eps2, 
                       fullFileName, NTotal, NFirst, NSecond, NThird);    
  #else
    int offset = rank*100000;
   read_dumbp_file(bodyPositions, bodyVelocities, bodiesIDs, eps2, fullFileName, offset);  
  #endif  
  
 
};




long long my_dev::base_mem::currentMemUsage;
long long my_dev::base_mem::maxMemUsage;

int main(int argc, char** argv)
{

  vector<real4> bodyPositions;
  vector<real4> bodyVelocities;
  vector<int>   bodyIDs;

  float eps      = 0.05;
  float theta    = 0.75;
  float timeStep = 1.0 / 16.0;
  int  tEnd      = 1000;
  int devID      = 0;

  string fileName       =  "";
  string logFileName    = "gpuLog.log";
  string snapshotFile   = "snapshot_";
  int snapshotIter      = -1;
  float  killDistance   = -1.0;
  float  remoDistance   = -1.0;
  int    snapShotAdd    =  0;


   if (argc <= 1) {
    cout << "Arguments: (in between [] are optional \n";
    cout << "\t-inputFile (dumbp format) \n";
    cout << "\t-[gpulogfile  (gpuLog.log is default)] \n";
    cout << "\t-[device id (0 is default, tries any other device if 0 fails)]\n";
    cout << "\t-[Timestep value  (1/16 is default)]\n";
    cout << "\t-[N-body end time (1000 is default)]\n";
    cout << "\t-[eps  (Will be squared) (0.05 is default)]\n";
    cout << "\t-[theta (0.75 is fefault)]\n";
    cout << "\t-[snapshot base filename (N-body time is appended in 000000 format) ('snapshot_' is default]\n";
    cout << "\t-[snapshot iteration (Nbody time)  (-1 to disable, is also default)]\n";
    cout << "\t-[Killlll distance  (-1 to disable, is also default)]\n";
    cout << "\t-[Particle removal distance  (-1 to disable, is also default)]\n";
    cout << "\t-[Value to add to the snapshot value (0 is default)] \n";

    exit(0);
  }

  if (argc > 1) {
    fileName = string(argv[1]);
  }
  if (argc > 2) {
    logFileName = string(argv[2]);
  }
  if (argc > 3) {
    devID = atoi(argv[3]);
  }
  if (argc > 4) {
    timeStep = atof(argv[4]);
  }
  if (argc > 5) {
    tEnd = atoi(argv[5]);
  }
  if (argc > 6) {
    eps = atof(argv[6]);
  }
  if (argc > 7) {
    theta = atof(argv[7]);
  }
  if(argc > 8)
  {
    snapshotFile = string(argv[8]);
  }
  if(argc > 9)
  {
    snapshotIter = atoi(argv[9]);
  }
  if (argc > 10) {
    killDistance = atof(argv[10]);
  }
  if (argc > 11) {
    remoDistance = atof(argv[11]);
  }
  if(argc > 12)
  {
    snapShotAdd = atoi(argv[12]);
  }
  
  cout << "Used settings: \n";
  cout << "Theta: \t\t"             << theta        << "\t\teps: \t\t"          << eps << endl;
  cout << "Timestep: \t"          << timeStep     << "\t\ttEnd: \t\t"         << tEnd << endl;
  cout << "snapshotFile: \t"      << snapshotFile << "\tsnapshotIter: \t" << snapshotIter << endl;
  cout << "Input file: \t"        << fileName     << "\t\tdevID: \t\t"        << devID << endl;
  cout << "Kill distance: \t"      << killDistance     << "\t\tRemove dist: \t"   << remoDistance << endl;
  cout << "Snapshot Addition: \t"  << snapShotAdd << endl;


  int NTotal, NFirst, NSecond, NThird;
  NTotal = NFirst = NSecond = NThird = 0;

  //Creat the octree class and set the properties
//   octree *tree = new octree(devID, theta, eps, snapshotFile, snapshotIter, timeStep, tEnd);      //device, theta, eps
  octree *tree = new octree(devID, theta, eps, snapshotFile, snapshotIter,  timeStep, tEnd, killDistance, remoDistance, snapShotAdd);
                            
                            
  //Get parallel processing information  
  int procId = tree->mpiGetRank();
  int nProcs = tree->mpiGetNProcs();
  
  //Used for profiler
  char *gpu_prof_log;
  gpu_prof_log=getenv("CUDA_PROFILE_LOG");
  if(gpu_prof_log){
    char tmp[50];
    sprintf(tmp,"process%d_%s",procId,gpu_prof_log);
    setenv("CUDA_PROFILE_LOG",tmp,1);
  }
      
  if(nProcs > 1)
  {
    logFileName.append("-");
    
    char buff[16];
    sprintf(buff,"%d", procId);
    logFileName.append(buff);
  }
  
  ofstream logFile(logFileName.c_str());
    
  tree->set_context(logFile, false); //Do logging to file and enable timing (false = enabled)
  
  if(procId == 0)
    read_dumbp_file_parallel(bodyPositions, bodyVelocities, bodyIDs, eps, fileName, procId, nProcs, NTotal, NFirst, NSecond, NThird, tree);
  else
    tree->ICRecv(0, bodyPositions, bodyVelocities,  bodyIDs);
  
  
  //Set the properties of the data set, it only is really used by process 0, which does the 
  //actual file I/O  
  tree->setDataSetProperties(NTotal, NFirst, NSecond, NThird);
  
  if(procId == 0)  
    cout << "Dataset particle information:\t" << NTotal << " " << NFirst << " " << NSecond << " " << NThird << std::endl;
  
  //Sanity check for standard plummer spheres
  double mass = 0, totalMass;
  for(unsigned int i=0; i < bodyPositions.size(); i++)
  {
    mass += bodyPositions[i].w;
  }
  
  tree->load_kernels();

  MPI_Reduce(&mass,&totalMass,1, MPI_DOUBLE, MPI_SUM,0, MPI_COMM_WORLD);
  
  if(procId == 0)   cerr << "Combined Mass: "  << totalMass << "\tNTotal: " << NTotal << std::endl;
  
  
  //Decide what kind of file type to read....
  //read_snapshot(bodyPositions, bodyVelocities, bodyIDs, eps, fileName);
//   read_dumbp_file(bodyPositions, bodyVelocities, bodyIDs, eps, fileName);
//   read_dumbp_bd_file(bodyPositions, bodyVelocities, bodyIDs, eps, fileName);
  //read_bd_file(bodyPositions, bodyVelocities, eps, fileName);
  //   read_bd_dumbp_file(bodyPositions, bodyVelocities, bodyIDs, eps, fileName, NTotal, NStars, NGas, NBH);
  
//   tree->write_dumbp_snapshot_parallel(&bodyPositions[0], &bodyVelocities[0], &bodyIDs[0], bodyPositions.size(), "test.ascii");

  
  //Parallel code version
  //read_dumbp_file(bodyPositions, bodyVelocities, bodyIDs, eps, fileName, procId, nProcs, NTotal, NStars, NGas, NBH);
  
  
  tree->createORB();
  
  //First distribute the initial particle distribution
  //over all available processes
  if(tree->nProcs > 1)
  {
    tree->createDistribution(&bodyPositions[0], bodyPositions.size());  
  }

  //Print the domain division
  if(tree->nProcs > 1)
  {
    if(tree->procId == 0)
      for(int i = 0;i< tree->nProcs;i++)     
      {      
        cerr << i << " " << tree->xlow[i].x << " " << tree->xlow[i].y << " " << tree->xlow[i].z << " " 
           << tree->xhigh[i].x << " " << tree->xhigh[i].y << " " << tree->xhigh[i].z <<endl;
      }
  }

  tree->mpiSync();    
  

  printf("Starting! \n");
  

  double t0 = tree->get_time();

  tree->localTree.setN(bodyPositions.size());
  tree->allocateParticleMemory(tree->localTree);

  //Load data onto the device
  for(uint i=0; i < bodyPositions.size(); i++)
  {
    tree->localTree.bodies_pos[i] = bodyPositions[i];
    tree->localTree.bodies_vel[i] = bodyVelocities[i];
    tree->localTree.bodies_ids[i] = bodyIDs[i];

    tree->localTree.bodies_Ppos[i] = bodyPositions[i];
    tree->localTree.bodies_Pvel[i] = bodyVelocities[i];
  }

  tree->localTree.bodies_pos.h2d();
  tree->localTree.bodies_vel.h2d();
  tree->localTree.bodies_Ppos.h2d();
  tree->localTree.bodies_Pvel.h2d();
  tree->localTree.bodies_ids.h2d();
   
  //Distribute the particles so each process has particles
  //assigned to his domain
  if(nProcs > 1)
  {    
    double ttemp = tree->get_time();
    cout << "Before exchange tree has : " << tree->localTree.n << " particles \n" ;
    while(tree->exchange_particles_with_overflow_check(tree->localTree));
    cout << "After exchange tree has : " << tree->localTree.n << " particles \n" ;    
    //Send the new and old particles to the device
    tree->localTree.bodies_pos.h2d();
    tree->localTree.bodies_vel.h2d();
    tree->localTree.bodies_ids.h2d();
    tree->localTree.bodies_acc0.h2d();
    tree->localTree.bodies_acc1.h2d();
    tree->localTree.bodies_time.h2d();
    
    //This is only required the first time since we have no predict call before we build the tree
    //Every next call this is not required since the predict call fills the predicted positions
    tree->localTree.bodies_Ppos.copy(tree->localTree.bodies_pos, tree->localTree.bodies_pos.get_size());
    tree->localTree.bodies_Pvel.copy(tree->localTree.bodies_vel, tree->localTree.bodies_vel.get_size());
    
    printf("Initial exchange Took in total: %lg sec\n", tree->get_time()-ttemp);
  }
  

  //Start construction of the tree
  tree->sort_bodies(tree->localTree);
  tree->build(tree->localTree);
  tree->allocateTreePropMemory(tree->localTree);
  tree->compute_properties(tree->localTree);

  
  //Start the integration
  tree->iterate();

  printf("Finished!!! Took in total: %lg sec\n", tree->get_time()-t0);

  tree->desort_bodies(tree->localTree);

  logFile.close();

  //   write_dumbp_file(&tree->localTree.bodies_pos[0], &tree->localTree.bodies_vel[0], tree->localTree.n, "outputFile.txt") ;
  MPI_Finalize();

  delete tree;
  tree = NULL;

  return 0;
}
