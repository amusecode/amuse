#include <iomanip>
#include <omp.h>

#include <types.h>
#include <my_errors.h>
#include <kernel.h>
#include <functions.h>
#include <utilis.h>

using namespace std;

HostError CudaInit(unsigned int *M, int NGPU, int rank, string gpu_name, const bool setdev, vector<unsigned int> &dev, string path){

	cudaDeviceProp *properties;
	int count;
	DeviceSafeCall(cudaGetDeviceCount(&count));
	properties = new cudaDeviceProp [count];
	ofstream hlog;
	unsigned int *to_use = new unsigned int [NGPU];
	int k = 0;
   string temp;
   char*  output_name;

	if((unsigned)NGPU != dev.size() && setdev) return HNoNumber;

	if(rank == 0)
	   temp = path + "HiGPUslog.dat";
      output_name = to_char(temp);
      hlog.open(output_name, ios::app);

	if(count < NGPU || count <= 0)
		return HNoGpus;

	for(int i = 0; i < count; i++){
		DeviceSafeCall(cudaGetDeviceProperties(&properties[i], i));
		if(rank == 0)
			hlog<<" Available : "<<properties[i].name<<" as device : "<<i<<endl;
	}

	if(rank == 0)
		hlog<<"============================================="<<endl;

    if(!setdev){

	for(int i = 0; i < count; i++){
	    if(gpu_name.length()>0 && to_string(properties[i].name) != gpu_name) { 
		continue;
	    } else {
		to_use[k] = i;
		k++;
		if(k >= NGPU) {break;} 
	    }
	}

   }else{

      for(unsigned int i = 0; i < dev.size(); i++){
          if(to_string(properties[dev[i]].name) != gpu_name)
             continue;
          else{
             to_use[k] = dev[i];
             k++;
          }
      }
   }

  if(k<NGPU)
     return HNoGpus;


	if(rank==0) {
		for(int i = 0; i < NGPU; i++)
			hlog<<" Using : "<<properties[to_use[i]].name<<" (device "<<to_use[i]<<")"<<endl;
	}

	if(rank == 0){
		if(properties[to_use[0]].major == 2)
			*M = properties[to_use[0]].multiProcessorCount * 1536;
		else if(properties[to_use[0]].major == 3)
			*M = properties[to_use[0]].multiProcessorCount * 2048;
		else if(properties[to_use[0]].major == 1){
			if(properties[to_use[0]].minor == 3)
				*M = properties[to_use[0]].multiProcessorCount * 1024;
			else
				return HNoDouble;
		}
		cout<<" Maximum number of parallel threads on the gpu : "<<*M<<endl;
	}

	MPISafeCall(MPI_Bcast(M, 1, MPI_INT, 0, MPI_COMM_WORLD));

	dev.resize(NGPU);

	for(int i = 0; i < NGPU; i++)
		dev[i] = to_use[i];


	delete [] properties;
	delete [] to_use;

	if(rank == 0)
		hlog.close();



	return HNoError;
}
