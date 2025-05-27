#ifndef _MY_OCL_H_
#define _MY_OCL_H_

typedef float real;
typedef unsigned int uint;

#define __align__(n) __attribute__((aligned(n)))

#define __builtin_align__(a) __align__(a)

#define __cuda_assign_operators(tag)

struct real4 {
  real x, y, z, w;
};

struct real3 {
  real x, y, z;
};

struct real2 {
  real x, y;
};

struct float2 {
  real x, y;
};
struct float3 {
  real x, y, z;
};

/*DEVICE_BUILTIN*/
struct __builtin_align__(16) float4
{
  float x, y, z, w;
  __cuda_assign_operators(float4)
};

/*DEVICE_BUILTIN*/
struct __builtin_align__(16) double2
{
  double x, y;
  __cuda_assign_operators(double2)
};



//////////////////

struct int4 {
  int x, y, z, w;
};

struct int3 {
  int x, y, z;
};

struct int2 {
  int x, y;
};

////////////////

struct uint4 {
  uint x, y, z, w;
};

struct uint3 {
  uint x, y, z;
};

struct uint2 {
  uint x, y;
};


#include <string>
#include <fstream>

#include <stdio.h>
#include <stdlib.h>
#include <CL/cl.h>
#include <CL/cl_platform.h>
#include <string.h>
#include <cassert>
#include <vector>
#include <cmath>
#include <sys/time.h>


using namespace std;


//Function made by NVIDIA
//////////////////////////////////////////////////////////////////////////////
//! Loads a Program file and prepends the cPreamble to the code.
//!
//! @return the source string if succeeded, 0 otherwise
//! @param cFilename        program filename
//! @param cPreamble        code that is prepended to the loaded file, typically a set of #defines or a header
//! @param szFinalLength    returned length of the code string
//////////////////////////////////////////////////////////////////////////////
inline char* oclLoadProgSource(const char* cFilename, const char* cPreamble, size_t* szFinalLength)
{
    // locals
    FILE* pFileStream = NULL;
    size_t szSourceLength;

    // open the OpenCL source code file
    #ifdef _WIN32   // Windows version
        if(fopen_s(&pFileStream, cFilename, "rb") != 0)
        {
            return NULL;
        }
    #else           // Linux version
        pFileStream = fopen(cFilename, "rb");
        if(pFileStream == 0)
        {
            return NULL;
        }
    #endif
    
    size_t szPreambleLength = strlen(cPreamble);

    // get the length of the source code
    fseek(pFileStream, 0, SEEK_END);
    szSourceLength = ftell(pFileStream);
    fseek(pFileStream, 0, SEEK_SET);

    // allocate a buffer for the source code string and read it in
    char* cSourceString = (char *)malloc(szSourceLength + szPreambleLength + 1);
    memcpy(cSourceString, cPreamble, szPreambleLength);
    if (fread((cSourceString) + szPreambleLength, szSourceLength, 1, pFileStream) != 1)
    {
        fclose(pFileStream);
        free(cSourceString);
        return 0;
    }

    // close the file and return the total length of the combined (preamble + source) string
    fclose(pFileStream);
    if(szFinalLength != 0)
    {
        *szFinalLength = szSourceLength + szPreambleLength;
    }
    cSourceString[szSourceLength + szPreambleLength] = '\0';

    return cSourceString;
}


inline const char* oclPrintError(cl_int err)
{
    switch (err) {
        case CL_SUCCESS:                           return "Success!\n";
        case CL_DEVICE_NOT_FOUND:                  return "Device not found.\n";
        case CL_DEVICE_NOT_AVAILABLE:              return "Device not available\n";
        case CL_COMPILER_NOT_AVAILABLE:            return "Compiler not available\n";
        case CL_MEM_OBJECT_ALLOCATION_FAILURE:     return "Memory object allocation failure\n";
        case CL_OUT_OF_RESOURCES:                  return "Out of resources\n";
        case CL_OUT_OF_HOST_MEMORY:                return "Out of host memory\n";
        case CL_PROFILING_INFO_NOT_AVAILABLE:      return "Profiling information not available\n";
        case CL_MEM_COPY_OVERLAP:                  return "Memory copy overlap\n";
        case CL_IMAGE_FORMAT_MISMATCH:             return "Image format mismatch\n";
        case CL_IMAGE_FORMAT_NOT_SUPPORTED:        return "Image format not supported\n";
        case CL_BUILD_PROGRAM_FAILURE:             return "Program build failure\n";
        case CL_MAP_FAILURE:                       return "Map failure\n";
        case CL_INVALID_VALUE:                     return "Invalid value\n";
        case CL_INVALID_DEVICE_TYPE:               return "Invalid device type\n";
        case CL_INVALID_PLATFORM:                  return "Invalid platform\n";
        case CL_INVALID_DEVICE:                    return "Invalid device\n";
        case CL_INVALID_CONTEXT:                   return "Invalid context\n";
        case CL_INVALID_QUEUE_PROPERTIES:          return "Invalid queue properties\n";
        case CL_INVALID_COMMAND_QUEUE:             return "Invalid command queue\n";
        case CL_INVALID_HOST_PTR:                  return "Invalid host pointer\n";
        case CL_INVALID_MEM_OBJECT:                return "Invalid memory object\n";
        case CL_INVALID_IMAGE_FORMAT_DESCRIPTOR:   return "Invalid image format descriptor\n";
        case CL_INVALID_IMAGE_SIZE:                return "Invalid image size\n";
        case CL_INVALID_SAMPLER:                   return "Invalid sampler\n";
        case CL_INVALID_BINARY:                    return "Invalid binary\n";
        case CL_INVALID_BUILD_OPTIONS:             return "Invalid build options\n";
        case CL_INVALID_PROGRAM:                   return "Invalid program\n";
        case CL_INVALID_PROGRAM_EXECUTABLE:        return "Invalid program executable\n";
        case CL_INVALID_KERNEL_NAME:               return "Invalid kernel name\n";
        case CL_INVALID_KERNEL_DEFINITION:         return "Invalid kernel definition\n";
        case CL_INVALID_KERNEL:                    return "Invalid kernel\n";
        case CL_INVALID_ARG_INDEX:                 return "Invalid argument index\n";
        case CL_INVALID_ARG_VALUE:                 return "Invalid argument value\n";
        case CL_INVALID_ARG_SIZE:                  return "Invalid argument size\n";
        case CL_INVALID_KERNEL_ARGS:               return "Invalid kernel arguments\n";
        case CL_INVALID_WORK_DIMENSION:            return "Invalid work dimension\n";
        case CL_INVALID_WORK_GROUP_SIZE:           return "Invalid work group size\n";
        case CL_INVALID_WORK_ITEM_SIZE:            return "Invalid work item size\n";
        case CL_INVALID_GLOBAL_OFFSET:             return "Invalid global offset\n";
        case CL_INVALID_EVENT_WAIT_LIST:           return "Invalid event wait list\n";
        case CL_INVALID_EVENT:                     return "Invalid event\n";
        case CL_INVALID_OPERATION:                 return "Invalid operation\n";
        case CL_INVALID_GL_OBJECT:                 return "Invalid OpenGL object\n";
        case CL_INVALID_BUFFER_SIZE:               return "Invalid buffer size\n";
        case CL_INVALID_MIP_LEVEL:                 return "Invalid mip-map level\n";
        default:return "Unknown\n";
    }
}


inline void __oclsafeCall(cl_int err, const char *file, const int line) {
  if(CL_SUCCESS != err) {
    fprintf(stderr, "oclSafeCall() Runtime API error in file <%s>, line %i : %s.\n",
	    file, line, oclPrintError(err));
    exit(-1);
  }
}
#define oclSafeCall(err)    __oclsafeCall(err, __FILE__, __LINE__)
#define oclCheckError(err)  __oclsafeCall(err, __FILE__, __LINE__)


namespace my_dev {

  class context {
  protected:
    size_t dev;

    cl_device_type device_type;
    cl_context hContext;

    cl_uint ciDeviceCount;   
    cl_device_id *hDevices;

    cl_command_queue hCmdQueue;

    cl_int ciErrNum;

    bool hContext_flag;
    bool hCmdQueue_flag;
    
    //Logging tools
    bool logfile_flag;
    bool disable_timing;
    
    ofstream *logFile;
    
    int logID;  //Unique ID to every log line
    
    
    //Events:
    cl_event start, stop;
    struct timeval Tvalue;
    struct timezone dummy;
    
  public:

    context() {
      hDevices = NULL;
      hContext_flag = false;
      hCmdQueue_flag = false;
      disable_timing    = false;
      
    }
    ~context() {      
      if (hDevices != 0) free(hDevices);
      if (hCmdQueue_flag) clReleaseCommandQueue(hCmdQueue);
      if (hContext_flag ) clReleaseContext     (hContext );
    }
    
    
    
    int create(std::ofstream &log, bool disableTiming = false, cl_device_type device_type = CL_DEVICE_TYPE_GPU)
    {
      logfile_flag      = true;
      logFile = &log;
      logID = 0;
      return create(disable_timing, device_type);     

    }
     
    int create(bool disableTiming = false, cl_device_type device_type = CL_DEVICE_TYPE_GPU) {
      assert(!hContext_flag);
      
      disable_timing    = disableTiming;

            
      printf("Creating OpenCL context \n");
      
      this->device_type = device_type;
      hContext = clCreateContextFromType(0, device_type, NULL, NULL, &ciErrNum);
      oclCheckError(ciErrNum);
      hContext_flag = true;
   
      //Get the number of devices and the deviceIDs and descriptions
      oclSafeCall(clGetDeviceIDs (NULL, this->device_type, 0, NULL, &ciDeviceCount));
     
      hDevices = (cl_device_id*) malloc(sizeof(cl_device_id)*ciDeviceCount);
      oclSafeCall(clGetDeviceIDs (NULL, this->device_type, ciDeviceCount, hDevices, &ciDeviceCount));

      printf("Found %d suitable devices: \n",ciDeviceCount);
      for(unsigned int i=0; i < ciDeviceCount; i++)
      {
        char device_string[1024];
        clGetDeviceInfo(hDevices[i], CL_DEVICE_NAME, sizeof(device_string), &device_string, NULL);
        printf(" %d: %s\n",i, device_string);
      }

      return ciDeviceCount;
    }
    
    void createQueue(size_t dev = 0) {
      assert(!hCmdQueue_flag);
      assert(hContext_flag);
      this->dev = dev;
      assert(dev < ciDeviceCount);
      hCmdQueue = clCreateCommandQueue(hContext, hDevices[dev],  0, &ciErrNum);
      oclCheckError(ciErrNum);
      hCmdQueue_flag = true;
      
      printf("Using device: %d\n", (int)dev);
      
//      printf("Using device: %ld dev[0]: %ld, dev[1]: %ld\n", (long int)hDevices[dev], (long int)hDevices[0], (long int)hDevices[1]);
//       cl_device_id device;  
//       ciErrNum = clGetCommandQueueInfo(hCmdQueue, CL_QUEUE_DEVICE, sizeof(cl_device_id), &device, NULL);
//       printf("Na de call gebruiken we dev: %d \n", device);
    }
    
    
    void startTiming(int stream=0)
    {
      if(disable_timing) return;
      
      //Save the current time
      gettimeofday(&Tvalue,&dummy);     
      
      //Events work different in openCL than in CUDA
      //So for now use a simple wallclock time timer
      
//       CU_SAFE_CALL(cuEventCreate(&start, CU_EVENT_DEFAULT));  
//       CU_SAFE_CALL(cuEventCreate(&stop, CU_EVENT_DEFAULT));
//       CU_SAFE_CALL(cuEventRecord(start, stream));
    }
    
    //Text and ID to be printed with the log message on screen / in the file
    void stopTiming(const char *text, int type = -1, int stream=0)
    {
      if(disable_timing) return;
      
      struct timeval Tvalue2;
      struct timezone dummy2;
      gettimeofday(&Tvalue2,&dummy2);

      //Calculate the elapsed time
      double startTime =  ((double) Tvalue.tv_sec +1.e-6*((double) Tvalue.tv_usec));
      double endTime =  ((double) Tvalue2.tv_sec +1.e-6*((double) Tvalue2.tv_usec));


      double time = endTime-startTime; 
      time *= 1000;
      
      
/*      CU_SAFE_CALL(cuEventRecord(stop, stream));
      CU_SAFE_CALL(cuEventSynchronize(stop));
      float time;
      CU_SAFE_CALL(cuEventElapsedTime(&time, start, stop));      
      CU_SAFE_CALL(cuEventDestroy(start));
      CU_SAFE_CALL(cuEventDestroy(stop));   */   
      
      printf("%s took:\t%f\t millisecond\n", text, time);
      
      if(logfile_flag)
      {
        (*logFile) << logID++ << "\t"  << type << "\t" << text << "\t" << time << endl;
      }
    }    
    
    /////////////
    
    cl_context&       get_context()       {return hContext;}
    cl_command_queue& get_command_queue() {return hCmdQueue;}
    cl_device_id      operator[](int i)   {return hDevices[i];}
    cl_device_id      get_device()        {return hDevices[dev];}
    
    //////////       
    
  };

  ///////////////////////
    
  class base_mem
  {
    public:
    //Memory usage counters
    static long long currentMemUsage;
    static long long maxMemUsage;  
    
    void increaseMemUsage(int bytes)
    {
      currentMemUsage +=  bytes;   
      
      if(currentMemUsage > maxMemUsage)
        maxMemUsage = currentMemUsage;
    }
    
    void decreaseMemUsage(int bytes)
    {
      currentMemUsage -=  bytes;
    }
    
    static void printMemUsage()
    {      
      printf("Current usage: %lld bytes ( %lld MB) \n", currentMemUsage, currentMemUsage / (1024*1024));
      printf("Maximum usage: %lld bytes ( %lld MB) \n", maxMemUsage, maxMemUsage / (1024*1024));     
      
/*      unsigned int free, total;
      cuMemGetInfo(&free, &total); 
      printf("Build-in usage: free: %d bytes ( %d MB , total: %d) \n", free, free / (1024*1024), total / (1024*1024));   */  
      
    }  
    
    static long long getMaxMemUsage()
    {      
      return maxMemUsage;
    }     
    
    
    
  };



  template<class T>
  class dev_mem : base_mem {
  protected:
    cl_context hContext;
    cl_command_queue hCmdQueue;

    int size;
    cl_mem hDeviceMem;
    std::vector<T> host_ptr;
    bool pinned_mem, context_flag, flags;
    bool hDeviceMem_flag;
    
    void ocl_free() {
      assert(context_flag);
      if (hDeviceMem_flag) {
	assert(size > 0);
	oclSafeCall(clReleaseMemObject(hDeviceMem));
      }
    }

  public:
    
    ///////// Constructors

    dev_mem() {
      size = 0;
      pinned_mem = false;
      hDeviceMem_flag = false;
      context_flag = false;
    }

    dev_mem(class context &c) {
      size = 0;
      pinned_mem = false;
      context_flag = false;
      hDeviceMem_flag = false;
      setContext(c);
    }
    
    dev_mem(class context &c, int n, bool zero = false,
	    cl_mem_flags flags = CL_MEM_READ_WRITE, bool pinned = false) {
      context_flag = false;
      hDeviceMem_flag = false;
      setContext(c);
      if (zero) this->ccalloc(n, pinned, flags);
      else      this->cmalloc(n, pinned, flags);
    }
    
    dev_mem(class context &c, std::vector<T> data,
	    cl_mem_flags flags = CL_MEM_READ_WRITE,  bool pinned = false) {
      context_flag    = false;
      hDeviceMem_flag = false;
      setContext(c);
      this->cmalloc(data, pinned, flags);
    }

    //////// Destructor
    
    ~dev_mem() {
      ocl_free();
    }
    
    ///////////

    void setContext(class context &c) {      
      assert(!context_flag);
      this->hContext   = c.get_context();
      this->hCmdQueue  = c.get_command_queue();      
      context_flag = true;
    }

    ///////////

    void cmalloc(int n, bool pinned = false,  cl_mem_flags flags = CL_MEM_READ_WRITE) {
      assert(pinned == false);  //Until we fix the mapped memory in openCL!
      assert(context_flag);
      assert(!hDeviceMem_flag);
      
      this->pinned_mem = pinned;      
      this->flags = flags;
      if (size > 0) ocl_free();
      size = n;
      host_ptr.resize(size);
      cl_int err;
      if (flags & (CL_MEM_USE_HOST_PTR == 0)) hDeviceMem = clCreateBuffer(hContext, flags, size*sizeof(T), NULL, &err);
      else                             	    hDeviceMem = clCreateBuffer(hContext, flags, size*sizeof(T), &host_ptr[0], &err);
      oclCheckError(err);
      hDeviceMem_flag = true;
    }

    void ccalloc(int n, bool pinned = false,  cl_mem_flags flags = CL_MEM_READ_WRITE) {
      assert(pinned == false);  //TODO Until we fix the mapped memory in openCL!
      assert(context_flag);
      assert(!hDeviceMem_flag);
      
      this->pinned_mem = pinned;
      this->flags = flags;
      if (size > 0) ocl_free();
      size = n;
      host_ptr.resize(size);      
      memset(&host_ptr[0], 0, size*sizeof(T));   

      cl_int err;
      if ((flags & CL_MEM_USE_HOST_PTR) == 0) hDeviceMem = clCreateBuffer(hContext, flags, size*sizeof(T), NULL, &err);
      else                             	    hDeviceMem = clCreateBuffer(hContext, flags, size*sizeof(T), &host_ptr[0], &err);
      oclCheckError(err);
      h2d();
      hDeviceMem_flag = true;
    }

    void cmalloc(std::vector<T> &data, bool pinned = false, cl_mem_flags flags = CL_MEM_READ_WRITE) {
      assert(pinned == false);  //Until we fix the mapped memory in openCL!
      assert(context_flag);
      assert(!hDeviceMem_flag);
      this->flags = flags;
      this->pinned_mem = pinned;
      if (size > 0) ocl_free();
      size = data.size();
      host_ptr = data;
      cl_int err;
      if ((flags & CL_MEM_USE_HOST_PTR) == 0) {
	hDeviceMem = clCreateBuffer(hContext, flags | CL_MEM_COPY_HOST_PTR,
				    size*sizeof(T), 
				    &host_ptr[0], 
				    &err);
      } else {
	hDeviceMem = clCreateBuffer(hContext, flags, size*sizeof(T), &host_ptr[0], &err);
      }
      oclCheckError(err);
      hDeviceMem_flag = true;
    }
    
    void zeroMem()
    {
      assert(context_flag);
      assert(hDeviceMem_flag);
      memset(&host_ptr[0], 0, size*sizeof(T));    //Set the vector to value 0 
      
      //Cant fidn a function to set memory directly so copy from host to device for now
      h2d();
      
    }
    
    ///////////

    //////////////

    void d2h(cl_bool OCL_BLOCKING = CL_TRUE)   {
      assert(context_flag);
      assert(hDeviceMem_flag);
      assert(size > 0);
      if ((flags & CL_MEM_USE_HOST_PTR) == 0) return;
      oclSafeCall(clEnqueueReadBuffer(hCmdQueue, hDeviceMem, OCL_BLOCKING, 0,
				      size*sizeof(T),
				      &host_ptr[0], 0, 0, 0));
    }
    
    void h2d(cl_bool OCL_BLOCKING = CL_TRUE)   {
      assert(context_flag);
      assert(hDeviceMem_flag);
      assert(size > 0);
      if ((flags & CL_MEM_USE_HOST_PTR) == 0) return;
      oclSafeCall(clEnqueueWriteBuffer(hCmdQueue, hDeviceMem, OCL_BLOCKING, 0,
				       size*sizeof(T),
				       &host_ptr[0], 0, 0, 0));
    }
    
    void copy(dev_mem &src_buffer, int n, cl_bool OCL_BLOCKING = CL_TRUE){
      assert(context_flag);
      assert(hDeviceMem_flag);
      if (size != n) {
        ocl_free();
        cmalloc(n, this->pinned_mem, flags);
        size = n;
      }     
       oclSafeCall(clEnqueueCopyBuffer(hCmdQueue,
                                      src_buffer.d(),
                                      hDeviceMem,
                                      0, 0, n*sizeof(T),
                                      0, NULL, NULL));    
      //Copy on the host
      memcpy (((void*) &host_ptr[0]), ((void*) &src_buffer[0]), n*sizeof(T));                                
    }
 
    /////////
    
    T& operator[] (int i){return host_ptr[i];}
    cl_mem& d() {return hDeviceMem;}
    void*   p() {return (void*)&hDeviceMem;}
    int get_size(){return size;}
    
  };     // end of class dev_mem

  ////////////////////


  class kernel {
  protected:
    char *hKernelSource;
    char *hKernelFilename;
    char *hKernelBinary;
    char *hKernelName;

    size_t hKernelLength;
    cl_program hProgram;
    cl_kernel  hKernel;
    cl_context hContext;
    cl_device_id hDeviceId;
    cl_device_id hDeviceIdTEST;
    cl_command_queue hCmdQueue;
    int ciErrNum;
    vector<size_t> hGlobalWork;
    vector<size_t> hLocalWork;

    bool context_flag;
    bool kernel_flag;
    bool program_flag;
    bool work_flag;   

  public:

    kernel() {
      hKernelSource = (char*)malloc(1024);
      hKernelBinary = (char*)malloc(1024);
      hKernelName   = (char*)malloc(256);
      hKernelFilename = (char*)malloc(1024);
      hGlobalWork.clear();
      hLocalWork.clear();

      context_flag = false;
      kernel_flag  = false;
      program_flag = false;
      work_flag    = false;
      
    }
    ~kernel() {
      free(hKernelSource);
      free(hKernelBinary);
      free(hKernelName);
      free(hKernelFilename);
      if (kernel_flag) clReleaseKernel(hKernel);
      if (program_flag) clReleaseProgram(hProgram);
    }

    kernel(class context &c) {
      hKernelSource = (char*)malloc(1024);
      hKernelBinary = (char*)malloc(1024);
      hKernelName   = (char*)malloc(256);
      hKernelFilename = (char*)malloc(1024);
      hGlobalWork.clear();
      hLocalWork.clear();

      context_flag = false;
      kernel_flag  = false;
      program_flag = false;
      work_flag    = false;

      setContext(c);
    }

    ////////////

    void setContext(class context &c) {
      assert(!context_flag);
      this->hContext   = c.get_context();
      this->hCmdQueue  = c.get_command_queue();
      this->hDeviceId  = c.get_device();      
      context_flag = true;
    }

    ////////////
    
    void load_source(const char *kernel_name, const char *subfolder,
                     const char *compilerOptions = "",
                     int maxrregcount = -1,
                     int architecture = 0) {
      assert(context_flag);
      assert(!program_flag);
      
      //TODO!!!!!!!!
      // create & compile program
      sprintf(hKernelFilename, "%s%s", subfolder, kernel_name);
      

      hKernelSource = oclLoadProgSource(hKernelFilename, "", &hKernelLength);
      hProgram = clCreateProgramWithSource(hContext, 1, (const char**)&hKernelSource,
 	//				   &hKernelLength, &ciErrNum);
                                           NULL, &ciErrNum);
                                         
      oclCheckError(ciErrNum);
      
#ifdef MAC
      const char* flags = "-cl-mad-enable -DMAC";
#else
      const char* flags = " -cl-mad-enable";
#endif
      
      //Combine flags and custom compile options
      char finalBuildOptions[1024];
      sprintf(finalBuildOptions, "%s %s ", flags, compilerOptions);
      
      if(maxrregcount >= 0)
      {
        sprintf(finalBuildOptions, "%s %s=%d", finalBuildOptions, "-cl-nv-maxrregcount", maxrregcount);      
      }

      ///the -cl-nv-maxrregcount=n build option when building the kernel (
      //#pragma OPENCL EXTENSION cl_khr_fp64 : enable

      ciErrNum = clBuildProgram(hProgram, 0, NULL, finalBuildOptions, NULL, NULL);
      
      if(ciErrNum != CL_SUCCESS)      
        print_compiler_output();

      program_flag = true;
    }

    void print_compiler_output() {
      printf("Compilation of the source file failed: %s \n", hKernelFilename);
      char buildLog[10240];
      clGetProgramBuildInfo(hProgram, hDeviceIdTEST, CL_PROGRAM_BUILD_LOG, sizeof(buildLog), buildLog, NULL);
      printf("Compiler output: \n %s \n", buildLog);
      exit(0);
    }


    void create(const char *kernel_name) {
      assert(program_flag);
      assert(!kernel_flag);
      sprintf(hKernelName, kernel_name,"");
      hKernel = clCreateKernel(hProgram, hKernelName, &ciErrNum);
      oclCheckError(ciErrNum);
      kernel_flag = true;
    }

    //'size'  is used for dynamic shared memory
    template<class T>
    void set_arg(int arg, void* ptr, int size = 1)  {
      assert(kernel_flag);
      ciErrNum = clSetKernelArg(hKernel, arg, size*sizeof(T), ptr);
      oclCheckError(ciErrNum);
    }
    
    
    void setWork(int items, int n_threads, int blocks = -1)
    {
      //Sets the number of blocks and threads based on the number of items
      //and number of threads per block.
      //TODO see if we can use the new kernel define for thread numbers?      
      vector<size_t> localWork(2), globalWork(2);
    
      int nx, ny;
      
      if(blocks == -1)
      {      
        //Calculate dynamic
        int ng = (items) / n_threads + 1;
        nx = (int)sqrt(ng);
        ny = (ng -1)/nx +  1; 
      }
      else
      {
        //Specified number of blocks and numbers of threads make it a
        //2D grid if nessecary        
        if(blocks > 65536)
        {
          nx = (int)sqrt(blocks);
          ny = (blocks -1)/nx +  1;           
        }
        else
        {
          nx = blocks;
          ny = 1;
        }
      }
    
      globalWork[0] = nx*n_threads;  globalWork[1] = ny*1;
      localWork [0] = n_threads;     localWork[1]  = 1;   
      setWork(globalWork, localWork);
    }    

    void setWork(vector<size_t> global_work, vector<size_t> local_work) {
      assert(kernel_flag);
      assert(global_work.size() == local_work.size());
      hGlobalWork = global_work;
      hLocalWork  = local_work;
    }
    
    void printWorkSize()
    {
      printf("Blocks: (%ld, %ld, %ld) Threads: (%ld, %ld, %ld) \n", 
              hGlobalWork[0], hGlobalWork[1], hGlobalWork[2],
              hLocalWork[0], hLocalWork[1], hLocalWork[2]);             
    }
        

    void execute(vector<size_t> global_work, vector<size_t> local_work, 
		 cl_event* event = NULL) {
      setWork(global_work, local_work);
      cl_uint work_dim = hGlobalWork.size();
      oclSafeCall(clEnqueueNDRangeKernel(hCmdQueue, hKernel, work_dim, 0,
			     &hGlobalWork[0],
			     &hLocalWork[0],
			     0, NULL, event));
    }
    
    void execute(cl_event* event = NULL) {
      assert(kernel_flag);
      cl_uint work_dim = hGlobalWork.size(); 
      oclSafeCall(clEnqueueNDRangeKernel(hCmdQueue, hKernel, work_dim, 0,
			     &hGlobalWork[0],
			     &hLocalWork[0],
			     0, NULL, event));
    }
    
    ////
    
    cl_kernel&  get_kernel() {return hKernel;}
    cl_program& get_program() {return hProgram;}
      
  };

}     // end of namespace my_ocl

#endif // _MY_OCL_H_


