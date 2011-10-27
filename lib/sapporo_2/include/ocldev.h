#ifndef __OCL_H__
#define __OCL_H__

#define __OPENCL_DEV__

#include <vector>
#include <cassert>
#include <iostream>
#include <cmath>
#include <cstring>
#include <cstdio>
#include <cstdlib>

#ifdef __MACOSX__
#include <OpenCL/opencl.h>
#else
#include <CL/cl.h>
#include <CL/cl_platform.h>
#endif

namespace dev {

  //Function made by NVIDIA
  //////////////////////////////////////////////////////////////////////////////
  //! Loads a Program file and prepends the cPreamble to the code.
  //!
  //! @return the source string if succeeded, 0 otherwise
  //! @param cFilename        program filename
  //! @param cPreamble        code that is prepended to the loaded file, typically a set of #defines or a header
  //! @param szFinalLength    returned length of the code string
  //////////////////////////////////////////////////////////////////////////////
  inline char* oclLoadProgSource(const char* cFilename, const char* cPreamble, size_t* szFinalLength) {
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
    
    const size_t szPreambleLength = strlen(cPreamble);
    
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
  
  inline const char* oclPrintError(const cl_int err) {
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
    default:return "Unknown OpenCL error\n";
    }
  }


  inline void __oclsafeCall(const cl_int err, const char *file, const int line) {
    if(CL_SUCCESS != err) {
      fprintf(stderr, "oclSafeCall() Runtime API error in file <%s>, line %i : %s.\n",
	      file, line, oclPrintError(err));
      assert(false); 
    }
  }
#define oclSafeCall(err)    __oclsafeCall(err, __FILE__, __LINE__)
#define oclCheckError(err)  __oclsafeCall(err, __FILE__, __LINE__)

  class context {
  protected:
    size_t                     devId;
    cl_device_type             DeviceType;
    cl_context                 Context;
    cl_uint                    DeviceCount;
    std::vector<cl_device_id>  Devices;
    cl_command_queue           CommandQueue;
    bool                       ContextFlag;
    bool                       CommandQueueFlag;
    
    int logId;
    std::vector<cl_platform_id> PlatformIDs;
    cl_platform_id              PlatformID;

  public:
    context() : ContextFlag(false), CommandQueueFlag(false) {};
    ~context() {
      if (ContextFlag     ) clReleaseContext(Context);
      if (CommandQueueFlag) clReleaseCommandQueue(CommandQueue);
    };
    
    const int getPlatformInfo(std::ostream &s = std::cerr) {
      
      s << "Getting list of OpenCL platforms ...\n";
      
      cl_uint numPlatforms;
      oclSafeCall(clGetPlatformIDs(0, NULL, &numPlatforms));
      assert(numPlatforms > 0);
      PlatformIDs.resize(numPlatforms+1);
      oclSafeCall(clGetPlatformIDs(numPlatforms, &PlatformIDs[0], NULL));
      for (cl_uint dev = 0; dev < numPlatforms; dev++) {
	char platform_string[1024];	
        oclSafeCall(clGetPlatformInfo(PlatformIDs[dev], CL_PLATFORM_NAME, sizeof(platform_string), &platform_string, NULL));
	std::cerr << dev << ": " << platform_string << "\n";
      }

      return numPlatforms;
    }

    const int getDeviceCount(const cl_device_type device_type = CL_DEVICE_TYPE_GPU, const int platform_id = 0) {
      assert(!ContextFlag);
      
      std::cerr << "Getting list of OpenCL devices ...\n";
      
      DeviceType = device_type;
      
      cl_uint numPlatforms;
      oclSafeCall(clGetPlatformIDs(0, NULL, &numPlatforms));
      assert(numPlatforms > 0);
      PlatformIDs.resize(numPlatforms+1);
      oclSafeCall(clGetPlatformIDs(numPlatforms, &PlatformIDs[0], NULL));
      for (cl_uint dev = 0; dev < numPlatforms; dev++) {
	char platform_string[1024];	
        oclSafeCall(clGetPlatformInfo(PlatformIDs[dev], CL_PLATFORM_NAME, sizeof(platform_string), &platform_string, NULL));
	std::cerr << " " << dev << ": " << platform_string << "\n";
      }
      fprintf(stderr, "Using platform %d \n", platform_id);
      PlatformID = PlatformIDs[platform_id];
      
      oclSafeCall(clGetDeviceIDs(PlatformID, DeviceType, 0, NULL, &DeviceCount));

      Devices.resize(DeviceCount);
      oclSafeCall(clGetDeviceIDs(PlatformID, DeviceType, DeviceCount, &Devices[0], &DeviceCount));

      std::cerr << "Found " << DeviceCount << " suitable devices: \n";
      for (cl_uint dev = 0; dev < DeviceCount; dev++) {
	char device_string[1024];	
        oclSafeCall(clGetDeviceInfo(Devices[dev], CL_DEVICE_NAME, sizeof(device_string), &device_string, NULL));
	std::cerr << " " << dev << ": " << device_string << "\n";
      }
    
      return DeviceCount;
    };

    void createQueue(const int dev = 0) {
      // assert(dev < DeviceCount);
      devId = dev;
      
      assert(!CommandQueueFlag);
      assert(!ContextFlag);
      
      cl_int ciErrNum;
      if (dev >= 0) {
	Context = clCreateContext(0, 1, &Devices[devId], NULL, NULL, &ciErrNum);
	oclCheckError(ciErrNum);
	
	CommandQueue = clCreateCommandQueue(Context, Devices[devId], 0, &ciErrNum);
	oclCheckError(ciErrNum);
	
	std::cerr << "Using device: " << devId << std::endl;
      } else {
	int dev = 0;
	while(1) {
	  std::cerr << "Trying device: " << dev << std::endl;
	  Context = clCreateContext(0, 1, &Devices[dev], NULL, NULL, &ciErrNum);
	  if (ciErrNum != CL_SUCCESS) {
	    dev = (dev + 1) % DeviceCount;
	  } else {
	    devId = dev;
	    CommandQueue = clCreateCommandQueue(Context, Devices[devId], 0, &ciErrNum);
	    oclCheckError(ciErrNum);
	    break;
	  }
	  
	}
      }
    }
    
    const cl_context&       get_context()       const {return Context;}
    const cl_command_queue& get_command_queue() const {return CommandQueue;}
    const cl_device_id&     operator[](const int i) const {return Devices[i];}
    const cl_device_id&     device()  const {return Devices[devId];}
    
  };
  
  
  template<class T>
  class memory {
  protected:
    cl_context       Context;
    cl_command_queue CommandQueue;
    bool ContextFlag;

    size_t n;
    cl_mem DeviceMem;
    cl_mem_flags DeviceMemFlags;
    std::vector<T> HostMem;

    void ocl_free() {
      if (n > 0) {
	assert(ContextFlag);
	oclSafeCall(clReleaseMemObject(DeviceMem));
	HostMem.clear();
	n = 0;
      }
    }
    
    void setContext(const cl_context &context, const cl_command_queue &command_queue) {
      assert(!ContextFlag);
      Context      = context;
      CommandQueue = command_queue;
      ContextFlag  = true;
    }

  public:
    memory() :  ContextFlag(false), n(0) {};
    memory(class context &c) :  ContextFlag(false), n(0){setContext(c);}
    memory(class context &c, const int _n, const cl_mem_flags flags = CL_MEM_READ_WRITE) :  ContextFlag(false), n(0) {
      setContext(c);
      allocate(_n, flags);
    };
    memory(class memory &x) :  ContextFlag(false), n(0) {   setContext(x.get_context(), x.get_command_queue());  }
    memory(class memory &x, const int _n, const cl_mem_flags flags = CL_MEM_READ_WRITE) :  ContextFlag(false), n(0) {
      setContext(x.get_context(), x.get_command_queue());
      allocate(_n, flags);
    }
    ~memory() {ocl_free();}
    

    void setContext(const context &c) { setContext(c.get_context(), c.get_command_queue()); }

    const std::vector<T> to_vector() const {return HostMem;}
    
    void allocate(const int _n, const cl_mem_flags flags = CL_MEM_READ_WRITE) {
      if (n > 0) ocl_free();

      assert(ContextFlag);
      DeviceMemFlags = flags;
      n = _n;
      
      HostMem.resize(n);
      memset(&HostMem[0], 0, n * sizeof(T));

      cl_int ciErrNum;
      if (!((flags & CL_MEM_USE_HOST_PTR) == 0)) 
 	DeviceMem = clCreateBuffer(Context, DeviceMemFlags, n*sizeof(T), &HostMem[0], &ciErrNum);
      else                             	      
	DeviceMem = clCreateBuffer(Context, DeviceMemFlags, n*sizeof(T), NULL, &ciErrNum);
      
      oclCheckError(ciErrNum);
    }

    void zeroMem() {
      assert(ContextFlag);
      assert(n > 0);
      memset(&HostMem[0], 0, n*sizeof(T));
      h2d();
    }

    void set(const std::vector<T> &in) {
      const int n = in.size();
      allocate(n, DeviceMemFlags);
      for (int i = 0; i < n; i++) 
	HostMem[i] = in[i];      
    }

    void device2host() {d2h();}
    void host2device() {h2d();}

    void d2h(const cl_bool OCL_BLOCKING = CL_TRUE) {
      assert(ContextFlag);
      assert(n > 0);
      if (!((DeviceMemFlags & CL_MEM_USE_HOST_PTR) == 0)) return;
      oclSafeCall(clEnqueueReadBuffer(CommandQueue, DeviceMem, OCL_BLOCKING, 0,
				      n*sizeof(T),
				      &HostMem[0], 0, 0, 0));
    }

    void h2d(const cl_bool OCL_BLOCKING = CL_TRUE) {
      assert(ContextFlag);
      assert(n > 0);
      if (!((DeviceMemFlags & CL_MEM_USE_HOST_PTR) == 0)) return;
      oclSafeCall(clEnqueueWriteBuffer(CommandQueue, DeviceMem, OCL_BLOCKING, 0,
				       n*sizeof(T),
				       &HostMem[0], 0, 0, 0));
    }

    void copy(const memory &src, const cl_bool OCL_BLOCKING = CL_TRUE) {
      assert(ContextFlag);
      if (n != src.n) {
	ocl_free();
	cmalloc(src.n, DeviceMemFlags);
      }
      oclSafeCall(clEnqueueCopyBuffer(CommandQueue,
                                      src.DeviceMem,
                                      DeviceMem,
                                      0, 0, n*sizeof(T),
                                      0, NULL, NULL));    
      d2h();
    }
    
    const T& operator[](const int i) const {return HostMem[i];};
    T& operator[](const int i) {return HostMem[i];}

    const cl_mem& get_device_mem() const {return DeviceMem;}
    void*   p() const {return (void*)&DeviceMem;}
    void*   ptr() const {return p();}
    const size_t size() const {return n;}
    
    const cl_context&       get_context()       const {return Context;}
    const cl_command_queue& get_command_queue() const {return CommandQueue;}

  };

  class kernel {
  protected:
    char *KernelSource;
    char *KernelFilename;
    char *KernelBinary;
    char *KernelName;

    cl_context       Context;
    cl_command_queue CommandQueue;

    size_t       KernelLength;
    cl_program   Program;
    cl_kernel    Kernel;
    cl_device_id DeviceId;
    std::vector<size_t> GlobalWork;
    std::vector<size_t> LocalWork;

    bool ContextFlag;
    bool KernelFlag;
    bool ProgramFlag;
    bool WorkFlag;   
    
    void clean() {
      KernelSource   = (char*)malloc(1024);
      KernelBinary   = (char*)malloc(1024);
      KernelName     = (char*)malloc(256);
      KernelFilename = (char*)malloc(1024);
      GlobalWork.clear();
      LocalWork.clear();

      ContextFlag = false;
      KernelFlag  = false;
      ProgramFlag = false;
      WorkFlag    = false;
    };
 
    void print_compiler_output() {
      fprintf(stderr, "Compilation of the source file failed: %s \n", KernelFilename);
      char buildLog[10240];
      clGetProgramBuildInfo(Program, DeviceId, CL_PROGRAM_BUILD_LOG, sizeof(buildLog), buildLog, NULL);
      fprintf(stderr,"Compiler output: \n %s \n", buildLog);
      assert(false);
    }

    void setContext(const cl_context &context, const cl_command_queue &command_queue) {
      assert(!ContextFlag);
      Context      = context;
      CommandQueue = command_queue;
      ContextFlag  = true;
    }

  public:

    kernel() {clean();}
    ~kernel() {
      free(KernelSource);
      free(KernelBinary);
      free(KernelName);
      free(KernelFilename);
      if (KernelFlag)  clReleaseKernel (Kernel);
      if (ProgramFlag) clReleaseProgram(Program);
    }

    kernel(class context &c) {
      clean();
      setContext(c);
      DeviceId = c.device();
    }
    kernel(class kernel &k) {
      clean();
      setContext(k.get_context(), k.get_command_queue());
      DeviceId = k.DeviceId;
    }
    void setContext(const context &c) { 
      setContext(c.get_context(), c.get_command_queue()); 
      DeviceId = c.device();
    }
    
    void load_source(const char *kernel_name, 
		     const char *subfolder = "",
                     const char *compilerOptions = "",
                     int maxrregcount = -1) {
      assert(ContextFlag);
      assert(!ProgramFlag);
      
      // create & compile program
      sprintf(KernelFilename, "%s%s", subfolder, kernel_name);
      

      KernelSource = oclLoadProgSource(KernelFilename, "", &KernelLength);
      
      cl_int ciErrNum;
      Program = clCreateProgramWithSource(Context, 1, (const char**)&KernelSource, NULL, &ciErrNum);
      oclCheckError(ciErrNum);
      
#ifdef MAC
      const char* flags = "-cl-mad-enable -DMAC";
#else
      const char* flags = " -cl-mad-enable";
#endif
      
      //Combine flags and custom compile options
      char finalBuildOptions[1024];
      sprintf(finalBuildOptions, "%s %s ", flags, compilerOptions);
      
      if(maxrregcount >= 0) {
        sprintf(finalBuildOptions, "%s %s=%d", finalBuildOptions, "-cl-nv-maxrregcount", maxrregcount);      
      }

      ///the -cl-nv-maxrregcount=n build option when building the kernel (
      //#pragma OPENCL EXTENSION cl_khr_fp64 : enable

      ciErrNum = clBuildProgram(Program, 0, NULL, finalBuildOptions, NULL, NULL);
      
      if(ciErrNum != CL_SUCCESS)      
        print_compiler_output();
      
      ProgramFlag = true;
    }

    void create(const char *kernel_name) {
      assert(ProgramFlag);
      assert(!KernelFlag);
      sprintf(KernelName, kernel_name,"");

      cl_int ciErrNum;
      fprintf(stderr, "Creating kernel %s \n", kernel_name);
      Kernel = clCreateKernel(Program, KernelName, &ciErrNum);
      oclCheckError(ciErrNum);
      KernelFlag = true;
    }

    /////////////

    void setWork(const std::vector<size_t> &global_work, const std::vector<size_t> &local_work) {
      assert(KernelFlag);
      assert(global_work.size() == local_work.size());
      GlobalWork = global_work;
      LocalWork  = local_work;
    }
    
    void setWork(const int nx_threads, const int nx_items,
		 const int ny_threads, const int ny_items) {
      std::vector<size_t> localWork(2), globalWork(2);
      const int ngx = (nx_items - 1) / nx_threads + 1;
      const int ngy = (ny_items - 1) / ny_threads + 1;
      globalWork[0] = ngx*nx_threads;  globalWork[1] = ngy*ny_threads;
      localWork [0] = nx_threads;      localWork [1] = ny_threads;   
      setWork(globalWork, localWork);
    }

    void setWork_block1D(const int n_threads, const int blocks) {
      std::vector<size_t> localWork(2), globalWork(2);
      const int nx = blocks;
      const int ny = 1;
      globalWork[0] = nx*n_threads;  globalWork[1] = ny;
      localWork [0] = n_threads;      localWork[1] = 1;   
      setWork(globalWork, localWork);
    }
    
    void setWork_block2D(const int n_threads, const int blocks) {
      std::vector<size_t> localWork(2), globalWork(2);
      const int nx = (int)std::sqrt(blocks);
      const int ny = (blocks -1)/nx +  1;           
      globalWork[0] = nx*n_threads;  globalWork[1] = ny;
      localWork [0] = n_threads;      localWork[1] = 1;   
      setWork(globalWork, localWork);
    }

    void setWork_1D(const int n_threads, const int items){
      std::vector<size_t> localWork(2), globalWork(2);
      const int ng = (items - 1) / n_threads + 1;
      const int nx = ng;
      const int ny = 1;
      globalWork[0] = nx*n_threads;  globalWork[1] = ny;
      localWork [0] = n_threads;      localWork[1] = 1;   
      setWork(globalWork, localWork);
    }    
    
    void setWork_2D(const int n_threads, const int items) {
      std::vector<size_t> localWork(2), globalWork(2);
      const int ng = (items - 1) / n_threads + 1;
      const int nx = (int)std::sqrt(ng);
      const int ny = (ng - 1)/nx +  1; 
      globalWork[0] = nx*n_threads;  globalWork[1] = ny;
      localWork [0] = n_threads;      localWork[1] = 1;   
      setWork(globalWork, localWork);
    }    
    
    void setWork_threadblock2D(const int nx_threads, const int ny_threads, 
                               const int nx_blocks,  const int ny_blocks) {
      std::vector<size_t> localWork(3), globalWork(3);
         
      globalWork[0] = nx_blocks*nx_threads;  globalWork[1] = ny_blocks*ny_threads;
      localWork [0] = nx_threads;      localWork[1] = ny_threads;   
      localWork [2] = globalWork [2]  = 1;
//       setWork(globalWork, localWork);
 
//       GlobalWork.resize(3);
//       LocalWork.resize(3);
//       
//       GlobalWork[0] = nx_blocks; GlobalWork[1] = ny_blocks; GlobalWork[2] = 1;
//       LocalWork[0] = nx_threads; LocalWork[1] = ny_threads; LocalWork[2] = 1;
      
      setWork(globalWork, localWork);   
      
      WorkFlag = true;      
      
    }
        
    void printWorkSize()
    {
      printf("Blocks: (%ld, %ld, %ld) Threads: (%ld, %ld, %ld) \n", 
              GlobalWork[0], GlobalWork[1], GlobalWork[2],
              LocalWork[0],  LocalWork[1],  LocalWork[2]);             
    }

    template<class T>
    void set_arg(const int arg, void* ptr, const int size = 1)  {
      assert(KernelFlag);
      oclSafeCall(clSetKernelArg(Kernel, arg, size*sizeof(T), ptr));
    }
    
    const int localDim()  const {return  LocalWork[0]* LocalWork[1];};
    const int globalDim() const {return GlobalWork[0]*GlobalWork[1];};
    const int num_groups() const {return globalDim()/localDim();};
    

    void execute(cl_event* event = NULL) {
      assert(KernelFlag);
      const cl_uint work_dim = GlobalWork.size(); 
//       printf("Executing on command queue: %d \n", CommandQueue);
      oclSafeCall(clEnqueueNDRangeKernel(CommandQueue, Kernel, work_dim, 0,
					 &GlobalWork[0],
					 &LocalWork[0],
					 0, NULL, event));
                                         
//       oclSafeCall(clFinish(CommandQueue));                                       
    }

    const cl_kernel&        get_kernel()        const {return Kernel;}
    const cl_program&       get_program()       const {return Program;}
    const cl_context&       get_context()       const {return Context;}
    const cl_command_queue& get_command_queue() const {return CommandQueue;}
    void wait() const {clFinish(CommandQueue);}
    

  };
}

#endif // __OCL_H__
