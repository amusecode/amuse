#ifndef _MY_CUDA_H_
#define _MY_CUDA_H_

#include <cmath>


#include <sys/time.h>

//typedef unsigned int uint;

// struct real4 {
//   real x, y, z, w;
// };
// struct __builtin_align__(16) float4
// {
//   float x, y, z, w;
//   __cuda_assign_operators(float4)
// };


// struct real3 {
//   real x, y, z;
// };
// 
// struct real2 {
//   real x, y;
// };


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <string>
#include <fstream>
#include <cassert>
#include <vector>
#include <cuda.h>
#include <vector_types.h>
#include <builtin_types.h>
#include <iostream>

//Some easy to use typedefs
typedef float4 real4;
typedef float real;

using namespace std;



/*
 #define int clFinish( param ) {  \
    return cuCtxSynchronize();   \
  }*/

#define cl_mem void*


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


  

inline const char* cuPrintError(int err)
{
    switch (err) {
        case CUDA_SUCCESS : return "CUDA_SUCCESS";
        case CUDA_ERROR_INVALID_VALUE : return "CUDA_ERROR_INVALID_VALUE";
        case CUDA_ERROR_OUT_OF_MEMORY : return "CUDA_ERROR_OUT_OF_MEMORY";
        case CUDA_ERROR_NOT_INITIALIZED : return "CUDA_ERROR_NOT_INITIALIZED";
        case CUDA_ERROR_DEINITIALIZED : return "CUDA_ERROR_DEINITIALIZED";
        case CUDA_ERROR_NO_DEVICE : return "CUDA_ERROR_NO_DEVICE";
        case CUDA_ERROR_INVALID_DEVICE : return "CUDA_ERROR_INVALID_DEVICE";
        case CUDA_ERROR_INVALID_IMAGE : return "CUDA_ERROR_INVALID_IMAGE";
        case CUDA_ERROR_INVALID_CONTEXT : return "CUDA_ERROR_INVALID_CONTEXT";
        case CUDA_ERROR_CONTEXT_ALREADY_CURRENT : return "CUDA_ERROR_CONTEXT_ALREADY_CURRENT";
        case CUDA_ERROR_MAP_FAILED : return "CUDA_ERROR_MAP_FAILED";
        case CUDA_ERROR_UNMAP_FAILED : return "CUDA_ERROR_UNMAP_FAILED";
        case CUDA_ERROR_ARRAY_IS_MAPPED : return "CUDA_ERROR_ARRAY_IS_MAPPED";
        case CUDA_ERROR_ALREADY_MAPPED : return "CUDA_ERROR_ALREADY_MAPPED";
        case CUDA_ERROR_NO_BINARY_FOR_GPU : return "CUDA_ERROR_NO_BINARY_FOR_GPU";
        case CUDA_ERROR_ALREADY_ACQUIRED : return "CUDA_ERROR_ALREADY_ACQUIRED";
        case CUDA_ERROR_NOT_MAPPED : return "CUDA_ERROR_NOT_MAPPED";
        case CUDA_ERROR_INVALID_SOURCE : return "CUDA_ERROR_INVALID SOURCE";
        case CUDA_ERROR_FILE_NOT_FOUND : return "CUDA_ERROR_FILE_NOT_FOUND";
        case CUDA_ERROR_INVALID_HANDLE : return "CASE_ERROR_INVALID_HANDLE";
        case CUDA_ERROR_NOT_FOUND : return "CUDA_ERROR_NOT_FOUND";
        case CUDA_ERROR_NOT_READY : return "CUDA_ERROR_NOT_READY";
        case CUDA_ERROR_LAUNCH_FAILED : return "CUDA_ERROR_LAUNCH_FAILED";
        case CUDA_ERROR_LAUNCH_OUT_OF_RESOURCES : return "CUDA_ERROR_LAUNCH_OUT_OF_RESOUCES";
        case CUDA_ERROR_LAUNCH_TIMEOUT : return "CUDA_ERROR_LAUNCH_TIMEOUT";
        case CUDA_ERROR_LAUNCH_INCOMPATIBLE_TEXTURING : return "CUDA_ERROR_LAUNCH_INCOMPATIBLE_TEXTURING";
        case CUDA_ERROR_UNKNOWN : return "CUDA_ERROR_UNKNOWN";
        default : 
          char buff[256];
          sprintf(buff,"%s : %d .", "Unknown CUresult", err);
          fprintf(stderr,"%s\n", buff);
          return "Unknown CUresult.";
          //return buff;
    }
}


#  define CU_SAFE_CALL_NO_SYNC( call ) {                                     \
    CUresult err = call;                                                     \
    if( CUDA_SUCCESS != err) {                                               \
        fprintf(stderr, "Cuda driver error <%s> in file '%s' in line %i.\n",   \
                cuPrintError(err), __FILE__, __LINE__ );                                   \
        exit(EXIT_FAILURE);                                                  \
    } }

#  define CU_SAFE_CALL( call )       CU_SAFE_CALL_NO_SYNC(call);

#  define CU_SAFE_CALL_NO_SYNC_KERNEL( call , kernel ) {                                     \
    CUresult err = call;                                                     \
    if( CUDA_SUCCESS != err) {                                               \
        fprintf(stderr, "Cuda driver error <%s> in file '%s' in line %i. (Kernel: %s)\n",   \
                cuPrintError(err), __FILE__, __LINE__, kernel );                                   \
        exit(EXIT_FAILURE);                                                  \
    } }

#  define CU_SAFE_CALL_KERNEL( call , kernel )       CU_SAFE_CALL_NO_SYNC_KERNEL(call, kernel);


//OpenCL to CUDA macro / functions
__inline__ int clFinish(int param)
{
  return cuCtxSynchronize();   
}



namespace my_dev {

  class context {
  protected:
    size_t dev;
  
    CUcontext hContext;
    CUdevice  hDevice;
    
    int ciDeviceCount;   
    int ciErrNum;

    bool hContext_flag;
    bool hInit_flag;
    bool logfile_flag;
    bool disable_timing;
    
    ofstream *logFile;
    
    int logID;  //Unique ID to every log line
    
    
    //Events:
    CUevent start, stop;
    
    //Compute capabilty, important for default compilation mode
    int ccMajor;
    int ccMinor;
    int defaultComputeMode;    
    
    
  public:

    context() {
      hContext_flag     = false;     
      hInit_flag        = false;
      logfile_flag      = false;
      disable_timing    = false;
      
      CU_SAFE_CALL(cuInit(0));        //Initialize driver API
             
      hInit_flag        = true;                 
    }
    ~context() {
      if (hContext_flag ) cuCtxDetach     (hContext );
    }
    
    //Set the default compute mode based on the current in use device
    void setDefaultComputeMode()
    {
      switch(ccMajor)
      {
        case 1:
          switch(ccMinor)
          {
            case 0:
              defaultComputeMode = CU_TARGET_COMPUTE_10;
              break;
            case 1:
              defaultComputeMode = CU_TARGET_COMPUTE_11;
              break;
            case 2:
              defaultComputeMode = CU_TARGET_COMPUTE_12;
              break;
            case 3:
              defaultComputeMode = CU_TARGET_COMPUTE_13;
              break;              
          }
          break;
       case 2:
          switch(ccMinor)
          {
            case 0:
              defaultComputeMode = CU_TARGET_COMPUTE_20;
              break;
            case 1:
              defaultComputeMode = CU_TARGET_COMPUTE_21;
              break;            
          }
          break;             
      }  
   }    
    
    
    int create(std::ofstream &log, bool disableTiming = false)
    {
      disable_timing = disableTiming;
      logfile_flag = true;
      logFile = &log;
      logID = 0;
      return create(disable_timing); 
    }

    int create(bool disableT = false) {
      assert(hInit_flag);     
      
      disable_timing = disableT;
      
      printf("Creating CUDA context \n");
            
      // Get number of devices supporting CUDA
      ciDeviceCount = 0;
      CU_SAFE_CALL(cuDeviceGetCount(&ciDeviceCount));
            
      printf("Found %d suitable devices: \n",ciDeviceCount);
      for(int i=0; i < ciDeviceCount; i++)
      {
        char device_string[1024];
        CU_SAFE_CALL(cuDeviceGetName(device_string, 1024, i));
        printf(" %d: %s\n",i, device_string);
      }

      return ciDeviceCount;
    }
    
    void createQueue(size_t dev = 0, int ctxCreateFlags = 0) 
    {
      //use CU_CTX_MAP_HOST as flag for zero-copy memory
      //Here we finally create and assign the context to a device
      assert(!hContext_flag);
      assert(hInit_flag);
      this->dev = dev;
      assert((int)dev < ciDeviceCount);
                  
      printf("Trying to use device: %d ...", (int)dev);
      //Get the device handle for dev
      CU_SAFE_CALL(cuDeviceGet(&hDevice, dev)); 

      //Faster and async kernel launches when using large size arrays of local memory
      ctxCreateFlags |= CU_CTX_LMEM_RESIZE_TO_MAX;

      //Create the context for this device handle
      //CU_SAFE_CALL(cuCtxCreate(&hContext, ctxCreateFlags, hDevice));
      
       //
      int res = cuCtxCreate(&hContext, ctxCreateFlags, hDevice);
      if(res != CUDA_SUCCESS)
      {
       printf("failed (error #: %d), now trying all devices starting at 0 \n", res);

	for(int i=0; i < ciDeviceCount; i++)
	{
		printf("Trying device: %d  ...", i);
		CU_SAFE_CALL(cuDeviceGet(&hDevice, i)); 
                if(cuCtxCreate(&hContext, ctxCreateFlags, hDevice) != CUDA_SUCCESS)
		{
			printf("failed!\n");
			if(i+1 == ciDeviceCount)
			{
  			  printf("All devices failed, exit! \n");
			  exit(0);
			}
		}
		else
		{
			printf("success! \n");
			break;
		}
	}
      }
      else
      {
	      printf("success!\n");
      }
      
     //Retrieve CC of the selected device
      cuDeviceComputeCapability(&ccMajor, &ccMinor, hDevice);
      setDefaultComputeMode();      
      
      hContext_flag = true;
    }
    
    
    void startTiming(CUstream stream=0)
    {
      if(disable_timing) return;
      
      CU_SAFE_CALL(cuEventCreate(&start, CU_EVENT_DEFAULT));  
      CU_SAFE_CALL(cuEventCreate(&stop, CU_EVENT_DEFAULT));
      CU_SAFE_CALL(cuEventRecord(start, stream));
    }
    
    //Text and ID to be printed with the log message on screen / in the file
    void stopTiming(const char *text, int type = -1, CUstream stream=0)
    {
      if(disable_timing) return;
      
      CU_SAFE_CALL(cuEventRecord(stop, stream));
      CU_SAFE_CALL(cuEventSynchronize(stop));
      float time;
      CU_SAFE_CALL(cuEventElapsedTime(&time, start, stop));      
      CU_SAFE_CALL(cuEventDestroy(start));
      CU_SAFE_CALL(cuEventDestroy(stop));      
      
      printf("%s took:\t%f\t millisecond\n", text, time);
      
      if(logfile_flag)
      {
        (*logFile) << logID++ << "\t"  << type << "\t" << text << "\t" << time << endl;
      }
    }
    
    /////////////
    
    CUcontext&       get_context()       {return hContext;}
    int              get_command_queue() {return 0;}
    //cl_device_id      operator[](int i)   {return hDevices[i];}
    CUdevice&      get_device()        {return hDevice;}

    const int&  getDefaultComputeMode() const {return defaultComputeMode;}
    //////////       
    
  };
  
  
  ////////////////////////////////////////
  
    //Class to handle streams / queues
  
  class dev_stream
  {
    private:
      CUstream stream;
      
    public:     
      dev_stream(unsigned int flags = 0)
      {
        createStream(flags);
      }
            
      
      void createStream(unsigned int flags = 0)
      {
         CU_SAFE_CALL(cuStreamCreate(&stream, flags));
      }   
      
      void destroyStream()
      {
        CU_SAFE_CALL(cuStreamDestroy(stream));
      }
      
      void sync()
      {
        CU_SAFE_CALL(cuStreamSynchronize(stream));        
      }
      
      CUstream s()
      {
        return stream;
      }
      
      
      ~dev_stream() {
      destroyStream();
    }
  };
  
  
  /*
  class dev_texture
  {
    private:
      CUtexref          texRef; //Used if memory is bind to a texture
      bool              texSet;
      CUdeviceptr       memPtr;
      size_t            size; //In bytes
      
    public:      
      CUtexref getTexture()             { return texRef; }
      void     setTexture(CUtexref tex) { texRef = tex; }
      
      CUdeviceptr getMemPtr()             { return memPtr; }
      void     setMemPtr(CUdeviceptr ptr, size_t sizeT)
      { 
        memPtr = ptr; 
        size   = sizeT;
        
        //Allocate / Assign the memory to the texture
        CU_SAFE_CALL(cuTexRefSetAddress(0, texRef, memPtr, size));   
      }      
    
      bool getIsTexSet()                { return texSet; }
      void setTexIsSet(bool val)        { texSet = val; } //If false texture has to be reset
      
      dev_texture()
      {
        texSet = false;        
      }
      
      ~dev_texture()
      {
        //Do we have to free resources?
      }      
  };
  */
  
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
      
      size_t free, total;
      cuMemGetInfo(&free, &total); 
      printf("Build-in usage: free: %ld bytes ( %ld MB , total: %ld) \n", free, free / (1024*1024), total / (1024*1024));     
      
    }  
    
    static long long getMaxMemUsage()
    {      
      return maxMemUsage;
    }   
    
  };


  template<class T>
  class dev_mem : base_mem {
  protected:
    CUcontext hContext;
    
    typedef struct textureInfo
    {
      CUtexref texture;
      int      texOffset; //The possible extra offset when using textures and combined memory
      int      texSize;
    } textureInfo;        
        
 
    vector<textureInfo> textures;
    
    int size;
    CUdeviceptr hDeviceMem;
    T           *host_ptr;
    void        *DeviceMemPtr;
    void        *tempDeviceMemPtr;
        
    bool pinned_mem, context_flag, flags;
    bool hDeviceMem_flag;
    bool childMemory; //Indicates that this is a shared buffer that will be freed by a parent
    
    void cuda_free() {      
      if(childMemory) //Only free if we are NOT a child
      {
        return;
      }

      assert(context_flag);
      if (hDeviceMem_flag)
      {
	assert(size > 0);
	CU_SAFE_CALL(cuMemFree(hDeviceMem));
        decreaseMemUsage(size*sizeof(T));
        
        if(pinned_mem){
          CU_SAFE_CALL(cuMemFreeHost((void*)host_ptr));}
        else{
          free(host_ptr);}
          hDeviceMem_flag = false;
      }
    } //cuda_free

  public:
    

    ///////// Constructors

    dev_mem() {
      size              = 0;
      pinned_mem        = false;
      hDeviceMem_flag   = false;
      context_flag      = false;
      host_ptr          = NULL;
      childMemory       = false;
    }

    dev_mem(class context &c) {
      size              = 0;      
      pinned_mem        = false;
      context_flag      = false;
      hDeviceMem_flag   = false;
      host_ptr          = NULL;
      childMemory       = false;
      setContext(c);
    }
    
    //CUDA has no memory flags like opencl
    //so just put it to 0 and keep function format same for 
    //compatability
    dev_mem(class context &c, int n, bool zero = false,
	    int flags = 0, bool pinned = false) {
      context_flag      = false;
      childMemory       = false;      
      hDeviceMem_flag   = false;
      pinned_mem        = pinned;
      size              = 0;
      setContext(c);
      if (zero) this->ccalloc(n, flags);
      else      this->cmalloc(n, flags);
    }
    
    dev_mem(class context &c, std::vector<T> data,
	    int flags = 0,  bool pinned = false) {
      context_flag    = false;
      hDeviceMem_flag = false;
      childMemory       = false;      
      pinned_mem      = pinned;
      size            = 0;
      setContext(c);
      this->cmalloc(data, flags);
    }
    
    void free_mem()
    {
      cuda_free();
    }

    //////// Destructor
    
    ~dev_mem() {      
        cuda_free();
    }
    
    ///////////

    void setContext(class context &c) {
      if(context_flag)
      {
	      //Check if the context is changed, that is not allowed!
	      assert(c.get_context() == this->hContext);
      }

      this->hContext   = c.get_context();     
      context_flag     = true;
    }


    ///////////
    
    //Get the reference of memory allocated by another piece of memory
    
    void cmalloc_copy(bool pinned, bool flags, CUdeviceptr cudaMem, 
                      void* ParentHost_ptr, int offset, int n,
                      int allignOffset)
    {
      assert(context_flag);
//       assert(!hDeviceMem_flag);
//       this->pinned_mem  = src_buffer.get_pinned();      
//       this->flags       = src_buffer.get_flags();     
      this->pinned_mem  = pinned;
      this->flags       = flags;
      this->childMemory = true;

      size = n;
        
      //Dont forget to add the allignment values
      host_ptr = (T*)ParentHost_ptr +  allignOffset*sizeof(uint);

      hDeviceMem   = cudaMem + offset*sizeof(uint) + allignOffset*sizeof(uint);
      DeviceMemPtr = (void*)(size_t)(hDeviceMem);

      hDeviceMem_flag = true;      
    /*  
    T& operator[] (int i){ return host_ptr[i]; }
    CUdeviceptr& d() {return hDeviceMem;}
   // void*   p() {return (void*)hDeviceMem;}
//     CUdeviceptr   p() {return hDeviceMem;}
    void*   p() {return &DeviceMemPtr;}
    void*   p(int offset)
    {      
      //Calculate the new memory offset 
      tempDeviceMemPtr = (void*)(size_t)(hDeviceMem + offset*sizeof(T));          
      return &tempDeviceMemPtr;   
    }      */    
      
    }

    void cmalloc(int n, bool pinned = false,  int flags = 0) {
      assert(context_flag);
//       assert(!hDeviceMem_flag);
      this->pinned_mem = pinned;      
      this->flags = flags;
      if (size > 0) cuda_free();
      size = n;
            
      if(pinned_mem){    
        CU_SAFE_CALL(cuMemAllocHost((void**)&host_ptr, size*sizeof(T)));}
      else{
        host_ptr = (T*)malloc(size*sizeof(T));}
        
      CU_SAFE_CALL(cuMemAlloc(&hDeviceMem, size*sizeof(T)));
      increaseMemUsage(size*sizeof(T));
      DeviceMemPtr = (void*)(size_t)hDeviceMem;    
     
//       if (flags & CL_MEM_USE_HOST_PTR == 0) hDeviceMem = clCreateBuffer(hContext, flags, size*sizeof(T), NULL, &err);
//       else                             	    hDeviceMem = clCreateBuffer(hContext, flags, size*sizeof(T), &host_ptr[0], &err);
//       oclCheckError(err);
      hDeviceMem_flag = true;
    }

    void ccalloc(int n, bool pinned = false, int flags = 0) {
      assert(context_flag);
     // assert(!hDeviceMem_flag);
      
      this->pinned_mem = pinned;      
      this->flags = flags;
      if (size > 0) cuda_free();
      size = n;
      
      if(pinned_mem)      
        cuMemAllocHost((void**)&host_ptr, size*sizeof(T));
      else
        host_ptr = (T*)calloc(size, sizeof(T));
      
      CU_SAFE_CALL(cuMemAlloc(&hDeviceMem, size*sizeof(T)));           
      
      CU_SAFE_CALL(cuMemsetD8(hDeviceMem, 0, size*sizeof(T)));     
      increaseMemUsage(size*sizeof(T));
      DeviceMemPtr = (void*)(size_t)hDeviceMem;    
      
//       if (flags & CL_MEM_USE_HOST_PTR == 0) hDeviceMem = clCreateBuffer(hContext, flags, size*sizeof(T), NULL, &err);
//       else                             	    hDeviceMem = clCreateBuffer(hContext, flags, size*sizeof(T), &host_ptr[0], &err);
//       oclCheckError(err);
//      h2d();
      hDeviceMem_flag = true;
    }
    
//     void cresize(int n, bool pinned = false)
    //Set reduce to false to not reduce the size, to speed up pinned memory buffers
    void cresize(int n, bool reduce = true)     
    {
      if(size == n)     //No need if we are already at the correct size
        return;
      
      if(size > n && reduce == false) //Do not make the memory size smaller
      {
        return;
      }
      
//       d2h();    //Get datafrom the device
   
      if(pinned_mem)
      {
        //No realloc function so do it by hand
        T *tmp_ptr;            
        CU_SAFE_CALL(cuMemAllocHost((void**)&tmp_ptr, n*sizeof(T)));        
        //Copy old content to newly allocated mem
        int tmpSize = min(size,n);
               
        //Copy the old data to the new pointer and free the old location
        memcpy (((void*) tmp_ptr), ((void*) host_ptr), tmpSize*sizeof(T)); 
        CU_SAFE_CALL(cuMemFreeHost((void*)host_ptr));
        host_ptr = tmp_ptr;
      }
      else
      {
        //Resizes the current array
        //New size is smaller, don't do anything with the allocated memory                
        host_ptr = (T*)realloc(host_ptr, n*sizeof(T));
      }

      //Free the memory and re-allocate it
      CU_SAFE_CALL(cuMemFree(hDeviceMem));
      decreaseMemUsage(size*sizeof(T));    
      CU_SAFE_CALL(cuMemAlloc(&hDeviceMem, n*sizeof(T)));        
      increaseMemUsage(n*sizeof(T));    
      DeviceMemPtr = (void*)(size_t)hDeviceMem;    
      size = n;
      
      //Rebind the textures
      for(unsigned int i = 0; i < textures.size(); i++)
      { 
        //Sometimes textures are only bound to a part of the total memory
        //So check this
        if(textures[i].texOffset < 0)
        {
          CU_SAFE_CALL(cuTexRefSetAddress(0, textures[i].texture, hDeviceMem, n*sizeof(T)));   
        }
        else
        {
          CUdeviceptr tempDeviceMemPtr = (CUdeviceptr) a(textures[i].texOffset);
          CU_SAFE_CALL(cuTexRefSetAddress(0, textures[i].texture, tempDeviceMemPtr, sizeof(T)*textures[i].texSize));             
        }        
      }      

//       h2d();    //Move data to the device again
    }

            
    
    //Set the memory to zero
    void zeroMem()
    {
      assert(context_flag);
      assert(hDeviceMem_flag);
            
      memset(host_ptr, 0, size*sizeof(T));      
      CU_SAFE_CALL(cuMemsetD8(hDeviceMem, 0, size*sizeof(T)));           
    }

    /*
      //Not implemented in CUDA
      void malloc(std::vector<T> &data, int flags = 0) {
      assert(context_flag);
      assert(!hDeviceMem_flag);
      this->flags = flags;
      if (size > 0) cuda_free();
      size = data.size();
      host_ptr = data;
      int err;
      
      printf("TODO: This function does NOT copy memory in CUDA!! \n");      
      CU_SAFE_CALL(cuMemAlloc(&hDeviceMem, size*sizeof(T)));
//       if (flags & CL_MEM_USE_HOST_PTR == 0) {
// 	hDeviceMem = clCreateBuffer(hContext, flags | CL_MEM_COPY_HOST_PTR,
// 				    size*sizeof(T), 
// 				    &host_ptr[0], 
// 				    &err);
//       } else {
// 	hDeviceMem = clCreateBuffer(hContext, flags, size*sizeof(T), &host_ptr[0], &err);
//       }
//      oclCheckError(err);
      hDeviceMem_flag = true;
    }*/
    
    ///////////

    //////////////

    void d2h(bool OCL_BLOCKING = true, CUstream stream = 0)   {      
      assert(context_flag);
      assert(hDeviceMem_flag);
      assert(size > 0);
      
      if(OCL_BLOCKING)
      {
        CU_SAFE_CALL(cuMemcpyDtoH(&host_ptr[0], hDeviceMem, size*sizeof(T)));
      }
      else
      {
        //Async copy, ONLY works for page-locked memory therefore default parameter
        //is blocking.
        assert(pinned_mem);
        CU_SAFE_CALL(cuMemcpyDtoHAsync(&host_ptr[0], hDeviceMem, size*sizeof(T), stream));          
      }
    }
    
    //D2h that only copies a certain number of items to the host
    void d2h(int number, bool OCL_BLOCKING = true, CUstream stream = 0)   {      
      assert(context_flag);
      assert(hDeviceMem_flag);
      assert(size > 0);
      
      if(OCL_BLOCKING)
      {
        CU_SAFE_CALL(cuMemcpyDtoH(&host_ptr[0], hDeviceMem, number*sizeof(T)));
      }
      else
      {
        //Async copy, ONLY works for page-locked memory therefore default parameter
        //is blocking.
        assert(pinned_mem);
        CU_SAFE_CALL(cuMemcpyDtoHAsync(&host_ptr[0], hDeviceMem, number*sizeof(T), stream));          
      }
    }    
    
    void h2d(bool OCL_BLOCKING  = true, CUstream stream = 0)   {
      assert(context_flag);
      assert(hDeviceMem_flag);
      assert(size > 0);
      //if (flags & CL_MEM_USE_HOST_PTR == 0) return;
      if(OCL_BLOCKING)
      {
        CU_SAFE_CALL(cuMemcpyHtoD(hDeviceMem, &host_ptr[0], size*sizeof(T)));
      }
      else
      {
        //Async copy, ONLY works for page-locked memory therefore default parameter
        //is blocking.
        assert(pinned_mem);
        CU_SAFE_CALL(cuMemcpyHtoDAsync(hDeviceMem, host_ptr, size*sizeof(T), stream));          
      }        
    }
    
    //JB: Modified this so that it copies a device buffer to an other device
    //buffer, and the host buffer to the other host buffer
    void copy(dev_mem &src_buffer, int n, bool OCL_BLOCKING = true)   {
      assert(context_flag);
      assert(hDeviceMem_flag);
      if (size < n) {
	cuda_free();
	cmalloc(n, flags);
	size = n;
        printf("Resize in copy \n");
      }
      
      //Copy on the device
      CU_SAFE_CALL(cuMemcpyDtoD(hDeviceMem, src_buffer.d(), n*sizeof(T)));
      //Copy on the host
      memcpy (((void*) &host_ptr[0]), ((void*) &src_buffer[0]), n*sizeof(T));                                          
    }
    
    /////////
    
    T& operator[] (int i){ return host_ptr[i]; }
    
    CUdeviceptr   get_devMem() {return hDeviceMem;}
    CUdeviceptr& d() {return hDeviceMem;}
   // void*   p() {return (void*)hDeviceMem;}
//     CUdeviceptr   p() {return hDeviceMem;}
    void*   p() {return &DeviceMemPtr;}
    void*   p(int offset)
    {      
      //Calculate the new memory offset 
      tempDeviceMemPtr = (void*)(size_t)(hDeviceMem + offset*sizeof(T));          
      return &tempDeviceMemPtr;   
    }    
    void*   a(int offset)
    {      
      //Calculate the new memory offset       
      return (void*)(size_t)(hDeviceMem + offset*sizeof(T));
    }     
    
    //Texture related functions
    int      addTexture(CUtexref tex, int offset, int texSize)
    {
      textureInfo temp;
      temp.texture      = tex;
      temp.texOffset    = offset;
      temp.texSize      = texSize;
      
      textures.push_back(temp);
      return textures.size()-1;
    }


    int  get_size(){return size;}
    bool get_pinned(){return pinned_mem;}
    bool get_flags(){return flags;}
  };     // end of class dev_mem

  ////////////////////
  
   
  class kernel {
  protected:    
    char *hKernelFilename;    
    char *hKernelName;

    //cl_program hProgram;       
    
    CUcontext   hContext;
    CUdevice    hDevice;
    CUmodule    cuModule;
    CUfunction  hKernel;

    vector<size_t> hGlobalWork;
    vector<size_t> hLocalWork;
    
    vector<void*> argumentList;
    vector<int> argumentOffset;
    
    //Kernel argument stuff
//     enum {MAXKERNELARGUMENTS = 128};
    #define MAXKERNELARGUMENTS 128
    typedef struct kernelArg
    {
      int alignment;    //Alignment of the variable type
      int sizeoftyp;    //Size of the variable type
      void* ptr;        //The pointer to the memory
      int size;         //The number of elements (incase of shared memory)
//       my_dev::dev_texture *texture; //The possible reference to a texture
      CUtexref texture;
      int      texOffset; //The possible extra offset when using textures and combined memory
      int      texSize;
      int      texIdx;
    } kernelArg;        
    
    std::vector<kernelArg> kernelArguments;        


    bool context_flag;
    bool kernel_flag;
    bool program_flag;
    bool work_flag;
    
    int sharedMemorySize;
    int paramOffset;
    
    int computeMode;

  public:

    kernel() {
      hKernelName     = (char*)malloc(256);
      hKernelFilename = (char*)malloc(1024);
      hGlobalWork.clear();
      hLocalWork.clear();

      context_flag = false;
      kernel_flag  = false;
      program_flag = false;
      work_flag    = false;
      
      sharedMemorySize = 0;
      paramOffset      = 0;

      //Kernel argument stuff
      kernelArguments.resize(MAXKERNELARGUMENTS);
      kernelArg argTemp; 
      argTemp.alignment = -1;   argTemp.sizeoftyp = -1; 
      argTemp.ptr       = NULL; argTemp.size      = -1;
      argTemp.texIdx    = -1;   argTemp.texOffset = 0;
      kernelArguments.assign(MAXKERNELARGUMENTS, argTemp);         
      
    }
    ~kernel() {      
      free(hKernelName);
      free(hKernelFilename);
     //if (kernel_flag) clReleaseKernel(hKernel);
      if (program_flag) cuModuleUnload(cuModule);
    }

    kernel(class context &c) {
      hKernelName   = (char*)malloc(256);
      hKernelFilename = (char*)malloc(1024);
      hGlobalWork.clear();
      hLocalWork.clear();

      context_flag = false;
      kernel_flag  = false;
      program_flag = false;
      work_flag    = false;
      
      sharedMemorySize = 0;
      paramOffset      = 0;
      setContext(c);
    }
    
    CUmodule getModule()
    {
      return cuModule;
    }

    ////////////

    void setContext(class context &c) {
      assert(!context_flag);
      this->hContext   = c.get_context();     
      this->hDevice    = c.get_device();
      
      computeMode      = c.getDefaultComputeMode();
      context_flag     = true;      
    }

    ////////////
    
    void load_source(const char *fileName, string &ptx_source)
    {
      FILE *fp;
      
      fp = fopen(fileName, "rb");
      
      if(fp == NULL)
      {
        printf("Cannot open source file: %s \n", fileName);
        exit(-1);
      }
      
      fseek(fp, 0, SEEK_END);
      int file_size = ftell(fp);      
      ptx_source.reserve(file_size + 512);     
      fseek(fp, 0, SEEK_SET);
      size_t read = fread(&ptx_source[0], sizeof(char), file_size, fp);
      
      //Set the last char to NULL else old values in the extra memory
      //buffer will crash the compiler....
      ptx_source[file_size] = '\0';

      if(read == 0)
      {
        printf("Cannot read source file: %s \n", fileName);
        exit(-1);        
      }
      
      fclose(fp);          
    }
    
    void load_source(const char *kernel_name, const char *subfolder,
                     const char *compilerOptions = "",
                     int maxrregcount = -1,
                     int architecture = CU_TARGET_COMPUTE_20) {
//                     int architecture = CU_TARGET_COMPUTE_12) {
      assert(context_flag);
      assert(!program_flag);
      
      //In cuda version we assume that the code is already compiled into ptx
      //so that the file loaded/specified is in fact a PTX file
      sprintf(hKernelFilename, "%s%s", subfolder, kernel_name);
      
      printf("Loading source: %s ...", hKernelFilename);
      
      string temp = string(hKernelFilename);
      if(temp.rfind("ptx") != string::npos)
      {  
        //Load PTX source
        const unsigned int maxJitOptions = 6;
        CUjit_option *jitOptions = new CUjit_option[maxJitOptions];
        void **jitOptVals        = new void*[maxJitOptions];
        
        
        int jitOptionCount = 0;
        //use JIT compiling to set the max register number
        if(maxrregcount > 0)
        {
          //Set the maximum number of registers option
          jitOptions[jitOptionCount] = CU_JIT_MAX_REGISTERS;
          int jitRegCount = maxrregcount;
          jitOptVals[jitOptionCount] = (void *)jitRegCount;        
          jitOptionCount++;        
        }
        
        //Fermi requires at least compute mode 2, so if device is compute mode 2
        //but the given option is not compute mode 2 (eg constant <= 3) we change 
        //it to compute mode 2.0
        //TODO should make this a bit more strict and usefull
        int maxArchitecture = architecture;
        if(architecture <= 3 && computeMode >= 4)
          maxArchitecture = 4;
        if(computeMode < architecture)
          maxArchitecture = computeMode;
        
            
       //Check for double precision, only use it if device supports it
        //Set the architecture
        {                
          jitOptions[jitOptionCount] = CU_JIT_TARGET;
          int arch = maxArchitecture;
//           int arch = architecture;
          jitOptVals[jitOptionCount] = (void *)arch;        
          jitOptionCount++;  
        }     
        
        
        // set up size of compilation log buffer                                                                                                    
        jitOptions[jitOptionCount] = CU_JIT_INFO_LOG_BUFFER_SIZE_BYTES;                                                                                          
        int jitLogBufferSize = 1024;                                                                                                                
        jitOptVals[jitOptionCount] = (void *)(size_t)jitLogBufferSize;  
        jitOptionCount++;
        
        
        // set up pointer to the compilation log buffer                                                                                             
        jitOptions[jitOptionCount] = CU_JIT_INFO_LOG_BUFFER;                                                                                                     
        char *jitLogBuffer = new char[jitLogBufferSize];                                                                                            
        jitOptVals[jitOptionCount] = jitLogBuffer;                                                                                                               
        jitOptionCount++;
                                                          
        
        string ptxSource;
        load_source(hKernelFilename, ptxSource);
        
  //       hier gebleven bij jit moeten we source code inladen
  //       sterker nog dat willen we altijd, dus altijd code inladen! Woo!
        
  //         CU_SAFE_CALL(cuModuleLoad(&cuModule, hKernelFilename));
    
  //         jitOptionCount = 0;
          CU_SAFE_CALL(cuModuleLoadDataEx(&cuModule, ptxSource.c_str(), jitOptionCount, jitOptions, (void **)jitOptVals));
          
          // printf("> PTX JIT log:\n%s\n", jitLogBuffer);        
          
          delete[] jitOptVals;
          delete[] jitOptions;
          delete[] jitLogBuffer;
      }
      else
      {
        //Load CUBIN source
        string ptxSource;
        load_source(hKernelFilename, ptxSource);
        CU_SAFE_CALL(cuModuleLoad(&cuModule, hKernelFilename));
      }
      
      printf("done!\n");


      program_flag = true;
    }

    void create(const char *kernel_name) {
      assert(program_flag);
      assert(!kernel_flag);
      sprintf(hKernelName, kernel_name,"");
      
      printf("%s \n", kernel_name);
      CU_SAFE_CALL(cuModuleGetFunction(&hKernel, cuModule, hKernelName));

      kernel_flag = true;
    }
    
    //NVIDIA macro
    #define ALIGN_UP(offset, alignment) (offset) = ((offset) + (alignment) - 1) & ~((alignment) -1)
    void completeArguments()
    {
      //Reset the parameter offset and amount of shared memory
      paramOffset       = 0;
      sharedMemorySize  = 0;
      
      //Loop over all set arguments and set them
      for(int i=0; i < MAXKERNELARGUMENTS; i++)
      {      
        //First of all check if this argument has to be set or that we've finished already
        if(kernelArguments[i].size == -1)
          continue;
        
        if(kernelArguments[i].size == -2)
        {
          //This is a texture
          CUtexref cu_texref = kernelArguments[i].texture;
          
          //Bind the texture and contiue to the next kernel parameter
          CU_SAFE_CALL(cuParamSetTexRef(hKernel, CU_PARAM_TR_DEFAULT, cu_texref));              
          
          continue;
        }
      
        //Now, check if this is a shared memory argument
        if(kernelArguments[i].ptr == NULL && kernelArguments[i].size > 1)
        { 
//           printf("%d is shared mem!\n", i);
          //Increase the shared memory size
          sharedMemorySize  +=  kernelArguments[i].size*kernelArguments[i].sizeoftyp;             
        }
        else
        {
//           printf("%d is real mem!\n", i);
          //This is an actual argument that we have to set
          ALIGN_UP(paramOffset, kernelArguments[i].alignment);  
          CU_SAFE_CALL(cuParamSetv(hKernel, paramOffset,
                                   kernelArguments[i].ptr, kernelArguments[i].sizeoftyp));     
          paramOffset += kernelArguments[i].sizeoftyp;                                        
        } //end if
      }//end for
    }//end completeArguments    

    //'size'  is used for dynamic shared memory
    //Cuda does not have a function like clSetKernelArg
    //therefore we keep track of a vector with arguments
    //that will be processed when we launch the kernel
    template<class T>
    void set_arg(unsigned int arg, void* ptr, int size = 1)  {
      assert(kernel_flag);
      
      //TODO have to check / think about if we want size default initialised
      //to 1 or to zero
      
      kernelArg tempArg;
      tempArg.alignment = __alignof(T);
      tempArg.sizeoftyp = sizeof(T);
      tempArg.ptr       = ptr;
      tempArg.size      = size;
      tempArg.texIdx    = -1;
      tempArg.texture   = 0;
      tempArg.texSize   = -1;
      tempArg.texOffset = -1;
      kernelArguments[arg] = tempArg;
 
      return;
    }

    //Offset and mem_size should both be set if you want only part of an array
    //bound to a texture
    template<class T>
    void set_arg(unsigned int arg, my_dev::dev_mem<T> &memobj, int adSize,
                 const char *textureName, int offset = -1, int mem_size = -1)  { //Texture 
      assert(kernel_flag);

      //Only configure the texture if this is the first call
      if(kernelArguments[arg].size != -2)
      { 
        CUtexref cu_texref;
        //Load the texture reference from the module, only required once
        CU_SAFE_CALL(cuModuleGetTexRef(&cu_texref, cuModule, textureName));
        
        //Texture format
        CU_SAFE_CALL(cuTexRefSetFormat(cu_texref, CU_AD_FORMAT_FLOAT, adSize));  
        
        //Assign memory
        if(offset < 0)
        {
          CU_SAFE_CALL(cuTexRefSetAddress(0, cu_texref, memobj.d(), sizeof(T)*memobj.get_size()));   
        }
        else
        {
          //Get the offsetted memory location
          //Fail if mem_size is not set!
          assert(mem_size > 0);
          CUdeviceptr tempDeviceMemPtr = (CUdeviceptr) memobj.a(offset);
          CU_SAFE_CALL(cuTexRefSetAddress(0, cu_texref, tempDeviceMemPtr, sizeof(T)*mem_size));
//           CU_SAFE_CALL(cuTexRefSetAddress(0, cu_texref, memobj.d(), sizeof(T)*mem_size));   
        }
        
        //Now store it in the calling memory object incase the memory is 
        //realloacted and the texture reference has to be updated        
        int idx = memobj.addTexture(cu_texref, offset, mem_size);
        
        kernelArg tempArg;
        tempArg.size          = -2;   //Texture
        tempArg.texIdx        = idx;
        tempArg.texture       = cu_texref;
        tempArg.texOffset     = offset;
        tempArg.texSize       = mem_size;
        tempArg.ptr           = NULL;
        kernelArguments[arg]  = tempArg;        
      }
      else
      {
        //Assign memory, has to be done EVERY kernel call otherwise things mess up!!
        CUtexref cu_texref = kernelArguments[arg].texture;     
        
        //Change the offsets if needed
        kernelArguments[arg].texOffset = offset;
        kernelArguments[arg].texSize   = mem_size;
        
        //Check for the texture offset and texture size, if they are both set
        if(kernelArguments[arg].texOffset >= 0) assert(kernelArguments[arg].texSize > 0);
        
        if(kernelArguments[arg].texOffset < 0)
        {
          CU_SAFE_CALL(cuTexRefSetAddress(0, cu_texref, memobj.d(), sizeof(T)*memobj.get_size()));   
        }
        else
        {          
          CUdeviceptr tempDeviceMemPtr = (CUdeviceptr) memobj.a(kernelArguments[arg].texOffset);
          CU_SAFE_CALL(cuTexRefSetAddress(0, cu_texref, tempDeviceMemPtr, sizeof(T)*kernelArguments[arg].texSize));             
        }
      }

      return;
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
        if(blocks >= 65536)
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
      
      hGlobalWork.resize(3);
      hLocalWork.resize(3);
      
      hLocalWork[0] = local_work[0];
      hGlobalWork[0] = global_work[0];
      
      hLocalWork[1]  = (local_work.size() > 1) ? local_work[1] : 1;
      hGlobalWork[1] = (global_work.size() > 1) ? global_work[1] : 1;
      
      hLocalWork[2]  = (local_work.size() > 2) ? local_work[2] : 1;
      
      //Since the values between CUDA and OpenCL differ:
      //Cuda is specific size of each block, while OpenCL
      //is the combined size of the lower blocks and this block
      //we have to divide the values
      
      hGlobalWork[0] /= hLocalWork[0];
      hGlobalWork[1] /= hLocalWork[1];
      hGlobalWork[2] /= hLocalWork[2];
      
      work_flag = true;
    }


    void execute(vector<size_t> global_work, vector<size_t> local_work, 
		 int* event = NULL) {
      setWork(global_work, local_work);
      
      completeArguments();
    

      CU_SAFE_CALL(cuParamSetSize(hKernel, paramOffset));
       
      if(sharedMemorySize > 0)
        CU_SAFE_CALL(cuFuncSetSharedSize(hKernel, sharedMemorySize));
      
      CU_SAFE_CALL(cuFuncSetBlockShape(hKernel, hLocalWork[0], hLocalWork[1], hLocalWork[2]));
      CU_SAFE_CALL(cuLaunchGridAsync(hKernel, hGlobalWork[0], hGlobalWork[1], 0));
     
      //TODO put this in a define sequence so we can compile it with and without debug
//       printf("Waiting on kernel: %s to finish...", hKernelName);
      CU_SAFE_CALL(cuCtxSynchronize());      
//       printf("finished \n");
      
    }
    
    void printWorkSize()
    {
      printf("Blocks: (%ld, %ld, %ld) Threads: (%ld, %ld, %ld) \n", 
              hGlobalWork[0], hGlobalWork[1], hGlobalWork[2],
              hLocalWork[0], hLocalWork[1], hLocalWork[2]);             
    }
    
    void execute(CUstream hStream = 0, int* event = NULL) {
      assert(kernel_flag);
      assert(work_flag);

      completeArguments();     
      CU_SAFE_CALL(cuParamSetSize(hKernel, paramOffset));
      
      if(sharedMemorySize > 0)
        CU_SAFE_CALL(cuFuncSetSharedSize(hKernel, sharedMemorySize));
      
      CU_SAFE_CALL(cuFuncSetBlockShape(hKernel, hLocalWork[0], hLocalWork[1], hLocalWork[2]));
      CU_SAFE_CALL(cuLaunchGridAsync(hKernel, hGlobalWork[0], hGlobalWork[1], hStream));   


//       cuFuncSetCacheConfig(hKernel, CU_FUNC_CACHE_PREFER_L1);
      //TODO put this in a define sequence so we can compile it with and without debug
//       fprintf(stderr,"Waiting on kernel: %s to finish...", hKernelName);
      CU_SAFE_CALL_KERNEL(cuCtxSynchronize(), hKernelName);   
//       fprintf(stderr,"finished \n");

    }
    ////

    CUfunction&  get_kernel() {return hKernel;}
    CUmodule& get_program() {return cuModule;}  //Program == module?
      
  };
  
  
   
}     // end of namespace my_cuda



#endif // _MY_CUDA_H_


