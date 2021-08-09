#include "evolve.h"

#ifdef EVOLVE_OPENCL

#include <stdio.h>
#include <stdlib.h>

#define CL_TARGET_OPENCL_VERSION 200
#ifdef __APPLE__
    #include "OpenCL/opencl.h"
#else
    #include "CL/cl.h"
#endif

#include "evolve_cl.h"
#include "evolve_kern.clh"

static cl_device_id device_id;             // compute device id 
static cl_context context;                 // compute context
static cl_command_queue queue;             // compute command queue
static cl_program program;                 // compute program
static cl_kernel kick_krn;                 // kernels
static cl_kernel timestep_krn;
static cl_kernel potential_krn;

static cl_mem _ipos, _ivel, _jpos, _jvel, _acc, _pot, _timestep; // opencl buffer objects

static UINT nalloc=0;

#define SELECT(x,a,b,c)   (x==0?a:(x==1?b:c))

void init_cl()
{
  cl_int err;
  char *name;
  
  err = clGetDeviceIDs(NULL, SELECT(opencl_device_type,CL_DEVICE_TYPE_DEFAULT,CL_DEVICE_TYPE_CPU,CL_DEVICE_TYPE_GPU), 1, &device_id, NULL);
  if (err != CL_SUCCESS) ENDRUN("OpenCL could not get device of requested type %s", SELECT(opencl_device_type,"default","cpu","gpu"));
  {
    size_t size;
    clGetDeviceInfo(device_id, CL_DEVICE_NAME, 0, NULL, &size);
    name=(char*) malloc(size);
    clGetDeviceInfo(device_id, CL_DEVICE_NAME, size, name, NULL);
  }

  context = clCreateContext(0, 1, &device_id, NULL, NULL, &err);
  if (!context || err!=CL_SUCCESS) ENDRUN("OpenCL failed to create context");

  queue = clCreateCommandQueueWithProperties(context, device_id, NULL, &err);
  queue = clCreateCommandQueueWithProperties(context, device_id, NULL, &err);
  if (!queue || err!=CL_SUCCESS) ENDRUN("OpenCL failed to create command queue");

  program = clCreateProgramWithSource(context, 1, (const char **) &srcstr, NULL, &err);
  if (!program || err!=CL_SUCCESS) ENDRUN("OpenCL failed to create program")

  err = clBuildProgram(program, 0, NULL, NULL, NULL, NULL);
  if (err != CL_SUCCESS)
  {
        size_t len;
        char *buffer;
        printf("OpenCL failed to build program executable!\n");
        clGetProgramBuildInfo(program, device_id, CL_PROGRAM_BUILD_LOG, 0, NULL, &len);
        buffer=(char*) malloc(len);
        clGetProgramBuildInfo(program, device_id, CL_PROGRAM_BUILD_LOG, len, buffer, NULL);
        printf("%s\n", buffer);
        free(buffer);
        ENDRUN("exiting")
  }
  
  kick_krn = clCreateKernel(program, "kick_kernel", &err);
  if(err!= CL_SUCCESS) ENDRUN("OpenCL failed to create kick_kernel"); 
  timestep_krn = clCreateKernel(program, "timestep_kernel", &err);
  if(err!= CL_SUCCESS) ENDRUN("OpenCL failed to create timestep_kernel"); 
  potential_krn = clCreateKernel(program, "potential_kernel", &err);
  if(err!= CL_SUCCESS) ENDRUN("OpenCL failed to create potential_kernel"); 

  printf("OpenCL initialized using device %s\n", name);
  free(name);
}

void release_cl_buffers()
{
    clReleaseMemObject(_ipos);
    clReleaseMemObject(_jpos);
    clReleaseMemObject(_ivel);
    clReleaseMemObject(_jvel);
    clReleaseMemObject(_acc);
    clReleaseMemObject(_timestep);
    clReleaseMemObject(_pot);
    nalloc=0;
}


void init_cl_buffers(UINT nthread)
{
  cl_int err;

  if(nalloc) release_cl_buffers();

  _ipos = clCreateBuffer(context,  CL_MEM_READ_ONLY,  nthread*sizeof(CLFLOAT4), NULL, &err);
  if(err!=CL_SUCCESS) ENDRUN("failure creating cl buffer")
  _ivel = clCreateBuffer(context,  CL_MEM_READ_ONLY,  nthread*sizeof(CLFLOAT4), NULL, &err);
  if(err!=CL_SUCCESS) ENDRUN("failure creating cl buffer")
  _jpos = clCreateBuffer(context,  CL_MEM_READ_ONLY,  nthread*sizeof(CLFLOAT4), NULL, &err);
  if(err!=CL_SUCCESS) ENDRUN("failure creating cl buffer")
  _jvel = clCreateBuffer(context,  CL_MEM_READ_ONLY,  nthread*sizeof(CLFLOAT4), NULL, &err);
  if(err!=CL_SUCCESS) ENDRUN("failure creating cl buffer")

  _acc = clCreateBuffer(context, CL_MEM_WRITE_ONLY, nthread*sizeof(CLFLOAT4), NULL, &err);
  if(err!=CL_SUCCESS) ENDRUN("failure creating cl buffer")
  _pot = clCreateBuffer(context, CL_MEM_WRITE_ONLY, nthread*sizeof(CLFLOAT), NULL, &err);
  if(err!=CL_SUCCESS) ENDRUN("failure creating cl buffer")
  _timestep = clCreateBuffer(context, CL_MEM_WRITE_ONLY, nthread*sizeof(CLFLOAT), NULL, &err);
  if(err!=CL_SUCCESS) ENDRUN("failure creating cl buffer")
  
  nalloc=nthread;
}

void close_cl()
{
  release_cl_buffers();
  clReleaseProgram(program);
  clReleaseKernel(kick_krn);
  clReleaseKernel(timestep_krn);
  clReleaseKernel(potential_krn);
  clReleaseCommandQueue(queue);
  clReleaseContext(context);
}

void kick_cl(struct sys s1, struct sys s2, DOUBLE dt)
{
  cl_int err;
  size_t global, local;
  int i, n2;
  int groupsize,nthread,blocksize;
  CLFLOAT cleps2=(CLFLOAT) eps2;
  struct particle *ipart;
  
  n2=s2.n-s2.nzero;
  if(s1.n==0 || n2==0) return;

  groupsize=NTHREAD;
  if(s1.n < NTHREAD) groupsize=s1.n;
  if(groupsize<=8) groupsize=8; // needed for some implementations (e.g. pocl)
  nthread=((s1.n-1)/groupsize+1)*groupsize;
  blocksize=BLOCKSIZE;
  if(n2<blocksize) blocksize=n2;

  if(nthread>nalloc) init_cl_buffers(nthread);

  CLFLOAT4* ipos = (CLFLOAT4*) clEnqueueMapBuffer(queue, _ipos, CL_TRUE,  CL_MAP_READ, 0, nthread * sizeof(CLFLOAT4), 0, NULL, NULL, &err);
  if(err!=CL_SUCCESS) ENDRUN("clEnqueueMapBuffer fail")
  CLFLOAT4* jpos = (CLFLOAT4*) clEnqueueMapBuffer(queue, _jpos, CL_TRUE,  CL_MAP_READ, 0, nthread * sizeof(CLFLOAT4), 0, NULL, NULL, &err);
  if(err!=CL_SUCCESS) ENDRUN("clEnqueueMapBuffer fail")

  for(i=0;i<s1.n;i++)
  {
    ipart=GETPART(s1,i);
    ipos[i].s0=ipart->pos[0];
    ipos[i].s1=ipart->pos[1];
    ipos[i].s2=ipart->pos[2];
    ipos[i].s3=ipart->mass;
  }  
  for(i=s1.n;i<nthread;i++) ipos[i]=(CLFLOAT4) {{0.0,0.0,0.0,0.0}};  
  for(i=0;i<n2;i++)
  {
    ipart=GETPART(s2,i);
    jpos[i].s0=ipart->pos[0];
    jpos[i].s1=ipart->pos[1];
    jpos[i].s2=ipart->pos[2];
    jpos[i].s3=ipart->mass;
  }
  //~ for(i=0;i<nthread;i++) acc[i]=(CLFLOAT4) {{0.0,0.0,0.0,0.0}};  
  
  err=clEnqueueUnmapMemObject(queue, _ipos, (void*)ipos, 0, NULL, NULL);
  if(err!=CL_SUCCESS) ENDRUN("clEnqueueUnmapBuffer fail")
  err=clEnqueueUnmapMemObject(queue, _jpos, (void*)jpos, 0, NULL, NULL);
  if(err!=CL_SUCCESS) ENDRUN("clEnqueueUnmapBuffer fail")
                        
  err=0;
  err |= clSetKernelArg(kick_krn, 0, sizeof(unsigned int), &n2);                      
  err |= clSetKernelArg(kick_krn, 1, sizeof(int), &blocksize);                      
  err |= clSetKernelArg(kick_krn, 2, sizeof(CLFLOAT), &cleps2);                      
  err |= clSetKernelArg(kick_krn, 3, sizeof(cl_mem), &_ipos);                      
  err |= clSetKernelArg(kick_krn, 4, sizeof(cl_mem), &_jpos);                      
  err |= clSetKernelArg(kick_krn, 5, sizeof(cl_mem), &_acc);                      
  err |= clSetKernelArg(kick_krn, 6, blocksize*sizeof(CLFLOAT4), NULL);                      

  if(err!=CL_SUCCESS) ENDRUN("clSetKernelArg fail")

  global=nthread;
  local=groupsize;
  err=clEnqueueNDRangeKernel(queue, kick_krn, 1, NULL, &global, &local, 0, NULL, NULL);

  if(err!=CL_SUCCESS) ENDRUN("clEnqueueNDRangeKernel fail")

  clFinish(queue);

  CLFLOAT4* acc = (CLFLOAT4*) clEnqueueMapBuffer(queue, _acc, CL_TRUE,  CL_MAP_WRITE, 0, nthread * sizeof(CLFLOAT4), 0, NULL, NULL, &err);
  if(err!=CL_SUCCESS) ENDRUN("clEnqueueMapBuffer fail")

  for(i=0;i<s1.n;i++)
  {
    ipart=GETPART(s1,i);
    COMPSUMV(ipart->vel[0],ipart->vel_e[0],dt*acc[i].s0);
    COMPSUMV(ipart->vel[1],ipart->vel_e[1],dt*acc[i].s1);
    COMPSUMV(ipart->vel[2],ipart->vel_e[2],dt*acc[i].s2);
  }

  err=clEnqueueUnmapMemObject(queue, _acc, (void*)acc, 0, NULL, NULL);
  if(err!=CL_SUCCESS) ENDRUN("clEnqueueUnmapBuffer fail")

}

void _timestep_cl(struct sys s1, struct sys s2,int dir)
{
  cl_int err;
  size_t global, local;
  int i, n2;
  int groupsize,nthread,blocksize;
  CLFLOAT cleps2=(CLFLOAT) eps2;
  CLFLOAT cldtparam=(CLFLOAT) dt_param;
  struct particle *ipart;
    
  n2=s2.n;
  if(s1.n==0 || n2==0) return;

  groupsize=NTHREAD;
  if(s1.n < NTHREAD) groupsize=s1.n;
  if(groupsize<=8) groupsize=8; // needed for some implementations (e.g. pocl)
  nthread=((s1.n-1)/groupsize+1)*groupsize;
  blocksize=BLOCKSIZE;
  if(n2<blocksize) blocksize=n2;

  if(nthread>nalloc) init_cl_buffers(nthread);

  CLFLOAT4* ipos = (CLFLOAT4*) clEnqueueMapBuffer(queue, _ipos, CL_TRUE,  CL_MAP_READ, 0, nthread * sizeof(CLFLOAT4), 0, NULL, NULL, &err);
  if(err!=CL_SUCCESS) ENDRUN("clEnqueueMapBuffer fail")
  CLFLOAT4* jpos = (CLFLOAT4*) clEnqueueMapBuffer(queue, _jpos, CL_TRUE,  CL_MAP_READ, 0, nthread * sizeof(CLFLOAT4), 0, NULL, NULL, &err);
  if(err!=CL_SUCCESS) ENDRUN("clEnqueueMapBuffer fail")
  CLFLOAT4* ivel = (CLFLOAT4*) clEnqueueMapBuffer(queue, _ivel, CL_TRUE,  CL_MAP_READ, 0, nthread * sizeof(CLFLOAT4), 0, NULL, NULL, &err);
  if(err!=CL_SUCCESS) ENDRUN("clEnqueueMapBuffer fail")
  CLFLOAT4* jvel = (CLFLOAT4*) clEnqueueMapBuffer(queue, _jvel, CL_TRUE,  CL_MAP_READ, 0, nthread * sizeof(CLFLOAT4), 0, NULL, NULL, &err);
  if(err!=CL_SUCCESS) ENDRUN("clEnqueueMapBuffer fail")

  for(i=0;i<s1.n;i++)
  {
    ipart=GETPART(s1,i);
    ipos[i].s0=ipart->pos[0]; ivel[i].s0=ipart->vel[0];
    ipos[i].s1=ipart->pos[1]; ivel[i].s1=ipart->vel[1];
    ipos[i].s2=ipart->pos[2]; ivel[i].s2=ipart->vel[2];
    ipos[i].s3=ipart->mass;   ivel[i].s3=0.;
  }  
  for(i=s1.n;i<nthread;i++) ipos[i]=(CLFLOAT4) {{0.0,0.0,0.0,0.0}};  
  for(i=s1.n;i<nthread;i++) ivel[i]=(CLFLOAT4) {{0.0,0.0,0.0,0.0}};  
  for(i=0;i<n2;i++)
  {
    ipart=GETPART(s2,i);
    jpos[i].s0= ipart->pos[0]; jvel[i].s0=ipart->vel[0];
    jpos[i].s1= ipart->pos[1]; jvel[i].s1=ipart->vel[1];
    jpos[i].s2= ipart->pos[2]; jvel[i].s2=ipart->vel[2];
    jpos[i].s3= ipart->mass;   jvel[i].s3=0.;
  }
  //~ for(i=0;i<nthread;i++) timestep[i]=(CLFLOAT) 0.;  

  err=clEnqueueUnmapMemObject(queue, _ipos, (void*)ipos, 0, NULL, NULL);
  if(err!=CL_SUCCESS) ENDRUN("clEnqueueUnmapBuffer fail")
  err=clEnqueueUnmapMemObject(queue, _jpos, (void*)jpos, 0, NULL, NULL);
  if(err!=CL_SUCCESS) ENDRUN("clEnqueueUnmapBuffer fail")
  err=clEnqueueUnmapMemObject(queue, _ivel, (void*)ivel, 0, NULL, NULL);
  if(err!=CL_SUCCESS) ENDRUN("clEnqueueUnmapBuffer fail")
  err=clEnqueueUnmapMemObject(queue, _jvel, (void*)jvel, 0, NULL, NULL);
  if(err!=CL_SUCCESS) ENDRUN("clEnqueueUnmapBuffer fail")
                        
  err=0;
  err |= clSetKernelArg(timestep_krn, 0, sizeof(unsigned int), &n2);                      
  err |= clSetKernelArg(timestep_krn, 1, sizeof(int), &blocksize);                      
  err |= clSetKernelArg(timestep_krn, 2, sizeof(CLFLOAT), &cleps2);                      
  err |= clSetKernelArg(timestep_krn, 3, sizeof(CLFLOAT), &cldtparam);                      
  err |= clSetKernelArg(timestep_krn, 4, sizeof(cl_mem), &_ipos);                      
  err |= clSetKernelArg(timestep_krn, 5, sizeof(cl_mem), &_ivel);                      
  err |= clSetKernelArg(timestep_krn, 6, sizeof(cl_mem), &_jpos);                      
  err |= clSetKernelArg(timestep_krn, 7, sizeof(cl_mem), &_jvel);                      
  err |= clSetKernelArg(timestep_krn, 8, sizeof(cl_mem), &_timestep);                      
  err |= clSetKernelArg(timestep_krn, 9, blocksize*sizeof(CLFLOAT4), NULL);                      
  err |= clSetKernelArg(timestep_krn, 10, blocksize*sizeof(CLFLOAT4), NULL);                      
  err |= clSetKernelArg(timestep_krn, 11, sizeof(int), &dir);                      

  if(err!=CL_SUCCESS) ENDRUN("clSetKernelArg fail")

  global=nthread;
  local=groupsize;
  err=clEnqueueNDRangeKernel(queue, timestep_krn, 1, NULL, &global, &local, 0, NULL, NULL);

  if(err!=CL_SUCCESS) ENDRUN("kernel run fail")

  clFinish(queue);
  
  CLFLOAT* timestep = (CLFLOAT*) clEnqueueMapBuffer(queue, _timestep, CL_TRUE,  CL_MAP_WRITE, 0, nthread * sizeof(CLFLOAT), 0, NULL, NULL, &err);
  if(err!=CL_SUCCESS) ENDRUN("clEnqueueMapBuffer fail")

  for(i=0;i<s1.n;i++) 
  {
    ipart=GETPART(s1,i);
    ipart->timestep=timestep[i]; //if(timestep[i]<ipart->timestep)
  }

  err=clEnqueueUnmapMemObject(queue, _timestep, (void*)timestep, 0, NULL, NULL);
  if(err!=CL_SUCCESS) ENDRUN("clEnqueueUnmapBuffer fail")

}

void timestep_cl(struct sys s1, struct sys s2,int dir)
{

  if(accel_zero_mass && s1.nzero*s2.nzero>CLWORKLIMIT/4)
  {
    struct sys s1m=zerosys, s1ml=zerosys, s2m=zerosys;
    
    s1m.n=s1.n-s1.nzero;
    if(s1m.n>0) s1m.part=s1.part;
    
    s1ml.n=s1.nzero;
    s1ml.nzero=s1.nzero;
    if(s1ml.n>0) s1ml.part=s1.zeropart;
    if(s1ml.n>0) s1ml.zeropart=s1.zeropart;
    
    s2m.n=s2.n-s2.nzero;
    if(s2m.n>0) s2m.part=s2.part;    
    _timestep_cl(s1m,  s2, dir);
    _timestep_cl(s1ml, s2m, dir);
  } else
  {
    _timestep_cl(s1,s2,dir);
  }
}


void potential_cl(struct sys s1, struct sys s2)
{
  cl_int err;
  size_t global, local;
  int i, n2;
  int groupsize,nthread,blocksize;
  CLFLOAT cleps2=(CLFLOAT) eps2;
  struct particle *ipart;
  
  n2=s2.n-s2.nzero;
  if(s1.n==0 || n2==0) return;

  groupsize=NTHREAD;
  if(s1.n < NTHREAD) groupsize=s1.n;
  if(groupsize<=8) groupsize=8; // needed for some implementations (e.g. pocl)
  nthread=((s1.n-1)/groupsize+1)*groupsize;
  blocksize=BLOCKSIZE;
  if(n2<blocksize) blocksize=n2;

  if(nthread>nalloc) init_cl_buffers(nthread);

  CLFLOAT4* ipos = (CLFLOAT4*) clEnqueueMapBuffer(queue, _ipos, CL_TRUE,  CL_MAP_READ, 0, nthread * sizeof(CLFLOAT4), 0, NULL, NULL, &err);
  if(err!=CL_SUCCESS) ENDRUN("clEnqueueMapBuffer fail")
  CLFLOAT4* jpos = (CLFLOAT4*) clEnqueueMapBuffer(queue, _jpos, CL_TRUE,  CL_MAP_READ, 0, nthread * sizeof(CLFLOAT4), 0, NULL, NULL, &err);
  if(err!=CL_SUCCESS) ENDRUN("clEnqueueMapBuffer fail")

  for(i=0;i<s1.n;i++)
  {
    ipart=GETPART(s1,i);
    ipos[i].s0=ipart->pos[0];
    ipos[i].s1=ipart->pos[1];
    ipos[i].s2=ipart->pos[2];
    ipos[i].s3=ipart->mass;
  }  
  for(i=s1.n;i<nthread;i++) ipos[i]=(CLFLOAT4) {{0.0,0.0,0.0,0.0}};  
  for(i=0;i<n2;i++)
  {
    ipart=GETPART(s2,i);
    jpos[i].s0=ipart->pos[0];
    jpos[i].s1=ipart->pos[1];
    jpos[i].s2=ipart->pos[2];
    jpos[i].s3=ipart->mass;
  }
  //~ for(i=0;i<nthread;i++) pot[i]=0.0;  

  err=clEnqueueUnmapMemObject(queue, _ipos, (void*)ipos, 0, NULL, NULL);  
  if(err!=CL_SUCCESS) ENDRUN("clEnqueueUnmapBuffer fail")
  err=clEnqueueUnmapMemObject(queue, _jpos, (void*)jpos, 0, NULL, NULL);
  if(err!=CL_SUCCESS) ENDRUN("clEnqueueUnmapBuffer fail")
                        
  err=0;
  err |= clSetKernelArg(potential_krn, 0, sizeof(unsigned int), &n2);                      
  err |= clSetKernelArg(potential_krn, 1, sizeof(int), &blocksize);                      
  err |= clSetKernelArg(potential_krn, 2, sizeof(CLFLOAT), &cleps2);                      
  err |= clSetKernelArg(potential_krn, 3, sizeof(cl_mem), &_ipos);                      
  err |= clSetKernelArg(potential_krn, 4, sizeof(cl_mem), &_jpos);                      
  err |= clSetKernelArg(potential_krn, 5, sizeof(cl_mem), &_pot);                      
  err |= clSetKernelArg(potential_krn, 6, blocksize*sizeof(CLFLOAT4), NULL);                      

  if(err!=CL_SUCCESS) ENDRUN("clSetKernelArg fail")

  global=nthread;
  local=groupsize;
  err=clEnqueueNDRangeKernel(queue, potential_krn, 1, NULL, &global, &local, 0, NULL, NULL);

  if(err!=CL_SUCCESS) ENDRUN("clEnqueueNDRangeKernel fail")

  err=clFinish(queue);
  if(err!=CL_SUCCESS) ENDRUN("clFinish fail")

  CLFLOAT* pot = (CLFLOAT*) clEnqueueMapBuffer(queue, _pot, CL_TRUE,  CL_MAP_WRITE, 0, nthread * sizeof(CLFLOAT), 0, NULL, NULL, &err);
  if(err!=CL_SUCCESS) ENDRUN("clEnqueueMapBuffer fail")

  for(i=0;i<s1.n;i++) GETPART(s1,i)->pot+=pot[i];

  err=clEnqueueUnmapMemObject(queue, _pot, (void*)pot, 0, NULL, NULL);
  if(err!=CL_SUCCESS) ENDRUN("clEnqueueUnmapBuffer fail")

}

#endif
