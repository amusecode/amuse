#include <assert.h>
//#include <cutil.h>
//

template <typename T, bool PINNED>
struct cuVector
{
  T *dev_pointer;
  T *host_pointer;
  int size;
  cuVector()
  {
    dev_pointer = NULL;
    host_pointer = NULL;
    size = 0;
  }
  void init(const T *data, const int _size)
  {
    assert(_size > 0);
    cufree();
    allocate(_size);
    assert(size == _size);
    CUDA_SAFE_CALL(cudaMemcpy(host_pointer, data, size*sizeof(T), cudaMemcpyHostToHost));
  }
  void copy(const T *data, const int _size)
  {
    assert(_size <= size);
    CUDA_SAFE_CALL(cudaMemcpy(host_pointer, data, _size*sizeof(T), cudaMemcpyHostToHost));
  }

  ~cuVector()
  {
    cufree();
  }
  void allocate(int _size)
  {
    assert(_size > 0);
    size = _size;
    void *p;
    CUDA_SAFE_CALL(cudaMalloc(&p, size * sizeof(T)));
    assert(p);
    dev_pointer = (T*)p;
    if (PINNED) 
    {
      CUDA_SAFE_CALL(cudaMallocHost(&p, size * sizeof(T)));
    }
    else 
      p = malloc(size*sizeof(T));
    assert(p);
    host_pointer = (T*)p;
  }
  void cufree()
  {
    if (size > 0)
    {
      CUDA_SAFE_CALL(cudaFree(dev_pointer));
      if (PINNED)   
      {
        CUDA_SAFE_CALL(cudaFreeHost(host_pointer));
      }
      else          
        free(host_pointer);
      dev_pointer = NULL;
      host_pointer = NULL;
      size = 0;
    }
  }
  void h2d(int count)
  {
    assert(count <= size);
    if (size > 0)
      CUDA_SAFE_CALL(cudaMemcpy(dev_pointer, host_pointer, count * sizeof(T), cudaMemcpyHostToDevice));
  }
  void h2d()
  {
    if (size > 0)
      this->h2d(size);
  }
  void d2h(int count)
  {
    assert(count <= size);
    if (size > 0)
      CUDA_SAFE_CALL(cudaMemcpy(host_pointer, dev_pointer, count * sizeof(T), cudaMemcpyDeviceToHost));
  }
  void d2h()
  {
    if (size > 0)
      this->d2h(size);
  }
  T &operator [] (int i)
  {
    return host_pointer[i];
  }
  operator T* ()
  {
    return dev_pointer;
  }

  T* d() {return dev_pointer;}
};

