#include <assert.h>
// #include <cutil.h>

template <typename T>
struct cudaPointer{
  T *dev_pointer;
  T *host_pointer;
  int size;
  cudaPointer(){
    dev_pointer = NULL;
    host_pointer = NULL;
    size = 0;
  }
  //  ~cudaPointer(){
  // free();
  //  }
  void allocate(int _size){
    size = _size;
    void *p;
    CUDA_SAFE_CALL(cudaMalloc(&p, size * sizeof(T)));
    assert(p);
    dev_pointer = (T*)p;
    CUDA_SAFE_CALL(cudaMallocHost(&p, size * sizeof(T)));
    assert(p);
    host_pointer = (T*)p;
  }
  void free(){
    CUDA_SAFE_CALL(cudaFree(dev_pointer));
    CUDA_SAFE_CALL(cudaFreeHost(host_pointer));
    dev_pointer = NULL;
    host_pointer = NULL;
    size = 0;
  }
  void htod(int count){
    CUDA_SAFE_CALL(cudaMemcpy(dev_pointer, host_pointer, count * sizeof(T), cudaMemcpyHostToDevice));
  }
  void htod(){
    this->htod(size);
  }
  void dtoh(int count){
    CUDA_SAFE_CALL(cudaMemcpy(host_pointer, dev_pointer, count * sizeof(T), cudaMemcpyDeviceToHost));
  }
  void dtoh(){
    this->dtoh(size);
  }
  T &operator [] (int i){
    return host_pointer[i];
  }
  operator T* (){
    return dev_pointer;
  }
};
