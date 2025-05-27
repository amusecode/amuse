#include <functions.h>
#include <omp.h>

HostError NextParticles(unsigned int N, unsigned int ompthreads, unsigned int* counter, unsigned int* vetint, double ATIME, double *local_time, double* step, int* next, unsigned long* nextsize){

   *nextsize = 0;

#pragma omp parallel
{
      unsigned int cpu_thread_id = omp_get_thread_num();
      unsigned int istart = cpu_thread_id*N/ompthreads;
      counter[cpu_thread_id] = 0;

#pragma omp for
      for(unsigned int i = 0; i < N; i++){
         if( (local_time[i] + step[i]) == ATIME){
            vetint[counter[cpu_thread_id] + istart] = i;
            counter[cpu_thread_id]++;
         }
      }

#pragma omp barrier
      unsigned int from = 0;
      unsigned int to = 0;

      for(unsigned int i = cpu_thread_id; i > 0; i--)
         from += counter[i-1];

      to = from + counter[cpu_thread_id];

      for(unsigned int i = from; i < to; i++)
         next[i] = vetint[istart + i - from];
}

      for(unsigned int i = 0; i < ompthreads; i++)
         *nextsize += counter[i];


   return HNoError;

   }

