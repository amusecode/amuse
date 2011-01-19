#include <pthread.h>
#include "Main_glsimpleX.h"

using namespace std;

int main(int argc, char** argv) {  

  pthread_t viz_thread;

  cerr.precision(5);
  
  MPI::Init( argc, argv );

  static SimpleX SimpleXGrid;

  pthread_create(&viz_thread,NULL,viewer, (void*) &SimpleXGrid);

  main_loop(argc, argv, &SimpleXGrid);

  pthread_join(viz_thread, NULL);
  
  MPI::Finalize();
  
  pthread_exit(NULL);

  return EXIT_SUCCESS;

}
