
#include "Main_simpleX.h"

using namespace std;

int main(int argc, char** argv) {  

  cerr.precision(5);
  
  MPI::Init( argc, argv );

  SimpleX SimpleXGrid;

  main_loop(argc, argv, &SimpleXGrid);

  MPI::Finalize();
  
  return EXIT_SUCCESS;

}
