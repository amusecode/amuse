#include <pthread.h>

#define __VIEWER 1
#include "viewer.cpp"

pthread_t viz_thread;

int start_viewer(){
pthread_t viz_thread;

pthread_create(&viz_thread,NULL,viewer, (void*) SimpleXGrid);

return 0;

}

int stop_viewer(){

pthread_join(viz_thread, NULL);

return 0;

}
