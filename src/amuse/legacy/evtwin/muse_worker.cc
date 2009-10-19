#include <mpi.h>
#include "muse_worker.h"


class message_header {
public:
  int tag;
  int len;
  int number_of_doubles;
  int number_of_ints;
  int number_of_floats;
  int number_of_chars;
  message_header(): tag(0), len(1), number_of_doubles(0), number_of_ints(0), number_of_floats(0), number_of_chars(0){}
  void send(MPI::Intercomm & intercomm, int rank){
    int header[6];
    header[0] = tag;
    header[1] = len;
    header[2] =number_of_doubles;
    header[3] =number_of_ints;
    header[4] =number_of_floats;
    header[5] =number_of_chars;
    intercomm.Send(header, 6, MPI_INT, 0, 999);
  }
  void recv(MPI::Intercomm & intercomm, int rank) {
    int header[6];
    intercomm.Bcast(header, 6, MPI_INT, 0);
    tag = header[0];
    len = header[1];
    number_of_doubles=header[2];
    number_of_ints=header[3];
    number_of_floats=header[4];
    number_of_chars=header[5];
  }
};

void run_loop() {
  
  MPI::Intercomm parent = MPI::COMM_WORLD.Get_parent();
  int rank = parent.Get_rank();
  
  bool must_run_loop = true;
  
  int max_len = 10;
  int * ints_in = new int[ max_len * 255];
  int * ints_out = new int[ max_len * 255];
  char * chars_in = new char[ max_len * 255];
  char * chars_out = new char[ max_len * 255];
  float * floats_in = new float[ max_len * 255];
  float * floats_out = new float[ max_len * 255];
  double * doubles_in = new double[ max_len * 255];
  double * doubles_out = new double[ max_len * 255];
  
  while(must_run_loop) {
    
    message_header request_header;
    message_header reply_header;
    
    request_header.recv(parent,rank);
    if (request_header.len > max_len) {
      max_len = request_header.len + 255;
      delete ints_in;
      delete ints_out;
      delete chars_in;
      delete chars_out;
      delete floats_in;
      delete floats_out;
      delete doubles_in;
      delete doubles_out;
      ints_in = new int[ max_len * 255];
      ints_out = new int[ max_len * 255];
      chars_in = new char[ max_len * 255];
      chars_out = new char[ max_len * 255];
      floats_in = new float[ max_len * 255];
      floats_out = new float[ max_len * 255];
      doubles_in = new double[ max_len * 255];
      doubles_out = new double[ max_len * 255];
    }
    if(request_header.number_of_doubles > 0) {
      parent.Bcast(doubles_in, request_header.number_of_doubles * request_header.len, MPI_DOUBLE, 0);
    }
    if(request_header.number_of_ints > 0) {
      parent.Bcast(ints_in, request_header.number_of_ints * request_header.len, MPI_INT, 0);
    }
    if(request_header.number_of_floats > 0) {
      parent.Bcast(floats_in, request_header.number_of_floats * request_header.len, MPI_FLOAT, 0);
    }
    if(request_header.number_of_chars > 0) {
      parent.Bcast(chars_in, request_header.number_of_chars * request_header.len, MPI_CHARACTER, 0);
    }
    
    reply_header.tag = request_header.tag;
    
    reply_header.len = request_header.len;
    
    switch(request_header.tag) {
      case 0:
        must_run_loop = false;
        break;
      case 92616527:
        doubles_out[0] = wrapper_get_age(
          ints_in[0]
        );
        reply_header.number_of_doubles = 1;
        break;
      case 179224572:
        doubles_out[0] = wrapper_get_mass(
          ints_in[0]
        );
        reply_header.number_of_doubles = 1;
        break;
      case 422430783:
        doubles_out[0] = wrapper_get_radius(
          ints_in[0]
        );
        reply_header.number_of_doubles = 1;
        break;
      case 543183209:
        wrapper_set_init_run_name(
          chars_in
        );
        break;
      case 871394194:
        ints_out[0] = wrapper_load_zams_star(
          doubles_in[0] ,
          doubles_in[1]
        );
        reply_header.number_of_ints = 1;
        break;
      case 1240315773:
        doubles_out[0] = wrapper_get_temperature(
          ints_in[0]
        );
        reply_header.number_of_doubles = 1;
        break;
      case 1350852878:
        ints_out[0] = wrapper_initialise_twin(
          chars_in ,
          ints_in[0] ,
          doubles_in[0]
        );
        reply_header.number_of_ints = 1;
        break;
      case 1419185947:
        wrapper_set_init_dat_name(
          chars_in
        );
        break;
      case 1572939278:
        wrapper_flush_star();
        break;
      case 1613624722:
        wrapper_swap_out_star(
          ints_in[0]
        );
        break;
      case 1620580413:
        wrapper_select_star(
          ints_in[0]
        );
        break;
      case 1741206111:
        ints_out[0] = wrapper_twin_evolve();
        reply_header.number_of_ints = 1;
        break;
      case 1918979957:
        wrapper_swap_in_star(
          ints_in[0]
        );
        break;
      case 1945501403:
        doubles_out[0] = wrapper_get_luminosity(
          ints_in[0]
        );
        reply_header.number_of_doubles = 1;
        break;
      case 2056389809:
        ints_out[0] = wrapper_get_stellar_type(
          ints_in[0]
        );
        reply_header.number_of_ints = 1;
        break;
      default:
        reply_header.tag = -1;
    }
    
    reply_header.send(parent, rank);
    if(reply_header.number_of_doubles > 0) {
      parent.Send(doubles_out, reply_header.number_of_doubles * request_header.len, MPI_DOUBLE, 0, 999);
    }
    if(reply_header.number_of_ints > 0) {
      parent.Send(ints_out, reply_header.number_of_ints * request_header.len, MPI_INT, 0, 999);
    }
    if(reply_header.number_of_floats > 0) {
      parent.Send(floats_out, reply_header.number_of_floats * request_header.len, MPI_FLOAT, 0, 999);
    }
    if(reply_header.number_of_chars > 0) {
      parent.Send(chars_out, reply_header.number_of_chars * request_header.len, MPI_CHARACTER, 0, 999);
    }
  }
  delete ints_in;
  delete ints_out;
  delete chars_in;
  delete chars_out;
  delete floats_in;
  delete floats_out;
  delete doubles_in;
  delete doubles_out;
  
  parent.Free();
}

int main(int argc, char *argv[])
{
  
  int provided;
  MPI::Init_thread(argc, argv, MPI_THREAD_MULTIPLE);
  
  run_loop();
  
  MPI_Finalize();
  return 0;
}
