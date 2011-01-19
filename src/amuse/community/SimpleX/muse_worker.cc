#include <mpi.h>
#include <iostream>
#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include "muse_worker.h"


class message_header {
public:
  int tag;
  int len;
  int number_of_doubles;
  int number_of_ints;
  int number_of_floats;
  int number_of_strings;
  int number_of_bools;
  int number_of_longs;
  message_header(): tag(0), len(1), number_of_doubles(0), number_of_ints(0), number_of_floats(0), number_of_strings(0), number_of_bools(0), number_of_longs(0){}
  void send(MPI::Intercomm & intercomm, int rank){
    int header[8];
    header[0] = tag;
    header[1] = len;
    header[2] =number_of_doubles;
    header[3] =number_of_ints;
    header[4] =number_of_floats;
    header[5] =number_of_strings;
    header[6] =number_of_bools;
    header[7] =number_of_longs;
    intercomm.Send(header, 8, MPI_INT, 0, 999);
  }
  void recv(MPI::Intercomm & intercomm, int rank) {
    int header[8];
    intercomm.Bcast(header, 8, MPI_INT, 0);
    tag = header[0];
    len = header[1];
    number_of_doubles=header[2];
    number_of_ints=header[3];
    number_of_floats=header[4];
    number_of_strings=header[5];
    number_of_bools=header[6];
    number_of_longs=header[7];
  }
};
void onexit(void) {
    int flag = 0;
    MPI_Finalized(&flag);
    
    if(!flag) {
        MPI::Intercomm parent = MPI::COMM_WORLD.Get_parent();
        int rank = parent.Get_rank();
        
        message_header reply_header;
        reply_header.tag = -2;
        reply_header.send(parent, rank);
        
        parent.Disconnect();
        MPI_Finalize();
    }
}



int internal__redirect_outputs(const char * stdoutfile, const char * stderrfile)
{
    char fullname[1024];
    int mpi_rank, mpi_err;
    
    fclose(stdin);
    
    mpi_err = MPI_Comm_rank( MPI_COMM_WORLD, &mpi_rank );
    if (strcmp(stdoutfile, "none"))
    {
        if (strcmp(stdoutfile, "/dev/null") ) 
        {
            snprintf( fullname, 1024, "%s.%.3d",stdoutfile, mpi_rank );
        }
        else
        {
            snprintf( fullname, 1024, "%s",stdoutfile);
        }
        freopen(fullname, "a+", stdout);
    }
    
    if (strcmp(stderrfile, "none"))
    {         
        if (strcmp(stderrfile, "/dev/null") ) 
        {
            snprintf( fullname, 1024, "%s.%.3d",stderrfile, mpi_rank ); 
        }
        else
        {
            snprintf( fullname, 1024, "%s",stderrfile);
        }
        freopen(fullname, "a+", stderr);
    }
        
    return 0;
}


void run_loop() {
  
  MPI::Intercomm parent = MPI::COMM_WORLD.Get_parent();
  int rank = parent.Get_rank();
  
  bool must_run_loop = true;
  char * characters = 0, * output_characters = 0;
  
  int max_len = 10;
  int * ints_in = new int[ max_len * 1];
  int * ints_out = new int[ max_len * 2];
  
  int * strings_in = new int[ max_len * 2];
  
  double * doubles_in = new double[ max_len * 6];
  double * doubles_out = new double[ max_len * 6];
  
  
  
  
  
  while(must_run_loop) {
    
    message_header request_header;
    message_header reply_header;
    
    request_header.recv(parent,rank);
    if (request_header.len > max_len) {
      max_len = request_header.len + 255;
      delete[] ints_in;
      delete[] ints_out;
      delete[] strings_in;
      delete[] doubles_in;
      delete[] doubles_out;
      ints_in = new int[ max_len * 1];
      ints_out = new int[ max_len * 2];
      
      strings_in = new int[ max_len * 2];
      
      doubles_in = new double[ max_len * 6];
      doubles_out = new double[ max_len * 6];
      
      
      
      
    }
    if(request_header.number_of_doubles > 0) {
      parent.Bcast(doubles_in, request_header.number_of_doubles * request_header.len, MPI_DOUBLE, 0);
    }
    if(request_header.number_of_ints > 0) {
      parent.Bcast(ints_in, request_header.number_of_ints * request_header.len, MPI_INT, 0);
    }
    
    if(request_header.number_of_strings > 0) {
      parent.Bcast(strings_in, request_header.number_of_strings * request_header.len, MPI_INTEGER, 0);
      characters = new char[strings_in[request_header.number_of_strings * request_header.len - 1] + 1];
      parent.Bcast(characters,  strings_in[request_header.number_of_strings * request_header.len- 1] + 1, MPI_CHARACTER, 0);
    }
    
    
    
    reply_header.tag = request_header.tag;
    
    reply_header.len = request_header.len;
    
    switch(request_header.tag) {
      case 0:
        must_run_loop = false;
        break;
      case 111364973:
        ints_out[0] = evolve(
          doubles_in[0] ,
          ints_in[0]
        );
        reply_header.number_of_ints = 1;
        break;
      
      case 572892820:
        for (int i = 0 ; i < request_header.len; i++){
          ints_out[i] = remove_particle(
            ints_in[i]
          );
        }
        reply_header.number_of_ints = 1;
        break;
      
      case 634707820:
        ints_out[0] = cleanup_module();
        reply_header.number_of_ints = 1;
        break;
      
      case 670175614:
        ints_out[0] = initialize(
          doubles_in[0]
        );
        reply_header.number_of_ints = 1;
        break;
      
      case 967950880:
        for (int i = 0 ; i < request_header.len; i++){
          get_state(
            ints_in[i] ,
            &doubles_out[i] ,
            &doubles_out[( 1 * request_header.len) + i] ,
            &doubles_out[( 2 * request_header.len) + i] ,
            &doubles_out[( 3 * request_header.len) + i] ,
            &doubles_out[( 4 * request_header.len) + i] ,
            &doubles_out[( 5 * request_header.len) + i]
          );
        }
        reply_header.number_of_doubles = 6;
        break;
      
      case 1141573512:
        ints_out[0] = internal__redirect_outputs(
          characters + ( 0- 1 < 0 ? 0 :strings_in[0 - 1] + 1) ,
          characters + ( 1- 1 < 0 ? 0 :strings_in[1 - 1] + 1)
        );
        reply_header.number_of_ints = 1;
        break;
      
      case 1449106038:
        ints_out[0] = reinitialize();
        reply_header.number_of_ints = 1;
        break;
      
      case 1859021207:
        for (int i = 0 ; i < request_header.len; i++){
          ints_out[i] = add_particle(
            &ints_out[( 1 * request_header.len) + i] ,
            doubles_in[i] ,
            doubles_in[( 1 * request_header.len) + i] ,
            doubles_in[( 2 * request_header.len) + i] ,
            doubles_in[( 3 * request_header.len) + i] ,
            doubles_in[( 4 * request_header.len) + i] ,
            doubles_in[( 5 * request_header.len) + i]
          );
        }
        reply_header.number_of_ints = 2;
        break;
      
      case 2009103572:
        ints_out[0] = setup_module();
        reply_header.number_of_ints = 1;
        break;
      
      default:
        reply_header.tag = -1;
    }
    
    MPI::COMM_WORLD.Barrier();
    
    
    if(rank == 0) {
      
      
      reply_header.send(parent, rank);
      
      
      if(reply_header.number_of_doubles > 0) {
        parent.Send(doubles_out, reply_header.number_of_doubles * request_header.len, MPI_DOUBLE, 0, 999);
      }
      if(reply_header.number_of_ints > 0) {
        parent.Send(ints_out, reply_header.number_of_ints * request_header.len, MPI_INT, 0, 999);
      }
    
    }
    if (characters) { delete[] characters; characters = 0;}
    if (output_characters) { delete[] output_characters; output_characters = 0;}
  }
  delete[] ints_in;
  delete[] ints_out;
  delete[] strings_in;
  delete[] doubles_in;
  delete[] doubles_out;
  
  parent.Disconnect();
}

int main(int argc, char *argv[])
{
  
  MPI::Init_thread(argc, argv, MPI_THREAD_MULTIPLE);
  atexit(onexit);
  
  run_loop();
  
  MPI_Finalize();
  return 0;
}