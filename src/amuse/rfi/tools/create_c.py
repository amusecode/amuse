from amuse.support.core import late
from amuse.support import exceptions
from amuse.rfi.tools.create_code import GenerateASourcecodeString
from amuse.rfi.tools.create_code import GenerateASourcecodeStringFromASpecificationClass
from amuse.rfi.tools.create_code import DTypeSpec
from amuse.rfi.tools.create_code import dtypes
from amuse.rfi.tools.create_code import DTypeToSpecDictionary
from amuse.rfi.tools import create_definition
from amuse.rfi.core import LegacyFunctionSpecification

dtype_to_spec = DTypeToSpecDictionary({
    'int32' : DTypeSpec('ints_in', 'ints_out',
                    'HEADER_INTEGER_COUNT', 'int', 'MPI_INT'),
    'int64' : DTypeSpec('longs_in', 'longs_out',
                    'HEADER_LONG_COUNT', 'long long int', 'MPI_LONG_LONG_INT'),
    'float32' : DTypeSpec('floats_in', 'floats_out',
                    'HEADER_FLOAT_COUNT', 'float', 'MPI_FLOAT'),
    'float64' : DTypeSpec('doubles_in', 'doubles_out',
                    'HEADER_DOUBLE_COUNT', 'double', 'MPI_DOUBLE'),
    'bool' : DTypeSpec('booleans_in', 'booleans_out',
                    'HEADER_BOOLEAN_COUNT', 'int', 'MPI_INTEGER'),
    'string' : DTypeSpec('strings_in', 'strings_out',
                    'HEADER_STRING_COUNT', 'int', 'MPI_INTEGER'),
})

HEADER_CODE_STRING = """
#ifndef NOMPI
    #include <mpi.h>
#endif
#include <iostream>
#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include <sys/socket.h>
#include <netinet/in.h>
#include <netdb.h>
#include <unistd.h>
"""

CONSTANTS_AND_GLOBAL_VARIABLES_STRING = """
static int ERROR_FLAG = 256;
static int HEADER_SIZE = 10; //integers

static int HEADER_FLAGS = 0;
static int HEADER_CALL_ID = 1;
static int HEADER_FUNCTION_ID = 2;
static int HEADER_CALL_COUNT = 3;
static int HEADER_INTEGER_COUNT = 4;
static int HEADER_LONG_COUNT = 5;
static int HEADER_FLOAT_COUNT = 6;
static int HEADER_DOUBLE_COUNT = 7;
static int HEADER_BOOLEAN_COUNT = 8;
static int HEADER_STRING_COUNT = 9;

static bool TRUE_BYTE = 1;
static bool FALSE_BYTE = 0;

static bool mpiIntercom = false;

static int socketfd = 0;

static int * header_in;
static int * header_out;

static int * ints_in;
static int * ints_out;

static long long int * longs_in;
static long long int * longs_out;

static float * floats_in;
static float * floats_out;

static double * doubles_in;
static double * doubles_out;

static int * booleans_in;
static int * booleans_out;

/* sizes of strings */
static int * string_sizes_in;
static int * string_sizes_out;

/* pointers to input and output strings (contents not stored here) */
static char * * strings_in;
static char * * strings_out;

/* actual string data */
static char * characters_in = 0;
static char * characters_out = 0;
"""


POLLING_FUNCTIONS_STRING = """
static int polling_interval = 0;

int internal__get_message_polling_interval(int * outval)
{
    *outval = polling_interval;
    
    return 0;
}

int internal__set_message_polling_interval(int inval)
{
    polling_interval = inval;
    
    return 0;
}
"""


RECV_HEADER_SLEEP_STRING = """
#include <unistd.h>

int mpi_recv_header(MPI_Comm & parent)
{
    MPI_Request header_request;
    MPI_Status request_status;
   
    MPI_Irecv(header_in, HEADER_SIZE, MPI_INT, 0, 989, parent, &header_request);
        
    if(polling_interval > 0)
    {
        int is_finished = 0;
        MPI_Test(&header_request, &is_finished, &request_status);
        while(!is_finished) {
            usleep(polling_interval);
            MPI_Test(&header_request, &is_finished, &request_status);
        }
        MPI_Wait(&header_request, &request_status);
    } else {
        MPI_Wait(&header_request, &request_status);
    }
    return 0;
}
"""

FOOTER_CODE_STRING = """
void onexit_mpi(void) {
#ifndef NOMPI
    int flag = 0;
    MPI_Finalized(&flag);
    
    if(!flag) {
        MPI_Comm parent;
        MPI_Comm_get_parent(&parent);
        
        int rank = 0;
        
        MPI_Comm_rank(parent, &rank);
        
        header_out[HEADER_FLAGS] = ERROR_FLAG;

        header_out[HEADER_CALL_ID] = 0;
        header_out[HEADER_FUNCTION_ID] = 0;
        header_out[HEADER_CALL_COUNT] = 0;
        header_out[HEADER_INTEGER_COUNT] = 0;
        header_out[HEADER_LONG_COUNT] = 0;
        header_out[HEADER_FLOAT_COUNT] = 0;
        header_out[HEADER_DOUBLE_COUNT] = 0;
        header_out[HEADER_BOOLEAN_COUNT] = 0;
        header_out[HEADER_STRING_COUNT] = 0;

        MPI_Send(header_out, HEADER_SIZE, MPI_INT, 0, 999, parent);
        
        MPI_Comm_disconnect(&parent);
        
        MPI_Finalize();
    }
#endif
}

void onexit_sockets(void) {
    close(socketfd);
}

void send_array_sockets(void *buffer, int length, int file_descriptor, int rank) {
    int total_written = 0;
    int bytes_written;

    if (rank != 0) {
        return;
    }

    while (total_written < length) {
        bytes_written = write(file_descriptor, ((char *) buffer) + total_written,
                        length - total_written);

        if (bytes_written == -1) {
            fprintf(stderr, "could not write data\\n");
            exit(1);
        }

        total_written = total_written + bytes_written;
    }
}

void receive_array_sockets(void *buffer, int length, int file_descriptor, int rank) {
    int total_read = 0;
    int bytes_read;

    if (rank != 0) {
        return;
    }

    while (total_read < length) {
        bytes_read = read(file_descriptor, ((char *) buffer) + total_read,
                        length - total_read);

        if (bytes_read == -1) {
            fprintf(stderr, "could not read data\\n");
            exit(1);
        }

        total_read = total_read + bytes_read;
    }
}

void new_arrays(int max_call_count) {
  ints_in = new int[ max_call_count * MAX_INTS_IN];
  ints_out = new int[ max_call_count * MAX_INTS_OUT];

  longs_in = new long long int[ max_call_count * MAX_LONGS_IN];
  longs_out = new long long int[ max_call_count * MAX_LONGS_OUT];

  floats_in = new float[ max_call_count * MAX_FLOATS_IN];
  floats_out = new float[ max_call_count * MAX_FLOATS_OUT];

  doubles_in = new double[ max_call_count * MAX_DOUBLES_IN];
  doubles_out = new double[ max_call_count * MAX_DOUBLES_OUT];
  
  booleans_in = new int[ max_call_count * MAX_BOOLEANS_IN];
  booleans_out = new int[ max_call_count * MAX_BOOLEANS_OUT];
  
  string_sizes_in = new int[ max_call_count * MAX_STRINGS_IN];
  string_sizes_out = new int[ max_call_count * MAX_STRINGS_OUT];

  strings_in = new char *[ max_call_count * MAX_STRINGS_IN];
  strings_out = new char *[ max_call_count * MAX_STRINGS_OUT];
}

void delete_arrays() {
  delete[] ints_in;
  delete[] ints_out;
  delete[] longs_in;
  delete[] longs_out;
  delete[] floats_in;
  delete[] floats_out;
  delete[] doubles_in;
  delete[] doubles_out;
  delete[] booleans_in;
  delete[] booleans_out;
  delete[] string_sizes_in;
  delete[] string_sizes_out;
  delete[] strings_in;
  delete[] strings_out;
}

void run_mpi(int argc, char *argv[]) {
#ifndef NOMPI
  int provided;
  MPI_Comm parent;
  int rank = 0;
  
  mpiIntercom = true;

  //fprintf(stderr, "C worker: running in mpi mode\\n");
  
  MPI_Init_thread(&argc, &argv, MPI_THREAD_MULTIPLE, &provided);
  MPI_Comm_get_parent(&parent);
  MPI_Comm_rank(parent, &rank);
  atexit(onexit_mpi);
  
  bool must_run_loop = true;
  
  int max_call_count = 10;
  
  header_in = new int[HEADER_SIZE];
  header_out = new int[HEADER_SIZE];

  new_arrays(max_call_count);  

  while(must_run_loop) {
    //fprintf(stderr, "receiving header\\n");
    mpi_recv_header(parent);
    
    //fprintf(stderr, "C worker code: got header %d %d %d %d %d %d %d %d %d %d\\n", header_in[0], header_in[1], header_in[2], header_in[3], header_in[4], header_in[5], header_in[6], header_in[7], header_in[8], header_in[9]);
    
    int call_count = header_in[HEADER_CALL_COUNT];

    if (call_count > max_call_count) {
      delete_arrays();
      max_call_count = call_count + 255;
      new_arrays(max_call_count);
    }
    
    if(header_in[HEADER_INTEGER_COUNT] > 0) {
      MPI_Bcast(ints_in, header_in[HEADER_INTEGER_COUNT] , MPI_INT, 0, parent);
    }
    
    if(header_in[HEADER_LONG_COUNT] > 0) {
      MPI_Bcast(longs_in, header_in[HEADER_LONG_COUNT], MPI_LONG_LONG_INT, 0, parent);
    }
    
    if(header_in[HEADER_FLOAT_COUNT] > 0) {
      MPI_Bcast(floats_in, header_in[HEADER_FLOAT_COUNT], MPI_FLOAT, 0, parent);
    }
    
    if(header_in[HEADER_DOUBLE_COUNT] > 0) {
      MPI_Bcast(doubles_in, header_in[HEADER_DOUBLE_COUNT], MPI_DOUBLE, 0, parent);
    }
    
    if(header_in[HEADER_BOOLEAN_COUNT] > 0) {
      MPI_Bcast(booleans_in, header_in[HEADER_BOOLEAN_COUNT], MPI_INTEGER, 0, parent);
    }
    
    if(header_in[HEADER_STRING_COUNT] > 0) {
      MPI_Bcast(string_sizes_in, header_in[HEADER_STRING_COUNT], MPI_INTEGER, 0, parent);
      
      int total_string_size = 0;
      for (int i = 0; i < header_in[HEADER_STRING_COUNT];i++) {
        total_string_size += string_sizes_in[i] + 1;
      }
      
      characters_in = new char[total_string_size];
      MPI_Bcast(characters_in, total_string_size, MPI_CHARACTER, 0, parent);

      int offset = 0;
      for (int i = 0 ; i <  header_in[HEADER_STRING_COUNT];i++) {
          strings_in[i] = characters_in + offset;
          offset += string_sizes_in[i] + 1;
      } 
    }

    header_out[HEADER_FLAGS] = 0;
    header_out[HEADER_CALL_ID] = header_in[HEADER_CALL_ID];
    header_out[HEADER_FUNCTION_ID] = header_in[HEADER_FUNCTION_ID];
    header_out[HEADER_CALL_COUNT] = call_count;
    header_out[HEADER_INTEGER_COUNT] = 0;
    header_out[HEADER_LONG_COUNT] = 0;
    header_out[HEADER_FLOAT_COUNT] = 0;
    header_out[HEADER_DOUBLE_COUNT] = 0;
    header_out[HEADER_BOOLEAN_COUNT] = 0;
    header_out[HEADER_STRING_COUNT] = 0;

    //fprintf(stderr, "c worker mpi: handling call\\n");
    
    must_run_loop = handle_call();
    
    //fprintf(stderr, "c worker mpi: call handled\\n");
    
    MPI_Barrier(MPI_COMM_WORLD);
    
    if(rank == 0) {
      MPI_Send(header_out, HEADER_SIZE, MPI_INT, 0, 999, parent);
      
      if(header_out[HEADER_INTEGER_COUNT] > 0) {
        MPI_Send(ints_out, header_out[HEADER_INTEGER_COUNT], MPI_INT, 0, 999, parent);
      }
      if(header_out[HEADER_LONG_COUNT] > 0) {
        MPI_Send(longs_out, header_out[HEADER_LONG_COUNT], MPI_LONG_LONG_INT, 0, 999, parent);
      }
      if(header_out[HEADER_FLOAT_COUNT] > 0) {
        MPI_Send(floats_out, header_out[HEADER_FLOAT_COUNT], MPI_FLOAT, 0, 999, parent);
      }
      if(header_out[HEADER_DOUBLE_COUNT] > 0) {
        MPI_Send(doubles_out, header_out[HEADER_DOUBLE_COUNT], MPI_DOUBLE, 0, 999, parent);
      }
      if(header_out[HEADER_BOOLEAN_COUNT] > 0) {
        MPI_Send(booleans_out, header_out[HEADER_BOOLEAN_COUNT], MPI_INTEGER, 0, 999, parent);
      }
      if(header_out[HEADER_STRING_COUNT] > 0) {
        int offset = 0;
        for( int i = 0; i < header_out[HEADER_STRING_COUNT] ; i++) {
          
          int length = strlen(strings_out[i]);
          string_sizes_out[i] = length;
          offset += length + 1;
        }
        
        characters_out = new char[offset + 1];
        offset = 0;
        
        for( int i = 0; i < header_out[HEADER_STRING_COUNT]  ; i++) {
          strcpy(characters_out+offset, strings_out[i]);
          offset += string_sizes_out[i] + 1;
        }
        MPI_Send(string_sizes_out, header_out[HEADER_STRING_COUNT], MPI_INTEGER, 0, 999, parent);
        MPI_Send(characters_out, offset, MPI_BYTE, 0, 999, parent);
      }
    
    }
    
    if (characters_in) { 
        delete[] characters_in;
        characters_in = 0;
    }
    
    if (characters_out) {
        delete[] characters_out;
        characters_out = 0;
    }
    //fprintf(stderr, "call done\\n");
  }
  delete_arrays();
  
  MPI_Comm_disconnect(&parent);
  MPI_Finalize();
  //fprintf(stderr, "mpi finalized\\n");
#else
  fprintf(stderr, "mpi support not compiled into worker\\n");
  exit(1);
#endif
}

void run_sockets_mpi(int argc, char *argv[], int port) {
#ifndef NOMPI
 bool must_run_loop = true;
  int max_call_count = 10;
  struct sockaddr_in serv_addr;
  struct hostent *server;
  
  mpiIntercom = false;

  MPI::Init_thread(argc, argv, MPI_THREAD_MULTIPLE);
  int rank = MPI::COMM_WORLD.Get_rank();
  
  if (rank == 0) {  
    //fprintf(stderr, "C worker: running in sockets+mpi mode\\n");
  
   
    socketfd = socket(AF_INET, SOCK_STREAM, 0);
    
    if (socketfd < 0) {
      fprintf(stderr, "cannot open socket\\n");
      exit(1);
    }
    
    server = gethostbyname("localhost");
    
    bzero((char *) &serv_addr, sizeof(serv_addr));
    serv_addr.sin_family = AF_INET;
    bcopy((char *) server->h_addr, (char *) &serv_addr.sin_addr.s_addr, server->h_length);
    serv_addr.sin_port = htons(port);
  
    if (connect(socketfd, (struct sockaddr *) &serv_addr, sizeof(serv_addr)) < 0) {
      fprintf(stderr, "cannot connect socket\\n");
      exit(1);
    }
    
    //fprintf(stderr, "sockets_mpi: finished initializing code\\n");
  
    atexit(onexit_sockets);
  
  }
  
  header_in = new int[HEADER_SIZE];
  header_out = new int[HEADER_SIZE];

  new_arrays(max_call_count);  
  
  while(must_run_loop) {
    //fprintf(stderr, "sockets_mpi: receiving header\\n");
    receive_array_sockets(header_in, HEADER_SIZE * sizeof(int), socketfd, rank);
    MPI::COMM_WORLD.Bcast(header_in, HEADER_SIZE, MPI_INT, 0);
    
    //fprintf(stderr, "C sockets_mpi worker code: got header %d %d %d %d %d %d %d %d %d %d\\n", header_in[0], header_in[1], header_in[2], header_in[3], header_in[4], header_in[5], header_in[6], header_in[7], header_in[8], header_in[9]);
    
    int call_count = header_in[HEADER_CALL_COUNT];

    if (call_count > max_call_count) {
      delete_arrays();
      max_call_count = call_count + 255;
      new_arrays(max_call_count);
    }
    
    if (header_in[HEADER_INTEGER_COUNT] > 0) {
      receive_array_sockets(ints_in, header_in[HEADER_INTEGER_COUNT] * sizeof(int), socketfd, rank);
      MPI::COMM_WORLD.Bcast(ints_in, header_in[HEADER_INTEGER_COUNT], MPI_INTEGER, 0);
    }
     
    if (header_in[HEADER_LONG_COUNT] > 0) {
      receive_array_sockets(longs_in, header_in[HEADER_LONG_COUNT] * sizeof(long long int), socketfd, rank);
      MPI::COMM_WORLD.Bcast(longs_in, header_in[HEADER_LONG_COUNT], MPI_LONG_LONG_INT, 0);
    }
    
    if(header_in[HEADER_FLOAT_COUNT] > 0) {
      receive_array_sockets(floats_in, header_in[HEADER_FLOAT_COUNT] * sizeof(float), socketfd, rank);
      MPI::COMM_WORLD.Bcast(floats_in, header_in[HEADER_FLOAT_COUNT], MPI_FLOAT, 0);
    }
    
    if(header_in[HEADER_DOUBLE_COUNT] > 0) {
      receive_array_sockets(doubles_in, header_in[HEADER_DOUBLE_COUNT] * sizeof(double), socketfd, rank);
      MPI::COMM_WORLD.Bcast(doubles_in, header_in[HEADER_DOUBLE_COUNT], MPI_DOUBLE, 0);
    }
    
    if(header_in[HEADER_BOOLEAN_COUNT] > 0) {
      //receive_array_sockets(booleans_in, header_in[HEADER_BOOLEAN_COUNT], socketfd , rank);
      for (int i = 0; i < header_in[HEADER_BOOLEAN_COUNT]; i++) {
        booleans_in[i] = 0;
        receive_array_sockets(&booleans_in[i], 1, socketfd , rank);
      }
      MPI::COMM_WORLD.Bcast(booleans_in, header_in[HEADER_BOOLEAN_COUNT], MPI_INTEGER, 0);
    }
    
    if(header_in[HEADER_STRING_COUNT] > 0) {
      receive_array_sockets(string_sizes_in, header_in[HEADER_STRING_COUNT] * sizeof(int), socketfd, rank);
      MPI::COMM_WORLD.Bcast(string_sizes_in, header_in[HEADER_STRING_COUNT], MPI_INT, 0);
      for (int i = 0; i < header_in[HEADER_STRING_COUNT]; i++) {
        strings_in[i] = new char[string_sizes_in[i] + 1];
        receive_array_sockets(strings_in[i], string_sizes_in[i], socketfd, rank);
        MPI::COMM_WORLD.Bcast(strings_in[i], string_sizes_in[i], MPI_CHARACTER, 0);
        strings_in[i][string_sizes_in[i]] = '\\0';
      }
    }
    
    header_out[HEADER_FLAGS] = 0;
    header_out[HEADER_CALL_ID] = header_in[HEADER_CALL_ID];
    header_out[HEADER_FUNCTION_ID] = header_in[HEADER_FUNCTION_ID];
    header_out[HEADER_CALL_COUNT] = call_count;
    header_out[HEADER_INTEGER_COUNT] = 0;
    header_out[HEADER_LONG_COUNT] = 0;
    header_out[HEADER_FLOAT_COUNT] = 0;
    header_out[HEADER_DOUBLE_COUNT] = 0;
    header_out[HEADER_BOOLEAN_COUNT] = 0;
    header_out[HEADER_STRING_COUNT] = 0;

    //fprintf(stderr, "c worker sockets_mpi: handling call\\n");
    
    must_run_loop = handle_call();
    
    //fprintf(stderr, "c worker sockets_mpi: call handled\\n");
    
    MPI_Barrier(MPI_COMM_WORLD);
    
    if (rank == 0) {

      send_array_sockets(header_out, HEADER_SIZE * sizeof(int), socketfd, 0);
          
      if(header_out[HEADER_INTEGER_COUNT] > 0) {
        send_array_sockets(ints_out, header_out[HEADER_INTEGER_COUNT] * sizeof(int), socketfd, 0);
      }
          
      if(header_out[HEADER_LONG_COUNT] > 0) {
        send_array_sockets(longs_out, header_out[HEADER_LONG_COUNT] * sizeof(long long int), socketfd, 0);
      }
          
      if(header_out[HEADER_FLOAT_COUNT] > 0) {
        send_array_sockets(floats_out, header_out[HEADER_FLOAT_COUNT] * sizeof(float), socketfd, 0);
      }
          
      if(header_out[HEADER_DOUBLE_COUNT] > 0) {
        send_array_sockets(doubles_out, header_out[HEADER_DOUBLE_COUNT] * sizeof(double), socketfd, 0);
      }
          
      if(header_out[HEADER_BOOLEAN_COUNT] > 0) {
          for (int i = 0; i < header_out[HEADER_BOOLEAN_COUNT]; i++) {
            if (booleans_out[i]) {
              send_array_sockets(&TRUE_BYTE, 1, socketfd, 0);
            } else {
              send_array_sockets(&FALSE_BYTE, 1, socketfd, 0);
            }
         }
      }
          
      if(header_out[HEADER_STRING_COUNT] > 0) {
        for (int i = 0; i < header_out[HEADER_STRING_COUNT]; i++) {
          string_sizes_out[i] = strlen(strings_out[i]);
        }
        send_array_sockets(string_sizes_out, header_out[HEADER_STRING_COUNT] * sizeof(int), socketfd, 0);
          
        for (int i = 0; i < header_out[HEADER_STRING_COUNT]; i++) {
          send_array_sockets(strings_out[i], string_sizes_out[i] * sizeof(char), socketfd, 0);
        }
      }
        
      //fprintf(stderr, "sockets_mpicall done\\n");
    }

  }
  delete_arrays();
  
  if (rank == 0) {
    close(socketfd);
  }
  
  MPI_Finalize();
  
  //fprintf(stderr, "sockets_mpi done\\n");
#else
  fprintf(stderr, "mpi support not compiled into worker\\n");
  exit(1);
#endif
}

void run_sockets(int port) {
  bool must_run_loop = true;
  int max_call_count = 10;
  struct sockaddr_in serv_addr;
  struct hostent *server;
  
  mpiIntercom = false;

  //fprintf(stderr, "C worker: running in sockets mode\\n");
   
  socketfd = socket(AF_INET, SOCK_STREAM, 0);
    
  if (socketfd < 0) {
    fprintf(stderr, "cannot open socket\\n");
    exit(1);
  }
    
  server = gethostbyname("localhost");
    
  bzero((char *) &serv_addr, sizeof(serv_addr));
  serv_addr.sin_family = AF_INET;
  bcopy((char *) server->h_addr, (char *) &serv_addr.sin_addr.s_addr, server->h_length);
  serv_addr.sin_port = htons(port);
  
  if (connect(socketfd, (struct sockaddr *) &serv_addr, sizeof(serv_addr)) < 0) {
    fprintf(stderr, "cannot connect socket\\n");
    exit(1);
  }
    
  //fprintf(stderr, "sockets: finished initializing code\\n");
  
  atexit(onexit_sockets);
  
  header_in = new int[HEADER_SIZE];
  header_out = new int[HEADER_SIZE];

  new_arrays(max_call_count);  
  
  while(must_run_loop) {
    //fprintf(stderr, "sockets: receiving header\\n");
    receive_array_sockets(header_in, HEADER_SIZE * sizeof(int), socketfd, 0);
    //fprintf(stderr, "C sockets worker code: got header %d %d %d %d %d %d %d %d %d %d\\n", header_in[0], header_in[1], header_in[2], header_in[3], header_in[4], header_in[5], header_in[6], header_in[7], header_in[8], header_in[9]);
    
    int call_count = header_in[HEADER_CALL_COUNT];

    if (call_count > max_call_count) {
      delete_arrays();
      max_call_count = call_count + 255;
      new_arrays(max_call_count);
    }
    
    if (header_in[HEADER_INTEGER_COUNT] > 0) {
      receive_array_sockets(ints_in, header_in[HEADER_INTEGER_COUNT] * sizeof(int), socketfd, 0);
    }
     
    if (header_in[HEADER_LONG_COUNT] > 0) {
      receive_array_sockets(longs_in, header_in[HEADER_LONG_COUNT] * sizeof(long long int), socketfd, 0);
    }
    
    if(header_in[HEADER_FLOAT_COUNT] > 0) {
      receive_array_sockets(floats_in, header_in[HEADER_FLOAT_COUNT] * sizeof(float), socketfd, 0);
    }
    
    if(header_in[HEADER_DOUBLE_COUNT] > 0) {
      receive_array_sockets(doubles_in, header_in[HEADER_DOUBLE_COUNT] * sizeof(double), socketfd, 0);
    }
    
    if(header_in[HEADER_BOOLEAN_COUNT] > 0) {
      //receive_array_sockets(booleans_in, header_in[HEADER_BOOLEAN_COUNT], socketfd , 0);
      for (int i = 0; i < header_in[HEADER_BOOLEAN_COUNT]; i++) {
        booleans_in[i] = 0;
        receive_array_sockets(&booleans_in[i], 1, socketfd , 0);
      }
    }
    
    if(header_in[HEADER_STRING_COUNT] > 0) {
      receive_array_sockets(string_sizes_in, header_in[HEADER_STRING_COUNT] * sizeof(int), socketfd, 0);
      for (int i = 0; i < header_in[HEADER_STRING_COUNT]; i++) {
        strings_in[i] = new char[string_sizes_in[i] + 1];
        receive_array_sockets(strings_in[i], string_sizes_in[i], socketfd, 0);
        strings_in[i][string_sizes_in[i]] = '\\0';
      }
    }
    
    header_out[HEADER_FLAGS] = 0;
    header_out[HEADER_CALL_ID] = header_in[HEADER_CALL_ID];
    header_out[HEADER_FUNCTION_ID] = header_in[HEADER_FUNCTION_ID];
    header_out[HEADER_CALL_COUNT] = call_count;
    header_out[HEADER_INTEGER_COUNT] = 0;
    header_out[HEADER_LONG_COUNT] = 0;
    header_out[HEADER_FLOAT_COUNT] = 0;
    header_out[HEADER_DOUBLE_COUNT] = 0;
    header_out[HEADER_BOOLEAN_COUNT] = 0;
    header_out[HEADER_STRING_COUNT] = 0;

    //fprintf(stderr, "c worker sockets: handling call\\n");
    
    must_run_loop = handle_call();
    
    //fprintf(stderr, "c worker sockets: call handled\\n");

    send_array_sockets(header_out, HEADER_SIZE * sizeof(int), socketfd, 0);
      
    if(header_out[HEADER_INTEGER_COUNT] > 0) {
      send_array_sockets(ints_out, header_out[HEADER_INTEGER_COUNT] * sizeof(int), socketfd, 0);
    }
      
    if(header_out[HEADER_LONG_COUNT] > 0) {
      send_array_sockets(longs_out, header_out[HEADER_LONG_COUNT] * sizeof(long long int), socketfd, 0);
    }
      
    if(header_out[HEADER_FLOAT_COUNT] > 0) {
      send_array_sockets(floats_out, header_out[HEADER_FLOAT_COUNT] * sizeof(float), socketfd, 0);
    }
      
    if(header_out[HEADER_DOUBLE_COUNT] > 0) {
      send_array_sockets(doubles_out, header_out[HEADER_DOUBLE_COUNT] * sizeof(double), socketfd, 0);
    }
      
    if(header_out[HEADER_BOOLEAN_COUNT] > 0) {
        for (int i = 0; i < header_out[HEADER_BOOLEAN_COUNT]; i++) {
          if (booleans_out[i]) {
            send_array_sockets(&TRUE_BYTE, 1, socketfd, 0);
          } else {
            send_array_sockets(&FALSE_BYTE, 1, socketfd, 0);
          }
       }
    }
      
    if(header_out[HEADER_STRING_COUNT] > 0) {
      for (int i = 0; i < header_out[HEADER_STRING_COUNT]; i++) {
        string_sizes_out[i] = strlen(strings_out[i]);
      }
      send_array_sockets(string_sizes_out, header_out[HEADER_STRING_COUNT] * sizeof(int), socketfd, 0);
      
      for (int i = 0; i < header_out[HEADER_STRING_COUNT]; i++) {
        send_array_sockets(strings_out[i], string_sizes_out[i] * sizeof(char), socketfd, 0);
      }
    }
    
    //fprintf(stderr, "call done\\n");
  }
  delete_arrays();
  
  close(socketfd);
  //fprintf(stderr, "sockets done\\n");
}
 
int main(int argc, char *argv[]) {
  int port;
  bool use_mpi = NEEDS_MPI;
  
  for(int i = 0 ; i < argc; i++) {
    //fprintf(stderr, "argument %d is %s\\n", i, argv[i]);
  }

  if (argc == 1) {
    run_mpi(argc, argv);
  } else {
    port = atoi(argv[1]);
    
    if (argc >= 3) {
       if (strcmp(argv[2], "no_mpi") == 0) {
         use_mpi = false;
       } else if (strcmp(argv[2], "mpi") == 0) {
         use_mpi = true;
       } else if (strcmp(argv[2], "auto") == 0) {
         use_mpi = NEEDS_MPI;
       }

    }
    
    if (use_mpi) {
      run_sockets_mpi(argc, argv, port);
    } else {
      run_sockets(port);
    }
  }  

  return 0;
}   


"""

class MakeCCodeString(GenerateASourcecodeString):
    @late
    def dtype_to_spec(self):
        return dtype_to_spec
       
         

class GenerateACStringOfAFunctionSpecification(MakeCCodeString):
    @late
    def specification(self):
        raise exceptions.AmuseException("No specification set, please set the specification first")
   
        
    def start(self):
        
        self.specification.prepare_output_parameters()
        self.output_casestmt_start()
        self.out.indent()
        
        if self.specification.must_handle_array:
            pass
        elif self.specification.can_handle_array:
            self.out.lf() + 'for (int i = 0 ; i < call_count; i++){'
            self.out.indent()
 
        self.output_copy_inout_variables()
        self.output_function_start()
        self.output_function_parameters()
        self.output_function_end()
        
        if self.specification.must_handle_array:
            if not self.specification.result_type is None:
                spec = self.dtype_to_spec[self.specification.result_type]
                self.out.lf() + 'for (int i = 1 ; i < call_count; i++){'
                self.out.indent()
                self.out.lf() + spec.output_var_name + '[i]' + ' = ' + spec.output_var_name + '[0]' + ';'
                self.out.dedent()
                self.out.lf() + '}'
        elif self.specification.can_handle_array:
            self.out.dedent()
            self.out.lf() + '}'
        
        self.output_lines_with_number_of_outputs()
        self.output_casestmt_end()
        self.out.dedent()
        self._result = self.out.string
    
    def index_string(self, index, must_copy_in_to_out=False):
        if self.specification.must_handle_array and not must_copy_in_to_out:
            if index == 0:
                return '0'
            else:
                return '( %d * call_count)' % index
        elif self.specification.can_handle_array or (self.specification.must_handle_array and must_copy_in_to_out):
            if index == 0:
                return 'i'
            else:
                return '( %d * call_count) + i' % index
        else:
            return index
    
    
    def input_var(self, name, index):
        if self.specification.must_handle_array:
            self.output_var(name, index)
        else:
            self.out.n() + name
            self.out + '[' + self.index_string(index) + ']'
        
    def output_var(self, name, index):
        self.out.n() + '&' + name
        self.out + '[' + self.index_string(index) + ']'
    
    def output_function_parameters(self):
        self.out.indent()
        
        first = True
        
        for parameter in self.specification.parameters:
            spec = self.dtype_to_spec[parameter.datatype]
            
            if first:
                first = False
            else:
                self.out + ' ,'
                
            if parameter.direction == LegacyFunctionSpecification.IN:
                    self.input_var(spec.input_var_name, parameter.input_index)
            if parameter.direction == LegacyFunctionSpecification.INOUT:
                    self.output_var(spec.output_var_name, parameter.output_index)
            elif parameter.direction == LegacyFunctionSpecification.OUT:
                    self.output_var(spec.output_var_name, parameter.output_index)
            elif parameter.direction == LegacyFunctionSpecification.LENGTH:
                self.out.n() + 'call_count'
    
        self.out.dedent()
        
        
    def output_copy_inout_variables(self):
        for parameter in self.specification.parameters:
            spec = self.dtype_to_spec[parameter.datatype]
            
            if parameter.direction == LegacyFunctionSpecification.INOUT:
                if self.specification.must_handle_array:
                    self.out.lf() + 'for (int i = 0 ; i < call_count; i++){'
                    self.out.indent()

                self.out.n() + spec.output_var_name
                self.out + '[' + self.index_string(parameter.output_index, must_copy_in_to_out=True) + ']'
                self.out + ' = '
                self.out + spec.input_var_name + '[' + self.index_string(parameter.input_index, must_copy_in_to_out=True) + ']' + ';'
            
                if self.specification.must_handle_array:
                    self.out.dedent()
                    self.out.lf() + '}'
                                
    def output_lines_with_number_of_outputs(self):
        dtype_to_count = {}
        
        for parameter in self.specification.output_parameters:
            count = dtype_to_count.get(parameter.datatype, 0)
            dtype_to_count[parameter.datatype] = count + 1
                
        if not self.specification.result_type is None:
            count = dtype_to_count.get(self.specification.result_type, 0)
            dtype_to_count[self.specification.result_type] = count + 1
            
        for dtype in dtype_to_count:       
            spec = self.dtype_to_spec[dtype]
            count = dtype_to_count[dtype]
            self.out.n() 
            self.out + 'header_out[' + spec.counter_name 
            self.out + '] = ' + count + ' * call_count;'
            pass
            
    def output_function_end(self):
        if len(self.specification.parameters) > 0:
            self.out.n()
            
        self.out + ')' + ';'
        
    def output_function_start(self):
        self.out.n() 
        if not self.specification.result_type is None:
            spec = self.dtype_to_spec[self.specification.result_type]
            self.out + spec.output_var_name
            self.out + '[' + self.index_string(0) + ']' + ' = '
        self.out + self.specification.name + '('
        
    def output_casestmt_start(self):
        self.out + 'case ' + self.specification.id + ':'
        
    def output_casestmt_end(self):
        self.out.n() + 'break;'
        
        

class GenerateACHeaderDefinitionStringFromAFunctionSpecification(MakeCCodeString):
   
        
    def start(self):
        self.output_function_start()
        self.output_function_parameters()
        self.output_function_end()
        self._result = self.out.string
            
    def output_function_parameters(self):        
        first = True
        
        for parameter in self.specification.parameters:
            spec = self.dtype_to_spec[parameter.datatype]
            
            if first:
                first = False
            else:
                self.out + ', '
                
            if parameter.datatype == 'string':
                self.out + 'char'
            else:
                self.out + spec.type
            self.out + ' '
            if parameter.is_output() or (parameter.is_input() and self.specification.must_handle_array):
                self.out + '*' + ' '
            if parameter.datatype == 'string':
                self.out + '*' + ' '
            self.out + parameter.name
                
            
    def output_function_end(self):
        self.out + ')' + ';'
        
    def output_function_start(self):
        self.out.n()
        if not self.specification.result_type is None:
            spec = self.dtype_to_spec[self.specification.result_type]
            self.out + spec.type
            self.out + ' '
        else:
            self.out + 'void' + ' '
        self.out + self.specification.name + '('
        
class GenerateACSourcecodeStringFromASpecificationClass\
    (GenerateASourcecodeStringFromASpecificationClass):

    @late
    def specification_class(self):
        raise exceptions.AmuseException("No specification_class set, please set the specification_class first")
    
    @late
    def dtype_to_spec(self):
        return dtype_to_spec

    def output_sourcecode_for_function(self):
        return GenerateACStringOfAFunctionSpecification()
    
    def start(self):
        self.out + HEADER_CODE_STRING

        self.output_local_includes()
        
        self.output_needs_mpi()
        
        self.output_code_constants()
        
        self.out.lf() + CONSTANTS_AND_GLOBAL_VARIABLES_STRING
        
        self.out.lf() + POLLING_FUNCTIONS_STRING
        
        if self.must_generate_mpi:
            self.out.lf() + RECV_HEADER_SLEEP_STRING
        
        self.output_handle_call()
        
        self.out.lf() + FOOTER_CODE_STRING
        
        self._result = self.out.string
        
    def output_local_includes(self):
        if hasattr(self.specification_class, 'include_headers'):
            for x in self.specification_class.include_headers:
                self.out.n() + '#include "' + x + '"'
        self.out.lf()

        
    def output_needs_mpi(self):
        if self.needs_mpi:
            self.out.lf() + 'static bool NEEDS_MPI = true;'
        else:
            self.out.lf() + 'static bool NEEDS_MPI = false;'
        self.out.lf().lf()
    
    def output_code_constants(self):
        for dtype in self.dtype_to_spec.keys():
            dtype_spec = self.dtype_to_spec[dtype]
            
            maxin = self.mapping_from_dtype_to_maximum_number_of_inputvariables.get(dtype, 0)
            self.out + 'static int MAX_' + dtype_spec.input_var_name.upper() + ' = ' + maxin + ";"
            self.out.lf()
            
            maxout = self.mapping_from_dtype_to_maximum_number_of_outputvariables.get(dtype, 0)
            self.out + 'static int MAX_' + dtype_spec.output_var_name.upper() + ' = ' + maxout + ";"
            self.out.lf()
            
    def output_handle_call(self):
        self.out.lf().lf() + 'bool handle_call() {'
        self.out.indent()
        
        self.out.lf() + 'int call_count = header_in[HEADER_CALL_COUNT];'
        
        self.out.lf().lf() + 'switch(header_in[HEADER_FUNCTION_ID]) {'
        self.out.indent()
        self.out.lf() + 'case 0:'
        self.out.indent().lf() + 'return false;'
        self.out.lf() + 'break;'
        self.out.dedent()
        
        self.output_sourcecode_for_functions()
        
        self.out.lf() + 'default:'
        self.out.indent()
        self.out.lf() + 'header_out[HEADER_FLAGS] = header_out[HEADER_FLAGS] | ERROR_FLAG;'
        self.out.lf() + 'strings_out[0] = new char[100];'
        self.out.lf() + 'sprintf(strings_out[0], "unknown function id: %d\\n", header_in[HEADER_FUNCTION_ID]);'
        self.out.lf() + 'fprintf(stderr, "unknown function id: %d\\n", header_in[HEADER_FUNCTION_ID]);'
        self.out.lf() + 'header_out[HEADER_STRING_COUNT] = 1;'
        self.out.lf() + 'return false;'
        self.out.dedent()
        
        self.out.dedent().lf() + '}'
        self.out.dedent()
        self.out.indent().lf() + 'return true;'
        self.out.dedent().lf() + '}'

class GenerateACHeaderStringFromASpecificationClass\
    (GenerateASourcecodeStringFromASpecificationClass):

    @late
    def ignore_functions_from_specification_classes(self):
        return []
        
    @late
    def underscore_functions_from_specification_classes(self):
        return []
        
    @late
    def dtype_to_spec(self):
        return dtype_to_spec
        
    @late
    def make_extern_c(self):
        return True
    
    def must_include_interface_function_in_output(self, x):
        if x.specification.name.startswith("internal__"):
            return False
            
        for cls in self.ignore_functions_from_specification_classes:
            if hasattr(cls, x.specification.name):
                return False
        
        return True
        
    def output_sourcecode_for_function(self):
        return GenerateACHeaderDefinitionStringFromAFunctionSpecification()
        
    def start(self):  
        if self.make_extern_c:
            self.out + 'extern "C" {'
            self.out.indent().lf()
            
        self.output_sourcecode_for_functions()
        
        if self.make_extern_c:
            self.out.dedent().lf() + '}'
        
        self.out.lf()
        
        self._result = self.out.string
        

class GenerateACStubStringFromASpecificationClass\
    (GenerateASourcecodeStringFromASpecificationClass):

    @late
    def dtype_to_spec(self):
        return dtype_to_spec
        
    @late
    def make_extern_c(self):
        return False
    
    def output_sourcecode_for_function(self):
        return create_definition.CreateCStub()

    def must_include_interface_function_in_output(self, x):
        return not x.specification.name.startswith("internal__")
     
    def start(self):  
    
        self.output_local_includes()
        
        self.out.lf()
        
        if self.make_extern_c:
            self.out + 'extern "C" {'
            self.out.indent().lf()
            
        self.output_sourcecode_for_functions()
        
        if self.make_extern_c:
            self.out.dedent().lf() + '}'
        
        self.out.lf()
        
        self._result = self.out.string
        
    
    def output_local_includes(self):
        self.out.n()
        if hasattr(self.specification_class, 'include_headers'):
            for x in self.specification_class.include_headers:
                self.out.n() + '#include "' + x + '"'
    
        
        
        
