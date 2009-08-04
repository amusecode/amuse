#include <mpi.h>
#include "muse_dynamics.h"
#include "parameters.h"
#include "local.h"

    int _add_particle(int id, double mass, double radius, double x, double y, double z, double vx, double vy, double vz) {
        dynamics_state state;
        state.id = id;
        state.mass = mass;
        state.radius = radius;
        state.x = x;
        state.y = y;
        state.z = z;
        state.vx = vx;
        state.vy = vy;
        state.vz = vz;
        return add_particle(state);
    }
    
    void _get_state(int id, int * id_out,  double * mass, double * radius, double * x, double * y, double * z, double * vx, double * vy, double * vz) {
        dynamics_state state = get_state(id);
        *id_out = state.id;
        *mass = state.mass;
        *radius = state.radius;
        *x = state.x;
        *y = state.y;
        *z = state.z;
        *vx = state.vx;
        *vy = state.vy;
        *vz = state.vz;
    }
    

class message_header {

public:
	int tag;
	int number_of_doubles;
	int number_of_ints;
	message_header(): tag(0), number_of_doubles(0), number_of_ints(0) {}

	void send(MPI::Intercomm & intercomm, int rank){
		int header[3];
		header[0] = tag;
		header[1] = number_of_doubles;
		header[2] = number_of_ints;		
		intercomm.Send(header, 3, MPI_INT, 0, 999);
	}

	void recv(MPI::Intercomm & intercom, int rank) {
		int header[6];

		intercom.Recv(header, 3, MPI_INT, 0, rank);
		tag = header[0];
		number_of_doubles = header[1];
		number_of_ints = header[2];
	}
};

void run_loop() {
  int rank = MPI::COMM_WORLD.Get_rank();
  
  MPI::Intercomm parent = MPI::COMM_WORLD.Get_parent();
  
  bool must_run_loop = true;
  
  while(must_run_loop) {
    int ints_in[255];
    int ints_out[255];
    double doubles_in[255];
    double doubles_out[255];
    
    message_header request_header;
    message_header reply_header;
    
    request_header.recv(parent,rank);
    if(request_header.number_of_doubles > 0) {
      parent.Recv(doubles_in, request_header.number_of_doubles, MPI_DOUBLE, 0, rank);
    }
    if(request_header.number_of_ints > 0) {
      parent.Recv(ints_in, request_header.number_of_ints, MPI_INT, 0, rank);
    }
    
    reply_header.tag = request_header.tag;
    
    switch(request_header.tag) {
      case 0:
        must_run_loop = false;
        break;
      case 1:
        ints_out[0] = setup_module();
        reply_header.number_of_ints = 1;
        break;
      case 2:
        ints_out[0] = cleanup_module();
        reply_header.number_of_ints = 1;
        break;
      case 3:
        ints_out[0] = initialize_particles(
          doubles_in[0]
        );
        reply_header.number_of_ints = 1;
        break;
      case 4:
        ints_out[0] = reinitialize_particles();
        reply_header.number_of_ints = 1;
        break;
      case 5:
        ints_out[0] = _add_particle(
          ints_in[0] ,
          doubles_in[0] ,
          doubles_in[1] ,
          doubles_in[2] ,
          doubles_in[3] ,
          doubles_in[4] ,
          doubles_in[5] ,
          doubles_in[6] ,
          doubles_in[7]
        );
        reply_header.number_of_ints = 1;
        break;
      case 6:
        ints_out[0] = evolve(
          doubles_in[0] ,
          ints_in[0]
        );
        reply_header.number_of_ints = 1;
        break;
      case 7:
        ints_out[0] = get_number();
        reply_header.number_of_ints = 1;
        break;
      case 8:
        _get_state(
          ints_in[0] ,
          &ints_out[0] ,
          &doubles_out[0] ,
          &doubles_out[1] ,
          &doubles_out[2] ,
          &doubles_out[3] ,
          &doubles_out[4] ,
          &doubles_out[5] ,
          &doubles_out[6] ,
          &doubles_out[7]
        );
        reply_header.number_of_ints = 1;
        reply_header.number_of_doubles = 8;
        break;
      case 20:
        if(request_header.number_of_doubles == 1){
          t = doubles_in[0];
        } else {
          reply_header.number_of_doubles = 1;
          doubles_out[0] = t;
        }
        break;
      case 21:
        if(request_header.number_of_doubles == 1){
          dt_param = doubles_in[0];
        } else {
          reply_header.number_of_doubles = 1;
          doubles_out[0] = dt_param;
        }
        break;
      case 22:
        if(request_header.number_of_doubles == 1){
          dt_dia = doubles_in[0];
        } else {
          reply_header.number_of_doubles = 1;
          doubles_out[0] = dt_dia;
        }
        break;
      case 23:
        if(request_header.number_of_doubles == 1){
          eps2 = doubles_in[0];
        } else {
          reply_header.number_of_doubles = 1;
          doubles_out[0] = eps2;
        }
        break;
      case 24:
        if(request_header.number_of_ints == 1){
          flag_collision = ints_in[0];
        } else {
          reply_header.number_of_ints = 1;
          ints_out[0] = flag_collision;
        }
        break;
      default:
        reply_header.tag = -1;
    }
    
    reply_header.send(parent, rank);
    if(reply_header.number_of_doubles > 0) {
      parent.Send(doubles_out, reply_header.number_of_doubles, MPI_DOUBLE, 0, 999);
    }
    if(reply_header.number_of_ints > 0) {
      parent.Send(ints_out, reply_header.number_of_ints, MPI_INT, 0, 999);
    }
  }
}

int main(int argc, char *argv[])
{
  MPI::Init(argc, argv);
  
  run_loop();
  
  MPI_Finalize();
  return 0;
}
