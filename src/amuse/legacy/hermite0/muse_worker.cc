#include <iostream>
#include <cstdlib>

#include "muse_dynamics.h"
#include "parameters.h"
#include "local.h"


#include <mpi.h>

using namespace std;

double E0 = 0;

class message_header {

public:
	int tag;
	int id;
	int arg_1;
	int arg_2;
	int number_of_doubles;
	int number_of_ints;
	message_header(): tag(0), id(0), arg_1(0), arg_2(0), number_of_doubles(0), number_of_ints(0) {}

	void send(MPI::Intercomm & intercomm, int rank){
		int header[6];
		header[0] = tag;
		header[1] = id;
		header[2] = arg_1;
		header[3] = arg_2;
		header[4] = number_of_doubles;
		header[5] = number_of_ints;		
		intercomm.Send(header, 6, MPI_INT, 0, 999);
	}

	void recv(MPI::Intercomm & intercom, int rank) {
		int header[6];

		intercom.Recv(header, 6, MPI_INT, 0, rank);
		tag = header[0];
		id  = header[1];
		arg_1 = header[2];
		arg_2 = header[3];
		number_of_doubles = header[4];
		number_of_ints = header[5];
	}
};



void run_loop() {
	int rank = MPI::COMM_WORLD.Get_rank();

	MPI::Intercomm parent = MPI::COMM_WORLD.Get_parent();
	
        bool must_run_loop = true;
	while(must_run_loop) { 
		dynamics_state state;
		double doubles_in[20];
		double doubles_out[20];
		int ints_in[20];
		int ints_out[20];

		message_header request_header;
		message_header reply_header;

		request_header.recv(parent,rank);
		std::cout<<"aaa"<<std::endl;
		if(request_header.number_of_doubles > 0) {
			parent.Recv(doubles_in, request_header.number_of_doubles, MPI_DOUBLE, 0, rank);
		}
		if(request_header.number_of_ints > 0) {
			parent.Recv(ints_in, request_header.number_of_ints, MPI_INT, 0, rank);
		}
		std::cout<<"bbb"<<std::endl;

		reply_header.tag = request_header.tag;

		switch(request_header.tag) {
			case 0:
				must_run_loop = false;
				break;
			case 1:
				reply_header.number_of_ints = 1;
				ints_out[0] = setup_module(0,0);
				break;
			case 2:
				reply_header.number_of_ints = 1;
				ints_out[0] = cleanup_module();
				break;
			case 3:
				reply_header.number_of_ints = 1;
				ints_out[0] = initialize_particles(doubles_in[0]);
				break;
			case 4:
				reply_header.number_of_ints = 1;
				ints_out[0] = reinitialize_particles();
				break;
			case 5:
				reply_header.number_of_ints = 1;
				state.id = ints_in[0];
				state.mass = doubles_in[0];
				state.radius = doubles_in[1];
				state.x = doubles_in[2];
				state.y = doubles_in[3];
				state.z = doubles_in[4];
				state.vx = doubles_in[5];
				state.vy = doubles_in[6];
				state.vz = doubles_in[7];
				ints_out[0] =add_particle(state);
				break;
			case 6:
				reply_header.number_of_ints = 1;
				ints_out[0] =evolve(doubles_in[0], ints_in[0]);
				break;
			case 7:
				reply_header.number_of_ints = 1;
				ints_out[0] =get_number();
				break;
			case 8:
				state = get_state(ints_in[0]);
				reply_header.number_of_doubles = 8;
				ints_out[0] = state.id;
				reply_header.number_of_ints = 1;
				doubles_out[0] = state.mass;
				doubles_out[1] = state.radius;
				doubles_out[2] = state.x;
				doubles_out[3] = state.y;
				doubles_out[4] = state.z;
				doubles_out[5] = state.vx;
				doubles_out[6] = state.vy;
				doubles_out[7] = state.vz;
				break;
			case 20:
				if(request_header.number_of_doubles == 1) {
				    t = doubles_in[0];
				} else {
				    reply_header.number_of_doubles = 1;
				    doubles_out[0]  = t;
				}
				break;
			case 21:
				if(request_header.number_of_doubles == 1) {
				    dt_param = doubles_in[0];
				} else {
				    reply_header.number_of_doubles = 1;
				    doubles_out[0]  = dt_param;
				}
				break;
			case 22:
				if(request_header.number_of_doubles == 1) {
				    dt_dia = doubles_in[0];
				} else {
				    reply_header.number_of_doubles = 1;
				    doubles_out[0]  = dt_dia;
				}
				break;
			case 23:
				if(request_header.number_of_doubles == 1) {
				    eps2 = doubles_in[0];
				} else {
				    reply_header.number_of_doubles = 1;
				    doubles_out[0]  = eps2;
				}
				break;
			case 24:
				if(request_header.number_of_ints == 1) {
				    flag_collision = ints_in[0];
				} else {
				    reply_header.number_of_ints = 1;
				    ints_out[0] =flag_collision;
				}
				break;
			default:
               			reply_header.tag = -1;
                
		}

		
		reply_header.send(parent, rank);
		if( reply_header.number_of_doubles > 0) {
			parent.Send(doubles_out, reply_header.number_of_doubles, MPI_DOUBLE, 0, 999);
		}
		if( reply_header.number_of_ints > 0) {
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
