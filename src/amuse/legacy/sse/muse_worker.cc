#include <iostream>
#include <cstdlib>

extern "C" {

void initialize_( double * z_in, double * neta_in, double * bwind_in,
	double * hewind_in, double * sigma_in, int * ifflag_in,
	int * wdflag_in, int * bhflag_in, int * nsflag_in, int * mxns_in,
	double * pts1_in, double * pts2_in, double * pts3_in, int * status);


void evolve_(int * kw, double * mass, double * mt, double * r, double * lum, 
	double * mc, double * rc, double * menv, double * renv,
	double * ospin, double * epoch, double * tm, double * tphys, double * tphysf,
	int * kw1, double * mass1, double * mt1, double * r1, double * lum1, 
	double * mc1, double * rc1, double * menv1, double * renv1,
	double * ospin1, double * epoch1, double * tm1, double * tphys1, double * tphysf1);

void get_time_step_(int * kw, double * mass, double * age, double * mt, double * tm, double * epoch, double * dt);
}


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
	
	std::cout<<"SSE starting the loop!"<<std::endl;
        bool must_run_loop = true;
	while(must_run_loop) { 
		double doubles_in[20];
		double doubles_out[20];
		int ints_in[20];
		int ints_out[20];

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
				std::cout<<"n i:"<<request_header.number_of_ints<<std::endl;
				initialize_(
					&doubles_in[0],
					&doubles_in[1],
					&doubles_in[2],
					&doubles_in[3],
					&doubles_in[4],
					&ints_in[0],
					&ints_in[1],
					&ints_in[2],
					&ints_in[3],
					&ints_in[4],
					&doubles_in[5],
					&doubles_in[6],
					&doubles_in[7],
					&reply_header.arg_1
				);
				break;
			case 2:
				reply_header.number_of_doubles = 1;
				get_time_step_(
					&ints_in[0],
					&doubles_in[0],
					&doubles_in[1],
					&doubles_in[2],
					&doubles_in[3],
					&doubles_in[4],
					&doubles_out[0]
				);
				break;
			case 3:
				reply_header.number_of_doubles = 13;
				reply_header.number_of_ints = 1;
				evolve_(
					&ints_in[0],
					&doubles_in[0],
					&doubles_in[1],
					&doubles_in[2],
					&doubles_in[3],
					&doubles_in[4],
					&doubles_in[5],
					&doubles_in[6],
					&doubles_in[7],
					&doubles_in[8],
					&doubles_in[9],
					&doubles_in[10],
					&doubles_in[11],
					&doubles_in[12],
					&ints_out[0],
					&doubles_out[0],
					&doubles_out[1],
					&doubles_out[2],
					&doubles_out[3],
					&doubles_out[4],
					&doubles_out[5],
					&doubles_out[6],
					&doubles_out[7],
					&doubles_out[8],
					&doubles_out[9],
					&doubles_out[10],
					&doubles_out[11],
					&doubles_out[12]
				);
				break;
			default:
				reply_header.tag = -1;
		}

		
		reply_header.send(parent, rank);
		if( reply_header.number_of_doubles > 0) {
			parent.Send(doubles_out, reply_header.number_of_doubles, MPI_DOUBLE, 0, 999);
		}
		if( reply_header.number_of_ints > 0) {
			parent.Send(ints_out, reply_header.number_of_ints, MPI_DOUBLE, 0, 999);
		}
	}
   
	std::cout<<"SSE loop ended!"<<std::endl;

}
 
int main(int argc, char *argv[])
{
	MPI::Init(argc, argv);

	run_loop();    

	MPI_Finalize();
	return 0;
}
