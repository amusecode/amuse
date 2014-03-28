#ifndef NOMPI
#include <mpi.h>
#endif

using namespace std;
#include <vector>
#include <string>

#ifndef __Communicator_h
#define __Communicator_h

class Communicator {
  int rank, size;
  string name;

  int numStar;
  int my_begin, my_end;
  vector<int> all_begins, all_numbers;

  int group, sendbuf_rank, recvbuf_sum;
  MPI_Group orig_group, new_group1, new_group2;  
  MPI_Comm new_comm, new_comm1, new_comm2;

  public:

  // constructors
  Communicator();

  // mpi
  void start_mpi();
  void stop_mpi();

  void make_two_groups();
  void free_groups();

  // work
  void divide_work(int numElements);

  // getters
  string get_name();

  int get_rank();
  int get_size();
  int get_my_begin();
  int get_my_end();

  int get_group();
  MPI_Comm get_comm_world();
  MPI_Comm get_comm_group1();
  MPI_Comm get_comm_group2();

  // int
  void bcast(int &x);
  void bcast(vector<int> &x);
  void gather(int &x, vector<int> &y);
  void gather(vector<int> &x, vector<int> &y);
  void join(vector<int> &x, vector<int> &y);

  void bcast(int &x, MPI_Comm comm);
  void bcast(vector<int> &x, MPI_Comm comm);
  void gather(int &x, vector<int> &y, MPI_Comm comm);
  void gather(vector<int> &x, vector<int> &y, MPI_Comm comm);
  void join(vector<int> &x, vector<int> &y, MPI_Comm comm);

  // float
  void bcast(float &x);
  void bcast(vector<float> &x);
  void gather(float &x, vector<float> &y);
  void gather(vector<float> &x, vector<float> &y);
  void join(vector<float> &x, vector<float> &y);

  void bcast(float &x, MPI_Comm comm);
  void bcast(vector<float> &x, MPI_Comm comm);
  void gather(float &x, vector<float> &y, MPI_Comm comm);
  void gather(vector<float> &x, vector<float> &y, MPI_Comm comm);
  void join(vector<float> &x, vector<float> &y, MPI_Comm comm);

  // double
  void bcast(double &x);
  void bcast(vector<double> &x);
  void gather(double &x, vector<double> &y);
  void gather(vector<double> &x, vector<double> &y);
  void join(vector<double> &x, vector<double> &y);

  void bcast(double &x, MPI_Comm comm);
  void bcast(vector<double> &x, MPI_Comm comm);
  void gather(double &x, vector<double> &y, MPI_Comm comm);
  void gather(vector<double> &x, vector<double> &y, MPI_Comm comm);
  void join(vector<double> &x, vector<double> &y, MPI_Comm comm);
};

#endif


