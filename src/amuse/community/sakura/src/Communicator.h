#ifndef NOMPI
#include <mpi.h>
#endif

#include <vector>
#include <string>

#ifndef __Communicator_h
#define __Communicator_h

class Communicator {
  int rank, size;
  std::string name;

  int numStar;
  int my_begin, my_end;
  std::vector<int> all_begins, all_numbers;

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
  std::string get_name();

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
  void bcast(std::vector<int> &x);
  void gather(int &x, std::vector<int> &y);
  void gather(std::vector<int> &x, std::vector<int> &y);
  void join(std::vector<int> &x, std::vector<int> &y);

  void bcast(int &x, MPI_Comm comm);
  void bcast(std::vector<int> &x, MPI_Comm comm);
  void gather(int &x, std::vector<int> &y, MPI_Comm comm);
  void gather(std::vector<int> &x, std::vector<int> &y, MPI_Comm comm);
  void join(std::vector<int> &x, std::vector<int> &y, MPI_Comm comm);

  // float
  void bcast(float &x);
  void bcast(std::vector<float> &x);
  void gather(float &x, std::vector<float> &y);
  void gather(std::vector<float> &x, std::vector<float> &y);
  void join(std::vector<float> &x, std::vector<float> &y);

  void bcast(float &x, MPI_Comm comm);
  void bcast(std::vector<float> &x, MPI_Comm comm);
  void gather(float &x, std::vector<float> &y, MPI_Comm comm);
  void gather(std::vector<float> &x, std::vector<float> &y, MPI_Comm comm);
  void join(std::vector<float> &x, std::vector<float> &y, MPI_Comm comm);

  // double
  void bcast(double &x);
  void bcast(std::vector<double> &x);
  void gather(double &x, std::vector<double> &y);
  void gather(std::vector<double> &x, std::vector<double> &y);
  void join(std::vector<double> &x, std::vector<double> &y);

  void bcast(double &x, MPI_Comm comm);
  void bcast(std::vector<double> &x, MPI_Comm comm);
  void gather(double &x, std::vector<double> &y, MPI_Comm comm);
  void gather(std::vector<double> &x, std::vector<double> &y, MPI_Comm comm);
  void join(std::vector<double> &x, std::vector<double> &y, MPI_Comm comm);
};

#endif


