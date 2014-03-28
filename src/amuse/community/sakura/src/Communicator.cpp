#include "Communicator.h"

//////////////////////////////////////////////////////////////////////////
// Constructors
//////////////////////////////////////////////////////////////////////////
Communicator::Communicator() {
  ;
}
//////////////////////////////////////////////////////////////////////////
// MPI
//////////////////////////////////////////////////////////////////////////
void Communicator::start_mpi() {
  //MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  char processor_name[MPI_MAX_PROCESSOR_NAME];
  int name_len;
  MPI_Get_processor_name(processor_name, &name_len);
  name = processor_name;  
  new_comm = MPI_COMM_WORLD;
}
void Communicator::stop_mpi() {
  //MPI_Finalize();
}
void Communicator::make_two_groups() {
  /* Extract the original group handle */
  MPI_Comm_group(MPI_COMM_WORLD, &orig_group);

  /* Divide tasks into two distinct groups based upon rank */
  vector<int> v1, v2;
  for(int i=0; i<size/2; i++) v1.push_back(i);
  for(int i=size/2; i<size; i++) v2.push_back(i);
  int* ranks1 = &v1[0];
  int* ranks2 = &v2[0];
  MPI_Group_incl(orig_group, size/2, ranks1, &new_group1);
  MPI_Group_incl(orig_group, size/2, ranks2, &new_group2);

  /* Create new communicators */
  MPI_Comm_create(MPI_COMM_WORLD, new_group1, &new_comm1);
  MPI_Comm_create(MPI_COMM_WORLD, new_group2, &new_comm2);

  /* Keep communicators different than MPI_COMM_NULL based on the rank */
  if(rank < size/2) {
    new_comm = new_comm1;
    group = 1;
  }
  else {
    new_comm = new_comm2;
    group = 2;
  }

  /* Get new ranks */
  MPI_Comm_rank(new_comm, &rank);
  MPI_Comm_size(new_comm, &size);
}
void Communicator::free_groups() {
  /* Free grouping */
  MPI_Comm_free(&new_comm);
  MPI_Group_free(&new_group1);
  MPI_Group_free(&new_group2);
}

//////////////////////////////////////////////////////////////////////////
// Work
//////////////////////////////////////////////////////////////////////////
void Communicator::divide_work(int numElements) {
  numStar = numElements;
  int myNumber = (int)(numStar/size);
  int restNumber = numStar - myNumber*size;
  for(int i=0; i<restNumber; i++) {
    if(i == rank) myNumber++;
  }
  gather(myNumber, all_numbers, new_comm);
  bcast(all_numbers, new_comm);
  my_begin = 0;
  for(int i=0; i<rank; i++) {
    my_begin += all_numbers[i];
  }
  my_end = my_begin + myNumber;
  gather(my_begin, all_begins, new_comm);
  bcast(all_begins, new_comm);
}
//////////////////////////////////////////////////////////////////////////
// Getters
//////////////////////////////////////////////////////////////////////////
string Communicator::get_name() {
  return name;
}

int Communicator::get_rank() {
  return rank;
}
int Communicator::get_size() {
  return size;
}
int Communicator::get_my_begin() {
  return my_begin;
}
int Communicator::get_my_end() {
  return my_end;
}

int Communicator::get_group() {
  return group;
}
MPI_Comm Communicator::get_comm_world() {
  return MPI_COMM_WORLD;
}
MPI_Comm Communicator::get_comm_group1() {
  return new_comm1;
}
MPI_Comm Communicator::get_comm_group2() {
  return new_comm2;
}
//////////////////////////////////////////////////////////////////////////
// INT
//////////////////////////////////////////////////////////////////////////
void Communicator::bcast(int &x) {
  MPI_Bcast(&x, 1, MPI_INT, 0, MPI_COMM_WORLD);   
}
void Communicator::bcast(vector<int> &x) {
  int N = x.size();
  bcast(N);
  x.resize(N);
  MPI_Bcast(&x.front(), N, MPI_INT, 0, MPI_COMM_WORLD);
}
void Communicator::gather(int &x, vector<int> &y) {
  y.resize(size);
  MPI_Gather(&x, 1, MPI_INT, &y.front(), 1, MPI_INT, 0, MPI_COMM_WORLD);  
}
void Communicator::gather(vector<int> &x, vector<int> &y) {
  int N = x.size();
  vector<int> N_col;
  gather(N, N_col);
  vector<int> i_col(size);
  for(int i=1; i<size; i++) {
    for(int j=i; j<size; j++) {
      i_col[j] += N_col[i-1];
    }
  }
  int N_sum = 0;
  for(int i=0; i<N_col.size(); i++) N_sum += N_col[i];
  y.resize(N_sum);
  MPI_Gatherv(&x.front(), N, MPI_INT, &y.front(), &N_col.front(), &i_col.front(), MPI_INT, 0, MPI_COMM_WORLD); 
}
void Communicator::join(vector<int> &x, vector<int> &y) {
  int M = x.size();
  y.resize(M);
  int N = M/numStar;
  int begin = my_begin*N;
  int number = all_numbers[rank]*N;
  vector<int> all_begin = all_begins;
  vector<int> all_number = all_numbers;
  for(int i=0; i<size; i++) {
    all_begin[i] *= N;
    all_number[i] *= N;
  }
  MPI_Gatherv(&x[begin], number, MPI_INT, &y.front(), &all_number.front(), &all_begin.front(), MPI_INT, 0, MPI_COMM_WORLD);
}

void Communicator::bcast(int &x, MPI_Comm comm) {
  MPI_Bcast(&x, 1, MPI_INT, 0, comm); 
}
void Communicator::bcast(vector<int> &x, MPI_Comm comm) {
  int N = x.size();
  bcast(N, comm);
  x.resize(N);
  MPI_Bcast(&x.front(), N, MPI_INT, 0, comm);
}
void Communicator::gather(int &x, vector<int> &y, MPI_Comm comm) {
  y.resize(size);
  MPI_Gather(&x, 1, MPI_INT, &y.front(), 1, MPI_INT, 0, comm);  
}
void Communicator::gather(vector<int> &x, vector<int> &y, MPI_Comm comm) {
  int N = x.size();
  vector<int> N_col;
  gather(N, N_col);
  vector<int> i_col(size);
  for(int i=1; i<size; i++) {
    for(int j=i; j<size; j++) {
      i_col[j] += N_col[i-1];
    }
  }
  int N_sum = 0;
  for(int i=0; i<N_col.size(); i++) N_sum += N_col[i];
  y.resize(N_sum);
  MPI_Gatherv(&x.front(), N, MPI_INT, &y.front(), &N_col.front(), &i_col.front(), MPI_INT, 0, comm); 
}
void Communicator::join(vector<int> &x, vector<int> &y, MPI_Comm comm) {
  int M = x.size();
  y.resize(M);
  int N = M/numStar;
  int begin = my_begin*N;
  int number = all_numbers[rank]*N;
  vector<int> all_begin = all_begins;
  vector<int> all_number = all_numbers;
  for(int i=0; i<size; i++) {
    all_begin[i] *= N;
    all_number[i] *= N;
  }
  MPI_Gatherv(&x[begin], number, MPI_INT, &y.front(), &all_number.front(), &all_begin.front(), MPI_INT, 0, comm);
}
/////////////////////////////////////////////////////////////////////////
// FLOAT
//////////////////////////////////////////////////////////////////////////
void Communicator::bcast(float &x) {
  MPI_Bcast(&x, 1, MPI_FLOAT, 0, MPI_COMM_WORLD);   
}
void Communicator::bcast(vector<float> &x) {
  int N = x.size();
  bcast(N);
  x.resize(N);
  MPI_Bcast(&x.front(), N, MPI_FLOAT, 0, MPI_COMM_WORLD);
}
void Communicator::gather(float &x, vector<float> &y) {
  y.resize(size);
  MPI_Gather(&x, 1, MPI_FLOAT, &y.front(), 1, MPI_FLOAT, 0, MPI_COMM_WORLD);  
}
void Communicator::gather(vector<float> &x, vector<float> &y) {
  int N = x.size();
  vector<int> N_col;
  gather(N, N_col);
  vector<int> i_col(size);
  for(int i=1; i<size; i++) {
    for(int j=i; j<size; j++) {
      i_col[j] += N_col[i-1];
    }
  }
  int N_sum = 0;
  for(int i=0; i<N_col.size(); i++) N_sum += N_col[i];
  y.resize(N_sum);
  MPI_Gatherv(&x.front(), N, MPI_FLOAT, &y.front(), &N_col.front(), &i_col.front(), MPI_FLOAT, 0, MPI_COMM_WORLD); 
}
void Communicator::join(vector<float> &x, vector<float> &y) {
  int M = x.size();
  y.resize(M);
  int N = M/numStar;
  int begin = my_begin*N;
  int number = all_numbers[rank]*N;
  vector<int> all_begin = all_begins;
  vector<int> all_number = all_numbers;
  for(int i=0; i<size; i++) {
    all_begin[i] *= N;
    all_number[i] *= N;
  }
  MPI_Gatherv(&x[begin], number, MPI_FLOAT, &y.front(), &all_number.front(), &all_begin.front(), MPI_FLOAT, 0, MPI_COMM_WORLD);
}

void Communicator::bcast(float &x, MPI_Comm comm) {
  MPI_Bcast(&x, 1, MPI_FLOAT, 0, comm);   
}
void Communicator::bcast(vector<float> &x, MPI_Comm comm) {
  int N = x.size();
  bcast(N);
  x.resize(N);
  MPI_Bcast(&x.front(), N, MPI_FLOAT, 0, comm);
}
void Communicator::gather(float &x, vector<float> &y, MPI_Comm comm) {
  y.resize(size);
  MPI_Gather(&x, 1, MPI_FLOAT, &y.front(), 1, MPI_FLOAT, 0, comm);  
}
void Communicator::gather(vector<float> &x, vector<float> &y, MPI_Comm comm) {
  int N = x.size();
  vector<int> N_col;
  gather(N, N_col);
  vector<int> i_col(size);
  for(int i=1; i<size; i++) {
    for(int j=i; j<size; j++) {
      i_col[j] += N_col[i-1];
    }
  }
  int N_sum = 0;
  for(int i=0; i<N_col.size(); i++) N_sum += N_col[i];
  y.resize(N_sum);
  MPI_Gatherv(&x.front(), N, MPI_FLOAT, &y.front(), &N_col.front(), &i_col.front(), MPI_FLOAT, 0, comm); 
}
void Communicator::join(vector<float> &x, vector<float> &y, MPI_Comm comm) {
  int M = x.size();
  y.resize(M);
  int N = M/numStar;
  int begin = my_begin*N;
  int number = all_numbers[rank]*N;
  vector<int> all_begin = all_begins;
  vector<int> all_number = all_numbers;
  for(int i=0; i<size; i++) {
    all_begin[i] *= N;
    all_number[i] *= N;
  }
  MPI_Gatherv(&x[begin], number, MPI_FLOAT, &y.front(), &all_number.front(), &all_begin.front(), MPI_FLOAT, 0, comm);
}
/////////////////////////////////////////////////////////////////////////
// DOUBLE
//////////////////////////////////////////////////////////////////////////
void Communicator::bcast(double &x) {
  MPI_Bcast(&x, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);   
}
void Communicator::bcast(vector<double> &x) {
  int N = x.size();
  bcast(N);
  x.resize(N);
  MPI_Bcast(&x.front(), N, MPI_DOUBLE, 0, MPI_COMM_WORLD);
}
void Communicator::gather(double &x, vector<double> &y) {
  y.resize(size);
  MPI_Gather(&x, 1, MPI_DOUBLE, &y.front(), 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);  
}
void Communicator::gather(vector<double> &x, vector<double> &y) {
  int N = x.size();
  vector<int> N_col;
  gather(N, N_col);
  vector<int> i_col(size);
  for(int i=1; i<size; i++) {
    for(int j=i; j<size; j++) {
      i_col[j] += N_col[i-1];
    }
  }
  int N_sum = 0;
  for(int i=0; i<N_col.size(); i++) N_sum += N_col[i];
  y.resize(N_sum);
  MPI_Gatherv(&x.front(), N, MPI_DOUBLE, &y.front(), &N_col.front(), &i_col.front(), MPI_DOUBLE, 0, MPI_COMM_WORLD); 
}
void Communicator::join(vector<double> &x, vector<double> &y) {
  int M = x.size();
  y.resize(M);
  int N = M/numStar;
  int begin = my_begin*N;
  int number = all_numbers[rank]*N;
  vector<int> all_begin = all_begins;
  vector<int> all_number = all_numbers;
  for(int i=0; i<size; i++) {
    all_begin[i] *= N;
    all_number[i] *= N;
  }
  MPI_Gatherv(&x[begin], number, MPI_DOUBLE, &y.front(), &all_number.front(), &all_begin.front(), MPI_DOUBLE, 0, MPI_COMM_WORLD);
}

void Communicator::bcast(double &x, MPI_Comm comm) {
  MPI_Bcast(&x, 1, MPI_DOUBLE, 0, comm);   
}
void Communicator::bcast(vector<double> &x, MPI_Comm comm) {
  int N = x.size();
  bcast(N);
  x.resize(N);
  MPI_Bcast(&x.front(), N, MPI_DOUBLE, 0, comm);
}
void Communicator::gather(double &x, vector<double> &y, MPI_Comm comm) {
  y.resize(size);
  MPI_Gather(&x, 1, MPI_DOUBLE, &y.front(), 1, MPI_DOUBLE, 0, comm);  
}
void Communicator::gather(vector<double> &x, vector<double> &y, MPI_Comm comm) {
  int N = x.size();
  vector<int> N_col;
  gather(N, N_col);
  vector<int> i_col(size);
  for(int i=1; i<size; i++) {
    for(int j=i; j<size; j++) {
      i_col[j] += N_col[i-1];
    }
  }
  int N_sum = 0;
  for(int i=0; i<N_col.size(); i++) N_sum += N_col[i];
  y.resize(N_sum);
  MPI_Gatherv(&x.front(), N, MPI_DOUBLE, &y.front(), &N_col.front(), &i_col.front(), MPI_DOUBLE, 0, comm); 
}
void Communicator::join(vector<double> &x, vector<double> &y, MPI_Comm comm) {
  int M = x.size();
  y.resize(M);
  int N = M/numStar;
  int begin = my_begin*N;
  int number = all_numbers[rank]*N;
  vector<int> all_begin = all_begins;
  vector<int> all_number = all_numbers;
  for(int i=0; i<size; i++) {
    all_begin[i] *= N;
    all_number[i] *= N;
  }
  MPI_Gatherv(&x[begin], number, MPI_DOUBLE, &y.front(), &all_number.front(), &all_begin.front(), MPI_DOUBLE, 0, comm);
}



