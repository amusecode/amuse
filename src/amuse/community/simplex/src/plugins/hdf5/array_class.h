/*
  Container class implementing a n-dimensional array as a 1d array
  for use with h5w HDF5 wrapper...
*/

#ifndef ARR_1D_H
#define ARR_1D_H

#include <iostream>
#include <stdlib.h>
#include <assert.h>
#include <string.h>
#include <typeinfo>
using namespace std;

template<typename T> void set_datatype(char *type, T var);

// Container class
template<typename T> class arr_1D
{
public:
  arr_1D() {arr=NULL; rank=0; dims=NULL;};
  arr_1D(int r, int *d);
  ~arr_1D();
  
  // initializes the object again, deleting previous data and shape
  void reinit(int r, int *d);

  // Indexing operators
  T& operator() ( int i );
  T& operator() ( int i, int j );
  T& operator() ( int i, int j, int k );
  // add more for arrays with n>3

  T* get_arr() {return arr;};  // returns a pointer to the 1D array - pass this to HDF5 reading/writing routines
  char* get_dtype() {return type;};
  int* get_dims() {return dims;};
  int get_rank() {return rank;};
  
private:
  T *arr;
  int *dims;
  int rank;
  char type[20];
  
};

// Constructor
template <typename T> arr_1D<T>::arr_1D (int r, int *d)
{
  int i;
  int tot_length;
  
  rank = r;
  if (dims)
    delete [] dims;
  dims = new int[rank];
  for (i=0; i< d[i]; i++)
    dims[i] = d[i];
    
  tot_length=1;
  for (i=0; i<rank; i++)
    tot_length*=dims[i];
  
  //allocate 1D array
  arr = new T [tot_length];
  set_datatype(type, arr[0]);

}

//Deconstructor
template <typename T> arr_1D<T>::~arr_1D ()
{
  if (dims)
    delete [] dims;
  dims = NULL; 

  if (arr)
    delete [] arr;
  arr = NULL;
}

// Reinitialize the instance...
template <typename T> void arr_1D<T>::reinit (int r, int *d)
{
  int i;
  int tot_length;
  
  rank = r;
  if (dims){
    delete [] dims;
    dims = NULL;
  }
  dims = new int[rank];
  for (i=0; i<rank; i++)
    dims[i] = d[i];
  
  tot_length=1;
  for (i=0; i<rank; i++)
    tot_length*=dims[i];
  
  if (arr)
    delete [] arr;
  //allocate 1D array
  arr = new T [tot_length];
  set_datatype(type, arr[0]);

}



// Overloading () operator - one for every dimesnion of array...
// 1D
template <typename T>
T& arr_1D<T>::operator() ( int i )
{
  if (rank == 1)
    {
      // check boundaries
      assert(i>=0 && i<dims[0]);
      return arr[i];
    }
  else
    {
      cout << "This array has " << rank << " dimensions!" << endl;
      exit(-1);
    }
}

// 2D
template <typename T>
T& arr_1D<T>::operator() ( int i , int j)
{
  if (rank == 2)
    {
      assert(i>=0 && i<dims[0]);
      assert(j>=0 && j<dims[1]);

      int count;
      
      count = i*dims[1] + j;
      return arr[count];
    }
  else
    {
      cout << "This array has " << rank << " dimenstions!" << endl;
      exit(-1);
    }

}

// 3D
template <typename T>
T& arr_1D<T>::operator() ( int i , int j, int k)
{
  if (rank == 3)
    {
      assert(i>=0 && i<dims[0]);
      assert(j>=0 && j<dims[1]);
      assert(k>=0 && k<dims[2]);

      int count;

      count = i*dims[1] + j*dims[2] + k;
      return arr[count];
    }
  else
    {
      cout << "This array has " << rank << " dimenstions!" << endl;
      exit(-1);
    }

}

// finds which hdf5 datatype should be used to store the variable passed as the argument !
template<typename T> void set_datatype(char *type, T var)
{
  char dtype[20];
  char t_ind;  
  strcpy(dtype, typeid(var).name());

  // check if the variable is a pointer to arr_1D class  
  if(strstr(dtype,"P6arr_1D"))
    {
      t_ind = 'P';
      strcpy(type,"arr_1D");
    }
  else // if not, continue normaly...
    t_ind = dtype[strlen(dtype)-1];
  
  if (t_ind=='i')
    strcpy(type,"int");
  else if (t_ind=='c')
    strcpy(type, "char");
  else if (t_ind=='f')
    strcpy(type, "float");
  else if (t_ind=='d')
    strcpy(type, "double");
  else if (t_ind=='s')   // short int is short
    strcpy(type, "short");
  else if (t_ind=='j')
    strcpy(type, "unsigned int");
  else if (t_ind=='t')
    strcpy(type, "unsigned short");
  else if (t_ind=='l')   // long int is long 
    strcpy(type, "long");
  else if (t_ind=='x')
    strcpy(type, "long long");
  else if (t_ind=='e')
    strcpy(type, "long double");
  else if (t_ind=='h')
    strcpy(type, "unsigned char");
  else if (t_ind=='y')
    strcpy(type, "unsigned long long");
  else if (t_ind=='m')
    strcpy(type, "unsigned long");
  
  else
    {
      cout << "Datatype string not recognized! Exiting..." << endl;
      exit(-1);
    }
}



#endif
