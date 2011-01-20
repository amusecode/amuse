/*
  This header file contains wrapper routines for reading and writing data using HDF5
*/
#ifndef H5W_H
#define H5W_H

#include <iostream>
#include <stdlib.h>
#include <typeinfo>
#include <string>
#include <stdlib.h>
#include "array_class.h"

//extern "C" {
#include "hdf5.h"
//}
using namespace std;


// routines for finding the HDF5 datatype of a given variable or from a character flag
int get_h5_datatype(const char *name);



/*
  Main class holding all the HDF5 routines...
*/


class h5w{

  public:

    h5w() {};
    h5w(char fname[]);
    h5w(char fname[],char mode);

    void open(char fname[]); 
    void close(void); 
    void make_dataset(const char *name, const char *type, int rank, int *dims);  // makes a dataset on atomic data
    void make_group(const char *name);                        // makes a group with absolute pathname "name"
    void make_group(const char *name, int size_hint);         // makes a group with absolute pathname "name"

    template<typename T> void write_attr(const char *dataset_name, const char *attr_name, T *data);  // all single values by reference
    template<typename T> void write_attr(const char *dataset_name, const char *attr_name, T data);  // all single values by value
    void write_attr(const char *dataset_name, const char *attr_name, char *data);        // for strings
    void write_attr(const char *dataset_name, const char *attr_name, const char *data);  // for strings
    template<typename T> void write_attr(const char *dataset_name, const char *attr_name, arr_1D<T> *data); // for arrays, arr_1D instance

    template<typename T> void read_attr(const char *dataset_name, const char *attr_name, T *data);  // all single values by reference
    void read_attr(const char *dataset_name, const char *attr_name, char *data);  // all single values by reference
    template<typename T> void read_attr(const char *dataset_name, const char *attr_name, arr_1D<T> *data); // for arrays, arr_1D instance

    /* Warning! Reading and writing of data is done through the arr_1D container class! */
    template<typename T> void write_data(const char *dataset_path, unsigned long long int *offset, arr_1D<T>* data);
    template<typename T> void read_data(const char *dataset_path, int *offset, arr_1D<T>* data);

    string get_data_type(char *dataset_name);
    // this routine calls a set of dummy routines and will read any type of numberic data to a
    // double container
    void read_data_as_double(char *dataset_name, int *offset, arr_1D<double>* dbl_arr);

  
  private:

    // wrapper routines for reading any type of numeric data to a double type container
    void read_int_as_double(char *dataset_name, int *offset, arr_1D<double>* dbl_arr);
    void read_uint_as_double(char *dataset_name, int *offset, arr_1D<double>* dbl_arr);
    void read_long_as_double(char *dataset_name, int *offset, arr_1D<double>* dbl_arr);
    void read_llong_as_double(char *dataset_name, int *offset, arr_1D<double>* dbl_arr);
    void read_ulong_as_double(char *dataset_name, int *offset, arr_1D<double>* dbl_arr);
    void read_ullong_as_double(char *dataset_name, int *offset, arr_1D<double>* dbl_arr);
    void read_float_as_double(char *dataset_name, int *offset, arr_1D<double>* dbl_arr);
    void read_short_as_double(char *dataset_name, int *offset, arr_1D<double>* dbl_arr);
    void read_ushort_as_double(char *dataset_name, int *offset, arr_1D<double>* dbl_arr);
    void read_char_as_double(char *dataset_name, int *offset, arr_1D<double>* dbl_arr);
    void read_uchar_as_double(char *dataset_name, int *offset, arr_1D<double>* dbl_arr);

    hid_t file;
    char filename[250];
  
};

// write_attr overload to write single variable by refenrece
template<typename T> void h5w::write_attr(const char *dataset_name, const char *attr_name, T *data)
{
  herr_t status;
  char type_tmp[20];
  hid_t group, dataset, object; // first 2 are temporary vars to see if the object is a group or dataset
  hid_t attr_datatype, attr_dataspace, attr;
  hid_t datatype_in;
  bool object_is_group=1;

  group = H5Gopen(file, dataset_name);
  if (group < 0)
    {
    dataset = H5Dopen(file, dataset_name);
    if (dataset < 0)
      {
        cout << "Path given for attribute writting not a group or a dataset!" << endl;
        cout << "Exiting!" << endl;
        exit(-1);
      }
    object = dataset;
    H5Gclose(group);
    object_is_group = 0;
    }
  else
    {
      object = group;
      object_is_group = 1;
    }

  // find datatype of input data
  set_datatype(type_tmp, data);
  datatype_in = H5Tcopy(get_h5_datatype(type_tmp));
  
  attr = H5Aopen_name(object, attr_name);
  if (attr > 0) // attribute already exists, delete it and write new stuff...
    H5Adelete(object,attr_name);
  
  attr_dataspace = H5Screate(H5S_SCALAR);
  // create new attribute
  attr = H5Acreate(object, attr_name, datatype_in, attr_dataspace, H5P_DEFAULT);

  //write 
  status = H5Awrite(attr, datatype_in, data);
  if (status < 0)
    {
      cout<< "Writing atribute: " << attr_name << " to object: " << dataset_name << " failed !" << endl << "Exiting!" << endl;
      exit(-1);
    }

  if(object_is_group){
    H5Gclose(group);
    H5Aclose(attr);
  }else{
    H5Dclose(dataset);
    H5Aclose(attr);
  }
  
}


// write_arr overload to pass single variable by value ...
template<typename T> void h5w::write_attr(const char *dataset_name, const char *attr_name, T data)
{
  herr_t status;
  char type_tmp[20];
  hid_t group, dataset, object; // first 2 are temporary vars to see if the object is a group or dataset
  hid_t attr_dataspace, attr;
  hid_t datatype_in;
  bool object_is_group=1;


  group = H5Gopen(file, dataset_name);
  if (group < 0)
    {
    dataset = H5Dopen(file, dataset_name);
    if (dataset < 0)
      {
        cout << "Path given for attribute writting not a group or a dataset!" << endl;
        cout << "Exiting!" << endl;
        exit(-1);
      }
    object = dataset;
    H5Gclose(group);
    object_is_group = 0;
    }
  else
    {
      object = group;
      object_is_group = 1;
    }

  // find datatype of input data
  set_datatype(type_tmp, data);
  datatype_in = H5Tcopy(get_h5_datatype(type_tmp));
  
  attr = H5Aopen_name(object, attr_name);
  if (attr > 0) // attribute already exists, delete it and write new stuff...
    H5Adelete(object,attr_name);
  
  attr_dataspace = H5Screate(H5S_SCALAR);
  // create new attribute
  attr = H5Acreate(object, attr_name, datatype_in, attr_dataspace, H5P_DEFAULT);

  //write 
  status = H5Awrite(attr, datatype_in, &data);
  if (status < 0)
    {
      cout<< "Writing atribute: " << attr_name << " to object: " << dataset_name << " failed !" << endl << "Exiting!" << endl;
      exit(-1);
    }

  if(object_is_group){
    H5Gclose(group);
    H5Aclose(attr);
  }else{
    H5Dclose(dataset);
    H5Aclose(attr);
  }
  
}

// read_attr overload to read single variable by refenrece
template<typename T> void h5w::read_attr(const char *dataset_name, const char *attr_name, T *data)
{
  herr_t status;
  char type_tmp[20];
  hid_t group, dataset, object; // first 2 are temporary vars to see if the object is a group or dataset
  hid_t attr_datatype, attr;
  hid_t datatype_in;

  group = H5Gopen(file, dataset_name);
  if (group < 0)
    {
    dataset = H5Dopen(file, dataset_name);
    if (dataset < 0)
      {
        cout << "Path given for attribute writting not a group or a dataset!" << endl;
        cout << "Exiting!" << endl;
        exit(-1);
      }
    object = dataset;
    H5Gclose(group);
    }
  else
    {
      object = group;
      //H5Dclose(dataset);
    }

  // find datatype of input data
  set_datatype(type_tmp, data);
  datatype_in = H5Tcopy(get_h5_datatype(type_tmp));
    
  attr = H5Aopen_name(object, attr_name);
  if (attr < 0) // attribute doesn't exist
    {
      cout << "Attribute " << attr_name << "could not be opened in object: " << dataset_name << endl;
      cout << "Exiting!" << endl;
      exit(-1);
    }
  
  attr_datatype = H5Aget_type(attr);
  if (attr_datatype < 0) // attribute doesn't exist
    {
      cout << "Could not get type of attribute: " << attr_name << endl;
      cout << "Exiting!" << endl;
      exit(-1);
    }
  if (H5Tequal(datatype_in, attr_datatype) <= 0)   // compare datatypes 
    {
      cout << "Type of variable passed not the same as type of data to be read in attribute: " << attr_name << endl;
      cout << "Exiting!" << endl;
      exit(-1);
    }
 

  //read 
  status = H5Aread(attr, datatype_in, data);
  if (status < 0)
    {
      cout<< "Reading atribute: " << attr_name << " from object: " << dataset_name << " failed !" << endl << "Exiting!" << endl;
      exit(-1);
    }

}

// read_attr overload to read an array of values into a arr_1D<T> instance reference
template<typename T> void h5w::read_attr(const char *dataset_name, const char *attr_name, arr_1D<T> *data)
{
  int i;
  herr_t status;
  hid_t group, dataset, object; // first 2 are temporary vars to see if the object is a group or dataset
  hid_t attr_datatype, attr_dataspace, attr;
  hsize_t *attr_dims;
  int attr_rank;
  hid_t datatype_in;
  
  group = H5Gopen(file, dataset_name);
  if (group < 0)
    {
    dataset = H5Dopen(file, dataset_name);
    if (dataset < 0)
      {
        cout << "Path given for attribute writting not a group or a dataset!" << endl;
        cout << "Exiting!" << endl;
        exit(-1);
      }
    object = dataset;
    H5Gclose(group);
    }
  else
    {
      object = group;
      H5Dclose(dataset);
    }

  // find datatype of input data
  datatype_in = H5Tcopy(get_h5_datatype(data->get_dtype()));
  
  attr = H5Aopen_name(object, attr_name);
  if (attr < 0) // no such attribute
    {
      cout << "Attribute " << attr_name << "could not be opened in object: " << dataset_name << endl;
      cout << "Exiting!" << endl;
      exit(-1);
    }

  attr_datatype = H5Aget_type(attr);
  if (attr_datatype < 0) // attribute doesn't exist
    {
      cout << "Could not get type of attribute: " << attr_name << endl;
      cout << "Exiting!" << endl;
      exit(-1);
    }
  if (H5Tequal(datatype_in, attr_datatype) <= 0)   // compare datatypes 
    {
      cout << "Type of variable passed not the same as type of data to be read in attribute: " << attr_name << endl;
      cout << "Exiting!" << endl;
      exit(-1);
    }
  

  attr_dataspace = H5Aget_space(attr);
  attr_rank = H5Sget_simple_extent_dims(attr_dataspace, attr_dims, NULL);
  if ((data->get_rank()) != (attr_rank))
    {
      cout << "Data container for attribute:" << attr_name << " does not have the same rank as the data" << endl;
      cout << "Exiting!" << endl;
      exit(-1);
    }
  
  int tmp_x, tmp_y;
  for (i=0; i< data->get_rank(); i++)
    {   
      tmp_x =  data->get_dims()[i];
      tmp_y =  attr_dims[i];
      //      if (attr_dims[i] != data->get_dims()[i])   !!! WHY IS THIS A BUG !?!?!
      if (tmp_x != tmp_y)
        {
          cout << "Data container for attribute:" << attr_name << " does not have the same dimensions as the data" << endl;
          cout << "Exiting!" << endl;
          exit(-1);
        }
    } 

  //read
  status = H5Aread(attr, datatype_in,data->get_arr());
  if (status < 0)
    {
      cout << "Reading atribute: " << attr_name << " from object: " << dataset_name << " failed !" << endl;
      cout << "Exiting!" << endl;
      exit(-1);
    }
  
}

/*
  Write data to chosen dataset - This routine assumes stride = 0, just write connected data chunks to set
 */
template<typename T> void h5w::write_data(const char *dataset_path, unsigned long long int *offset, arr_1D<T>* data)
{
/*
  dataset_path - absolute path of the dataset to which to write to
  offset       - where in the dataset to write the data
  data         - container class holding an n-dimensional array as a 1D array,
                 also, has the rank, dimension lengths and datatype of the array!               
*/

  // general
  herr_t status;
  int i;
 
  // variables of the set already in file
  hid_t dataspace, datatype, dataset;
  hsize_t *dim, *wr_count, *wr_offset;
  int rank;
  
  // dataspace and other variables for data to be written
  hid_t chunk_dataspace, chunk_datatype;
  hsize_t *chunk_count, *chunk_offset, *chunk_dim;
  
  chunk_datatype = H5Tcopy(get_h5_datatype(data->get_dtype()));  // hdf5 datatype of input data
  // dataset properties
  dataset = H5Dopen(file, dataset_path);
  datatype = H5Dget_type(dataset);
  dataspace = H5Dget_space(dataset);
  rank = H5Sget_simple_extent_ndims(dataspace);
  
  if (H5Tequal(datatype, chunk_datatype) <= 0)   // compare datatypes
    {
      cout << "write_data(): Wrong input datatype for writing to dataset: " << dataset_path << endl;
      cout << "Exiting!" << endl;
      exit(-1);
    }
  
  
  if (rank < data->get_rank())
    {
      cout << "write_data(): Rank of dataset smaller than rank of input data for dataset: " << dataset_path << endl;
      cout << "Exiting!" << endl;
      exit(-1);
    }
  dim = new hsize_t[rank];
  status = H5Sget_simple_extent_dims(dataspace, dim, NULL);
 
  for (i=0; i<rank; i++)
    if ( (int) dim[i] < data->get_dims()[i])
      {
        cout << "write_data(): Size of input data chunk larger than the size of dataset along one or more dimensions" << endl;
        cout << i << " " << dim[i] << " " << data->get_dims()[i] << endl;
        cout << "Exiting!" << endl;
        exit(-1);
      }
  
  // where in the dataset to write
  wr_count = new hsize_t[rank];
  wr_offset = new hsize_t[rank];

  // write all data from input chunk - count = input_dim !!
  chunk_count = new hsize_t[data->get_rank()];
  chunk_offset = new hsize_t[data->get_rank()];
  chunk_dim = new hsize_t[data->get_rank()];
  

  for (i=0; i<rank; i++)
    {
      wr_count[i] = 1;
      wr_offset[i] = offset[i];    // where in the dataset to write the data chunk
    }
  for (i=0; i< data->get_rank(); i++)
    {
      wr_count[i] = data->get_dims()[i];
      chunk_count[i]=data->get_dims()[i];
      chunk_offset[i] = 0;         // write the whole chunk
      chunk_dim[i] = data->get_dims()[i];
    }   
        
  // hyperslab in the  file
  status = H5Sselect_hyperslab(dataspace, H5S_SELECT_SET, wr_offset, NULL, wr_count, NULL);
  
  // hyperslab of the input data
  chunk_dataspace = H5Screate_simple(data->get_rank(),chunk_dim,NULL);
  status = H5Sselect_hyperslab(chunk_dataspace,  H5S_SELECT_SET, chunk_offset, NULL, chunk_count, NULL);

  // writing
  status = H5Dwrite(dataset, chunk_datatype, chunk_dataspace, dataspace, H5P_DEFAULT, data->get_arr());


  if(dim){
    delete [] dim;
    dim = NULL;
  }

  if(wr_count){
    delete [] wr_count;
    wr_count = NULL;
  }

  if(wr_offset){
    delete [] wr_offset;
    wr_offset = NULL;
  }

  if(chunk_count){
    delete [] chunk_count;
    chunk_count = NULL;
  }  

  if(chunk_offset){
    delete [] chunk_offset;
    chunk_offset = NULL;
  }
 
  if(chunk_dim){
    delete [] chunk_dim;
    chunk_dim = NULL;
  }

  // close everything!
  H5Tclose(chunk_datatype);
  H5Tclose(datatype);
  H5Sclose(chunk_dataspace);
  H5Sclose(dataspace);
  H5Dclose(dataset);

}

/*
  Read data chunk from a chosen dataset
 */
template<typename T> void h5w::read_data(const char *dataset_path, int *offset, arr_1D<T>* data)
{
/*
  dataset_path - absolute path of the dataset to which to write to
  offset       - where in the dataset to write the data
  data         - container class holding an n-dimensional array as a 1D array,
                 also, has the rank, dimension lengths and datatype of the array!               
*/

  // general
  herr_t status;
  int i;
 
  // variables of the set already in file
  hid_t dataspace, datatype, dataset;
  hsize_t *dim, *wr_count, *wr_offset;
  int rank;
  
  // dataspace and other variables for data to be written
  hid_t chunk_dataspace, chunk_datatype;
  hsize_t *chunk_count, *chunk_offset, *chunk_dim;
  
  
  chunk_datatype = H5Tcopy(get_h5_datatype(data->get_dtype()));  // hdf5 datatype of input data
  // dataset properties
  dataset = H5Dopen(file, dataset_path);
  datatype = H5Dget_type(dataset);
  dataspace = H5Dget_space(dataset);
  rank = H5Sget_simple_extent_ndims(dataspace);
  
  if (H5Tequal(datatype, chunk_datatype) <= 0)   // compare datatypes
    {
      cout << "Wrong input datatype for writing to dataset: " << dataset_path << endl;
      cout << "Exiting!" << endl;
      exit(-1);
    }
  
  
  if (rank < data->get_rank())
    {
      cout << "read_data(): Rank of dataset smaller than rank of input data for dataset: " << dataset_path << endl;
      cout << "Exiting!" << endl;
      exit(-1);
    }
  dim = new hsize_t[rank];
  status = H5Sget_simple_extent_dims(dataspace, dim, NULL);
 
  for (i=0; i<rank; i++)
    if ( (int) dim[i] < data->get_dims()[i])
      {
        cout << "read_data(): Size of input data chunk larger than the size of dataset along one or more dimensions" << endl;
        cout << i << " " << dim[i] << " " << data->get_dims()[i] << endl;
        cout << "Exiting!" << endl;
        exit(-1);
      }
  
  // where in the dataset to write
  wr_count = new hsize_t[rank];
  wr_offset = new hsize_t[rank];

  // write all data from input chunk - count = input_dim !!
  chunk_count = new hsize_t[data->get_rank()];
  chunk_offset = new hsize_t[data->get_rank()];
  chunk_dim = new hsize_t[data->get_rank()];
  

  for (i=0; i<rank; i++)
    {
      wr_count[i] = 1;
      wr_offset[i] = offset[i];    // where in the dataset to write the data chunk
    }
  for (i=0; i< data->get_rank(); i++)
    {
      wr_count[i] = data->get_dims()[i];
      chunk_count[i]=data->get_dims()[i];
      chunk_offset[i] = 0;         // write the whole chunk
      chunk_dim[i] = data->get_dims()[i];
    }   
        
  // hyperslab in the  file
  status = H5Sselect_hyperslab(dataspace, H5S_SELECT_SET, wr_offset, NULL, wr_count, NULL);
  
  // hyperslab of the input data
  chunk_dataspace = H5Screate_simple(data->get_rank(),chunk_dim,NULL);
  status = H5Sselect_hyperslab(chunk_dataspace,  H5S_SELECT_SET, chunk_offset, NULL, chunk_count, NULL);

  // writing
  status = H5Dread(dataset, chunk_datatype, chunk_dataspace,dataspace,  H5P_DEFAULT, data->get_arr());

  // close everything!
  H5Tclose(chunk_datatype);
  H5Tclose(datatype);
  H5Sclose(chunk_dataspace);
  H5Sclose(dataspace);
  H5Dclose(dataset);

  if(dim){
    delete [] dim;
    dim = NULL;
  }
  if(wr_count){
    delete [] wr_count;
    wr_count = NULL;
  }
  if(wr_offset){
    delete [] wr_offset;
    wr_offset = NULL;
  }
  if(chunk_count){
    delete [] chunk_count;
    chunk_count = NULL;
  }
  if(chunk_offset){
    delete [] chunk_offset;
    chunk_offset = NULL;
  }
  if(chunk_dim){
    delete [] chunk_dim;
    chunk_dim = NULL;
  }

}


#endif
