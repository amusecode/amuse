/*
  methods for HDF5 wrapper class
*/
#include <string>
#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include "h5w.h"

using namespace std;

// Constructor with the file name - if the file does not exist create it, else just open it...
h5w::h5w(char fname[])
{

  // do not call HDF5 error stream by default...
  H5Eset_auto(NULL, NULL);

  strcpy(filename, fname);
  // First, try to open a file
  cout << "Trying to open file: " << filename << endl;
  file = H5Fopen(filename,H5F_ACC_RDWR, H5P_DEFAULT);  

  // if not successful, create file...
  if (file < 0 )
    {
      cout << "Could not open file: " << filename << endl;
      cout << "Creating it now!" << endl << endl;

      file = H5Fcreate(filename,H5F_ACC_TRUNC, H5P_DEFAULT,H5P_DEFAULT); // this will overwrite the file ...
      
      if (file < 0)
        {
          cout << "Could not create file: " << filename << endl;
          cout << "Exiting!" << endl;
          exit(-1);
        }
    }
  cout << "File " << filename << " ready for use" << endl;

}

// Constructor with the file name and mode - will do what the mode says...
h5w::h5w(char fname[], const char mode)
{
  // do not call HDF5 error stream by default...
  H5Eset_auto(NULL, NULL);
  
  strcpy(filename, fname);
  // open
  if (mode == 'o')
    {
      //cerr << "Opening file: " << filename << endl;
      file = H5Fopen(filename,H5F_ACC_RDWR, H5P_DEFAULT);  
      if (file < 0)
        {
          cout << "Could not open file: " << filename << endl;
          cout << "Exiting!" << endl << endl;
          exit(-1);
        }
    }
  else if (mode =='n') 
    {
      //cout << "Creating file: " << filename << endl;
      file = H5Fcreate(filename,H5F_ACC_TRUNC, H5P_DEFAULT,H5P_DEFAULT); // this will overwrite the file ...
      if (file < 0)
        {
          cout << "Could not create file: " << filename << endl;
          cout << "Exiting!" << endl;
          exit(-1);
        }
    }
  else 
    {
      cout << "Mode for h5w object constructor not known!" << endl << "Exiting!" << endl;
      exit(-1);
    }
  
  //cout << "File " << filename << " ready for use" << endl;
}

void h5w::open(char fname[])
{

  // do not call HDF5 error stream by default...
  H5Eset_auto(NULL, NULL);

  strcpy(filename, fname);
  // First, try to open a file
  cout << "Trying to open file: " << filename << endl;
  file = H5Fopen(filename,H5F_ACC_RDWR, H5P_DEFAULT);  
  if (file < 0)
    {
      cout << "Could not open file: " << filename << endl;
      cout << "Exiting!" << endl << endl;
      exit(-1);
    }
  cout << "File " << filename << " ready for use" << endl;

}


void h5w::close(void)
{
  herr_t status;
  
  status = H5Fclose(file);
  if (status < 0)
    {
      cout << "Could not close the file: " << filename << endl << "Exiting!" << endl;
      exit(-1);
    }
}

// makes a dataset on atomic data
void h5w::make_dataset(const char *name, const char *type, const int rank, int *dims)
{
  /*
    name - full name of dataset
    type - datatype of dataset  (see get_h5_datatype for string names of HDF5 datatypes
    rank - number of dimensions of dataset
    dims - dimensions of dataset
  */

  hid_t datatype, dataset, dataspace;
  hsize_t *dim_set;
  herr_t status;

  int i;

  dim_set = new hsize_t [rank];

  for (i=0;i<rank;i++)
    dim_set[i] = dims[i];
  
  dataspace = H5Screate_simple(rank, dim_set, NULL);
  datatype = get_h5_datatype(type);
  status = H5Tset_order(datatype, H5T_ORDER_LE);      // order is little endian by default
  
  dataset = H5Dcreate(file, name, datatype, dataspace, H5P_DEFAULT);
  if (dataset < 0)
    {
      cout << "Dataset " << name << " could not be created in file " << filename <<endl;
      cout << "Exiting!" << endl;
      exit(-1);
    }
      
  if(dim_set){
    delete [] dim_set;
    dim_set = NULL;
  }


  H5Tclose(datatype);
  H5Dclose(dataset);
  H5Sclose(dataspace);

}
// makes a group with absolute pathname "name"  
void h5w::make_group(const char *name)
{
  hid_t grp;
  
  grp = H5Gcreate(file, name, 0);

  if (grp<0)
    {
      cout << "Could not create group: " << name << endl;
      cout << "Exiting!" << endl;
      exit(-1);
    }
}
// makes a group with absolute pathname "name"   - WITH SIZE_HINT
void h5w::make_group(const char *name, int size_hint = 0)
{
  hid_t grp;
  
  grp = H5Gcreate(file, name, size_hint);

  if (grp<0)
    {
      cout << "Could not create group: " << name << endl;
      cout << "Exiting!" << endl;
      exit(-1);
    }
}
    


// get the HF5 datatype for a given var name (passed as a string)
int get_h5_datatype(const char *name)
{
  if (strcmp(name,"int")==0)
    return H5T_NATIVE_INT;
  if (strcmp(name,"char")==0)
    // return H5T_NATIVE_CHAR;
    return H5T_C_S1;
  if (strcmp(name,"short char")==0)
    return H5T_NATIVE_SCHAR;
  if (strcmp(name,"unsigned char")==0)
    return H5T_NATIVE_UCHAR;
  if (strcmp(name,"short")==0)
    return H5T_NATIVE_SHORT;
  if (strcmp(name,"unsigned short")==0)
    return H5T_NATIVE_USHORT;
  if (strcmp(name,"unsigned int")==0)
    return H5T_NATIVE_UINT;
  if (strcmp(name,"unsigned long")==0)
    return H5T_NATIVE_ULONG;
  if (strcmp(name,"long long")==0)
    return H5T_NATIVE_LLONG;
  if (strcmp(name,"float")==0)
    return H5T_NATIVE_FLOAT;
  if (strcmp(name,"double")==0)
    return H5T_NATIVE_DOUBLE;
  if (strcmp(name,"long double")==0)   // No 128-bit doubles on Intel machines ?!? (at least none in H5Tpublic.h)
    return H5T_NATIVE_LDOUBLE;
  if (strcmp(name,"long")==0)
    return H5T_NATIVE_LONG;
  if (strcmp(name,"unsigned long long")==0)
    return H5T_NATIVE_ULLONG;

  cout << "Datatype string not recognized! Exiting..." << endl;
  exit(-1);
}


// write_attr overload to write single string by refenrece
void h5w::write_attr(const char *dataset_name, const char *attr_name, char *data)
{
  herr_t status;
  char type_tmp[20];
  hid_t group, dataset, object; // first 2 are temporary vars to see if the object is a group or dataset
  hid_t attr_dataspace, attr;
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
  if (strcmp(type_tmp,"char")==0)
    {
      int len;
      len = strlen(data);
      H5Tset_size(datatype_in, len);
    }
  
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
  
}


// write_attr overload to write single string by value
void h5w::write_attr(const char *dataset_name, const char *attr_name, const char *data)
{
  herr_t status;
  char type_tmp[20];
  hid_t group, dataset, object; // first 2 are temporary vars to see if the object is a group or dataset
  hid_t attr_dataspace, attr;
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
  if (strcmp(type_tmp,"char")==0)
    {
      int len;
      len = strlen(data);
      H5Tset_size(datatype_in, len);
    }
  
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
  
}

// read_attr overload to read single variable by refenrece
void h5w::read_attr(const char *dataset_name, const char *attr_name, char *data)
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

  hsize_t attr_len;
  attr_len = H5Tget_size(attr_datatype);  // size is in bytes!
  H5Tset_size(datatype_in, attr_len);

  if (H5Tequal(datatype_in, attr_datatype) <= 0)   // compare datatypes 
    {
      cout << "Type of variable passed not the same as type of data to be read in attribute: " << attr_name << endl;
      cout << "Exiting!" << endl;
      exit(-1);
    }
  
  char *tmp;
  tmp = new char[attr_len];

  //read 
  status = H5Aread(attr, attr_datatype, tmp);
  if (status < 0)
    {
      cout<< "Reading atribute: " << attr_name << " from object: " << dataset_name << " failed !" << endl << "Exiting!" << endl;
      exit(-1);
    }
  strcpy(data,tmp);

  H5Aclose(attr);
  H5Tclose(attr_datatype);
  H5Tclose(datatype_in);
}


// returns a string holding the name of the datatype of the dataset
string h5w::get_data_type(char *dataset_name)
{
  hid_t dset, dtype;
  string result;

  dset = H5Dopen(file, dataset_name);
  if (dset < 0)
    {
      cout << "Could not open dataset: " << dataset_name << " for obtaining its type!" << endl;
      exit(-1);
    }
  dtype = H5Dget_type(dset);

  if (H5Tequal(dtype,H5T_NATIVE_INT) > 0)
    return "int";
  if (H5Tequal(dtype,H5T_NATIVE_DOUBLE) > 0)
    return "double";
  if (H5Tequal(dtype,H5T_NATIVE_CHAR) > 0)
    return "char";
  if (H5Tequal(dtype,H5T_NATIVE_UCHAR) > 0)
    return "unsigned char";
  if (H5Tequal(dtype,H5T_NATIVE_SHORT) > 0)
    return "short";
  if (H5Tequal(dtype,H5T_NATIVE_USHORT) > 0)
    return "unsigned short";
  if (H5Tequal(dtype,H5T_NATIVE_UINT) > 0)
    return "unsigned int";
  if (H5Tequal(dtype,H5T_NATIVE_ULONG) > 0)
    return "unsigned long";
  if (H5Tequal(dtype,H5T_NATIVE_LLONG) > 0)
    return "long long";
  if (H5Tequal(dtype,H5T_NATIVE_FLOAT) > 0)
    return "float";
  if (H5Tequal(dtype,H5T_NATIVE_LONG) > 0)
    return "long";
  if (H5Tequal(dtype,H5T_NATIVE_ULLONG) > 0)
    return "unsigned long long";

  cout << "Datatype string of the dataset not recognized! Exiting..." << endl;
  exit(-1);
  
}


//                                      
// read one variable from the HDF5 file and return it as a double array
// (used for generalizing data input...)
//

void h5w::read_data_as_double(char *dataset_name, int *offset, arr_1D<double>* dbl_arr)
{
  string dtype;

  dtype = get_data_type(dataset_name);

  // now that the type is know call the appropriate wrapper function
  if (dtype == "double")
    {
      read_data(dataset_name, offset, dbl_arr);
      return;
    }
  if (dtype == "float")
    {
      read_float_as_double(dataset_name, offset, dbl_arr);
      return;
    }
  if (dtype == "int")
    {
      read_int_as_double(dataset_name, offset, dbl_arr);
      return;
    }
  if (dtype == "unsigned int")
    {
      read_uint_as_double(dataset_name, offset, dbl_arr);
      return;
    }
  if (dtype == "long")
    {
      read_long_as_double(dataset_name, offset, dbl_arr);
      return;
    }
  if (dtype == "long long")
    {
      read_llong_as_double(dataset_name, offset, dbl_arr);
      return;
    }
  if (dtype == "unsigned long")
    {
      read_ulong_as_double(dataset_name, offset, dbl_arr);
      return;
    }
  if (dtype == "unsigned long long")
    {
      read_ullong_as_double(dataset_name, offset, dbl_arr);
      return;
    }
  if (dtype == "char")
    {
      read_char_as_double(dataset_name, offset, dbl_arr);
      return;
    }
  if (dtype == "unsigned char")
    {
      read_uchar_as_double(dataset_name, offset, dbl_arr);
      return;
    }
  if (dtype == "short")
    {
      read_short_as_double(dataset_name, offset, dbl_arr);
      return;
    }
  if (dtype == "unsigned short")
    {
      read_ushort_as_double(dataset_name, offset, dbl_arr);
      return;
    }
  
  cout << "Datatype not known?" << endl;
  exit(0);


 
}
                         
/*
  Dummy routines to reading different types of data...
*/
void h5w::read_int_as_double(char *dataset_name, int *offset, arr_1D<double>* dbl_arr)
{
  arr_1D<int> data;
  data.reinit(dbl_arr->get_rank(), dbl_arr->get_dims());
  
  int total_len=1;
  read_data(dataset_name, offset, &data);
  for (int i=0; i <dbl_arr->get_rank(); i++)
    total_len *=dbl_arr->get_dims()[i];
  for (int i=0; i < total_len; i++)
    dbl_arr->get_arr()[i] = (double) data.get_arr()[i];

}


void h5w::read_uint_as_double(char *dataset_name, int *offset, arr_1D<double>* dbl_arr)
{
  arr_1D<unsigned int> data;
  data.reinit(dbl_arr->get_rank(), dbl_arr->get_dims());
  
  int total_len=1;
  read_data(dataset_name, offset, &data);
  for (int i=0; i <dbl_arr->get_rank(); i++)
    total_len *=dbl_arr->get_dims()[i];
  for (int i=0; i < total_len; i++)
    dbl_arr->get_arr()[i] = (double) data.get_arr()[i];
}

void h5w::read_long_as_double(char *dataset_name, int *offset, arr_1D<double>* dbl_arr)
{
  arr_1D<long> data;
  data.reinit(dbl_arr->get_rank(), dbl_arr->get_dims());
  
  int total_len=1;
  read_data(dataset_name, offset, &data);
  for (int i=0; i <dbl_arr->get_rank(); i++)
    total_len *=dbl_arr->get_dims()[i];
  for (int i=0; i < total_len; i++)
    dbl_arr->get_arr()[i] = (double) data.get_arr()[i];
}

void h5w::read_ulong_as_double(char *dataset_name, int *offset, arr_1D<double>* dbl_arr)
{
  arr_1D<unsigned long> data;
  data.reinit(dbl_arr->get_rank(), dbl_arr->get_dims());
  
  int total_len=1;
  read_data(dataset_name, offset, &data);
  for (int i=0; i <dbl_arr->get_rank(); i++)
    total_len *=dbl_arr->get_dims()[i];
  for (int i=0; i < total_len; i++)
    dbl_arr->get_arr()[i] = (double) data.get_arr()[i];
}
void h5w::read_ullong_as_double(char *dataset_name, int *offset, arr_1D<double>* dbl_arr)
{
  arr_1D<unsigned long long> data;
  data.reinit(dbl_arr->get_rank(), dbl_arr->get_dims());
  
  int total_len=1;
  read_data(dataset_name, offset, &data);
  for (int i=0; i <dbl_arr->get_rank(); i++)
    total_len *=dbl_arr->get_dims()[i];
  for (int i=0; i < total_len; i++)
    dbl_arr->get_arr()[i] = (double) data.get_arr()[i];
}

void h5w::read_float_as_double(char *dataset_name, int *offset, arr_1D<double>* dbl_arr)
{
  arr_1D<float> data;
  data.reinit(dbl_arr->get_rank(), dbl_arr->get_dims());
  
  int total_len=1;
  read_data(dataset_name, offset, &data);
  for (int i=0; i <dbl_arr->get_rank(); i++)
    total_len *=dbl_arr->get_dims()[i];
  for (int i=0; i < total_len; i++)
    dbl_arr->get_arr()[i] = (double) data.get_arr()[i];
}
  


void h5w::read_short_as_double(char *dataset_name, int *offset, arr_1D<double>* dbl_arr)
{
  arr_1D<short> data;
  data.reinit(dbl_arr->get_rank(), dbl_arr->get_dims());
  
  int total_len=1;
  read_data(dataset_name, offset, &data);
  for (int i=0; i <dbl_arr->get_rank(); i++)
    total_len *=dbl_arr->get_dims()[i];
  for (int i=0; i < total_len; i++)
    dbl_arr->get_arr()[i] = (double) data.get_arr()[i];
}
  

void h5w::read_ushort_as_double(char *dataset_name, int *offset, arr_1D<double>* dbl_arr)
{
  arr_1D<unsigned short> data;
  data.reinit(dbl_arr->get_rank(), dbl_arr->get_dims());
  
  int total_len=1;
  read_data(dataset_name, offset, &data);
  for (int i=0; i <dbl_arr->get_rank(); i++)
    total_len *=dbl_arr->get_dims()[i];
  for (int i=0; i < total_len; i++)
    dbl_arr->get_arr()[i] = (double) data.get_arr()[i];
}

void h5w::read_char_as_double(char *dataset_name, int *offset, arr_1D<double>* dbl_arr)
{
  arr_1D<char> data;
  data.reinit(dbl_arr->get_rank(), dbl_arr->get_dims());
  
  int total_len=1;
  read_data(dataset_name, offset, &data);
  for (int i=0; i <dbl_arr->get_rank(); i++)
    total_len *=dbl_arr->get_dims()[i];
  for (int i=0; i < total_len; i++)
    dbl_arr->get_arr()[i] = (double) data.get_arr()[i];
}
  

void h5w::read_uchar_as_double(char *dataset_name, int *offset, arr_1D<double>* dbl_arr)
{
  arr_1D<unsigned char> data;
  data.reinit(dbl_arr->get_rank(), dbl_arr->get_dims());
  
  int total_len=1;
  read_data(dataset_name, offset, &data);
  for (int i=0; i <dbl_arr->get_rank(); i++)
    total_len *=dbl_arr->get_dims()[i];
  for (int i=0; i < total_len; i++)
    dbl_arr->get_arr()[i] = (double) data.get_arr()[i];
}
  
void h5w::read_llong_as_double(char *dataset_name, int *offset, arr_1D<double>* dbl_arr)
{
  arr_1D<long long> data;
  data.reinit(dbl_arr->get_rank(), dbl_arr->get_dims());
  
  int total_len=1;
  read_data(dataset_name, offset, &data);
  for (int i=0; i <dbl_arr->get_rank(); i++)
    total_len *=dbl_arr->get_dims()[i];
  for (int i=0; i < total_len; i++)
    dbl_arr->get_arr()[i] = (double) data.get_arr()[i];
}

