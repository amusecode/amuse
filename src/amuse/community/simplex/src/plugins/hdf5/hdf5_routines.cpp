#include "h5w.h"

vector<int> source_indices;
void smear_out_sources(void);
/*
  This routine writes the simplex output in HDF5 format - DO NOT USE OPENMP UNLESS YOU MAKE CALLS TO HDF5 CRITICAL
*/
void write_hdf5_output(char *fname)
{
  int i,j,k;
  arr_1D<double> double_arr;
  arr_1D<unsigned long long> int_arr;
  int dims[2];
  int offset[2];

  char fullname[200];
  sprintf(fullname,"%s/%s",createOps.output_dir,fname);
  
  
  h5w file(fullname,'n');

  // write structure of the file
  file.make_group("/Header");
  file.make_group("/Vertices");
  dims[0] = createOps.numSites;
  dims[1] = 3;                              // only write 3D data for now...
  file.make_dataset("/Vertices/Coordinates","double",2,dims);
  file.write_attr("/Vertices/Coordinates","var_name","coords");
  dims[1] = 1;
  file.make_dataset("/Vertices/H_neutral_fraction","double",1,dims);
  file.write_attr("/Vertices/H_neutral_fraction","var_name","HI");
  file.make_dataset("/Vertices/number_density","double",1,dims);
  file.write_attr("/Vertices/number_density","var_name","n");
  file.make_dataset("/Vertices/Volume","double",1,dims);
  file.write_attr("/Vertices/Volume","var_name","volume");
  file.make_dataset("/Vertices/Luminosity","double",1,dims);
  file.write_attr("/Vertices/Luminosity","var_name","lum");
  
#ifdef OUTPUT_GADGET_ID
  file.make_dataset("/Vertices/Gadget_ID","unsigned long long",1,dims);
  file.write_attr("/Vertices/Gadget_ID","var_name","gadget_id");
#endif

  // write header

  file.write_attr("/Header","Number_of_Sites", createOps.numSites);
  file.write_attr("/Header","Number_of_sweeps", createOps.numSweeps);
  file.write_attr("/Header","BoxSizeInPc", createOps.sizeBox);
  file.write_attr("/Header","redshift", createOps.redShift);
  file.write_attr("/Header","Simulation_time", createOps.simTime);
  if (createOps.diffusePhotons)
    file.write_attr("/Header","difusePhotons_on","true");
  else
    file.write_attr("/Header","difusePhotons_on","false");
  if (createOps.recombination)
    file.write_attr("/Header","recombination_on","true");
  else
    file.write_attr("/Header","recombination_on","false");
#ifdef CONVERT_TO_PHYSICAL
  file.write_attr("/Header","converted_to_physical", "true");
#else
  file.write_attr("/Header","converted_to_physical", "false");
#endif
  

  // writing is done one variable at a time, through the arr_1D instance of size chunk_size
  // (its not equal to numSites in order to preserve memory....
  int chunk_size = 100000;  // MAKE THIS A PARAMETER
  if (chunk_size > createOps.numSites)
    chunk_size = createOps.numSites;
  
  // do the writing !!!
  offset[1] = 0;
  for (i=0; i <createOps.numSites; i+=chunk_size)
    {
      offset[0] = i;
      if (i+chunk_size >= createOps.numSites) // make sure not to write outside of data range
        dims[0] = createOps.numSites-i;
      else
        dims[0] = chunk_size;
      
      // writing coordinates
      dims[1] = 3; 
      
      double_arr.reinit(2,dims);
      for (j=0;j<dims[0];j++)
        for (k=0;k<3;k++)
          double_arr(j,k) = sites[i+j].x[k];
      file.write_data("/Vertices/Coordinates",offset, &double_arr);
      
      // writing neutral fraction
      
      double_arr.reinit(1,dims);
      for (j=0;j<dims[0];j++)
        double_arr(j) = sites[i+j].numAtoms / (sites[i+j].dens*sites[i+j].volume); // HI
      file.write_data("/Vertices/H_neutral_fraction",offset, &double_arr);
  
      // writing number density of gas
      
      double_arr.reinit(1,dims);
      for (j=0;j<dims[0];j++)
        double_arr(j) = sites[i+j].dens;
      file.write_data("/Vertices/number_density",offset, &double_arr);

      // writing cell volumes
      double_arr.reinit(1,dims);
      for (j=0;j<dims[0];j++)
        double_arr(j) = sites[i+j].volume;
      file.write_data("/Vertices/Volume",offset, &double_arr);

      // writing Luminositites
      double_arr.reinit(1,dims);
      for (j=0;j<dims[0];j++)
        double_arr(j) = sites[i+j].flux;
      file.write_data("/Vertices/Luminosity",offset, &double_arr);


      // writing Gadget ID of particles
#ifdef OUTPUT_GADGET_ID
      
      int_arr.reinit(1,dims);
      for (j=0;j<dims[0];j++)
        int_arr(j) = sites[i+j].gadget_ID;

      file.write_data("/Vertices/Gadget_ID",offset, &int_arr);
#endif
    }
  
  file.close();
  cout << "HDF5 file writen: " << fname << endl;
    
}

// Reads the IC from a gadget HDF5 file 
void read_vertex_list_from_gadget_HDF5()
{
  int sources=0;
  
  double unitLengthConv, unitMassConv;
  double UnitLength_in_cm, UnitMass_in_g, UnitLuminosity;
  double boxSize;
  double Omega0, OmegaBaryon, OmegaLambda;
  
  h5w file(createOps.nameInput,'o');        // open file for reading


  file.read_attr("/Header","HubbleParam", &createOps.hubbleParam); 
  cout << "\nHubble Parameter: h = " << createOps.hubbleParam << endl;

  file.read_attr("/Header","Omega0", &createOps.omegaMatter); 
  cout << "Omega0 = " << createOps.omegaMatter << endl;
  file.read_attr("/Header","OmegaLambda", &createOps.omegaLambda); 
  cout << "OmegaLambda = " << createOps.omegaLambda << endl;
  file.read_attr("/Header","OmegaBaryon", &createOps.omegaBaryon); 
  cout << "OmegaBaryon = " << createOps.omegaBaryon << endl;

  file.read_attr("/Header","Redshift", &createOps.redShift); 
  cout << endl << "Redshift: z = " << createOps.redShift << endl;
  
  unsigned int tmp;       // Stupid bug - if I read directly to numSites it "overwrites" perPadding for some reason...
  file.read_attr("/Header","NumPart_Total", &tmp);
  createOps.numSites = tmp;

  file.read_attr("/Header","UnitLength_in_cm", &UnitLength_in_cm);
  file.read_attr("/Header","UnitMass_in_g", &UnitMass_in_g);
  file.read_attr("/Header","UnitLuminosity", &UnitLuminosity);

  file.read_attr("/Header","BoxSize", &boxSize);
  
  boxSize = boxSize*UnitLength_in_cm/parsecToCm;  // box size read from
                                                  // gadget file in pc

  cout << "Box size: " << boxSize <<" pc" <<  endl;

  cout << endl;
  
  // conversion from comoving to physical:
  unitLengthConv = 1.e6;       // converts from Mpc to pc
  unitMassConv   = UnitMass_in_g/1.67262e-24;   // Converts from grams to atom numbers
  createOps.sizeBox = boxSize;
#ifdef CONVERT_TO_PHYSICAL
  unitLengthConv = unitLengthConv / (1.0+createOps.redShift) / createOps.hubbleParam;
  unitMassConv   *= (createOps.omegaBaryon/createOps.omegaMatter*0.73);   // Converts DM mass to HYDROGEN (x = 0.73) tracing the DM
  unitMassConv /= createOps.hubbleParam;
  createOps.sizeBox = createOps.sizeBox / (1.0+createOps.redShift) / createOps.hubbleParam;
  UnitLuminosity = UnitLuminosity / pow(createOps.hubbleParam, 2);
  
  cout << "Units converted from comoving to physical " << endl;
  cout << "New box size is: " << createOps.sizeBox << " pc\n " << endl;

#endif

  // Allocate memory
  allocSitesArray();
  cout << endl;

  // now read actual data
  arr_1D<float> float_arr;
  arr_1D<unsigned long long> int_arr;
  
  int dims[2], offset[2];
  int i,j,k,site, coord;

  int chunk_size = 1000;  // MAKE THIS A PARAMETER
  if (chunk_size > createOps.numSites)
    chunk_size = createOps.numSites;
  
  // do the reading !!!
  offset[1] = 0;
  for (i=0; i <createOps.numSites; i+=chunk_size)
    {
      offset[0] = i;
      if (i+chunk_size >= createOps.numSites) // make sure not to write outside of data range
        dims[0] = createOps.numSites-i;
      else
        dims[0] = chunk_size;
      
      // coordinates
      dims[1] = 3; 
      
      float_arr.reinit(2,dims);

      file.read_data("/PartType0/Coordinates",offset, &float_arr);
      for (site=0; site<dims[0]; site++)
        for (coord=0; coord<dims[1]; coord++)
          sites[site+i].x[coord] = float_arr(site,coord)*unitLengthConv/createOps.sizeBox; // coordinates in box size units

      // number of atoms (stored in density, to be calculated when volume of cells 
      // is obtained by the triangulation...)
      dims[1] = 0;
      float_arr.reinit(1,dims);

      file.read_data("/PartType0/Masses",offset, &float_arr);
      for (site=0; site<dims[0]; site++)
          sites[site+i].dens = float_arr(site) * unitMassConv; 
      // CONVERT DM TO GAS PROPERLY !!!!!
      
      // ionised fraction (1-neutral fraction)
      dims[1] = 0;
      float_arr.reinit(1,dims);

      file.read_data("/PartType0/NeutralFraction",offset, &float_arr);
      for (site=0; site<dims[0]; site++)
        {
          sites[site+i].ionised = 1.0 - float_arr(site); 
          if (sites[site+i].ionised < 0)
            sites[site+i].ionised = 0.0;
        }

      // luminosities
      dims[1] = 0;
      float_arr.reinit(1,dims);

      file.read_data("/PartType0/Luminosity",offset, &float_arr);
      for (site=0; site<dims[0]; site++)
        {
          if (float_arr(site) > 0)
            {
              sites[site+i].flux = float_arr(site)*UnitLuminosity;
              source_indices.push_back(site+i);

//               printf("source: x=%f, y=%f, z=%f -> f = %e\n", 
//                      sites[site+i].x[0],sites[site+i].x[1],sites[site+i].x[2],
//                      sites[site+i].flux);
              sources++;
            }
          else 
            sites[site+i].flux = 0.0;
        }
      
      // gadget ID
      dims[1] = 0;
      int_arr.reinit(1,dims);

      file.read_data("/PartType0/ParticleIDs",offset, &int_arr);
      for (site=0; site<dims[0]; site++)
        sites[site+i].gadget_ID = int_arr(site);

      // set borders
      for (site=0; site<dims[0]; site++)
        sites[site+i].border = 0;

    }
  
  file.close();
  cout << endl;
  cout << "Input HDF5 file read! Number of source vertices: " << sources << endl;
  cout << endl;

}


// During the reading of HDF5 input, the indices of vertices holding the sources will be stored in an array
// Now, use this array to find the sources and smear out their luminosities EQUALLY (should be more clever...) among
// their neighbours
void smear_out_sources(void)
{

  int i,j,index;
  double new_flux;

  cout << " -- Smearing out "<< source_indices.size() << " sources -- " << endl;
  for (i=0; i<source_indices.size(); i++)
    {
      index = source_indices[i];
      new_flux = sites[index].flux / ((double)sites[index].numNeigh+1.0); // all vertices get equal flux

      for (j=0; j<sites[index].numNeigh; j++)
        if (sites[sites[index].neighId[j]].border == 0)
          sites[sites[index].neighId[j]].flux = new_flux;
      sites[index].flux = new_flux;
    }


}

