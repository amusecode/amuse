#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>
#include <cstdio>

#include <types.h>
#include <functions.h>
#include <my_errors.h>
#include <utilis.h>
#include <mpi_types.h>

#include "sys/time.h"

#include <stdarg.h>

using namespace std;

HostError CPU_memcheck(const char *file, const int line, string path){
	
	string str1, str2, str3;
   int memtot, memfree, membuf, cached;
	int local_rank;
	MPISafeCall(MPI_Comm_rank(MPI_COMM_WORLD, &local_rank));

	if(local_rank == 0){

		ofstream meminfo;
		char *output_name;
		string temp;
		temp = path + "cpu_memory.dat";
		output_name = to_char(temp);
		meminfo.open(output_name);
		if(!meminfo)
         return HNoFile;

	   ifstream proc("/proc/meminfo");
		if(!proc)
			return HNoFile;

   	proc>>str1>>str2>>str3;
		if(str1!="MemTotal:")
			return HInvalidMeminfo;
		isNumeric_int(to_char(str2),&memtot);

   	proc>>str1>>str2>>str3;
		if(str1!="MemFree:")
			return HInvalidMeminfo;
   	isNumeric_int(to_char(str2),&memfree);

   	proc>>str1>>str2>>str3;
		if(str1!="Buffers:")
			return HInvalidMeminfo;
   	isNumeric_int(to_char(str2),&membuf);

   	proc>>str1>>str2>>str3;
		if(str1!="Cached:")
			return HInvalidMeminfo;
   	isNumeric_int(to_char(str2),&cached);


		meminfo<<fixed<<setprecision(1);
		float used = memtot/1024.0 - memfree/1024.0 - membuf/1024.0 - cached/1024.0;
		float total = memtot/1024.0; 
   	//cout<<" Buffered CPU mem : "<< membuf/1024.0<<" MB "<<endl;
   	//cout<<" Cached CPU mem   : "<< cached/1024.0<<" MB \n"<<endl;
   	meminfo<<" Used CPU memory @ file/line "<<file<<"/"<<line<<" : "<<used/total*100.0  <<" % \n"<<endl;
   	//cout<<" Free  CPU memory : "<< memfree/1024.0 + membuf/1024.0 + cached/1024.0 <<" MB "<<endl;
   	//cout<<" Total CPU memory : "<< memtot/1024.0 <<" MB "<<endl;

		if(total-used < 20.0)
			return HNoMemory;

		meminfo.close();
	}

	return HNoError;

}

HostError open_file(string argument, string path){

      ofstream generic;
      string temp;
      char *output_name;

      temp = path + "gpu_memory.dat";
		output_name = to_char(temp);
		generic.open(output_name, ios::out);
      generic.close();

      temp = path + "cpu_memory.dat";
      output_name = to_char(temp);
      generic.open(output_name, ios::out);
      generic.close();

      temp = path + "HiGPUslog.dat";
      output_name = to_char(temp);
      generic.open(output_name, ios::out);
      generic<<" Launch string used   : "<< argument <<endl;
      generic.close();


      temp = path + "energy.dat";
      output_name = to_char(temp);
      generic.open(output_name, ios::out);
      generic.close();

#ifdef CHECK_TIMES

      temp = path + "active_particles.dat";
      output_name = to_char(temp);
      generic.open(output_name, ios::out);
      generic.close();
#endif


   return HNoError;
}

HostError append_file(string argument, string path){

      ofstream generic;
    
      string temp;
      char *output_name;

      temp = path + "gpu_memory.dat";
      output_name = to_char(temp);
      generic.open(output_name, ios::app);
		generic<<" RESTART AT THIS POINT ***************************************************************"<<endl;
      generic.close();

      temp = path + "cpu_memory.dat";
      output_name = to_char(temp);
      generic.open(output_name, ios::app);
		generic<<" RESTART AT THIS POINT ***************************************************************"<<endl;
      generic.close();

      temp = path + "HiGPUslog.dat";
      output_name = to_char(temp);
      generic.open(output_name, ios::app);
		generic<<" RESTART AT THIS POINT ***************************************************************"<<endl;
      generic<<" Launch string used   : "<< argument <<endl;
		generic.close();


      temp = path + "energy.dat";
      output_name = to_char(temp);
      generic.open(output_name, ios::app);
      generic<<" \n \n"<<endl;
      generic.close();

#ifdef CHECK_TIMES

   temp = path + "active_particles.dat";
   output_name = to_char(temp);
	generic.open(output_name, ios::app);
   generic.close();

#endif

  return HNoError;

}

HostError print_info(double plc __attribute__((unused)), double plm __attribute__((unused)), double rs __attribute__((unused)), double ms __attribute__((unused)), string path){


#ifdef CHECK_ERRORS
   cout<<" Check errors option ENABLED "<<endl;
#else
   cout<<" Check errors option DISABLED "<<endl;
#endif

#ifdef CHECK_TIMES
   ofstream generic;
	string temp;
	char *output_name;
   temp = path + "times.dat";
   output_name = to_char(temp);
   generic.open(output_name, ios::out);
   generic.close();
   cout<<" Check times option ENABLED "<<endl;
#else
   cout<<" Check times option DISABLED "<<endl;
#endif

#ifdef PLUMMER
      cout<<" Plummer Potential option ENABLED "<<endl;
      ofstream hlog;
		string temp;
		char *output_name;
      temp = path + "HiGPUslog.dat";
      output_name = to_char(temp);
      hlog.open(output_name, ios::app);
      hlog<<" Plummer (core) : "<<plc<<endl;
      hlog<<" Plummer (mass) : "<<plm<<endl;
      hlog<<" =============================================== "<<endl;
      hlog.close();
#else
      cout<<" Plummer Potential option DISABLED "<<endl;
#endif

#ifdef GALAXY
      cout<<" Milky Way Potential option ENABLED "<<endl;
      ofstream hlog;
		temp = path + "HiGPUslog.dat";
		output_name = to_char(temp);
		hlog.open(output_name, ios::app);
      hlog<<" Scale radius : "<<rs<<endl;
      hlog<<" Scale mass   : "<<ms<<endl;
      hlog<<" =============================================== "<<endl;
      hlog.close();
#else
      cout<<" Milky Way Potential option DISABLED "<<endl;
#endif



   return HNoError;

}

HostError adjust_param_ifwarmstart(bool *CDM, bool *CDV, unsigned *FMAX, string warm_start_file, string path, string *FINP, double *GTIME_WARM, double DTPRINT){

   *FINP = path + warm_start_file;
	*CDM = 0;
   *CDV = 0;
   unsigned FMAX_old = *FMAX;
   string loc = to_string(*FMAX);
   loc = (warm_start_file).substr(0, loc.length());
   *FMAX = to_uint(loc);
   cout<<" FMAX NEW : "<<*FMAX<<endl;
   *GTIME_WARM = (*FMAX-FMAX_old)*DTPRINT;
   cout<<" GTIME WARM : "<<*GTIME_WARM<<endl;

   return HNoError;

}


HostError check_argv(int ac, char *av[], string *param, bool *warm_start, string *warm_start_file, bool *vir, double *ratio, bool *setdev, vector<unsigned int>& dev, double *plummer_core, double *plummer_mass, double *mscale, double *rscale, string *argument, bool *cleanstop){

   bool pl = 0;
	bool gl = 0;

   for(int i = 0; i < ac; i++){
      if(av[i][0]=='-'){
         if(to_string(&av[i][1])=="h" || to_string(&av[i][1])=="help"){
            cout<<endl;
            cout<<" ------------------------------------------------------- "<<endl;
            cout<<" ----------------------Help of HiGPUs------------------- "<<endl;
            cout<<" ------------------------------------------------------- "<<endl;
            cout<<" Usage : mpirun -np [#nodes] ./HiGPUs.x [options] "<<endl;
            cout<<" Options : "<<endl;
            cout<<"           -f [file]  (input_param.txt) = it specifies the file containing the simulation parameters"<<endl;
            cout<<"           -h = it shows this help screen "<<endl;
            cout<<"           -r [file] = restart the simulation from the specified file (HiGPUslog.dat is necessary !!)"<<endl;
            cout<<"           -v [value] = rescale to a specific ""value"" of virial ratio"<<endl;
				cout<<"           -d [id]:[id]: ... = enable manual selection of the available gpu devices "<<endl;
            cout<<"           -p b=[value]:M=[value] = if you compiled with -DPLUMMER, you need to specify the parameters b and M of the Plummer potential"<<endl;
				cout<<"           -gal R=[value]:M=[value] =  if you compiled with -DGALAXY, you need to specify scale radius R (pc) and scale mass M (M_sun)"<<endl;
				cout<<" ------------------------------------------------------- "<<endl;
				cout<<" ------------------------------------------------------- "<<endl;
				cout<<" ------------------------------------------------------- "<<endl;
            cout<<endl;

				*cleanstop = 1;
				return HNoError;

         }

         else if(to_string(&av[i][1])=="r" || to_string(&av[i][1])=="restart"){
            *warm_start = 1;
            *warm_start_file = to_string(&av[i+1][0]);
         }

         else if(to_string(&av[i][1])=="f"){
            *param = to_string(&av[i+1][0]);
            if(param -> size() == 0) return HInvalidArgv;
         }

         else if(to_string(&av[i][1])=="v"){
            *vir = 1;
            *ratio = to_double(to_string(&av[i+1][0]));
            if(*ratio < 0.0) return HInvalidArgv;
         }

         else if(to_string(&av[i][1])=="d"){
            *setdev = 1;
            string dev_param = to_string(&av[i+1][0]);
            unsigned int size = dev_param.size();
            for( unsigned int t = 0; t < size ; t++){
                if(dev_param[t] != ':') dev.push_back(to_uint(&dev_param[t]));
            }

            if(dev.size() == 0) return HInvalidArgv;
         }

         else if(to_string(&av[i][1])=="plummer" || to_string(&av[i][1])=="p"){

            pl = 1;
            string plum_param = to_string(&av[i+1][0]);
            size_t found1 = plum_param.find('=',0);
            size_t found2 = plum_param.find(':',found1+1);

            if (found1==string::npos || found2==string::npos)
                     return HInvalidArgv;

            string b_plum = plum_param.substr(int(found1+1), int(found2)-int(found1+1));
            *plummer_core = to_double(b_plum);

            if(*plummer_core <= 0.0)
               return HInvalidArgv;

            found1 = plum_param.find('=',found2+1);

            if (found1==string::npos)
                     return HInvalidArgv;

            b_plum = plum_param.substr(int(found1+1), plum_param.length()-int(found1+1));
            *plummer_mass = to_double(b_plum);

            if(*plummer_mass <= 0.0)
               return HInvalidArgv;
         }

			else if(to_string(&av[i][1])=="galaxy" || to_string(&av[i][1])=="gal"){

            gl = 1;
            string gala_param = to_string(&av[i+1][0]);
            size_t found1 = gala_param.find('=',0);
            size_t found2 = gala_param.find(':',found1+1);

            if (found1==string::npos || found2==string::npos)
                     return HInvalidArgv;

            string r_gala = gala_param.substr(int(found1+1), int(found2)-int(found1+1));
            *rscale = to_double(r_gala);

            if(*rscale <= 0.0)
               return HInvalidArgv;

            found1 = gala_param.find('=',found2+1);

            if (found1==string::npos)
                     return HInvalidArgv;

            r_gala = gala_param.substr(int(found1+1), gala_param.length()-int(found1+1));
            *mscale = to_double(r_gala);

            if(*mscale <= 0.0)
               return HInvalidArgv;
         }


         else
            return HInvalidArgv;

      }

   }

   if((*param).length() == 0){
      cout<<" Assuming input file for parameters input_param.txt"<<endl;
     *param = "input_param.txt";
   }

   if((*warm_start_file).length() == 0 && *warm_start == 1)
      return HInvalidArgv;



#ifdef PLUMMER
   if(!pl)
      return HNoDefPlummer2;
#else
   if(pl)
      return HNoDefPlummer;
#endif

#ifdef GALAXY
   if(!gl)
      return HNoDefGalaxy2;
#else
   if(gl)
      return HNoDefGalaxy;
#endif


   for(int i=0; i < ac; i++)
      *argument += to_string(&av[i][0]) + " ";


   return(HNoError);
}

HostError isDivisible(unsigned int *N, unsigned int *M, int size, unsigned int NGPU, unsigned int TPB, unsigned int *BFMAX){

   int product = size*NGPU*TPB;

	double esp = ceil(log(*N) / log (2.0));

   *M = *N;
   *N = (unsigned int) pow(2.0,esp);

   if(*N < product) *N = product;


   if(*N%product != 0){
      cout<<" Minimum N allowed per computational node : "<< product<<endl;
      cout<<" N must be an integer multiple of : "<<product<<endl;
      return HInvalidN;
   }

   *BFMAX = *N/product;

   return HNoError;

}

HostError Max_dimension(unsigned int TPB, unsigned int BFMAX, unsigned int N, unsigned int *dim_max, unsigned int GPUMINTHREADS){

   *dim_max = N;

   for(unsigned int nextsize = 1; nextsize < N; nextsize++){

      int BL = ceil((double)nextsize/TPB);
      int dim = TPB*BL;
      unsigned int threads = TPB;
      unsigned int bfmax = BFMAX;

      unsigned int Bfactor = 1;

      while(threads*BL < GPUMINTHREADS && threads > 32){
       threads /= 2;
       bfmax *= 2;
       BL = ceil((double)nextsize/threads);
      }
      dim = threads*BL;

      while(threads*BL < GPUMINTHREADS && Bfactor < bfmax){
       BL *= 2;
       Bfactor *= 2;
      }

    unsigned int sz = dim*Bfactor;
      if(sz > *dim_max)
         *dim_max = sz;

      }

   int local_rank;
   MPI_Comm_rank(MPI_COMM_WORLD, &local_rank);

   if(local_rank == 0)
      cout<<" Maximum dimension chosen for the accelerations array : "<<*dim_max<<endl;

   MPI_Barrier(MPI_COMM_WORLD);

   return HNoError;

}

HostError cpu_read_external(const string file_name, double4 *pos, float4 *vel, double4 *veldb, const unsigned int N, const unsigned int M, const bool CDM, const bool CDV, const bool warm){

   ifstream fin_data;

   char* name = to_char(file_name);
   fin_data.open(name);

   if(warm)
      cout<<" WARM START FROM FILE : "<<file_name<<endl;

   if(!fin_data)
      return HNoFile;

   unsigned int lines = countlines(fin_data);

   if(M > lines)
      return HNoLines;

   for(unsigned int i = 0; i < M; i++)
      fin_data >> pos[i].x >> pos[i].y >> pos[i].z  >> veldb[i].x >> veldb[i].y >> veldb[i].z >> pos[i].w;

   fin_data.close();

   for(unsigned int i = 0; i < M; i++)  veldb[i].w = pos[i].w;

   cout<<"\n Read initial conditions file : "<<name<<"\n";

   Utilis use;

   if(CDM && !warm){
      vector <double> cm = use.CdM(pos,M);
      cout<<" Previous center of mass coordinates : "<<cm[0]<<"  "<<cm[1]<<"  "<<cm[2]<<"\n";
   }
   else
      cout<<" WARNING : Center of mass different from the system of reference origin \n";

   if(CDV && !warm){
      vector <double> cv = use.CdM(veldb,M);
      cout<<" Previous center of velocity coordinates : "<<cv[0]<<"  "<<cv[1]<<"  "<<cv[2]<<"\n";
   }
else
      cout<<" WARNING : Center of velocity different from the system of reference origin \n";

   for(unsigned int i = 0; i < M; i++){
      vel[i].x = veldb[i].x;
      vel[i].y = veldb[i].y;
      vel[i].z = veldb[i].z;
      vel[i].w = veldb[i].w;
	}

    for(unsigned int i = M; i < N; i++){
      pos[i].x = rand()/RAND_MAX*100000.;
      pos[i].y = rand()/RAND_MAX*100000.;
      pos[i].z = rand()/RAND_MAX*100000.;
      pos[i].w = 0.0;
      veldb[i].x = rand()/RAND_MAX*100000.;
      veldb[i].y = rand()/RAND_MAX*100000.;
      veldb[i].z = rand()/RAND_MAX*100000.;
      veldb[i].w = 0.0;
      vel[i].x = rand()/RAND_MAX*100000.;
      vel[i].y = rand()/RAND_MAX*100000.;
      vel[i].z = rand()/RAND_MAX*100000.;
      vel[i].w = 0.0;
   }



   return HNoError;

}


HostError cpu_read_params(

      const string file_to_read,
      unsigned int *N,
      unsigned int *gpus,
      unsigned int *threads,
      double       *totaltime,
      double       *dtmax,
      double       *dtmin,
      double       *epsilon,
      double       *eta6,
      double       *eta4,
      double       *dtprint,
      unsigned int *filemax,
      bool         *cdm,
      bool         *cdv,
      string       *file,
      string       *gpu_name,
		string       path) {

      ifstream fin_data;
      string temp;
		char *file_name, *output_name;

     file_name = to_char(file_to_read);
     fin_data.open(file_name);

     if(!fin_data)
       return HNoFile;

      ofstream hlog;
		temp = path + "HiGPUslog.dat";
      output_name = to_char(temp);
      hlog.open(output_name, ios::app);

     fin_data >> *N;
     fin_data.ignore(300,'!');

     fin_data >> *gpus;
     fin_data.ignore(300,'!');

     fin_data >> *threads;
     fin_data.ignore(300,'!');

     fin_data >> *totaltime;
     fin_data.ignore(300,'!');

     fin_data >> *dtmax;
     fin_data.ignore(300,'!');

     fin_data >> *dtmin;
     fin_data.ignore(300,'!');

     fin_data >> *epsilon;
     fin_data.ignore(300,'!');

     fin_data >> *eta6;
     fin_data.ignore(300,'!');

     fin_data >> *eta4;
     fin_data.ignore(300,'!');

     fin_data >> *dtprint;
     fin_data.ignore(300,'!');

     fin_data >> *filemax;
     fin_data.ignore(300,'!');

     fin_data >> *cdm;
     fin_data.ignore(300,'!');

     fin_data >> *cdv;
     fin_data.ignore(300,'!');

     fin_data >> *file;
     fin_data.ignore(300,'!');

     getline(fin_data, *gpu_name);
     getline(fin_data, *gpu_name);

     size_t found = gpu_name->find('#');
     if (found==string::npos)
        return HNotFound;

     gpu_name->resize(int(found));

     *dtmax = pow(2.0, *dtmax);
     *dtmin = pow(2.0, *dtmin);


     hlog<<" ================================================== "<<endl;
     hlog<<" Read parameters file : "<<file_to_read<<endl;
     hlog<<" N                   : "<<*N<<endl;
     hlog<<" Gpus                : "<<*gpus<<endl;
     hlog<<" Threads per block   : "<<*threads<<endl;
     hlog<<" Time of integration : "<<*totaltime<<endl;
     hlog<<" Max time step       : "<<*dtmax<<endl;
     hlog<<" Min time step       : "<<*dtmin<<endl;
     hlog<<" Softening           : "<<*epsilon<<endl;
     hlog<<" eta 6th order       : "<<*eta6<<endl;
     hlog<<" eta 4th order       : "<<*eta4<<endl;
     hlog<<" time for printing   : "<<*dtprint<<endl;
     hlog<<" Max output files    : "<<*filemax<<endl;
     hlog<<" CDM scale           : "<<*cdm<<endl;
     hlog<<" CDV scale           : "<<*cdv<<endl;
     hlog<<" Input file          : "<<*file<<endl;
     hlog<<" Process id          : "<<(int) getpid()<<endl;
     hlog<<" ================================================== "<<endl;

     fin_data.close();
     hlog.close();

     return (HNoError);

}


HostError __MPIbcast_otherparams(bool *warm_start, bool *VIR, bool *setdev, double *plummer_core, double *plummer_mass, double *rscale, double *mscale, double *ratio, vector<unsigned> &dev){

   unsigned size_dev = dev.size();
   int local_rank;
   MPI_Comm_rank(MPI_COMM_WORLD, &local_rank);


   MPISafeCall(MPI_Barrier(MPI_COMM_WORLD));
   MPISafeCall(MPI_Bcast(warm_start, 1, MPI::BOOL, 0, MPI_COMM_WORLD));
   MPISafeCall(MPI_Bcast(VIR, 1, MPI::BOOL, 0, MPI_COMM_WORLD));
   MPISafeCall(MPI_Bcast(setdev, 1, MPI::BOOL, 0, MPI_COMM_WORLD));
   MPISafeCall(MPI_Bcast(plummer_core, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD));
   MPISafeCall(MPI_Bcast(plummer_mass, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD));
	MPISafeCall(MPI_Bcast(rscale, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD));
   MPISafeCall(MPI_Bcast(mscale, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD));
   MPISafeCall(MPI_Bcast(ratio, 1, MPI::DOUBLE, 0, MPI_COMM_WORLD));
   MPISafeCall(MPI_Bcast(&size_dev, 1, MPI::INT, 0, MPI_COMM_WORLD));

   if(local_rank != 0){
      for(unsigned int i = 0; i < size_dev; i++) dev.push_back(0);
   }

   for(unsigned int i = 0; i < size_dev; i++){
      MPISafeCall(MPI_Bcast(&dev[i], 1, MPI::INT, 0, MPI_COMM_WORLD));
   }

   return HNoError;

}
HostError __MPIstart(int *rank, int *size, MPI_Datatype *mpi_float4, MPI_Datatype *mpi_double4){

   MPISafeCall(MPI_Comm_set_errhandler(MPI_COMM_WORLD, MPI_ERRORS_RETURN));

   MPISafeCall(MPI_Comm_rank(MPI_COMM_WORLD, rank));
   MPISafeCall(MPI_Comm_size(MPI_COMM_WORLD, size));

   MPISafeCall(MPI_Type_contiguous (4,MPI_FLOAT,mpi_float4));
   MPISafeCall(MPI_Type_commit(mpi_float4));
   MPISafeCall(MPI_Type_contiguous (4,MPI_DOUBLE, mpi_double4));
   MPISafeCall(MPI_Type_commit(mpi_double4));


   return HNoError;
}
HostError __MPIbcastParam(unsigned int *N, unsigned int *M, unsigned int *NGPU, unsigned int *TPB, double *TTIME,  double *DTMAX,
                           double *DTMIN, double *EPS, double *ETA6, double *ETA4, double *DTPRINT, unsigned int *BFMAX, string *gpu_name, int rank, int size){

   char *value = new char [256];
   if(rank==0)
      value = to_char(*gpu_name);
  
   MPISafeCall(MPI_Bcast( N,     1, MPI_INT,    0, MPI_COMM_WORLD));
   MPISafeCall(MPI_Bcast( M,     1, MPI_INT,    0, MPI_COMM_WORLD));
	MPISafeCall(MPI_Bcast( NGPU,  1, MPI_INT,    0, MPI_COMM_WORLD));
   MPISafeCall(MPI_Bcast( TPB,   1, MPI_INT,    0, MPI_COMM_WORLD));
   MPISafeCall(MPI_Bcast( BFMAX, 1, MPI_INT,    0, MPI_COMM_WORLD));
   MPISafeCall(MPI_Bcast( TTIME, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD));
   MPISafeCall(MPI_Bcast( DTMAX, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD));
   MPISafeCall(MPI_Bcast( DTMIN, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD));
   MPISafeCall(MPI_Bcast( DTPRINT, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD));
   MPISafeCall(MPI_Bcast( EPS,   1, MPI_DOUBLE, 0, MPI_COMM_WORLD));
   MPISafeCall(MPI_Bcast( ETA6,  1, MPI_DOUBLE, 0, MPI_COMM_WORLD));
   MPISafeCall(MPI_Bcast( ETA4,  1, MPI_DOUBLE, 0, MPI_COMM_WORLD));
   MPISafeCall(MPI_Bcast( value,  256, MPI_CHAR, 0, MPI_COMM_WORLD));

   *gpu_name = to_string(value);

   return random_check_Bcast(rank, size, "4i-6d", *N, *NGPU, *TPB, *BFMAX, *TTIME, *DTMAX, *DTMIN, *EPS, *ETA6, *ETA4, *DTPRINT);

}
HostError __MPIbcastPosVel(double4 *pos, float4 *vel, double4 *veldb, unsigned int N, int rank, int size){

   MPISafeCall(MPI_Bcast( pos, N, mpi_double4, 0, MPI_COMM_WORLD));
   MPISafeCall(MPI_Bcast( vel, N, mpi_float4, 0, MPI_COMM_WORLD));
   MPISafeCall(MPI_Bcast( veldb, N, mpi_double4, 0, MPI_COMM_WORLD));

   srand ( time(NULL) );
   int __rand = N*(double)rand()/RAND_MAX;

   return random_check_Bcast(rank, size, "2d", pos[__rand].y, veldb[__rand].z);

}

HostError random_check_Bcast( int rank, int size, const char* Format ... )
{

	if(size < 2)
		return HNoError;

	int to = (rank+1) % size;
	int from = (rank + size -1) % size;

	MPI_Status status;
      va_list Arguments;
      va_start(Arguments, Format);
		int number;
		int intArg0, intArg1;
		double doubleArg0, doubleArg1;

		if(Format[0]=='\0')
			return HInvalidFormat;

      for(int i = 0; Format[i] != '\0'; ++i )
      {
			if(Format[i]=='-')
				continue;

			isNumeric_int(&Format[i], &number);

            if (number > 0 )
            {
					i++;

					if(Format[i]=='\0')
						return HInvalidFormat;

					for(int j = 0; j < number; j++){

					if (Format[i] == 'i'){
                  intArg0=va_arg(Arguments, int);

						MPISafeCall(MPI_Sendrecv(&intArg0, 1, MPI_INT, to, 100, &intArg1, 1, MPI_INT, from, 100, MPI_COMM_WORLD, &status));

						if(intArg0 != intArg1)
							return HBcastFailed;
						}

					else if (Format[i] == 'd'){
						doubleArg0 = va_arg(Arguments, double);
						
						MPISafeCall(MPI_Sendrecv(&doubleArg0, 1, MPI_DOUBLE, to, 100, &doubleArg1, 1, MPI_DOUBLE, from, 100, MPI_COMM_WORLD, &status));

						if(doubleArg0 != doubleArg1)
							return HBcastFailed;
					}

					else
						return HInvalidFormat;
					}
				}

					else if (number == -1){
						if (Format[i] == 'i'){
							intArg0=va_arg(Arguments, int);

							MPISafeCall(MPI_Sendrecv(&intArg0, 1, MPI_INT, to, 100, &intArg1, 1, MPI_INT, from, 100, MPI_COMM_WORLD, &status));

							if(intArg0 != intArg1)
								return HBcastFailed;
						}

						else if (Format[i] == 'd'){
							doubleArg0 = va_arg(Arguments, double);

							MPISafeCall(MPI_Sendrecv(&doubleArg0, 1, MPI_DOUBLE, to, 100, &doubleArg1, 1, MPI_DOUBLE, from, 100, MPI_COMM_WORLD, &status));

							if(doubleArg0 != doubleArg1)
								return HBcastFailed;
						}
						else
							return HInvalidFormat;
					}

					else
						return HInvalidFormat;

				}

      va_end(Arguments);
		return HNoError;
}
