#include <string>
#include <functions.h>
#include <sys/statvfs.h>
#include <iomanip>
#include <fstream>
#include <utilis.h>
#include <my_errors.h>

HostError CheckHDDMemory(bool *cleanstop __attribute__((unused)), string path){

   struct statvfs *buffer = new struct statvfs [1];
   statvfs(getenv("PWD"), buffer);
   double conv = (1024. * 1024.) * 1024.; //GB

   double total = (buffer->f_blocks*buffer->f_frsize) / conv;
   double sys_space = 0.05*total;
   double available = (buffer->f_bfree*buffer->f_frsize)/conv - sys_space;
   double used = total - sys_space - available;
   double perc = (used/(total-sys_space))*100.;

	ofstream out;
	string temp;
	char* output_name;
	temp = path + "HD_state.dat";
   output_name = to_char(temp);
   out.open(output_name, ios::app);
	out<<fixed<<setprecision(1);

	out<<" This file shows the state of your Hard Disk while HiGPUs is running. "<<endl;
	out<<" This file is updated whenever a snapshot is written. "<<endl;
	out<<" If the available space is less than 5 GB the program will abort its execution. "<<endl;
	out<<" --------------------------------------------------------------------------"<<endl;
	out<<" --------------------------------------------------------------------------"<<endl;
   out<<" Available  : "<<available<<" GB "<<endl;
   out<<" Total      : "<<total<<" GB "<<endl;
   out<<" Used       : "<<used<<" GB "<<endl;
   out<<" Percentage : "<<perc<<" % "<<endl;
	out<<" --------------------------------------------------------------------------"<<endl;
   out<<" --------------------------------------------------------------------------"<<endl;


   if(available < 5.0){ //less than 5 GB
      out<<" NO HD SPACE !!!! PLEASE FREE YOUR HD "<<endl;
      out<<" This program is going to terminate its execution now ! "<<endl;
		out.close();

      delete [] buffer;

#ifdef CHECK_ERRORS
      return HNoSpace;
#else
      *cleanstop = 1;
		return HNoError;
#endif
   }

	out.close();

   return HNoError;
}


HostError getEnergyLine_HiGPUslog(string *lastLine, string path)
{

    ifstream data;
	 string temp;
    char *input_name;
    temp = path + "HiGPUslog.dat";
    input_name = to_char(temp);
    data.open(input_name);

    if(!data)
       return HNoFile;

    char *name;

   getline(data, *lastLine);
   name = to_char(*lastLine);
      while(name[0]!='#'){
         getline(data, *lastLine);
         name = to_char(*lastLine);
      }

    data.close();

    return HNoError;
}



HostError AcquireEnergy(double *E, string path){

   string last_line;
   HostSafeCall(getEnergyLine_HiGPUslog(&last_line, path));

   size_t found  = last_line.find('#',1);
   size_t found2 = last_line.find('#',found+1);

   if (found==string::npos)
      return HNotFound;

   string energy = last_line.substr(int(found+1), int(found2)-int(found+1));
   *E = to_double(energy);

   return HNoError;

}
