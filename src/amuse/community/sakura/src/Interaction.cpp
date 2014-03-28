#include "Interaction.h"

Interaction::Interaction() {
  ;
}

void Interaction::calc_a(vector<Particle> &particle, Communicator &communicator) {
  int N = particle.size();
  for(int i=0; i<N; i++) {
    particle[i].set_ax(0.0);	  
    particle[i].set_ay(0.0);	  
    particle[i].set_az(0.0);	  
  }

  vector<double> myacc(3*N, 0);

  for(int i=communicator.get_my_begin(); i<communicator.get_my_end(); i++) {
    for(int j=0; j<N; j++) {
      if(i != j) {
	double rjix = particle[j].get_x()-particle[i].get_x();
	double rjiy = particle[j].get_y()-particle[i].get_y();
	double rjiz = particle[j].get_z()-particle[i].get_z();	      
	double dr2 = rjix*rjix + rjiy*rjiy + rjiz*rjiz;
	double dr = sqrt(dr2);
	double dr3 = dr*dr2;
	
	double dax = particle[j].get_mass() / dr3 * rjix;
	double day = particle[j].get_mass() / dr3 * rjiy;
	double daz = particle[j].get_mass() / dr3 * rjiz;	      
	      
	myacc[i*3+0] += dax;
	myacc[i*3+1] += day;
	myacc[i*3+2] += daz;
      }	      
    }
  }

  vector<double> acc(3*N, 0);
  communicator.join(myacc, acc);
  communicator.bcast(acc);
  
  for(int i=0; i<N; i++) {
    particle[i].set_ax(acc[i*3+0]);
    particle[i].set_ay(acc[i*3+1]);
    particle[i].set_az(acc[i*3+2]);
  }  
}






