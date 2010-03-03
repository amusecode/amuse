#include "worker_code.h"
#include "src/code.h"
#include <map>

class Particle{

public:
    int index;
    double x,y,z;
    int n0, n1, n2;
    
    Particle(int index, double x, double y, double z):index(index), x(x), y(y), z(z), n0(0), n1(0), n2(0) {
    }
    Particle(const Particle & original):index(original.index), x(original.x), y(original.y), z(original.z) {
    }
};

typedef std::map<int, Particle *> ParticlesMap;
typedef std::map<int, Particle *>::iterator ParticlesMapIterator;

int highest_index = 0;
ParticlesMap particlesMap;

int find_nearest_neighbors(){
  std::size_t n = particlesMap.size();
  
  double * x = new double[n];
  double * y = new double[n];
  double * z = new double[n];
  int * n0 = new int[n];
  int * n1 = new int[n];
  int * n2 = new int[n];
  Particle ** particles = new Particle*[n];
  ParticlesMapIterator i;
  int c = 0;
  
  for(i = particlesMap.begin(); i != particlesMap.end(); i++) {
    Particle * p = (*i).second;
    particles[c] = p;
    x[c] = p->x;
    y[c] = p->y;
    z[c] = p->z;
    
    c++;
  }
  
  int errorcode = find_nearest_neighbors(n, x, y, z, n0, n1, n2);
  if(errorcode) {
    return errorcode;
  }
  
  for(std::size_t j = 0 ; j < n; j++) {
    Particle * p = particles[j];
    
    p->n0 = n0[j] >= 0 ? particles[n0[j]]->index : -1;
    p->n1 = n1[j] >= 0 ? particles[n1[j]]->index : -1;
    p->n2 = n1[j] >= 0 ? particles[n2[j]]->index : -1;
  }
  
  return 0;
}

int get_nearest_neighbor(int index_of_the_particle, 
  double * index_of_the_neighbor, double * distance){
  
  if(index_of_the_particle > highest_index) {
    return -1;
  }
  
  Particle * p0 = particlesMap[index_of_the_particle];
  Particle * p1 = particlesMap[p0->n0];
  *index_of_the_neighbor  = p0->n0;
  *distance = distance_between_points(p0->x, p0->y, p0->z, p1->x, p1->y, p1->z);
  
  return 0;
}

int new_particle(int * index_of_the_particle, double x, double y, 
  double z){
  
  *index_of_the_particle = highest_index;
  
  Particle * p = new Particle(highest_index, x, y, z);
  particlesMap[highest_index] = p;
  
  highest_index++;
  
  return 0;
}

int get_close_neighbors(
  int index_of_the_particle, 
  double * index_of_first_neighbor, 
  double * index_of_second_neighbor, 
  double * index_of_third_neighbor
  ){

  if(index_of_the_particle > highest_index) {
    return -1;
  }
  
  Particle * p = particlesMap[index_of_the_particle];
  *index_of_first_neighbor  = p->n0;
  *index_of_second_neighbor = p->n1;
  *index_of_third_neighbor  = p->n2;
  
  return 0;
}

int delete_particle(int index_of_the_particle){
  
  if(index_of_the_particle > highest_index) {
    return -1;
  }
  
  particlesMap.erase(index_of_the_particle);
  
  return 0;
}

int set_state(int index_of_the_particle, double x, double y, double z){
    
  if(index_of_the_particle > highest_index) {
    return -1;
  }
  
  Particle * p = particlesMap[index_of_the_particle];
  
  p->x = x;
  p->y = y;
  p->z = z;
  
  return 0;
}

int get_state(int index_of_the_particle, double * x, double * y, 
  double * z){
      
  if(index_of_the_particle > highest_index) {
    return -1;
  }
  
  Particle * p = particlesMap[index_of_the_particle];
  
  *x = p->x;
  *y = p->y;
  *z = p->z;
  
  return 0;
}

int get_number_of_particles(int * value){
  *value = (int) particlesMap.size();
  return 0;
}
