#include "Cluster.h"

Cluster::Cluster(vector<double> data) {
  int N = data.size()/7;
  s.resize(N);
  mpreal m;
  vector<mpreal> r(3), v(3);
  for(int i=0; i<N; i++) {
    m    = (mpreal)data[i*7+0];
    r[0] = (mpreal)data[i*7+1];
    r[1] = (mpreal)data[i*7+2];
    r[2] = (mpreal)data[i*7+3];
    v[0] = (mpreal)data[i*7+4];
    v[1] = (mpreal)data[i*7+5];
    v[2] = (mpreal)data[i*7+6];
    s[i] = Star(m, r, v);
  }
  this->time = 0;
}
Cluster::Cluster(vector<mpreal> data) {
  int N = data.size()/7;
  s.resize(N);
  mpreal m;
  vector<mpreal> r(3), v(3);
  for(int i=0; i<N; i++) {
    m    = data[i*7+0];
    r[0] = data[i*7+1];
    r[1] = data[i*7+2];
    r[2] = data[i*7+3];
    v[0] = data[i*7+4];
    v[1] = data[i*7+5];
    v[2] = data[i*7+6];
    s[i] = Star(m, r, v);
  }
  this->time = 0;
}

void Cluster::calcAcceleration_dt() {
  int N = s.size();

  for(int i=0; i<N; i++) {
    s[i].a.assign(3, "0");
  }

  dt = "1e100";

  mpreal dx = "0";
  mpreal dy = "0";
  mpreal dz = "0";
  mpreal RdotR = "0";
  mpreal apre = "0";
  mpreal fm = "0";
  mpreal daix = "0";
  mpreal daiy = "0";
  mpreal daiz = "0";
  mpreal a2i = "0";
  mpreal mydti = "0";

  for(int i=0; i<N-1; i++) {
    for(int j=i+1; j<N; j++) {

      dx = s[j].r[0]-s[i].r[0];
      dy = s[j].r[1]-s[i].r[1];
      dz = s[j].r[2]-s[i].r[2];
      RdotR = dx*dx + dy*dy + dz*dz + eps2;
      apre = "1"/sqrt(RdotR*RdotR*RdotR);      

      fm = s[j].m*apre;
      daix = fm*dx;
      daiy = fm*dy;
      daiz = fm*dz;

      s[i].a[0] += daix;
      s[i].a[1] += daiy;
      s[i].a[2] += daiz;

      a2i = daix*daix+daiy*daiy+daiz*daiz;
      mydti = RdotR/a2i;
      if(mydti < dt) dt = mydti;

      fm = s[i].m*apre;
      daix = fm*dx;
      daiy = fm*dy;
      daiz = fm*dz;
      s[j].a[0] -= daix;
      s[j].a[1] -= daiy;
      s[j].a[2] -= daiz;

      a2i = daix*daix+daiy*daiy+daiz*daiz;
      mydti = RdotR/a2i;
      if(mydti < dt) dt = mydti;
    }
  }

  dt = pow(dt, "0.25");
}
void Cluster::calcAcceleration() {
  int N = s.size();

  for(int i=0; i<N; i++) {
    s[i].a.assign(3, "0");
  }

  mpreal dx = "0";
  mpreal dy = "0";
  mpreal dz = "0";
  mpreal RdotR = "0";
  mpreal apre = "0";
  mpreal fm = "0";
  mpreal daix = "0";
  mpreal daiy = "0";
  mpreal daiz = "0";
  mpreal a2i = "0";
  mpreal mydti = "0";

  for(int i=0; i<N-1; i++) {
    for(int j=i+1; j<N; j++) {
      dx = s[j].r[0]-s[i].r[0];
      dy = s[j].r[1]-s[i].r[1];
      dz = s[j].r[2]-s[i].r[2];
      RdotR = dx*dx + dy*dy + dz*dz + eps2;
      apre = "1"/sqrt(RdotR*RdotR*RdotR);      

      fm = s[j].m*apre;
      daix = fm*dx;
      daiy = fm*dy;
      daiz = fm*dz;
      s[i].a[0] += daix;
      s[i].a[1] += daiy;
      s[i].a[2] += daiz;

      fm = s[i].m*apre;
      daix = fm*dx;
      daiy = fm*dy;
      daiz = fm*dz;
      s[j].a[0] -= daix;
      s[j].a[1] -= daiy;
      s[j].a[2] -= daiz;
    }
  }
}
  
void Cluster::updatePositions(mpreal dt) {
  int N = s.size();
  for(int i=0; i<N; i++) {
    s[i].a0 = s[i].a;
    for(int k = 0; k != 3; ++k)
      s[i].r[k] += dt*s[i].v[k] + "0.5"*dt*dt*s[i].a0[k];
  }
}

void Cluster::updateVelocities(mpreal dt) {
  int N = s.size();
  for(int i=0; i<N; i++) {
    for(int k = 0; k != 3; ++k) {
       s[i].v[k] += "0.5"*dt*(s[i].a0[k]+s[i].a[k]);
    }
  }  
}
void Cluster::step(mpreal &dt) {
  updatePositions(dt);
  calcAcceleration();
  updateVelocities(dt);
}
    
vector<mpreal> Cluster::energies() {
  mpreal init = "0";
  vector<mpreal> E(3), rij(3);
  E.assign(3,"0");
    
  for (vector<Star>::iterator si = s.begin(); si != s.end(); ++si) {
    E[1] += "0.5"*si->m*inner_product(si->v.begin(), si->v.end(), si->v.begin(), init); 
  }
    
  for (vector<Star>::iterator si = s.begin(); si != s.end(); ++si) {
    for (vector<Star>::iterator sj = si+1; sj != s.end(); ++sj) {
      for (int i = 0; i != 3; ++i) 
        rij[i] = si->r[i]-sj->r[i];
      E[2] -= si->m*sj->m/sqrt(inner_product(rij.begin(), rij.end(),  rij.begin(), init)); 
    }
  }
  E[0] = E[1] + E[2];
  return E;
}

vector<double> Cluster::get_data_double() {
  vector<double> ddata;  
  for (vector<Star>::iterator si = s.begin(); si != s.end(); ++si) {
    ddata.push_back(si->m.toDouble());
    ddata.push_back(si->r[0].toDouble());
    ddata.push_back(si->r[1].toDouble());
    ddata.push_back(si->r[2].toDouble());
    ddata.push_back(si->v[0].toDouble());
    ddata.push_back(si->v[1].toDouble());
    ddata.push_back(si->v[2].toDouble());	  
  }	  
  return ddata;
}
vector<mpreal> Cluster::get_data() {
  vector<mpreal> ddata;  
  for (vector<Star>::iterator si = s.begin(); si != s.end(); ++si) {
    ddata.push_back(si->m);
    ddata.push_back(si->r[0]);
    ddata.push_back(si->r[1]);
    ddata.push_back(si->r[2]);
    ddata.push_back(si->v[0]);
    ddata.push_back(si->v[1]);
    ddata.push_back(si->v[2]);	  
  }	  
  return ddata;
}

