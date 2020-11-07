#include "Cluster.h"

Cluster::Cluster(vector<double> data) {
  int N = data.size()/arg_cnt;
  s.resize(N);
  mpreal m;
  vector<mpreal> r(3), v(3);
  #ifdef use_additional_acc
  vector<mpreal> a_step(3);
  #endif // use_additional_acc
  for(int i=0; i<N; i++) {
    m    = (mpreal)data[i*arg_cnt+arg_m];
    r[0] = (mpreal)data[i*arg_cnt+arg_r_0];
    r[1] = (mpreal)data[i*arg_cnt+arg_r_1];
    r[2] = (mpreal)data[i*arg_cnt+arg_r_2];
    v[0] = (mpreal)data[i*arg_cnt+arg_v_0];
    v[1] = (mpreal)data[i*arg_cnt+arg_v_1];
    v[2] = (mpreal)data[i*arg_cnt+arg_v_2];
    #ifdef use_additional_acc
    a_step[0] = (mpreal)data[i*arg_cnt+arg_a_step_0];
    a_step[1] = (mpreal)data[i*arg_cnt+arg_a_step_1];
    a_step[2] = (mpreal)data[i*arg_cnt+arg_a_step_2];
    s[i] = Star(m, r, v, a_step);
    #else
    s[i] = Star(m, r, v);
    #endif // use_additional_acc
  }
//  this->time = 0;
}
Cluster::Cluster(vector<mpreal> data) {
  int N = data.size()/arg_cnt;
  s.resize(N);
  mpreal m;
  vector<mpreal> r(3), v(3);
  #ifdef use_additional_acc
  vector<mpreal> a_step(3);
  #endif // use_additional_acc
  for(int i=0; i<N; i++) {
    m    = data[i*arg_cnt+arg_m];
    r[0] = data[i*arg_cnt+arg_r_0];
    r[1] = data[i*arg_cnt+arg_r_1];
    r[2] = data[i*arg_cnt+arg_r_2];
    v[0] = data[i*arg_cnt+arg_v_0];
    v[1] = data[i*arg_cnt+arg_v_1];
    v[2] = data[i*arg_cnt+arg_v_2];
    #ifdef use_additional_acc
    a_step[0] = data[i*arg_cnt+arg_a_step_0];
    a_step[1] = data[i*arg_cnt+arg_a_step_1];
    a_step[2] = data[i*arg_cnt+arg_a_step_2];
    s[i] = Star(m, r, v, a_step);
    #else
    s[i] = Star(m, r, v);
    #endif // use_additional_acc
  }
//  this->time = 0;
}

  mpreal dx; //constructor sets to zero by default = "0";
  mpreal dy; //constructor sets to zero by default = "0";
  mpreal dz; //constructor sets to zero by default = "0";
  mpreal RdotR; //constructor sets to zero by default = "0";
  mpreal apre; //constructor sets to zero by default = "0";
  mpreal fm; //constructor sets to zero by default = "0";
  mpreal daix; //constructor sets to zero by default = "0";
  mpreal daiy; //constructor sets to zero by default = "0";
  mpreal daiz; //constructor sets to zero by default = "0";
  mpreal a2i; //constructor sets to zero by default = "0";
  mpreal mydti; //constructor sets to zero by default = "0";


void Cluster::calcAcceleration_dt() {
  int N = s.size();
  #ifdef use_additional_acc
     for(int i=0; i<N; i++) {
        for(int j=0; j<3; j++){
    DEBUG::DEB << s[i].a_step[j].toString() << "  "<< N<<"  "<<i<< "  "<< j<< "Cluster\n"; DEBUG::DEB.flush();
           s[i].a[j]=s[i].a_step[j];
        }
     }
  #else
  for(int i=0; i<N; i++) {
    s[i].a.assign(3, "0");
  }
  #endif // use_additional_acc
  dt.setInf(+1); // = "1e100";
/*
  mpreal dx; //constructor sets to zero by default = "0";
  mpreal dy; //constructor sets to zero by default = "0";
  mpreal dz; //constructor sets to zero by default = "0";
  mpreal RdotR; //constructor sets to zero by default = "0";
  mpreal apre; //constructor sets to zero by default = "0";
  mpreal fm; //constructor sets to zero by default = "0";
  mpreal daix; //constructor sets to zero by default = "0";
  mpreal daiy; //constructor sets to zero by default = "0";
  mpreal daiz; //constructor sets to zero by default = "0";
  mpreal a2i; //constructor sets to zero by default = "0";
  mpreal mydti; //constructor sets to zero by default = "0";
*/
  for(int i=0; i<N-1; i++) {
    for(int j=i+1; j<N; j++) {

      dx = s[j].x[0]-s[i].x[0];
      dy = s[j].x[1]-s[i].x[1];
      dz = s[j].x[2]-s[i].x[2];
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
     for(int i=0; i<N; i++) {
        for(int j=0; j<3; j++){
    DEBUG::DEB << s[i].a[j].toString() << "  "<< N<<"  "<<i<< "  "<< j<< "acc\n"; DEBUG::DEB.flush();
        }
     }
  dt = pow(dt, "0.25");
}
void Cluster::calcAcceleration() {
  int N = s.size();
  #ifdef use_additional_acc
     for(int i=0; i<N; i++) {
        for(int j=0; j<3; j++){
           s[i].a[j]=s[i].a_step[j];
        }
     }
  #else
  for(int i=0; i<N; i++) {
    s[i].a.assign(3, "0");
  }
  #endif // use_additional_acc
/*
  mpreal dx; //constructor sets to zero by default = "0";
  mpreal dy; //constructor sets to zero by default = "0";
  mpreal dz; //constructor sets to zero by default = "0";
  mpreal RdotR; //constructor sets to zero by default = "0";
  mpreal apre; //constructor sets to zero by default = "0";
  mpreal fm; //constructor sets to zero by default = "0";
  mpreal daix; //constructor sets to zero by default = "0";
  mpreal daiy; //constructor sets to zero by default = "0";
  mpreal daiz; //constructor sets to zero by default = "0";
  mpreal a2i; //constructor sets to zero by default = "0";
  mpreal mydti; //constructor sets to zero by default = "0";
*/
  for(int i=0; i<N-1; i++) {
    for(int j=i+1; j<N; j++) {
      dx = s[j].x[0]-s[i].x[0];
      dy = s[j].x[1]-s[i].x[1];
      dz = s[j].x[2]-s[i].x[2];
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
      s[i].x[k] += dt*s[i].v[k] + "0.5"*dt*dt*s[i].a0[k];
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
/*
vector<mpreal> Cluster::energies() {
  mpreal init; //constructor sets to zero by default = "0";
  vector<mpreal> E(3), rij(3);
  E.assign(3,"0");

  for (vector<Star>::iterator si = s.begin(); si != s.end(); ++si) {
    E[1] += "0.5"*si->m*inner_product(si->v.begin(), si->v.end(), si->v.begin(), init);
  }

  for (vector<Star>::iterator si = s.begin(); si != s.end(); ++si) {
    for (vector<Star>::iterator sj = si+1; sj != s.end(); ++sj) {
      for (int i = 0; i != 3; ++i)
        rij[i] = si->x[i]-sj->x[i];
      E[2] -= si->m*sj->m/sqrt(inner_product(rij.begin(), rij.end(),  rij.begin(), init));
    }
  }
  E[0] = E[1] + E[2];
  return E;
}
*/

/*vector<double> Cluster::get_data_double() {
  vector<double> ddata;
  for (vector<Star>::iterator si = s.begin(); si != s.end(); ++si) {
    ddata.push_back(si->m.toDouble());
    ddata.push_back(si->x[0].toDouble());
    ddata.push_back(si->x1].toDouble());
    ddata.push_back(si->x[2].toDouble());
    ddata.push_back(si->v[0].toDouble());
    ddata.push_back(si->v[1].toDouble());
    ddata.push_back(si->v[2].toDouble());
    #ifdef use_additional_acc
    ddata.push_back(si->a_step[0].toDouble());
    ddata.push_back(si->a_step[1].toDouble());
    ddata.push_back(si->a_step[2].toDouble());
    #endif // use_additional_acc
  }
  return ddata;
}*/
vector<mpreal> Cluster::get_data() {
  vector<mpreal> ddata;
  for (vector<Star>::iterator si = s.begin(); si != s.end(); ++si) {
    ddata.push_back(si->m);
    ddata.push_back(si->x[0]);
    ddata.push_back(si->x[1]);
    ddata.push_back(si->x[2]);
    ddata.push_back(si->v[0]);
    ddata.push_back(si->v[1]);
    ddata.push_back(si->v[2]);
    #ifdef use_additional_acc
    ddata.push_back(si->a_step[0]);
    ddata.push_back(si->a_step[1]);
    ddata.push_back(si->a_step[2]);
    #endif // use_additional_acc
  }
  return ddata;
}

#ifdef use_additional_acc
void Cluster::remove_step_Acceleration() {
  int N = s.size();
  for(int i=0; i<N; i++) {
    for(int j=0; j<3;j++)
        s[i].a_step[j].setZero() ;
  }
}
#endif // use_additional_acc
