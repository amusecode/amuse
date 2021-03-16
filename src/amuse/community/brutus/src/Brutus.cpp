#include "Brutus.h"

Brutus::Brutus() {
//DEBUG::DEB <<  "Brutus::Brutus()\n"; DEBUG::DEB.flush();
  t = "0";
  data.clear();
  N = 0;

  tolerance = "1e-6";
//  numBits = 56;
  eta= new mpreal();
  *eta = get_eta(tolerance);

  setup();
}
Brutus::Brutus(vector<mpreal> &data) {
//DEBUG::DEB <<  "Brutus::Brutus(vector<mpreal> &data)\n"; DEBUG::DEB.flush();
  t = "0";
  this->data = data;
  N = data.size()/arg_cnt;

  tolerance = "1e-6";
//  numBits = 56;
  eta= new mpreal();
  *eta = get_eta(tolerance);

  setup();
}
Brutus::Brutus(mpreal &t, vector<mpreal> &data, mpreal &tolerance) {
//DEBUG::DEB <<  "Brutus::Brutus(mpreal &t, vector<mpreal> &data, mpreal &tolerance)\n"; DEBUG::DEB.flush();
  this->t = t;
  this->data = data;
  N = data.size()/arg_cnt;

  this->tolerance = tolerance;
//  numBits = get_numBits(tolerance);
  eta= new mpreal();
  *eta = get_eta(tolerance);

  setup();
}
Brutus::Brutus(mpreal &t, vector<mpreal> &data, mpreal &tolerance, int &numBits) {
//DEBUG::DEB <<  "Brutus::Brutus(mpreal &t, vector<mpreal> &data, mpreal &tolerance, int &numBits)\n"; DEBUG::DEB.flush();
  this->t = t;
  this->data = data;
  N = data.size()/arg_cnt;

  this->tolerance = tolerance;
//  this->numBits = numBits;
  eta= new mpreal();
  *eta = get_eta(tolerance);

  setup();
}

void Brutus::set_data(vector<mpreal> &data) {
//DEBUG::DEB <<  "Brutus::set_data(vector<mpreal> &data)\n"; DEBUG::DEB.flush();
  this->data = data;

  N = data.size()/arg_cnt;


//Star::DEBUG <<"\n"; Star::DEBUG.flush();
  Cluster c(data);
  cl = c;
}
void Brutus::set_eta(mpreal p_eta) {
  eta_was_set=true;
  *eta = p_eta;
  //DEBUG::DEB <<  "Brutus::set_eta  "<< *eta <<"\n"; DEBUG::DEB.flush();
}
void Brutus::set_tolerance(mpreal &tolerance) {
//DEBUG::DEB <<  "Brutus::set_tolerance(mpreal &tolerance)\n"; DEBUG::DEB.flush();
  this->tolerance = tolerance;
  bs.set_tolerance(tolerance);
}
/*void Brutus::set_numBits(int &numBits) {
  this->numBits = numBits;
}*/
void Brutus::set_t_begin(mpreal &t_begin) {
//DEBUG::DEB <<  "Brutus::set_t_begin(mpreal &t_begin)\n"; DEBUG::DEB.flush();
  this->t = t_begin;
}

mpreal Brutus::get_eta(mpreal tolerance) {
  mpreal a = "0", b = "0";
//DEBUG::DEB <<  "Brutus::get_eta(mpreal tolerance)" << tolerance.toString()<<"  "<<this->eta->toString() << "\n"; DEBUG::DEB.flush();
  if(eta_was_set)
    return *eta;
  if(tolerance > "1e-50") {
    a = "-0.012";
    b = "-0.40";
  }
  else {
    a = "-0.029";
    b = "0.45";
  }

  mpreal abslogtol = abs( log10(tolerance) );
  mpreal logeta = a*abslogtol+b;
  mpreal eta = pow("10", logeta);

//DEBUG::DEB <<  "Brutus::get_eta(mpreal tolerance)" << tolerance.toString()<<"  "<<eta.toString() << "\n"; DEBUG::DEB.flush();
  return eta;
}
mpreal Brutus::get_tolerance() {
//DEBUG::DEB <<  "Brutus::get_tolerance()\n"; DEBUG::DEB.flush();
  return tolerance;
}
/*int Brutus::get_numBits() {
  return numBits;
}*/
/*int Brutus::get_numBits(mpreal tolerance) {
  mpreal absloge = abs( log10(tolerance) );
  return 4*(int)absloge.toLong()+32;
}*/
/*
mpreal Brutus::fit_slope(vector<mpreal> &x, vector<mpreal> &y) {
  mpreal a = "0";
  int N = x.size();

  if(N < 2) {
    a = "-1e100";
  }
  else {
    mpreal SX="0", SY="0", SX2="0", SXY="0", Xav="0", Yav="0";
    for(int i=0; i<N; i++) {
      SX += x[i];
      SY += y[i];
      SX2 += x[i]*x[i];
      SXY += x[i]*y[i];
    }
    Xav = SX / (mpreal)N;
    Yav = SY / (mpreal)N;
    a = (SXY - SX * Yav) / (SX2 - SX * Xav);
  }

  return a;
}
*/
void Brutus::setup() {
//DEBUG::DEB <<  "Brutus::setup()\n"; DEBUG::DEB.flush();
  Cluster c(data);
  c.eps2 = "0";
  cl = c;

  Bulirsch_Stoer b(tolerance);
  bs = b;

  *eta = get_eta(tolerance);
}

void Brutus::evolve(mpreal t_end) {
  while (t<t_end) {
    cl.calcAcceleration_dt();
    dt = *eta*cl.dt;
//DEBUG::DEB <<  "Brutus::evolve  "<< eta->toString() <<"\n"; DEBUG::DEB.flush();
    if(t+dt > t_end) dt = t_end-t;

    bool converged = bs.integrate(cl, dt);

    if(!converged) {
      cerr << "Not converged at " << t << "!" << endl;
      exit(1);
    }

    t += dt;
  }
  this->data = cl.get_data();
}

mpreal Brutus::get_t() {
//DEBUG::DEB <<  "Brutus::get_t()\n"; DEBUG::DEB.flush();
  return t;
}
vector<mpreal> Brutus::get_data() {
//DEBUG::DEB <<  "Brutus::get_data()\n"; DEBUG::DEB.flush();
  return data;
}
vector<double> Brutus::get_data_double() {
//DEBUG::DEB <<  "Brutus::get_data_double()\n"; DEBUG::DEB.flush();
  int N = data.size()/arg_cnt;
  vector<double> v(arg_cnt*N, 0);
  for(int i=0; i<N; i++) {
    for(int j=0; j<arg_cnt; j++) {
      v[i*arg_cnt+j] = data[i*arg_cnt+j].toDouble();
    }
  }
  return v;
}
vector<string> Brutus::get_data_string() {
//DEBUG::DEB <<  "Brutus::get_data_string()\n"; DEBUG::DEB.flush();
  int N = data.size()/arg_cnt;
  vector<string> v(arg_cnt*N, "0");
  for(int i=0; i<N; i++) {
    for(int j=0; j<arg_cnt; j++) {
      v[i*arg_cnt+j] = data[i*arg_cnt+j].toString();
    }
  }
  return v;
}




