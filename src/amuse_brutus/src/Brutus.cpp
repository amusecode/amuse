#include "Brutus.h"

Brutus::Brutus() {
  t = "0";
  data.clear();
  N = 0;  

  tolerance = "1e-6";
  numBits = 56;

  eta = get_eta(tolerance);

  setup();
}
Brutus::Brutus(std::vector<mpreal> &data) {
  t = "0";
  this->data = data;
  N = data.size()/7;  

  tolerance = "1e-6";
  numBits = 56;

  eta = get_eta(tolerance);

  setup();
}
Brutus::Brutus(mpreal &t, std::vector<mpreal> &data, mpreal &tolerance) {
  this->t = t;
  this->data = data;
  N = data.size()/7;  

  this->tolerance = tolerance;
  numBits = get_numBits(tolerance);

  eta = get_eta(tolerance);

  setup();
}
Brutus::Brutus(mpreal &t, std::vector<mpreal> &data, mpreal &tolerance, int &numBits) {
  this->t = t;
  this->data = data;
  N = data.size()/7;  

  this->tolerance = tolerance;
  this->numBits = numBits;

  eta = get_eta(tolerance);

  setup();
}

void Brutus::set_data(std::vector<mpreal> &data) {
  this->data = data;
  N = data.size()/7; 

  Cluster c(data);
  cl = c;
}
void Brutus::set_eta(mpreal &eta) {
  this->eta = eta;
}
void Brutus::set_tolerance(mpreal &tolerance) {
  this->tolerance = tolerance;
  bs.set_tolerance(tolerance);
}
void Brutus::set_numBits(int &numBits) {
  this->numBits = numBits;
}
void Brutus::set_t_begin(mpreal &t_begin) {
  this->t = t_begin;
}

mpreal Brutus::get_eta(mpreal tolerance) {
  mpreal a = "0", b = "0";

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
  
  return eta;
}
mpreal Brutus::get_tolerance() {
  return tolerance;
}
int Brutus::get_numBits() {
  return numBits;
}  
int Brutus::get_numBits(mpreal tolerance) {
  mpreal absloge = abs( log10(tolerance) );
  return 4*(int)absloge.toLong()+32;
}
mpreal Brutus::fit_slope(std::vector<mpreal> &x, std::vector<mpreal> &y) {
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

void Brutus::setup() {
  Cluster c(data);
  c.eps2 = "0";
  cl = c;

  Bulirsch_Stoer b(tolerance);
  bs = b;

  eta = get_eta(tolerance);  
}

void Brutus::evolve(mpreal t_end) {
  while (t<t_end) {
    cl.calcAcceleration_dt();
    dt = eta*cl.dt;

    if(t+dt > t_end) dt = t_end-t;

    bool converged = bs.integrate(cl, dt);

    if(!converged) {
      std::cerr << "Not converged at " << t << "!" << std::endl;
      exit(1);
    }

    t += dt;
  }
  this->data = cl.get_data();
}
  
mpreal Brutus::get_t() {
  return t;
}
std::vector<mpreal> Brutus::get_data() {
  return data;
}
std::vector<double> Brutus::get_data_double() {
  int N = data.size()/7;
  std::vector<double> v(7*N, 0);
  for(int i=0; i<N; i++) {
    for(int j=0; j<7; j++) {
      v[i*7+j] = data[i*7+j].toDouble();
    }
  }
  return v;
}
std::vector<std::string> Brutus::get_data_string() {
  int N = data.size()/7;
  std::vector<std::string> v(7*N, "0");
  for(int i=0; i<N; i++) {
    for(int j=0; j<7; j++) {
      v[i*7+j] = data[i*7+j].toString();
    }
  }
  return v;
}




