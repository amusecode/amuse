#ifndef MATRIX3_H
#define MATRIX3_H
#include<iostream>
#include<iomanip>
#include<cmath>
#include"Vector3.h"

class Matrix3
{
 private:
  double m[3][3];

 public:
  Matrix3(double c=0.0){
    m[0][0] = c; m[0][1] = c; m[0][2] = c;
    m[1][0] = c; m[1][1] = c; m[1][2] = c;
    m[2][0] = c; m[2][1] = c; m[2][2] = c;
  }

  Matrix3(Vector3 &v0, Vector3 &v1, Vector3 &v2){
    m[0][0] = v0[0]; m[0][1] = v1[0]; m[0][2] = v2[0];
    m[1][0] = v0[1]; m[1][1] = v1[1]; m[1][2] = v2[1];
    m[2][0] = v0[2]; m[2][1] = v1[2]; m[2][2] = v2[2];
  }

  Matrix3(const Matrix3& c){
    m[0][0] = c.m[0][0]; m[0][1] = c.m[0][1]; m[0][2] = c.m[0][2];
    m[1][0] = c.m[1][0]; m[1][1] = c.m[1][1]; m[1][2] = c.m[1][2];
    m[2][0] = c.m[2][0]; m[2][1] = c.m[2][1]; m[2][2] = c.m[2][2];
  }

  ~Matrix3(){}

  double* operator[](int i){return m[i];}

  Matrix3& operator=(const Matrix3 &c){
    m[0][0] = c.m[0][0]; m[0][1] = c.m[0][1]; m[0][2] = c.m[0][2];
    m[1][0] = c.m[1][0]; m[1][1] = c.m[1][1]; m[1][2] = c.m[1][2];
    m[2][0] = c.m[2][0]; m[2][1] = c.m[2][1]; m[2][2] = c.m[2][2];
    return *this;
  }

  Matrix3 operator+(const Matrix3 &a){
    Matrix3 x;
    x.m[0][0]=m[0][0]+a.m[0][0]; x.m[0][1]=m[0][1]+a.m[0][1];  x.m[0][2]=m[0][2]+a.m[0][2];
    x.m[1][0]=m[1][0]+a.m[1][0]; x.m[1][1]=m[1][1]+a.m[1][1];  x.m[1][2]=m[1][2]+a.m[1][2];
    x.m[2][0]=m[2][0]+a.m[2][0]; x.m[2][1]=m[2][1]+a.m[2][1];  x.m[2][2]=m[2][2]+a.m[2][2];
    return x;
  }

  Matrix3 operator-(const Matrix3 &a){
    Matrix3 x;
    x.m[0][0]=m[0][0]-a.m[0][0]; x.m[0][1]=m[0][1]-a.m[0][1];  x.m[0][2]=m[0][2]-a.m[0][2];
    x.m[1][0]=m[1][0]-a.m[1][0]; x.m[1][1]=m[1][1]-a.m[1][1];  x.m[1][2]=m[1][2]-a.m[1][2];
    x.m[2][0]=m[2][0]-a.m[2][0]; x.m[2][1]=m[2][1]-a.m[2][1];  x.m[2][2]=m[2][2]-a.m[2][2];
    return x;
  }

  Matrix3 operator*(const Matrix3 &a){
    Matrix3 x;
    x.m[0][0] = m[0][0]*a.m[0][0] + m[0][1]*a.m[1][0] + m[0][2]*a.m[2][0];
    x.m[0][1] = m[0][0]*a.m[0][1] + m[0][1]*a.m[1][1] + m[0][2]*a.m[2][1];
    x.m[0][2] = m[0][0]*a.m[0][2] + m[0][1]*a.m[1][2] + m[0][2]*a.m[2][2];

    x.m[1][0] = m[1][0]*a.m[0][0] + m[1][1]*a.m[1][0] + m[1][2]*a.m[2][0];
    x.m[1][1] = m[1][0]*a.m[0][1] + m[1][1]*a.m[1][1] + m[1][2]*a.m[2][1];
    x.m[1][2] = m[1][0]*a.m[0][2] + m[1][1]*a.m[1][2] + m[1][2]*a.m[2][2];

    x.m[2][0] = m[2][0]*a.m[0][0] + m[2][1]*a.m[1][0] + m[2][2]*a.m[2][0];
    x.m[2][1] = m[2][0]*a.m[0][1] + m[2][1]*a.m[1][1] + m[2][2]*a.m[2][1];
    x.m[2][2] = m[2][0]*a.m[0][2] + m[2][1]*a.m[1][2] + m[2][2]*a.m[2][2];
    return x;
  }


  Vector3 operator * (Vector3 &a){
    Vector3 x;
    x[0] = m[0][0]*a[0] + m[0][1]*a[1] + m[0][2]*a[2];
    x[1] = m[1][0]*a[0] + m[1][1]*a[1] + m[1][2]*a[2];
    x[2] = m[2][0]*a[0] + m[2][1]*a[1] + m[2][2]*a[2];
    return x;
  }

  Matrix3 operator * (const double s) const{
    Matrix3 x;
    x.m[0][0]=m[0][0]*s; x.m[0][1]=m[0][1]*s; x.m[0][2]=m[0][2]*s;
    x.m[1][0]=m[1][0]*s; x.m[1][1]=m[1][1]*s; x.m[1][2]=m[1][2]*s;
    x.m[2][0]=m[2][0]*s; x.m[2][1]=m[2][1]*s; x.m[2][2]=m[2][2]*s;
    return x;
  }

  Matrix3 transposed(){
    Matrix3 mtmp;
    for(int i=0; i<3; i++){
      for(int j=0; j<3; j++){
	mtmp.m[j][i] = m[i][j];
      }
    }
    return mtmp;
  }

  void unit(){
    for(int i=0; i<3; i++){
      for(int j=0; j<3; j++){
	if(i != j){m[i][j]=0.0;}
	else{m[i][i] = 1.0;}
      }
    }
  }

  Matrix3& operator += (const Matrix3 &a){
    Matrix3 x;
    m[0][0] += a.m[0][0]; m[0][1] += a.m[0][1];  m[0][2] += a.m[0][2];
    m[1][0] += a.m[1][0]; m[1][1] += a.m[1][1];  m[1][2] += a.m[1][2];
    m[2][0] += a.m[2][0]; m[2][1] += a.m[2][1];  m[2][2] +=a.m[2][2];
    return *this;
  }

  Matrix3& operator -= (const Matrix3 &a){
    Matrix3 x;
    m[0][0] -= a.m[0][0]; m[0][1] -= a.m[0][1];  m[0][2] -= a.m[0][2];
    m[1][0] -= a.m[1][0]; m[1][1] -= a.m[1][1];  m[1][2] -= a.m[1][2];
    m[2][0] -= a.m[2][0]; m[2][1] -= a.m[2][1];  m[2][2] -=a.m[2][2];
    return *this;
  }

  Matrix3& operator *= (const double &a){
    Matrix3 x;
    m[0][0] *= a; m[0][1] *= a;  m[0][2] *= a;
    m[1][0] *= a; m[1][1] *= a;  m[1][2] *= a;
    m[2][0] *= a; m[2][1] *= a;  m[2][2] *= a;
    return *this;
  }

  void rotation(double I, double OMEGA, double omega){
    m[0][0]=cos(omega)*cos(OMEGA) - sin(omega)*sin(OMEGA)*cos(I);
    m[0][1]=-sin(omega)*cos(OMEGA) - cos(omega)*sin(OMEGA)*cos(I);
    m[0][2]=sin(OMEGA)*sin(I);

    m[1][0]=cos(omega)*sin(OMEGA) + sin(omega)*cos(OMEGA)*cos(I);
    m[1][1]=-sin(omega)*sin(OMEGA) + cos(omega)*cos(OMEGA)*cos(I);
    m[1][2]=-cos(OMEGA)*sin(I);

    m[2][0]=sin(omega)*sin(I);
    m[2][1]=cos(omega)*sin(I);
    m[2][2]=cos(I);
  }

  double determinant(){
    double det = 0.0;
    det = m[0][0] * m[1][1] * m[2][2] 
      + m[0][1] * m[1][2] * m[2][0] 
      + m[0][2] * m[1][0] * m[2][1] 
      - m[0][2] * m[1][1] * m[2][0] 
      - m[0][1] * m[1][0] * m[2][2] 
      - m[0][0] * m[1][2] * m[2][1];
    return det;
  }

  friend Matrix3 outer_product(const Vector3 &a, const Vector3 &b);

  friend Matrix3 quadrupole(const Vector3 &a);

  friend Matrix3 inertiamoment(const double m, const Vector3 &a);

  friend Vector3 operator * (const Vector3 &a, const Matrix3 &c);

  friend Matrix3 operator * (const double &s, const Matrix3 &c);

  friend std::ostream& operator << (std::ostream& c, Matrix3& mtmp);
};

inline std::ostream& operator << (std::ostream& c, Matrix3& mtmp){
  c<<std::setprecision(15)<<mtmp.m[0][0]<<"   "<<std::setprecision(15)<<mtmp.m[0][1]<<"    "<<std::setprecision(15)<<mtmp.m[0][2]<<std::endl;
  c<<std::setprecision(15)<<mtmp.m[1][0]<<"   "<<std::setprecision(15)<<mtmp.m[1][1]<<"    "<<std::setprecision(15)<<mtmp.m[1][2]<<std::endl;
  c<<std::setprecision(15)<<mtmp.m[2][0]<<"   "<<std::setprecision(15)<<mtmp.m[2][1]<<"    "<<std::setprecision(15)<<mtmp.m[2][2]<<std::endl;
  //c<<u.v[0]<<"   "<<u.v[1]<<"    "<<u.v[2];
  return c;
}

inline Matrix3 outer_product(const Vector3 &a, const Vector3 &b){
  Matrix3 x;
  x[0][0]=a[0]*b[0]; x[0][1]=a[0]*b[1]; x[0][2]=a[0]*b[2];
  x[1][0]=a[1]*b[0]; x[1][1]=a[1]*b[1]; x[1][2]=a[1]*b[2];
  x[2][0]=a[2]*b[0]; x[2][1]=a[2]*b[1]; x[2][2]=a[2]*b[2];
  return x;
}


inline Matrix3 quadrupole(const Vector3 &a){
  Matrix3 x=3.0*outer_product(a,a);
  double a_squared=a*a;
  Matrix3 delta;
  delta.unit();
  return x-a_squared*delta;
}

inline Matrix3 inertiamoment(const double mass, const Vector3 &pos){
  Matrix3 x = mass * outer_product(pos, pos);
  return x;
}

inline Vector3 operator*(const Vector3 &a, const Matrix3 &c){
  Vector3 x;
  x[0]=a[0]*c.m[0][0] + a[1]*c.m[1][0] + a[2]*c.m[2][0];
  x[1]=a[0]*c.m[0][1] + a[1]*c.m[1][1] + a[2]*c.m[2][1];
  x[2]=a[0]*c.m[0][2] + a[1]*c.m[1][2] + a[2]*c.m[2][2];
  return x;
}

inline Matrix3 operator * (const double &s, const Matrix3 &c){
  Matrix3 x;
  x[0][0]=c.m[0][0]*s; x[0][1]=c.m[0][1]*s; x[0][2]=c.m[0][2]*s;
  x[1][0]=c.m[1][0]*s; x[1][1]=c.m[1][1]*s; x[1][2]=c.m[1][2]*s;
  x[2][0]=c.m[2][0]*s; x[2][1]=c.m[2][1]*s; x[2][2]=c.m[2][2]*s;
  return x;
}

/*
inline void inertial_moment(const double &mass, const Vector3 &pos, Matrix3 &I){
  I[0][0] += mass*pos[0]*pos[0]; I[0][1] += mass*pos[0]*pos[1]; I[0][2] += mass*pos[0]*pos[2];
  I[1][1] += mass*pos[1]*pos[1]; I[1][2] += mass*pos[1]*pos[2];
  I[2][2] += mass*pos[2]*pos[2];
  I[1][0] = I[0][1];
  I[2][0] = I[0][2];
  I[2][1] = I[1][2];
}
*/

#endif
