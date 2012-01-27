#ifndef VECTOR3_H
#define VECTOR3_H

#include<iostream>
#include<iomanip>
#include<cstdlib>
#include<cstdio>

class Vector3
{
 public:
  Vector3(double x, double y, double z){
    v[0]=x; v[1]=y; v[2]=z;
  }

  Vector3(const Vector3& u){
    v[0]=u.v[0]; v[1]=u.v[1]; v[2]=u.v[2];
  }

  Vector3(double c[3]){
    v[0] = c[0]; v[1] = c[1]; v[2] = c[2];
  }

  Vector3(double c=0.0){
    v[0] = c; v[1] = c; v[2] = c;
  }

  ~Vector3(){}

  void def(double x=0.0, double y=0.0, double z=0.0){
    v[0]=x; v[1]=y; v[2]=z;
  }

  void set(double x, double y, double z){
    v[0]=x; v[1]=y; v[2]=z;
  }

  const double& operator[](int i) const {return v[i];}

  double& operator[](int i){return v[i];}

  //double& get_value(int i) const {return v[i];}
  
  /******************************
    arithmetic operator
  ********************************/
  Vector3& operator = (const Vector3& u){
    v[0]=u.v[0]; v[1]=u.v[1]; v[2]=u.v[2];
    return *this;
  }

  Vector3 operator + (const Vector3 &u) const{
    Vector3 x=0.0;
    x.v[0]=v[0]+u.v[0];
    x.v[1]=v[1]+u.v[1];
    x.v[2]=v[2]+u.v[2];
    return x;
  }

  Vector3 operator - (const Vector3 &u) const{
    Vector3 x=0.0;
    x.v[0]=v[0]-u.v[0];
    x.v[1]=v[1]-u.v[1];
    x.v[2]=v[2]-u.v[2];
    return x;
  }

  double operator * (const Vector3 u) const{
    double s=0.0;
    s=v[0]*u.v[0]+v[1]*u.v[1]+v[2]*u.v[2];
    return s;
  }

  Vector3 operator ^ (const Vector3 u) const{
    Vector3 x;
    x.v[0]=v[1]*u.v[2]-v[2]*u.v[1];
    x.v[1]=v[2]*u.v[0]-v[0]*u.v[2];
    x.v[2]=v[0]*u.v[1]-v[1]*u.v[0];
    return x;
  }

  Vector3 operator * (const double s) const{
    Vector3 x;
    x.v[0]=v[0]*s;
    x.v[1]=v[1]*s;
    x.v[2]=v[2]*s;
    return x;
  }
  
  Vector3 operator / (const double s) const{
    Vector3 x;
    x.v[0]=v[0]/s;
    x.v[1]=v[1]/s;
    x.v[2]=v[2]/s;
    return x;
  }


//	Vector +=, -=, *=, /=
  
  inline Vector3& operator += (const Vector3 &u){
    v[0]+=u.v[0]; v[1]+=u.v[1]; v[2]+=u.v[2];
    return *this;
  }

  inline Vector3& operator -= (const Vector3 &u){
    v[0]-=u.v[0]; v[1]-=u.v[1]; v[2]-=u.v[2];
    return *this;
  }
  
  inline Vector3& operator *= (const double b){
    v[0] *= b; v[1] *= b; v[2] *= b;
    return *this;
  }
  
  inline Vector3& operator /= (const double b){
    register double binv = 1.0/b;
    v[0] *= binv; v[1] *= binv; v[2] *= binv;
    return *this;
  }
  
  //friend Vector3 operator + (const Vector3 &u0, const Vector3 &u1);
  //friend Vector3 operator - (const Vector3 &u0, const Vector3 &u1);
  friend Vector3 operator * (const double &s, const Vector3 &u);
  //friend double operator * (const Vector3 &u0, const Vector3 &u1);
  
  
  /******************************
   stream operator
  *****************************/
  
  friend std::ostream& operator << (std::ostream& c, const Vector3& u);
  //don't forget &
  friend std::istream& operator >> (std::istream& c, Vector3& u);
  private:
  double v[3];

  //friend Vector3 atov(char *str);

};

/*
inline Vector3 operator + (const Vector3 &u0, const Vector3 &u1){
  Vector3 x;
  x.v[0]=u0.v[0]+u1.v[0];
  x.v[1]=u0.v[1]+u1.v[1];
  x.v[2]=u0.v[2]+u1.v[2];
  return x;
}
*/
/*
inline Vector3 operator - (const Vector3 &u0, const Vector3 &u1){
  Vector3 x;
  x.v[0]=u0.v[0]-u1.v[0];
  x.v[1]=u0.v[1]-u1.v[1];
  x.v[2]=u0.v[2]-u1.v[2];
  return x;
}
*/

inline Vector3 operator * (const double &s, const Vector3 &u){
  Vector3 x;
  x.v[0]=u.v[0]*s;
  x.v[1]=u.v[1]*s;
  x.v[2]=u.v[2]*s;
  return x;
}

/*
inline double operator * (const Vector3 &u0, const Vector3 &u1){
  double s=0.0;
  s=u0.v[0]*u1.v[0] + u0.v[1]*u1.v[1] + u0.v[2]*u1.v[2];
  return s;
}
*/

inline std::ostream& operator << (std::ostream& c, const Vector3& u){
  c<<std::setprecision(15)<<u.v[0]<<"   "<<std::setprecision(15)<<u.v[1]<<"    "<<std::setprecision(15)<<u.v[2];
  //c<<u.v[0]<<"   "<<u.v[1]<<"    "<<u.v[2];
  return c;
}

inline std::istream& operator >> (std::istream& c, Vector3& u){
  c>>u.v[0]; c>>u.v[1]; c>>u.v[2];
  return c;

}

inline Vector3 atov(char *str){
  Vector3 v=0.0;
  char buf1[1024], buf2[1024], buf3[1024];
  if(sscanf(str,"%s %s %s",buf1,buf2,buf3)<3){std::cerr<<"error atov"<<std::endl;}
  else{
    v[0]=atof(buf1);
    v[1]=atof(buf2);
    v[2]=atof(buf3);
  }
  return v;
}

inline Vector3 atov(char *str1, char *str2, char *str3){
  Vector3 v=0.0;
  v[0]=atof(str1);
  v[1]=atof(str2);
  v[2]=atof(str3);
  return v;
}

#endif
