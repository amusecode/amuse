#ifndef STELLAR_EVOLUTION_H
#define STELLAR_EVOLUTION_H

/*
Tout et al MNRAS 281, 257 (1996) 
*/

static const double L[7][5] = {{3.970417e-01, -3.2913574e-01, 3.4776688e-01, 3.7470851e-01, 9.011915e-02},
			       {8.527626e+00, -2.441225973e+01, 5.643597107e+01, 3.706152575e+01, 5.4562406e+00},
			       {2.5546e-04, -1.23461e-03, -2.3246e-04,  4.5519e-04, 1.6176e-04},
			       {5.432889e+00, -8.62157806e+00, 1.344202049e+01, 1.451584135e+01, 3.39793084e+00},
			       {5.563579e+00,-1.032345224e+01, 1.944322980e+01, 1.897361347e+01, 4.16903097e+00},
			       {7.8866060e-01, -2.90870942e+00,  6.54713531e+00, 4.05606657e+00, 5.3287322e-01},
			       {5.86685e-03, -1.704237e-02, 3.872348e-02, 2.570041e-02, 3.83376e-03}};

static const double R[9][5] = {{1.715359e+00, 6.2246212e-01, -9.2557761e-01, -1.16996966e+00, -3.0631491e-01},
			       {6.597788e+00, -4.2450044e-01,-1.213339427e+01,-1.073509484e+01, -2.51487077e+00},
			       {1.008855000e+01, -7.11727086e+00,-3.167119479e+01, -2.424848322e+01,-5.33608972e+00},
			       {1.012495e+00, 3.2699690e-01, -9.23418e-03, -3.876858e-02, -4.12750e-03},
			       {7.490166e-02, 2.410413e-02, 7.233664e-02, 3.040467e-02, 1.97741e-03}, 
			       {1.077422e-02, 0.0, 0.0, 0.0, 0.0},
			       {3.082234e+00, 9.447205e-01, -2.15200882e+00, -2.49219496e+00, -6.3848738e-01},
			       {1.784778e+01, -7.4534569e+00,-4.896066856e+01,-4.005386135e+01, -9.09331816e+00},
			       {2.2582e-04, -1.86899e-03, 3.88783e-03, 1.42402e-03,-7.671e-05}};


/* sloar metalicity functions
double get_radius_zams0(const double &m)
double get_luminosity_zams0(const double &m)
*/
inline double get_radius_zams0(const double &m){
  double msqrt = sqrt(m);
  double m2 = m * m;  
  double m2_5 = m2 * msqrt;  
  double m4 = m2 * m2;  
  double m6 = m4 * m2;
  double m6_5 = m6 * msqrt;
  double m8_5 = m6_5 * m2;
  double m11 = m8_5 * m2_5;
  double m18_5 = m8_5 * m6 * m4;
  double m19 = m18_5 * msqrt;
  double m19_5 = m19 * msqrt;
  return (R[0][0]*m2_5 + R[1][0]*m6_5 + R[2][0]*m11 + R[3][0]*m19 + R[4][0]*m19_5) / (R[5][0] + R[6][0]*m2 + R[7][0]*m8_5 + m18_5 + R[8][0]*m19_5);
}

inline double get_luminosity_zams0(const double &m){
  double msqrt = sqrt(m);
  double m2 = m * m;  
  double m3 = m * m2;  
  double m5 = m3 * m2;
  double m5_5 = m5 * msqrt;
  double m7 = m5 * m2;
  double m8 = m5 * m3;
  double m9_5 = m5_5 * m2 * m2; 
  double m11 = m8 * m3;
  return (L[0][0]*m5_5 + L[1][0]*m11) / (L[2][0] + m3 + L[3][0]*m5 + L[4][0]*m7 + L[5][0]*m8 + L[6][0]*m9_5);
}



/*
below functions have not be checked
*/

inline double get_coefficient_zams(const double &z, double c[5]){
  double zlog = log10(z*50); // = z / z_sun,  z_sun = 0.02
  return (((c[4]*zlog + c[3]) * zlog + c[2]) * zlog + c[1]) * zlog +c[0];
} 

inline void get_L_coef_zams_all(const double &z, double Lcoef[7]){
  double zlog = log10(z*50); // = z / z_sun,  z_sun = 0.02
  for(int i=0; i<7; i++){
    Lcoef[i] = (((L[i][4]*zlog + L[i][3]) * zlog + L[i][2]) * zlog + L[i][1]) * zlog + L[i][0];
  }
} 

inline void get_R_coef_zams_all(const double &z, double Rcoef[9]){
  double zlog = log10(z*50); // = z / z_sun,  z_sun = 0.02
  for(int i=0; i<9; i++){
    Rcoef[i] = (((R[i][4]*zlog + R[i][3]) * zlog + R[i][2]) * zlog + R[i][1]) * zlog + R[i][0];
  }
} 

inline void get_coef_zams_all(const double &z,  double Lcoef[7],  double Rcoef[9]){
  double zlog = log10(z*50); // = z / z_sun,  z_sun = 0.02
  for(int i=0; i<7; i++){
    Lcoef[i] = (((L[i][4]*zlog + L[i][3]) * zlog + L[i][2]) * zlog + L[i][1]) * zlog + L[i][0];
  }
  for(int i=0; i<9; i++){
    Rcoef[i] = (((R[i][4]*zlog + R[i][3]) * zlog + R[i][2]) * zlog + R[i][1]) * zlog + R[i][0];
  }
} 

inline double get_luminosity_zams(const double &m, const double Lcoef[7]){
  double msqrt = sqrt(m);
  double m2 = m * m;  
  double m3 = m * m2;  
  double m5 = m3 * m2;
  double m5_5 = m5 * msqrt;
  double m7 = m5 * m2;
  double m8 = m5 * m3;
  double m9_5 = m5_5 * m2 * m2; 
  double m11 = m3 * m8;
  return (Lcoef[0]*m5_5 + Lcoef[1]*m11) / (Lcoef[2] + m3 + Lcoef[3]*m5 + Lcoef[4]*m7 + Lcoef[5]*m8 + Lcoef[6]*m9_5);
}

inline double get_radius_zams(const double &m, const double Rcoef[9]){
  double msqrt = sqrt(m);
  double m2 = m * m;  
  double m2_5 = m2 * msqrt;  
  double m4 = m2 * m2;  
  double m6 = m4 * m2;
  double m6_5 = m6 * msqrt;
  double m8_5 = m6_5 * m2;
  double m11 = m8_5 * m2_5;
  double m18_5 = m8_5 * m6 * m4;
  double m19 = m18_5 * msqrt;
  double m19_5 = m19 * msqrt;
  return (Rcoef[0]*m2_5 + Rcoef[1]*m6_5 + Rcoef[2]*m11 + Rcoef[3]*m19 + Rcoef[4]*m19_5) / (Rcoef[5] + Rcoef[6]*m2 + Rcoef[7]*m8_5 + m18_5 + Rcoef[8]*m19_5);
}


#endif //STELLAR_EVOLUTION_H
