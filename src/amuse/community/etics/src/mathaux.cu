#include "common.hpp"
// #include "mathaux.hpp"

double FactorialSCF_tmpname(int x) {
    return (x <= 1) ? 1 : x*FactorialSCF_tmpname(x - 1);
}

double GammaPlusHalf(int x) {
    return pow(0.25, x) * sqrt(M_PI) * FactorialSCF_tmpname(2*x) / FactorialSCF_tmpname(x);
}

void RadialCoefficients(Real *Result) {
    for (int n = 0; n <= NMAX; n++) {
        for (int l = 0; l <= LMAX; l++) {
            double AnlTilde =  - pow(2.0, 8*l + 6) * FactorialSCF_tmpname(n) * (n + 2*l + 1.5);
            AnlTilde *= GammaPlusHalf(2*l + 1) * GammaPlusHalf(2*l + 1);
            AnlTilde /= (4*M_PI * (n*(n + 3 + 4*l)/2 + (l+1)*(2*l+1)) * FactorialSCF_tmpname(n+4*l+2));
            Result[(LMAX+1)*n+l] = (Real)AnlTilde;
        }
    }
}

__host__ __device__ Real Pl(int l, Real x) {
    // No range check or anything... The first 10 polynomials are hardcoded
    switch (l) {
        case 0  : return 1;
        case 1  : return x;
        case 2  : return (1.5*x*x - 0.5);
        case 3  : return (2.5*x*x*x - 1.5*x);
        case 4  : {Real x2=x*x; return (4.375*x2*x2  - 3.75*x2 + 0.375);}
        case 5  : {Real x2=x*x, x3=x2*x; return (7.875*x3*x2 - 8.75*x3 + 1.875*x);}
        case 6  : {Real x2=x*x, x4=x2*x2; return (14.4375*x4*x2 - 19.6875*x4 + 6.5625*x2 - 0.3125);}
        case 7  : {Real x2=x*x, x3=x2*x, x5=x3*x2; return (26.8125*x5*x2 - 43.3125*x5 + 19.6875*x3 - 2.1875*x);}
        case 8  : {Real x2=x*x, x4=x2*x2, x6=x4*x2, x8=x4*x4; return (50.2734375*x8 - 93.84375*x6 + 54.140625*x4 - 9.84375*x2 + 0.2734375);}
        case 9  : {Real x2=x*x, x3=x2*x, x5=x3*x2, x7=x5*x2, x9=x7*x2;
                   return (94.9609375*x9 - 201.09375*x7 + 140.765625*x5 - 36.09375*x3 + 2.4609375*x);}
        case 10 : {Real x2=x*x, x4=x2*x2, x6=x4*x2, x8=x4*x4, x10=x6*x4;
                   return (180.42578125*x10 - 427.32421875*x8 + 351.9140625*x6 - 117.3046875*x4 + 13.53515625*x2 - 0.24609375);}
        case 11 : {Real x2=x*x, x3=x2*x, x5=x3*x2, x7=x5*x2, x9=x7*x2, x11=x9*x2;
                   return (344.44921875*x11 - 902.12890625*x9 + 854.6484375*x7 - 351.9140625*x5 + 58.65234375*x3 - 2.70703125*x);}
        case 12 : {Real x2=x*x, x4=x2*x2, x6=x4*x2, x8=x4*x4, x10=x6*x4, x12=x6*x6;
                   return (660.1943359375*x12 - 1894.470703125*x10 + 2029.7900390625*x8 - 997.08984375*x6 + 219.9462890625*x4 - 17.595703125*x2 + 0.2255859375);}
        case 13 : {Real x2=x*x, x3=x2*x, x5=x3*x2, x7=x5*x2, x9=x7*x2, x11=x9*x2, x13=x9*x2;
                   return (1269.6044921875*x13 - 3961.166015625*x11 + 4736.1767578125*x9 - 2706.38671875*x7 + 747.8173828125*x5 -87.978515625*x3 + 2.9326171875*x);}
        case 14 : {Real x2=x*x, x4=x2*x2, x6=x4*x2, x8=x4*x4, x10=x6*x4, x12=x6*x6, x14=x8*x6;
                   return (2448.52294921875*x14 - 8252.42919921875*x12 + 10893.2065429688*x10 - 7104.26513671875*x8 + 2368.08837890625*x6 - 373.90869140625*x4 + 21.99462890625*x2 - 0.20947265625);}
        default : return -1e30;
    }
}

void AngularCoefficients(Real *Result) {
    for (int l = 0; l <= LMAX; l++) {
        for (int m = 0; m <= l; m++) {
            Result[(l+1)*l/2+m] = (sqrt((2*l+1)/(4*M_PI)) * sqrt(FactorialSCF_tmpname(l-m)/FactorialSCF_tmpname(l+m)));
        }
    }
}
