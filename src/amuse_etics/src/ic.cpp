/**
 * @file    ic.cpp
 * @author  Yohai Meiron <ymeiron@pku.edu.cn>
 * @brief   Functions to generate initial conditions.
 *
 * These are functions to generate a Hernquist or a Plummer sphere. In both
 * cases the model has a total mass of 1 and the total energy is -0.25 Henon
 * energy units. The Hernquist model generator was adapted from a program by
 * John Dubinski (1999, original version) from the NEMO toolbox; the accessory
 * functions are from Numerical Recipes in C (Press et al. 1992). The Plummer
 * generator is adapted from the gen-plum program by Peter Berczik.
 */

#include <iostream>
#include <cstdlib>
#include <cstdio>
#include <cmath>
#include "common.hpp"
#include "ic.hpp"

// Internal functions
namespace etics {
    namespace ic {
        double E;

        void polint(double xa[], double ya[], int n, double x, double *y, double *dy);
        double trapzd(double (*func)(double), double a, double b, int n);
        double qromb(double (*func)(double), double a, double b);
        double xrandom(double a, double b);
        double xrandom();
        double f(double E2), drho2(double u);
    }
}

double etics::ic::trapzd(double (*func)(double), double a, double b, int n) {
    double x, tnm, sum, del;
    static double s;
    int it, j;

    if (n == 1) {
        return (s=0.5*(b-a)*((*func)(a) + (*func)(b)));
    } else {
        for (it = 1, j = 1; j < n-1; j++) it <<= 1;
        tnm = it;
        del = (b-a)/tnm;
        x = a + 0.5*del;
        for (sum = 0, j = 1; j <= it; j++, x += del) sum += (*func)(x);
        s = 0.5*(s+(b-a)*sum/tnm);
        return s;
    }
}

void etics::ic::polint(double xa[], double ya[], int n, double x, double *y, double *dy) {
    int ns = 1;
    double den, dif, dift, ho, hp, w;
    double *c, *d;

    dif = fabs(x-xa[1]);
    c = new double[n+1];
    d = new double[n+1];
    for (int i = 1; i <= n; i++) {
        if ((dift=fabs(x-xa[i])) < dif) {
            ns = i;
            dif = dift;
        }
        c[i] = ya[i];
        d[i] = ya[i];
    }
    *y = ya[ns--];
    for (int m = 1; m < n; m++) {
        for (int i = 1; i <= n-m; i++) {
            ho = xa[i]-x;
            hp = xa[i+m]-x;
            w = c[i+1]-d[i];
            if ((den=ho-hp) == 0.0) {
                std::cerr << "Error in routine polint" << std::endl;
                exit(1);
            }
            den = w/den;
            d[i] = hp*den;
            c[i] = ho*den;
        }
        *y += (*dy=(2*ns < (n-m) ? c[ns+1] : d[ns--]));
    }
    delete[] c;
    delete[] d;
}

#define EPS 1.0e-6
#define JMAX 20
#define JMAXP (JMAX+1)
#define K 5
double etics::ic::qromb(double (*func)(double), double a, double b) {
    double ss, dss;
    double s[JMAXP+1], h[JMAXP+1];

    h[1] = 1.0;
    for (int j = 1; j <= JMAX; j++) {
        s[j] = trapzd(func, a, b, j);
        if (j >= K) {
            polint(&h[j-K], &s[j-K], K, 0, &ss, &dss);
            if (fabs(dss) < EPS*fabs(ss)) return ss;
        }
        s[j+1] = s[j];
        h[j+1] = 0.25*h[j];
    }
    std::cerr << "Too many steps in routine qromb" << std::endl;
    exit(1);
    return 0;
}
#undef EPS
#undef JMAX
#undef JMAXP
#undef K

double etics::ic::xrandom(double a, double b) {
    double Result = ((double)rand())/((double)RAND_MAX);
    Result *= (b - a);
    Result += a;
    return Result;
}

double etics::ic::xrandom() {
    return xrandom(0, 1);
}

/* Distribution function of a Hernquist model */
double etics::ic::f(double E2) {
    E = E2;
    return qromb(drho2, 0.0, 2.0*sqrt(E))*(sqrt(32.)*(M_PI*M_PI*M_PI));
}

double etics::ic::drho2(double u) {
    double p = E - 0.25*u*u;
    double r;
    if (p <= 0) return 0;
    else {
        double p1 = 1/p;
        double b0 = 0.5*(2-p1);
        double c0 = (1-p1);
        r = -b0 + sqrt(b0*b0 - c0);
    }

    double r2 = r*r, r3 = r2*r;

    double c0 = 1.0 + r;
    double c02 = c0*c0;
    double c03 = c02*c0;
    double c04 = c02*c02;
    double dp1 = -1/c02;
    double dp2 = 2/c03;

    double drho1r = -(1+4*r)/(r2*c04);
    double drho2r =  2.0*(14*r2 + 5*r + 1)/(r3*c04*c0);

    return drho2r/(dp1*dp1) - drho1r/(dp1*dp1*dp1)*dp2;
}

void etics::ic::hernquist(int N, int Seed, Particle **ParticleList) {
    using namespace etics::ic;
    double x, y, z, vx, vy, vz;
    double fmax, f0, f1, v2, vmax, vmax2, Ep=0, Ek=0;

    srand(Seed);

    *ParticleList = new Particle[N];

    for(int i = 0; i < N; i++) {
        double eta = xrandom(0.0, 1.0);
        double radius = sqrt(eta)/(1-sqrt(eta));
        double phi = xrandom(0.0, 2*M_PI);
        double cth = xrandom(-1.0, 1.0);
        double sth = sqrt(1.0 - cth*cth);
        x = radius*sth*cos(phi);
        y = radius*sth*sin(phi);
        z = radius*cth;
        double psi0 = 1/(1 + radius);
        Ep += -psi0;
        vmax2 = 2.0*psi0;
        vmax = sqrt(vmax2);
        fmax = f(psi0);
        f0 = fmax; f1 = 1.1*fmax; // just some initial values
        while(f1 > f0) {
            // pick a velocity
            v2 = 2.0*vmax2;
            while(v2 > vmax2) { // rejection technique
                vx = vmax*xrandom(-1.0, 1.0);
                vy = vmax*xrandom(-1.0, 1.0);
                vz = vmax*xrandom(-1.0, 1.0);
                v2 = vx*vx + vy*vy + vz*vz;
            }
            f0 = f(psi0 - 0.5*v2);
            f1 = fmax*xrandom(0.0, 1.0);
        }
        Ek += 0.5*v2;
        Particle p;
        p.m = 1.0/N;
        p.pos = vec3(x, y, z);
        p.vel = vec3(vx, vy, vz);
        p.ID = i;
        p.Status = 0;
        p.CalculateR2();
        (*ParticleList)[i] = p;
    }
    Ep /= 2*N; Ek /= N;
    double E = Ep + Ek;
    double RadiusFactor = -4*E; // should be close to 1/3
    double VelocityFactor = 1/sqrt(-4*E); // should be close to sqrt(3)
    for(int i = 0; i < N; i++) {
        (*ParticleList)[i].pos = (*ParticleList)[i].pos * RadiusFactor;
        (*ParticleList)[i].vel = (*ParticleList)[i].vel * VelocityFactor;
    }
#warning you must centralize yourself!!!
}

#define RMAX 10.0
void etics::ic::plummer(int N, int Seed, Particle **ParticleList) {
    double X,  Y,  Z, Vx, Vy, Vz,
           X1, X2, X3, X4, X5, X6, X7,
           R, Ve, V;

    srand(Seed);

    *ParticleList = new Particle[N];

    int i = 0;
    while(i < N) {
        X1 = xrandom(); X2 = xrandom(); X3 = xrandom();
        R = 1.0/sqrt( (pow(X1,-2.0/3.0) - 1.0) );
        if (R > RMAX) continue;
        Z = (1.0 - 2.0*X2)*R;
        X = sqrt(R*R - Z*Z) * cos(2.0*M_PI*X3);
        Y = sqrt(R*R - Z*Z) * sin(2.0*M_PI*X3);

        Ve = sqrt(2.0)*pow( (1.0 + R*R), -0.25 );

        X4 = 0.0;
        X5 = 0.0;
        while( 0.1*X5 >= X4*X4*pow( (1.0-X4*X4), 3.5) ) {
            X4 = xrandom();
            X5 = xrandom();
        }

        V = Ve*X4;

        X6 = xrandom();
        X7 = xrandom();

        Vz = (1.0 - 2.0*X6)*V;
        Vx = sqrt(V*V - Vz*Vz) * cos(2.0*M_PI*X7);
        Vy = sqrt(V*V - Vz*Vz) * sin(2.0*M_PI*X7);

#define CONV (3.0*M_PI/16.0)
        X *= CONV; Y *= CONV; Z *= CONV;
        Vx /= sqrt(CONV); Vy /= sqrt(CONV); Vz /= sqrt(CONV);
#undef CONV

        Particle p;
        p.m = 1.0/N;
        p.pos = vec3(X, Y, Z);
        p.vel = vec3(Vx, Vy, Vz);
        p.ID = i;
        p.Status = 0;
        p.CalculateR2();
        (*ParticleList)[i] = p;

        i++;
    } /* while(i < N) */
}
#undef RMAX
