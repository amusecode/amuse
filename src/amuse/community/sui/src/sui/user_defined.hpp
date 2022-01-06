// Include the standard C++ headers
#include <cmath>
#include <numbers>
#include <math.h>
#include <cfloat>
#include <cstdio>
#include <cassert>
#include <cstdint>
#include <iostream>
#include <fstream>
#include <sstream>
#include <array>
#include <vector>
#include <sys/stat.h>
#include <time.h>
// Include the header file of FDPS
#include <particle_simulator.hpp>
// Include the header file of Phantom-GRAPE library
#if defined(ENABLE_PHANTOM_GRAPE_X86)
#include <gp5util.h>
#endif

//#include "mathematical_constants.h"
//#include "physical_constants.h"


//const PS::F64 specific_heat_ratio = 5.0/3.0;  // FIXME needs to be AMUSE settable
const PS::F64 specific_heat_ratio = 1.0; //5.0/3.0;  // FIXME needs to be AMUSE settable
PS::F64 SCF_smth = 1.0; // scale factor for smoothing length
PS::F64 mass_avg;
const PS::F64 CFL_hydro = 0.3; // CFL(Courant-Friedrichs-Lewy) number for hydrodynamics
const PS::F64 alpha_AV = 1.0; // coefficient of artificial viscosity
const PS::S32 N_neighbor = 50; // number of neighbor particles
const PS::F64 CFL_dyn = 0.3; // coefficient used to limit a timestep in terms of dynamics
PS::F64 eps_grav = 0.01;

/* Kernel Function */
PS::F64 W(const PS::F64 r, const PS::F64 h){
    // M4 Cubic spline kernel
    // (see Eq. (4) in Springel (2005)[MNRAS,364,1105])
    const PS::F64 u = r/h;
    const PS::F64 cc=8.0/(std::numbers::pi*h*h*h);
    if (u <= 0.5) {
        const PS::F64 u2 = u*u;
        return cc*(1.0+u2*6.0*(-1.0+u));
    }
    else if ((0.5 < u) && (u <= 1.0)) {
        const PS::F64 s = 1.0-u;
        return cc*2.0*s*s*s;
    }
    else {
        return 0.0;
    }
}

PS::F64vec gradW(const PS::F64vec dr, const PS::F64 h){
    // This subroutine gives \nabla W(r,h), i.e.,
    // \dfrac{\partial W(r,h)}{\partial r}\dfrac{dr}{r}.
    const PS::F64 r = std::sqrt(dr * dr);
    const PS::F64 u=r/h;
    const PS::F64 cc = -48.0/(std::numbers::pi*h*h*h*h);
#if defined(USE_PRESCR_OF_THOMAS_COUCHMAN_1992)
    if (u <= 1.0/3.0) {
        return dr * cc*(1.0/3.0)/(r);
    }
    else if ((1.0/3.0 < u) && (u <= 0.5)) {
        return dr * cc*u*(2.0-3.0*u)/(r);
    }
    else if ((0.5 < u) && (u < 1.0)) {
        return dr * cc*(1.0-u)*(1.0-u)/(r);
    }
    else {
        return PS::F64vec(0.0, 0.0, 0.0);
    }
#else
    if ((0.0 < u) && (u <= 0.5)) {
        return dr * cc*u*(2.0-3.0*u)/(r);
    }
    else if ((0.5 < u) && (u < 1.0)) {
        return dr * cc*(1.0-u)*(1.0-u)/(r);
    }
    else {
        // r=0 case is included in this branch
        return PS::F64vec(0.0, 0.0, 0.0);
    }
#endif
}


PS::F64 dWdh(const PS::F64 r, const PS::F64 h){
   // This subroutine gives dW(r,h)/dh, i.e.,
   // \dfrac{\partial W(r,h)}{\partial h}.
   const PS::F64 u=r/h;
   const PS::F64 cc=-24.0/(std::numbers::pi*h*h*h*h);
   if (u <= 0.5) {
      const PS::F64 u2 = u*u;
      return cc*(1.0
                +u2*(-10.0
                     +12.0*u));
   }
   else if ((0.5 < u) && (u < 1.0)) {
      const PS::F64 s = 1.0-u;
      return cc*2.0*s*s*(1.0-2.0*u);
   }
   else {
      return 0.0;
   }

}

//** Force class for gravity calculation
class Force_grav {
public:
    PS::F64vec acc;
    PS::F64 pot; 
    void clear() {
        acc = 0.0;
        pot = 0.0;
    }
};
//** Force classes for SPH calculation
class Force_dens{
public:
   PS::S32 flag;
   PS::F64 dens;
   PS::F64 smth;
   PS::F64 gradh;
   PS::F64 divv;
   PS::F64vec rotv;
   void clear(){
      flag  = 0;
      dens  = 0.0;
      gradh = 0.0;
      divv  = 0.0;
      rotv  = 0.0;
   }
};
class Force_hydro{
public:
   PS::F64vec acc;
   PS::F64 eng_dot;
   PS::F64 ent_dot;
   PS::F64 dt;
   void clear(){
      acc = 0.0;
      eng_dot = 0.0;
      ent_dot = 0.0;
   }
};

class FP_test {
public:
    PS::S64    id;
    PS::F64vec pos;
    PS::F64vec vel;
    PS::F64vec acc_grav; // gravitational acceleration
    PS::F64    pot_grav; // gravitational potential
    PS::F64vec acc_hydro; // acceleration due to pressure-gradient

    PS::F64    dens; // mass density
    PS::F64    eng; // specific internal energy
    //PS::F64    pres; // pressure
    PS::F64    smth; // smoothing length
    PS::F64    gradh; // grad-h term
    PS::F64    divv; // divergence of velocity
    PS::F64vec rotv; // rotation of velocity
    //PS::F64    snds; // sound speed

    PS::F64vec getPos() const {
        return pos;
    }

    void setPos(const PS::F64vec& pos){
        this->pos = pos;
    }

    PS::F64 getCharge() const {
        return 0;
    }

    void copyFromForce(const Force_grav& f) {
        this->acc_grav = f.acc;
        this->pot_grav = f.pot;
    }

    void copyFromForce(const Force_dens& f){
        //this->flag  = f.flag;
        this->dens  = f.dens;
        this->smth  = f.smth;
        this->gradh = f.gradh;
        this->divv  = f.divv;
        this->rotv  = f.rotv;
        
    }

    void copyFromForce(const Force_hydro& f){
        this->acc_hydro = f.acc;
        //this->eng_dot   = f.eng_dot;
        //this->ent_dot   = f.ent_dot;
        //this->dt        = f.dt;
    }

};

class FP_sph {
public:
    PS::S64    id;
    PS::F64    mass;
    PS::F64vec pos;
    PS::F64vec vel;
    PS::F64vec acc_grav; // gravitational acceleration
    PS::F64    pot_grav; // gravitational potential
    PS::F64vec acc_hydro; // acceleration due to pressure-gradient
    PS::S32    flag;
    PS::F64    dens; // mass density
    PS::F64    eng; // specific internal energy
    PS::F64    ent; // entropy
    PS::F64    pres; // pressure
    PS::F64    smth; // smoothing length
    PS::F64    gradh; // grad-h term
    PS::F64    divv; // divergence of velocity
    PS::F64vec rotv; // rotation of velocity
    PS::F64    BalSW; // Balsara switch
    PS::F64    snds; // sound speed
    PS::F64    eng_dot; // time rate of change of `eng`
    PS::F64    ent_dot; // time rate of change of `ent`
    PS::F64    dt; // hydrodynamic time step for this particle
    PS::F64vec vel_half;
    PS::F64    eng_half;
    PS::F64    ent_half;

    void copyFromForce(const Force_grav& f) {
        this->acc_grav = f.acc;
        this->pot_grav = f.pot;
    }
    void copyFromForce(const Force_dens& f){
        this->flag  = f.flag;
        this->dens  = f.dens;
        this->smth  = f.smth;
        this->gradh = f.gradh;
        this->divv  = f.divv;
        this->rotv  = f.rotv;
        
    }
    void copyFromForce(const Force_hydro& f){
        this->acc_hydro = f.acc;
        this->eng_dot   = f.eng_dot;
        this->ent_dot   = f.ent_dot;
        this->dt        = f.dt;
    }
    PS::F64 getCharge() const{
        return this->mass;
    }
    PS::F64vec getPos() const{
        return this->pos;
    }
    PS::F64 getRSearch() const{
        return this->smth;
    }
    void setPos(const PS::F64vec& pos){
       this->pos = pos;
    }
    void writeAscii(FILE* fp) const{
        fprintf(fp,
                "%lld\t%lf\t%lf\t%lf\t%lf\t%lf\t"
                "%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n",
                this->id, this->mass,
                this->pos.x, this->pos.y, this->pos.z,
                this->vel.x, this->vel.y, this->vel.z,
                this->dens, this->eng, this->ent, this->pres);
    }
    void readAscii(FILE* fp){
        fscanf(fp,
               "%lld\t%lf\t%lf\t%lf\t%lf\t%lf\t"
               "%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n",
               &this->id, &this->mass,
               &this->pos.x, &this->pos.y, &this->pos.z,
               &this->vel.x, &this->vel.y, &this->vel.z,
               &this->dens, &this->eng, &this->ent, &this->pres);
    }
    void writeBinaryPos(FILE* fp) const {
        fwrite(&this->pos, sizeof(this->pos), 1, fp);
    }
    void setEntropy(){
        ent = (specific_heat_ratio - 1.0) * eng / std::pow(dens, specific_heat_ratio - 1.0);
    }
    void setPressure(){
#if defined(ISOTHERMAL_EOS)
        // In this case, eng = const.
        pres = (specific_heat_ratio - 1.0) * dens * eng;
        ent  = pres / std::pow(dens, specific_heat_ratio);
#else
#if defined(USE_ENTROPY)
        pres = ent * std::pow(dens, specific_heat_ratio);
        eng  = pres / ((specific_heat_ratio - 1.0) * dens);
#else
        pres = (specific_heat_ratio - 1.0) * dens * eng;
        ent  = pres / std::pow(dens, specific_heat_ratio);
#endif
#endif
        snds = std::sqrt(specific_heat_ratio * pres / dens);
#if defined(USE_BALSARA_SWITCH)
        BalSW = std::fabs(divv) / (std::fabs(divv) + std::sqrt(rotv * rotv) + 1.0e-4 * snds / smth); 
#else
        BalSW = 1.0;
#endif
    }

};

//** Full Particle Classes
class FP_nbody{
public:
    PS::S64    id;
    PS::F64    mass;
    PS::F64vec pos;
    PS::F64vec vel;
    PS::F64vec acc;
    PS::F64    pot;

    PS::F64vec getPos() const {
        return pos;
    }
    PS::F64 getCharge() const {
        return mass;
    }
    void setPos(const PS::F64vec& pos){
       this->pos = pos;
    }
    void copyFromForce(const Force_grav & f) {
        this->acc = f.acc;
        this->pot = f.pot;
    }
    void writeAscii(FILE* fp) const {
        fprintf(fp, "%lld\t%g\t%g\t%g\t%g\t%g\t%g\t%g\n", 
                this->id, this->mass,
                this->pos.x, this->pos.y, this->pos.z,
                this->vel.x, this->vel.y, this->vel.z);
    }
    void readAscii(FILE* fp) {
        fscanf(fp, "%lld\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n", 
               &this->id, &this->mass,
               &this->pos.x, &this->pos.y, &this->pos.z,
               &this->vel.x, &this->vel.y, &this->vel.z);
    }
};

//** Essential Particle Class
class EP_grav {
public:
    PS::S64 id;
    PS::F64 mass;
    PS::F64vec pos;

    void copyFromFP(const FP_nbody& fp) {
        this->id   = fp.id;
        this->mass = fp.mass;
        this->pos  = fp.pos;
    }
    void copyFromFP(const FP_sph& fp) {
        this->id   = fp.id;
        this->mass = fp.mass;
        this->pos  = fp.pos;
    }
    PS::F64 getCharge() const {
        return this->mass;
    }
    PS::F64vec getPos() const {
        return this->pos;
    }
};

class EP_hydro {
public:
   PS::S64    id;
   PS::F64vec pos;
   PS::F64vec vel;
   PS::F64    mass;
   PS::F64    smth;
   PS::F64    dens;
   PS::F64    pres;
   PS::F64    gradh;
   PS::F64    snds;
   PS::F64    BalSW;

   void copyFromFP(const FP_sph& fp){
      this->id    = fp.id;
      this->pos   = fp.pos;
      this->vel   = fp.vel;
      this->mass  = fp.mass;
      this->smth  = fp.smth;
      this->dens  = fp.dens;
      this->pres  = fp.pres;
      this->gradh = fp.gradh;
      this->snds  = fp.snds;
      this->BalSW = fp.BalSW;
   }
   PS::F64vec getPos() const{
      return this->pos;
   }
   PS::F64 getRSearch() const{
      return SCF_smth * this->smth;
   }
   void setPos(const PS::F64vec& pos){
      this->pos = pos;
   }
};

class CalcDensity{
public:
    void operator () (const EP_hydro * ep_i,
                      const PS::S32 n_ip,
                      const EP_hydro * ep_j,
                      const PS::S32 n_jp,
                      Force_dens * force){
#if defined(ENABLE_VARIABLE_SMOOTHING_LENGTH)
        const PS::F64 eps = 1.0e-6;
        const PS::F64 M_trgt = mass_avg * N_neighbor;
        for (PS::S32 i = 0; i < n_ip; i++) {
            PS::F64 dens = 0.0;
            PS::F64 h = ep_i[i].smth;
            const PS::F64 h_max_alw = SCF_smth * h; // maximum allowance
            PS::F64 h_L = 0.0;
            PS::F64 h_U = h_max_alw;
            PS::F64 dh_prev = 0.0;
            PS::S32 n_unchanged = 0;
            // Software caches
            PS::F64 * mj  = (PS::F64 *)malloc(sizeof(PS::F64) * n_jp);
            PS::F64 * rij = (PS::F64 *)malloc(sizeof(PS::F64) * n_jp);
            for (PS::S32 j = 0; j < n_jp; j++) {
                mj[j] = ep_j[j].mass;
                const PS::F64vec dr = ep_j[j].pos - ep_i[i].pos;
                rij[j] = std::sqrt(dr * dr);
            }
            for (;;) {
                // Calculate density
                dens = 0.0;
                for (PS::S32 j = 0; j < n_jp; j++) {
                   dens += mj[j] * W(rij[j], h);
                }
                // Check if the current value of the smoohting length satisfies 
                // Eq.(5) in Springel (2005).
                const PS::F64 M = 4.0 * std::numbers::pi * h * h * h * dens / 3.0;
                if ((h < h_max_alw) && (std::abs(M/M_trgt - 1.0) < eps)) {
                    // In this case, Eq.(5) holds within a specified accuracy.
                    force[i].flag = 1;
                    force[i].dens = dens;
                    force[i].smth = h;
                    break;
                }
                if (((h == h_max_alw) && (M < M_trgt)) || (n_unchanged == 4)) {
                    // In this case, we skip this particle forcibly.
                    // In order to determine consistently the density
                    // and the smoohting length for this particle,
                    // we must re-perform calcForceAllAndWriteBack().
                    force[i].flag = 0;
                    force[i].dens = dens;
                    force[i].smth = h_max_alw;
                    break;
                }
                // Update h_L & h_U
                if (M < M_trgt) {
                   if (h_L < h) h_L = h;
                }
                else if (M_trgt < M) {
                   if (h < h_U) h_U = h;
                }
                const PS::F64 dh = h_U - h_L;
                if (dh == dh_prev) {
                   n_unchanged++;
                }
                else {
                   dh_prev = dh;
                   n_unchanged = 0;
                }
                // Update smoothing length
                h = std::pow((3.0 * M_trgt)/(4.0 * std::numbers::pi * dens), 1.0/3.0);
                if ((h <= h_L) || (h == h_U)) {
                   // In this case, we switch to the bisection search.
                   // The inclusion of '=' in the if statement is very
                   // important to escape a limit cycle.
                   h = 0.5 * (h_L + h_U);
                }
                else if (h_U < h) {
                   h = h_U;
                }
            }
            // Calculate grad-h term
            if (force[i].flag == 1) {
                PS::F64 drho_dh = 0.0;
                for (PS::S32 j = 0; j < n_jp; j++) {
                   drho_dh += mj[j] * dWdh(rij[j], h);
                }
                force[i].gradh = 1.0 / (1.0 + (h * drho_dh) / (3.0 * dens));
            } 
            else {
                force[i].gradh = 1.0; // dummy value
            }
#if defined(USE_BALSARA_SWITCH)
            // Compute \div v & \rot v for Balsara switch
            for (PS::S32 j = 0; j < n_jp; j++) {
               const PS::F64vec dr = ep_i[i].pos - ep_j[j].pos;
               const PS::F64vec dv = ep_i[i].vel - ep_j[j].vel;
               force[i].divv -= mj[j] * dv * gradW(dr, force[i].smth);
               force[i].rotv -= mj[j] * dv ^ gradW(dr, force[i].smth);
            }
            force[i].divv /= force[i].dens;
            force[i].rotv /= force[i].dens;
#endif
            // Release memory
            free(mj);
            free(rij);
        }
#else
        for (PS::S32 i = 0; i < n_ip ; i++){
            // Compute density
            for (PS::S32 j = 0; j < n_jp; j++){
               const PS::F64vec dr = ep_j[j].pos - ep_i[i].pos;
               const PS::F64 rij = std::sqrt(dr * dr);
               force[i].dens += ep_j[j].mass * W(rij, ep_i[i].smth);
            }
            force[i].smth = ep_i[i].smth;
            force[i].gradh = 1.0;
#if defined(USE_BALSARA_SWITCH)
            // Compute \div v & \rot v for Balsara switch
            for (PS::S32 j = 0; j < n_jp; j++) {
               const PS::F64vec dr = ep_i[i].pos - ep_j[j].pos;
               const PS::F64vec dv = ep_i[i].vel - ep_j[j].vel;
               force[i].divv -= ep_j[j].mass * dv * gradW(dr, force[i].smth);
               force[i].rotv -= ep_j[j].mass * dv ^ gradW(dr, force[i].smth);
            }
            force[i].divv /= force[i].dens;
            force[i].rotv /= force[i].dens;
#endif
        }
#endif
    }
};

class CalcHydroForce{
public:
    void operator () (const EP_hydro * ep_i,
                      const PS::S32 n_ip,
                      const EP_hydro * ep_j,
                      const PS::S32 n_jp,
                      Force_hydro * force){
        for (PS::S32 i = 0; i < n_ip; i++){
           const PS::F64vec pos_i = ep_i[i].pos;
           const PS::F64vec vel_i = ep_i[i].vel;
           const PS::F64 smth_i   = ep_i[i].smth;
           const PS::F64 dens_i   = ep_i[i].dens;
           const PS::F64 pres_i   = ep_i[i].pres;
           const PS::F64 f_i      = ep_i[i].gradh;
           const PS::F64 snds_i   = ep_i[i].snds;
           const PS::F64 povrho2_i = pres_i / (dens_i * dens_i);
           PS::F64 v_sig_max = 0.0;
           for (PS::S32 j = 0; j < n_jp; j++){
              const PS::F64vec dr = pos_i - ep_j[j].pos;
              const PS::F64vec dv = vel_i - ep_j[j].vel;
              const PS::F64 w_ij = (dv * dr < 0) ? dv * dr / std::sqrt(dr * dr) : 0;
              const PS::F64 v_sig = snds_i + ep_j[j].snds - 3.0 * w_ij;
              v_sig_max = std::max(v_sig_max, v_sig);
              const PS::F64 AV = - 0.5 * alpha_AV * v_sig * w_ij / (0.5 * (dens_i + ep_j[j].dens))
                                 * 0.5 * (ep_i[i].BalSW + ep_j[j].BalSW);
              const PS::F64vec gradW_i  = gradW(dr, smth_i);
              const PS::F64vec gradW_j  = gradW(dr, ep_j[j].smth);
              const PS::F64vec gradW_ij = 0.5 * (gradW_i + gradW_j);
              const PS::F64 povrho2_j = ep_j[j].pres / (ep_j[j].dens * ep_j[j].dens);
              const PS::F64 f_j = ep_j[j].gradh;
              force[i].acc     -= ep_j[j].mass * (f_i * povrho2_i * gradW_i
                                                 +f_j * povrho2_j * gradW_j
                                                 +AV * gradW_ij);
              force[i].eng_dot += ep_j[j].mass * (f_i * povrho2_i * gradW_i
                                                 +0.5 * AV * gradW_ij) * dv;
              force[i].ent_dot += 0.5 * ep_j[j].mass * AV * gradW_ij * dv;
           }
           const PS::F64 p = specific_heat_ratio - 1.0;
           force[i].ent_dot *= p/std::pow(dens_i, p);
           force[i].dt = CFL_hydro * 2.0 * ep_i[i].smth / v_sig_max;
        }
    }
};

/* Interaction functions */
#if defined(ENABLE_PHANTOM_GRAPE_X86)
template <class TParticleJ>
void CalcGravity(const EP_grav * iptcl,
                 const PS::S32 ni,
                 const TParticleJ * jptcl,
                 const PS::S32 nj,
                 Force_grav * force) {
    const PS::S32 nipipe = ni;
    const PS::S32 njpipe = nj;
    PS::F64 (*xi)[3] = (PS::F64 (*)[3])malloc(sizeof(PS::F64) * nipipe * PS::DIMENSION);
    PS::F64 (*ai)[3] = (PS::F64 (*)[3])malloc(sizeof(PS::F64) * nipipe * PS::DIMENSION);
    PS::F64  *pi     = (PS::F64  *    )malloc(sizeof(PS::F64) * nipipe);
    PS::F64 (*xj)[3] = (PS::F64 (*)[3])malloc(sizeof(PS::F64) * njpipe * PS::DIMENSION);
    PS::F64  *mj     = (PS::F64  *    )malloc(sizeof(PS::F64) * njpipe);
    for(PS::S32 i = 0; i < ni; i++) {
        xi[i][0] = iptcl[i].getPos()[0];
        xi[i][1] = iptcl[i].getPos()[1];
        xi[i][2] = iptcl[i].getPos()[2];
        ai[i][0] = 0.0;
        ai[i][1] = 0.0;
        ai[i][2] = 0.0;
        pi[i]    = 0.0;
    }
    for(PS::S32 j = 0; j < nj; j++) {
        xj[j][0] = jptcl[j].getPos()[0];
        xj[j][1] = jptcl[j].getPos()[1];
        xj[j][2] = jptcl[j].getPos()[2];
        mj[j]    = jptcl[j].getCharge();
    }
    PS::S32 devid = PS::Comm::getThreadNum();
    g5_set_xmjMC(devid, 0, nj, xj, mj);
    g5_set_nMC(devid, nj);
    g5_calculate_force_on_xMC(devid, xi, ai, pi, ni);
    for(PS::S32 i = 0; i < ni; i++) {
        force[i].acc[0] += ai[i][0];
        force[i].acc[1] += ai[i][1];
        force[i].acc[2] += ai[i][2];
        force[i].pot    -= pi[i];
    }
    free(xi);
    free(ai);
    free(pi);
    free(xj);
    free(mj);
}
#else
template <class TParticleJ>
void CalcGravity (const EP_grav * ep_i,
                  const PS::S32 n_ip,
                  const TParticleJ * ep_j,
                  const PS::S32 n_jp,
                  Force_grav * force) {
    const PS::F64 eps2 = eps_grav * eps_grav;
    for(PS::S32 i = 0; i < n_ip; i++){
        PS::F64vec xi = ep_i[i].getPos();
        PS::F64vec ai = 0.0;
        PS::F64 poti = 0.0;
        for(PS::S32 j = 0; j < n_jp; j++){
            PS::F64vec rij    = xi - ep_j[j].getPos();
            PS::F64    r3_inv = rij * rij + eps2;
            PS::F64    r_inv  = 1.0/sqrt(r3_inv);
            r3_inv  = r_inv * r_inv;
            r_inv  *= ep_j[j].getCharge();
            r3_inv *= r_inv;
            ai     -= r3_inv * rij;
            poti   -= r_inv;
        }
        force[i].acc += ai;
        force[i].pot += poti;
    }
}
#endif
