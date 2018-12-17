
/* Leapfrog integrators */
void InitialKick(PS::ParticleSystem<FP_nbody> & psys, const PS::F64 dt){
    for (PS::S32 i = 0; i < psys.getNumberOfParticleLocal(); i++){
        psys[i].vel += 0.5 * dt * psys[i].acc;
    }
}
void InitialKick(PS::ParticleSystem<FP_sph> & psys, const PS::F64 dt){
    for (PS::S32 i = 0; i < psys.getNumberOfParticleLocal(); i++){
        psys[i].vel_half = psys[i].vel + 0.5 * dt * (psys[i].acc_grav + psys[i].acc_hydro);
#if !defined(ISOTHERMAL_EOS)
        psys[i].eng_half = psys[i].eng + 0.5 * dt * psys[i].eng_dot;
        psys[i].ent_half = psys[i].ent + 0.5 * dt * psys[i].ent_dot;
#endif
    }
}


void FullDrift(PS::ParticleSystem<FP_nbody>& psys, const PS::F64 dt){
    for (PS::S32 i = 0; i < psys.getNumberOfParticleLocal(); i++){
        psys[i].pos += dt * psys[i].vel;
    }
}
void FullDrift(PS::ParticleSystem<FP_sph>& psys, const PS::F64 dt){
    for (PS::S32 i = 0; i < psys.getNumberOfParticleLocal(); i++){
        psys[i].pos += dt * psys[i].vel_half;
    }
}


void Predict(PS::ParticleSystem<FP_sph>& psys, const PS::F64 dt){
    for(PS::S32 i = 0; i < psys.getNumberOfParticleLocal(); i++){
        psys[i].vel += dt * (psys[i].acc_grav + psys[i].acc_hydro);
#if !defined(ISOTHERMAL_EOS)
        psys[i].eng += dt * psys[i].eng_dot;
        psys[i].ent += dt * psys[i].ent_dot;
#endif
    }
}


void FinalKick(PS::ParticleSystem<FP_nbody>& psys, const PS::F64 dt){
    for (PS::S32 i = 0; i < psys.getNumberOfParticleLocal(); i++){
        psys[i].vel += 0.5 * dt * psys[i].acc;
    }
}
void FinalKick(PS::ParticleSystem<FP_sph>& psys, const PS::F64 dt){
    for (PS::S32 i = 0; i < psys.getNumberOfParticleLocal(); i++){
        psys[i].vel = psys[i].vel_half + 0.5 * dt * (psys[i].acc_grav + psys[i].acc_hydro);
#if !defined(ISOTHERMAL_EOS)
        psys[i].eng = psys[i].eng_half + 0.5 * dt * psys[i].eng_dot;
        psys[i].ent = psys[i].ent_half + 0.5 * dt * psys[i].ent_dot;
#endif
    }
#if defined(DUMP_VELOCITY_OF_SPH_PARTICLE)
    const PS::F64 coeff = 0.1;
    for (PS::S32 i = 0; i < psys.getNumberOfParticleLocal(); i++) {
        psys[i].vel *= std::exp(- coeff * (CFL_hydro/0.1) * psys[i].snds * dt / psys[i].smth);
    }
#endif
}
