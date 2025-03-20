#ifndef GALACTICS_RIGIDPOT_H
#define GALACTICS_RIGIDPOT_H

void InitRigidPotential(char * potname);

void getforce(float x, float y, float z, float *ax, float *ay, float *az, float *pot);

#endif

