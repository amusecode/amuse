struct particle{
	double x,y,z;
	double vx,vy,vz;
}particle;

void operator_sei(double dt, struct particle *p);
void operator_H0(double dt, struct particle* p);
void operator_phi(double dt, struct particle* p);
void operator_HKin(double dt, struct particle* p);
void operator_kick(double dt, struct particle* p);
void operator_lf(double dt, struct particle* p);
void operator_quinn(double dt, struct particle* p);
