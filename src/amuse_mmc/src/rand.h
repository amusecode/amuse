#if !defined __RAND_H
#define __RAND_H

#define INV_RAND_MAX2 4.6566130595601737e-10
#define RAND_MAX1     2147483646
#define RAND_MAX2     2147483562

/* two random number generators according to NRC */
int rand1();

/* rand2() is a bit slower, but it has longer period */
int rand2();

/* generate unit vector with random orientation */
void isorand(double *r);

void init_rand(int newseed);

void init_rand2(int newseed1, int newseed2);

int get_rand_seed1();

int get_rand_seed2();

#endif
