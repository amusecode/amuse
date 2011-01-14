#include <math.h>
#include <rand.h>

const int _RAND_M0=2147483647, _RAND_A0=16807, _RAND_Q0=127773, _RAND_R0=2836;
const int _RAND_M1=2147483563, _RAND_A1=40014, _RAND_Q1=53668, _RAND_R1=12211;
const int _RAND_M2=2147483399, _RAND_A2=40692, _RAND_Q2=52774, _RAND_R2=3791;
int _RAND_SEED0=1, _RAND_SEED1=1, _RAND_SEED2=1, _ISET=0;

int rand1()
{
	int k = _RAND_SEED0 / _RAND_Q0;
	_RAND_SEED0 = _RAND_A0 * (_RAND_SEED0 - k * _RAND_Q0) - _RAND_R0 * k;
	if (_RAND_SEED0 < 0) _RAND_SEED0 += _RAND_M0;
	return _RAND_SEED0;
};

int rand2()
{
	int k = _RAND_SEED1 / _RAND_Q1;
	_RAND_SEED1 = _RAND_A1 * (_RAND_SEED1 - k * _RAND_Q1) - _RAND_R1 * k;
	if (_RAND_SEED1 < 0) _RAND_SEED1 += _RAND_M1;
	k = _RAND_SEED2 / _RAND_Q2;
	_RAND_SEED2 = _RAND_A2 * (_RAND_SEED2 - k * _RAND_Q2) - _RAND_R2 * k;
	if (_RAND_SEED2 < 0) _RAND_SEED2 += _RAND_M2;
	k = _RAND_SEED1 - _RAND_SEED2;
	if (k < 1) k += RAND_MAX2;
	return k;
};

/* generate unit vector with random orientation */
void isorand(double *r)
{
	double rad;

	do
	 {
		r[0] = 2. * rand2() * INV_RAND_MAX2 - 1.;
		r[1] = 2. * rand2() * INV_RAND_MAX2 - 1.;
		r[2] = 2. * rand2() * INV_RAND_MAX2 - 1.;
		rad = r[0]*r[0] + r[1]*r[1] + r[2]*r[2];
	 }
	while ((rad > 1.) || (rad < 0.0001));

	rad = sqrt(rad);
	r[0] /= rad;
	r[1] /= rad;
	r[2] /= rad;
};

void init_rand(int newseed)
{
	if ((newseed <1) || (newseed > RAND_MAX1)) _RAND_SEED0 = 1;
	else _RAND_SEED0 = newseed;
	if ((newseed <1) || (newseed > RAND_MAX2)) _RAND_SEED1 = 1;
	else _RAND_SEED1 = newseed;
};

void init_rand2(int newseed1, int newseed2)
{
	if ((newseed1 <1) || (newseed1 > RAND_MAX1)) _RAND_SEED0 = 1;
	else _RAND_SEED0 = newseed1;
	if ((newseed1 <1) || (newseed1 > RAND_MAX2)) _RAND_SEED1 = 1;
	else _RAND_SEED1 = newseed1;
	if ((newseed2 <1) || (newseed2 > RAND_MAX2)) _RAND_SEED2 = 1;
	else _RAND_SEED2 = newseed2;
};

int get_rand_seed1()
{
	return _RAND_SEED1;
};

int get_rand_seed2()
{
	return _RAND_SEED2;
};
