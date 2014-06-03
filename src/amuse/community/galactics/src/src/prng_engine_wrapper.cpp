#include <random>
#include <functional>

#include "prng_engine.hpp"

#define SEED  1234321
#define SKIP  100000

sitmo::prng_engine & global_urng()
{
  static sitmo::prng_engine u;
  return u;
}

extern "C" void ran_seed(long long i)
{
  global_urng().seed(SEED);
  global_urng().discard(i*SKIP);
}

extern "C" float ran1(int * idum)
{
     static std::uniform_real_distribution<float> d(0.,1.);
     return d( global_urng() );
}
 
