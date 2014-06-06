//#include <random>
//#include <functional>
#include <limits>
#include <cmath>
#include <iostream>

#include "prng_engine.hpp"

#define SEED  1234321

/*
 * copied from standard c++ library
 */
template<typename _RealType,
        typename _UniformRandomNumberGenerator>
    _RealType
    generate_canonical(_UniformRandomNumberGenerator& __urng)
    {
      const size_t __b = static_cast<size_t>(std::numeric_limits<_RealType>::digits);
      const long double __r = static_cast<long double>(__urng.max())
                - static_cast<long double>(__urng.min()) + 1.0L;
      const size_t __log2r = std::log(__r) / std::log(2.0L);
      size_t __k = std::max<size_t>(1UL, (__b + __log2r - 1UL) / __log2r);
      _RealType __sum = _RealType(0);
      _RealType __tmp = _RealType(1);
      for (; __k != 0; --__k)
      {
          __sum += _RealType(__urng() - __urng.min()) * __tmp;
          __tmp *= __r;
      }
      return __sum / __tmp;
    }

sitmo::prng_engine & global_urng()
{
  static sitmo::prng_engine u;
  return u;
}

extern "C" void ran_seed(long long i)
{
  global_urng().seed(SEED);
  global_urng().discard(i);
}

extern "C" float ran1(int * idum)
{
     return generate_canonical<float, sitmo::prng_engine>( global_urng() );
}
 
// fortran callable
extern "C" void ran_seed_(long long *i)
{
  ran_seed(*i);
}
extern "C" float ran1_(int * idum)
{
  return ran1(idum);
}



