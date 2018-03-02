#ifndef THRUST_COMMON_H
#define THRUST_COMMON_H

#include <thrust/device_vector.h>
#include <thrust/transform.h>
#include <thrust/sequence.h>
#include <thrust/copy.h>
#include <thrust/fill.h>
#include <thrust/replace.h>
#include <thrust/functional.h>
#include <thrust/iterator/counting_iterator.h>
#include <thrust/random.h>

#include <cmath>
#include <ctime>

namespace edda{


///
/// \brief Simply return a random value between 0..1
///
/// Need to pass a unique index value
///
struct Rand{
  __host__ __device__
  inline float operator() (int index)
  {
    thrust::default_random_engine rng(index<<1);
    return rng();
  }
};

namespace detail {
  ///
  /// \brief Create a Thrust default random engine with a given seed
  ///
  /// The caller should ensure that the seeds provided are not repetative in patterns.
  ///
  struct GenRand: public thrust::unary_function<int, thrust::default_random_engine>
  {
    __device__
    thrust::default_random_engine operator () (int idx)
    {
        thrust::default_random_engine randEng;
        randEng.discard(idx << 1);
        return randEng;
    }
  };
}

///
/// \brief randomEngineIterator creates a Thrust default random engine with a given seed
/// \return A Thrust iterator
///
/// Note: The caller should ensure that the seeds provided are not repetative in patterns
///
inline thrust::transform_iterator<detail::GenRand, thrust::counting_iterator<int> >
randomEngineIterator(int seed) {
  return thrust::make_transform_iterator(thrust::make_counting_iterator(seed), detail::GenRand()) ;
}

// Calling Cuda math functions
namespace Math{

} // Math

} // edda

#endif // THRUST_COMMON_H
