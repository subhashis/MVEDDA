#ifndef NDARRAY_ITERATOR
#define NDARRAY_ITERATOR

#include <thrust/iterator/iterator_adaptor.h>
#include "common.h"
#include "core/thrust_common.h"
#include "core/ndarray.h"


namespace edda {

template<typename >
class ndarray_iterator
    : public thrust::iterator_adaptor<
        ndarray_iterator<T>, // the first template parameter is the name of the iterator we're creating
        Iterator                   // the second template parameter is the name of the iterator we're adapting
                                   // we can use the default for the additional template parameters
      >
{
  public:
    // shorthand for the name of the iterator_adaptor we're deriving from
    typedef thrust::iterator_adaptor<
      repeat_iterator<Iterator>,
      Iterator
    > super_t;
    __host__ __device__
    repeat_iterator(const Iterator &x, int n) : super_t(x), begin(x), n(n) {}
    // befriend thrust::iterator_core_access to allow it access to the private interface below
    friend class thrust::iterator_core_access;
  private:
    // repeat each element of the adapted range n times
    unsigned int n;
    // used to keep track of where we began
    const Iterator begin;
    // it is private because only thrust::iterator_core_access needs access to it
    __host__ __device__
    typename super_t::reference dereference() const
    {
      return *(begin + (this->base() - begin) / n);
    }
};




class ndarray_iterator {

private:
    NdArray &ndarray;

public:
    ndarray_iterator(const NdArray &ndarray_) {
         ndarray = ndarray_;
    }

    Type operator() ()
};


struct GetSample_functor {
  template <class Dist> //, ENABLE_IF_BASE_OF(Dist, dist::Distribution)>
  __host__ __device__
  double operator() (thrust::tuple<Dist, thrust::default_random_engine> &tuple)
  {
    Dist &dist = thrust::get<0>(tuple);
    return dist::getSample(dist, thrust::get<1>(tuple));
  }
};

///
/// Output iterator must point to elements in type 'Real'
///
template <class InputIterator, class OutputIterator>
void randomSampleField(InputIterator begin, InputIterator end, OutputIterator out)
{
  static int seed ;
  seed += time(NULL) % 10000;
  int n = end-begin;
  thrust::transform( thrust::make_zip_iterator(thrust::make_tuple(begin, randomEngineIterator(seed)) ) ,
                     thrust::make_zip_iterator(thrust::make_tuple(end, randomEngineIterator(seed+n)) ),
                     out, GetSample_functor()) ;
  seed += n;
}

} // edda

#endif // NDARRAY_ITERATOR
