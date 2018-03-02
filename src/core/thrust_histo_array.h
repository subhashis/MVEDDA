#ifndef GMM_ARRAY_H
#define GMM_ARRAY_H

#include <memory>
#include <thrust/iterator/constant_iterator.h>
#include <thrust/host_vector.h>

#include <core/ndarray.h>
#include <core/thrust_common.h>
#include <distributions/histogram.h>

namespace edda{

  ///
  ///  Histogram array class: integrate separate arrays to Thrust-readable arrays
  ///

namespace detail {
  struct MakeHistogram: public thrust::unary_function<int, dist::Histogram > {

    MakeHistogram() {}

    __host__ __device__
    dist::Histogram operator() (int idx) const
    {
      dist::Histogram histo;
      // TODO
      return histo;
    }

  };

} // detail


class HistoArray{


  detail::MakeHistoArray getMakeHistogram() {
    return detail::MakeHistogram();
  }
public:

  int narrays;
  int num_of_elems;

  HistoArray() {}

  HistoArray(std::vector<NdArray<Real> > &data_) {

  }


  inline thrust::transform_iterator<detail::MakeHistogram, thrust::counting_iterator<int> >
  begin() {
    return thrust::make_transform_iterator( thrust::make_counting_iterator(0),
                                            detail::MakeHistogram(dDataPtrArray, narrays)
          );
  }

  inline thrust::transform_iterator<detail::MakeHistogram, thrust::counting_iterator<int> >
  end() {
    return thrust::make_transform_iterator( thrust::make_counting_iterator(num_of_elems),
                                            detail::MakeHistogram(dDataPtrArray, narrays)
        );
  }


};

} // edda

#endif // #ifndef GMM_ARRAY_H
