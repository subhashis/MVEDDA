#ifndef GMM_ARRAY_H
#define GMM_ARRAY_H

#include <memory>
#include <thrust/iterator/constant_iterator.h>
#include <thrust/host_vector.h>

#include <core/ndarray.h>
#include <core/thrust_common.h>
#include <distributions/gaussian_mixture.h>

namespace edda{

typedef Tuple<NdArray<Real>::SelfDevicePtr, MAX_GMs*3> DeviceGMMArray;

namespace detail {
    struct MakeStridedGmm: public thrust::unary_function<int, dist::DefaultGaussianMixture > {
      // [GMMs][n] array
      DeviceGMMArray dDataPtrArray;
      int narrays;

      MakeStridedGmm() {}
      MakeStridedGmm(const DeviceGMMArray &dDataPtrArray_, int narrays_) : dDataPtrArray(dDataPtrArray_), narrays(narrays_) {}

      __host__ __device__
      dist::DefaultGaussianMixture operator() (int idx) const
      {
        dist::DefaultGaussianMixture gmm;
        int narray = narrays;
        for (int i=0; i<narray; i++) {
          gmm.models[i/3].p[i%3] = dDataPtrArray[i]->get_val(idx);
        }
        if (narray==2)
          gmm.models[0].w = 1.;
        return gmm;
      }

    };

} // detail

class GmmArray{
  ///
  /// This is a collection of 1D Real arrays in the order of mean0, var0, weight0, mean1, var1, weight1...
  /// The minimum numbers of arrays is 2 (mean and var)
  ///
  Tuple<NdArray<Real>, MAX_GMs*3> dataArray;
  Tuple<NdArray<Real>::SelfDevicePtr, MAX_GMs*3> dDataPtrArray;

  detail::MakeStridedGmm getMakeStridedGmm() {
    return detail::MakeStridedGmm(dDataPtrArray, narrays);
  }
public:

  int narrays;
  int num_of_elems;

  GmmArray() {}

  GmmArray(std::vector<NdArray<Real> > &data_) {
    narrays = data_.size();
    if ( narrays > MAX_GMs*3 ) {
      printf ( "Warning: The provided GMM models are larger than the maximum size.  Extra models will be discarded!");
      narrays = MAX_GMs*3;
    }

    assert(narrays >= 2);

    num_of_elems = data_[0].get_num_of_elems();

    for (int i=0; i<narrays; i++)
    {
      dataArray[i].take( data_[i] );
      dDataPtrArray[i] = dataArray[i].get_selft_ptr();
    }
  }


  inline thrust::transform_iterator<detail::MakeStridedGmm, thrust::counting_iterator<int> >
  begin() {
    return thrust::make_transform_iterator( thrust::make_counting_iterator(0),
                                            detail::MakeStridedGmm(dDataPtrArray, narrays)
          );
  }

  inline thrust::transform_iterator<detail::MakeStridedGmm, thrust::counting_iterator<int> >
  end() {
    return thrust::make_transform_iterator( thrust::make_counting_iterator(num_of_elems),
                                            detail::MakeStridedGmm(dDataPtrArray, narrays)
        );
  }


};

} // edda

#endif // #ifndef GMM_ARRAY_H
