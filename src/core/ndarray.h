// Copyright 2015 The Edda Authors. All rights reserved.
// Use of this source code is governed by a MIT-style license that can be
// found in the LICENSE file.

#ifndef NDARRAY_H_
#define NDARRAY_H_

#include <cstdio>
#include <cassert>
#include <iostream>
#include <memory>
#include <thrust/random/uniform_real_distribution.h>
#ifdef __CUDACC__
#include <cuda.h>
#endif

#include "thrust_common.h"

namespace edda {

const int kNdArrayMaxDims = 8;

///
/// \brief A general n-dimensional array compatible to Cuda
///
template <typename Type>
class NdArray {

  thrust::device_ptr<Type> dData;
  NdArray<Type> *dSelfPtr = 0;
  int num_of_dims = 0;
  int num_of_elems = 0;
  int dims[kNdArrayMaxDims];
  int strides[kNdArrayMaxDims];
  bool ownership = false;

public:
  typedef NdArray<Type> *SelfDevicePtr;

  ///
  /// \brief Create an empty array
  ///
  __host__ __device__
  NdArray(): ownership(false) { }

  ///
  /// \brief Create an empty array
  ///
  NdArray(int num_of_dims_, int* dims_) {

    init_shape(num_of_dims_, dims_);

    dData = thrust::device_malloc<Type>(num_of_elems);

    ownership = true;
    printf("NdArray created\n");
  }

  ///
  /// \brief Create an empty array
  ///
  NdArray(const std::initializer_list<int>& dims_) {

    init_shape(dims_);

    dData = thrust::device_malloc<Type>(num_of_elems);

    ownership = true;
    printf("NdArray created\n");
  }

  ///
  /// \brief Create a device array and copy host content to it.
  ///
  NdArray(Type* data, int num_of_dims_, int* dims_)
    : NdArray(num_of_dims_, dims_)
  {
    thrust::copy(data, data+num_of_elems, dData );
  }

  ///
  /// \brief Create a device array and copy host content to it.
  ///
  NdArray(Type* data, const std::initializer_list<int>& dims_)
    : NdArray( dims_ )
  {
    thrust::copy(data, data+num_of_elems, dData );

  }

  ///
  /// \brief pass a device pointer.
  ///
  /// Note: NdArray now owns the array
  ///
  __host__ __device__
  NdArray(thrust::device_ptr<Type> dData, int num_of_dims_, int* dims_)
    : NdArray( num_of_dims_, dims_ )
  {
    this->dData  = dData;

  }

  NdArray(const NdArray &obj) {
    this->operator=( obj );
  }

  ///
  /// \brief Assignment will be a shallow copy of the input
  ///
  NdArray<Type> &operator=(const NdArray<Type> &obj) {
    printf("opeartor=\n");
    this->dData = obj.dData;
    this->num_of_dims = obj.num_of_dims;
    this->num_of_elems = obj.num_of_elems;
    for (int i=0; i<num_of_dims; i++) {
      this->dims[i] = obj.dims[i];
      this->strides[i] = obj.strides[i];
    }
    this->ownership = false;
    return *this;
  }

  ///
  /// \brief Take over the ownership of the input
  ///
  void take(NdArray<Type> &obj) {
    free(); // free the current data
    *this = obj;
    this->ownership = true;
    obj.ownership = false;
  }

  ~NdArray() {
    free();
  }

  void free() {
    if (this->ownership ) {
      printf("Release NdArray\n");
      thrust::device_free( dData );
    }

#ifdef __CUDACC__
      if ( dSelfPtr )
        cudaFree( dSelfPtr );
#endif
  }

  template <typename OutputIterator>
  void copy_to_host(OutputIterator out) {

    thrust::copy(begin(), end(), out);
  }

  __host__ __device__ int get_num_of_dims() const {
    return num_of_dims;
  }

  __host__ __device__ int get_num_of_elems() const {
    return num_of_elems;
  }

  __host__ __device__ const int* get_dims() const {
    return dims;
  }

  __host__ __device__
  void set_ownership(bool own) {
    this->ownership = own;
    printf("Set Ownership\n");
  }

  __host__ __device__
  thrust::device_ptr<Type> & data() {
    return dData;
  }

  const thrust::device_ptr<Type> begin() const {
    return dData;
  }

  const thrust::device_ptr<Type> end() const {
    return dData + num_of_elems;
  }

  __host__ __device__ thrust::device_ptr<Type> get_ptr(const std::initializer_list<int>& ind) const {
    int nd = num_of_dims;
    const int* s = strides;
    thrust::device_ptr<Type> dptr = dData;
    auto it = ind.begin();
    while (nd--) {
      dptr += (*s++) * (*it++);
    }
    return dptr;
  }

  __device__ __host__ Type get_val(const std::initializer_list<int>& ind) const {
    return *(get_ptr(ind));
  }

  __device__ __host__ Type get_val(int ind) const {
      //printf("val[%d] = %f\n", ind, (float) dData[ind]);
    return dData[ind];
  }

  __host__ __device__ void set_val(
      const std::initializer_list<int>& ind, const Type& val) {
    *(get_ptr(ind)) = val;
  }

  __host__ __device__ void Reshape(
      const std::initializer_list<int>& newshape) {

    // total size check
    int newsize=1;
    auto it = newshape.begin();
    for (; it != newshape.end(); ++it)
      newsize *= *it;
    assert (newsize == num_of_elems);

    num_of_dims = newshape.size();

    it = newshape.begin();
    for (int i = 0; i < num_of_dims; ++i, ++it)
      dims[i] = *it;

    UpdateStrides();
  }

  // NdArray* Slice(const std::initializer_list<int>& from,
  //                const std::initializer_list<int>& to) {
  //   int* dims = new int[num_of_dims];
  //   auto itf = from.begin(), itt = to.begin();
  //   for (int i = 0; i < num_of_dims; ++i) {
  //     dims[i] = *(itt + i) - *(itf + i) + 1;
  //   }

  //   int* strides = new int[num_of_dims];
  //   for (int i = 0; i < num_of_dims; ++i) {
  //     strides[i] = strides_[i];
  //   }

  //   return new NdArray(get_ptr(from), num_of_dims, dims, strides);
  // }

  __host__
  SelfDevicePtr get_selft_ptr() {
#ifdef __CUDACC__
    if (!dSelfPtr) {

      cudaMalloc(&dSelfPtr, sizeof(NdArray<Type>));
      cudaMemcpy(dSelfPtr, this, sizeof(NdArray<Type>), cudaMemcpyHostToDevice);

    }
    return dSelfPtr;
#else
    return this;
#endif
  }

 private:
  void init_shape(const std::initializer_list<int>& dims_) {

    assert(dims_.size() <= kNdArrayMaxDims);

    num_of_dims = dims_.size();

    auto it = dims_.begin();
    for (int i = 0; i < num_of_dims; ++i)
      dims[i] = *(it + i);

    UpdateStrides();

    num_of_elems = 1;
    for (int i = 0; i < num_of_dims; ++i)
      num_of_elems *= dims[i];
  }

  void init_shape(int num_of_dims_, int* dims_) {

    assert(num_of_dims_ <= kNdArrayMaxDims);

    num_of_dims = num_of_dims_;

    for (int i = 0; i < num_of_dims; ++i)
      dims[i] = dims_[i];

    UpdateStrides();

    num_of_elems = 1;
    for (int i = 0; i < num_of_dims; ++i)
      num_of_elems *= dims[i];

  }

  __host__ __device__ void UpdateStrides() {
    int stride = 1;
    for (int i = num_of_dims - 1; i >= 0; --i) {
      strides[i] = stride;
      stride *= dims[i] ? dims[i] : 1;
    }
  }

};



}  // namespace edda

#endif  // NDARRAY_H_
