// Copyright 2015 The Edda Authors. All rights reserved.
// Use of this source code is governed by a MIT-style license that can be
// found in the LICENSE file.

#ifndef DISTR_ARRAY_H_
#define DISTR_ARRAY_H_

#include <cassert>
#include <iostream>
#include <memory>
#include <stdexcept>
#include <vector>
#include <boost/any.hpp>

#include "distributions/variant.h"
#include "core/vector_matrix.h"
#include "core/shared_ary.h"

namespace edda {

///
/// \brief The class interface of distrbution arrays used by the Dataset class
///
class DistrArray
{
public:

  virtual ~DistrArray() {}

  ///
  /// Get the number of elements
  ///
  virtual size_t getLength() =0;

  ///
  /// Get the number of element components.  This is used for elements in variable-length vector or matrix.
  ///
  virtual int getNumComponents() = 0;

  ///
  /// Set the target component index of vector data when calling getScalar()
  ///
  virtual void SetTargetComponent(int idx) = 0;

  ///
  /// Get the target component index of vector data
  ///
  virtual int GetTargetComponent() = 0;

  ///
  /// Get random sampling of scalar distribution
  ///
  virtual Real getScalar(size_t idx) =0;

  ///
  /// Get random sampling of vector distribution
  ///
  virtual std::vector<Real> getVector(size_t idx)=0;

  virtual dist::Variant getDistr(size_t idx) =0;

  virtual std::vector<dist::Variant> getDistrVector(size_t idx)=0;  //!!! should be deprecated !!!

  //virtual void setDistr(size_t idx, dist::Variant) =0;

  //virtual void setDistrVector(size_t idx, std::vector<dist::Variant>)=0;

  ///
  /// Get the array
  ///
  virtual boost::any getRawArray() = 0;

  ///
  /// Get gistribution name for data writer
  ///
  virtual std::string getDistrName() = 0;

};

//---------------------------------------------------------------------------------------
/// \brief Array of scalar distributions.
///
/// This is the class that holds the actual array, with a smart pointer.
/// \param T The element type of an array.
///
template<typename Distr, ENABLE_IF_BASE_OF(Distr, dist::DistributionTag)>
class ScalarDistrArray: public DistrArray
{
protected:
  shared_ary<Distr> array;
public:
  ScalarDistrArray(shared_ary<Distr> array) { this->array = array; }

  virtual ~ScalarDistrArray() { }

  virtual size_t getLength() { return array.getLength(); }

  virtual int getNumComponents() { return 1; }

  virtual void SetTargetComponent(int idx) { /*do nothing*/ }

  virtual int GetTargetComponent() {return 0;}

  virtual dist::Variant getDistr(size_t idx) { return array[idx]; }

  virtual std::vector<dist::Variant> getDistrVector(size_t idx) { return std::vector<dist::Variant> (1, array[idx] ); }

  virtual Real getScalar(size_t idx) { return getSample(array[idx]); }

  virtual std::vector<Real> getVector(size_t idx) { return std::vector<Real> (1, Real(getSample(array[idx])) ); }

  //virtual boost::any getItem(size_t idx) { return boost::any( array[idx] );  }

  //virtual void setItem(size_t idx, int component, const boost::any &item) { array[idx] = boost::any_cast<T>( item );  }

  virtual boost::any getRawArray() { return boost::any(array); }

  virtual std::string getDistrName() { return getName(Distr()); }
};


//---------------------------------------------------------------------------------------

///
/// \brief Array of vector distributions
///
/// Note the distributions among vector components are independent,
/// which is in the contrast of JointDistrArray
///


//!!! should be deprecated !!!
template<typename Distr, int N, ENABLE_IF_BASE_OF(Distr, dist::DistributionTag)>
class VectorDistrArray: public DistrArray
{
protected:
  shared_ary<Vector<Distr,N> > array;
  int target_comp;
public:
  VectorDistrArray(shared_ary<Vector<Distr,N> > array) { this->array = array; target_comp = 0;}

  virtual ~VectorDistrArray() { }

  virtual size_t getLength() { return array.getLength(); }

  virtual int getNumComponents() { return N; }

  virtual void SetTargetComponent(int idx) {target_comp = idx;}

  virtual int GetTargetComponent() {return target_comp;}

  virtual dist::Variant getDistr(size_t idx) {
      throw std::runtime_error("Requesting scalar in a VectorArray");
    }

  virtual std::vector<dist::Variant> getDistrVector(size_t idx) {
    std::vector<dist::Variant> v(N);
    for (int i=0; i<N; i++) {
      v[i] = array[idx][i];
    }
    return v;
  }

  virtual Real getScalar(size_t idx) {
    std::vector<Real> vec = getVector(idx);
    return vec[target_comp];
  }

  virtual std::vector<Real> getVector(size_t idx) {
    std::vector<Real> v(N);
    for (int i=0; i<N; i++) {
      v[i] = getSample(array[idx][i]);
    }
    return v;
  }

  //virtual boost::any getItem(size_t idx) { return boost::any( array[idx] );  }

  //virtual void setItem(size_t idx, int component, const boost::any &item) { array[idx] = boost::any_cast<Vector<T,N> >( item );  }

  virtual boost::any getRawArray() { return boost::any(array); }

  virtual std::string getDistrName() { return getName(Distr()); }
};


//---------------------------------------------------------------------------------------
///
/// \brief An array of joint distributions
///
template<typename Distr, ENABLE_IF_BASE_OF(Distr, dist::DistributionTag)>
class JointDistrArray: public DistrArray
{
protected:
  shared_ary<Distr> array;
  int num_comps;  // number of components
  int target_comp;
public:
  JointDistrArray(shared_ary<Distr> array, int num_comps) { this->array = array; this->num_comps = num_comps; }

  // If the user does not provide the number of components, obtain from the first element of the array
  JointDistrArray(shared_ary<Distr> array) {
    this->array = array;
    this->num_comps = 0;
    if (array.getLength()>0)
      this->num_comps = getJointSample(array[0]).size();
  }

  virtual ~JointDistrArray() { }

  virtual size_t getLength() { return array.getLength(); }

  virtual int getNumComponents() { return num_comps; }

  virtual void SetTargetComponent(int idx) {target_comp = idx;}

  virtual int GetTargetComponent() {return target_comp;}

  virtual dist::Variant getDistr(size_t idx) { return array[idx]; }


  /// !!! should be deprecated, since we do not use VECTOR3 any more, but restore the vector type into three independent distr_array
  virtual std::vector<dist::Variant> getDistrVector(size_t idx) {
    throw std::runtime_error("Requesting a vector of distributions in a joint distribution array.");
  }  

  virtual Real getScalar(size_t idx) { return getJointSample(array[idx])[target_comp]; }

  virtual std::vector<Real> getVector(size_t idx) { return getJointSample(array[idx]); }

  //virtual boost::any getItem(size_t idx) { return boost::any( array[idx] );  }

  //virtual void setItem(size_t idx, int component, const boost::any &item) { array[idx] = boost::any_cast<T>( item );  }

  virtual boost::any getRawArray() { return boost::any(array); }

  virtual std::string getDistrName() { return getName(Distr()); }
};


} // namespace edda

#endif // DISTR_ARRAY_H_
