// Copyright 2015 The Edda Authors. All rights reserved.
// Use of this source code is governed by a MIT-style license that can be
// found in the LICENSE file.

#ifndef GMM_VTK_ARRAY_H
#define GMM_VTK_ARRAY_H

#include <vector>
#include <stdexcept>
#include <string>

#include <vtkImageData.h>
#include <vtkSmartPointer.h>
#include <vtkDataArray.h>
#include <vtkFloatArray.h>
#include <vtkPointData.h>
#include <vtkCellData.h>
#include <vtkFieldData.h>
#include <vtkStructuredGrid.h>

#include <distributions/gaussian_mixture.h>
#include <dataset/distr_array.h>
#include <core/thrust_gmm_array.h>

namespace edda {

//---------------------------------------------------------------------------------------

/// \brief GmmVtkDataArray implements AbstractDataArray.  It holds vtkDataArrys and returns GMMs with interleaved memory accessing.
class GmmVtkDataArray: public DistrArray
{
protected:
  std::vector<vtkSmartPointer<vtkDataArray> > arrays;

  size_t length = 0;
  int components = 1;
  int target_comp;
public:
  GmmVtkDataArray(vtkFieldData *fieldData, const char *arrayNamePrefix="")  ;

  /// \brief Based on the input array order, assign mean0, var0, weight0, mean1, var1, weight1,...
  /// The number of arrays should be multiples of 3
  GmmVtkDataArray(std::vector<vtkSmartPointer<vtkDataArray> > arrays_);

  typedef std::vector<vtkSmartPointer<vtkDataArray> > RawArrayType;

  virtual ~GmmVtkDataArray();

  virtual size_t getLength();

  virtual int getNumComponents() ;

  virtual void SetTargetComponent(int idx) {target_comp = idx;}

  virtual int GetTargetComponent() {return target_comp;}

  virtual dist::Variant getDistr(size_t idx);

  virtual std::vector<dist::Variant> getDistrVector(size_t idx);

  virtual Real getScalar(size_t idx);

  virtual std::vector<Real> getVector(size_t idx);

  virtual boost::any getRawArray() { return boost::any(arrays); }

  std::shared_ptr<GmmArray> genNdArray() ;

  virtual std::string getDistrName() {
    std::stringstream ss;
    ss << "GaussianMixture" << MAX_GMs;
    return ss.str();
  }
};

} // namespace edda

#endif // GMM_VTK_ARRAY_H
