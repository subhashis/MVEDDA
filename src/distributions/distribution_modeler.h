// Copyright 2015 The Edda Authors. All rights reserved.
// Use of this source code is governed by a MIT-style license that can be
// found in the LICENSE file.

#ifndef DISTRIBUTION_MODELER_H
#define DISTRIBUTION_MODELER_H

#include <iostream>
#include <string>
#include <vector>
#include <stdlib.h>
#include "edda.h"
#include "distributions/gaussian.h"
#include "distributions/distribution.h"
#include <distributions/gaussian_mixture.h>
#include <distributions/gmm.h>
#include "distributions/histogram.h"
#include <distributions/estimate_gmm.h>
#include "distributions/variant.h"

#include "distributions/joint_GMM.h"
#include "distributions/joint_gaussian.h"
#include "distributions/joint_histogram.h"

#include "dataset/dataset.h"
#include "dataset/distr_array.h"
#include "core/interpolator.h"
#include "io/edda_vtk_writer.h"

using namespace std;

namespace edda{

/// \brief Define of Distribution Modeler class

/// Distributional Modeler class is a high level abstraction of a collection of distributions.
/// Essentially it holds the collection of different distribution objects (GMM, Histograms etc).

class DistributionModeler {
protected:
	size_t len;
	dist::Variant* distArray;
public:
	/// \brief Default constructor
    DistributionModeler()
   	{
        len = 0;
	}
	/// \brief Constructor with specified length
	DistributionModeler(int l)
   	{
        len = l;
        distArray = new dist::Variant[len];
	}
	/// \brief Evalautes the GMM for the specified data
	/// \param data pointer to the array of initial sample(data)
	/// \param size number of sample
	/// \param nGmm number of Gaussian components
	/// \param index corresponding index of the distribution in the array of distributions
	void computeGMM(float *data, size_t size, size_t nGmm, size_t index)
	{
		//TODO: check for out of bound errors.
		double *dataD;		
  		dataD = new double[size];

  		for(size_t i=0; i<size; i++)
	 		dataD[i] = (double) data[i];

	 	dist::GMM new_distr;
		new_distr = eddaComputeGMM(dataD, size, nGmm);

		distArray[index] = new_distr;
	}
	/// \brief Evalautes the Histogram for the specified data
	/// \param data pointer to the array of initial sample(data)
	/// \param size number of sample
	/// \param index corresponding index of the distribution in the array of distributions
	/// \param binCount number of desired bins in the histogram
	void computeHistogram(float *data, size_t size, size_t index, int binCount)
	{
		dist::Histogram new_distr;
		new_distr = eddaComputeHistogram(data, size, binCount);

		distArray[index] = new_distr;
	}

	void computeHistogram(float *data, size_t size, size_t index, int binCount, dist::Variant &ret_distr)
	{
		dist::Histogram new_distr;
		new_distr = eddaComputeHistogram(data, size, binCount);

		ret_distr = new_distr;
		distArray[index] = new_distr;
	}

	/// \brief Evalautes the JointGMM for the specified multivariate data
	/// \param data pointer to the array of initial sample(data)
	/// \param size number of sample
	/// \param nGmm number of Gaussian components
	/// \param index corresponding index of the distribution in the array of distributions
	void computeJointGMM(std::vector<Real*>& data, size_t size, size_t nGmm, size_t index)
	{
		dist::JointGMM new_distr;
		new_distr = eddaComputeJointGMM(data, size, nGmm);
		
		distArray[index] = new_distr;
	}
	/// \brief Evalautes the JointHistogram for the specified multivariate data
	/// \param data pointer to the array of initial sample(data)
	/// \param size number of sample
	/// \param mins vector of minimum values for the different variables.
	/// \param maxs vector of maximum values for the different variables.
	/// \param nBins number of desired bins in the histogram
	/// \param index corresponding index of the distribution in the array of distributions
	void computeJointHistogram(std::vector<Real*>& data, size_t size, const std::vector<Real>& mins, const std::vector<Real>& maxs, const std::vector<int>& nBins, size_t index)
	{
		dist::JointHistogram new_distr;
		new_distr = dist::eddaComputeJointHistogram(data, size, mins, maxs, nBins);
		
		distArray[index] = new_distr;
	}
	/// \brief Returns the array of distributions
	DistrArray * getDistrArray()
	{
		shared_ary<dist::Variant> pArray (new dist::Variant[len], len);
		for(int i=0; i<len; i++)
		{
			pArray[i] = distArray[i];
		}
		DistrArray * abs_array = new ScalarDistrArray<dist::Variant>(pArray);
		return abs_array;
	}
	/// \brief Returns the array of joint distributions
	/// \param numComp number of variables involved
	DistrArray * getMVDistrArray(int numComp)
	{
		shared_ary<dist::Variant> pArray (new dist::Variant[len], len);
		for(int i=0; i<len; i++)
		{
			pArray[i] = distArray[i];
		}
		DistrArray * abs_array = new JointDistrArray<dist::Variant>(pArray,numComp);
		return abs_array;
	}
	
};

}
#endif // DISTRIBUTION_MODELER_H
