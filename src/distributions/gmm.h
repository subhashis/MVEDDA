// Copyright 2015 The Edda Authors. All rights reserved.
// Use of this source code is governed by a MIT-style license that can be
// found in the LICENSE file.

#ifndef GMM_H
#define GMM_H

#include <cmath>
#include <cstdlib>
#include <iostream>
#include <vector>

#include <boost/math/distributions.hpp>
#include "common.h"
#include "distribution_tag.h"
#include "gaussian.h"
#include "core/statistics.h"

namespace edda {
namespace dist {

//struct GMMTuple {
//  union{
//    struct{Real m,v,w;};  // mean, variance, weight
//    Real p[3];
//  };
//};

// ------------------------------------------------------------------------------
///
/// \brief Defines a Gaussian Mixture class
///
class EDDA_EXPORT GMM: public ContinuousDistributionTag {

  void modelReduction(const std::vector<GMMTuple> &models_) {
    if (models_.size() <=  GMs) {
      init();
      for (size_t i=0; i<models_.size(); i++)
        models[i] = models_[i];

    } else {
      // Too many Gaussian models
      throw NotImplementedException();
    }
  }

  ///
  /// \brief Reset all model weights to 0
  ///
  __host__ __device__
  void init () {
    //memset(this, 0, sizeof(GaussianMixture<GMs>));
    models[0].v = 1;
    models[0].m = 0;
    models[0].w = 1;
    for (int i=1; i<GMs; i++)
      models[i].w = 0;
  }

public:
  //Tuple<GMMTuple, GMs> models;
  std::vector<GMMTuple> models;
  int GMs;

  ///
  /// \brief Constructor
  ///
  __host__ __device__
	  GMM() {
		  GMs = 1;
		  models.resize(GMs);
		  init();
	  }

  ///
  /// \brief Set GMM by a given GMM
  /// \param gms number of Gaussians of GMM
  ///
  __host__ __device__
  GMM(int gms) { 
	  GMs = gms;
	  models.resize(GMs);
	  init(); 
  }

  ///
  /// \brief Set GMM by a given GMM
  /// \param models a given GMM
  ///
  void assign (const std::vector<GMMTuple> &models) {
    std::vector<GMMTuple> vmodels_;
    for (int i=0; i<GMs; i++)
    {
      GMMTuple t = models[i];
      vmodels_.push_back(t);
    }
    modelReduction(vmodels_);
    normalizeWeights();
  }

  ///
  /// \brief Constructor
  /// \param models a given GMM
  ///
  GMM(const std::vector<GMMTuple> &models_) {
    modelReduction(models_);
    normalizeWeights();
  }

  ///
  /// \brief Return number of Gaussian components
  ///
  int getNumComponenets() const
  {
	  return GMs;
  }

  ///
  /// \brief Scale sum of weights to 1.
  ///
  __host__ __device__
  void normalizeWeights() {
    double sum = 0;
    int i;
    for (i=0; i<GMs; i++) {
      sum += models[i].w;
    }
    if (sum == 0) return;
    for (i=0; i<GMs; i++) {
      models[i].w /= sum;
    }
  }
};

// ------------------------------------------------------------------------------
// Below defines GaussianMixture related generic functions

///
/// \brief Return mean
/// \param dist a distribution
///
inline double getMean(const GMM &dist)
{
  double mean = 0;
  for (int i = 0; i<dist.getNumComponenets(); i++)
  {
    mean += dist.models[i].m * dist.models[i].w;
  }
  return mean;
}

///
/// \brief Return variance.
/// \param dist a distribution
///
/// if f(x) = sum_i( w_i * f_i(x) ), v_i = variance of f_i, and m_i = mean of f_i, then
/// var(f) = sum_i( w_i * v_i + w_i * m_i^2 ) - (sum_i( w_i * m_i) )^2.
///
/// ref: http://stats.stackexchange.com/questions/16608/what-is-the-variance-of-the-weighted-mixture-of-two-gaussians
/// (code not verified)
///
inline double getVar(const GMM &dist)
{
  // Let the first summation as term1 and second as term2
  double term1=0, term2=0;
  for (int i = 0; i<dist.getNumComponenets(); i++)
  {
    const GMMTuple & model = dist.models[i];
    term1 += (double)model.w * model.v + (double)model.w * model.m * model.m ;
    term2 += (double)model.w * model.m;
  }
  return term1 - term2 * term2;
}

///
/// \brief Return PDF of x
/// \param dist a distribution
/// \param x a sample (value)
///
inline double getPdf(const GMM &dist, const double x)
{
  double p=0;
  for (int i = 0; i<dist.getNumComponenets(); i++)
  {
    p += getPdf( Gaussian(dist.models[i].m, dist.models[i].v), x ) * dist.models[i].w;
  }
  return p;
}

///
/// \brief Return a sample
/// \param dist a distribution
///
/// Note: To ensure correct sampling distribution, the sum of weights
/// should be normalized to 1 before calling this function.
///
inline double getSample(const GMM &dist)
{
  float ratio = rand() / (float)RAND_MAX;
  float accumulated = 0;
  for (int i = dist.getNumComponenets() - 1; i >= 0; i--)
  {
    accumulated += dist.models[i].w;
    if (ratio < accumulated) {
      return getSample( Gaussian(dist.models[i].m, dist.models[i].v) );
    }
  }
  // return sample from the last model
  return getSample( Gaussian(dist.models[0].m, dist.models[0].v) );
}


inline double getInverseCDF_GMM(const GMM &dist, const double x)
{
  //float ratio = rand() / (float)RAND_MAX;
  float accumulated = 0;
  for (int i = dist.getNumComponenets() - 1; i >= 0; i--)
  {
    accumulated += dist.models[i].w;
    if (x < accumulated) {
      return getSample( Gaussian(dist.models[i].m, dist.models[i].v) );
    }
  }
  // return sample from the last model
  return getSample( Gaussian(dist.models[0].m, dist.models[0].v) );
}

///
/// \brief Return a sample
/// \param dist a distribution
/// \param rng a random engine
///
/// Note: To ensure correct sampling distribution, the sum of weights
/// should be normalized to 1 before calling this function.
///
__host__ __device__
inline double getSample(const GMM &dist, thrust::default_random_engine &rng)
{
  thrust::uniform_real_distribution<Real> uniform;
  float ratio = uniform(rng);
  float accumulated = 0;
  for (int i = dist.getNumComponenets() - 1; i >= 0; i--)
  {
    accumulated += dist.models[i].w;
    if (ratio < accumulated) {
      return getSample( Gaussian(dist.models[i].m, dist.models[i].v), rng );
    }
  }
  // return sample from the last model
  return getSample( Gaussian(dist.models[0].m, dist.models[0].v), rng );
}


///
/// \brief Return CDF of x
/// \param dist a distribution
/// \param x a sample (value)
///
__host__ __device__
inline double getCdf(const GMM &dist, double x)
{
  double cdf=0;
  for (int i = 0; i<dist.getNumComponenets(); i++)
  {
    if (dist.models[i].w >0)
      cdf += getCdf(Gaussian(dist.models[i].m, dist.models[i].v), x) * dist.models[i].w;
  }
  return cdf;
}

///
/// \brief Print itself
/// \param os outstream
/// \param dist a distribution
///
inline std::ostream& operator<<(std::ostream& os, const GMM &dist)
{
  os <<  "<UnivariateGMM(m,v,w):";
  for (int i = 0; i<dist.getNumComponenets(); i++)
    os << " (" << dist.models[i].m << "," << dist.models[i].v << "," << dist.models[i].w << ")";
  os << ">";
  return os;
}

///
/// \brief Return distribution information
/// \param x a distribution
///
__host__ __device__
inline std::string getName(const GMM &x) {
  std::stringstream ss;
  ss << "GaussianMixture" << x.getNumComponenets();
  return ss.str();
}

// ------------------------------------------------------------------------------
// Below defines Gaussian related arithmetics

///
/// \brief random variable with unary -
/// \param x a distribution
///
inline GMM operator-(const GMM &x)
{
  throw NotImplementedException();

}

///
/// \brief random variable +=
/// \param x distribution1
/// \param rhs distribution2
///
inline GMM& operator+=(GMM &x, const GMM& rhs) {
  throw NotImplementedException();
}

///
/// \brief random variable += with scalar
/// \param x distribution1
/// \param r a value which is applied to means of GMM
///
inline GMM& operator+=(GMM &x, const double r) {
	for (int i = 0; i<x.getNumComponenets(); i++)
    x.models[i].m += r;
  return x;
}

///
/// \brief random variable *= with scalar
/// \param x distribution1
/// \param r a value which is applied to means and variances of GMM
///
inline GMM& operator*=(GMM &x, const double r) {
	for (int i = 0; i<x.getNumComponenets(); i++) {
    x.models[i].m *= r;
    x.models[i].v *= r;
  }
  return x;
}

//typedef GaussianMixture<MAX_GMs> DefaultGaussianMixture;

}  // namespace dist



}  // namespace edda



#endif  // GAUSSIAN_MIXTURE_H
