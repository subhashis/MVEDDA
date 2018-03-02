// Copyright 2015 The Edda Authors. All rights reserved.
// Use of this source code is governed by a MIT-style license that can be
// found in the LICENSE file.

#ifndef GAUSSIAN_MIXTURE_H
#define GAUSSIAN_MIXTURE_H

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

struct GMMTuple {
  union{
    struct{Real m,v,w;};  // mean, variance, weight
    Real p[3];
  };
};

// ------------------------------------------------------------------------------
///
/// \brief Defines a Gaussian Mixture class
/// GMs: number of Gaussian Models
///
template <int GMs>
class EDDA_EXPORT GaussianMixture: public ContinuousDistributionTag {

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
  Tuple<GMMTuple, GMs> models;

  ///
  /// \brief Constructor
  ///
  __host__ __device__
  GaussianMixture() { init(); }

  ///
  /// \brief Set GMM by a given GMM
  /// \param models a given GMM
  ///
  template <int GMs_>
  void assign (const Tuple<GMMTuple, GMs_> &models) {
    std::vector<GMMTuple> vmodels_;
    for (int i=0; i<GMs_; i++)
    {
      GMMTuple t = models[i];
      vmodels_.push_back(t);
    }
    modelReduction(vmodels_);
    normalizeWeights();
  }

  ///
  /// \brief Constructor :by a given GMM
  /// \param models a given GMM
  ///
  GaussianMixture(const std::vector<GMMTuple> &models_) {
    modelReduction(models_);
    normalizeWeights();
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
template <int GMs>
inline double getMean(const GaussianMixture<GMs> &dist)
{
  double mean = 0;
  for (int i=0; i<GMs; i++)
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
template <int GMs>
inline double getVar(const GaussianMixture<GMs> &dist)
{
  // Let the first summation as term1 and second as term2
  double term1=0, term2=0;
  for (int i=0; i<GMs; i++)
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
template <int GMs>
inline double getPdf(const GaussianMixture<GMs> &dist, const double x)
{
  double p=0;
  for (int i=0; i<GMs; i++)
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
template <int GMs>
inline double getSample(const GaussianMixture<GMs> &dist)
{
  float ratio = rand() / (float)RAND_MAX;
  float accumulated = 0;
  for (int i=GMs-1; i >= 0; i--)
  {
    accumulated += dist.models[i].w;
    if (ratio < accumulated) {
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
template <int GMs>
__host__ __device__
inline double getSample(const GaussianMixture<GMs> &dist, thrust::default_random_engine &rng)
{
  thrust::uniform_real_distribution<Real> uniform;
  float ratio = uniform(rng);
  float accumulated = 0;
  for (int i=GMs-1; i>=0; i--)
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
template <int GMs>
__host__ __device__
inline double getCdf(const GaussianMixture<GMs> &dist, double x)
{
  double cdf=0;
  for (int i=0; i<GMs; i++)
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
template <int GMs>
inline std::ostream& operator<<(std::ostream& os, const GaussianMixture<GMs> &dist)
{
  os <<  "<GaussianMixture(m,v,w):";
  for (int i=0; i<GMs; i++)
    os << " (" << dist.models[i].m << "," << dist.models[i].v << "," << dist.models[i].w << ")";
  os << ">";
  return os;
}

///
/// \brief Return distribution information
/// \param x a distribution
///
__host__ __device__
template <int GMs>
inline std::string getName(const GaussianMixture<GMs> &x) {
  std::stringstream ss;
  ss << "GaussianMixture" << GMs;
  return ss.str();
}

// ------------------------------------------------------------------------------
// Below defines Gaussian related arithmetics

///
/// \brief random variable with unary -
/// \param x a distribution
///
template <int GMs>
inline GaussianMixture<GMs> operator-(const GaussianMixture<GMs> &x)
{
  throw NotImplementedException();

}

///
/// \brief random variable +=
/// \param x distribution1
/// \param rhs distribution2
///
template <int GMs>
inline GaussianMixture<GMs>& operator+=(GaussianMixture<GMs> &x, const GaussianMixture<GMs>& rhs) {
  throw NotImplementedException();
}

///
/// \brief random variable += with scalar
/// \param x distribution1
/// \param r a value which is applied to means of GMM
///
template <int GMs>
inline GaussianMixture<GMs>& operator+=(GaussianMixture<GMs> &x, const double r) {
  for (int i=0; i<GMs; i++)
    x.models[i].m += r;
  return x;
}

///
/// \brief random variable *= with scalar
/// \param x distribution1
/// \param r a value which is applied to means and variances of GMM
///
template <int GMs>
inline GaussianMixture<GMs>& operator*=(GaussianMixture<GMs> &x, const double r) {
  for (int i=0; i<GMs; i++) {
    x.models[i].m *= r;
    x.models[i].v *= r;
  }
  return x;
}

typedef GaussianMixture<MAX_GMs> DefaultGaussianMixture;

}  // namespace dist
}  // namespace edda



#endif  // GAUSSIAN_MIXTURE_H
