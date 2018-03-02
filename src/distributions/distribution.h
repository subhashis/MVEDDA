// Copyright 2015 The Edda Authors. All rights reserved.
// Use of this source code is governed by a MIT-style license that can be
// found in the LICENSE file.

#ifndef DISTRIBUTION_H_
#define DISTRIBUTION_H_

#include <iostream>
#include <string>

#include "gaussian.h"
#include "gaussian_mixture.h"
#include "histogram.h"

/// Distributions:
/// Define distribution interfaces
///

namespace edda{
namespace dist{

///
/// \brief The Distribution class is an interface class with virtual classes
///
class Distribution
{
public:
  virtual double getPdf(double x)=0;
  virtual double getCdf(double x)=0;
  virtual double getMean()=0;
  virtual double getSample()=0;
  virtual std::string getName() = 0;
  virtual std::ostream& operator<<(std::ostream& os)=0;
};

template <class T>
class DistributionWrapper: public Distribution
{
public:
  T dist;

  DistributionWrapper (const T &dist_) :dist(dist_) { }

  virtual double getPdf(double x) {
    return dist::getPdf(dist, x);
  }

  virtual double getCdf(double x) {
    return dist::getCdf(dist, x);
  }

  virtual double getMean() {
    return dist::getMean(dist);
  }

  virtual double getSample() {
    return dist::getSample(dist);
  }

  virtual std::ostream& operator<<(std::ostream& os) {
    os << dist;
    return os;
  }

  virtual std::string getName() {
    return dist::getName(dist);
  }
};

typedef DistributionWrapper<Gaussian> GaussianWrapper;
typedef DistributionWrapper<DefaultGaussianMixture> GaussianMixtureWrapper;
typedef DistributionWrapper<Real> RealValueWrapper;
//typedef VirtualDistribution<Histogram> HistogramWrapper;

} // dist
} // edda
#endif // DISTRIBUTION_H_
