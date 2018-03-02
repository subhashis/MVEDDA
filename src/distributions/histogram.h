// Copyright 2015 The Edda Authors. All rights reserved.
// Use of this source code is governed by a MIT-style license that can be
// found in the LICENSE file.

#ifndef HISTOGRAM_H
#define HISTOGRAM_H

/// Histogram: TODO
///

#include <iostream>
#include <typeinfo>
#include <iostream>

#include "common.h"
#include "edda_export.h"
#include "core/vector_matrix.h"
#include "core/tuple.h"
#include "core/thrust_common.h"
#include "distribution_tag.h"
#include "core/shared_ary.h"

namespace edda {
namespace dist {

///
/// \brief The Distribution class is a root class for all distribution-type classes.
///
class EDDA_EXPORT Histogram : public DiscreteDistributionTag {
public:

  ///
  /// \brief Empty constructor
  ///
  Histogram(){}

  ///
  /// \brief Constructor
  /// \param hist input histogram array
  ///
  Histogram(double* histData){//use when load from vtk
	  m_nBins = histData[0];
	  m_minValue = histData[1];
	  m_maxValue = histData[2];
	  m_binWidth = (m_maxValue - m_minValue) / (Real)(m_nBins);

	  m_cdf.resize(m_nBins);
	  for (int b = 0; b < m_nBins; b++){
		  m_cdf[b] = histData[b+3];
	  }
  }

  ///
  /// \brief Constructor
  /// \param hist input histogram array
  /// \param _nBins number of bins
  ///
  Histogram(float* histData, int _nBins){//use when load from .edda file
	  m_nBins = _nBins;
	  m_minValue = histData[0];
	  m_maxValue = histData[1];
	  m_binWidth = (m_maxValue - m_minValue) / (Real)(m_nBins);

	  m_cdf.resize(m_nBins);
	  for (int b = 0; b < m_nBins; b++){
		  m_cdf[b] = histData[b + 2];
	  }
  }

  ///
  /// \brief Constructor by given min max value of the histogram
  /// \param dataArray contain samples to construct a hisotgram
  /// \param nElements number of samples in dataArray
  /// \param _nBins number of bins
  /// \param _minValue data range (min)
  /// \param _maxValue data range (max)
  ///
  Histogram(Real *dataArray, int nElements, const int _nBins, const Real _minValue = 1, const Real _maxValue = -1){
	  if (_minValue > _maxValue){// no input min/max value
		  m_minValue = dataArray[0];
		  m_maxValue = dataArray[0];
		  for (int i = 1; i < nElements; i++){
			  if (m_minValue > dataArray[i])
				  m_minValue = dataArray[i];
			  if (m_maxValue < dataArray[i])
				  m_maxValue = dataArray[i];
		  }
	  }
	  else{
		  m_minValue = _minValue;
		  m_maxValue = _maxValue;
	  }

	  m_nBins = _nBins;

	  m_binWidth = (m_maxValue - m_minValue) / (Real)(m_nBins);

	  m_cdf.resize(m_nBins);

	  //modeling and convert to cdf
	  for (int i = 0; i < nElements; i++){
		  int b = valueToBinsIndex(dataArray[i]);
		  m_cdf[b] ++;
	  }

	  for (int b = 1; b < m_nBins; b++)
		  m_cdf[b] += m_cdf[b - 1];

	  for (int b = 0; b < m_nBins; b++){
		  m_cdf[b] /= (Real)nElements;
	  }
  }

  ///
  /// \brief Constructor by given a histogram
  /// \param hist the given histogram
  ///
  Histogram(const Histogram &hist)
  : m_nBins(hist.m_nBins), 
    m_minValue(hist.m_minValue),
    m_maxValue(hist.m_maxValue),
    m_binWidth(hist.m_binWidth)
  {
    m_cdf = hist.m_cdf;
  }

  ///
  /// \brief Return mean of histogram
  ///
  Real getMean() const{
    Real mean = 0;
    Real cValue = m_minValue + m_binWidth / 2.0;
    for (int i = 0; i < m_nBins; i++){
      if (i == 0) mean += m_cdf[0] * cValue;
      else mean += (m_cdf[i] - m_cdf[i - 1]) * cValue;
      cValue += m_binWidth;
    }
    return mean;
  }

  ///
  /// \brief Return variance of the histogram
  ///
  Real getVar() const{
    //using Var(X) = E[X^2] - (E[X])^2  to compute Var
    Real mean = 0;
    Real mSqrt = 0;
    Real cValue = m_minValue + m_binWidth / 2.0;
    for (int i = 0; i < m_nBins; i++){
      if (i == 0) mean += m_cdf[0] * cValue;
      else mean += (m_cdf[i] - m_cdf[i - 1]) * cValue;

      if (i == 0) mSqrt += m_cdf[0] * cValue * cValue;
      else mSqrt += (m_cdf[i] - m_cdf[i - 1]) * cValue * cValue;

      cValue += m_binWidth;
    }
    return mSqrt - mean*mean;
  }

  ///
  /// \brief Return probability of a given value
  /// \param x given value
  ///
  Real getPdf(const double x) const{
    int b = valueToBinsIndex(x);

    if (b == 0) return m_cdf[0];
    else return (m_cdf[b] - m_cdf[b - 1]);
  }

  ///
  /// \brief Return cumulative probability of a given value
  /// \param x given value
  ///
  Real getCdf(const double x) const{
    int b = valueToBinsIndex(x);

    return m_cdf[b];
  }

  ///
  /// \brief Return a sample dran from a hisotgram
  ///
  Real getSample() const{
    Real sample;
    Real r = rand() / (float)RAND_MAX;

    int low = 0;
    int high = m_nBins - 1;

    //bineary searh in cdf
    //do{
    //  int mid = floor((low + high) / 2.0);

    //  Real prevCdf;
    //  if (mid == 0) prevCdf = -0.1;
    //  else prevCdf = m_cdf[mid - 1];

    //  if (prevCdf < r && r <= m_cdf[mid]){
    //    //sample here
    //    Real subr = rand() / (float)RAND_MAX;
    //    Real s, t;
    //    rangeOfBin(s, t, mid);
    //    sample = s + m_binWidth * subr;
    //    break;
    //  }
    //  else if (m_cdf[mid] > r){
    //    high = mid - 1;
    //  }
    //  else{
    //    low = mid + 1;
    //  }
    //} while (1);

	//sequence search
	for (int b = 0; b < m_nBins; b++){
		if (m_cdf[b] >= r){
			Real subr = rand() / (float)RAND_MAX;
			Real s, t;
			rangeOfBin(s, t, b);
			sample = s + m_binWidth * subr;
			break;
		}
	}

    return sample;
  }


  Real getInverseCdf(const double x) const{
    Real sample;
    //Real r = rand() / (float)RAND_MAX;

    int low = 0;
    int high = m_nBins - 1;

    //sequence search
    for (int b = 0; b < m_nBins; b++){
      if (m_cdf[b] >= x){
        Real subr = rand() / (float)RAND_MAX;
        Real s, t;
        rangeOfBin(s, t, b);
        sample = s + m_binWidth * subr;
        break;
      }
    }

    return sample;

    //std::cout << "her you go! x = " << x << std::endl;
    //return 1.0;
  }

  ///
  /// \brief Output histogram information
  /// \param os outstream
  ///
  void output(std::ostream& os) const{
    os << "<Histogram Min:" << this->m_minValue << " Max:" << this->m_maxValue;
    for (int b = 0; b < m_nBins; b++)
      os << ", Bin " << b << ": " << m_cdf[b];
    os << ">";
  }

  ///
  /// \brief Return number of bins
  ///
  int getBins(){
	  return m_nBins;
  }

  ///
  /// \brief Return max value of data value range
  ///
  float getMaxValue(){
	  return m_maxValue;
  }

  ///
  /// \brief Return min value of data value range
  ///
  float getMinValue(){
	  return m_minValue;
  }

  ///
  /// \brief Return get probability of a bin
  /// \param b bin ID
  ///
  float getBinValue(int b){
	  //This usually return accumlative prob
	  return m_cdf[b];
  }

protected:
  ///
  /// \brief Return bin ID of a given value
  /// \param v given scalar value
  ///
  int valueToBinsIndex(Real v) const {
	if (v < m_minValue)
		  return 0;
	else if (v >= m_maxValue)
		  return (m_nBins - 1);
		else
		  return floor((v - m_minValue) / m_binWidth);
  }

  ///
  /// \brief Return data value range of a bin
  /// \param s Return data value low bound of the bin
  /// \param t Return data value upper bound of the bin
  /// \param b given bin ID
  ///
  void rangeOfBin(Real& s, Real &t, int b) const{
    s = m_minValue + b * m_binWidth;
    t = m_minValue + (b + 1) * m_binWidth;
  }

  int m_nBins;      // number of bins
  Real m_minValue;  // min of data value range 
  Real m_maxValue;  // max of data value range 
  Real m_binWidth;  // bin width of bins
  std::vector<Real> m_cdf;  // histogra's cdf vector
};

///
/// \brief Compute the mean of the distribution
/// \param dist a distribution (histogram)
///
inline double getMean(const Histogram &dist)
{
  return (double)dist.getMean();
}

///
/// \brief Compute Variance
/// \param dist a distribution (histogram)
///
inline double getVar(const Histogram &dist)
{
  return (double)dist.getVar();
}

///
/// \brief Return PDF of x
/// \param dist a distribution (histogram)
/// \param x a input scalar value 
///
inline double getPdf(const Histogram &dist, const double x)
{
  return (double)dist.getPdf(x);
}

///
/// \brief Get a Monte-Carlo sample of the distribution. We rely on the specific implementation from each distribution.
/// \param dist a distribution (histogram)
///
inline double getSample(const Histogram &dist)
{
  return (double)dist.getSample();
}


///
/// \brief Return CDF of x
/// \param dist a distribution (histogram)
/// \param x a input scalar value
///
//__host__ __device__
inline double getCdf(const Histogram &dist, double x)
{
  return (double)dist.getCdf(x);
}

///
/// \brief Print histogram information
/// \param os outstream
/// \param dist a distribution (histogram)
///
inline std::ostream& operator<<(std::ostream& os, const Histogram &dist)
{
  dist.output(os);
  return os;
}

///
/// \brief Print distribution type
/// \param os outstream
///
__host__ __device__
inline std::string getName(const Histogram &x) {
  return "Histogram";
}

}  // namespace dist

///
/// \brief Compute a histogram
/// \param dataArray contain samples to construct a hisotgram
/// \param nElements number of samples in dataArray
/// \param _nBins number of bins
/// \param _minValue data range (min)
/// \param _maxValue data range (max)
///
inline dist::Histogram eddaComputeHistogram(float *dataArray, int nElement, const int _nBins, const float _minValue = 1, const float _maxValue = -1)
{
	//convert to internal data type
	Real* _dataArray = (Real*)malloc(sizeof(Real)*nElement);
	for (int i = 0; i < nElement; i++)_dataArray[i] = (Real)dataArray[i];
	dist::Histogram h = dist::Histogram(_dataArray, nElement, _nBins, _minValue, _maxValue);
	free(_dataArray);
	return h;
}

///
/// \brief Compute a histogram
/// \param dataArray contain samples to construct a hisotgram
/// \param nElements number of samples in dataArray
/// \param _nBins number of bins
/// \param _minValue data range (min)
/// \param _maxValue data range (max)
///
inline dist::Histogram eddaComputeHistogram(double *dataArray, int nElement, const int _nBins, const double _minValue = 1, const double _maxValue = -1)
{
	//convert to internal data type
	Real* _dataArray = (Real*)malloc(sizeof(Real)*nElement);
	for (int i = 0; i < nElement; i++)_dataArray[i] = (Real)dataArray[i];
	dist::Histogram h = dist::Histogram(_dataArray, nElement, _nBins, _minValue, _maxValue);
	free(_dataArray);
	return h;
}

}  // namespace edda

#endif // HISTOGRAM_H
