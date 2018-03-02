// Copyright 2015 The Edda Authors. All rights reserved.
// Use of this source code is governed by a MIT-style license that can be
// found in the LICENSE file.

#ifndef DIST_JOINT_HISTOGRAM_H_
#define DIST_JOINT_HISTOGRAM_H_

#include <cstdlib>
#include <iostream>
#include <vector>
#include <unordered_set>
#define _USE_MATH_DEFINES  // For Visual Studio
#include <math.h>

#include <boost/limits.hpp>
#include <boost/assert.hpp>
#include <boost/static_assert.hpp>
#include <boost/bind.hpp>

#include "common.h"
#include "distribution_tag.h"
#include "core/statistics.h"
#include <boost/unordered_map.hpp>

namespace edda {
namespace dist {

/// \brief Define a Joint histogram class

/// Joint histogram is used for datasets containing multiple variables.
/// Inferences, such as marginalization, conditional probability, can be 
/// done with a joint histogram.
struct EDDA_EXPORT JointHistogram: public DiscreteDistributionTag, public JointDistributionTag {
	ublas_vector mean; ///<  joint mean of the distribution, a boost's vector
	
	/// \brief Default constructor
	__host__ __device__
	JointHistogram() {
		// initialize a default joint histogram for 3 variables
		num_vars = 3;
	
		min_vals = std::vector<Real>(num_vars, 0);
		max_vals = std::vector<Real>(num_vars, 0);
		num_bins = std::vector<int>(num_vars, 0);
		mean = ublas::zero_vector<Real>(3);
		setCovariance( ublas::identity_matrix<Real>(3) );
	}

	/// \brief Constructor with data and distribution parameters
	/// \param dataAry the input data which consists of multiple variables
	/// \param nElements number of data elements
	/// \param mins an array specifies the minimum values of each variable
	/// \param maxs an array specifies the maximum values of each variable
	/// \param nBins an array specifies the number of bins for each 
	///	variable in the joint histogram
	__host__ __device__
	JointHistogram(std::vector<Real*>& dataAry, int nElements, const std::vector<Real> &mins, 
					const std::vector<Real> &maxs, const std::vector<int> &nBins) {
		assert(dataAry.size()==mins.size() && mins.size()==maxs.size() && mins.size()==nBins.size());
		for(int i=0; i<mins.size(); i++) assert(maxs[i]>=mins[i]);

		// initialize a joint histogram from given data
		num_vars = dataAry.size();
		min_vals = mins;
		max_vals = maxs;
		num_bins = nBins;
		// compute the pdf, mean and cov 
		_update(dataAry, nElements);
	}

	/// \brief Constructor that constructs a given jontHistogram parameters 
	/// \param nVars number of components/variables
	///	\param mins an array specifies the minimum values of each variable
	/// \param maxs an array specifies the maximum values of each variable
	/// \param binWidths an array specifies the bin width in each variable
	///	\param nBins an array specifies the number of bins in each variable
	/// \param new_pdf the pdf of the given jointHistogram
	/// \param new_mean the joint mean of the given jointHistogram
	/// \param new_cov the covariance matrix of the given jointHistogram 
	__host__ __device__
	JointHistogram(int nVars, std::vector<Real> mins, std::vector<Real> maxs, std::vector<Real> binWidths, 
					std::vector<int> nBins, boost::unordered_map<std::vector<int>, Real> new_pdf, 
					ublas_vector new_mean, ublas_matrix new_cov) {
		// construct from a given JointHistogram parameters
		num_vars = nVars;
		min_vals = mins;
		max_vals = maxs;
		num_bins = nBins;
		bin_widths = binWidths;
		pdf = new_pdf;
		mean = new_mean;
		setCovariance(new_cov);
		_computeJointCdf();
	}

	/// \brief Set the min value of each variable
	/// \param mins each element of the vector represents the min 
	/// of the corresponding variable
	__host__ __device__
	void setMinVals(const std::vector<Real>& mins) { min_vals = mins; } 

	/// \brief Get the min value of each variable
	/// \return each element of the vector represents the min 
	/// of the corresponding variable	
	__host__ __device__
	std::vector<Real> getMinVals() const { return  min_vals; } 
	
	/// \brief Set the max value of each variable
	/// \param maxs each element of the vector represents 
	/// the max of the corresponding variable
	__host__ __device__
	void setMaxVals(const std::vector<Real>& maxs) { max_vals = maxs; } 
	
	/// \brief Get the max value of each variable
	/// \return each element of the vector represents 
	/// the max of the corresponding variable
	__host__ __device__
	std::vector<Real> getMaxVals() const { return  max_vals; } 

	/// \brief Set the bin width used for each variable
	/// \param widths each element of the vector represents 
	/// the bin width of the corresponding variable
	__host__ __device__
	void setBinWidths(const std::vector<Real>& widths) { bin_widths = widths; } 

	/// \brief Get the bin width used for each variable
	/// \return each element of the vector represents 
	/// the bin width of the corresponding variable
	__host__ __device__
	std::vector<Real> getBinWidths() const { return  bin_widths; } 

	/// \brief Set the number of bins for each variable
	/// \param bins each element of the vector represents 
	/// the number of bins for the corresponding variable
	__host__ __device__
	void setNumBins(const std::vector<int>& bins) { num_bins = bins; } 

	/// \brief Get the number of bins for each variable
	/// \return each element of the vector represents 
	/// the number of bins for the corresponding variable
	__host__ __device__
	std::vector<int> getNumBins() const { return  num_bins; } 

	/// \brief Set the number of variables of the joint histogram 
	/// \param nVars the number of variables 
	__host__ __device__
	void setNumVars(const int nVars) { num_vars = nVars; } 
	
	/// \brief Get the number of variables of the joint histogram 
	/// \return the number of variables 
	__host__ __device__
	int getNumVars() const { return  num_vars; } 

	/// \brief Set the distribution (pdf) of the joint histogram 
	/// \param p an unordered_map that stores a pdf. In the map structure,
	/// the key (std:vector<int>) is the index of a non-empty bin 
	/// of the joint histogram; the value is the probability of that bin.
	void setDistr(boost::unordered_map<std::vector<int>, Real> p) { pdf = p; }
	
	/// \brief Get the distribution (pdf) of the joint histogram 
	/// \return un unordered_map that stores a pdf. In the map structure,
	/// the key (std:vector<int>) is the index of a non-empty bin 
	/// of the joint histogram; the value is the probability of that bin.
	boost::unordered_map<std::vector<int>, Real> getDistr() const { return pdf; }

	/// \brief Compute a cumulative distribution function, which is similar to the CDF,
	/// but only for the non-empty bins and these non-empty bins are sorted in
	/// decreasing order (based on their probability) before computing this CDF. 
	/// \return a sparse cumulative distribution function, which contains
	/// only non-empty bins' cumulative distributions. Note that this is not the
	///	real CDF for the distribution, as the order of the bins	was sorted in 
	/// decreasing order first. This function is only used for sampling.
	std::vector<std::pair<std::vector<int>, Real>> getJointCdf() const { 
		assert(cdf.size()>0); 
		assert(cdf.size()==pdf.size()); 
		return cdf; 
	}
	
	/// \brief Set the covariance matrix of the joint histogram 
	/// \param cov the covariance matrix,
	__host__ __device__
	void setCovariance(const ublas_matrix &cov) {
		this->cov = cov;
		this->det = determinant(this->cov);
	}
	
	/// \brief Get the covariance matrix of the joint histogram 
	/// \return the covariance matrix,
	__host__ __device__
	const ublas_matrix &getCovariance() const {  return this->cov ; }

	/// \brief Get the determinant of the covariance matrix 
	/// \return the determinant of the covariance matrix,
	__host__ __device__
	double getDet() const {  return this->det ; }

	/// \brief Marginalization of the jointHistogram
	/// \param vars variables want to be left after marginalization, 
	///	i.e. project onto which variables
	/// \return the marginalized distribution, which is also a jointHistogram
	JointHistogram marginalization(std::unordered_set<int>& vars) const {
		assert(pdf.size()>0);
		assert(vars.size()>0);

		std::vector<int> new_num_bins;
		std::vector<Real> new_bin_widths;
		std::vector<Real> new_min_vals;
		std::vector<Real> new_max_vals;
		boost::unordered_map<std::vector<int>, Real> new_pdf;
		ublas_vector new_mean = ublas::zero_vector<Real>(vars.size()); // boost's vector
		ublas_matrix new_cov(ublas::identity_matrix<Real>(vars.size())); // boost's vector
		int varid = 0;
		for(auto itv=vars.begin(); itv!=vars.end(); ++itv) {
			new_num_bins.push_back(num_bins[*itv]);
			new_bin_widths.push_back(bin_widths[*itv]);
			new_min_vals.push_back(min_vals[*itv]);
			new_max_vals.push_back(max_vals[*itv]);
			new_mean[varid] = mean[*itv];
			varid++;
		}
		// marginalize the pdf
		for (auto itp=pdf.begin(); itp!=pdf.end(); ++itp) {
			std::vector<int> bin = itp->first;
			std::vector<int> key;
			for(auto itv=vars.begin(); itv!=vars.end(); ++itv) {
				key.push_back(bin[*itv]);
			}
			new_pdf[key] += itp->second;
		}
		// compute new mean and cov
		int iitr=0;	
		for(auto itv=vars.begin(); itv!=vars.end(); ++itv) {
			int jitr = 0;
			for(auto jtv=vars.begin(); jtv!=vars.end(); ++jtv) {
				new_cov(iitr,jitr) = cov(*itv, *jtv);
				jitr++;
			}
			iitr++;
		}
		JointHistogram jhist(vars.size(), new_min_vals, new_max_vals, new_bin_widths, new_num_bins, new_pdf, new_mean, new_cov);
		return jhist; 
	}
	
	/// \brief Compute the joint histogram of certain variables (vars) given the 
	/// rest of variables(cond_var) within certain ranges (bin_range)
	/// \param vars compute the joint histogram of these variables
	/// \param cond_var given the range of these variables
	/// \param bin_range the bin range of conditional variables (cond_var)
	///	\return return the conditional distribution, which is also a jointHistogram
	JointHistogram conditionalHist(std::unordered_set<int>& vars, std::vector<int>& cond_var, std::vector<std::pair<int,int>>& bin_range) { 
		// step 0: check input parameter
		assert(pdf.size()>0);
		assert(vars.size()>0);
		assert(cond_var.size()>0 && cond_var.size()==bin_range.size());
		for(int i=0; i<bin_range.size(); i++) {// validate bin range
			assert(bin_range[i].first<=bin_range[i].second);
		}
		for(int i=0; i<cond_var.size(); i++) {// check for duplication
			for(int j=i+1; j<cond_var.size(); j++) {
				assert(cond_var[i]!=cond_var[j]);
			}
		}
		// step 1: computing P(cond_var in bin_range)
		Real pB = 0;
		bool flag_valid = false;
		for (auto itp=pdf.begin(); itp!=pdf.end(); ++itp) {
			std::vector<int> bin = itp->first;
			bool flag_in_range = true;
			for(int i=0; i<cond_var.size(); i++) {
				if(bin[cond_var[i]]<bin_range[i].first || bin[cond_var[i]]>bin_range[i].second) {
					flag_in_range = false;
					break;
				}
			}
			if(flag_in_range){
				flag_valid = true;
				pB += itp->second;
			}
		}
		if(pB==0 || !flag_valid) {// no point in conditional range
			throw std::runtime_error("no data item falls in the given bin ranges!");
		}
	
		// step 2: computing P(vars AND cond_vars)/P(cond_var in bin_range) 
		boost::unordered_map<std::vector<int>, Real> new_pdf;
		for (auto itp=pdf.begin(); itp!=pdf.end(); ++itp) {
			std::vector<int> bin = itp->first;
			std::vector<int> key;
			for(auto itv=vars.begin(); itv!=vars.end(); ++itv) {
				key.push_back(bin[*itv]);
			}
			bool flag_in_range = true;
			for(int i=0; i<cond_var.size(); i++) {
				if(bin[cond_var[i]]<bin_range[i].first || bin[cond_var[i]]>bin_range[i].second) {
					flag_in_range = false;
					break;
				}
			}
			if(flag_in_range)
				new_pdf[key] += itp->second/pB;
		}
		
		// step 3: construct a new JointHistogram
		std::vector<int> new_num_bins;
		std::vector<Real> new_bin_widths;
		std::vector<Real> new_min_vals;
		std::vector<Real> new_max_vals;
		ublas_vector new_mean = ublas::zero_vector<Real>(vars.size()); // boost's vector
		ublas_matrix new_cov(ublas::identity_matrix<Real>(vars.size())); // boost's vector
		int varid = 0;
		for(auto itv=vars.begin(); itv!=vars.end(); ++itv) {
			new_num_bins.push_back(num_bins[*itv]);
			new_bin_widths.push_back(bin_widths[*itv]);
			new_min_vals.push_back(min_vals[*itv]);
			new_max_vals.push_back(max_vals[*itv]);
			new_mean[varid] = mean[*itv];
			varid++;
		}
		int iitr=0;	
		for(auto itv=vars.begin(); itv!=vars.end(); ++itv) {
			int jitr = 0;
			for(auto jtv=vars.begin(); jtv!=vars.end(); ++jtv) {
				new_cov(iitr,jitr) = cov(*itv, *jtv);
				jitr++;
			}
			iitr++;
		}
		JointHistogram jhist(vars.size(), new_min_vals, new_max_vals, new_bin_widths, new_num_bins, new_pdf, new_mean, new_cov);
		return jhist; 
	}
	
private:
	int num_vars;					// number of variables
	std::vector<Real> min_vals;		// min of each variable
	std::vector<Real> max_vals;		// max of each variable
	std::vector<int> num_bins;		// number of bins for each variable
	std::vector<Real> bin_widths;	// bin width of each variable
	boost::unordered_map<std::vector<int>, Real> pdf;	// probability distribution function
	std::vector<std::pair<std::vector<int>, Real>> cdf;	// cumulative distribution function
	ublas_matrix cov;				// covariance matrix 
	double det;						// determinant of the covariance matrix

	/// \brief Compute joint PDF, mean, and cov 
	/// \param dataAry raw multivariate array
	/// \param nElements number of elements in each variable
	__host__ __device__
	void _update(std::vector<Real*>& dataAry, int nElements) {
		// 1. compute the PDF (i.e. a sparse matrix)
		boost::unordered_map<std::vector<int>, int> cntMap;
		std::vector<int> key(num_vars);
		bin_widths.clear();	
		for(int i=0; i<num_vars; i++)	
			bin_widths.push_back((max_vals[i]-min_vals[i])/num_bins[i]);
		for(int j=0; j<nElements; j++) {
			for(int i=0; i<num_vars; i++) {// value to bin index
				key[i] = std::floor((dataAry[i][j]-min_vals[i])/bin_widths[i]);
				if(key[i]>=0 && key[i]<num_bins[i])
					continue;	
				else
					key[i] = num_bins[i]-1;
			}
			cntMap[key]++;// rely on the default initialization
		}
		// convert count to frequency
		for (auto it=cntMap.begin(); it!=cntMap.end(); ++it) {
			pdf[it->first] = it->second*1.0/nElements;	
		}
		_computeJointCdf();

		// 2. update the mean 
		mean = ublas::zero_vector<Real>(num_vars);
		for (auto it=pdf.begin(); it!=pdf.end(); ++it) {
			std::vector<int> bin = it->first;
			for(int i=0; i<bin.size(); i++) {
				mean[i] += (min_vals[i]+(bin[i]+0.5)*bin_widths[i])*it->second; 
			}
		}

		// 3. update the cov matrix
		setCovariance( ublas::identity_matrix<Real>(num_vars) );
		for(int i=0; i<num_vars; i++) {
			for(int j=0; j<num_vars; j++) {
				for(int k=0; k<nElements; k++) {
					cov(i,j) += (dataAry[i][k]-mean[i])*(dataAry[j][k]-mean[j]);	
				}
			}
		}
		cov *= 1.0/nElements;
		setCovariance(cov);
	}
	
	/// \brief Compute the pseudo CDF, and store into the member cdf
	__host__ __device__
	void _computeJointCdf() {
		assert(pdf.size()>0); 
		if(cdf.size()==pdf.size()) return; // already constructed before
		
		// sort PDF in decreasing order, i.e. move high frequency bins to front
		std::vector<std::pair<std::vector<int>, Real>> sorted_pdf(pdf.begin(), pdf.end());
		std::sort(sorted_pdf.begin(), sorted_pdf.end(), 
			boost::bind(&std::pair<std::vector<int>, Real>::second, _1) > 
			boost::bind(&std::pair<std::vector<int>, Real>::second, _2));
	
		// construct a CDF on the sparse JointHistogram
		Real accFreq = sorted_pdf[0].second;
		cdf.push_back(sorted_pdf[0]);	
		for(int i=1; i<sorted_pdf.size(); i++) {
			accFreq += sorted_pdf[i].second;
			std::pair<std::vector<int>, Real> tmp(sorted_pdf[i].first, accFreq);
			cdf.push_back(tmp);
		}
		//const float epsilon = std::numeric_limits<float>::epsilon();	
		//assert(abs(accFreq-1)<epsilon);
	}
};

// ------------------------------------------------------------------------------
// Below defines JointHistogram related generic functions
///
/// \brief Return joint mean of the joint histogram
/// \param dist the joint histogram
__host__ __device__
inline std::vector<Real> getJointMean(const JointHistogram &dist) {
	std::vector<Real> m(dist.mean.size());
	std::copy(dist.mean.begin(), dist.mean.end(), m.begin());
	return m;
}

/// \brief Return the probability of a multivariate value 
///	given a JointHistogram object (distribution)
///	\param dist the JointHistogram distribution
/// \param x the multivariate value 
/// \return the probability of the given value 
//TODO: need to change the function name, getJointFrequency?
__host__ __device__
inline double getJointPdf(const JointHistogram &dist, const std::vector<Real> &x_) { 
	int nvars = dist.getNumVars();
	std::vector<Real> mins = dist.getMinVals();
	std::vector<Real> widths = dist.getBinWidths();
	std::vector<int> bins = dist.getNumBins();
	std::vector<int> key(nvars);
	for(int i=0; i<nvars; i++) {// value to bin index
		key[i] = std::floor((x_[i]-mins[i])/widths[i]);
		if(key[i]>=0 && key[i]<bins[i])
			continue;	
		else
			key[i] = bins[i]-1;
	}
	return (dist.getDistr())[key]; 
}

/// \brief Using binary search to find the bin that a particular 
/// value is in. This function is used exclusively in getJointSample(). 
///	\param vec a sorted vector of pairs according to the second 
///	value of each pair. This vector is the cumulative density function
/// based on the increasing order of the pdf. 
/// \param p the searched value
/// \return the index of the element of vec that value p falls into
__host__
static inline int biSearchNextLarge(std::vector<std::pair<std::vector<int>, Real>>& vec, Real p) {
	assert(vec.size()>0);
	if(vec.size()==1) return 0;
	// vec is sorted, increasing order
	int low = 0;
	int high = vec.size()-1;

	while(low<=high) {
		int mid = low+ (high-low)/2;
		if(vec[mid].second<p) low = mid+1;
		else if(vec[mid].second>p) high = mid-1;
		else return mid+1;
	}

	if(high<0) return 0; //p<vec[0];
	else{
		if( low > (vec.size()-1))
			return vec.size()-1; // p>=vec[vec.size()-1]
		else
			return (low<high) ? low+1 : high+1;
	}
}

/// \brief Generate a random sample from a distribution
///	\param dist the JointHistogram object (the distribution)
///	\return a multivariate sample 
__host__
inline std::vector<Real> getJointSample(const JointHistogram &dist) {
	if(dist.getDistr().size()==0) {
		printf("ERROR: distribution size is 0.\n");
		return std::vector<Real>(dist.getNumVars(), 0);
	}
	std::vector<std::pair<std::vector<int>, Real>> cdf;
	cdf = dist.getJointCdf();
	Real x = static_cast<Real>(rand())/RAND_MAX;
	int i=0;
	bool useBiSearch = true;
	if(useBiSearch) {
		i = biSearchNextLarge(cdf, x);
	}
	else {
		while(x > cdf[i].second)	
			i++;
	}

	// if reach here, we find the bin
	std::vector<int> bin = cdf[i].first;
	std::vector<Real> smp;
	std::vector<Real> mins = dist.getMinVals();
	std::vector<Real> widths = dist.getBinWidths();
	bool binCenter = true;
	for(int i=0; i<bin.size(); i++) { 
		Real smpPos;
		if(binCenter) smpPos = 0.5;// use bin center
		else smpPos = static_cast<Real>(rand())/RAND_MAX;
		smp.push_back(mins[i]+(bin[i]+smpPos)*widths[i]);
	}
	return smp;
}

/// \brief Return a random sample using random engine. This function 
/// leaves the interface for parallel platforms
/// \param dist the JointHistogram object
/// \param rng random engine from thrust
/// \return a random float number
__host__ __device__
inline std::vector<Real> getJointSample(const JointHistogram &dist, thrust::default_random_engine &rng) {
  // TODO version used on device side
	return getJointSample(dist);
}

/// \brief Print the joint histogram itself
///	\param os destination stream 
///	\dist the joint histogram object
/// \return destination stream after inserting information 
///	of the JointHistogram object
__host__
inline std::ostream& operator<<(std::ostream& os, const JointHistogram &dist) {
    os <<  "<JointHistogram: mean=" << dist.mean << ", covariance=" << dist.getCovariance() << ">" ;
    return os;
}

/// \brief Return the distribution type in string, i.e. "JointHistogram"
///	\dist the joint histogram object
/// \return the type of the object in string, i.e. "JointHistogram"
__host__ __device__
inline std::string getName(const JointHistogram &dist) {
    return "JointHistogram";
}

/// \brief This function is the interface used to construct 
///	a JointHistogram object outside of the class
///	\param dataAry the raw multivariate data array. Each vector 
///	element is a pointer pointing to the start location of a 
///	block of memory that stores one variable
/// \param nElements number of elements in each variable. All variables
///	are assumed to have the same number of elements
///	\param mins minimum value of each variable
///	\param maxs maximum value of each variable
/// \param nBins number of bins in each variable
///	\return the constructed JointHistogram object
__host__ __device__
inline JointHistogram eddaComputeJointHistogram(std::vector<Real*>& dataAry, int nElements, const std::vector<Real> &mins, 
				const std::vector<Real> &maxs, const std::vector<int> &nBins) {
	JointHistogram tmp(dataAry, nElements, mins, maxs, nBins); 
	return tmp;
}

}  // namespace dist
}  // namespace edda

#endif  // DIST_HISTOGRAM_H_
