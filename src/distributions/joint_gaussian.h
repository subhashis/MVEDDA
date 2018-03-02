// Copyright 2015 The Edda Authors. All rights reserved.
// Use of this source code is governed by a MIT-style license that can be
// found in the LICENSE file.

#ifndef DIST_JOINT_GAUSSIAN_H_
#define DIST_JOINT_GAUSSIAN_H_

#include <cstdlib>
#include <iostream>
#include <vector>
#define _USE_MATH_DEFINES  // For Visual Studio
#include <math.h>

#include <boost/limits.hpp>
#include <boost/assert.hpp>
#include <boost/static_assert.hpp>
#include <boost/random.hpp>

#include "common.h"
#include "distribution_tag.h"
#include "core/statistics.h"
#include "invert_matrix.h"

#include "Eigen/Dense"

using namespace std;
using namespace Eigen;

namespace edda {
	namespace dist {

		// ------------------------------------------------------------------------------
		///
		/// \brief Defines a Gaussian class
		///
		struct EDDA_EXPORT JointGaussian : public ContinuousDistributionTag, public JointDistributionTag {

			// constructor
			__host__ __device__
				JointGaussian() {
					mean = ublas::zero_vector<Real>(3);
					setMatrices(ublas::identity_matrix<Real>(3));
				}

			///
			/// \brief Constructor
			///
			__host__ __device__
				JointGaussian(const ublas_vector &mean, const ublas_matrix &cov) {
					assert(mean.size() == cov.size1() && mean.size() == cov.size2());
					this->mean = mean;
					setMatrices(cov);

				}

			///
			/// \brief set and use eigen decomposition to compute the necessary matrixes, eigenMat: (for sampling) and "uMat" (for probability computation)
			/// \param cov the input covariance matrix
			///
			__host__ __device__
				void setMatrices(const ublas_matrix &cov) {
					this->cov = cov;
	
					int covSize = this->cov.size1(); //cov.size1 == cov.size2

					//BEGIN eigenDecomposition using Eigen Library********************
					MatrixXf rawCov(covSize, covSize);
					for (int i = 0; i < covSize; i++){
						for (int j = 0; j < covSize; j++){
							rawCov(i,j) = this->cov(i, j);
						}
					}
					SelfAdjointEigenSolver<MatrixXf> eigensolver(rawCov);
					if (eigensolver.info() != Success){ std::cout << "eigen decomposition fails" << std::endl; }
					ublas_matrix eigenVectors(covSize, covSize);
					ublas_vector eigenValues(covSize);
					for (int i = 0; i < covSize; i++)eigenValues(i) = eigensolver.eigenvalues()[i];
					for (int i = 0; i < covSize; i++){
						for (int j = 0; j < covSize; j++){
							eigenVectors(i, j) = eigensolver.eigenvectors()(i, j);
						}
					}
					//END eigenDecomposition and eigen value vector is in boost data structure now********************

					//compute eigenMat by eigenvector and eigenvalue (for sampling)
					eigenMat.resize(covSize, covSize);
					for (int j = 0; j < covSize; j++){
						double sqrtEvalue = sqrt(eigenValues(j));
						if (sqrtEvalue != sqrtEvalue)sqrtEvalue = 0; //discard the small negative Eigen value
						column(eigenMat, j) = column(eigenVectors, j) *sqrtEvalue;
					}
					
					//compute the U matrix and logPDet (for probability estimation)
					double cond = 0.0000000000222044604;//just a eps value
					double maxAbsEv = -1;
					for (int i = 0; i < covSize; i++){
						if (maxAbsEv < fabs(eigenValues(i))){
							maxAbsEv = fabs(eigenValues(i));
						}
					}
					double eps = cond * maxAbsEv;

					logPDet = 0;
					for (int i = 0; i < covSize; i++){
						if (eigenValues(i) > eps)
							logPDet += log(eigenValues(i));
					}

					uMat.resize(covSize, covSize);
					for (int j = 0; j < covSize; j++){
						double spinv = eigenValues(j);
						if (spinv > eps) spinv = 1.0 / spinv;
						column(uMat, j) = column(eigenVectors, j) *sqrt(spinv);
					}
				}

			///
			/// \brief Return a sample drawn from this joint Gaussian
			/// \param rng random engine
			///
			__host__ __device__
				std::vector<Real> getJointSample(thrust::default_random_engine &rng) const {
					//draw nVar samples from independent standard normal variant
					thrust::random::normal_distribution<double> ndist(0,1);
					ublas_vector r(cov.size1());
					for (int j = 0; j < r.size(); j++){
						r(j) = ndist(rng);
					}

					//transform the above sample by the covariance matrix
					ublas_vector s(cov.size1());
					boost::numeric::ublas::axpy_prod(eigenMat, r, s, true);
					std::vector<Real> retS(cov.size1());
					for (int j = 0; j < r.size(); j++){
						retS[j] = s(j) + mean(j);
					}

					return retS;
				}

			///
			/// \brief Return log probability of x (This function is used by EM)
			/// \param x_ vector of a sample for probability estimation
			///
			__host__ __device__
				inline double getJointLogPdf(const std::vector<Real> x_) const
			{
					int k = this->mean.size();
					assert(x_.size() == k);
					ublas_vector x(k);
					std::copy(x_.begin(), x_.end(), x.begin());
					x = x - this->mean;
					ublas_vector tmp;
					tmp = ublas::prod(ublas::trans(x), getUMat());
					double density = ublas::inner_prod(tmp, tmp);

					return -0.5 * (k * 1.83787706641 + getLogDet() + density);
			}

			///
			/// \brief Return mean vector of this Gaussian
			///
			__host__ __device__
				const ublas_vector &getMean() const { return this->mean; }

			///
			/// \brief Return covariance matrix of this Gaussian
			///
			__host__ __device__
				const ublas_matrix &getCovariance() const { return this->cov; }

			///
			/// \brief Return eigen matrix of this Gaussian's covariance matrix
			///
			__host__ __device__
				const ublas_matrix &getEigenMat() const { return this->eigenMat; }

			///
			/// \brief Return eigen matrix of this Gaussian's covariance matrix
			///
			__host__ __device__
				const ublas_matrix &getUMat() const { return this->uMat; }

			///
			/// \brief Return U matrix of this Gaussian's covariance matrix
			///
			__host__ __device__
				double getLogDet() const { return this->logPDet; }

		private:
			ublas_vector mean; 		// boost's vector for mean vector
			ublas_matrix cov;		//covariance matrix
			ublas_matrix eigenMat; //Eigen matrix: use when draw sample, sample = eigenMat * z + mean
			ublas_matrix uMat;	   //U-Matrix: use compute the point probability from Gaussian
			double logPDet;		   //log determinate of covariance matrix
		};

		// ------------------------------------------------------------------------------
		// Below defines JointGaussian related generic functions
		__host__ __device__
			inline std::vector<Real> getJointMean(const JointGaussian &dist)
		{
			ublas_vector mean = dist.getMean();
			std::vector<Real> m(mean.size());
			std::copy(mean.begin(), mean.end(), m.begin());
			return m;
		}

		///
		/// \brief Return PDF of x
		/// \param dist a distribution (Gaussian)
		/// \param x_ vector of a sample for probability estimation
		///
		__host__ __device__
			inline double getJointPdf(const JointGaussian &dist, const std::vector<Real> x_)
		{
			ublas_vector mean = dist.getMean();
			int k = mean.size();
			assert(x_.size() == k);
			ublas_vector x(k);  
			std::copy(x_.begin(), x_.end(), x.begin());
			x = x - mean;
			ublas_vector tmp;
			tmp = ublas::prod(ublas::trans(x), dist.getUMat());
			double density = ublas::inner_prod(tmp, tmp);

			return exp(-0.5 * (k * 1.83787706641 + dist.getLogDet() + density));
		}

		///
		/// \brief Return a random sample using random engine
		/// \param dist a distribution (Gaussian)
		/// \param rng random engine
		///
		__host__ __device__
			inline std::vector<Real> getJointSample(const JointGaussian &dist, thrust::default_random_engine &rng)
		{
				return getJointSample(dist, rng);
		}

		///
		/// \brief Print itself
		/// \param os outstream
		/// \param dist a distribution
		///
		__host__
			inline std::ostream& operator<<(std::ostream& os, const JointGaussian &dist)
		{
				os << "<JointGaussian: mean=" << dist.getMean() << ", covariance=" << dist.getCovariance() << ", eigenMatrix=" << dist.getEigenMat() << ">";
				return os;
		}

		///
		/// \brief Print the distribution type
		/// \param dist a distribution
		///
		__host__ __device__
		inline std::string getName(const JointGaussian &dist) {
			return "JointGaussian";
		}


	}  // namespace dist
}  // namespace edda

#endif  // DIST_GAUSSIAN_H_
