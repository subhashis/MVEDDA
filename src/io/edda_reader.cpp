#include <memory>
#include <vector>
#include <stdexcept>
#include <string>
#include <fstream>

#include "edda_reader.h"
#include "dataset/distr_array.h"
#include "dataset/dataset.h"
#include "core/interpolator.h"
#include "../test/bmp_image.h"

#include "edda.h"
#include "distributions/gaussian.h"
#include "distributions/distribution.h"
#include "distributions/variant.h"

#include "distributions/gmm.h"
#include "distributions/gaussian_mixture.h""
#include "distributions/joint_gaussian.h"
#include "distributions/joint_GMM.h"
#include "distributions/histogram.h"
#include "distributions/joint_histogram.h"

using namespace std;
using namespace edda;
using namespace dist;

#define LC_GENERIC_ARRAY

namespace edda {
#ifdef LC_GENERIC_ARRAY
	DistrArray * readMixArray(ifstream & myfile, int n)
	{
		dist::Variant* distArray;
		distArray = new dist::Variant[n];

		float* distData = new float[333 * 3]; // !!! this array is designed to avoid the time cost by multiple memory allocation and deletion. better check if the size is big enough before each time using it !!!

		for (int index = 0; index < n; index++){
			int distrTypeNumber = 1;
			myfile.read((char*)(&distrTypeNumber), sizeof (int));
			if (distrTypeNumber == 1){ //gaussian
				int nGM;
				myfile.read((char*)(&nGM), sizeof (int));

				myfile.read((char*)(distData), sizeof (float)*nGM * 3);

				dist::GMM new_gmm = dist::GMM(nGM);
				for (int m = 0; m<nGM; m++)
				{
					new_gmm.models[m].m = distData[3 * m];
					new_gmm.models[m].v = distData[3 * m + 1];
					new_gmm.models[m].w = distData[3 * m + 2];
				}
				distArray[index] = new_gmm;
			}
			else if (distrTypeNumber == 2){ //histogram
				int nbins;
				myfile.read((char*)(&nbins), sizeof (int));

				myfile.read((char*)(distData), sizeof (float));
				myfile.read((char*)(distData + 1), sizeof (float));
				myfile.read((char*)(distData + 2), sizeof (float)*nbins);

				distArray[index] = Histogram(distData, nbins);
			}
			else if (distrTypeNumber == 3){ //joint gmm
				int nVar, nComp;
				myfile.read((char*)(&nVar), sizeof (int));
				myfile.read((char*)(&nComp), sizeof (int));

				int sizeEachComp = 1 + nVar + (nVar + 1)*nVar / 2;
				myfile.read((char*)(distData), sizeof (float)* sizeEachComp * nComp);
				
				
				std::vector<Real> weights(nComp);
				std::vector<JointGaussian> gaus(nComp);
				for (int c = 0; c < nComp; c++)
				{
					weights[c] = distData[sizeEachComp * c];

					ublas_vector mean = ublas_vector(nVar, 0);
					ublas_matrix cov = ublas::zero_matrix<Real>(nVar, nVar);;
					for (int jj = 0; jj < nVar; jj++){
						mean(jj) = distData[sizeEachComp * c + 1 + jj];
					}
					int count = 0;
					for (int jj = 0; jj < nVar; jj++){
						for (int i = jj; i < nVar; i++){
							cov(jj, i) = distData[c*sizeEachComp + 1 + nVar + count];//cov row major
							count++;
						}
					}
					for (int jj = 0; jj < nVar; jj++){
						for (int i = 0; i < jj; i++){
							cov(jj, i) = cov(i, jj);
						}
					}
					gaus[c] = JointGaussian(mean, cov);
				}
				dist::JointGMM new_distr(weights, gaus, nVar, nComp);
				distArray[index] = new_distr;
			}
			else if (distrTypeNumber == 4) { //joint hist

				int* intData = new int[333 * 3]; // !!! this array is designed to avoid the time cost by multiple memory allocation and deletion. better check if the size is big enough before each time using it !!!

				int num_vars;
				myfile.read((char*)(&num_vars), sizeof(int));
				
				myfile.read((char*)(distData), sizeof(float)* num_vars);
				std::vector<Real> min_vals(distData, distData + num_vars);

				myfile.read((char*)(distData), sizeof(float)* num_vars);
				std::vector<Real> max_vals(distData, distData + num_vars);

				myfile.read((char*)(distData), sizeof(float)* num_vars);
				std::vector<Real> bin_widths(distData, distData + num_vars);

				myfile.read((char*)(intData), sizeof(int)* num_vars);
				std::vector<int> num_bins(intData, intData + num_vars);

				//read pdf
				boost::unordered_map<std::vector<int>, Real> pdf;
				int n_pdf;
				myfile.read((char*)(&n_pdf), sizeof(int));
				int count = 0;
				while (count < n_pdf) {
					int n_vint;
					myfile.read((char*)(&n_vint), sizeof(int));
					myfile.read((char*)(intData), sizeof(int)* n_vint);
					std::vector<int> vint(intData, intData + n_vint);
					float rfloat;
					myfile.read((char*)(&rfloat), sizeof(float));
					Real r = rfloat;
					pdf[vint] = r;					
					count++;
				}

				//read mean
				ublas_vector mean(num_vars, 0);
				myfile.read((char*)(distData), sizeof(float)* num_vars);
				for (int jj = 0; jj < num_vars; jj++) {
					mean(jj) = distData[jj];
				}

				//read cov	
				ublas_matrix cov = ublas::zero_matrix<Real>(num_vars, num_vars);;
				int sizeHalfCovMat = (num_vars + 1)*num_vars / 2;
				myfile.read((char*)(distData), sizeof(float)* sizeHalfCovMat);
				count = 0;
				for (int jj = 0; jj < num_vars; jj++) {
					for (int i = jj; i < num_vars; i++) {
						cov(jj, i) = distData[count];//cov row major
						count++;
					}
				}
				for (int jj = 0; jj < num_vars; jj++) {
					for (int i = 0; i < jj; i++) {
						cov(jj, i) = cov(i, jj);
					}
				}
				
				dist::JointHistogram new_distr(num_vars, min_vals, max_vals, bin_widths, num_bins, pdf, mean, cov);
				distArray[index] = new_distr;
			}
			else{
				myfile.close();
				throw NotImplementedException();
			}
		}

		delete[] distData;

		shared_ary<dist::Variant> pArray(new dist::Variant[n], n);
		for (int i = 0; i < n; i++)
		{
			pArray[i] = distArray[i];
		}
		DistrArray * abs_array = new ScalarDistrArray<dist::Variant>(pArray);
		return abs_array;		
	}

#elif
	DistrArray * readMixArray(ifstream & myfile, int n)
	{
		if (n == 0) { 
			dist::Variant* distArray = new dist::Variant[n];
			shared_ary<dist::Variant> pArray(distArray, n);
			DistrArray * abs_array = new ScalarDistrArray<dist::Variant>(pArray);
			return abs_array;
		}

		int distrTypeNumber = 1;
		myfile.read((char*)(&distrTypeNumber), sizeof(int));

		//seperate the process of univariate and multivariate. currently not support mixing of univariate and multivariate
		if (distrTypeNumber == 1 || distrTypeNumber == 2) {
			dist::Variant* distArray;
			distArray = new dist::Variant[n];

			float* distData = new float[333 * 3]; // !!! this array is designed to avoid the time cost by multiple memory allocation and deletion. better check if the size is big enough before each time using it !!!

			for (int index = 0; index < n; index++) {
				if (index > 0) {
					myfile.read((char*)(&distrTypeNumber), sizeof(int));
				}
				if (distrTypeNumber == 1) { //gaussian
					int nGM;
					myfile.read((char*)(&nGM), sizeof(int));

					myfile.read((char*)(distData), sizeof(float)*nGM * 3);

					dist::GMM new_gmm = dist::GMM(nGM);
					for (int m = 0; m<nGM; m++)
					{
						new_gmm.models[m].m = distData[3 * m];
						new_gmm.models[m].v = distData[3 * m + 1];
						new_gmm.models[m].w = distData[3 * m + 2];
					}
					distArray[index] = new_gmm;
				}
				else if (distrTypeNumber == 2) { //histogram
					int nbins;
					myfile.read((char*)(&nbins), sizeof(int));

					myfile.read((char*)(distData), sizeof(float));
					myfile.read((char*)(distData + 1), sizeof(float));
					myfile.read((char*)(distData + 2), sizeof(float)*nbins);

					distArray[index] = Histogram(distData, nbins);
				}
				else if (distrTypeNumber == 3 || distrTypeNumber == 4) { //joint distr
					//cannot mix univariate and multivariate
					myfile.close();
					throw NotSupportException();
				}
				else {
					myfile.close();
					throw NotImplementedException();
				}
			}

			delete[] distData;

			//shared_ary<dist::Variant> pArray(new dist::Variant[n], n);
			//for (int i = 0; i < n; i++)
			//{
			//	pArray[i] = distArray[i];
			//}
			shared_ary<dist::Variant> pArray(distArray, n);

			DistrArray * abs_array = new ScalarDistrArray<dist::Variant>(pArray);
			return abs_array;
		}
		else if (distrTypeNumber == 3) {
			dist::JointGMM* distArray;
			distArray = new dist::JointGMM[n];

			float* distData = new float[333 * 3]; // !!! this array is designed to avoid the time cost by multiple memory allocation and deletion. better check if the size is big enough before each time using it !!!

			for (int index = 0; index < n; index++) {
				if (index > 0) {
					myfile.read((char*)(&distrTypeNumber), sizeof(int));
				}
				if (distrTypeNumber == 1 || distrTypeNumber == 2 || distrTypeNumber == 4) { //other typesof distr
					//cannot mix other distr with a type of multivariate
					//in the future may support mix of JointGMM and JointHistogram
					myfile.close();
					throw NotSupportException();
				}
				else if (distrTypeNumber == 3) { //joint gmm
					int nVar, nComp;
					myfile.read((char*)(&nVar), sizeof(int));
					myfile.read((char*)(&nComp), sizeof(int));

					int sizeEachComp = 1 + nVar + (nVar + 1)*nVar / 2;
					myfile.read((char*)(distData), sizeof(float)* sizeEachComp * nComp);

					ublas_vector w = ublas_vector(nComp, 0);
					ublas_matrix m(nComp, nVar);
					ublas_matrix covs(nComp* nVar, nVar);
					for (int c = 0; c < nComp; c++)
					{
						w[c] = distData[sizeEachComp * c];
						ublas_vector mean = ublas_vector(nVar, 0);
						ublas_matrix cov = ublas::zero_matrix<Real>(nVar, nVar);;
						for (int jj = 0; jj < nVar; jj++) {
							mean(jj) = distData[sizeEachComp * c + 1 + jj];
						}
						int count = 0;
						for (int jj = 0; jj < nVar; jj++) {
							for (int i = jj; i < nVar; i++) {
								cov(jj, i) = distData[c*sizeEachComp + 1 + nVar + count];//cov row major
								count++;
							}
						}
						for (int jj = 0; jj < nVar; jj++) {
							for (int i = 0; i < jj; i++) {
								cov(jj, i) = cov(i, jj);
							}
						}
						//copy single GM to full GMM
						for (int jj = 0; jj < nVar; jj++) {
							m(c, jj) = mean(jj);
						}
						subrange(covs, c*nVar, (c + 1)*nVar, 0, nVar) = cov;
					}

					dist::JointGMM gmm;
					gmm.setGMM(nVar, nComp, w, m, covs);
					distArray[index] = gmm;
				}
				else {
					myfile.close();
					throw NotImplementedException();
				}
			}

			delete[] distData;

			////method1
			////have error in the future when using the DistrArray * stored in eddaDataset. indexing at 0 will be ok, but larger indexing will break the program
			//dist::Variant* test  = new dist::Variant[n];
			//for (int i = 0; i < n; i++)
			//{
			//	test[i] = distArray[i];
			//}
			//shared_ary<dist::Variant> pArray(test, n);
			//DistrArray * abs_array = new ScalarDistrArray<dist::Variant>(pArray);

			//method2
			//have error at the end of the program when recycling memory
			shared_ary<dist::JointGMM> pArray(distArray, n);
			DistrArray * abs_array = new JointDistrArray<dist::JointGMM>(pArray, 2);

			return abs_array;
		}
		
	}
#endif

	template <typename T>
	shared_ptr<Dataset<T> > loadEddaDatasetTemplate(const string &edda_file)
	{
		ifstream myfile(edda_file.c_str(), ios::binary);

		streampos begin = myfile.tellg();
		myfile.seekg(0, ios::end);
		streampos end = myfile.tellg();
		int fileByteSize = end - begin;

		myfile.seekg(0, ios::beg);

		char eddaFileMark[4];
		myfile.read(eddaFileMark, sizeof(char)*4);
		if (eddaFileMark[0] != 'E' || eddaFileMark[1] != 'D' || eddaFileMark[2] != 'D' || eddaFileMark[3] != 'A'){
			myfile.close();
			//TODO: throw a proper exception
			exit(0);
		}

		char majorVersion, minorVersion;
		myfile.read(&majorVersion, sizeof(char));
		myfile.read(&minorVersion, sizeof(char));

		if (majorVersion == 0 && minorVersion == 1){
			/*
			//not supported any more
			*/
		}
		else if (majorVersion == 0 && minorVersion == 2){
			int gridTypeNumber;
			myfile.read((char*)(&gridTypeNumber), sizeof(int));
			if (gridTypeNumber == 1){
				int dims[3];
				myfile.read((char*)(&dims), sizeof (int)* 3);
				float spacing[3];
				myfile.read((char*)(&spacing), sizeof (float)* 3);

				
				int numDistrArray = 1;
				myfile.read((char*)(&numDistrArray), sizeof (int));
				
				std::vector<DistrArray *> dVec(numDistrArray);
				for (int i = 0; i < numDistrArray; i++){
					dVec[i] = readMixArray(myfile, dims[0] * dims[1] * dims[2]);
				}
				
				myfile.close();

				return make_shared<Dataset<T>>(new RegularCartesianGrid(dims[0], dims[1], dims[2]), dVec);
			}
			else{
				myfile.close();
				throw NotImplementedException();
			}
		}
		else{
			myfile.close();
			throw NotImplementedException();
		}

		printf("Read file from %s.\n", edda_file.c_str());
	}



	shared_ptr<Dataset<Real> > loadEddaScalarDataset_noneVTK(const string &edda_file)
	{
		return loadEddaDatasetTemplate<Real>(edda_file);
	}
}