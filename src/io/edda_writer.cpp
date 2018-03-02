#include <vtkImageData.h>
#include <vtkXMLImageDataWriter.h>
#include <vtkSmartPointer.h>
#include <vtkPointData.h>
#include <vtkFloatArray.h>
#include <vtkStringArray.h>
#include <vtkFieldData.h>

#include "distributions/gaussian_mixture.h"
#include "distributions/gmm.h"
#include "distributions/histogram.h"
#include "distributions/joint_gaussian.h"
#include "distributions/joint_GMM.h"
#include "distributions/joint_histogram.h"

#include "edda_writer.h"

using namespace std;

namespace edda{

	const dist::GMMTuple getGmmModels(dist::Variant &distr, int GMs, int model)
	{
		switch (GMs)
		{
		case 2:
			return boost::get<dist::GaussianMixture<2> >(distr).models[model];
		case 3:
			return boost::get<dist::GaussianMixture<3> >(distr).models[model];
		case 4:
			return boost::get<dist::GaussianMixture<4> >(distr).models[model];
		case 5:
			return boost::get<dist::GaussianMixture<5> >(distr).models[model];
		default:
			throw runtime_error("The Gaussian mixture models exceeds default size");
		}
	}


	void writeMixArrays(ofstream & myFile, DistrArray *array)
	{
		float* gmData = new float[1000]; // !!! this array is designed to avoid the time cost by multiple memory allocation and deletion. better check if the size is big enough before each time using it !!!

		//write data
		int n = array->getLength();
		for (int j = 0; j < n; j++){
			dist::Variant curDist = array->getDistr(j);

			string s = getName(curDist);
			if (s.compare(0, 15, "GaussianMixture") == 0) {
				// Only compare the first 15 chars because this string ends with the number of Gaussian models

				//	8. distr type. 1: GaussianMixture. 2: Histogram
				//	if GaussianMixture:
				int distrTypeNumber = 1;
				myFile.write((char*)(&distrTypeNumber), sizeof (int));

				//	9.1. number of gaussian components
				int GMs = stoi(s.substr(15));
				myFile.write((char*)(&GMs), sizeof (int));

				//	9.2. GMM parameters
				std::vector<dist::GMMTuple> gmmtv = boost::get<dist::GMM >(curDist).models;
				for (int i = 0; i < GMs * 3; i++){
					if (i % 3 == 0)
						gmData[i] = gmmtv[i / 3].m;
					else if (i % 3 == 1)
						gmData[i] = gmmtv[i / 3].v;
					else
						gmData[i] = gmmtv[i / 3].w;
				}
				myFile.write((char*)(gmData), sizeof (float)*GMs * 3);
			}
			else if (s.compare(0, 15, "Histogram") == 0) {
				//	8. distr type. 1: GaussianMixture. 2: Histogram
				//	if Histogram:
				int distrTypeNumber = 2;
				myFile.write((char*)(&distrTypeNumber), sizeof (int));

				dist::Histogram curHist = boost::get<dist::Histogram>(curDist);

				int nbins = curHist.getBins();
				float minv = curHist.getMinValue();
				float maxv = curHist.getMaxValue();
				//	9.1. number of bins
				myFile.write((char*)&nbins, sizeof(int));
				//	9.2. min/max values
				myFile.write((char*)&minv, sizeof(float));
				myFile.write((char*)&maxv, sizeof(float));

				//	9.3. bin values
				float * values = (float*)malloc(sizeof(float)*nbins);
				for (int b = 0; b < nbins; b++)
				{
					values[b] = curHist.getBinValue(b);
				}
				myFile.write((char*)(values), sizeof(float)*nbins);
				free(values);
			}
			else if (s.compare(0, 15, "JointGMM") == 0) {
				int distrTypeNumber = 3;
				myFile.write((char*)(&distrTypeNumber), sizeof (int));

				dist::JointGMM curJGMM = boost::get<dist::JointGMM>(curDist);

				//	9.1. number of variales and number of gaussian components
				int nVar = curJGMM.getNumVariables();
				int nComp = curJGMM.getNumComponents();
				myFile.write((char*)(&nVar), sizeof(int));
				myFile.write((char*)(&nComp), sizeof(int));

				//	9.2. weights, mean, covariance matrix
				int sizeEachComp = 1 + nVar + (nVar + 1)*nVar / 2;
				float * values = (float*)malloc(sizeof(float)*sizeEachComp*nComp);
				for (int c = 0; c < nComp; c++)
				{
					values[c*sizeEachComp] = curJGMM.getWeight(c);

					dist::JointGaussian jg = curJGMM.getJointGaussian(c);
					const ublas_vector curMean = jg.getMean();
					const ublas_matrix curCov = jg.getCovariance();
					for (int i = 0; i < nVar; i++){
						values[c*sizeEachComp + 1 + i] = curMean(i);
					}
					//since the cov matrix must be symmetric, only save half
					int count = 0;
					for (int j = 0; j < nVar; j++){
						for (int i = j; i < nVar; i++){
							values[c*sizeEachComp + 1 + nVar + count] = curCov(j, i);//row major
							count++;
						}
					}
				}
				myFile.write((char*)(values), sizeof(float)*sizeEachComp*nComp);
				free(values);
			}
			else if (s.compare(0, 15, "JointHistogram") == 0) {
				int distrTypeNumber = 4;
				myFile.write((char*)(&distrTypeNumber), sizeof(int));

				dist::JointHistogram curJHist = boost::get<dist::JointHistogram>(curDist);

				//	number of variales
				int num_vars = curJHist.getNumVars();
				myFile.write((char*)(&num_vars), sizeof(int));

				std::vector<float> temp(num_vars);//convert to float for the current version

				std::vector<Real> min_vals = curJHist.getMinVals();
				for (int v = 0; v < num_vars; v++) {
					temp[v] = min_vals[v];
				}
				myFile.write((char*)(&(temp[0])), sizeof(float)*num_vars);

				std::vector<Real> max_vals = curJHist.getMaxVals();
				for (int v = 0; v < num_vars; v++) {
					temp[v] = max_vals[v];
				}
				myFile.write((char*)(&(temp[0])), sizeof(float)*num_vars); 
				
				std::vector<Real> bin_widths = curJHist.getBinWidths();
				for (int v = 0; v < num_vars; v++) {
					temp[v] = bin_widths[v];
				}
				myFile.write((char*)(&(temp[0])), sizeof(float)*num_vars);

				std::vector<int> num_bins = curJHist.getNumBins();
				myFile.write((char*)(&(num_bins[0])), sizeof(int)*num_vars);

				//next store the pdf
				boost::unordered_map<std::vector<int>, Real> pdf = curJHist.getDistr();
				//size of the map
				int n_pdf = pdf.size();
				myFile.write((char*)(&n_pdf), sizeof(int));
				//each pair of the map
				for (auto it = pdf.begin(); it != pdf.end(); ++it) {
					std::vector<int> vint = it->first;
					float r = (float)(it->second);
					int n_vint = vint.size();
					myFile.write((char*)(&n_vint), sizeof(int));
					myFile.write((char*)(&(vint[0])), sizeof(int)*n_vint);
					myFile.write((char*)(&r), sizeof(float));
				}

				//mean
				std::vector<Real> mean = getJointMean(curJHist);
				for (int v = 0; v < num_vars; v++) {
					temp[v] = mean[v];
				}
				myFile.write((char*)(&(temp[0])), sizeof(float)*num_vars);

				//cov matrix
				ublas_matrix cov = curJHist.getCovariance();
				int sizeHalfCovMat = (num_vars + 1)*num_vars / 2;
				float * values = (float*)malloc(sizeof(float)*sizeHalfCovMat);
				//since the cov matrix must be symmetric, only save half
				int count = 0;
				for (int j = 0; j < num_vars; j++) {
					for (int i = j; i < num_vars; i++) {
						values[count] = cov(j, i);//row major
						count++;
					}
				}
				myFile.write((char*)(values), sizeof(float)*sizeHalfCovMat);
				free(values);
			}
			else{
				throw NotImplementedException();
			}
		}
		delete[] gmData;
	}


	template <typename T>
	void writeEddaDatasetTemplate(shared_ptr<Dataset<T> > dataset, const string &edda_file)
	{
		const string array_name_prefix = "";

		ofstream myFile(edda_file.c_str(), ios::out | ios::binary);

		//1. "EDDA" as a marker
		//2. version number	

		char eddaFileMark[4] = { 'E', 'D', 'D', 'A' };
		myFile.write(eddaFileMark, sizeof (char)* 4);
		char majorVersion = 0, minorVersion = 2;
		myFile.write(&majorVersion, sizeof (char));
		myFile.write(&minorVersion, sizeof (char));

		CartesianGrid *cartesianGrid = dynamic_cast<CartesianGrid *>(dataset->getGrid());
		if (cartesianGrid) {
			//3. gridType. 1: Regular CartesianGrid.
			//if Regular CartesianGrid:
			//	4. dimension
			//	5. spacing

			int gridTypeNumber = 1;
			myFile.write((char*)(&gridTypeNumber), sizeof (int));
			int *dims = dataset->getDimension();
			myFile.write((char*)(dims), sizeof (int)* 3);
			float spacing[3];
			dataset->getSpacing(spacing[0], spacing[1], spacing[2]);
			myFile.write((char*)(spacing), sizeof (float)* 3);
		}
		else {
			// TODO for other grid types
			throw NotImplementedException();
		}

		//	6. number of distribution arrays
		int numDistrArray = dataset->getNumDistrArray();
		myFile.write((char*)(&numDistrArray), sizeof (int));

		for (int ida = 0; ida < numDistrArray; ida++){
			DistrArray *array = dataset->getArray(ida);
			writeMixArrays(myFile, array);
		}

		printf("Saved file to %s.\n", edda_file.c_str());
		myFile.close();
	}

	void writeEddaDataset(shared_ptr<Dataset<Real> > dataset, const string &edda_file)
	{
		writeEddaDatasetTemplate(dataset, edda_file);
	}
}