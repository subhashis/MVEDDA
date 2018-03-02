#include <string>
#include <vector>

#include "dataset/distr_array.h"
#include "distributions/distribution_modeler.h"
#include "distributions/joint_GMM.h"
#include "distributions/joint_gaussian.h"
#include "io/edda_reader.h"
#include "io/edda_writer.h"


#include "distributions/joint_histogram.h"


using namespace edda;
using namespace std;
using namespace edda::dist;

int main(int argc, char* argv[])
{

	if (argc < 10)
	{
		cout << "USAGE: check the parameter list!!\n";
		exit(11);
	}
	string filename1 = argv[1];
	string filename2 = argv[2];
	string filename3 = argv[3];
	int xdim = atoi(argv[4]);
	int ydim = atoi(argv[5]);
	int zdim = atoi(argv[6]);
	int blockXdim = atoi(argv[7]);
	int blockYdim = atoi(argv[8]);
	int blockZdim = atoi(argv[9]);


	//////for quicker testing
	//cout << "USAGE: default settings\n";
	//string filename1 = "E:/Data/edda/PRECIPf25.binLE.raw_corrected_2_subsampled";
	//string filename2 = "E:/Data/edda/QGRAUPf25.binLE.raw_corrected_2_subsampled";
	//string filename3 = "E:/Data/edda/QRAINf25.binLE.raw_corrected_2_subsampled";
	//int xdim =30;
	//int ydim = 30;
	//int zdim = 30;
	//int blockXdim = 10;
	//int blockYdim = 10;
	//int blockZdim = 10;


	int nVar = 3;
	int nComp = 2;


	//from here to the creation of "Dataset<Real> *ds" is the same with the example in MVDistrModelerExample_Hist.cpp
	float *inData1;
	inData1 = new float[xdim*ydim*zdim];

	FILE *fIn1 = fopen(filename1.c_str(), "rb");
	if (!fIn1)
	{
		fprintf(stderr, "Error opening file\n");
		exit(13);
	}
	size_t readStatus1 = fread(inData1, sizeof(float), xdim*ydim*zdim, fIn1);
	fclose(fIn1);

	float *inData2;
	inData2 = new float[xdim*ydim*zdim];

	FILE *fIn2 = fopen(filename2.c_str(), "rb");
	if (!fIn2)
	{
		fprintf(stderr, "Error opening file\n");
		exit(13);
	}
	size_t readStatus2 = fread(inData2, sizeof(float), xdim*ydim*zdim, fIn2);
	fclose(fIn2);

	float *inData3;
	inData3 = new float[xdim*ydim*zdim];

	FILE *fIn3 = fopen(filename3.c_str(), "rb");
	if (!fIn3)
	{
		fprintf(stderr, "Error opening file\n");
		exit(13);
	}
	size_t readStatus3 = fread(inData3, sizeof(float), xdim*ydim*zdim, fIn3);
	fclose(fIn3);

	int newW(0), newH(0), newD(0);

	if (xdim % blockXdim == 0)
		newW = xdim / blockXdim;
	else
	{
		fprintf(stderr, "Wrong dimension combination\n");
		exit(14);
	}
	if (ydim % blockYdim == 0)
		newH = ydim / blockYdim;
	else
	{
		fprintf(stderr, "Wrong dimension combination\n");
		exit(14);
	}
	if (zdim % blockZdim == 0)
		newD = zdim / blockZdim;
	else
	{
		fprintf(stderr, "Wrong dimension combination\n");
		exit(14);
	}

	DistributionModeler mv_dm(newW*newH*newD);

	cout << "starting partitioning\n";
	int counter = 0;

	for (size_t z = 0; z<zdim; z += blockZdim)
	{
		for (size_t y = 0; y<ydim; y += blockYdim)
		{
			for (size_t x = 0; x<xdim; x += blockXdim)
			{

				Real* data1 = (Real*)malloc(sizeof(Real)*blockXdim*blockYdim*blockZdim);
				Real* data2 = (Real*)malloc(sizeof(Real)*blockXdim*blockYdim*blockZdim);
				Real* data3 = (Real*)malloc(sizeof(Real)*blockXdim*blockYdim*blockZdim);
				std::vector<Real> min_val(3, FLT_MAX);
				std::vector<Real> max_val(3, -FLT_MAX);
				std::vector<int> bins(3, 30);


				int i = 0;
				for (size_t zz = z; zz<z + blockZdim; zz++)
				{
					for (size_t yy = y; yy<y + blockYdim; yy++)
					{
						for (size_t xx = x; xx<x + blockXdim; xx++)
						{
							data1[i] = (Real)(inData1[zz*xdim*ydim + yy*xdim + xx]);
							data2[i] = (Real)(inData2[zz*xdim*ydim + yy*xdim + xx]);
							data3[i] = (Real)(inData3[zz*xdim*ydim + yy*xdim + xx]);
							if (data1[i]<min_val[0]) min_val[0] = data1[i];
							if (data1[i]>max_val[0]) max_val[0] = data1[i];
							if (data2[i]<min_val[1]) min_val[1] = data2[i];
							if (data2[i]>max_val[1]) max_val[1] = data2[i];
							if (data3[i]<min_val[2]) min_val[2] = data3[i];
							if (data3[i]>max_val[2]) max_val[2] = data3[i];
							i++;
						}
					}
				}

				std::vector<Real*> trainSamples;
				trainSamples.push_back(data1);
				trainSamples.push_back(data2);
				trainSamples.push_back(data3);

				std::cout << "dimensions: [" << z << "][" << y << "][" << x << "]\n";

				mv_dm.computeJointHistogram(trainSamples, blockXdim*blockYdim*blockZdim, min_val, max_val, bins, counter);

				counter++;
			}
		}
	}
	std::vector<DistrArray *> dVec;
	dVec.push_back(mv_dm.getDistrArray());

	Dataset<Real> *ds = new Dataset<Real>(new RegularCartesianGrid(newW, newH, newD), dVec);
	shared_ptr<Dataset<Real>> shr_ds(ds);

	//write the dataset using the writer
	writeEddaDataset(shr_ds, "testDataJointHist.edda");

	//read the dataset using the reader
	shared_ptr<Dataset<Real>> shr_ds2 = loadEddaScalarDataset_noneVTK("testDataJointHist.edda");
	

	//basic compare
	int numDistrArray1 = shr_ds->getNumDistrArray();
	int numDistrArray2 = shr_ds2->getNumDistrArray();
	if (numDistrArray1 != numDistrArray2){
		cout << "Joint Histogram IO failed! number of arrays changed! " << endl;
		return 0;
	}

	DistrArray *array1 = shr_ds->getArray(0); //we know there is only 1 array
	DistrArray *array2 = shr_ds2->getArray(0); //we know there is only 1 array
	int n1 = array1->getLength();
	int n2 = array2->getLength();
	if (n1 != n2){
		cout << "Joint Histogram IO failed! length of arrays changed! " << endl;
		return 0;
	}

	//print and compare one arbitrary candidate
	cout << endl << "Compare one single Histogram from original dataset, and the dataset after IO: " << endl << endl;

	dist::Variant curDist1 = array1->getDistr(n1 / 2);
	dist::JointHistogram curJHist1 = boost::get<dist::JointHistogram>(curDist1);
	cout << "joint Histogram No. " << n1 / 2 << " from the original dataset:" << endl;
	cout << curJHist1 << endl << endl;

	dist::Variant curDist2 = array2->getDistr(n1 / 2);
	dist::JointHistogram curJHist2 = boost::get<dist::JointHistogram>(curDist2);
	cout << "joint Histogram No. " << n1 / 2 << " from the dataset after IO:" << endl;
	cout << curJHist2 << endl;
	
	//compare every parameter one by one
	double dif = 0.0;

	for (int j = 0; j < n1; j++){
		dist::Variant curDist1 = array1->getDistr(j);
		dist::Variant curDist2 = array2->getDistr(j);

		string s = getName(curDist2);
		if (s.compare(0, 15, "JointHistogram") != 0) {
			cout << "Joint Histogram IO failed! arrays element not joint Histogram! " << endl;
			return 0;
		}

		dist::JointHistogram curJHist1 = boost::get<dist::JointHistogram>(curDist1);
		dist::JointHistogram curJHist2 = boost::get<dist::JointHistogram>(curDist2);
		
		int num_vars1 = curJHist1.getNumVars();
		int num_vars2 = curJHist2.getNumVars();
		if (num_vars1 != num_vars2){
			cout << "Joint Histogram IO failed! num_vars of at least one joint Histogram changed! " << endl;
			return 0;
		}

		std::vector<Real> min_vals1 = curJHist1.getMinVals();
		std::vector<Real> min_vals2 = curJHist2.getMinVals();
		std::vector<Real> max_vals1 = curJHist1.getMaxVals();
		std::vector<Real> max_vals2 = curJHist2.getMaxVals();
		std::vector<Real> bin_widths1 = curJHist1.getBinWidths();
		std::vector<Real> bin_widths2 = curJHist2.getBinWidths();
		std::vector<int> num_bins1 = curJHist1.getNumBins();
		std::vector<int> num_bins2 = curJHist2.getNumBins();

		for (int v = 0; v < num_vars1; v++) {
			dif += abs(min_vals1[v] - min_vals2[v]);
			dif += abs(max_vals1[v] - max_vals2[v]);
			dif += abs(bin_widths1[v] - bin_widths2[v]);
			dif += abs(num_bins1[v] - num_bins2[v]);
		}

		//compare pdf
		boost::unordered_map<std::vector<int>, Real> pdf1 = curJHist1.getDistr();
		boost::unordered_map<std::vector<int>, Real> pdf2 = curJHist2.getDistr();
		int n_pdf1 = pdf1.size();
		int n_pdf2 = pdf2.size();
		if (n_pdf1 != n_pdf2) {
			cout << "Joint Histogram IO failed! pdf size of at least one joint Histogram changed! " << endl;
			return 0;
		}
		for (auto it1 = pdf1.begin(); it1 != pdf1.end(); ++it1) {
			std::vector<int> vint1 = it1->first;
			if (pdf2.find(vint1) == pdf2.end()) {
				cout << "Joint Histogram IO failed! at least one key of at least one pdf of at least one joint Histogram changed! " << endl;
				return 0;
			}
			Real r1 = it1->second;
			Real r2 = pdf2[vint1];
			int n_vint1 = vint1.size();
			dif += abs(r1 - r2);
		}
		
		//compare mean and cov
		std::vector<Real> curMean1 = getJointMean(curJHist1);
		ublas_matrix curCov1 = curJHist1.getCovariance();
		std::vector<Real> curMean2 = getJointMean(curJHist2);
		ublas_matrix curCov2 = curJHist2.getCovariance();
		for (int i = 0; i < num_vars1; i++) {
			dif += abs(curMean1[i] - curMean2[i]);
		}
		for (int j = 0; j < num_vars1; j++) {
			for (int i = j; i < num_vars1; i++) {
				dif += abs(curCov1(j, i) - curCov2(j, i));
			}
		}		
	}

	cout << "the total difference between parameters of the modeler result and the IO result is: " << dif << endl;

	return 1;
}