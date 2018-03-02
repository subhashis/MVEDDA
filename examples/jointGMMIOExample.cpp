#include <string>
#include <vector>

#include "dataset/distr_array.h"
#include "distributions/distribution_modeler.h"
#include "distributions/joint_GMM.h"
#include "distributions/joint_gaussian.h"
#include "io/edda_reader.h"
#include "io/edda_writer.h"

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


	//from here to the creation of "Dataset<Real> *ds" is the same with the example in MVDistrModelerExample_GMM.cpp
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


	//edda ensemble data modeling
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
							i++;
						}
					}
				}

				std::vector<Real*> trainSamples;
				trainSamples.push_back(data1);
				trainSamples.push_back(data2);
				trainSamples.push_back(data3);

				std::cout << "dimensions: [" << z << "][" << y << "][" << x << "]\n";

				mv_dm.computeJointGMM(trainSamples, blockXdim*blockYdim*blockZdim, nComp, counter);

				counter++;
			}
		}
	}
	std::vector<DistrArray *> dVec;
	dVec.push_back(mv_dm.getDistrArray());

	Dataset<Real> *ds = new Dataset<Real>(new RegularCartesianGrid(newW, newH, newD), dVec);
	shared_ptr<Dataset<Real>> shr_ds(ds);

	//write the dataset using the writer
	writeEddaDataset(shr_ds, "testDataJointGMM.edda");

	//read the dataset using the reader
	shared_ptr<Dataset<Real>> shr_ds2 = loadEddaScalarDataset_noneVTK("testDataJointGMM.edda");


	//basic compare
	int numDistrArray1 = shr_ds->getNumDistrArray();
	int numDistrArray2 = shr_ds2->getNumDistrArray();
	if (numDistrArray1 != numDistrArray2){
		cout << "Joint GMM IO failed! number of arrays changed! " << endl;
		return 0;
	}

	DistrArray *array1 = shr_ds->getArray(0); //we know there is only 1 array
	DistrArray *array2 = shr_ds2->getArray(0); //we know there is only 1 array
	int n1 = array1->getLength();
	int n2 = array2->getLength();
	if (n1 != n2){
		cout << "Joint GMM IO failed! length of arrays changed! " << endl;
		return 0;
	}

	//print and compare one arbitrary candidate
	cout << endl << "Compare one single GMM from original dataset, and the dataset after IO: " << endl << endl;
	dist::Variant curDist1 = array1->getDistr(n1 / 2);
	dist::JointGMM curJGMM1 = boost::get<dist::JointGMM>(curDist1);
	cout << "joint GMM No. " << n1 / 2 << " of the original dataset:" << endl;
	cout << curJGMM1 << endl << endl;

	dist::Variant curDist2 = array2->getDistr(n1 / 2);
	dist::JointGMM curJGMM2 = boost::get<dist::JointGMM>(curDist2);
	cout << "joint GMM No. " << n1 / 2 << " of the dataset after IO:" << endl;
	cout << curJGMM2 << endl;

	////compare sample. currently not supported
	//vector<Real> sample1 = curJGMM1.getJointSample();
	//vector<Real> sample2 = curJGMM2.getJointSample();
	//cout << "A sample from joint GMM No. " << n1 / 2 << " of the original dataset:" << endl;
	//for (int i = 0; i < sample1.size(); i++) {
	//	cout << sample1[i] << " ";
	//}
	//cout << endl;
	//cout << "A sample from joint GMM No. " << n1 / 2 << " of the dataset after IO:" << endl;
	//for (int i = 0; i < sample2.size(); i++) {
	//	cout << sample2[i] << " ";
	//}
	//cout << endl;


	//compare every parameter one by one
	double dif = 0.0;

	for (int j = 0; j < n1; j++){
		dist::Variant curDist1 = array1->getDistr(j);
		dist::Variant curDist2 = array2->getDistr(j);

		string s = getName(curDist2);
		if (s.compare(0, 15, "JointGMM") != 0) {
			cout << "Joint GMM IO failed! arrays element not joint GMM! " << endl;
			return 0;
		}

		dist::JointGMM curJGMM1 = boost::get<dist::JointGMM>(curDist1);
		dist::JointGMM curJGMM2 = boost::get<dist::JointGMM>(curDist2);
		
		int nVar1 = curJGMM1.getNumVariables();
		int nComp1 = curJGMM1.getNumComponents();
		int nVar2 = curJGMM2.getNumVariables();
		int nComp2 = curJGMM2.getNumComponents();
		if (nVar1 != nVar2 || nComp1 != nComp2){
			cout << "Joint GMM IO failed! nVar or nComp of at least one joint GMM changed! " << endl;
			return 0;
		}

		for (int c = 0; c < nComp1; c++)
		{
			dif += abs(curJGMM1.getWeight(c) - curJGMM2.getWeight(c));

			dist::JointGaussian jg1 = curJGMM1.getJointGaussian(c);
			const ublas_vector curMean1 = jg1.getMean();
			const ublas_matrix curCov1 = jg1.getCovariance();
			dist::JointGaussian jg2 = curJGMM2.getJointGaussian(c);
			const ublas_vector curMean2 = jg2.getMean();
			const ublas_matrix curCov2 = jg2.getCovariance();
			for (int i = 0; i < nVar; i++){
				dif += abs(curMean1(i) - curMean2(i));
			}
			for (int j = 0; j < nVar; j++){
				for (int i = j; i < nVar; i++){
					dif += abs(curCov1(j, i) - curCov2(j, i));
				}
			}
		}
	}

	cout << "the total difference between parameters of the modeler result and the IO result is: " << dif << endl;

	return 1;
}