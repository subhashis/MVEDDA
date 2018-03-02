#include <iostream>
#include <fstream>
#include <stdio.h>
#include <stdlib.h>
#include <cmath>
#include <random>
#include <algorithm>
#include <vector>
#include <tuple>

#include "distributions/distribution.h"
#include "distributions/distribution_modeler.h"

#include "dataset/distr_array.h"
#include "edda_export.h"

#include "edda.h"
#include "distributions/gaussian.h"
#include "distributions/distribution.h"
#include "distributions/variant.h"

#include "distributions/gmm.h"
#include "distributions/gaussian_mixture.h"
#include "distributions/joint_gaussian.h"
#include "distributions/joint_GMM.h"
#include "distributions/histogram.h"
#include "distributions/joint_histogram.h"


using namespace std;
using namespace edda;
using namespace dist;



void distrArrayWriter(const string &fname, DistrArray *array)
{

	std::ofstream outFile;
	outFile.open(fname.c_str(), ios::out|ios::binary);
	int n = array->getLength();
	int numComp = array->getNumComponents();
	cout << "array length = " << n << endl;
	cout << "array num components = " << numComp << endl;
	for(int i=0; i<n; i++)
	{
		dist::Variant curDist = array->getDistr(i);
		//string s = getName(curDist);
		//cout << s << endl;
		//TODO:can select the distrbution ID when using mixed distribution types
		dist::Histogram curHist = boost::get<dist::Histogram>(curDist);
		int nBins = curHist.getBins();
		float minV = curHist.getMinValue();
		float maxV = curHist.getMaxValue();

		outFile.write((char*)&nBins, sizeof(int));
		outFile.write((char*)&minV, sizeof(float));
		outFile.write((char*)&maxV, sizeof(float));
		

		float *values = (float*)malloc(sizeof(float)*nBins);
		for(int b=0; b<nBins; b++)
		{
			values[b] = curHist.getBinValue(b);
		}
		outFile.write((char*)values, sizeof(float)*nBins);
		free(values);
		
	}
	cout << "Saved file to " << fname << endl;
	outFile.close();

	

}

void pccParameterWriter(const string &fname, float *array, int size)
{

	std::ofstream outFile;
	outFile.open(fname.c_str(), ios::out|ios::binary);

	outFile.write((char*)array, sizeof(float)*size);
	cout << "Saved file to " << fname << endl;
	outFile.close();
}

void dataWriter(const string &fname, float *array, int size)
{

	std::ofstream outFile;
	outFile.open(fname.c_str(), ios::out|ios::binary);

	outFile.write((char*)array, sizeof(float)*size);
	cout << "Saved file to " << fname << endl;
	outFile.close();

}

void vecDataWriter(const string &fname, std::vector<float>* dataVec)
{
	std::ofstream outFile;
	outFile.open(fname.c_str(), ios::out|ios::binary);
	int n = dataVec->size();
	float *values = (float*)malloc(sizeof(float)*n);
	for(int i=0; i<n; i++)
	{
		values[i] = (*dataVec)[i];
	}
	outFile.write((char*)values, sizeof(float)*n);
	free(values);
	cout << "Vector data saved to " << fname << endl;
	outFile.close();
}


float * pccParameterReader(const string &fname, int size)
{
	FILE *fIn;
	float *retArray = new float[size];
    fIn = fopen(fname.c_str(),"rb");
    if(!fIn)
    {
        fprintf(stderr, "Error opening file\n");
        exit(13);
    }
    fread(retArray, sizeof(float), size, fIn);
    fclose(fIn);

    return retArray;

}
float * pccParameterReader_multifile(string * fnames, int fn_size, int size)
{
	
	float *retArray = new float[size];
	int per_filesize = size/fn_size;

	for(int f=0; f<fn_size; f++){
		std::ifstream inFile;
		inFile.open(fnames[f].c_str(), ios::binary);
		
		streampos begin = inFile.tellg();
		inFile.seekg(0, ios::end);
		streampos end = inFile.tellg();
		int fileByteSize = end - begin;
		inFile.seekg(0, ios::beg);
		//cout << "fname[" << i << "] = " << fileByteSize/sizeof(float) << endl;
		inFile.read((char*)(retArray + f*per_filesize), sizeof (float)*per_filesize);

		inFile.close();
	}
    

    return retArray;

}

DistrArray * distrArrayReader(const string &fname, int size)
{
	//int distrType = 9 // currently hardcode it to select a distribution type

	std::ifstream inFile;
	inFile.open(fname.c_str(), ios::binary);
	
	streampos begin = inFile.tellg();
	inFile.seekg(0, ios::end);
	streampos end = inFile.tellg();
	int fileByteSize = end - begin;
	inFile.seekg(0, ios::beg);

	dist::Variant* dArray;
	dArray = new dist::Variant[size];

	for(int i=0; i<size; i++)
	{
		int nBins;
		inFile.read((char*)(&nBins), sizeof (int));

		float *distrData = (float*)malloc(sizeof(float)*(nBins+2));

		inFile.read((char*)(distrData), sizeof (float));
		inFile.read((char*)(distrData + 1), sizeof (float));
		inFile.read((char*)(distrData + 2), sizeof (float)*nBins);

		dArray[i] = Histogram(distrData, nBins);

		free(distrData);
		
	}
	shared_ary<dist::Variant> pArray(new dist::Variant[size], size);
	for (int i = 0; i < size; i++)
	{
		pArray[i] = dArray[i];
	}
	DistrArray * abs_array = new ScalarDistrArray<dist::Variant>(pArray);

	cout << "Loaded file " << fname  << " successfully."<< endl;

	return abs_array;

}

DistrArray * distrArrayReader_new(const string &fname, int size, int distrType)
{
	//int distrType = 9 // currently hardcode it to select a distribution type

	std::ifstream inFile;
	inFile.open(fname.c_str(), ios::binary);
	
	streampos begin = inFile.tellg();
	inFile.seekg(0, ios::end);
	streampos end = inFile.tellg();
	int fileByteSize = end - begin;
	inFile.seekg(0, ios::beg);

	dist::Variant* dArray;
	dArray = new dist::Variant[size];

	for(int i=0; i<size; i++)
	{
		if(distrType > 30)
		{
			int nBins;
			//inFile.read((char*)(&nBins), sizeof (int));
			nBins = distrType;

			float *distrData = (float*)malloc(sizeof(float)*(nBins+2));

			inFile.read((char*)(distrData), sizeof (float));
			inFile.read((char*)(distrData + 1), sizeof (float));
			inFile.read((char*)(distrData + 2), sizeof (float)*nBins);

			dArray[i] = Histogram(distrData, nBins);

			free(distrData);
		}
		else
		{
			int nGM = distrType/3;
			float *distrData = (float*)malloc(sizeof(float)*nGM*3);
			inFile.read((char*)(distrData), sizeof(float)*nGM*3);
			dist::GMM new_gmm = dist::GMM(nGM);
			for (int m = 0; m<nGM; m++)
			{
				new_gmm.models[m].m = distrData[3 * m];
				new_gmm.models[m].v = distrData[3 * m + 1];
				new_gmm.models[m].w = distrData[3 * m + 2];
			}
			dArray[i] = new_gmm;
		}

	}
	shared_ary<dist::Variant> pArray(new dist::Variant[size], size);
	for (int i = 0; i < size; i++)
	{
		pArray[i] = dArray[i];
	}
	DistrArray * abs_array = new ScalarDistrArray<dist::Variant>(pArray);

	cout << "Loaded file " << fname  << " successfully."<< endl;

	return abs_array;

}

DistrArray * distrArrayReader_new_multifile(string * fnames, int fn_size, int size, int num_distr_component, int start, int end, int distrType)
{
	dist::Variant* dArray;
	dArray = new dist::Variant[size];

	float * full_distr_data;
	int full_size = num_distr_component*size;
	full_distr_data = new float[full_size];

	int per_filesize = full_size/fn_size;

	for(int i=0; i<fn_size; i++)
	{
		std::ifstream inFile;
		inFile.open(fnames[i].c_str(), ios::binary);
		
		streampos begin = inFile.tellg();
		inFile.seekg(0, ios::end);
		streampos end = inFile.tellg();
		int fileByteSize = end - begin;
		inFile.seekg(0, ios::beg);
		//cout << "fname[" << i << "] = " << fileByteSize/sizeof(float) << endl;
		inFile.read((char*)(full_distr_data + i*per_filesize), sizeof (float)*per_filesize);

		inFile.close();


	}

	int idx=0;
	for(int f=0; f<fn_size; f++){
		int realStart = start + (f*per_filesize);
		int realEnd = end + (f*per_filesize);

		cout << "f[" << f << "]=(" << realStart << "," << realEnd << ")" <<endl;
		int testC = 0; 
 		for(int i=realStart; i<realEnd; i+=distrType){
			if(distrType > 30){
				int nBins;
				nBins = distrType-2;
				dist::Histogram hg = Histogram(&full_distr_data[i], nBins);
				dArray[idx] = hg;
				testC++;
				//cout << hg << endl;
				idx++;
			}
			else{
				int nGM = distrType/3;
				dist::GMM new_gmm = dist::GMM(nGM);
				for (int m = 0; m<nGM; m++)
				{
					new_gmm.models[m].m = full_distr_data[i + (3 * m)];
					new_gmm.models[m].v = full_distr_data[i + (3 * m + 1)];
					new_gmm.models[m].w = full_distr_data[i + (3 * m + 2)];
					//cout << "(" << full_distr_data[i + (3 * m)] << "," << full_distr_data[i + (3 * m + 1)] << "," << full_distr_data[i + (3 * m + 2)] << ")";
				}
				//cout << "test = " << testC << endl;
				testC++;
				//cout << "idx = " << idx << endl;
				dArray[idx] = new_gmm;
				//cout << new_gmm << endl;
				//cout << "_test1\n";
				idx++;
			}
		}
		cout << "idx = " << idx << endl;
		cout << "test = " << testC << endl;
	}
	delete[] full_distr_data;

	shared_ary<dist::Variant> pArray(new dist::Variant[size], size);
	for (int i = 0; i < size; i++)
	{
		pArray[i] = dArray[i];
	}
	DistrArray * abs_array = new ScalarDistrArray<dist::Variant>(pArray);

	for(int f=0; f<fn_size; f++){
		cout << "Loaded file " << fnames[f]  << " successfully."<< endl;
	}
	delete[] dArray;

	return abs_array;

}



void vectorWriter_int(const string &fname, std::vector<std::vector<int>>* vec_vec)
{
	std::ofstream outFile;
	outFile.open(fname.c_str(), ios::out|ios::binary);
	int n = vec_vec->size();
	//int numComp = array->getNumComponents();
	cout << "array length = " << n << endl;
	//cout << "array num components = " << numComp << endl;
	for(int i=0; i<n; i++)
	{
		int curVecSize = (*vec_vec)[i].size();
		outFile.write((char*)&curVecSize, sizeof(int));
		int *values = (int*)malloc(sizeof(int)*curVecSize);
		for(int b=0; b<curVecSize; b++)
		{
			values[b] = (*vec_vec)[i][b];
		}
		outFile.write((char*)values, sizeof(int)*curVecSize);
		free(values);		
	}

	cout << "Saved file to " << fname << endl;
	outFile.close();
}

float * dataReader(const string &fname, int size)
{
	FILE *fIn;
	float *retArray = new float[size];
    fIn = fopen(fname.c_str(),"rb");
    if(!fIn)
    {
        fprintf(stderr, "Error opening file\n");
        exit(13);
    }
    fread(retArray, sizeof(float), size, fIn);
    fclose(fIn);

    return retArray;

}

void vecDataReader(const string &fname, std::vector<std::vector<int> >* vec_vec)
{
	std::ifstream inFile;
	inFile.open(fname.c_str(), ios::binary);
	
	streampos begin = inFile.tellg();
	inFile.seekg(0, ios::end);
	streampos end = inFile.tellg();
	int fileByteSize = end - begin;
	inFile.seekg(0, ios::beg);

	int size = vec_vec->size();
	for(int i=0; i<size; i++){
		int curSize;
		inFile.read((char*)(&curSize), sizeof (int));
		int *ids = (int*)malloc(sizeof(int)*(curSize));
		inFile.read((char*)(ids), sizeof(int)*curSize);
		for(int j=0; j<curSize; j++){
			(*vec_vec)[i].push_back(ids[j]);
		}
		free(ids);

	}

}

