#include <iostream>
#include <stdlib.h>
#include <cmath>
#include <vector>
#include <sstream>

#include "distributions/distribution_modeler.h"
#include "dataset/grid.h"

#define IOTEST

#ifdef IOTEST
#include "io/edda_writer.h"
#include "io/edda_reader.h"
#endif

using namespace std;
using namespace edda;



int main(int argc, char* argv[])
{
	string filename;
	int xdim;
	int ydim;
	int zdim;
	int blockXdim;
	int blockYdim;
	int blockZdim;

		cout << "USAGE: default settings\n";
		filename = string(SAMPLE_DATA_PATH) + "/subIsabel.raw";
		xdim = 100;
		ydim = 100;
		zdim = 100;
		blockXdim = 10;
		blockYdim = 10;
		blockZdim = 10;


	float *inData;
	inData = new float[xdim*ydim*zdim];	
    
    FILE *fIn = fopen(filename.c_str(),"rb");
    if(!fIn)
    {
		fprintf(stderr, "Error opening file %s\n", filename.c_str());
        exit(13);
    }
    size_t readStatus = fread(inData, sizeof(float), xdim*ydim*zdim, fIn);
    fclose(fIn);

    int newW(0),newH(0),newD(0);
	newW = xdim / blockXdim;
	newH = ydim / blockYdim;
	newD = zdim / blockZdim;

    //first creata an edda data using DistributionModeler
	DistributionModeler dm1(newW*newH*newD);
	DistributionModeler dm2(newW*newH*newD);
	DistributionModeler dm3(newW*newH*newD);
	float *data1 = new float[blockXdim*blockYdim*blockZdim];
	float *data2 = new float[blockXdim*blockYdim*blockZdim];
	float *data3 = new float[blockXdim*blockYdim*blockZdim];
	int counter =0;
    for(size_t z=0; z<zdim; z += blockZdim)
    {
    	for(size_t y=0; y<ydim; y += blockYdim)
    	{
    		for(size_t x=0; x<xdim; x += blockXdim)
    		{	
    			int i = 0;
    			for(size_t zz=z; zz<z+blockZdim; zz++)
    			{
    				for(size_t yy=y; yy<y+blockYdim; yy++)
    				{
    					for(size_t xx=x; xx<x+blockXdim; xx++)
    					{
							//duplicate the same number 3 times with a little change, to fake a vector3 data
    						data1[i] = inData[zz*xdim*ydim + yy*xdim + xx];
							data2[i] = inData[zz*xdim*ydim + yy*xdim + xx]+1;
							data3[i] = inData[zz*xdim*ydim + yy*xdim + xx]+2;
							i++;
    					}
    				}
    			}

				int re = counter % 2;
				if (re == 0){
					dm1.computeGMM(data1, blockXdim*blockYdim*blockZdim, 2, counter);
					dm2.computeHistogram(data2, blockXdim*blockYdim*blockZdim, counter, 4);
					dm3.computeGMM(data3, blockXdim*blockYdim*blockZdim, 4, counter);
				}
				else if (re == 1){
					dm1.computeHistogram(data1, blockXdim*blockYdim*blockZdim, counter, 3);
					dm2.computeGMM(data2, blockXdim*blockYdim*blockZdim, 3, counter);
					dm3.computeHistogram(data3, blockXdim*blockYdim*blockZdim, counter, 5);
				}
				counter++;
    		}
    	}
    }
	delete[] inData;
	delete[] data1;
	delete[] data2;
	delete[] data3;

	std::vector<DistrArray *> dVec;
	dVec.push_back(dm1.getDistrArray());
	dVec.push_back(dm2.getDistrArray());
	dVec.push_back(dm3.getDistrArray());    
	Dataset<Real> *ds = new Dataset<Real>(new RegularCartesianGrid(newW, newH, newD), dVec);
	shared_ptr<Dataset<Real>> shr_ds(ds);

	//write the dataset using the writer
	writeEddaDataset(shr_ds, "testData.edda");

	//read the dataset using the reader
	shared_ptr<Dataset<Real>> shr_ds2 = loadEddaScalarDataset_noneVTK("testData.edda");
	

	//single test of the result of IO functions
	std::vector<dist::Variant> distrs;

	cout << endl;
	distrs = shr_ds->at_comp_distr_new(6, 5, 5);
	cout << "modeler result at_comp(6,5,5) : " << endl;
	for (int i = 0; i < distrs.size(); i++){
		cout << distrs[i] << endl;
	}
	cout << endl;
	distrs = shr_ds2->at_comp_distr_new(6, 5, 5);
	cout << "IO reuslt at_comp(6,5,5) : " << endl;
	for (int i = 0; i < distrs.size(); i++){
		cout << distrs[i] << endl;
	}
	cout << endl;

	distrs = shr_ds->at_comp_distr_new(7, 5, 5);
	cout << "modeler result at_comp(7,5,5) : " << endl;
	for (int i = 0; i < distrs.size(); i++){
		cout << distrs[i] << endl;
	}
	cout << endl;
	distrs = shr_ds2->at_comp_distr_new(7, 5, 5);
	cout << "IO reuslt at_comp(7,5,5) : " << endl;
	for (int i = 0; i < distrs.size(); i++){
		cout << distrs[i] << endl;
	}
	cout << endl;


	//one-by-one test of the result of IO functions
	int* dims = shr_ds2->getDimension();
	double dif = 0;
	for (int k = 0; k < dims[2]; k++){
		for (int j = 0; j < dims[1]; j++){
			for (int i = 0; i < dims[0]; i++){
				std::vector<dist::Variant> distrs = shr_ds->at_comp_distr_new(i, j, k);
				std::vector<dist::Variant> distrs2 = shr_ds2->at_comp_distr_new(i, j, k);
				if (distrs.size() == distrs2.size()){
					for (int i = 0; i < distrs.size(); i++){
						if (getName(distrs[i]).compare(0, 15, "GaussianMixture") == 0 && getName(distrs2[i]).compare(0, 15, "GaussianMixture") == 0){
							dist::GMM curDist1 = boost::get<dist::GMM >(distrs[i]);
							dist::GMM curDist2 = boost::get<dist::GMM >(distrs2[i]);
							for (int model = 0; model < curDist2.getNumComponenets(); model++){
								dif = dif + abs(curDist1.models[model].m - curDist2.models[model].m)
									+ abs(curDist1.models[model].v - curDist2.models[model].v)
									+ abs(curDist1.models[model].w - curDist2.models[model].w);
							}
						}
						else if (getName(distrs[i]).compare(0, 15, "Histogram") == 0 && getName(distrs2[i]).compare(0, 15, "Histogram") == 0){
							dist::Histogram curDist1 = boost::get<dist::Histogram >(distrs[i]);
							dist::Histogram curDist2 = boost::get<dist::Histogram >(distrs2[i]);
							int nbins1 = curDist1.getBins();
							int nbins2 = curDist2.getBins();
							float minv1 = curDist1.getMinValue();
							float maxv1 = curDist1.getMaxValue();
							float minv2 = curDist2.getMinValue();
							float maxv2 = curDist2.getMaxValue();
							dif = dif + abs(nbins1 - nbins2) + abs(minv1 - minv2) + abs(maxv1 - maxv2);
							for (int b = 0; b < min(nbins1, nbins2); b++){
								dif = dif + abs(curDist1.getBinValue(b) - curDist2.getBinValue(b));
							}
						}
						else{
							dif = dif + 1000;
						}
					}
				}
				else{
					dif = dif + 1000;
				}
			}
		}
	}
	cout << "the total difference between the modeler result and the IO result is: " << dif << endl;

    
  	return 0;
}
