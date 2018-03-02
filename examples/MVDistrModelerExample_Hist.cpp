#include <iostream>
#include <stdlib.h>
#include <cmath>
#include <vector>
#include <sstream>

#include "distributions/distribution_modeler.h"

using namespace std;
using namespace edda;



int main(int argc, char* argv[])
{
	if(argc < 8)
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

    float *inData1;
    inData1 = new float[xdim*ydim*zdim];    
    
    FILE *fIn1 = fopen(filename1.c_str(),"rb");
    if(!fIn1)
    {
        fprintf(stderr, "Error opening file\n");
        exit(13);
    }
    size_t readStatus1 = fread(inData1, sizeof(float), xdim*ydim*zdim, fIn1);
    fclose(fIn1);

    float *inData2;
    inData2 = new float[xdim*ydim*zdim];    
    
    FILE *fIn2 = fopen(filename2.c_str(),"rb");
    if(!fIn2)
    {
        fprintf(stderr, "Error opening file\n");
        exit(13);
    }
    size_t readStatus2 = fread(inData2, sizeof(float), xdim*ydim*zdim, fIn2);
    fclose(fIn2);

    float *inData3;
    inData3 = new float[xdim*ydim*zdim];    
    
    FILE *fIn3 = fopen(filename3.c_str(),"rb");
    if(!fIn3)
    {
        fprintf(stderr, "Error opening file\n");
        exit(13);
    }
    size_t readStatus3 = fread(inData3, sizeof(float), xdim*ydim*zdim, fIn3);
    fclose(fIn3);

    int newW(0),newH(0),newD(0);

    if(xdim % blockXdim == 0)
        newW = xdim / blockXdim;
    else
    {
        fprintf(stderr, "Wrong dimension combination\n");
        exit(14);
    }
    if(ydim % blockYdim == 0)
        newH = ydim / blockYdim;
    else
    {
        fprintf(stderr, "Wrong dimension combination\n");
        exit(14);
    }
    if(zdim % blockZdim == 0)
        newD = zdim / blockZdim;
    else
    {
        fprintf(stderr, "Wrong dimension combination\n");
        exit(14);
    }


    //edda ensemble data modeling
    DistributionModeler mv_dm(newW*newH*newD);      

    cout << "starting partitioning\n";
    int counter =0;

    for(size_t z=0; z<zdim; z += blockZdim)
    {
        for(size_t y=0; y<ydim; y += blockYdim)
        {
            for(size_t x=0; x<xdim; x += blockXdim)
            {

                Real* data1 = (Real*)malloc(sizeof(Real)*blockXdim*blockYdim*blockZdim);
		        Real* data2 = (Real*)malloc(sizeof(Real)*blockXdim*blockYdim*blockZdim);
		        Real* data3 = (Real*)malloc(sizeof(Real)*blockXdim*blockYdim*blockZdim);
		        std::vector<Real> min_val(3, FLT_MAX);
		        std::vector<Real> max_val(3, -FLT_MAX);
		        std::vector<int> bins(3, 30);		
		

                int i = 0;
                for(size_t zz=z; zz<z+blockZdim; zz++)
                {
                    for(size_t yy=y; yy<y+blockYdim; yy++)
                    {
                        for(size_t xx=x; xx<x+blockXdim; xx++)
                        {
                            data1[i] = (Real)(inData1[zz*xdim*ydim + yy*xdim + xx]);
                            data2[i] = (Real)(inData2[zz*xdim*ydim + yy*xdim + xx]);
                            data3[i] = (Real)(inData3[zz*xdim*ydim + yy*xdim + xx]);
			                if(data1[i]<min_val[0]) min_val[0] = data1[i];
			                if(data1[i]>max_val[0]) max_val[0] = data1[i];
			                if(data2[i]<min_val[1]) min_val[1] = data2[i];
			                if(data2[i]>max_val[1]) max_val[1] = data2[i];
			                if(data3[i]<min_val[2]) min_val[2] = data3[i];
			                if(data3[i]>max_val[2]) max_val[2] = data3[i];
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
	dVec.push_back(mv_dm.getMVDistrArray(3));
	
	Dataset<Real> *ds = new Dataset<Real> (new RegularCartesianGrid(newW, newH, newD), dVec);

    //Testing: to see if the edda dataset was properly created with the joint distribution.
    /*shared_ptr<Dataset<Real>> shr_ds(ds);
    std::cout<< "number of arrays :\n " << shr_ds->getNumDistrArray() << std::endl;

    DistrArray *testArray = shr_ds->getArray(0);
    std::cout << "testArray = " << testArray->getDistr(10) << std::endl;
    std::cout << "=================================\n";*/
    
	return 0;
}
