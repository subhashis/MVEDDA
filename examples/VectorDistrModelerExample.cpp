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
    DistributionModeler dm1(newW*newH*newD);   
    DistributionModeler dm2(newW*newH*newD);   
    DistributionModeler dm3(newW*newH*newD);   

    cout << "starting partitioning\n";
    int counter =0;

    for(size_t z=0; z<zdim; z += blockZdim)
    {
        for(size_t y=0; y<ydim; y += blockYdim)
        {
            for(size_t x=0; x<xdim; x += blockXdim)
            {

                float *data1 = new float[blockXdim*blockYdim*blockZdim];
                float *data2 = new float[blockXdim*blockYdim*blockZdim];
                float *data3 = new float[blockXdim*blockYdim*blockZdim];

                int i = 0;
                for(size_t zz=z; zz<z+blockZdim; zz++)
                {
                    for(size_t yy=y; yy<y+blockYdim; yy++)
                    {
                        for(size_t xx=x; xx<x+blockXdim; xx++)
                        {
                            data1[i] = inData1[zz*xdim*ydim + yy*xdim + xx];
                            data2[i] = inData2[zz*xdim*ydim + yy*xdim + xx];
                            data3[i] = inData3[zz*xdim*ydim + yy*xdim + xx];
                            i++;
                        }
                    }
                }
                std::cout << "dimensions: [" << z << "][" << y << "][" << x << "]\n";
                dm1.computeGMM(data1, blockXdim*blockYdim*blockZdim, 2, counter);
                dm2.computeGMM(data2, blockXdim*blockYdim*blockZdim, 2, counter);
                dm3.computeGMM(data3, blockXdim*blockYdim*blockZdim, 2, counter);
                counter++;
            }
        }
    }
    std::vector<DistrArray *> dVec;
    dVec.push_back(dm1.getDistrArray());
    dVec.push_back(dm2.getDistrArray());
    dVec.push_back(dm3.getDistrArray());

    Dataset<Real> *ds = new Dataset<Real> (new RegularCartesianGrid(newW, newH, newD), dVec);
    
  	return 0;
}
