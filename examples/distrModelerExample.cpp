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
	string filename;
	int xdim;
	int ydim;
	int zdim;
	int blockXdim;
	int blockYdim;
	int blockZdim;

	if(argc < 8)
	{
		cout << "USAGE: check the parameter list!!\n";
		exit(11);
	}
	else
	{
		filename = argv[1];
		xdim = atoi(argv[2]);
		ydim = atoi(argv[3]);
		zdim = atoi(argv[4]);
		blockXdim = atoi(argv[5]);
		blockYdim = atoi(argv[6]);
		blockZdim = atoi(argv[7]);
	}

	float *inData;
	inData = new float[xdim*ydim*zdim];	
    
    FILE *fIn = fopen(filename.c_str(),"rb");
    if(!fIn)
    {
        fprintf(stderr, "Error opening file\n");
        exit(13);
    }
    size_t readStatus = fread(inData, sizeof(float), xdim*ydim*zdim, fIn);
    fclose(fIn);

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
	std::cout << "successfully read the file\n";
    //edda data modeling
	DistributionModeler dm(newW*newH*newD);   
	
	int counter =0;
    for(size_t z=0; z<zdim; z += blockZdim)
    {
    	for(size_t y=0; y<ydim; y += blockYdim)
    	{
    		for(size_t x=0; x<xdim; x += blockXdim)
    		{

				float *data;
				data = new float[blockXdim*blockYdim*blockZdim];
    			int i = 0;
    			for(size_t zz=z; zz<z+blockZdim; zz++)
    			{
    				for(size_t yy=y; yy<y+blockYdim; yy++)
    				{
    					for(size_t xx=x; xx<x+blockXdim; xx++)
    					{
    						data[i] = inData[zz*xdim*ydim + yy*xdim + xx];
    						i++;
    					}
    				}
    			}

                std::cout << "dimensions: [" << z << "][" << y << "][" << x << "]\n";

    			
				dm.computeGMM(data, blockXdim*blockYdim*blockZdim, 3, counter);

				counter++;
    		}
    	}
    }
    
    std::vector<DistrArray *> dVec;
    dVec.push_back(dm.getDistrArray());

    Dataset<Real> *ds = new Dataset<Real> (new RegularCartesianGrid(newW, newH, newD), dVec);
	
    /*//edda ensemble data modeling
    DistributionModelerNew dm1(newW*newH*newD);   
    DistributionModelerNew dm2(newW*newH*newD);   
    DistributionModelerNew dm3(newW*newH*newD);   

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
                dm1.computeGMM(data1, blockXdim*blockYdim*blockZdim, counter);
                dm2.computeGMM(data1, blockXdim*blockYdim*blockZdim, counter);
                dm3.computeGMM(data1, blockXdim*blockYdim*blockZdim, counter);
                counter++;
            }
        }
    }
    std::vector<DistrArray *> dVec;
    dVec.push_back(dm1.getDistrArray());
    dVec.push_back(dm2.getDistrArray());
    dVec.push_back(dm3.getDistrArray());

    Dataset<Real> *ds = new Dataset<Real> (new RegularCartesianGrid(newW, newH, newD), dVec);*/

    
  	return 0;
}
