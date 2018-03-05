#include <iostream>
#include <stdlib.h>
#include <cmath>
#include <vector>
#include <sstream>

#include "distributions/distribution.h"
#include "distributions/distribution_modeler.h"

//using namespace vinecopulib;

#include <sys/time.h>

#include "dataset/distr_array.h"
#include "utility_functions.h"
#include "mveddaIO.h"

using namespace std;
using namespace edda;


int main(int argc, char* argv[])
{
		FILE *fIn;
		float *retArray = new float[50*50*10*9];
		string fname = "/media/hazarika/WD2/ResearchDataStore2/codda/isabel/block_5x5x5/distrFiles/allGMM3/Pressure_GMM3_250x250x50_5x5x5.mvedda";
	    fIn = fopen(fname.c_str(),"rb");
	    if(!fIn)
	    {
	        fprintf(stderr, "Error opening file\n");
	        exit(13);
	    }
	    fread(retArray, sizeof(float), 50*50*10*9, fIn);
	    fclose(fIn);

	    for(int i=0; i< 50*50*10; i++)
	    {
	    	for(int j=0; j<9; j++)
	    	{
	    		cout << retArray[i*9+j] << ", ";
	    	}
	    	cout << endl;
	    }
}