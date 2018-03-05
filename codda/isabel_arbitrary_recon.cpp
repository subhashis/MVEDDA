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

typedef unsigned long long timestamp_t;

#define XDIM 250
#define YDIM 250
#define ZDIM 50

#define BLOCKX 5
#define BLOCKY 5
#define BLOCKZ 5

#define nVar 2
#define spatialVar 3

#define hist_nbin 32
#define numGaussian 3

#define numSample 500

#define DISTR_TYPE 1 //0:histogram, 1:GMM

//location to read the grid-map values from 
const string grid_output_path = "/home/hazarika/visData/mvcopula_output_data/isabel_MV_t20/res_250x250x50_5x5x5/grid_map/";

const string output_path = "/home/hazarika/visData/mvcopula_output_data/isabel_MV_t20/res_250x250x50_5x5x5/";
const string var0_filename = "Pressure_GMM3_250x250x50_5x5x5.mvedda";
const string var1_filename = "Velocity_GMM3_250x250x50_5x5x5.mvedda";
const string parameter_filename = output_path + "Pressure_Velocity_250x250x50_5x5x5.spcc";


static timestamp_t get_timestamp ()
{
  struct timeval now;
  gettimeofday (&now, NULL);
  return  now.tv_usec + (timestamp_t)now.tv_sec * 1000000;
}

void computeWeightedValues(float *mv_samples, std::vector<int> *grid_id_list, float *gridMap, float *pressureField, float *velocityField, float *sumOfWeightsField)
{
	int listSize = grid_id_list->size();
	for(int s=0; s<numSample; s++){
		float v1 = mv_samples[s*(nVar+spatialVar) + 0];
		float v2 = mv_samples[s*(nVar+spatialVar) + 1];

		float xs = mv_samples[s*(nVar+spatialVar) + nVar + 0];
		float ys = mv_samples[s*(nVar+spatialVar) + nVar + 1];
		float zs = mv_samples[s*(nVar+spatialVar) + nVar + 2];
		for(int l=0; l<listSize; l++){
			int grid_id = (*grid_id_list)[l];
			float px = gridMap[grid_id*3+0];
			float py = gridMap[grid_id*3+1];
			float pz = gridMap[grid_id*3+2];

			float dist2 = (xs-px)*(xs-px) + (ys-py)*(ys-py) + (zs-pz)*(zs-pz);
			float dist = sqrt(dist2);
			if(dist == 0.0)
			{
				cout << "zero encountered at idx = " << grid_id << endl;
				exit(11);
			}
			float weight = 1.0/dist;

			pressureField[grid_id] += (v1*weight);
			velocityField[grid_id] += (v2*weight);
			sumOfWeightsField[grid_id] += weight;
		}


	}
}

int main(int argc, char* argv[])
{

	if(argc < 4) {
        printf("[Invalid Arguments] : Check Usage!!\n");
        exit(13);
    }

	istringstream ss_x(argv[1]);
	int target_x;
	if (!(ss_x >> target_x))
	    cerr << "Invalid number " << argv[1] << '\n';

	istringstream ss_y(argv[2]);
	int target_y;
	if (!(ss_y >> target_y))
	    cerr << "Invalid number " << argv[2] << '\n';

	istringstream ss_z(argv[3]);
	int target_z;
	if (!(ss_z >> target_z))
	    cerr << "Invalid number " << argv[3] << '\n';

	string s1 = argv[1], s2 = argv[2], s3 = argv[3];
	string grid_res = s1 + "x" + s2 + "x" + s3;

	//generate the grid-map and block-map filenames for the user provided resolution
	string fname = "grid_map_" + grid_res + ".raw";
	fname = grid_output_path + fname;

	string fname1 = "block_map_" + grid_res + ".raw";
	fname1 = grid_output_path + fname1;
	//cout << fname << endl;

	string filename[nVar];
	//get the user provided filepath/name for distribution files for each variable
	filename[0] = output_path + var0_filename;
	//filename[0] = "/media/hazarika/WD2/ResearchDataStore2/codda/isabel/block_5x5x5/distrFiles/allGMM3/Pressure_GMM3_250x250x50_5x5x5.mvedda";
	filename[1] = output_path + var1_filename;


	int reduced_XDIM(0),reduced_YDIM(0),reduced_ZDIM(0);
	reduced_XDIM = XDIM/BLOCKX;
	reduced_YDIM = YDIM/BLOCKY;
	reduced_ZDIM = ZDIM/BLOCKZ;

	int perBlock_xsize = target_x/reduced_XDIM;
	int perBlock_ysize = target_y/reduced_YDIM;
	int perBlock_zsize = target_z/reduced_ZDIM;


	std:vector<DistrArray *> dVec;
	dVec.push_back(distrArrayReader_new(filename[0], reduced_ZDIM*reduced_YDIM*reduced_XDIM, 9));
	//dVec.push_back(distrArrayReader_new(filename[1], reduced_ZDIM*reduced_YDIM*reduced_XDIM, 9));



	/*int numPar = ((nVar+spatialVar)*(nVar+spatialVar-1))/2;
	float *pcc_parameters = new float[numPar*reduced_ZDIM*reduced_YDIM*reduced_XDIM];
	pcc_parameters  = pccParameterReader(parameter_filename, numPar*reduced_ZDIM*reduced_YDIM*reduced_XDIM);

	float *gridMap = dataReader(fname, target_x*target_y*target_z*3);
	std::vector<std::vector<int> > list_gridPoints(reduced_ZDIM*reduced_YDIM*reduced_XDIM);
	vecDataReader(fname1, &list_gridPoints);

	//target reconstructed fields + initialization
	float *pressureField = new float[target_x*target_y*target_z];
	float *velocityField = new float[target_x*target_y*target_z];
	float *sumOfWeightsField = new float[target_x*target_y*target_z];

	for(int i=0; i<target_x*target_y*target_z; i++){
		pressureField[i] = 0.0f;
		velocityField[i] = 0.0f;
		sumOfWeightsField[i] = 0.0f;
	}

	timestamp_t t0 = get_timestamp();

	for(int z=0; z<reduced_ZDIM; z++){
		for(int y=0; y<reduced_YDIM; y++){
			for(int x=0; x<reduced_XDIM; x++){

				int block_idx = z*reduced_YDIM*reduced_XDIM + y*reduced_XDIM + x;

				//determine the x,y,z range before transformation
		        int x1 = x*BLOCKX, x2 = x*BLOCKX + BLOCKX - 1;
		        int y1 = y*BLOCKY, y2 = y*BLOCKY + BLOCKY - 1;
		        int z1 = z*BLOCKZ, z2 = z*BLOCKZ + BLOCKZ - 1;

		        std::vector<dist::Variant> local_dVec(nVar);
				for(int v=0; v<nVar; v++)
				{
					local_dVec[v] =  dVec[v]->getDistr(block_idx);
				}
				float *mv_samples;
				mv_samples = getMVSamples(&pcc_parameters[block_idx*numPar], nVar, spatialVar, numSample, &local_dVec, x1, x2, y1, y2, z1, z2);
				local_dVec.clear();

				computeWeightedValues(mv_samples, &list_gridPoints[block_idx], gridMap, pressureField, velocityField, sumOfWeightsField);

				std::cout << "[" << z << "][" << y << "][" << x << "]" << "idx=" << block_idx << endl ;

				delete[] mv_samples;

			}
		}
	}

	timestamp_t t1 = get_timestamp();


	for(int i=0; i<target_x*target_y*target_z; i++){
		pressureField[i] /= sumOfWeightsField[i];
		velocityField[i] /= sumOfWeightsField[i];
	}

	double secs = (t1 - t0) / 1000000.0L;
	cout << "time taken = " << secs << endl;


	
	dataWriter(output_path + "output/recon_pressure_"+ grid_res +"_s500.raw", pressureField, target_x*target_y*target_z);
	dataWriter(output_path + "output/recon_velocity_"+ grid_res +"_s500.raw", velocityField, target_x*target_y*target_z);

	delete[] gridMap;
	list_gridPoints.clear();
	delete[] pressureField;
	delete[] velocityField;*/

	return 0;
}