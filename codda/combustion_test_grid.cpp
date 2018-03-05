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

#define XDIM 480
#define YDIM 720
#define ZDIM 120

#define BLOCKX 5
#define BLOCKY 5
#define BLOCKZ 5

#define nVar 3
#define spatialVar 3

#define hist_nbin 32
#define numGaussian 3

#define numSample 500

#define DISTR_TYPE 1 //0:histogram, 1:GMM

const string grid_output_path = "/media/hazarika/WD2/ResearchDataStore2/codda/combustion/block_res_5x5x5/grid_map/";

static timestamp_t get_timestamp ()
{
  struct timeval now;
  gettimeofday (&now, NULL);
  return  now.tv_usec + (timestamp_t)now.tv_sec * 1000000;
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


	string new_fname = "grid_map_" + grid_res + "_new.raw";
	new_fname = grid_output_path + new_fname;

	string new_fname1 = "block_map_" + grid_res + "_new.raw";
	new_fname1 = grid_output_path + new_fname1;


	int reduced_XDIM(0),reduced_YDIM(0),reduced_ZDIM(0);
	reduced_XDIM = XDIM/BLOCKX;
	reduced_YDIM = YDIM/BLOCKY;
	reduced_ZDIM = ZDIM/BLOCKZ;
	

	float *gridMap = dataReader(fname, target_x*target_y*target_z*3);
	std::vector<std::vector<int> > list_gridPoints(reduced_ZDIM*reduced_YDIM*reduced_XDIM);
	vecDataReader(fname1, &list_gridPoints);

	float *new_gridMap = dataReader(new_fname, target_x*target_y*target_z*3);
	std::vector<std::vector<int> > new_list_gridPoints(reduced_ZDIM*reduced_YDIM*reduced_XDIM);
	vecDataReader(new_fname1, &new_list_gridPoints);


	//check grid_map
	int grid_map_mismatch = 0;
	for(int i=0; i<target_x*target_y*target_z*3; i++)
	{
		if(gridMap[i] != new_gridMap[i])
		{
			grid_map_mismatch++;
		}
	}

	cout << "grid_map_mismatch = " << grid_map_mismatch << endl;

	int vect_size_missmatch = 0;
	int vect_value_missmatch = 0;
	for(int i=0; i<reduced_ZDIM*reduced_YDIM*reduced_XDIM; i++)
	{
		if(list_gridPoints[i].size() != new_list_gridPoints[i].size())
			vect_size_missmatch++;
		else
		{
			for(int j=0; j<list_gridPoints[i].size(); j++)
			{
				if(list_gridPoints[i][j] != new_list_gridPoints[i][j])
					vect_value_missmatch++;
			}
		}
	}

	cout << "vect_size_missmatch = " << vect_size_missmatch << endl;
	cout << "vect_value_missmatch = " << vect_value_missmatch << endl;




	

	return 0;
}