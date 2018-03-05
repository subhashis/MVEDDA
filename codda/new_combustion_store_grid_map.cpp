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

float map_coordinates(int o_min, int o_max, int t_min, int t_max, int t_query){
	float o_diff = float(o_max - o_min);
	float t_diff = float(t_max - t_min);

	float ratio = float(t_query - t_min)/t_diff;

	float a = (ratio*o_diff) + o_min;

	return a;
}

bool isInside(float px, float py, float pz, float new_xmin, float new_xmax, float new_ymin, float new_ymax, float new_zmin, float new_zmax){
	bool flag(true);
	if(px < new_xmin || px > new_xmax)
	{
		flag = false;
		return flag;
	}
	if(py < new_ymin || py > new_ymax)
	{
		flag = false;
		return flag;
	}
	if(pz < new_zmin || pz > new_zmax)
	{
		flag = false;
		return flag;
	}
	return flag;
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
	string fname = "grid_map_" + grid_res + "_new.raw";
	fname = grid_output_path + fname;

	string fname1 = "block_map_" + grid_res + "_new.raw";
	fname1 = grid_output_path + fname1;
	//cout << fname << endl;


	int reduced_XDIM(0),reduced_YDIM(0),reduced_ZDIM(0);
	reduced_XDIM = XDIM/BLOCKX;
	reduced_YDIM = YDIM/BLOCKY;
	reduced_ZDIM = ZDIM/BLOCKZ;

	int perBlock_xsize = target_x/reduced_XDIM;
	int perBlock_ysize = target_y/reduced_YDIM;
	int perBlock_zsize = target_z/reduced_ZDIM;

	float *gridMap = new float[target_x*target_y*target_z*3];

	timestamp_t t0 = get_timestamp();
	for(int z=0; z<target_z; z++){
		for(int y=0; y<target_y; y++){
			for(int x=0; x<target_x; x++){
				float px = map_coordinates(0, XDIM-1, 0, target_x-1, x);
				float py = map_coordinates(0, YDIM-1, 0, target_y-1, y);
				float pz = map_coordinates(0, ZDIM-1, 0, target_z-1, z);

				gridMap[z*target_y*target_x*3 + y*target_x*3 + x*3 + 0] = px;
				gridMap[z*target_y*target_x*3 + y*target_x*3 + x*3 + 1] = py;
				gridMap[z*target_y*target_x*3 + y*target_x*3 + x*3 + 2] = pz;


				//cout << "<" << x << "," << y << "," << z << "> ---> " << "<" << px << "," << py << "," << pz << ">" << endl;
					
			}
		}
	}
	//now assign the grid location under influnce of a block (reduced data)
	//int max_count = -100, min_count = 10000;

	std::vector<std::vector<int> > list_gridPoints(reduced_ZDIM*reduced_YDIM*reduced_XDIM); 

	for(int z=0; z<reduced_ZDIM; z++){
		for(int y=0; y<reduced_YDIM; y++){
			for(int x=0; x<reduced_XDIM; x++){
				int idx = z*reduced_YDIM*reduced_XDIM + y*reduced_XDIM + x;

				//determine the x,y,z range of the block
		        int xmin = x*BLOCKX, xmax = x*BLOCKX + BLOCKX - 1;
		        int ymin = y*BLOCKY, ymax = y*BLOCKY + BLOCKY - 1;
		        int zmin = z*BLOCKZ, zmax = z*BLOCKZ + BLOCKZ - 1;

		        int new_xmin = xmin - 2, new_xmax = xmax + 2;
		        int new_ymin = ymin - 2, new_ymax = ymax + 2;
		        int new_zmin = zmin - 2, new_zmax = zmax + 2;

		        float t_xmin = map_coordinates(0, target_x-1, 0, XDIM-1, xmin);
				float t_xmax = map_coordinates(0, target_x-1, 0, XDIM-1, xmax);

				float t_ymin = map_coordinates(0, target_y-1, 0, YDIM-1, ymin);
				float t_ymax = map_coordinates(0, target_y-1, 0, YDIM-1, ymax);
				
				float t_zmin = map_coordinates(0, target_z-1, 0, ZDIM-1, zmin);
				float t_zmax = map_coordinates(0, target_z-1, 0, ZDIM-1, zmax);

				int t_xstart = int(floor(t_xmin)) - 3;
				int t_xend = int(ceil(t_xmax)) + 3;

				int t_ystart = int(floor(t_ymin)) - 3;
				int t_yend = int(ceil(t_ymax)) + 3;

				int t_zstart = int(floor(t_zmin)) - 3;
				int t_zend = int(ceil(t_zmax)) + 3;

				for(int zz=t_zstart; zz<=t_zend; zz++){
					for(int yy=t_ystart; yy<=t_yend; yy++){
						for(int xx=t_xstart; xx<=t_xend; xx++){

							if(isInside(xx,yy,zz,0,target_x-1,0,target_y-1,0,target_z-1))
							{
								float px = gridMap[zz*target_y*target_x*3 + yy*target_x*3 + xx*3 + 0];
								float py = gridMap[zz*target_y*target_x*3 + yy*target_x*3 + xx*3 + 1];
								float pz = gridMap[zz*target_y*target_x*3 + yy*target_x*3 + xx*3 + 2];

								if(isInside(px,py,pz,new_xmin,new_xmax,new_ymin,new_ymax,new_zmin,new_zmax)){
									//count_points++;
									int grid_id = zz*target_y*target_x + yy*target_x + xx;
									list_gridPoints[idx].push_back(grid_id); 
								}

							}
							
							
							//cout << "<" << x << "," << y << "," << z << "> ---> " << "<" << px << "," << py << "," << pz << ">" << endl;
								
						}
					}
				}
		        
		        //std::cout << "[" << z << "][" << y << "][" << x << "]" << "idx=" << idx <<endl ;
			}
		}
	} 

	timestamp_t t1 = get_timestamp();

	double secs = (t1 - t0) / 1000000.0L;
	cout << "time taken = " << secs << endl;
	dataWriter(fname, gridMap, target_x*target_y*target_z*3);
	vectorWriter_int(fname1, &list_gridPoints);

	delete[] gridMap;
	list_gridPoints.clear();

	return 0;
}