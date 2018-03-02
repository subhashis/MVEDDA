#include <iostream>
#include <stdlib.h>
#include <cmath>
#include <vector>
#include <sstream>

#include <vtkSmartPointer.h>
#include <vtkPolyDataMapper.h>
#include <vtkXMLStructuredGridReader.h>
#include <vtkXMLStructuredGridWriter.h>
#include <vtkStructuredGridGeometryFilter.h>
#include <vtkDataSet.h>
#include <vtkStructuredGrid.h>
#include <vtkFloatArray.h>
#include <vtkPointData.h>
#include <vtkDataArray.h>

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

#define BLOCKX 5
#define BLOCKY 5
#define BLOCKZ 5

#define nVar 3
#define spatialVar 3

#define hist_nbin 32
#define numGaussian 3

#define numSample 500


#define DISTR_TYPE 1 //0:histogram, 1:GMM

//const string grid_output_path = "/home/hazarika/visData/mvcopula_output_data/isabel_MV_t20/res_250x250x50_5x5x5/grid_map/";

static timestamp_t get_timestamp ()
{
  struct timeval now;
  gettimeofday (&now, NULL);
  return  now.tv_usec + (timestamp_t)now.tv_sec * 1000000;
}

void write_vts(vtkSmartPointer<vtkStructuredGrid> newSgrid, const string &vts_output_path)
{
    
    //string filename = "/media/hazarika/WD3/turbine_insitu/test_passage_18/test_passage_18_0_0_recon_test4_single_mveddafile_read.vts";
    //string filename = vts_output_path + "out_"+procId+"_15_recon.vts";
    
    //Write file: grid point centered data            
    vtkSmartPointer<vtkXMLStructuredGridWriter> writer = vtkSmartPointer<vtkXMLStructuredGridWriter>::New();
    writer->SetFileName(vts_output_path.c_str());
    writer->SetInputData(newSgrid);
    writer->SetEncodeAppendedData(0);
    writer->Write();     
}

float map_coordinates(int o_min, int o_max, int t_min, int t_max, int t_query){
  float o_diff = float(o_max - o_min);
  float t_diff = float(t_max - t_min);

  float ratio = float(t_query - t_min)/t_diff;

  float a = (ratio*o_diff) + o_min;

  return a;
}

bool isInside(float px, float py, float pz, float new_xmin, float new_xmax, float new_ymin, float new_ymax, float new_zmin, float new_zmax)
{
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

void computeWeightedValues(float *mv_samples, std::vector<int> *grid_id_list, float *gridMap, float *entropyField, float *uvelField, float *temperatureField, float *sumOfWeightsField)
{
	int listSize = grid_id_list->size();
	//cout << "listSize=" << listSize << endl;
	for(int s=0; s<numSample; s++){
		float v1 = mv_samples[s*(nVar+spatialVar) + 0];
		float v2 = mv_samples[s*(nVar+spatialVar) + 1];
		float v3 = mv_samples[s*(nVar+spatialVar) + 2];

		
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

			


			entropyField[grid_id] += (v1*weight);
			uvelField[grid_id] += (v2*weight);
			temperatureField[grid_id] += (v3*weight);
			sumOfWeightsField[grid_id] += weight;

			// if(isNaN(entropyField[grid_id]))
			// {
			// 	cout << "at v1" << endl;
			// 	//cout << "{" << v1 << "," << v2 << "," << v3 << "," << xs <<  ", " << ys << "," << zs << "}" << endl ;
			// 	//cout << "weight=" << weight << " dist=" << dist << " entropyField[" << grid_id << "]=" << entropyField[grid_id] << " v1*weight=" << v1*weight << endl;
			// 	exit(14);
			// }
			// if(isNaN(uvelField[grid_id]))
			// {
			// 	cout << "at v2" << endl;
			// 	exit(14);
			// }
			// if(isNaN(uvelField[grid_id]))
			// {
			// 	cout << "at v3" << endl;
			// 	exit(14);
			// }
		}


	}
}

void computeConditionProb(float *mv_samples, std::vector<int> *grid_id_list, float *gridMap, float *conditionField, float *sumOfTotalSamples)
{
	int listSize = grid_id_list->size();
	//cout << "listSize=" << listSize << endl;
	int count = 0;
	for(int s=0; s<numSample; s++){
		float v1 = mv_samples[s*(nVar+spatialVar) + 0];
		float v2 = mv_samples[s*(nVar+spatialVar) + 1];
		float v3 = mv_samples[s*(nVar+spatialVar) + 2];
		if(v1 >= 0.8 && v2 <= -0.05)
		{
			count++;
		}
	}
	for(int l=0; l<listSize; l++){
		int grid_id = (*grid_id_list)[l];
		
		conditionField[grid_id] += float(count);
		sumOfTotalSamples[grid_id] += float(numSample);
	}
}


int main(int argc, char* argv[])
{
	if(argc < 2) {
        printf("[Invalid Arguments] : Check Usage!!\n");
        exit(13);
    }

	string procId = argv[1];
	string timeId = "20602";

	vtkSmartPointer<vtkXMLStructuredGridReader> reader1 = vtkSmartPointer<vtkXMLStructuredGridReader>::New();
	string full_vts_file = "/media/hazarika/WD2/ResearchDataStore2/codda/turbine/turbine_full_annulus_grids/full_grid_data/time20602/out_"+procId+"_20602.vts";
	reader1->SetFileName(full_vts_file.c_str());
	reader1->Update();

	vtkSmartPointer<vtkStructuredGrid> full_sgrid = vtkStructuredGrid::SafeDownCast(reader1->GetOutput());
	//int numPoints1 = full_sgrid->GetNumberOfPoints();
	//cout << numPoints1 << endl;
	int full_extent[6];
	full_sgrid->GetExtent(full_extent);
	int XDIM = full_extent[1] + 1;
	int YDIM = full_extent[3] + 1;
	int ZDIM = full_extent[5] + 1;

	vtkSmartPointer<vtkXMLStructuredGridReader> reader2 = vtkSmartPointer<vtkXMLStructuredGridReader>::New();
	string target_vts_file = "/media/hazarika/WD2/ResearchDataStore2/codda/turbine/turbine_full_annulus_grids/grid_data/time20602/out_"+procId+"_20602.vts";
	//string target_vts_file = "/media/hazarika/WD2/ResearchDataStore2/codda/turbine/turbine_full_annulus_grids/full_grid_data/time20602/out_"+procId+"_20602.vts";
	reader2->SetFileName(target_vts_file.c_str());
	reader2->Update();

	vtkSmartPointer<vtkStructuredGrid> target_sgrid = vtkStructuredGrid::SafeDownCast(reader2->GetOutput());
	//int numPoints2 = target_sgrid->GetNumberOfPoints();
	//cout << numPoints2 << endl;
	int target_extent[6];
	target_sgrid->GetExtent(target_extent);
	int target_x = target_extent[1] + 1;
	int target_y = target_extent[3] + 1;
	int target_z = target_extent[5] + 1;

	cout << "full res : [" << XDIM << "," << YDIM << "," << ZDIM << "]" << endl;
	cout << "target res : [" << target_x << "," << target_y << "," << target_z << "]" << endl;
	

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
	timestamp_t t1 = get_timestamp();

	//partitioning scheme for turbine data: last block/dimension needs to be handled separatetly
	int newX = XDIM/BLOCKX;
	int newY = YDIM/BLOCKY;
	int newZ = ZDIM/BLOCKZ;

	int XDIM1 = newX*BLOCKX;
	int YDIM1 = newY*BLOCKY;
	int ZDIM1 = newZ*BLOCKZ;

	

	std::vector<std::vector<int> > list_gridPoints(newX*newY*newZ);
	int blockId = 0;
	timestamp_t t2 = get_timestamp();
	for(size_t z=0; z<ZDIM1; z+=BLOCKZ){
		for(size_t y=0; y<YDIM1; y+=BLOCKY){
	    	for(size_t x=0; x<XDIM1; x+=BLOCKX){
	    		int x_start = x;
	    		int y_start = y;
	    		int z_start = z;

	    		int x_end = x+BLOCKX-1;
	    		int y_end = y+BLOCKY-1;
	    		int z_end = z+BLOCKZ-1;

	    		if(x/BLOCKX == newX-1){
	    		  if(x_end < XDIM-1)
	    		  	x_end = XDIM-1;
	    		}
	    		if(y/BLOCKY == newY-1){
	    		  if(y_end < YDIM-1)
	    		    y_end = YDIM-1;
	    		}
	    		if(z/BLOCKZ == newZ-1){
	    		  if(z_end < ZDIM-1)
	    		    z_end = ZDIM-1;
	    		}

	    		int xB_size = x_end - x_start + 1;
	    		int yB_size = y_end - y_start + 1;
	    		int zB_size = z_end - z_start + 1;

	    		int new_xmin = x_start - 2, new_xmax = x_end + 2;
	    		int new_ymin = y_start - 2, new_ymax = y_end + 2;
	    		int new_zmin = z_start - 2, new_zmax = z_end + 2;

	    		//int count_points = 0;

	    		for(int zz=0; zz<target_z; zz++){
	    		  for(int yy=0; yy<target_y; yy++){
	    		    for(int xx=0; xx<target_x; xx++){
	    		      float px = gridMap[zz*target_y*target_x*3 + yy*target_x*3 + xx*3 + 0];
	    		      float py = gridMap[zz*target_y*target_x*3 + yy*target_x*3 + xx*3 + 1];
	    		      float pz = gridMap[zz*target_y*target_x*3 + yy*target_x*3 + xx*3 + 2];

	    		      if(isInside(px,py,pz,new_xmin,new_xmax,new_ymin,new_ymax,new_zmin,new_zmax)){
	    		        //count_points++;
	    		        int grid_id = zz*target_y*target_x + yy*target_x + xx;
	    		        list_gridPoints[blockId].push_back(grid_id); 
	    		      }
	    		      
	    		      //cout << "<" << x << "," << y << "," << z << "> ---> " << "<" << px << "," << py << "," << pz << ">" << endl;
	    		        
	    		    }
	    		  }
	    		}
	    		blockId++;
	    	}
	    }
	}

	timestamp_t t3 = get_timestamp();

	

	//respective filenames for different processors
	string distr_filename[1];
	string para_filename[1];

	distr_filename[0] = "/media/hazarika/WD2/ResearchDataStore2/codda/turbine/turbine_full_annulus_grids/distrFiles/time"+timeId+"/allvar_GMM3.single_mvedda_"+procId+"_" + timeId;
	// distr_filename[1] = filepath + "allvar_GMM3.single_mvedda_1_15";
	// distr_filename[2] = filepath + "allvar_GMM3.single_mvedda_2_15";
	// distr_filename[3] = filepath + "allvar_GMM3.single_mvedda_3_15";

	para_filename[0] = "/media/hazarika/WD2/ResearchDataStore2/codda/turbine/turbine_full_annulus_grids/distrFiles/time"+timeId+"/allvar_corr.spcc_"+procId+"_" + timeId;
	// para_filename[1] = filepath + "allvar_corr.spcc_1_15";
	// para_filename[2] = filepath + "allvar_corr.spcc_2_15";
	// para_filename[3] = filepath + "allvar_corr.spcc_3_15";


	//read the distribution files and parameter files
	std:vector<DistrArray *> dVec;
	int num_distr_component = (nVar*9);
	int per_filecount = newX*newY*newZ;
	int start = 0;
	int end = start + per_filecount*9;
	int fn_size = 1;
	dVec.push_back(distrArrayReader_new_multifile(distr_filename, fn_size, per_filecount, num_distr_component, start, end, 9));
	//distrArrayReader_new_multifile(distr_filename, fn_size, reduced_ZDIM*reduced_YDIM*reduced_XDIM, num_distr_component, start, end, 9);
	start = end;
	end += per_filecount*9;
	dVec.push_back(distrArrayReader_new_multifile(distr_filename, fn_size, per_filecount, num_distr_component, start, end, 9));
	//distrArrayReader_new_multifile(distr_filename, fn_size, reduced_ZDIM*reduced_YDIM*reduced_XDIM, num_distr_component, start, end, 9);
	start = end;
	end += per_filecount*9;
	dVec.push_back(distrArrayReader_new_multifile(distr_filename, fn_size, per_filecount, num_distr_component, start, end, 9));
	

	int numPar = ((nVar+spatialVar)*(nVar+spatialVar-1))/2;
	float *pcc_parameters = new float[numPar*per_filecount];
	pcc_parameters  = pccParameterReader_multifile(para_filename, fn_size, numPar*per_filecount);


	float *conditional_array = new float[target_x*target_y*target_z];
	float *sumOfTotalSamplesField = new float[target_x*target_y*target_z];

	for(int i=0; i<target_x*target_y*target_z; i++)
	{
		conditional_array[i] = 0.0;
		sumOfTotalSamplesField[i] = 0.0;
	}


	blockId = 0;
	timestamp_t t4 = get_timestamp();
	for(size_t z=0; z<ZDIM1; z+=BLOCKZ){
		for(size_t y=0; y<YDIM1; y+=BLOCKY){
	    	for(size_t x=0; x<XDIM1; x+=BLOCKX){
	    		int x_start = x;
	    		int y_start = y;
	    		int z_start = z;

	    		int x_end = x+BLOCKX-1;
	    		int y_end = y+BLOCKY-1;
	    		int z_end = z+BLOCKZ-1;

	    		if(x/BLOCKX == newX-1){
	    		  if(x_end < XDIM-1)
	    		  	x_end = XDIM-1;
	    		}
	    		if(y/BLOCKY == newY-1){
	    		  if(y_end < YDIM-1)
	    		    y_end = YDIM-1;
	    		}
	    		if(z/BLOCKZ == newZ-1){
	    		  if(z_end < ZDIM-1)
	    		    z_end = ZDIM-1;
	    		}

	    		int xB_size = x_end - x_start + 1;
	    		int yB_size = y_end - y_start + 1;
	    		int zB_size = z_end - z_start + 1;

	    		std::vector<dist::Variant> local_dVec(nVar);
	    		local_dVec[0] = dVec[0]->getDistr(blockId);
	    		local_dVec[1] = dVec[1]->getDistr(blockId);
	    		local_dVec[2] = dVec[2]->getDistr(blockId);

	    		float *mv_samples;
	    		mv_samples = getMVSamples(&pcc_parameters[blockId*numPar], nVar, spatialVar, numSample, &local_dVec, x_start, x_end, y_start, y_end, z_start, z_end);
	    		local_dVec.clear();

	    		
	    		//cout << "<" << x << "," << y << "," << z << ">"<<endl;
	    		//computeWeightedValues(mv_samples, &list_gridPoints[blockId], gridMap, entropy_recon_array, uvel_recon_array, temperature_recon_array, sumOfWeightsField);
				computeConditionProb(mv_samples, &list_gridPoints[blockId], gridMap, conditional_array, sumOfTotalSamplesField);
				
				blockId++;
				delete[] mv_samples;
	    	}
	    }
	}
	timestamp_t t5 = get_timestamp();

	delete[] pcc_parameters;

	for(int i=0; i<target_x*target_y*target_z; i++)
	{		
		conditional_array[i] /= sumOfTotalSamplesField[i];		
	}

	vtkSmartPointer<vtkFloatArray> condArray = vtkSmartPointer<vtkFloatArray>::New();
	condArray->SetName("conditionField");
	condArray->SetNumberOfComponents(1);
	condArray->SetNumberOfValues(target_x*target_y*target_z);

	for(vtkIdType i=0; i<target_x*target_y*target_z; i++){
	  float f1 = conditional_array[int(i)];
	  condArray->SetValue(i,f1);
	}

	target_sgrid->GetPointData()->AddArray(condArray);


	delete[] conditional_array;
	delete[] sumOfTotalSamplesField;

	string vts_output_path = "/media/hazarika/WD2/ResearchDataStore2/codda/turbine/turbine_full_annulus_grids/grid_data_output_conditional/time"+timeId+"/out_"+procId+"_"+timeId+"_condition.vts";
	//string vts_output_path = "/media/hazarika/WD2/ResearchDataStore2/codda/turbine/turbine_full_annulus_grids/full_grid_data_output/time"+timeId+"/out_"+procId+"_"+timeId+"_recon.vts";
	write_vts(target_sgrid, vts_output_path);

	double secs0 = (t1 - t0) / 1000000.0L;
	double secs1 = (t3 - t2) / 1000000.0L;
	double secs2 = (t5 - t4) / 1000000.0L;
	cout << "grid map time = " << secs0 << endl;
	cout << "block map time = " << secs1 << endl;
	cout << "recon time = " << secs2 << endl;
	

	return 0;
}

