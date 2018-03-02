#include <iostream>
#include <fstream>
#include <string>
#include <cstdio>
#include <cfloat>

#include "edda.h"
#include "distributions/gaussian.h"
#include "distributions/distribution.h"
#include "io/edda_vtk_reader.h"
#include "dataset/dataset.h"

//this file use old VTK IO type, which is not supported anymore. This file may be deleted in the future

using namespace std;
using namespace edda;

int main(int argc, char **argv)
{
  srand(time(NULL));  // random seeding

  VECTOR3 pos, pos2, pos3;
  pos = VECTOR3(2, 2, 2);
  pos2 = VECTOR3(2.1, 2.1, 2.1);
  pos3 = VECTOR3(0.0, 0.0, 0.0);

  cout << "Loading sample file" << endl;

  //for testing curvilinear grid
  string filename = string(SAMPLE_DATA_PATH) + "/out_92651_0.vts";
  pos2 = VECTOR3(-0.02, -0.4312, 0.006); //for out_92651_0.vts //vert id 1474 (4,7,3)

  //string filename = string(SAMPLE_DATA_PATH) + "/out_0_9602.vts";
  //pos2 = VECTOR3(-0.03, 0.416, 0.208);
  
  // load data with random sampling
  shared_ptr<Dataset<Real> > dataset = loadEddaScalarDataset(filename, "");

  
  Real value;
  int i,j,k;
  i = (int)pos[0];
  j = (int)pos[1];
  k = (int)pos[2];
  
  value = dataset->at_comp(i,j,k);
  cout << "at_comp " << i << "," << j << "," << k << " : " << value << endl;

  if (dataset->at_phys(pos2, value) == SUCCESS)
	  cout << pos2 << ": " << value << endl;
  else
	  cout << pos2 << ": " << "fail to get value at given position" << endl;

  
  if (dataset->at_phys(pos3, value) == SUCCESS)
	  cout << pos3 << ": " << value << endl;
  else
	  cout << pos3 << ": " << "fail to get value at given position" << endl;

  dist::Variant distr = dataset->at_comp_distr(i,j,k);
  cout << "at_comp_distr " << i << "," << j << "," << k << " : " << distr << endl;

  VECTOR3 minB, maxB;
  dataset->getGrid()->boundary(minB, maxB);

  const int numSampleX = 500, numSampleZ = 500;
  bool* successfullySampled = new bool[numSampleZ*numSampleX];
  float* sampleResults = new float[numSampleZ*numSampleX];
  
  float vmax = -FLT_MAX, vmin = FLT_MAX;
  float selectedY = (maxB[1] - minB[1]) / 2 + minB[1];
  for (int j = 0; j < numSampleZ; j++){
	  for (int i = 0; i < numSampleX; i++){
		VECTOR3 pos = VECTOR3((maxB[0] - minB[0]) / (numSampleX - 1)*i + minB[0], selectedY, (maxB[2] - minB[2]) / (numSampleZ - 1)*j + minB[2]);
		  Real value;
		  if (dataset->at_phys(pos, value) == SUCCESS){
			  sampleResults[j*numSampleX+i] = value;
			  successfullySampled[j*numSampleX + i] = true;
			  vmax = max(vmax, value);
			  vmin = min(vmin, value);
		  }
		  else{
			  successfullySampled[j*numSampleX + i] = false;
		  }
	  }
  }

  //at places where sampling is failed, the value is set to vmin - (vmax - vmin) as an indicator 
  for (int j = 0; j < numSampleZ; j++){
	  for (int i = 0; i < numSampleX; i++){
		  if (!successfullySampled[j*numSampleX + i]){
			  sampleResults[j*numSampleX + i] = vmin - (vmax - vmin);
		  }
	  }
  }

  const char* outputFilename = "crossSectionAtCenterY.raw";
  ofstream myfile;
  myfile.open(outputFilename, ios::out | ios::binary);
  myfile.write((const char*)sampleResults, sizeof(float)*numSampleZ*numSampleX);
  myfile.close();

  cout << "uniformly sampled at the cross section at y coordinate " << selectedY << endl;
  cout << "the result is saved to " << outputFilename << endl;
  
  delete[] successfullySampled;
  delete[] sampleResults;

  return 0;
}
