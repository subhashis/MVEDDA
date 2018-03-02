#include <iostream>
#include <string>
#include <cstdio>
#include "edda.h"
#include "distributions/gaussian.h"
#include "distributions/distribution.h"
#include "io/edda_vtk_reader.h"
#include "dataset/dataset.h"

#include "io/edda_writer.h"
#include "io/edda_reader.h"

using namespace std;
using namespace edda;

//this file use old VTK IO type, which is not supported anymore. This file may be deleted in the future

int main(int argc, char **argv)
{
  srand(time(NULL));  // random seeding

  string filename;
  string arrayNamePrefix;
  if (argc>1) {
    filename = string(argv[1]);
    arrayNamePrefix = string(argv[2]);
  } else {
    cout << "Loading sample file" << endl;
    filename = string(SAMPLE_DATA_PATH) + "/isabel_pressure_small.vti";
  }

  //filename = string(SAMPLE_DATA_PATH) + "/isabel_pressure_small.vti";

  //original VTK loader
  shared_ptr<Dataset<Real> > dataset = loadEddaScalarDataset(filename, arrayNamePrefix);

  //write edda dataset with our writer
  writeEddaDataset(dataset, "testData.edda");
  
  //load edda dataset with our reader
  shared_ptr<Dataset<Real> > dataset2 = loadEddaScalarDataset_noneVTK("testData.edda");

  cout << "data loaded by VTK reader: " << endl;

  {
	  VECTOR3 pos;
	  Real value;
	  dist::Variant distr;
	  int i;

	  pos = VECTOR3(2.1, 2.1, 2.1);
	  dataset->at_phys(pos, value);
	  cout << pos << ": " << value << endl;

	  value = dataset->at_comp(5, 5, 5);
	  cout << "at_comp(5,5,5) : " << value << endl;

	  distr = dataset->at_comp_distr(5, 5, 5);
	  cout << "at_comp(5,5,5) : " << distr << endl;

	  distr = dataset->at_comp_distr(2, 2, 2);
	  cout << "at_comp(2,2,2) : " << distr << endl;
  }

  cout << endl << "data loaded by our reader: " << endl;

  {
	  VECTOR3 pos;
	  Real value;
	  dist::Variant distr;
	  int i;

	  pos = VECTOR3(2.1, 2.1, 2.1);
	  dataset2->at_phys(pos, value);
	  cout << pos << ": " << value << endl;

	  value = dataset2->at_comp(5, 5, 5);
	  cout << "at_comp(5,5,5) : " << value << endl;

	  distr = dataset2->at_comp_distr(5, 5, 5);
	  cout << "at_comp(5,5,5) : " << distr << endl;

	  distr = dataset2->at_comp_distr(2, 2, 2);
	  cout << "at_comp(2,2,2) : " << distr << endl;
  }

  //one by one check only for GMM5!!
  if (dataset->getArray()->getDistrName() == "GaussianMixture5"){
	  int* dims = dataset2->getDimension();
	  double dif = 0;
	  for (int k = 0; k < dims[2]; k++){
		  for (int j = 0; j < dims[1]; j++){
			  for (int i = 0; i < dims[0]; i++){
				  dist::Variant distr = dataset->at_comp_distr(i, j, k);
				  dist::Variant distr2 = dataset2->at_comp_distr(i, j, k);
				  dist::GaussianMixture<5> curDist1 = boost::get<dist::GaussianMixture<5> >(distr);
				  dist::GaussianMixture<5> curDist2 = boost::get<dist::GaussianMixture<5> >(distr2);
				  for (int model = 0; model < 5; model++){
					  dif = dif + abs(curDist1.models[model].m - curDist2.models[model].m)
						  + abs(curDist1.models[model].v - curDist2.models[model].v)
						  + abs(curDist1.models[model].w - curDist2.models[model].w);
				  }
			  }
		  }
	  }
	  cout << "the differences between the old vtk format and the new format is: " << dif << endl;
  }
  else{
	  cout << "one by one check is not performed " << endl;
  }
  return 0;
}
