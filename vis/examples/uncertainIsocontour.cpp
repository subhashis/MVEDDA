// Copyright 2015 The Edda Authors. All rights reserved.
// Use of this source code is governed by a MIT-style license that can be
// found in the LICENSE file.

#include <string>
#include <cstdlib>
#include <cstdio>
#include <iostream>

#include "distributions/gaussian.h"
#include "distributions/distribution.h"
#include "io/path.h"
#include "filters/level_crossing_prob.h"
#include "io/edda_reader.h"

using namespace std;
using namespace edda;

void process_edda_file(string edda_file, float isov)
{
  std::vector<dist::GMM> gmmArray;
  int dim[3];
  int* dims;
  if (getFileExtension(edda_file).compare("edda")==0) {
	shared_ptr<Dataset<Real>> shr_ds2 = loadEddaScalarDataset_noneVTK(edda_file.c_str());
	dims = shr_ds2->getDimension();
	dim[0] = dims[0];
	dim[1] = dims[1];
	dim[2] = dims[2];

	for(int z=0; z<dims[2]; z++) {
	  for(int y=0; y<dims[1]; y++) {
		for(int x=0; x<dims[0]; x++) {
		  dist::Variant distr = shr_ds2->at_comp_distr_new(x, y, z)[0];
		  dist::GMM curDist = boost::get<dist::GMM >(distr);
		  gmmArray.push_back(curDist);
		}
   	  }
	}
  } else {
    printf("File format not supported\n");
    exit(1);
  }

  // compute the prob volume using level-crossing
  shared_ptr<NdArray<float>> out_ndarray;
  ReturnStatus r = levelCrossingProb(gmmArray.begin(), dim, isov, out_ndarray);
  if (r!=ReturnStatus::SUCCESS) {
    return ;
  }
  // save prob volume to a vector
  std::vector<float> buffer;
  std::vector<float> mean;
  for(int i=0; i<out_ndarray->get_num_of_elems(); i++) {
	float tmp = out_ndarray->get_val(i);
	if(tmp>1) tmp=1;
	if(tmp<0) tmp=0;
	buffer.push_back(tmp);
  }
  // compute the mean volume
  for(int i=0; i<dim[0]*dim[1]*dim[2]; i++) {
	float curmean = 0;
	float curweight = 0;
	for(int model=0; model<gmmArray[i].getNumComponenets(); model++) {
	  curmean += gmmArray[i].models[model].w*gmmArray[i].models[model].m;
	  curweight += gmmArray[i].models[model].w; 
	}
	mean.push_back(curmean/curweight);
  }

  // save as raw file
  printf("writing prob.raw with %d elements\n", buffer.size());
  FILE* fp = fopen("prob.raw", "w");
  fwrite(&buffer[0], sizeof(float), buffer.size(), fp);
  fclose(fp);
  printf("writing mean.raw with %d elements\n", mean.size());
  fp = fopen("mean.raw", "w");
  fwrite(&mean[0], sizeof(float), mean.size(), fp);
  fclose(fp);
}

int main(int argc, char** argv) {
  cout << "isoProbField <edda file> <iso-value>" << endl;
  if (argc<=2)
    return -1;
  string input_file = argv[1];
  float isov = atof(argv[2]);
  
  cout <<"processing the edda file now"<<endl;
  if (getFileExtension(input_file).compare("edda")==0) {
    process_edda_file(input_file, isov);
  }
  else {
	cout << "please provide files with .edda extension" <<endl;
  } 
}
