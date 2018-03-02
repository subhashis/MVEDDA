#include <iostream>

#include <thrust/device_vector.h>
#include <thrust/transform.h>
#include <thrust/sequence.h>
#include <thrust/copy.h>
#include <thrust/fill.h>
#include <thrust/replace.h>
#include <thrust/functional.h>
#include <iostream>

#include "core/shared_ary.h"
#include "distributions/distribution.h"
#include "distributions/gaussian.h"
#include "core/interpolator.h"
#include "test_common.h"
#include "core/thrust_random_sample.h"
#include "distributions/gaussian_mixture.h"


using namespace std;
using namespace edda;

int main(int argc, char* argv[]) {
  dist::Gaussian g1(1.,1.);
  dist::Gaussian g2(2.,2.);
  dist::Gaussian g3=g1+g2;

  {
    thrust::host_vector<dist::Gaussian> h_array(100);
    int i;
    for (i=0; i<100; i++)
      h_array[i] = dist::Gaussian(0,1);

    thrust::device_vector<dist::Gaussian > d_array = h_array;

    thrust::device_vector<Real > d_out(h_array.size());

    randomSampleField(d_array.begin(), d_array.end(), d_out.begin());

    thrust::host_vector<Real > h_out = d_out;

    for (i=0; i<h_out.size(); i++)
      cout << "Gaussian sampling output: " << h_out[i] << endl;
  }
  {
    vector<dist::GMMTuple> models;
    dist::GMMTuple tuple;
    tuple.m = 0;
    tuple.v = 1;
    tuple.w = 1;
    models.push_back(tuple);
    tuple.m = 1;
    tuple.v = 1;
    tuple.w = 1;
    models.push_back(tuple);
    tuple.m = 2;
    tuple.v = 1;
    tuple.w = 1;
    models.push_back(tuple);
    thrust::host_vector<dist::GaussianMixture<3> > h_array(100);
    int i;
    for (i=0; i<100; i++)
      h_array[i] = dist::GaussianMixture<3>(models);

    thrust::device_vector<dist::GaussianMixture<3> > d_array = h_array;

    thrust::device_vector<Real > d_out(h_array.size());

    randomSampleField(d_array.begin(), d_array.end(), d_out.begin());

    thrust::host_vector<Real > h_out = d_out;

    for (i=0; i<h_out.size(); i++)
      cout << "GMM sampling output: " << h_out[i] << endl;
  }

  return 0;
}
