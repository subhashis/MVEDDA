#include <cstdlib>
#include <ctime>
#include <iostream>

#include "distributions/distribution.h"
#include "distributions/variant.h"
#include "dataset/distr_array.h"

using namespace edda;
using namespace std;
using namespace edda::dist;

DefaultGaussianMixture gen_gmm()
{
  DefaultGaussianMixture gmm;
  for (int i=0; i<MAX_GMs; i++)
  {
    gmm.models[i].m = i;  gmm.models[i].v = 1; gmm.models[i].w=1;
  }
  gmm.normalizeWeights();

  return gmm;
}

GaussianMixture<2> gen_gmm2()
{
  GaussianMixture<2> gmm;
  for (int i=0; i<2; i++)
  {
    gmm.models[i].m = i;  gmm.models[i].v = 1; gmm.models[i].w=1;
  }
  gmm.normalizeWeights();

  return gmm;
}


DistrArray *make_Gaussian_array() {
  shared_ary<Gaussian> array (new Gaussian[10], 10);
  DistrArray * abstract_array = new ScalarDistrArray<Gaussian>(array);
  return abstract_array;
}


DistrArray * make_GMM_array() {
  shared_ary<DefaultGaussianMixture> array (new DefaultGaussianMixture[10], 10);
  for (int i=0; i<10; i++)
    array[i] = gen_gmm();
  DistrArray * abstract_array = new ScalarDistrArray<DefaultGaussianMixture>(array);
  return abstract_array;
}

DistrArray * make_GMM2_array() {
  shared_ary<GaussianMixture<2> > array (new GaussianMixture<2>[10], 10);
  for (int i=0; i<10; i++)
    array[i] = gen_gmm2();
  DistrArray * abstract_array = new ScalarDistrArray<GaussianMixture<2> >(array);
  return abstract_array;
}

DistrArray * make_hybrid_array() {
  shared_ary<Variant> array (new Variant[10], 10);
  // add Gaussian
  array[0] = Gaussian(1,2);
  // assign GaussianMixture
  array[1] = gen_gmm();

  DistrArray * abstract_array = new ScalarDistrArray<Variant>(array);
  return abstract_array;
}

DistrArray * make_JointGaussian_array() {
  shared_ary<JointGaussian> array(new JointGaussian[10], 10);
  
  DistrArray * abstract_array = new JointDistrArray<JointGaussian>(array);
  return abstract_array;
}

DistrArray * make_JointHistogram_array() {
  shared_ary<JointHistogram> array(new JointHistogram[10], 10);
  DistrArray * abstract_array = new JointDistrArray<JointHistogram>(array);
  return abstract_array;
}

DistrArray * make_JointGMM_array() {
	shared_ary<JointGMM> array(new JointGMM[10], 10);
	DistrArray * abstract_array = new JointDistrArray<JointGMM>(array);
	return abstract_array;
}

DistrArray * make_NEWGMM_array() {
	shared_ary<GMM> array(new GMM[10], 10);
	DistrArray * abstract_array = new ScalarDistrArray<GMM>(array);
	return abstract_array;
}

int main()
{

  DistrArray * array1 = make_Gaussian_array();
  DistrArray * array2 = make_GMM_array();
  DistrArray * array22 = make_GMM2_array();
  DistrArray * array3 = make_hybrid_array();
  DistrArray * array4 = make_JointGaussian_array();
  DistrArray * array5 = make_JointHistogram_array();
  DistrArray * array6 = make_JointGMM_array();

  int i;
  cout << "array1: " << endl;
  for (i=0; i<10; i++)
  {
    cout << i << ": " << array1->getDistr(i) <<
            ": sample = " << array1->getScalar(i) << endl;
  }
  cout << endl ;

  cout << "array2: "<< endl;
  for (i=0; i<10; i++)
  {
    cout << i << ": " << array2->getDistr(i) <<
            ": sample = " << array2->getScalar(i) << endl;
  }
  cout << endl;

  cout << "array3: "<< endl;
  for (i=0; i<2; i++)
  {
    cout << i << ": " << array3->getDistr(i) <<
            ": sample = " << array3->getScalar(i) << endl;
  }
  cout << endl;

  cout << "array4: " << endl;
  for (i=0; i<10; i++)
  {
    vector<Real> v = array4->getVector(i);
    cout << i << ": " << array4->getDistr(i) <<
            ": sample = " << v[0] << " " << v[1] << " " << v[2] << endl;
  }

  cout << "array5: " << endl;
  for (i=0; i<10; i++)
  {
    vector<Real> v = array5->getVector(i);
    cout << i << ": " << array5->getDistr(i) <<
            ": sample = " << v[0] << " " << v[1] << " " << v[2] << endl;
  }

  cout << "array6: " << endl;
  for (i = 0; i<10; i++)
  {
	  vector<Real> v = array6->getVector(i);
	  cout << i << ": " << array6->getDistr(i) <<
		  ": sample = " << v[0] << " " << v[1] << " " << v[2] << endl;
  }
}
