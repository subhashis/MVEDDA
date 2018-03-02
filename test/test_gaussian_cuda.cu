
#include <iostream>

#include "distributions/distribution.h"
#include "distributions/gaussian.h"
#include "core/interpolator.h"
#include "test_common.h"

using namespace std;
using namespace edda;

int main(int argc, char* argv[]) {
  dist::Gaussian g1(1.,1.);
  dist::Gaussian g2(2.,2.);
  dist::Gaussian g3;

  cout << "g1= "<< g1 << endl;
  cout << "g1= "<< g2 << endl;

  // adding two random variables
  g3=g1+g2;

  cout << "g1+g2= " << g3 << endl;
  TEST( g1.mean == dist::getMean(g1) );
  TEST( approx_equal(dist::getMean(g3), 3. ) );
  TEST( approx_equal(dist::getVar(g3), 3. ) );

  // linear interpolation
  g3 =  edda::lerp(g1, g2, .1);

  cout << "lerp(g1, g2, .1) = " << g3 << endl;
  TEST( approx_equal(dist::getMean(g3), 1.1));
  TEST( approx_equal(dist::getVar(g3), 2.*0.1*0.1+1.*0.9*0.9));

  cout << "CDF of g1 at 1: " << getCdf(g1, 1) << endl;
  TEST( approx_equal(getCdf(g1, 1), 0.5));

  // 1.96 std from the mean (one side) will cover 97.5% of the distribution
  double cdf = dist::getCdf(g1, getMean(g1)+1.96);
  cout << "cdf(g1,2.96)=" << cdf << endl;
  TEST( approx_equal(cdf, 0.975002) );

  cout << "A random sample of g3: " << getSample(g3) << endl;

  cout << "size of Gaussian<double>: " << sizeof(g1) << endl;
  TEST( sizeof(g1) == sizeof(Real)*2 );


  return 0;
}
