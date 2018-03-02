#include "test_common.h"
#include <core/ndarray.h>

using namespace edda;

int main ()
{
  NdArray<Real> p({1});
  p.set_val({0}, 3);
  {
    printf("Assigning\n");
    NdArray<Real> q = p;
    //TEST(q.get_val({0})==3);
  }
  {
    printf("with vector\n");
    std::vector<NdArray<Real> > vec(10, p);
    vec[0] = NdArray<Real> ({2,2});
  }
  printf("Test ends\n");
  return 0;

}
