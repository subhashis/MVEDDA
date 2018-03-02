#include <cstdlib>
#include <ctime>
#include <iostream>
#include <iostream>

#include "io/edda_vtk_reader.h"
#include "io/edda_vtk_writer.h"
using namespace edda;
using namespace std;


int main()
{
  shared_ptr<Dataset<Real> > dataset1 = loadEddaScalarDataset(SAMPLE_DATA_PATH "isabel_pressure_small.vti");
  writeEddaVtkDataset(dataset1, "test.vti", "test_");

  shared_ptr<Dataset<Real> > dataset2 = loadEddaScalarDataset("test.vti", "test_");
  cout << dataset2->getArray()->getDistrName() << endl;

  return 0;
}
