#include <cstdlib>
#include <ctime>
#include <iostream>
#include <iostream>

#include "io/edda_vtk_reader.h"
#include "io/edda_vtk_writer.h"
#include "io/edda_reader.h"
#include "io/edda_writer.h"
#include "distributions/histogram.h"
#include "dataset/distr_array.h"

using namespace edda;
using namespace std;
using namespace edda::dist;


int main()
{
	shared_ary<Histogram> array(new Histogram[8], 8);
	for (int i = 0; i < 8; i++){
		float data[1];
		data[0] = 50 + i;
		//array[i] = Histogram(data, 1, 50, 60, 20);
		//array[i] = eddaComputeHistogram(data, 1, 50, 60, 20);
		array[i] = eddaComputeHistogram(data, 1, 20);
	}
	DistrArray * abstract_array = new ScalarDistrArray<Histogram>(array);

	for (int i = 0; i<8; i++)
	{
		printf("%d\n", i);
		cout << i << ": " << abstract_array->getScalar(i) <<
			": sample = " << getSample(abstract_array->getScalar(i)) << endl;
	}

	shared_ptr<Dataset<Real> > dataset = make_Dataset<Real>(
		new RegularCartesianGrid(2, 2, 2),
		abstract_array
		);
	writeEddaVtkDataset(dataset, "testHist.vti", "test_");
	cout << "finish histogram modeling and write into file, press any key to continue!" << endl;
	//getchar();

	shared_ptr<Dataset<Real> > dataset2 = loadEddaScalarDataset("testHist.vti", "test_");
	cout << dataset2->getArray()->getDistrName() << endl;
	int* dim;
	dim = dataset2->getDimension();
	cout << dim[0] << dim[1] << dim[2] << endl;

	for (int i = 0; i < dim[0]; i++)
	for (int j = 0; j < dim[1]; j++)
	for (int k = 0; k < dim[2]; k++){
		cout << "at_comp(" << i << "," << j << "," << k << ") : " << dataset2->at_comp(i, j, k) << endl;
		cout << "at_comp_distr : " << dataset2->at_comp_distr(i, j, k) << endl;
	}
	cout << "load histogram from file and sampling, press any key to finish!" << endl;
	//getchar();


	{
		cout << endl << endl;

		writeEddaDataset(dataset, "testHist.edda");
		shared_ptr<Dataset<Real> > dataset3 = loadEddaScalarDataset_noneVTK("testHist.edda");
		cout << dataset3->getArray()->getDistrName() << endl;
		int* dim;
		dim = dataset3->getDimension();
		cout << dim[0] << dim[1] << dim[2] << endl;

		for (int i = 0; i < dim[0]; i++)
		for (int j = 0; j < dim[1]; j++)
		for (int k = 0; k < dim[2]; k++){
			cout << "at_comp(" << i << "," << j << "," << k << ") : " << dataset3->at_comp(i, j, k) << endl;
			cout << "at_comp_distr : " << dataset3->at_comp_distr(i, j, k) << endl;
		}
		cout << "load histogram from file and sampling, press any key to finish!" << endl;


		int* dims = dataset3->getDimension();
		double dif = 0;
		for (int k = 0; k < dims[2]; k++){
			for (int j = 0; j < dims[1]; j++){
				for (int i = 0; i < dims[0]; i++){
					dist::Variant distr = dataset->at_comp_distr(i, j, k);
					dist::Variant distr2 = dataset2->at_comp_distr(i, j, k);

					dist::Histogram curHist1 = boost::get<dist::Histogram>(distr);
					dist::Histogram curHist2 = boost::get<dist::Histogram>(distr2);
					
					dif = dif + abs(curHist1.getBins() - curHist2.getBins())
						+ abs(curHist1.getMinValue() - curHist2.getMinValue())
						+ abs(curHist1.getMaxValue() - curHist2.getMaxValue());

					for (int b = 0; b < curHist1.getBins(); b++)
					{
						dif = dif + abs(curHist1.getBinValue(b) - curHist2.getBinValue(b));
					}
				}
			}
		}
		cout << "the differences between the old vtk format and the new format is: " << dif << endl;

	}
	return 0;
}