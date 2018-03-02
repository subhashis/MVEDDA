#include <memory>
#include <vtkDataSet.h>
#include <vtkXMLImageDataReader.h>
#include <vtkXMLStructuredGridReader.h>
#include <vtkImageData.h>
#include <vtkNew.h>
#include <vtkFieldData.h>
#include <vtkStringArray.h>
#include <vtkStdString.h>
#include <vtkPoints.h>

#include "io/path.h"

#include "io/gmm_vtk_data_array.h"
#include "io/edda_vtk_reader.h"
#include "dataset/distr_array.h"
#include "distributions/histogram.h"

using namespace std;
using namespace edda;
using namespace dist;

namespace edda {


string getDistrType(vtkFieldData* vtk_data, const string &array_name_prefix)
{
  string array_name = array_name_prefix + "distr_type";
  
  vtkStringArray *vtk_array = vtkStringArray::SafeDownCast(vtk_data->GetAbstractArray(array_name.c_str()));
  if (vtk_array==NULL) {
    cout << "Distribution data type not defind in field " << array_name << ".  Use default GaussianMixture." << endl;
    stringstream ss;
    ss << "GaussianMixture" << MAX_GMs; // default distribution type
    return ss.str();
  }
  return vtk_array->GetValue(0);
}

DistrArray *genHistoArray(vtkPointData *vtk_point_data)
{
  char arrayName[1024];
  sprintf(arrayName, "Variable_0");
  vtkSmartPointer<vtkDataArray> array = vtk_point_data->GetArray(arrayName);
  int c = array->GetNumberOfComponents();
  int n = array->GetNumberOfTuples();

  double* histData = (double*)malloc(sizeof(double)*c);
  shared_ary<Histogram> histAry(new Histogram[n], n);
  for (int nn = 0; nn < n; nn++){
          array->GetTuple(nn, histData);
          histAry[nn] = Histogram(histData);
  }
  free(histData);
  return new ScalarDistrArray<Histogram>(histAry);
}

template <typename T>
shared_ptr<Dataset<T> > loadEddaDataset(const string &edda_file, const string &array_name_prefix)
{
  string ext = getFileExtension(edda_file);

  if (ext.compare("vti")==0) {
    vtkNew<vtkXMLImageDataReader> reader;
    reader->SetFileName(edda_file.c_str());
    reader->Update();
    vtkImageData *vtkdata = reader->GetOutput();

    int *dim = vtkdata->GetDimensions();

    // check dataset type
    string type = getDistrType(vtkdata->GetFieldData(), array_name_prefix);

    if (type.compare(0, 15, "GaussianMixture")==0) {
        shared_ptr<Dataset<T> > dataset = make_Dataset<T>(
                         new RegularCartesianGrid(dim[0], dim[1], dim[2]),
                         new GmmVtkDataArray( vtkdata->GetPointData(), array_name_prefix.c_str() )
       );
      return dataset;


    } else if (type.compare("Histogram")==0)
    {
      shared_ptr<Dataset<T> > dataset = make_Dataset<T>(
                          new RegularCartesianGrid(dim[0], dim[1], dim[2]),
                          genHistoArray(vtkdata->GetPointData())
                          );
      return dataset;
    }
    else {
      cout << "Unknown distribution type: " << type << endl;
      exit(1);
    }


  } // end of if (ext.compare("vti")==0)
  else if (ext.compare("vts")==0){ // structured grids

    vtkNew<vtkXMLStructuredGridReader> reader;
    reader->SetFileName(edda_file.c_str());
    reader->Update();
    vtkStructuredGrid *vtkdata = reader->GetOutput();

    int *dim = vtkdata->GetDimensions();

    float *point_ary = (float *)vtkdata->GetPoints()->GetVoidPointer(0);

    // check dataset type
    string type = getDistrType(vtkdata->GetFieldData(), array_name_prefix);

    if (type.compare(0, 15, "GaussianMixture")==0) {
      shared_ptr<Dataset<T> > dataset = make_Dataset<T>(
                         new CurvilinearGrid(dim, point_ary),
                         new GmmVtkDataArray( vtkdata->GetPointData(), array_name_prefix.c_str() )
       );
      return dataset;

    } else if (type.compare("Histogram")==0)
    {
      shared_ptr<Dataset<T> > dataset = make_Dataset<T>(
                          new CurvilinearGrid(dim, point_ary),
                          genHistoArray(vtkdata->GetPointData() )
                          );
      return dataset;
    }
    else {
      cout << "Unknown distribution type: " << type << endl;
      exit(1);
    }


  } else {
    printf("File format of %s not supported\n", edda_file.c_str());
    exit(1);
  }
}

shared_ptr<Dataset<Real> > loadEddaScalarDataset(const string &edda_file, const string &array_name_prefix)
{
  return loadEddaDataset<Real>(edda_file, array_name_prefix);
}

shared_ptr<Dataset<VECTOR3> > loadEddaVector3Dataset(const string &edda_file, const string &array_name_prefix)
{
  return loadEddaDataset<VECTOR3>(edda_file, array_name_prefix);
}


} // edda
