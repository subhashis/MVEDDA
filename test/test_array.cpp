#include "vtkSmartPointer.h"
#include "vtkXMLPolyDataReader.h"
#include "vtkXMLPolyDataWriter.h"
#include "vtkPoints.h"
#include "vtkPolyData.h"
#include "vtkDoubleArray.h"
#include "vtkPointData.h"
#include "vtkVariant.h"
#include "vtkVariantArray.h"

int main( int argc, char *argv[])
{
  if( argc < 3 )
  {
    std::cerr << "Missing Parameters " << std::endl;
    std::cerr << "Usage:   " << argv[0] << " inputfile outputfile" <<
std::endl;
    std::cerr << "Example: \"" << argv[0] << " input.vtp output.vtp\"" <<
std::endl;
    return EXIT_FAILURE;
  }

  vtkSmartPointer<vtkXMLPolyDataReader> inputReader =
vtkSmartPointer<vtkXMLPolyDataReader>::New();
  inputReader->SetFileName(argv[1]);
  try
  {
      inputReader->Update();
  }
  catch(...)
  {
      std::cerr << "Error occurs when reading " << argv[2] << std::endl;
      return EXIT_FAILURE;
  }

  vtkSmartPointer<vtkPolyData> poly = vtkSmartPointer<vtkPolyData>::New();
  poly->DeepCopy(inputReader->GetOutput());

  vtkSmartPointer<vtkVariantArray> vararray =
vtkSmartPointer<vtkVariantArray>::New();
  vararray->SetName("Variant");
  vararray->SetNumberOfValues(poly->GetNumberOfPoints());
  for(vtkIdType id=0; id<poly->GetNumberOfPoints(); id++)
  {
      vtkSmartPointer<vtkDoubleArray> darray =
vtkSmartPointer<vtkDoubleArray>::New();
      darray->SetName("Double");
      darray->SetNumberOfValues(id%5+1);
      for(vtkIdType j=0; j<darray->GetNumberOfTuples(); j++)
darray->SetValue(j, j*0.5);
      vtkVariant variant(darray);

      vararray->SetValue(id, variant);
  }
  poly->GetPointData()->AddArray(vararray);

  vtkSmartPointer<vtkXMLPolyDataWriter> outputWriter =
vtkSmartPointer<vtkXMLPolyDataWriter>::New();
  outputWriter->SetInput(poly);
  outputWriter->SetFileName(argv[2]);
  outputWriter->SetDataModeToBinary();
  try
  {
      outputWriter->Write();
  }
  catch(...)
  {
      std::cerr << "Error occurs when writing " << argv[2] << std::endl;
      return EXIT_FAILURE;
  }

  return EXIT_SUCCESS;
}
