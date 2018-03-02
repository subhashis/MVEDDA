// Copyright 2015 The Edda Authors. All rights reserved.
// Use of this source code is governed by a MIT-style license that can be
// found in the LICENSE file.
#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <list>
#include <iterator>

#include "vtkDataSet.h"
#include "vtkSmartPointer.h"
#include "vtkLineSource.h"
#include "vtkRenderWindowInteractor.h"
#include "vtkInteractorStyleTrackballCamera.h"
#include "vtkRenderWindow.h"
#include "vtkRenderer.h"
#include "vtkActor.h"
#include "vtkPolyDataMapper.h"
#include "vtkOutlineFilter.h"
#include "vtkProperty.h"
#include "vtkLineWidget.h"
#include "vtkCommand.h"
#include "vtkCallbackCommand.h"
#include "vtkPointData.h"
#include "vtkFloatArray.h"

#include "vtkCubeSource.h"
#include "vtkIdList.h"
#include "vtkPoints.h"
#include "vtkCellArray.h"

#include "edda.h"
#include "distributions/gaussian.h"
#include "distributions/distribution.h"
#include "io/edda_vtk_reader.h"
#include "filters/stream_tracer.h"
#include "dataset/dataset.h"

#include "vtkTransformFilter.h"
#include "vtkTransform.h"

#include "io/edda_reader.h"
#include "io/edda_writer.h"


int SAMPLES = 100;

using namespace std;
using namespace edda;

#define vsp_new(cls,x) vtkSmartPointer<cls> x = vtkSmartPointer<cls>::New()

vtkSmartPointer<vtkLineWidget> lineWidget;
vtkSmartPointer<vtkRenderWindow> renWin;
vtkSmartPointer<vtkPolyData> streamlines;

shared_ptr<StreamTracer<VECTOR3> > streamTracer;

void computeStreamlines(vtkSmartPointer<vtkPolyData> vtk_seeds)
{
  int i;
  // convert seeds from vtk
  list<VECTOR3> seeds;
  for (i=0; i<vtk_seeds->GetNumberOfPoints(); i++)
  {
    double p[3];
    vtk_seeds->GetPoint(i, p);
    seeds.push_back(VECTOR3(p[0],p[1],p[2]));
  }

  list<list<VECTOR3> >traces;
  streamTracer->step_size = 0.05;
  streamTracer->max_steps = 1000;
  streamTracer->compute(seeds, traces);


  // convert traces to vtk
  vsp_new(vtkPoints, points);
  vsp_new(vtkCellArray, cells);
  vsp_new(vtkFloatArray, fieldData);
  fieldData->SetNumberOfComponents(1);

  int count = 0;
  for (auto traceItr = traces.begin(); traceItr != traces.end(); ++traceItr)
  {
    auto &singleTrace = *traceItr;
    vtkIdList *ids = vtkIdList::New();
    for (auto posItr = singleTrace.begin(); posItr!=singleTrace.end(); ++posItr)
    {
      VECTOR3 &p = *posItr;
      int id = points->InsertNextPoint(p[0], p[1], p[2]);
      ids->InsertNextId(id);

      fieldData->InsertTuple1(id, count);
    }
    cells->InsertNextCell(ids);
    count ++;
  }

  streamlines->SetPoints(points);
  streamlines->SetLines(cells);
  int array_idx = streamlines->GetPointData()->AddArray(fieldData);
  streamlines->GetPointData()->SetActiveAttribute(array_idx,  vtkDataSetAttributes::SCALARS);

}

void updateSeeds(vtkObject* caller, unsigned long eventId, void *clientdata, void *calldata)
{
  vsp_new(vtkPolyData, seeds) ;
  lineWidget->GetPolyData(seeds);
  computeStreamlines(seeds);
  renWin->Render();
}


int main(int argc, char **argv)
{
  cout << "Input arguments: <info file>" << endl;
  //if (argc<2)
  //    return -1;
  string filename;

  if (argc>1) {
	  filename = string(argv[1]);
  }
  else {
	  cout << "Loading sample file" << endl;
	  filename = string(SAMPLE_DATA_PATH) + "/SW6_PIV.vti";
  }

  cout << "Press 'i' to change the rake\n";

  // load data with random sampling
  shared_ptr<Dataset<VECTOR3> > datasetvtk = edda::loadEddaVector3Dataset(filename, "");

  //writeEddaDataset(datasetvtk, "testData_vector3.edda");

  //load edda dataset with our reader
  shared_ptr<Dataset<VECTOR3> > dataset = loadEddaVector3Dataset_noneVTK("testData_vector3.edda");

  {
	  //one by one check only for GMM5 and vector3 !!
	  if (dataset->getArray()->getDistrName() == "GaussianMixture5" && dataset->getArray()->getNumComponents()==3){
		  int* dims = dataset->getDimension();
		  double dif = 0;
		  for (int k = 0; k < dims[2]; k++){
			  for (int j = 0; j < dims[1]; j++){
				  for (int i = 0; i < dims[0]; i++){
					  std::vector<dist::Variant> distrVec = datasetvtk->at_comp_distr_vector(i, j, k);
					  std::vector<dist::Variant> distrVec2 = dataset->at_comp_distr_vector(i, j, k);
					  
					  for (int c = 0; c < 3; c++){
						  dist::GaussianMixture<5> curDist1 = boost::get<dist::GaussianMixture<5> >(distrVec[c]);
						  dist::GaussianMixture<5> curDist2 = boost::get<dist::GaussianMixture<5> >(distrVec2[c]);
						  for (int model = 0; model < 5; model++){
							  dif = dif + abs(curDist1.models[model].m - curDist2.models[model].m)
								  + abs(curDist1.models[model].v - curDist2.models[model].v)
								  + abs(curDist1.models[model].w - curDist2.models[model].w);
						  }
					  }
				  }
			  }
		  }
		  cout << "the differences between the old vtk format and the new format is: " << dif << endl;
	  }
	  else{
		  cout << "one by one check is not performed " << endl;
	  }
  }



  //shared_ptr<Dataset<VECTOR3> > dataset;
  VECTOR3 minB, maxB;
  dataset->getGrid()->boundary(minB, maxB);

  // create stream tracer
  streamTracer = std::make_shared<StreamTracer<VECTOR3> >(dataset);

  // create a dummy vtk bounding box
  vsp_new(vtkCubeSource, cube);
  cube->SetBounds(minB[0], maxB[0], minB[1], maxB[1], minB[2], maxB[2]);
  cube->Update();

  // user can change the seeding rake
  lineWidget = vtkSmartPointer<vtkLineWidget>::New();
  lineWidget->SetInputData(cube->GetOutput());
  lineWidget->SetResolution(19); // 20 seeds along the line
  lineWidget->SetAlignToYAxis();
  lineWidget->PlaceWidget();
  lineWidget->ClampToBoundsOn();

  vsp_new(vtkPolyData, seeds);
  lineWidget->GetPolyData(seeds);

  streamlines = vtkSmartPointer<vtkPolyData>::New();

  //computeStreamlines(seeds);

  vsp_new(vtkPolyDataMapper, mapper);
  mapper->SetScalarRange(0, 20); // for 20 seeds
  mapper->SetInputData(streamlines);
  mapper->SetColorModeToMapScalars();
  vsp_new(vtkActor, actor);
  actor->SetMapper(mapper);

  //
  // outline
  //
  vsp_new(vtkOutlineFilter, outline);
  outline->SetInputData(cube->GetOutput());

  vtkPolyDataMapper *outlineMapper = vtkPolyDataMapper::New();
  outlineMapper->SetInputConnection(outline->GetOutputPort());

  vtkActor *outlineActor = vtkActor::New();
  outlineActor->SetMapper(outlineMapper);
  outlineActor->GetProperty()->SetColor(.5,.5,.5);

  //
  // renderer
  //
  vsp_new(vtkRenderer, ren);
  renWin = vtkSmartPointer<vtkRenderWindow>::New();
  renWin->AddRenderer(ren);
  vsp_new(vtkRenderWindowInteractor, iren);
  iren->SetRenderWindow(renWin);
  vsp_new(vtkInteractorStyleTrackballCamera, style);
  iren->SetInteractorStyle(style);

  // line widget interactor
  lineWidget->SetInteractor(iren);
  lineWidget->SetDefaultRenderer(ren);
  vsp_new(vtkCallbackCommand, callback );
  callback->SetCallback(updateSeeds);
  lineWidget->AddObserver(vtkCommand::EndInteractionEvent, callback);


  ren->AddActor(actor);
  ren->AddActor(outlineActor);
  ren->SetBackground(0,0,0);
  ren->ResetCamera();

  renWin->SetSize(500,500);

  iren->Initialize();

  // Sign up to receive TimerEvent
  iren->AddObserver(vtkCommand::TimerEvent, callback);
  int timerId = iren->CreateRepeatingTimer(100);

  renWin->Render();
  iren->Start();

  return 0;
}



