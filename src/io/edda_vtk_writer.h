#ifndef EDDA_VTK_WRITER
#define EDDA_VTK_WRITER

#include <string>
#include "edda_export.h"
#include "dataset/dataset.h"

namespace edda{

  void EDDA_EXPORT writeEddaVtkDataset(std::shared_ptr<Dataset<Real> > , const std::string &edda_file, const std::string &array_name_prefix = "");

}

#endif
