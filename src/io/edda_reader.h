#ifndef EDDA_READER
#define EDDA_READER

#include <string>
#include "edda_export.h"
#include "dataset/dataset.h"

namespace edda{
	std::shared_ptr<Dataset<Real> > EDDA_EXPORT loadEddaScalarDataset_noneVTK(const std::string &edda_file);
}

#endif
