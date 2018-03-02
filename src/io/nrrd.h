// Copyright 2015 The Edda Authors. All rights reserved.
// Use of this source code is governed by a MIT-style license that can be
// found in the LICENSE file.

#ifndef NRRD_H_
#define NRRD_H_

#include "common.h"
#include "edda_export.h"

namespace edda {

ReturnStatus EDDA_EXPORT write_nrrd_3d(const char *nrrd_fname, const char *raw_fname,
                           int w, int h, int d, const char *type);

}  // namespace edda

#endif  // NRRD_H_
