// Copyright 2015 The Edda Authors. All rights reserved.
// Use of this source code is governed by a MIT-style license that can be
// found in the LICENSE file.

#ifndef EDDA_EXPORT_H_
#define EDDA_EXPORT_H_

#include "edda.h"

// Defines EDDA_EXPORT so that functionality implemented by
// the edda library can be exported to consumers.

#if defined(OS_WIN)
#define EDDA_EXPORT __declspec(dllexport)

#else  // defined(OS_WIN)
#define EDDA_EXPORT __attribute__((visibility("default")))
#endif

#endif  // EDDA_EXPORT_H_
