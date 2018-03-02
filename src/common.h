// Copyright 2015 The Edda Authors. All rights reserved.
// Use of this source code is governed by a MIT-style license that can be
// found in the LICENSE file.

#ifndef COMMON_H_
#define COMMON_H_

#include <cassert>
#include <stdexcept>
#include <limits>

#include <cstdint>
#include <boost/config.hpp>
#include <type_traits>

#include "edda.h"

namespace edda {

const double DEG_TO_RAD = 0.0174532925199432957692;  // PI / 180
const double RAD_TO_DEG = 57.2957795130823208768;  // 180 / PI
const double PI_BY_2    = 1.57079632679489661923;  // PI / 2

const double EPS        = 1.0E-6;

// The distribution precision is set here.
// Using float instead of double for smaller storage and because the nature of uncertainty has less precision.
typedef float Real;

// numeric limits:
//template <typename T> using limits = std::numeric_limits<T>;

// You can add more for needed return status
enum ReturnStatus { SUCCESS = 0, FAIL, OUT_OF_BOUND };


#ifdef OS_WIN
typedef long long int64_t;
#endif

// May not be needed in the release version but just keep here for now.
struct NotImplementedException: public std::runtime_error {
  NotImplementedException() : std::runtime_error("Sorry, method not implemented.") {}
};

class OutOfBoundException{};

struct NotSupportException : public std::runtime_error {
	NotSupportException() : std::runtime_error("Sorry, method not supported.") {}
};

///
/// \brief This is useful to constrain what types are applicable for generic functions or classes.
/// @param T The class to be tested
/// @param B The supposed base class
///
#define ENABLE_IF_BASE_OF(T, B)  typename std::enable_if<std::is_base_of<B, T>::value>::type* = nullptr

#define MAX_GMs 5

}  // namespace edda

#endif  // COMMON_H_
