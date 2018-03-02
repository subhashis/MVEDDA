#ifndef TEST_COMMON_H
#define TEST_COMMON_H

#include <iostream>

#include "common.h"

#define TEST(b) { if (!(b)) { std::cerr << "TEST FAIL: " << #b << std::endl;  return 1; } }

template<typename T>
bool approx_equal(T a, T b) { return std::abs(a-b) < edda::EPS ; }


#endif // TEST_COMMON_H
