# Copyright 2015 The Edda Authors. All rights reserved.
# Use of this source code is governed by a MIT-style license that can be
# found in the LICENSE file.

cmake_minimum_required (VERSION 2.8)

# Set the project's name.
project (edda)

### Provide options that the user can optionally select.
option (EDDA_BUILD_VIS "Build distribution vis tools." OFF)

option (EDDA_BUILD_TESTS "Build the test programs." OFF)
option (EDDA_BUILD_EXAMPLES "Build the example programs." ON)
option (EDDA_BUILD_MV "Build the example programs." ON)
option (EDDA_BUILD_SHARED_LIBS "Build into shared libraries." OFF)
#option (EDDA_WITH_OPENMP "Build with OpenMP parallelization." ON)

set (THRUST_DEVICE_SYSTEM "CPP" CACHE STRING "Thrust device backend (CUDA/TBB/OPENMP/CPP).")
set_property(CACHE THRUST_DEVICE_SYSTEM PROPERTY STRINGS "CUDA;TBB;OPENMP;CPP")
###

# VTK integration
find_package(VTK)
include (${VTK_USE_FILE})
list(APPEND LINK_LIBS ${VTK_LIBRARIES})

# Find vinecopulib package and dependencies
find_package(vinecopulib                  REQUIRED)
find_package(Boost 1.56                   REQUIRED)
include(cmake/findEigen3.cmake            REQUIRED)
find_package(Threads                      REQUIRED)
find_package(OpenCV     REQUIRED)


# Set required variables for includes and libraries
# In the second line
#   * VINECOPULIB_LIBRARIES is needed if vinecopulib has been built as a
#     shared lib (does nothing otherwise).
#   * CMAKE_THREAD_LIBS_INIT is needed for some linux systems
#     (but does nothing on OSX/Windows).
set(external_includes ${VINECOPULIB_INCLUDE_DIR} ${EIGEN3_INCLUDE_DIR} ${Boost_INCLUDE_DIRS} ${OpenCV_INCLUDE_DIRS})
set(external_libs ${VINECOPULIB_LIBRARIES} ${CMAKE_THREAD_LIBS_INIT} ${OpenCV_LIBS})



# Tell the compiler to use C++11, if not using Visual Studio.
if (NOT MSVC)
    include(CheckCXXCompilerFlag)
    CHECK_CXX_COMPILER_FLAG("-std=c++11" COMPILER_SUPPORTS_CXX11)
    CHECK_CXX_COMPILER_FLAG("-std=c++0x" COMPILER_SUPPORTS_CXX0X)
    if(COMPILER_SUPPORTS_CXX11)
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
    elseif(COMPILER_SUPPORTS_CXX0X)
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++0x")
    else()
        message(WARNING "The compiler ${CMAKE_CXX_COMPILER} has no C++11 support. Please use a different C++ compiler.")
    endif()
endif()


# All warnings on
#set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}  -Wall)

set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}  -fext-numeric-literals")


# Load settings for the boost project.
find_package (Boost REQUIRED)
if (Boost_FOUND)
  include_directories (${Boost_INCLUDE_DIRS})
  link_libraries (${Boost_LIBRARIES})
  
else ()
set ( BOOST_ROOT "/usr/local" CACHE PATH "the root directory containing Boost")
endif ()
set ( BOOST_ROOT "/usr/local" CACHE PATH "the root directory containing Boost")

# default CUDA off
set(EDDA_BUILD_WITH_CUDA OFF)

# To fix compilation problem: relocation R_X86_64_32 against `a local symbol' can not be
# used when making a shared object; recompile with -fPIC
# See http://www.cmake.org/pipermail/cmake/2007-May/014350.html
IF( CMAKE_SYSTEM_PROCESSOR STREQUAL "x86_64" )
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fPIC")
ENDIF()

## Thrust with CUDA backend
if (${THRUST_DEVICE_SYSTEM} STREQUAL "CUDA")
  find_package(CUDA REQUIRED)

  #list(APPEND CUDA_NVCC_FLAGS "-ffast-math")
  #list(APPEND CUDA_NVCC_FLAGS --compiler-options -fno-strict-aliasing -lineinfo -use_fast_math -Xptxas -dlcm=cg)
  #list(APPEND CUDA_NVCC_FLAGS  "--std=c++11") #cuda 7 feature
  #list(APPEND CUDA_NVCC_FLAGS -gencode arch=compute_20,code=sm_20)
  #list(APPEND CUDA_NVCC_FLAGS -gencode arch=compute_30,code=sm_30)
  #list(APPEND CUDA_NVCC_FLAGS -gencode arch=compute_35,code=sm_35)

  # To successfully build a library, this is required (rdc=true)
  # See http://stackoverflow.com/questions/13683575/cuda-5-0-separate-compilation-of-library-with-cmake
  set(CUDA_ATTACH_VS_BUILD_RULE_TO_CUDA_FILE ON)
  set(CUDA_SEPARABLE_COMPILATION ON)

  list(APPEND LINK_LIBS ${CUDA_LIBRARIES} ${CUDA_CUDART_LIBRARY})

  set(EDDA_BUILD_WITH_CUDA ON)

## Thrust with TBB backend
elseif (${THRUST_DEVICE_SYSTEM} STREQUAL "TBB")

  add_definitions(-DTHRUST_DEVICE_SYSTEM=THRUST_DEVICE_SYSTEM_TBB)

  set ( THRUST_INCLUDE_DIR "/usr/local/include" CACHE PATH "Thrust root dir")
  include_directories ( ${THRUST_INCLUDE_DIR} )

  set ( TBB_INSTALL_PATH "/usr/local" CACHE PATH "TBB root dir")
  include_directories(${TBB_INSTALL_PATH}/include)
  link_directories( ${TBB_INSTALL_PATH}/lib)
  list(APPEND LINK_LIBS tbb)

## Thrust with OpenMP backend
elseif (${THRUST_DEVICE_SYSTEM} STREQUAL "OPENMP")
  add_definitions ( -DTHRUST_DEVICE_SYSTEM=THRUST_DEVICE_SYSTEM_OMP)

  set ( THRUST_INCLUDE_DIR "/usr/local/include" CACHE PATH "Thrust root dir")
  include_directories ( ${THRUST_INCLUDE_DIR} )

  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fopenmp")

elseif (${THRUST_DEVICE_SYSTEM} STREQUAL "CPP")
  add_definitions ( -DTHRUST_DEVICE_SYSTEM=THRUST_DEVICE_SYSTEM_CPP)

  set ( THRUST_INCLUDE_DIR "/usr/local/include" CACHE PATH "Thrust root dir")
  include_directories ( ${THRUST_INCLUDE_DIR} )

else ()
  message(ERROR "Invalid setting of THRUST_DEVICE_SYSTEM")
endif ()


# Set libraries' output path to lib.
set (LIBRARY_OUTPUT_PATH "${PROJECT_BINARY_DIR}/lib")

# Set variables for platform detection.
if (${CMAKE_SYSTEM_NAME} MATCHES "Linux")
  set (OS_LINUX ON)
endif ()

if (${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
  set (OS_MACOSX ON)
endif ()

if (WIN32)
  set (OS_WIN ON)
endif ()

if (${EDDA_BUILD_SHARED_LIBS})
    set(LIB_TYPE SHARED)
else ()
    set(LIB_TYPE STATIC)
endif ()

# for examples
add_definitions(-DSAMPLE_DATA_PATH="${CMAKE_SOURCE_DIR}/sample_data/")

##############################################
# sources

# Make the configuration file.
configure_file ("${PROJECT_SOURCE_DIR}/edda.h.in"
                "${PROJECT_SOURCE_DIR}/src/edda.h")

# Add include directories to the build.
include_directories ("${PROJECT_SOURCE_DIR}/src")

# Add the src subdirectory to the build.
add_subdirectory (src)

##############################################
# tools

if (EDDA_BUILD_VIS)
  add_subdirectory (vis)
endif ()

# If EDDA_BUILD_TESTS is on, enable testing.
if (EDDA_BUILD_TESTS)
  enable_testing ()
  add_subdirectory (test)
endif ()

# If EDDA_BUILD_EXAMPLES is on, add the examples subdirectory to the build.
if (EDDA_BUILD_EXAMPLES)
  add_subdirectory (examples)
endif ()

# If EDDA_BUILD_MV is on, add the examples subdirectory to the build.
if (EDDA_BUILD_MV)
  #add_subdirectory (multivariate)
  add_subdirectory (codda)
endif ()
