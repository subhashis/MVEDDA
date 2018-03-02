# Copyright 2015 The Edda Authors. All rights reserved.
# Use of this source code is governed by a MIT-style license that can be
# found in the LICENSE file.

option (EDDA_BUILD_PARAVIEW_PLUGINS "Build Paraview Plugins. If ON, EDDA_BUILD_VTK_FILTERS must also be ON." OFF)

set(VIS_LINK_LIBS ${LINK_LIBS})
list(APPEND VIS_LINK_LIBS ${PROJECT_NAME}) # edda

# Choose either Paraview or VTK library
if (EDDA_BUILD_PARAVIEW_PLUGINS)
  find_package(ParaView REQUIRED)
  include(${PARAVIEW_USE_FILE})
  list(APPEND VIS_LINK_LIBS ${VTK_LIBRARIES})

else ()
  # VTK integration
  find_package(VTK)
  include (${VTK_USE_FILE})
  list(APPEND VIS_LINK_LIBS ${VTK_LIBRARIES})

endif ()


#################### VTK fix for Cuda #######################
if (EDDA_BUILD_WITH_CUDA)
    macro(REMOVE_VTK_DEFINITIONS)
        get_directory_property(_dir_defs DIRECTORY ${CMAKE_SOURCE_DIR} COMPILE_DEFINITIONS)
        set(_vtk_definitions)
        foreach(_item ${_dir_defs})
            if(_item MATCHES "vtk*")
                list(APPEND _vtk_definitions -D${_item})
            endif()
        endforeach()
        remove_definitions(${_vtk_definitions})
    endmacro(REMOVE_VTK_DEFINITIONS)

    REMOVE_VTK_DEFINITIONS()
endif()


# Add include directories to the build.
include_directories (${CMAKE_SOURCE_DIR}/src)
include_directories (${CMAKE_SOURCE_DIR}/vis/src)


# Add the src subdirectory to the build.
add_subdirectory (src)

add_subdirectory(examples)

