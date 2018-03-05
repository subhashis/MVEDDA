# Copyright 2015 The Edda Authors. All rights reserved.
# Use of this source code is governed by a MIT-style license that can be
# found in the LICENSE file.

project(codda)

set(sources
#  MV_hist_regular.cpp
#  createDensityField.cpp
#  createConditionalField.cpp
#  regular_partitioning.cpp
#   regular_partition_loader.cpp
#   regular_partition_loader1.cpp
#   regular_partition_recon.cpp
#  regular_partition_vine.cpp
#  biCop_main.cpp
#   vineCop_main.cpp
#   TOY_2D_test1.cpp
#   TOY_2D_test2.cpp
#    turbine_create_distribution.cpp
#    turbine_create_parameters.cpp
#    turbine_conditional_field.cpp
#    turbine_conditional_field_PWG.cpp
#    simple_GMM_test.cpp
#    isabel_distribution_analysis.cpp
#    isabel_recon.cpp
#    isabel_store_grid_map.cpp
#    isabel_arbitrary_recon.cpp
#    turbine_distribution_analysis.cpp
#    turbine_onePsg_store_grid_map.cpp
#    turbine_onePsg_store_vts_grid_map.cpp
#    turbine_onePsg_arbitrary_recon.cpp
#    turbine_onePsg_vts_arbitrary_recon.cpp
#    turbine_store_grid_map.cpp
#    turbine_arbitrary_recon.cpp
#    isabel_multi_res.cpp
#    rmse_compute.cpp
#    test_turbine_onePsg.cpp


#    turbine_perProcessor_onePsg_recon.cpp
#    turbine_perProcessor_fullAnnulus_recon.cpp
#    turbine_perProcessor_fullAnnulus_conditionalQuery.cpp
#    isabel_store_grid_map.cpp
#    isabel_arbitrary_recon.cpp
#    isabel_arbitrary_recon_backup.cpp
#    test_read.cpp
#	  combustion_store_grid_map.cpp
    combustion_arbitrary_recon.cpp
    new_isabel_store_grid_map.cpp
    isabel_test_grid.cpp
    new_combustion_store_grid_map.cpp
    combustion_test_grid.cpp
    

)


include_directories(${CMAKE_SOURCE_DIR}/src)
include_directories(${external_includes})

# Force compile .cu files
#add_definitions("-x c++")


foreach (file ${sources})

  get_filename_component(name ${file} NAME_WE)

  if (EDDA_BUILD_WITH_CUDA)
    cuda_add_executable ( ${name} ${file} )
  else ()
    SET_SOURCE_FILES_PROPERTIES(${file} PROPERTIES LANGUAGE CXX)   
    add_executable ( ${name} ${file} )
  endif ()

  target_link_libraries ( ${name} edda ${LINK_LIBS} ${VINECOPULIB_LIBRARIES} ${CMAKE_THREAD_LIBS_INIT} ${OpenCV_LIBS})
endforeach ()


