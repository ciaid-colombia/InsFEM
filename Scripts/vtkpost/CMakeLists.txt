cmake_minimum_required(VERSION 3.3 FATAL_ERROR)

project(vtkpost)

find_package(VTK COMPONENTS 
  vtkCommonCore
  vtkCommonDataModel
  vtkFiltersGeometry
  vtkIOXML QUIET)
if (NOT VTK_FOUND)
  message("Skipping vtkpost: ${VTK_NOT_FOUND_MESSAGE}")
  return ()
endif()
message (STATUS "VTK_VERSION: ${VTK_VERSION}")

find_package(MPI REQUIRED)
include_directories(SYSTEM ${MPI_C_INCLUDE_PATH})

if (VTK_VERSION VERSION_LESS "8.90.0")
  # old system
  include(${VTK_USE_FILE})
  add_library(vtkpost OBJECT vtkc.cpp)
  #add_executable(WriteVTI MACOSX_BUNDLE WriteVTI.cxx )
  target_link_libraries(vtkpost PRIVATE ${VTK_LIBRARIES})
else ()
  # include all components
  #add_executable(WriteVTI MACOSX_BUNDLE WriteVTI.cxx )
  add_library(vtkpost OBJECT vtkc.cpp)
  target_link_libraries(vtkpost PRIVATE ${VTK_LIBRARIES})
  # vtk_module_autoinit is needed
  vtk_module_autoinit(
    TARGETS vtkpost
    MODULES ${VTK_LIBRARIES}
    )
endif () 
