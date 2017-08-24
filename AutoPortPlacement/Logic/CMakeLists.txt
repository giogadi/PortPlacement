cmake_minimum_required(VERSION 3.5)

if(NOT DEFINED CMAKE_MACOSX_RPATH)
  set(CMAKE_MACOSX_RPATH OFF)
endif()

#-----------------------------------------------------------------------------
project(vtkSlicer${MODULE_NAME}ModuleLogic)

set(CMAKE_MODULE_PATH
  ${CMAKE_CURRENT_SOURCE_DIR}/CMake
  ${CMAKE_MODULE_PATH})

find_package(Eigen3 REQUIRED CONFIG)

find_package(NLopt REQUIRED)

set(KIT ${PROJECT_NAME})

set(${KIT}_EXPORT_DIRECTIVE "VTK_SLICER_${MODULE_NAME_UPPER}_MODULE_LOGIC_EXPORT")

set(${KIT}_INCLUDE_DIRECTORIES
  ${NLopt_INCLUDE_DIRS}
  )
include_directories(${${KIT}_INCLUDE_DIRECTORIES})

set(${KIT}_SRCS
  vtkSlicer${MODULE_NAME}Logic.cxx
  vtkSlicer${MODULE_NAME}Logic.h
  )

set(${KIT}_TARGET_LIBRARIES
  collisions
  davinci
  ${ITK_LIBRARIES}
  Eigen3::Eigen
  ${NLopt_LIBRARIES}
  optim
  vtkSlicerAnnotationsModuleMRML
  vtkSlicerMarkupsModuleMRML
  ${VTK_LIBRARIES}
  )

include_directories(${CMAKE_CURRENT_SOURCE_DIR})
include_directories(${CMAKE_CURRENT_BINARY_DIR})

add_subdirectory(collisions)
add_subdirectory(davinci-kinematics)
add_subdirectory(optim)

if(BUILD_TESTING)
  add_subdirectory(tests)
endif()

#-----------------------------------------------------------------------------
SlicerMacroBuildModuleLogic(
  NAME ${KIT}
  EXPORT_DIRECTIVE ${${KIT}_EXPORT_DIRECTIVE}
  INCLUDE_DIRECTORIES ${${KIT}_INCLUDE_DIRECTORIES}
  SRCS ${${KIT}_SRCS}
  TARGET_LIBRARIES ${${KIT}_TARGET_LIBRARIES}
  )
