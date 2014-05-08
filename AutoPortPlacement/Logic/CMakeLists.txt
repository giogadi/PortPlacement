project(vtkSlicer${MODULE_NAME}ModuleLogic)

set(CMAKE_MODULE_PATH
  ${CMAKE_CURRENT_SOURCE_DIR}/CMake
  ${CMAKE_MODULE_PATH})

# we need Eigen for parts of this logic
find_package(Eigen3 REQUIRED)
include_directories(${EIGEN3_INCLUDE_DIR})

set(nlopt_include_directory "" CACHE PATH "directory where the nlopt.h file is located")
  if(nlopt_include_directory STREQUAL "")
    message( FATAL_ERROR "please set the nlopt_include_directory variable value.")
  else()
    include_directories( ${nlopt_include_directory} )
  endif()

set(nlopt_library_directory "" CACHE PATH "directory where the nlopt library is located")
  if(nlopt_library_directory STREQUAL "")
    message( FATAL_ERROR "please set the nlopt_library_directory variable value.")
  else()
    link_directories( ${nlopt_library_directory} )
  endif()

set(KIT ${PROJECT_NAME})

set(${KIT}_EXPORT_DIRECTIVE "VTK_SLICER_${MODULE_NAME_UPPER}_MODULE_LOGIC_EXPORT")

set(${KIT}_INCLUDE_DIRECTORIES
  )

set(${KIT}_SRCS
  vtkSlicer${MODULE_NAME}Logic.cxx
  vtkSlicer${MODULE_NAME}Logic.h
  )

set(${KIT}_TARGET_LIBRARIES
  ${ITK_LIBRARIES}
  davinci
  optim
  vtkSlicerMarkupsModuleMRML
  vtkSlicerAnnotationsModuleMRML
  )

include_directories(${CMAKE_CURRENT_SOURCE_DIR})
include_directories(${CMAKE_CURRENT_BINARY_DIR})

add_subdirectory(collisions)
add_subdirectory(davinci-kinematics)
add_subdirectory(tests)
add_subdirectory(optim)

#-----------------------------------------------------------------------------
SlicerMacroBuildModuleLogic(
  NAME ${KIT}
  EXPORT_DIRECTIVE ${${KIT}_EXPORT_DIRECTIVE}
  INCLUDE_DIRECTORIES ${${KIT}_INCLUDE_DIRECTORIES}
  SRCS ${${KIT}_SRCS}
  TARGET_LIBRARIES ${${KIT}_TARGET_LIBRARIES}
  )
