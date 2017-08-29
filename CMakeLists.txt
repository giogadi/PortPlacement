cmake_minimum_required(VERSION 3.5)

project(PortPlacement)

#-----------------------------------------------------------------------------
set(EXTENSION_HOMEPAGE "http://www.slicer.org/slicerWiki/index.php/Documentation/Nightly/Extensions/PortPlacement")
set(EXTENSION_CATEGORY "IGT")
set(EXTENSION_CONTRIBUTORS "Andinet Enquobahrie (Kitware), Luis G. Torres (UNC), Jean-Christophe Fillion-Robin (Kitware)")
set(EXTENSION_DESCRIPTION "Assists in the planning of surgical port placement in a laparoscopic procedure.")
set(EXTENSION_ICONURL "http://www.slicer.org/slicerWiki/images/a/a7/Portplacement_icon.png")
set(EXTENSION_SCREENSHOTURLS "http://www.slicer.org/slicerWiki/images/3/32/Portplacement.png")
set(EXTENSION_DEPENDS "Eigen3")
set(EXTENSION_BUILD_SUBDIRECTORY inner-build)
set(EXTENSION_STATUS "This extension is intended for public use and we expect to fully maintain it.")

set(SUPERBUILD_TOPLEVEL_PROJECT inner)

#
# Implementation agnostic fPIC
#
set( CMAKE_POSITION_INDEPENDENT_CODE ON )

#-----------------------------------------------------------------------------
# Extension dependencies
find_package(Slicer REQUIRED)
include(${Slicer_USE_FILE})
mark_as_superbuild(Slicer_DIR)

find_package(Git REQUIRED)
mark_as_superbuild(GIT_EXECUTABLE)

find_package(Eigen3 REQUIRED CONFIG)
mark_as_superbuild(Eigen3_DIR)

#-----------------------------------------------------------------------------
# SuperBuild setup
option(${EXTENSION_NAME}_SUPERBUILD "Build ${EXTENSION_NAME} and the projects it depends on." ON)
mark_as_advanced(${EXTENSION_NAME}_SUPERBUILD)
if(${EXTENSION_NAME}_SUPERBUILD)
  include("${CMAKE_CURRENT_SOURCE_DIR}/SuperBuild.cmake")
  return()
endif()

#-----------------------------------------------------------------------------
add_subdirectory(PortPlacement)
add_subdirectory(AutoPortPlacement)

#-----------------------------------------------------------------------------
include(${Slicer_EXTENSION_GENERATE_CONFIG})
include(${Slicer_EXTENSION_CPACK})
