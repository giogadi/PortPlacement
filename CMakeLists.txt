cmake_minimum_required(VERSION 3.5)

project(PortPlacement)

#-----------------------------------------------------------------------------
set(EXTENSION_HOMEPAGE "http://www.slicer.org/slicerWiki/index.php/Documentation/Nightly/Extensions/PortPlacement")
set(EXTENSION_CATEGORY "Port Placement")
set(EXTENSION_CONTRIBUTORS "Luis G. Torres (UNC)")
set(EXTENSION_DESCRIPTION "Assists in the planning of surgical port placement in a laparoscopic procedure.")
set(EXTENSION_ICONURL "http://viewvc.slicer.org/viewvc.cgi/Slicer4/trunk/Extensions/Testing/PortPlacement/PortPlacement.png?revision=21747&view=co")
set(EXTENSION_SCREENSHOTURLS "http://wiki.slicer.org/slicerWiki/images/e/e2/Slicer-r19441-PortPlacement-screenshot.png")

#
# Implementation agnostic fPIC
#
set( CMAKE_POSITION_INDEPENDENT_CODE ON )

#-----------------------------------------------------------------------------
# Extension dependencies
find_package(Slicer REQUIRED)
include(${Slicer_USE_FILE})

find_package(Eigen3 REQUIRED CONFIG)
mark_as_superbuild(Eigen3_DIR)

#-----------------------------------------------------------------------------
add_subdirectory(PortPlacement)
add_subdirectory(AutoPortPlacement)

#-----------------------------------------------------------------------------
include(${Slicer_EXTENSION_CPACK})
