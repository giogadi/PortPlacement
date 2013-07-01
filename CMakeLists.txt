cmake_minimum_required(VERSION 2.8.9)

project(PortPlacement)

#-----------------------------------------------------------------------------
set(EXTENSION_HOMEPAGE "http://www.slicer.org/slicerWiki/index.php/Documentation/Nightly/Extensions/PortPlacement")
set(EXTENSION_CATEGORY "Port Placement")
set(EXTENSION_CONTRIBUTORS "Luis G. Torres (UNC)")
set(EXTENSION_DESCRIPTION "Assists in the planning of surgical port placement in a laparoscopic procedure.")
set(EXTENSION_ICONURL "http://viewvc.slicer.org/viewvc.cgi/Slicer4/trunk/Extensions/Testing/PortPlacement/PortPlacement.png?revision=21747&view=co")
set(EXTENSION_SCREENSHOTURLS "http://wiki.slicer.org/slicerWiki/images/e/e2/Slicer-r19441-PortPlacement-screenshot.png")

#-----------------------------------------------------------------------------
find_package(Slicer REQUIRED)
include(${Slicer_USE_FILE})

#-----------------------------------------------------------------------------
add_subdirectory(PortPlacement)

#-----------------------------------------------------------------------------
include(${Slicer_EXTENSION_CPACK})
