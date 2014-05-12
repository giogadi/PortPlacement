find_package(PkgConfig QUIET)
PKG_CHECK_MODULES(PC_NLopt QUIET nlopt)

find_path(NLopt_INCLUDE_DIR
  NAMES nlopt.h
  HINTS
    ${NLopt_DIR}
    ${PC_NLopt_INCLUDEDIR}
    ${PC_NLopt_INCLUDE_DIRS}
  PATH_SUFFIXES include
  )

find_library(NLopt_LIBRARY
  NAMES
    nlopt
    nlopt-0
    libnlopt
    libnlopt-0
  HINTS
    ${NLopt_DIR}
    ${PC_NLopt_LIBDIR}
    ${PC_NLopt_LIBRARY_DIRS}
  PATH_SUFFIXES
    lib
    lib64
  )

set(NLopt_INCLUDE_DIRS ${NLopt_INCLUDE_DIR})
set(NLopt_LIBRARIES ${NLopt_LIBRARY})

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(NLopt DEFAULT_MSG NLopt_INCLUDE_DIRS NLopt_LIBRARIES)
mark_as_advanced(NLopt_INCLUDE_DIR NLopt_LIBRARY)
