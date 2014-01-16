# build NLOpt
include(ExternalProject)
ExternalProject_Add(nlopt
  PREFIX nlopt
  URL http://ab-initio.mit.edu/nlopt/nlopt-2.4.1.tar.gz
  URL_MD5 "ed72dd533b641d7e1f93eb790dee828c"
  CONFIGURE_COMMAND ../nlopt/configure --with-cxx --enable-shared --without-python --without-guile --without-octave --without-matlab --prefix=${CMAKE_CURRENT_BINARY_DIR}/nlopt)

set(NLopt_INCLUDE_DIRS "${CMAKE_CURRENT_BINARY_DIR}/nlopt/include")
set(NLopt_LIBRARIES "${CMAKE_CURRENT_BINARY_DIR}/nlopt/lib/${CMAKE_SHARED_LIBRARY_PREFIX}nlopt_cxx${CMAKE_SHARED_LIBRARY_SUFFIX}")

include_directories(${NLopt_INCLUDE_DIRS})

add_library(optim optim.cxx)

add_dependencies(optim nlopt)

target_link_libraries(optim
  davinci
  ${NLopt_LIBRARIES})
