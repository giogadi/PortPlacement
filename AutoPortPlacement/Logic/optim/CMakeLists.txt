# build NLOpt
include(ExternalProject)
ExternalProject_Add(nlopt
  PREFIX nlopt
  URL http://ab-initio.mit.edu/nlopt/nlopt-2.3.tar.gz
  URL_MD5 "811a9f1c7a7f879c7d7b4caa059eb8d6"
  CONFIGURE_COMMAND ../nlopt/configure --with-cxx --enable-shared --without-python --without-guile --without-octave --without-matlab --prefix=${CMAKE_CURRENT_BINARY_DIR}/nlopt)

set(NLopt_INCLUDE_DIRS "${CMAKE_CURRENT_BINARY_DIR}/nlopt/include")
set(NLopt_LIBRARIES "${CMAKE_CURRENT_BINARY_DIR}/nlopt/lib/${CMAKE_SHARED_LIBRARY_PREFIX}nlopt_cxx${CMAKE_SHARED_LIBRARY_SUFFIX}")

include_directories(${NLopt_INCLUDE_DIRS})

add_library(optim optim.cxx)

add_dependencies(optim nlopt)

target_link_libraries(optim 
  davinci
  ${NLopt_LIBRARIES})
  
