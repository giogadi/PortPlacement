# build NLOpt
add_subdirectory(nlopt)
include_directories(${nlopt_INCLUDE_DIRS})

add_library(optim optim.cxx)

target_link_libraries(optim
  davinci
  nlopt)
