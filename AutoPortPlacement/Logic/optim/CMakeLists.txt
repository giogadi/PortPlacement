add_library(optim STATIC optim.cxx)

target_link_libraries(optim
  ${VTK_LIBRARIES}
  davinci
  nlopt)
