add_library(optim STATIC
  optim.cxx
  optim.h
  )
set_property(TARGET optim PROPERTY POSITION_INDEPENDENT_CODE ON)
target_link_libraries(optim
  davinci
  ${NLopt_LIBRARIES}
  ${VTK_LIBRARIES}
  )
