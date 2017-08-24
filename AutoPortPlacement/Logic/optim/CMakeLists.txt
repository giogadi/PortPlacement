
find_package(NLopt REQUIRED CONFIG)

add_library(optim STATIC
  optim.cxx
  optim.h
  )
set_property(TARGET optim PROPERTY POSITION_INDEPENDENT_CODE ON)
target_link_libraries(optim
  davinci
  ${NLOPT_LIBRARIES}
  ${VTK_LIBRARIES}
  )

target_include_directories(optim PUBLIC
  ${NLOPT_INCLUDE_DIRS}
  )

set_property(GLOBAL APPEND PROPERTY Slicer_TARGETS optim)
