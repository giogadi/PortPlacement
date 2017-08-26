add_library(collisions STATIC
  collisions.cxx
  collisions.h
  )
target_link_libraries(collisions PUBLIC Eigen3::Eigen)

set_property(TARGET collisions PROPERTY POSITION_INDEPENDENT_CODE ON)

set_property(GLOBAL APPEND PROPERTY Slicer_TARGETS collisions)
