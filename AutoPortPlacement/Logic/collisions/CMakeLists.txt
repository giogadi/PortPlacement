add_library(collisions STATIC
  collisions.cxx
  collisions.h
  )
set_property(TARGET collisions PROPERTY POSITION_INDEPENDENT_CODE ON)
