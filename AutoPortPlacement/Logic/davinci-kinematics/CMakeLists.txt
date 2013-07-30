include_directories(../pugixml-1.2/src/)

add_library(davinci davinci.cxx ../pugixml-1.2/src/pugixml.cpp)

target_link_libraries(davinci collisions)

