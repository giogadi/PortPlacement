include_directories(../pugixml-1.2/src/)

set(DAVINCI_PARAMS_FILE "${CMAKE_CURRENT_SOURCE_DIR}/../tests/data/davinci-parameters.xml")

configure_file(davinci.h.in davinci.h)

add_library(davinci davinci.cxx ../pugixml-1.2/src/pugixml.cpp)

target_link_libraries(davinci collisions)

