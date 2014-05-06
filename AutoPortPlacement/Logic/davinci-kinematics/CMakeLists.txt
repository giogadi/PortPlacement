set(DAVINCI_PARAMS_FILE "${CMAKE_CURRENT_SOURCE_DIR}/../data/davinci-parameters.xml")

configure_file(davinci.h.in davinci.h)

add_library(davinci STATIC davinci.cxx)

target_link_libraries(davinci collisions ${VTK_LIBRARIES} ${ITK_LIBRARIES})

