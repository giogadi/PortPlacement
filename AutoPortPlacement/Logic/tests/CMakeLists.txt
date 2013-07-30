add_executable(da-vinci-params-load-test 
  da-vinci-params-load-test.cxx)
target_link_libraries(da-vinci-params-load-test davinci)

add_executable(jacobian-test
  jacobian-test.cxx)
target_link_libraries(jacobian-test davinci)

add_executable(da-vinci-passive-ik-test
  da-vinci-passive-ik-test.cxx)
target_link_libraries(da-vinci-passive-ik-test davinci)

add_executable(min-distance-test
  min-distance-test.cxx)
target_link_libraries(min-distance-test collisions)

set(DAVINCI_PARAMS_FILE "${CMAKE_CURRENT_SOURCE_DIR}/data/davinci-parameters.xml")

add_test(daVinci_load_params_test da-vinci-params-load-test ${DAVINCI_PARAMS_FILE})
add_test(daVinci_jacobian_test jacobian-test ${DAVINCI_PARAMS_FILE})
add_test(daVinci_passive_ik_test da-vinci-passive-ik-test ${DAVINCI_PARAMS_FILE})
add_test(collisions_distance_test min-distance-test)
