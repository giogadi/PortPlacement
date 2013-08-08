add_executable(da-vinci-params-load-test 
  da-vinci-params-load-test.cxx)
target_link_libraries(da-vinci-params-load-test davinci)

add_executable(intra-ik-test
  intra-ik-test.cxx)
target_link_libraries(intra-ik-test davinci)

add_executable(jacobian-test
  jacobian-test.cxx)
target_link_libraries(jacobian-test davinci)

add_executable(da-vinci-passive-ik-test
  da-vinci-passive-ik-test.cxx)
target_link_libraries(da-vinci-passive-ik-test davinci)

add_executable(min-distance-test
  min-distance-test.cxx)
target_link_libraries(min-distance-test collisions)

add_executable(passive-lr-optim-test
  passive-lr-optim-test.cxx)
target_link_libraries(passive-lr-optim-test optim)

add_executable(config-ut-test
  config-ut-test.cxx)
target_link_libraries(config-ut-test davinci)

add_executable(find-feasible-plan-test
  find-feasible-plan-test.cxx)
target_link_libraries(find-feasible-plan-test optim)

add_test(daVinci_load_params_test da-vinci-params-load-test)
add_test(daVinci_intra_ik_test intra-ik-test)
add_test(daVinci_jacobian_test jacobian-test)
add_test(daVinci_passive_ik_test da-vinci-passive-ik-test)
add_test(collisions_distance_test min-distance-test)
add_test(passive_lr_optim_test passive-lr-optim-test)
add_test(config_ut_test config-ut-test)
add_test(find_feasible_plan_test find-feasible-plan-test)
