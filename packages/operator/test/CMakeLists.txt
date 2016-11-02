macro(include_test test_exec test_src)
    include_directories(${global_include_directories})
    set(test_dependencies utilities angular_discretization energy_discretization spatial_discretization data operator)
    
    add_executable(${test_exec} ${test_src})
    target_link_libraries(${test_exec} ${test_dependencies})
    add_test(${test_exec} ${test_exec})
    install(TARGETS ${test_exec} DESTINATION ${CMAKE_INSTALL_PREFIX}/test)
endmacro()

include_test(tst_moment_discrete tst_Moment_Discrete.cc)
