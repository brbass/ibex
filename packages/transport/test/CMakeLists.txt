macro(include_test test_exec test_src test_args)
    include_directories(${global_include_directories})
    set(test_dependencies utilities angular_discretization energy_discretization data solid_geometry spatial_discretization operator transport)
    
    add_executable(${test_exec} ${test_src})
    target_link_libraries(${test_exec} ${test_dependencies})
    add_test(${test_exec} ${test_exec} ${test_args})
    install(TARGETS ${test_exec} DESTINATION ${CMAKE_INSTALL_PREFIX}/test)
endmacro()

include_test(tst_rbf_collocation tst_RBF_Collocation.cc ${CMAKE_CURRENT_SOURCE_DIR}/input)

add_subdirectory(input)
