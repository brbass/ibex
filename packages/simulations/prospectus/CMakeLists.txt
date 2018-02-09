macro(include_test test_exec test_src test_args)
    include_directories(${global_include_directories})
    set(test_dependencies utilities angular_discretization energy_discretization data solid_geometry spatial_discretization operator transport solver)
    
    add_executable(${test_exec} ${test_src})
    target_link_libraries(${test_exec} ${test_dependencies})
    if(${is_test})
        add_test(${test_exec} ${test_exec} ${test_args})
    endif()
    install(TARGETS ${test_exec} DESTINATION ${CMAKE_INSTALL_PREFIX}/test)
endmacro()

include_test(tst_prospectus tst_Prospectus.cc "" false)
