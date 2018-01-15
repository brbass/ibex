macro(include_test test_exec test_src test_args)
    include_directories(${global_include_directories})
    set(test_dependencies utilities solid_geometry manufactured)
    
    add_executable(${test_exec} ${test_src})
    target_link_libraries(${test_exec} ${test_dependencies})
    add_test(${test_exec} ${test_exec} ${test_args})
    install(TARGETS ${test_exec} DESTINATION ${CMAKE_INSTALL_PREFIX}/test)
endmacro()

include_test(tst_manufactured_geometry tst_Manufactured_Geometry.cc "")
