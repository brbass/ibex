macro(include_test test_exec test_src)
    include_directories(${global_include_directories})
    set(test_dependencies utilities angular_discretization)
    
    add_executable(${test_exec} ${test_src})
    target_link_libraries(${test_exec} ${test_dependencies})
    add_test(${test_exec} ${test_exec})
    install(TARGETS ${test_exec} DESTINATION ${CMAKE_INSTALL_PREFIX}/test)
endmacro()

include_test(tst_gauss tst_Gauss.cc)
include_test(tst_ldfe tst_LDFE.cc)
