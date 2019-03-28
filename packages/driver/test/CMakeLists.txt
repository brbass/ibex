macro(include_executable test_exec test_src)
    include_directories(${global_include_directories})
    set(test_dependencies utilities angular_discretization energy_discretization data solid_geometry spatial_discretization operator transport solver driver)
    
    add_executable(${test_exec} ${test_src})
    target_link_libraries(${test_exec} ${test_dependencies})
    install(TARGETS ${test_exec} DESTINATION ${CMAKE_INSTALL_PREFIX}/test)
endmacro()

macro(include_test test_name test_exec test_args)
    add_test(${test_name} ${test_exec} ${test_args})
endmacro()

include_executable(tst_driver tst_Driver.cc)

include_test(tst_manufactured_constant tst_driver ${CMAKE_CURRENT_SOURCE_DIR}/input/manufactured_constant.xml)

add_subdirectory(input)
