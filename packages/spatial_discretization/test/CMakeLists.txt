macro(include_test test_exec test_src test_args)
    include_directories(${global_include_directories})
    set(test_dependencies utilities solid_geometry spatial_discretization)
    
    add_executable(${test_exec} ${test_src})
    target_link_libraries(${test_exec} ${test_dependencies})
    add_test(${test_exec} ${test_exec} ${test_args})
    install(TARGETS ${test_exec} DESTINATION ${CMAKE_INSTALL_PREFIX}/test)
endmacro()

include_test(tst_rbf tst_RBF.cc "")
include_test(tst_distance tst_Distance.cc ${CMAKE_CURRENT_SOURCE_DIR}/../../solid_geometry/test/input)
include_test(tst_meshless_function tst_Meshless_Function.cc "")
# include_test(tst_rbf_interpolation tst_RBF_Interpolation.cc ${CMAKE_CURRENT_SOURCE_DIR}/input)

add_subdirectory(input)
