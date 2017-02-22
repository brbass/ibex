macro(include_test test_exec test_src test_args)
    include_directories(${global_include_directories})
    set(test_dependencies utilities solid_geometry)
    
    add_executable(${test_exec} ${test_src})
    target_link_libraries(${test_exec} ${test_dependencies})
    add_test(${test_exec} ${test_exec} ${test_args})
    install(TARGETS ${test_exec} DESTINATION ${CMAKE_INSTALL_PREFIX}/test)
endmacro()

include_test(tst_surfaces tst_Surfaces.cc "")
include_test(tst_regions tst_Regions.cc "")
include_test(tst_solid_geometry tst_Constructive_Solid_Geometry.cc ${CMAKE_CURRENT_SOURCE_DIR}/input)

add_subdirectory(input)
