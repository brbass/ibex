macro(include_executable test_exec test_src)
    set(test_dependencies utilities angular_discretization energy_discretization data solid_geometry spatial_discretization operator transport solver heat)
    
    include_directories(${global_include_directories})
    message(${test_src})
    add_executable(${test_exec} ${test_src})
    target_link_libraries(${test_exec} ${test_dependencies})
    install(TARGETS ${test_exec} DESTINATION ${CMAKE_INSTALL_PREFIX}/simulations/ray_effects)
endmacro()

file(GLOB src *.cc *.hh)

include_executable(ray_effects "${src}")
