set(executable_name ibex)
set(executable_dependencies driver)

include_directories(${global_include_directories})

file(GLOB src *.cc *.hh)
add_executable(${executable_name} ${src})
target_link_libraries(${executable_name} ${executable_dependencies})
install(TARGETS ${executable_name} DESTINATION ${CMAKE_INSTALL_PREFIX})
