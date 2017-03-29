set(library_name solver)

set(library_include_directories ${global_include_directories} ${global_trilinos_include_directories})
set(library_link_directories ${global_trilinos_link_directories})
set(library_dependencies external utilities angular_discretization energy_discretization spatial_discretization data operator transport ${global_trilinos_link_libraries})

file(GLOB src *.cc *.hh)

include_directories(${library_include_directories})
link_directories(${library_link_directories})

add_library(${library_name} ${src})

target_link_libraries(${library_name} ${library_dependencies})

install(TARGETS ${library_name} DESTINATION ${CMAKE_INSTALL_PREFIX}/lib)

add_subdirectory(test)
