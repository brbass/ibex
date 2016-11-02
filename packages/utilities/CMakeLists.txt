set(library_name utilities)
set(library_dependencies external)

include_directories(${global_include_directories})
include_directories(${global_trilinos_include_directories})
link_directories(${global_trilinos_link_directories})

file(GLOB src *.cc *.hh)
add_library(${library_name} ${src})
target_link_libraries(${library_name} ${library_dependencies})
target_link_libraries(${library_name} ${global_trilinos_link_libraries})

install(TARGETS ${library_name} DESTINATION ${CMAKE_INSTALL_PREFIX}/lib)

add_subdirectory(test)
