set(library_name external)

file(GLOB src *.cc *hh)
add_library(${library_name} ${src}) 
install(TARGETS ${library_name} DESTINATION ${CMAKE_INSTALL_PREFIX}/lib)
