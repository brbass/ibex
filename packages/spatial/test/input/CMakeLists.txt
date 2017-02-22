file(GLOB input_files *.xml)

install(FILES ${input_files} DESTINATION ${CMAKE_INSTALL_PREFIX}/test/input)

