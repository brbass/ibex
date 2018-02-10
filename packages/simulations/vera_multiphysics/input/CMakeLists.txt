file(GLOB input_files *.xml)

foreach(input_path ${input_files})
    get_filename_component(input_file ${input_path} NAME)
    configure_file(${input_path} ${CMAKE_CURRENT_BINARY_DIR}/${input_file} COPYONLY)
endforeach(input_path)

install(FILES ${input_files} DESTINATION ${CMAKE_INSTALL_PREFIX}/simulations/vera_multiphysics/input)
