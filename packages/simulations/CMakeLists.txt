set(subdirectories vera_multiphysics) # prospectus

foreach(subdirectory ${subdirectories})
    add_subdirectory(${subdirectory})
endforeach()
