set(subdirectories vera_multiphysics slab_heat manufactured_heat) # prospectus

foreach(subdirectory ${subdirectories})
    add_subdirectory(${subdirectory})
endforeach()
