set(subdirectories vera_multiphysics slab_heat manufactured_heat output_points quadratures ray_effects purely_absorbing) # prospectus

foreach(subdirectory ${subdirectories})
    add_subdirectory(${subdirectory})
endforeach()
