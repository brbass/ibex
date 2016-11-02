set(subdirectories external utilities energy_discretization angular_discretization data solid_geometry spatial_discretization operator transport preconditioner solver driver ibex)

foreach(subdirectory ${subdirectories})
    add_subdirectory(${subdirectory})
endforeach()
