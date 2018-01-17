set(subdirectories external utilities energy angular data solid spatial operator transport solver manufactured driver ibex)

foreach(subdirectory ${subdirectories})
    add_subdirectory(${subdirectory})
endforeach()
