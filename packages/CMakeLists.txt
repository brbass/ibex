set(subdirectories external utilities energy angular data solid spatial operator transport solver manufactured heat driver ibex simulations)

foreach(subdirectory ${subdirectories})
    add_subdirectory(${subdirectory})
endforeach()
