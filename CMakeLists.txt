cmake_minimum_required (VERSION 3.1)

####################
# Compiler options #
####################

set(CMAKE_BUILD_TYPE Debug)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -O3")

#################
# Find Trilinos #
#################

set(global_trilinos_include_directories "")
set(global_trilinos_link_directories "")
set(global_trilinos_link_libraries "")

if (TRUE)
   # Find Trilinos
   find_package(Trilinos REQUIRED)
   message("Trilinos found in ${Trilinos_DIR}")

   # Set compiler options
   set(CMAKE_CXX_COMPILER ${Trilinos_CXX_COMPILER} )
   set(CMAKE_C_COMPILER ${Trilinos_C_COMPILER} )
   set(CMAKE_CXX_FLAGS  "${Trilinos_CXX_COMPILER_FLAGS} ${CMAKE_CXX_FLAGS}")
   set(CMAKE_C_FLAGS  "${Trilinos_C_COMPILER_FLAGS} ${CMAKE_C_FLAGS}")

   # Get list of include and link values
   set(global_trilinos_include_directories ${global_trilinos_include_directories} ${Trilinos_INCLUDE_DIRS} ${Trilinos_TPL_INCLUDE_DIRS})
   set(global_trilinos_link_directories ${global_trilinos_link_directories} ${Trilinos_LIBRARY_DIRS} ${Trilinos_TPL_LIBRARY_DIRS})
   set(global_trilinos_link_libraries ${global_trilinos_link_libraries} ${Trilinos_LIBRARIES} ${Trilinos_TPL_LIBRARIES})
   
else()
    set(trilinosPackages Amesos AztecOO Epetra Ifpack Teuchos Belos)
    foreach(trilinosPackage ${trilinosPackages})
        # Find package
        find_package(${trilinosPackage} REQUIRED)
        message("${trilinosPackage} found in ${${trilinosPackage}_DIR}")
        
        # Set compiler
        set(CMAKE_CXX_COMPILER ${${trilinosPackage}_CXX_COMPILER} )
        set(CMAKE_C_COMPILER ${${trilinosPackage}_C_COMPILER} )
        set(CMAKE_CXX_FLAGS  "${${trilinosPackage}_CXX_COMPILER_FLAGS} ${CMAKE_CXX_FLAGS}")
        set(CMAKE_C_FLAGS  "${${trilinosPackage}_C_COMPILER_FLAGS} ${CMAKE_C_FLAGS}")
        
        # Get lists of include and link values
        set(global_trilinos_include_directories ${global_trilinos_include_directories} ${${trilinosPackage}_INCLUDE_DIRS} ${${trilinosPackage}_TPL_INCLUDE_DIRS})
        set(global_trilinos_link_directories ${global_trilinos_link_directories} ${${trilinosPackage}_LIBRARY_DIRS} ${${trilinosPackage}_TPL_LIBRARY_DIRS})
        set(global_trilinos_link_libraries ${global_trilinos_link_libraries} ${${trilinosPackage}_LIBRARIES} ${${trilinosPackage}_TPL_LIBRARIES})
    endforeach()
endif()

###########
# Project #
###########

project(Ibex)

set(CMAKE_INSTALL_PREFIX ${CMAKE_CURRENT_SOURCE_DIR}/bin CACHE PATH "Installation directory")

if(TRUE)
    message("\nCMake CXX Compiler:\t${CMAKE_CXX_COMPILER}")
    message("\nCMake C Compiler:\t${CMAKE_C_COMPILER}")
    message("\nCMake CXX Flags:\t${CMAKE_CXX_FLAGS}")
    message("\nCMake C Flags:\t${CMAKE_C_FLAGS}")
    message("\nTrilinos include directories:\n${global_trilinos_include_directories}")
    message("\nTrilinos link directories:\n${global_trilinos_link_directories}")
    message("\nTrilinos link libraries:\n${global_trilinos_link_libraries}")
endif()

##########
# OpenMP #
##########

include(FindOpenMP)
if(OPENMP_FOUND)
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS} -DENABLE_OPENMP")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS} -DENABLE_OPENMP")
  set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
  message("OpenMP found")
else(OPENMP_FOUND)
  message("WARNING: OpenMP could not be found.")
endif(OPENMP_FOUND)

##########################
# Get header directories #
##########################

# Finds all header files recursively in given directories
macro(header_directories return_list)
    file(GLOB_RECURSE header_files packages/*.hh)
    set(header_directories "")
    foreach(header_file ${header_files})
        get_filename_component(header_directory ${header_file} DIRECTORY)
        set(header_directories ${header_directories} ${header_directory})
    endforeach()
    list(REMOVE_DUPLICATES header_directories)
    set(${return_list} ${header_directories})
endmacro()

# Lists of include and link values
set(global_include_directories "")

header_directories(global_include_directories)

########################
# Include source files #
########################

# Make lists global
set(global_include_directories ${global_include_directories} CACHE INTERNAL "global_include_directories")
set(global_trilinos_directories ${global_trilinos_directories} CACHE INTERNAL "global_trilinos_directories")
set(global_trilinos_libraries ${global_trilinos_libraries} CACHE INTERNAL "global_trilinos_libraries")
set(global_trilinos_link_libraries ${global_trilinos_link_libraries} CACHE INTERNAL "global_trilinos_link_libraries")

enable_testing()

set(subdirectories packages)
foreach(subdirectory ${subdirectories})
    add_subdirectory(${subdirectory})
endforeach()
