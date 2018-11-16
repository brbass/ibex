--------------------
Install Trilinos
--------------------

The only external dependency of ibex is Trilinos. If Trilinos is installed on your system with the appropriate packages (Epetra, Anasazi, AztecOO, Ifpack and Belos), this section can be skipped.

To install Trilinos with only the needed packages, do the following:

- Create a folder, e.g., /software/trilinos

- Clone the Trilinos repository from https://github.com/trilinos/Trilinos.git into, e.g., /software/trilinos/Trilinos

- Make a build directory, e.g., /software/trilinos/build

- Make sure MPI_BASE_DIR is set in your environment so Trilinos can find it

- Run cmake in the build directory using the following script (modify build and install paths appropriately):

::

   cmake \
   -DTPL_ENABLE_MPI=ON \
   -DTrilinos_ENABLE_ALL_PACKAGES=OFF \
   -DTrilinos_ENABLE_Epetra=ON \
   -DTrilinos_ENABLE_Anasazi=ON \
   -DTrilinos_ENABLE_AztecOO=ON \
   -DTrilinos_ENABLE_Ifpack=ON \
   -DTrilinos_ENABLE_Belos=ON \
   -DCMAKE_INSTALL_PREFIX=/software/trilinos/bin \
   /software/trilinos/Trilinos

- Run make from the build directory,

::

   make install -j32

- Add the following to your bashrc (with the correct paths):

::

   export TRILINOS_INC=/software/trilinos/bin/include
   export TRILINOS_LIB=/software/trilinos/bin/lib 
   export TRILINOS_ROOT=/software/trilinos/bin 
   export TRILINOS_LINK=/software/trilinos/bin/lib 
   export TRILINOS_INCLUDE=/software/trilinos/bin/include 
   export Trilinos_VERSION=12.8.1 
   export Trilinos_INCLUDE_DIRS="-I/software/trilinos/bin/include" 
   export Trilinos_LIBRARY_DIRS="-L/software/trilinos/bin/lib" 
   export CMAKE_INCLUDE_PATH="/software/trilinos/bin/include:$CMAKE_INCLUDE_PATH" 
   export CMAKE_LIBRARY_PATH="/software/trilinos/bin/lib:$CMAKE_LIBRARY_PATH" 
   export CMAKE_MODULE_PATH="/software/trilinos/bin/lib/cmake:$CMAKE_MODULE_PATH"

------------
Install ibex
------------

The following assumes that IBEX_SOURCE is set to the directory containing this file.

To build ibex, first set the release type (Debug or Release) in the CMakeLists.txt file in IBEX_SOURCE. For debugging, remove the -O3 symbol.

To build ibex in place, run the following (with an appropriate number of processors) in IBEX_SOURCE:

::
  
   mkdir build; cd build; cmake ..; cmake ..; make install -j32

This will install to the IBEX_SOURCE/bin directory

To make and install in a different directory,

- Create the build directory

- From the build directory, run

::

   cmake IBEX_SOURCE -DCMAKE_INSTALL_PREFIX=IBEX_BIN; cmake IBEX_SOURCE -DCMAKE_INSTALL_PREFIX=IBEX_BIN; make install -j32

Optionally, add the following scripts (with modified paths) to .bashrc to add the "ibex" executable and the Python scripts

::

   export PATH="/home/brbass/code/ibex/bin:$PATH"
   export PYTHONPATH="/home/brbass/code/ibex/scripts:$PYTHONPATH"

#############
# Run tests #
#############

- To test ibex, run the following from the build directory:

::

   ctest -j32

