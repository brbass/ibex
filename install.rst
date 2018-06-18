--------------------
Install Trilinos
--------------------

The only external dependency of ibex is Trilinos. If Trilinos is installed on your system with the appropriate packages (Epetra, Anasazi, AztecOO, Ifpack and Belos), this section can be skipped.

To install Trilinos with only the needed packages, do the following:

- Create a folder, e.g., /software/trilinos

- Clone the Trilinos repository from https://github.com/trilinos/Trilinos.git into, e.g., /software/trilinos/Trilinos

- Make a build directory, e.g., /software/trilinos/build

- Run cmake in the build directory using the following script (modify open-mpi, build, and install paths if needed):

::

   cmake \
   -DTPL_ENABLE_MPI=ON \
   -DMPI_BASE_DIR=/usr/local/Cellar/open-mpi/3.0.1 \
   -DTrilinos_ENABLE_ALL_PACKAGES=OFF \
   -DTrilinos_ENABLE_Epetra=ON \
   -DTrilinos_ENABLE_Anasazi=ON \
   -DTrilinos_ENABLE_AztecOO=ON \
   -DTrilinos_ENABLE_Ifpack=ON \
   -DTrilinos_ENABLE_Belos=ON \
   -DCMAKE_INSTALL_PREFIX=/software/trilinos/bin \
   /software/trilinos/Trilinos

- Make the program, e.g. from the build directory,

make install -j4

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

ibex is designed to be built in place. For a build outside of the source directory, the install directory will need to be changed in the CMakeLists.txt file in this directory.

- Choose release type (Debug or Release) in the CMakeLists.txt file in this directory (./CMakeLists.txt)
  - For debugging, remove the -O3 symbol

- Make build directory, e.g. ./build

- From build directory, run the following:

::

   cmake ..; cmake ..; make install -j4

- This should install ibex into the directory ./bin

- Optionally, add the following scripts (with modified paths) to .bashrc to add the "ibex" executable and the Python scripts

::

   export PATH="/home/brbass/code/ibex/bin:$PATH"
   export PYTHONPATH="/home/brbass/code/ibex/scripts:$PYTHONPATH"

#############
# Run tests #
#############

- To test ibex, run the following from the build directory:

::

   ctest -j4

