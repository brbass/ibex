#include <iostream>
#include <string>

#include <mpi.h>

#include "Driver.hh"

using namespace std;

int main(int argc, char **argv)
{
    MPI_Init(&argc, &argv);

    if (argc != 2)
    {
        cerr << "usage: ibex [input.xml]" << endl;
        return 1;
    }
    
    string filename = argv[1];
    
    Driver driver(filename);

    MPI_Finalize();
}
