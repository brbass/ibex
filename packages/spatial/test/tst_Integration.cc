#include <cmath>
#include <iomanip>
#include <iostream>
#include <mpi.h>

#include "Angular_Discretization.hh"
#include "Angular_Discretization_Factory.hh"
#include "Angular_Discretization_Parser.hh"
#include "Boundary_Source.hh"
#include "Boundary_Source_Parser.hh"
#include "Cartesian_Plane.hh"
#include "Check_Equality.hh"
#include "Constructive_Solid_Geometry.hh"
#include "Constructive_Solid_Geometry_Parser.hh"
#include "Cross_Section.hh"
#include "Cylinder_2D.hh"
#include "Energy_Discretization.hh"
#include "Energy_Discretization_Parser.hh"
#include "Material.hh"
#include "Material_Factory.hh"
#include "Material_Parser.hh"
#include "Region.hh"
#include "Weak_Spatial_Discretization.hh"
#include "Weak_Spatial_Discretization_Factory.hh"
#include "Weak_Spatial_Discretization_Parser.hh"

namespace ce = Check_Equality;
using namespace std;

int main()
{
    int checksum = 0;
    
    // Initialize MPI
    MPI_Init(&argc, &argv);

    // Run tests
    
    // Close MPI
    MPI_Finalize();
    
    return checksum;
}
