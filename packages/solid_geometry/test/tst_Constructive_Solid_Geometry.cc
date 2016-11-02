#include <iostream>
#include <memory>

#include "Angular_Discretization_Parser.hh"
#include "Boundary_Source_Parser.hh"
#include "Check_Equality.hh"
#include "Energy_Discretization_Parser.hh"
#include "Material_Parser.hh"
#include "Constructive_Solid_Geometry_Parser.hh"
#include "Vector_Functions.hh"

using namespace std;

namespace ce = Check_Equality;
namespace vf = Vector_Functions;

// Only checks 2D region find and optical depth

shared_ptr<Constructive_Solid_Geometry> get_solid_geometry(string xml_input_filename)
{
    pugi::xml_document input_document;
    
    if (!input_document.load_file(xml_input_filename.c_str()))
    {
        cout << "Could not open xml input file \"" + xml_input_filename + "\"" << endl;
        return shared_ptr<Constructive_Solid_Geometry>();
    }
    
    pugi::xml_node input_file = input_document.child("input");
    
    Angular_Discretization_Parser angular_parser(input_file);
    shared_ptr<Angular_Discretization> angular = angular_parser.get_ptr();
    
    Energy_Discretization_Parser energy_parser(input_file);
    shared_ptr<Energy_Discretization> energy = energy_parser.get_ptr();
    
    Material_Parser material_parser(input_file,
                                    angular,
                                    energy);
    vector<shared_ptr<Material> > materials = material_parser.get_ptr();
    
    Boundary_Source_Parser boundary_parser(input_file,
                                           angular,
                                           energy);
    vector<shared_ptr<Boundary_Source> > boundary_sources = boundary_parser.get_ptr();
    
    Constructive_Solid_Geometry_Parser solid_parser(input_file,
                                                     materials,
                                                     boundary_sources);
    return solid_parser.get_ptr();
}

int test_cylindrical_pincell(string input_folder)
{
    int checksum = 0;
    
    string xml_input_filename = input_folder + "/cylindrical_pincell.xml";
    
    shared_ptr<Constructive_Solid_Geometry> solid = get_solid_geometry(xml_input_filename);
    
    if (!solid)
    {
        checksum += 1;
        return checksum;
    }

    int const fuel_index = 0;
    int const moderator_index = 1;
    double const fuel_diameter = 2;
    double const pincell_length = 4;
    double const sigma_t_fuel = 1.0;
    double const sigma_t_moderator = 2.0;
    vector<double> const origin = {0, 0};

    // Test for Cartesian boundaries
    
    if (!solid->cartesian_boundaries())
    {
        checksum += 1;
        cout << "cylindrical_pincell, cartesian boundaries not found." << endl;
    }

    // Test finding of region

    {
        vector<vector<double> > const locations_to_check
            = {{0, 0},
               {1.5, 0},
               {0, 1.5},
               {-1.5, 0},
               {0, -1.5},
               {0, 2},
               {-2, 0},
               {-2.1, 0},
               {0, 2.1}};
    
        vector<int> const expected_region = {fuel_index,
                                             moderator_index,
                                             moderator_index,
                                             moderator_index,
                                             moderator_index,
                                             moderator_index,
                                             moderator_index,
                                             Constructive_Solid_Geometry::NO_REGION,
                                             Constructive_Solid_Geometry::NO_REGION};

        int num_locations = locations_to_check.size();

        for (int i = 0; i < num_locations; ++i)
        {
            int region = solid->find_region_including_surface(locations_to_check[i]);

            if (expected_region[i] != region)
            {
                cout << "cylindrical pincell: region incorrect for case " << i << endl;
                checksum += 1;
            }
        }
    }

    // Test intersection (not completed)
    
    // Test boundary (not completed)
    
    // Test optical distance

    {
        vector<vector<double> > const initial_positions
            = {{0, 2},
               {2, -2},
               {-2, -2},
               {1./sqrt(2), 1./sqrt(2)}};
        vector<vector<double> > const final_positions
            = {{0, -2},
               {-2, 2},
               {2, -2},
               {-1./sqrt(2), -1./sqrt(2)}};

        vector<double> const expected_optical_distance
            = {2 * sigma_t_fuel + 2 * sigma_t_moderator,
               2 * sigma_t_fuel + 2 * (2 * sqrt(2) - 1) * sigma_t_moderator,
               4 * sigma_t_moderator,
               2 * sigma_t_fuel};

        int num_cases = initial_positions.size();

        for (int i = 0; i < num_cases; ++i)
        {
            vector<double> optical_distance;

            solid->optical_distance(initial_positions[i],
                                    final_positions[i],
                                    optical_distance);
            
            if (!ce::approx(optical_distance[0], expected_optical_distance[i], 1e-11))
            {
                cout << "cylindrical pincell: optical distance failed for case " << i << endl;
                cout << "\texpected: " << expected_optical_distance[i];
                cout << "\tcalculated: " << optical_distance[0] << endl;
                checksum += 1;
            }
        }
    }
    
    return checksum;
}

int main(int argc, char **argv)
{
    int checksum = 0;

    if (argc != 2)
    {
        cerr << "usage: tst_Constructive_Solid_Geometry [input_folder]" << endl;
        return 1;
    }

    string input_folder = argv[1];
    
    checksum += test_cylindrical_pincell(input_folder);

    return checksum;
}
