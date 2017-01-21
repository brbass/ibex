#include <iostream>
#include <memory>

#include "Cross_Section.hh"
#include "Cylinder_3D.hh"
#include "Energy_Discretization.hh"
#include "LDFE_Quadrature.hh"
#include "Material.hh"
#include "Plane_3D.hh"
#include "Random_Number_Generator.hh"
#include "Region.hh"
#include "Sphere_3D.hh"
#include "Vector_Functions.hh"

namespace vf = Vector_Functions;

using namespace std;

int test_spherical_region()
{
    int dimension = 3;
    
    int checksum = 0;
    
    shared_ptr<Angular_Discretization> angular
        = make_shared<LDFE_Quadrature>(dimension,
                                       1,
                                       1);
    shared_ptr<Energy_Discretization> energy
        = make_shared<Energy_Discretization>(1);

    Cross_Section::Dependencies none_group;
    none_group.energy = Cross_Section::Dependencies::Energy::GROUP;
    Cross_Section::Dependencies scattering_group2;
    scattering_group2.angular = Cross_Section::Dependencies::Angular::SCATTERING_MOMENTS;
    scattering_group2.energy = Cross_Section::Dependencies::Energy::GROUP_TO_GROUP;
    
    shared_ptr<Cross_Section> sigma_t
        = make_shared<Cross_Section>(none_group,
                                     angular,
                                     energy,
                                     vector<double>({1}));
    shared_ptr<Cross_Section> sigma_s
        = make_shared<Cross_Section>(scattering_group2,
                                     angular,
                                     energy,
                                     vector<double>({0}));
    shared_ptr<Cross_Section> nu
        = make_shared<Cross_Section>(none_group,
                                     angular,
                                     energy,
                                     vector<double>({0}));
    shared_ptr<Cross_Section> sigma_f
        = make_shared<Cross_Section>(none_group,
                                     angular,
                                     energy,
                                     vector<double>({0}));
    shared_ptr<Cross_Section> chi
        = make_shared<Cross_Section>(none_group,
                                     angular,
                                     energy,
                                     vector<double>({0}));
    shared_ptr<Cross_Section> internal_source
        = make_shared<Cross_Section>(none_group,
                                     angular,
                                     energy,
                                     vector<double>({0}));

    shared_ptr<Material> material
        = make_shared<Material>(0,
                                angular,
                                energy,
                                sigma_t,
                                sigma_s,
                                nu,
                                sigma_f,
                                chi,
                                internal_source);
    
    int index = 0;
    Surface::Surface_Type surface_type = Surface::Surface_Type::BOUNDARY;
    double radius = 5;
    vector<double> origin = {0, 0, 0};

    vector<Surface::Relation> surface_relations
        = {Surface::Relation::INSIDE};
    vector<shared_ptr<Surface> > surfaces
        = {make_shared<Sphere_3D>(index,
                                  surface_type,
                                  radius,
                                  origin)};
    
    shared_ptr<Region> region
        = make_shared<Region>(index,
                              material,
                              surface_relations,
                              surfaces);

    Random_Number_Generator<double> rng(-radius, // min
                                        radius, // max
                                        7); // seed
    int num_tests = 1000;

    for (int i = 0; i < num_tests; ++i)
    {
        vector<double> position = rng.vector(dimension);

        Region::Relation expected_result;
        
        if (vf::magnitude(vf::subtract(position, origin)) < 5)
        {
            expected_result = Region::Relation::INSIDE;
        }
        else
        {
            expected_result = Region::Relation::OUTSIDE;
        }

        if (region->relation(position) != expected_result)
        {
            cout << "spherical region relation failed" << endl;
            checksum += 1;
        }
    }
    
    return checksum;
}

int main()
{
    int checksum = 0;
    
    checksum += test_spherical_region();

    return checksum;
}
