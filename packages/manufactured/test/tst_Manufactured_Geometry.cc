#include <iostream>
#include <memory>
#include <vector>

#include "Angular_Discretization.hh"
#include "Angular_Discretization_Factory.hh"
#include "Check_Equality.hh"
#include "Cross_Section.hh"
#include "Energy_Discretization.hh"
#include "Manufactured_Constant_Cross_Sections.hh"
#include "Manufactured_Constant_Solution.hh"
#include "Manufactured_Cross_Sections.hh"
#include "Manufactured_Solid_Geometry.hh"
#include "Manufactured_Solution.hh"
#include "Material.hh"
#include "Boundary_Source.hh"

using namespace std;
namespace ce = Check_Equality;

class Constant_Cross_Sections : public Manufactured_Cross_Sections
{
    Constant_Cross_Sections(shared_ptr<Angular_Discretization> angular,
                            shared_ptr<Energy_Discretization> energy):
        Manufactured_Cross_Sections(angular,
                                   energy)
    {
    }
    
    virtual void get_cross_sections(vector<double> const &position,
                                    vector<double> &sigma_t,
                                    vector<double> &sigma_s) const override
    {
        sigma_t = {1.0, 2.0};
        sigma_s = {0.5, 0.0, 0.0, 0.5};
        sigma_s.resize(angular_->number_of_scattering_moments() * energy_->number_of_groups() * energy_->number_of_groups(), 0);
    }
};

class Constant_Solution : public Manufactured_Solution
{
public:

    Constant_Solution(shared_ptr<Angular_Discretization> angular,
                      shared_ptr<Energy_Discretization> energy):
        Manufactured_Solution(angular,
                              energy)
    {
    }

    virtual vector<double> get_solution(vector<double> const &position) const override
    {
        vector<double> solution = {2., 1.};
        solution.resize(angular_->number_of_moments() * energy_->number_of_groups(), 0);
        return solution;
    }

    virtual vector<double> get_grad_solution(vector<double> const &position) const override
    {
        return vector<double>(angular_->number_of_moments() * energy_->number_of_groups() * angular_->dimension(), 0);
    }

};

int test_constant(int dimension)
{
    int checksum = 0;
    
    // Get angular discretization
    int angular_rule = dimension == 1 ? 256 : 3;
    int number_of_scattering_moments = 2;
    Angular_Discretization_Factory angular_factory;
    shared_ptr<Angular_Discretization> angular
        = angular_factory.get_angular_discretization(dimension,
                                                     number_of_scattering_moments,
                                                     angular_rule);
    int number_of_moments = angular->number_of_moments();
    int number_of_ordinates = angular->number_of_ordinates();
    double angular_normalization = angular->angular_normalization();

    // Get energy discretization
    int number_of_groups = 2;
    shared_ptr<Energy_Discretization> energy
        = make_shared<Energy_Discretization>(number_of_groups);
    
    // Get solution
    vector<double> solution_data = {2, 1};
    solution_data.resize(number_of_groups * number_of_moments, 0);
    shared_ptr<Manufactured_Solution> solution
        = make_shared<Manufactured_Constant_Solution>(angular,
                                                      energy,
                                                      solution_data);
    
    // Get cross sections
    vector<double> sigma_t_data = {1.0, 2.0};
    vector<double> sigma_s_data = {0.5, 0.0, 0.0, 0.5};
    sigma_s_data.resize(number_of_scattering_moments * number_of_groups * number_of_groups);
    shared_ptr<Manufactured_Cross_Sections> cross_sections
        = make_shared<Manufactured_Constant_Cross_Sections>(angular,
                                                            energy,
                                                            sigma_t_data,
                                                            sigma_s_data);
    
    // Get geometry
    shared_ptr<Solid_Geometry> solid
        = make_shared<Manufactured_Solid_Geometry>(angular,
                                                   energy,
                                                   solution,
                                                   cross_sections);

    // Setup data
    double tolerance = 1e-14;
    vector<double> position(dimension, 0.5);
    vector<double> expected; 
    
    // Check material for errors
    shared_ptr<Material> material = solid->material(position);
    shared_ptr<Cross_Section> internal_source = material->internal_source();
    expected.assign(number_of_groups * number_of_moments, 0);
    expected[0] = 2. * (1.0 - 0.5);
    expected[1] = 1. * (2.0 - 0.5);
    
    if (!ce::approx(internal_source->data(), expected, tolerance))
    {
        checksum += 1;
        cout << "internal source for dimension ";
        cout << dimension;
        cout << " incorrect" << endl;
    }
    
    // Check boundary source for errors
    shared_ptr<Boundary_Source> boundary_source = solid->boundary_source(position);
    expected.assign(number_of_groups * number_of_ordinates, 0);
    for (int o = 0; o < number_of_ordinates; ++o)
    {
        for (int g = 0; g < number_of_groups; ++g)
        {
            double val = g == 0 ? 2. : 1.;
            expected[g + number_of_groups * o] = val / angular_normalization;
        }
    }
    if (!ce::approx(boundary_source->data(), expected, tolerance))
    {
        checksum += 1;
        cout << "internal source for dimension ";
        cout << dimension;
        cout << "incorrect" << endl;
    }

    return checksum;
}

int main()
{
    int checksum = 0;

    for (int i : {1, 2, 3})
    {
        checksum += test_constant(i);
    }
    
    return checksum;
}
