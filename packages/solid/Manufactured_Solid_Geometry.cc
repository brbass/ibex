#include "Manufactured_Solid_Geometry.hh"

#include "Angular_Discretization.hh"
#include "Boundary_Source.hh"
#include "Energy_Discretization.hh"
#include "Manufactured_Solution.hh"
#include "Material.hh"
#include "Material_Factory.hh"
#include "XML_Node.hh"

using namespace std;

Manufactured_Solid_Geometry::
Manufactured_Solid_Geometry(shared_ptr<Angular_Discretization> angular,
                            shared_ptr<Energy_Discretization> energy,
                            shared_ptr<Manufactured_Solution> solution):
    angular_(angular),
    energy_(energy),
    solution_(solution),
    material_factory_(make_shared<Material_Factory>(angular,
                                                    energy))
{
    check_class_invariants();
}

shared_ptr<Material> Manufactured_Solid_Geometry::
material(vector<double> const &position) const
{
    // Get solution and gradient solution
    vector<double> solution = solution_->get_solution(position);
    vector<double> grad_solution = solution_->get_grad_solution(position);
    
    // Get cross sections
    vector<double> sigma_t;
    vector<double> sigma_s;
    solution_->get_cross_sections(position,
                                  sigma_t,
                                  sigma_s);

    // Get source
    vector<double> source = solution_->get_source(position,
                                                  grad_solution,
                                                  sigma_t,
                                                  sigma_s);

    // Create material
    int number_of_groups = energy_->number_of_groups();
    vector<double> nu(number_of_groups, 0);
    vector<double> sigma_f(number_of_groups, 0);
    vector<double> chi(number_of_groups, 0);
    return material_factory_->get_standard_material(0, // index
                                                    sigma_t,
                                                    sigma_s,
                                                    nu,
                                                    sigma_f,
                                                    chi,
                                                    source);
}

shared_ptr<Boundary_Source> Manufactured_Solid_Geometry::
boundary_source(vector<double> const &position) const
{
    // Get solution
    vector<double> solution = solution_->get_solution(position);

    // Convert solution from moment to discrete form
    int number_of_groups = energy_->number_of_groups();
    int number_of_moments = angular_->number_of_moments();
    int number_of_ordinates = angular_->number_of_ordinates();
    vector<double> discrete_solution(number_of_groups * number_of_ordinates);
    for (int g = 0; g < number_of_groups; ++g)
    {
        // Get moment-only form
        vector<double> local_data(number_of_moments);
        for (int m = 0; m < number_of_moments; ++m)
        {
            int k_sol = g + number_of_groups * m;
            local_data[m] = solution[k_sol];
        }

        // Convert to discrete
        angular_->moment_to_discrete(local_data);

        // Put data into discrete solution
        for (int o = 0; o < number_of_ordinates; ++o)
        {
            int k_sol = g + number_of_groups * o;
            discrete_solution[k_sol] = local_data[o];
        }
    }

    Boundary_Source::Dependencies deps;
    deps.angular = Boundary_Source::Dependencies::Angular::ORDINATES;
    vector<double> alpha(number_of_groups, 0);
    return make_shared<Boundary_Source>(0, // index
                                        deps,
                                        angular_,
                                        energy_,
                                        discrete_solution,
                                        alpha);
}

void Manufactured_Solid_Geometry::
check_class_invariants() const
{
    Assert(angular_);
    Assert(energy_);
    Assert(solution_);
    Assert(material_factory_);
}

void Manufactured_Solid_Geometry::
output(XML_Node output_node) const
{
}
