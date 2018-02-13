#include <memory>
#include <mpi.h>
#include <string>
#include <vector>

#include "Angular_Discretization_Parser.hh"
#include "Boundary_Source_Parser.hh"
#include "Cartesian_Plane.hh"
#include "Cross_Section.hh"
#include "Energy_Discretization.hh"
#include "Energy_Discretization_Parser.hh"
#include "Krylov_Eigenvalue.hh"
#include "LDFE_Quadrature.hh"
#include "Material.hh"
#include "Material_Parser.hh"
#include "Quadrature_Rule.hh"
#include "Solver.hh"
#include "Solver_Parser.hh"
#include "Transport_Discretization.hh"
#include "VERA_Solid_Geometry.hh"
#include "Weak_RBF_Sweep.hh"
#include "Weak_Spatial_Discretization.hh"
#include "Weak_Spatial_Discretization_Parser.hh"
#include "Weak_Sweep_Parser.hh"
#include "XML_Document.hh"
#include "XML_Node.hh"

using namespace std;

class VERA_Result
{
public:
    
    VERA_Result(shared_ptr<Solid_Geometry> solid,
                shared_ptr<Angular_Discretization> angular,
                shared_ptr<Energy_Discretization> energy,
                shared_ptr<Weak_Spatial_Discretization> spatial,
                shared_ptr<Solver::Result> result):
        solid_(solid),
        angular_(angular),
        energy_(energy),
        spatial_(spatial),
        result_(result)
    {
        number_of_ordinates_ = 4;
        Quadrature_Rule::cartesian_1d(Quadrature_Rule::Quadrature_Type::GAUSS_LEGENDRE,
                                      number_of_ordinates_,
                                      0,
                                      M_PI / 4,
                                      ordinates_,
                                      weights_);
    }

    double get_radial_fission_energy(double radius)
    {
        // Check whether position is inside the fuel
        if (radius > 0.4096)
        {
            return 0;
        }
        
        // Get size information
        int number_of_groups = energy_->number_of_groups();
        int number_of_moments = angular_->number_of_moments();

        // Get angle-independent material at this radius
        double const mev_to_joule = 1.6021766e-13;
        shared_ptr<Material> const material
            = solid_->material({radius, 0});
        vector<double> const chi_nu_sigma_f
            = material->sigma_f()->data();
        vector<double> const nu = {2.66457, 2.44118};
        vector<double> const kappa = {195.8, 193.4};
        vector<double> kappa_sigma_f(number_of_groups, 0);
        for (int gf = 0; gf < number_of_groups; ++gf)
        {
            for (int gt = 0; gt < number_of_groups; ++gt)
            {
                kappa_sigma_f[gf] += chi_nu_sigma_f[gf + number_of_groups * gt];
            }
            
            kappa_sigma_f[gf] *= kappa[gf] / nu[gf] * mev_to_joule;
        }
        
        // Integrate source radially from 0 to pi/4
        double source = 0;
        for (int q = 0; q < number_of_ordinates_; ++q)
        {
            // Get position
            double const theta = ordinates_[q];
            vector<double> const position
                = {radius * cos(theta),
                   radius * sin(theta)};
            
            // Get flux values
            vector<double> const flux
                = spatial_->expansion_values(number_of_groups * number_of_moments,
                                             position,
                                             result_->coefficients);
            
            // Get fission source values
            int const m = 0;
            for (int g = 0; g < number_of_groups; ++g)
            {
                int const k_flux = g + number_of_groups * m;
                source += flux[k_flux] * kappa_sigma_f[g] * weights_[q];
            }
        }
        
        // Normalize integral
        source *= 4 / M_PI;
        
        return source;
    }

private:

    shared_ptr<Solid_Geometry> solid_;
    shared_ptr<Angular_Discretization> angular_;
    shared_ptr<Energy_Discretization> energy_;
    shared_ptr<Weak_Spatial_Discretization> spatial_;
    shared_ptr<Solver::Result> result_;
    int number_of_ordinates_;
    vector<double> ordinates_;
    vector<double> weights_;
};

shared_ptr<VERA_Result>
run_transport(string filename,
              shared_ptr<VERA_Temperature> temperature)
{
    XML_Document input_file(filename);
    XML_Node input_node = input_file.get_child("input");
    
    // Get energy discretization
    Energy_Discretization_Parser energy_parser;
    shared_ptr<Energy_Discretization> energy
        = energy_parser.parse_from_xml(input_node.get_child("energy_discretization"));

    // Get angular discretization
    Angular_Discretization_Parser angular_parser;
    shared_ptr<Angular_Discretization> angular
        = angular_parser.parse_from_xml(input_node.get_child("angular_discretization"));

    // Get materials
    Material_Parser material_parser(angular,
                                    energy);
    bool include_ifba = input_node.get_child("materials").get_attribute<bool>("include_ifba");
    vector<shared_ptr<Material> > materials
        = material_parser.parse_from_xml(input_node.get_child("materials"));
    
    // Get boundary source
    Boundary_Source_Parser boundary_parser(angular,
                                           energy);
    vector<shared_ptr<Boundary_Source> > boundary_sources
        = boundary_parser.parse_from_xml(input_node.get_child("boundary_sources"));
    Assert(boundary_sources.size() == 1);
    
    // Get solid geometry
    shared_ptr<VERA_Solid_Geometry> solid
        = make_shared<VERA_Solid_Geometry>(include_ifba,
                                           temperature,
                                           angular,
                                           energy,
                                           materials,
                                           boundary_sources[0]);

    // Get boundary surfaces
    vector<shared_ptr<Cartesian_Plane> > boundary_surfaces
        = solid->cartesian_boundary_surfaces();
    
    // Get spatial discretization
    Weak_Spatial_Discretization_Parser spatial_parser(solid,
                                                      boundary_surfaces);
    shared_ptr<Weak_Spatial_Discretization> spatial
        = spatial_parser.get_weak_discretization(input_node.get_child("spatial_discretization"));
    
    // Get transport discretization
    shared_ptr<Transport_Discretization> transport
        = make_shared<Transport_Discretization>(spatial,
                                                angular,
                                                energy);

    // Get sweep
    Weak_Sweep_Parser sweep_parser(spatial,
                                   angular,
                                   energy,
                                   transport);
    shared_ptr<Weak_RBF_Sweep> sweep
        = sweep_parser.get_weak_rbf_sweep(input_node.get_child("transport"));
    
    // Get solver
    Solver_Parser solver_parser(spatial,
                                angular,
                                energy,
                                transport);
    shared_ptr<Solver> solver
        = solver_parser.get_krylov_eigenvalue(input_node.get_child("solver"),
                                              sweep);
    solver->solve();
    
    return make_shared<VERA_Result>(solid,
                                    angular,
                                    energy,
                                    spatial,
                                    solver->result());
}

shared_ptr<VERA_Temperature> 
run_heat(string filename,
         shared_ptr<VERA_Result> result)
{
    
}

int main(int argc, char **argv)
{
    // Initialize MPI
    MPI_Init(&argc, &argv);

    if (argc != 2)
    {
        cerr << "need input file" << endl;
        return 1;
    }

    string filename = argv[1];

    shared_ptr<VERA_Temperature> temperature
        =                         make_shared<VERA_Temperature>([](vector<double> const &){return 600;});
    shared_ptr<VERA_Result> result
        = run_transport(filename,
                        temperature);
    

    
    // Close MPI
    MPI_Finalize();

    return 0;
}
