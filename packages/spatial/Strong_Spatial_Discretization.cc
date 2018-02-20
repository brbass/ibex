#include "Strong_Spatial_Discretization.hh"

#include "Boundary_Source.hh"
#include "Check.hh"
#include "Material.hh"
#include "Solid_Geometry.hh"
#include "Weight_Function.hh"

using namespace std;

Strong_Spatial_Discretization::
Strong_Spatial_Discretization(vector<shared_ptr<Basis_Function> > &bases,
                              vector<shared_ptr<Weight_Function> > &weights,
                              shared_ptr<Dimensional_Moments> dimensional_moments,
                              shared_ptr<Weak_Spatial_Discretization_Options> options,
                              shared_ptr<KD_Tree> kd_tree):
    Weak_Spatial_Discretization(bases,
                                weights,
                                dimensional_moments,
                                options,
                                kd_tree)
{
    switch (options_->weighting)
    {
    case Weak_Spatial_Discretization_Options::Weighting::BASIS:
        perform_basis_integration();
        break;
    case Weak_Spatial_Discretization_Options::Weighting::POINT:
        perform_point_integration();
        break;
    default:
        AssertMsg(false, "strong discretization not compatible with given weighting");
        break;
    }

    options_->normalized = true;
}

void Strong_Spatial_Discretization::
perform_point_integration()
{
    shared_ptr<Solid_Geometry> solid = options_->solid;
    Assert(solid);
    
    for (int i = 0; i < number_of_points_; ++i)
    {
        // Get weight function information
        shared_ptr<Weight_Function> weight = weights_[i];
        vector<double> const position = weight->position();
        int const number_of_boundary_surfaces = weight->number_of_boundary_surfaces();
        int const number_of_basis_functions = weight->number_of_basis_functions();
        
        // Set dummy values for integrals
        Weight_Function::Integrals integrals;
        integrals.is_w.assign(number_of_boundary_surfaces, 0);
        integrals.is_b_w.assign(number_of_boundary_surfaces * number_of_basis_functions, 0);
        integrals.iv_w.assign(1, 0);
        integrals.iv_dw.assign(dimension_, 0);
        integrals.iv_b_w.assign(number_of_basis_functions, 0);
        integrals.iv_b_dw.assign(number_of_basis_functions * dimension_, 0);
        integrals.iv_db_w.assign(number_of_basis_functions * dimension_, 0);
        integrals.iv_db_dw.assign(number_of_basis_functions * dimension_ * dimension_, 0);
        
        // Get material and boundary source at this point
        shared_ptr<Material> material = solid->material(position);;
        vector<shared_ptr<Boundary_Source> > boundary_sources;
        if (number_of_boundary_surfaces > 0)
        {
            Assert(number_of_boundary_surfaces == 1);
            boundary_sources.push_back(solid->boundary_source(position));
        }
        
        // Put information into weight function
        weight->set_integrals(integrals,
                              material,
                              boundary_sources);
    }
}

void Strong_Spatial_Discretization::
perform_basis_integration()
{
    AssertMsg(false, "not implemented");
}
