#include "Weak_Spatial_Discretization.hh"

#include "Basis_Function.hh"
#include "Boundary_Source.hh"
#include "Cartesian_Plane.hh"
#include "Check.hh"
#include "Conversion.hh"
#include "Dimensional_Moments.hh"
#include "KD_Tree.hh"
#include "Meshless_Function.hh"
#include "Meshless_Normalization.hh"
#include "Weight_Function_Integration.hh"
#include "XML_Node.hh"

using namespace std;

Weak_Spatial_Discretization::
Weak_Spatial_Discretization(vector<shared_ptr<Basis_Function> > &bases,
                            vector<shared_ptr<Weight_Function> > &weights,
                            shared_ptr<Dimensional_Moments> dimensional_moments,
                            shared_ptr<Weak_Spatial_Discretization_Options> options,
                            shared_ptr<KD_Tree> kd_tree):
    basis_depends_on_neighbors_(bases[0]->function()->depends_on_neighbors()),
    weight_depends_on_neighbors_(weights[0]->function()->depends_on_neighbors()),
    number_of_points_(weights.size()),
    dimension_(weights[0]->dimension()),
    number_of_nodes_(weights[0]->number_of_nodes()),
    weights_(weights),
    bases_(bases),
    options_(options),
    dimensional_moments_(dimensional_moments),
    kd_tree_(kd_tree)
{
    // Check whether weak parameter options are finalized
    if (!options_->input_finalized)
    {
        options_->finalize_input();
    }

    // Get number of basis functions
    number_of_basis_functions_.resize(number_of_points_);
    for (int i = 0; i < number_of_points_; ++i)
    {
        number_of_basis_functions_[i] = weights_[i]->number_of_basis_functions();
    }

    // Get boundary weights
    int j = 0;
    boundary_weights_.resize(number_of_points_);
    for (int i = 0; i < number_of_points_; ++i)
    {
        if (weights_[i]->point_type() == Weight_Function::Point_Type::BOUNDARY)
        {
            boundary_weights_[j] = weights_[i];
            j += 1;
        }
    }
    number_of_boundary_weights_ = j;
    boundary_weights_.resize(number_of_boundary_weights_);

    // Get boundary bases
    j = 0;
    boundary_bases_.resize(number_of_points_);
    has_reflection_ = false;
    for (int i = 0; i < number_of_points_; ++i)
    {
        int number_of_boundary_surfaces = bases_[i]->number_of_boundary_surfaces();
        if (number_of_boundary_surfaces > 0)
        {
            boundary_bases_[j] = bases[i];
            bases[i]->set_boundary_index(j);
            j += 1;

            for (int s = 0; s < number_of_boundary_surfaces; ++s)
            {
                if (bases[i]->boundary_surface(s)->boundary_source()->has_reflection())
                {
                    has_reflection_ = true;
                }
            }
        }
    }
    number_of_boundary_bases_ = j;
    boundary_bases_.resize(number_of_boundary_bases_);

    // Get KD tree
    if (!kd_tree_)
    {
        vector<vector<double> > points(number_of_points_);
        for (int i = 0; i < number_of_points_; ++i)
        {
            vector<double> const point = weights_[i]->position();
        
            points[i] = point;
        }
        kd_tree_ = make_shared<KD_Tree>(dimension_,
                                        number_of_points_,
                                        points);
    }

    // Perform external integrals, if applicable
    if (options_->external_integral_calculation)
    {
        Weight_Function_Integration integrator(number_of_points_,
                                               options_,
                                               bases,
                                               weights);
        integrator.perform_integration();
    }
    
    check_class_invariants();
}

int Weak_Spatial_Discretization::
nearest_point(vector<double> const &position) const
{
    vector<int> index;
    vector<double> squared_distance;

    kd_tree_->find_neighbors(1,
                             position,
                             index,
                             squared_distance);
    
    return index[0];
}

double Weak_Spatial_Discretization::
collocation_value(int i,
                  vector<double> const &coefficients) const
{
    shared_ptr<Weight_Function> weight = weights_[i];
    int number_of_basis_functions = weight->number_of_basis_functions();
    vector<int> const basis_indices = weight->basis_function_indices();
    vector<double> const &v_b = weight->values().v_b;

    // Sum over coefficients
    double sum = 0;
    for (int j = 0; j < number_of_basis_functions; ++j)
    {
        int index = basis_indices[j];
        sum += v_b[j] * coefficients[index];
    }
    return sum;
}

double Weak_Spatial_Discretization::
weighted_collocation_value(int i,
                           vector<double> const &coefficients) const
{
    shared_ptr<Weight_Function> weight = weights_[i];
    Weight_Function::Integrals const integrals = weight->integrals();
    int number_of_basis_functions = weight->number_of_basis_functions();
    vector<int> const basis_indices = weight->basis_function_indices();
    vector<double> const &iv_b_w = integrals.iv_b_w;
    double const iv_w = integrals.iv_w[0];
    
    // Sum over coefficients
    double sum = 0;
    for (int j = 0; j < number_of_basis_functions; ++j)
    {
        int index = basis_indices[j];
        sum += iv_b_w[j] * coefficients[index];
    }
    return sum / iv_w;
}

void Weak_Spatial_Discretization::
collocation_values(vector<double> const &coefficients,
                   vector<double> &values) const
{
    // Get values at each point
    values.resize(number_of_points_);
    for (int i = 0; i < number_of_points_; ++i)
    {
        values[i] = collocation_value(i,
                                      coefficients);
    }
}

void Weak_Spatial_Discretization::
weighted_collocation_values(vector<double> const &coefficients,
                            vector<double> &values) const
{
    // Get weighted values at each point
    values.resize(number_of_points_);
    for (int i = 0; i < number_of_points_; ++i)
    {
        values[i] = weighted_collocation_value(i,
                                               coefficients);
    }
}

double Weak_Spatial_Discretization::
expansion_value(int i,
                vector<double> const &position,
                vector<double> const &coefficients) const
{
    Assert(coefficients.size() == number_of_points_);
    
    shared_ptr<Weight_Function> weight = weights_[i];
    int number_of_basis_functions = weight->number_of_basis_functions();
    vector<int> const basis_indices = weight->basis_function_indices();

    // Get value of function at a specific point
    double sum = 0;
    for (int j = 0; j < number_of_basis_functions; ++j)
    {
        int index = basis_indices[j];
        double val = weight->basis_function(j)->function()->value(position);
        sum += val * coefficients[index];
    }
    return sum;
}

double Weak_Spatial_Discretization::
expansion_value(vector<double> const &position,
                vector<double> const &coefficients) const
{
    int index = nearest_point(position);
    
    return expansion_value(index,
                           position,
                           coefficients);
}

vector<double> Weak_Spatial_Discretization::
expansion_values(int i,
                 int number_of_groups,
                 vector<double> const &position,
                 vector<double> const &coefficients) const
{
    Assert(coefficients.size() == number_of_points_ * number_of_groups);
    
    shared_ptr<Weight_Function> weight = weights_[i];
    int number_of_basis_functions = weight->number_of_basis_functions();
    vector<int> const basis_indices = weight->basis_function_indices();

    // Get the basis function values at the position
    vector<double> basis_vals(number_of_basis_functions);
    for (int j = 0; j < number_of_basis_functions; ++j)
    {
        basis_vals[j] = weight->basis_function(j)->function()->base_function()->value(position);
    }

    // Normalize if applicable
    if (basis_depends_on_neighbors_)
    {
        vector<vector<double> > center_positions(number_of_basis_functions);
        for (int j = 0; j < number_of_basis_functions; ++j)
        {
            center_positions[j] = weight->basis_function(j)->position();
        }
        shared_ptr<Meshless_Normalization> norm
            = weight->basis_function(0)->function()->normalization();
        norm->get_values(position,
                         center_positions,
                         basis_vals,
                         basis_vals);
    }
    
    // Get value of function at a specific point
    vector<double> vals(number_of_groups, 0.);
    for (int j = 0; j < number_of_basis_functions; ++j)
    {
        int basis_index = basis_indices[j];
        double val = basis_vals[j];
        for (int g = 0; g < number_of_groups; ++g)
        {
            int index = g + number_of_groups * basis_index;
            vals[g] += val * coefficients[index];
        }
    }
    
    return vals;
}

vector<double> Weak_Spatial_Discretization::
expansion_values(int number_of_groups,
                 vector<double> const &position,
                 vector<double> const &coefficients) const
{
    int index = nearest_point(position);
    
    return expansion_values(index,
                            number_of_groups,
                            position,
                            coefficients);
}

void Weak_Spatial_Discretization::
check_class_invariants() const
{
    Assert(weights_.size() == number_of_points_);
    Assert(boundary_weights_.size() == number_of_boundary_weights_);
    Assert(bases_.size() == number_of_points_);
    Assert(boundary_bases_.size() == number_of_boundary_bases_);
    
    for (int i = 0; i < number_of_points_; ++i)
    {
        Assert(weights_[i]);
        Assert(bases_[i]);
    }
    for (int i = 0; i < number_of_boundary_weights_; ++i)
    {
        Assert(boundary_weights_[i]);
    }
    for (int i = 0; i < number_of_boundary_bases_; ++i)
    {
        Assert(boundary_bases_[i]);
    }
}

void Weak_Spatial_Discretization::
output(XML_Node output_node) const
{
    // Output scalar values
    output_node.set_child_value(has_reflection_, "has_reflection");
    output_node.set_child_value(number_of_points_, "number_of_points");
    output_node.set_child_value(number_of_boundary_weights_, "number_of_boundary_weights");
    output_node.set_child_value(number_of_boundary_bases_, "number_of_boundary_bases");
    output_node.set_child_value(dimension_, "dimension");
    output_node.set_child_value(number_of_nodes_, "number_of_nodes");

    // Output options
    options_->output(output_node.append_child("options"));
    
    // Output point positions
    vector<double> points(dimension_ * number_of_points_);
    for (int i = 0; i < number_of_points_; ++i)
    {
        vector<double> const point = weights_[i]->position();
        for (int d = 0; d < dimension_; ++d)
        {
            points[d + dimension_ * i] = point[d];
        }
    }
    output_node.set_child_vector(points, "points", "dimension-point");

    // Output weights
    XML_Node weights_node = output_node.append_child("weights");
    for (int i = 0; i < number_of_points_; ++i)
    {
        weights_[i]->output(weights_node.append_child("weight"));
    }

    // Output bases
    XML_Node bases_node = output_node.append_child("bases");
    for (int i = 0; i < number_of_points_; ++i)
    {
        bases_[i]->output(bases_node.append_child("basis"));
    }
}

void Weak_Spatial_Discretization_Options::
output(XML_Node output_node) const
{
    output_node.set_child_value(external_integral_calculation, "external_integral_calculation");
    output_node.set_child_value(integration_ordinates, "integration_ordinates");
    output_node.set_child_value(include_supg, "include_supg");
    output_node.set_child_value(identical_basis_functions_conversion()->convert(identical_basis_functions), "identical_basis_functions");
    output_node.set_child_value(weighting_conversion()->convert(weighting), "weighting");
    output_node.set_child_value(tau_scaling_conversion()->convert(tau_scaling), "tau_scaling");
    output_node.set_child_value(total_conversion()->convert(total), "total");
}

shared_ptr<Conversion<Weak_Spatial_Discretization_Options::Weighting, string> > Weak_Spatial_Discretization_Options::
weighting_conversion() const
{
    vector<pair<Weighting, string> > conversions
        = {{Weighting::POINT, "point"},
           {Weighting::FLAT, "flat"},
           {Weighting::FLUX, "flux"},
           {Weighting::FULL, "full"},
           {Weighting::BASIS, "basis"}};
    return make_shared<Conversion<Weighting, string> >(conversions);
}

shared_ptr<Conversion<Weak_Spatial_Discretization_Options::Total, string> > Weak_Spatial_Discretization_Options::
total_conversion() const
{
    vector<pair<Total, string> > conversions
        = {{Total::ISOTROPIC, "isotropic"},
           {Total::MOMENT, "moment"}};
    return make_shared<Conversion<Total, string> >(conversions);
}

shared_ptr<Conversion<Weak_Spatial_Discretization_Options::Tau_Scaling, string> > Weak_Spatial_Discretization_Options::
tau_scaling_conversion() const
{
    vector<pair<Tau_Scaling, string> > conversions
        = {{Tau_Scaling::NONE, "none"},
           {Tau_Scaling::FUNCTIONAL, "functional"},
           {Tau_Scaling::LINEAR, "linear"},
           {Tau_Scaling::ABSOLUTE, "absolute"},
           {Tau_Scaling::CONSTANT, "constant"}};
    return make_shared<Conversion<Tau_Scaling, string> >(conversions);
}

shared_ptr<Conversion<Weak_Spatial_Discretization_Options::Identical_Basis_Functions, string> > Weak_Spatial_Discretization_Options::
identical_basis_functions_conversion() const
{
    vector<pair<Identical_Basis_Functions, string> > conversions
        = {{Identical_Basis_Functions::AUTO, "auto"},
           {Identical_Basis_Functions::TRUE, "true"},
           {Identical_Basis_Functions::FALSE, "false"}};
    return make_shared<Conversion<Identical_Basis_Functions, string> >(conversions);
}

void Weak_Spatial_Discretization_Options::
finalize_input()
{
    if (include_supg)
    {
        switch (weighting)
        {
        case Weighting::FULL:
            // Fallthrough intentional
        case Weighting::BASIS:
            normalized = true;
            break;
        default:
            normalized = false;
            break;
        }
    }
    else
    {
        normalized = true;
    }
    
    input_finalized = true;
}
