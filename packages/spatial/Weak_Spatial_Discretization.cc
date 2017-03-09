#include "Weak_Spatial_Discretization.hh"

#include "Basis_Function.hh"
#include "Boundary_Source.hh"
#include "Cartesian_Plane.hh"
#include "Check.hh"
#include "KD_Tree.hh"
#include "Meshless_Function.hh"
#include "XML_Node.hh"

using namespace std;

Weak_Spatial_Discretization::
Weak_Spatial_Discretization(vector<shared_ptr<Basis_Function> > &bases,
                            vector<shared_ptr<Weight_Function> > &weights):
    number_of_points_(weights.size()),
    dimension_(weights[0]->dimension()),
    number_of_nodes_(weights[0]->number_of_nodes()),
    weights_(weights),
    bases_(bases)
{
    // Get number of basis functions
    number_of_basis_functions_.resize(number_of_points_);
    for (int i = 0; i < number_of_points_; ++i)
    {
        number_of_basis_functions_[i] = weights_[i]->number_of_basis_functions();
    }

    // Get number of dimensional moments
    number_of_dimensional_moments_ = weights_[0]->number_of_dimensional_moments();
    
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
    vector<vector<double> > points(number_of_points_);
    for (int i = 0; i < number_of_points_; ++i)
    {
        vector<double> const point = weights_[i]->position();
        
        points[i] = point;
    }
    kd_tree_ = make_shared<KD_Tree>(dimension_,
                                    number_of_points_,
                                    points);
    
    check_class_invariants();
}

int Weak_Spatial_Discretization::
nearest_point(vector<double> const &position) const
{
    vector<int> index;
    vector<double> distance;

    kd_tree_->find_neighbors(1,
                             position,
                             index,
                             distance);

    return index[0];
}

double Weak_Spatial_Discretization::
collocation_value(int i,
                  vector<double> const &coefficients) const
{
    shared_ptr<Weight_Function> weight = weights_[i];
    int number_of_basis_functions = weight->number_of_basis_functions();
    vector<int> const basis_indices = weight->basis_function_indices();
    vector<double> const v_b = weight->v_b();

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
    int number_of_basis_functions = weight->number_of_basis_functions();
    vector<int> const basis_indices = weight->basis_function_indices();
    vector<double> const iv_b_w = weight->iv_b_w();
    double iv_w = weight->iv_w()[0];
    
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
        Assert(weights_[i]->number_of_dimensional_moments() == number_of_dimensional_moments_);
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
    output_node.set_child_value(number_of_points_, "number_of_points");
    output_node.set_child_value(number_of_boundary_weights_, "number_of_boundary_weights");
    output_node.set_child_value(number_of_boundary_bases_, "number_of_boundary_bases");
    output_node.set_child_value(dimension_, "dimension");
    output_node.set_child_value(number_of_nodes_, "number_of_nodes");
    
    XML_Node weights_node = output_node.append_child("weights");
    
    for (int i = 0; i < number_of_points_; ++i)
    {
        weights_[i]->output(weights_node.append_child("weight"));
    }
    for (int i = 0; i < number_of_points_; ++i)
    {
        bases_[i]->output(weights_node.append_child("basis"));
    }
}
