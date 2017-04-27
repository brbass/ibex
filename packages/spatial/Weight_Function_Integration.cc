#include "Weight_Function_Integration.hh"

#include <algorithm>
#include <cmath>

#include "Basis_Function.hh"
#include "Check.hh"
#include "KD_Tree.hh"
#include "Solid_Geometry.hh"
#include "Weight_Function.hh"

using namespace std;

Weight_Function_Integration::
Weight_Function_Integration(int number_of_points,
                            vector<shared_ptr<Basis_Function> > const &bases,
                            vector<shared_ptr<Weight_Function> > const &weights,
                            shared_ptr<Solid_Geometry> solid,
                            vector<vector<double> > limits,
                            vector<int> dimensional_cells):
    options_(weights[0]->options()),
    number_of_points_(number_of_points),
    bases_(bases),
    weights_(weights),
    solid_(solid)
{
    Assert(bases.size() == number_of_points_);
    Assert(weights.size() == number_of_points_);
    Assert(solid);
    
    mesh_ = make_shared<Mesh>(*this,
                              solid->dimension(),
                              limits,
                              dimensional_cells);
}

void Weight_Function_Integration::
initialize_integrals(vector<Weight_Function::Integrals> integrals) const
{
    integrals.resize(number_of_points_);
    for (int i = 0; i < number_of_points_; ++i)
    {
        shared_ptr<Weight_Function> weight = weights_[i];
        int number_of_boundary_surfaces = weight->number_of_boundary_surfaces();
        int number_of_basis_functions = weight->number_of_basis_functions();
        Weight_Function::Integrals &local_integrals = integrals[i];
        local_integrals.is_w.assign(number_of_boundary_surfaces, 0.);
        local_integrals.is_b_w.assign(number_of_boundary_surfaces * number_of_basis_functions);
        local_integrals.iv_w.assign(1, 0.);
        local_integrals.iv_dw.assign(dimension_, 0);
        local_integrals.iv_b_w.assign(number_of_basis_functions, 0);
        local_integrals.iv_b_dw.assign(number_of_basis_functions * dimension_, 0);
        local_integrals.iv_db_w.assign(number_of_basis_functions * dimension_, 0);
        local_integrals.iv_db_dw.assign(number_of_basis_functions * dimension_ * dimension_, 0);
    }
}

void Weight_Function_Integration::
initialize_materials(vector<Material_Data> materials) const
{
    // Get material data
    shared_ptr<Material> test_material = solid_geometry_->material(weights_[0]->position());
    shared_ptr<Angular_Discretization> angular_discretization = test_material->angular_discretization();
    shared_ptr<Energy_Discretization> energy_discretization = test_material->energy_discretization();
    int number_of_groups = energy_discretization->number_of_groups();
    int number_of_scattering_moments = angular_discretization->number_of_scattering_moments();
    int number_of_moments = angular_discretization->number_of_moments();

    // Initialize materials
    materials.resize(number_of_points_);
    for (int i = 0; i < number_of_points_; ++i)
    {
        Material_Data &material = materials[i];
        
        
    }
    
}

void Weight_Function_Integration::
perform_integration()
{
    // Initialize integrals to zero
    vector<Weight_Function::Integrals> integrals;
    initialize_integrals(integrals);
    
    // Initialize materials to zero
    vector<Material_Data> materials;
    initialize_materials(materials);
    
    for (int i = 0; i < mesh_->number_of_background_cells_; ++i)
    {
        
    }
}

Weight_Function_Integration::Mesh::
Mesh(Weight_Function_Integration const &wfi,
     int dimension,
     vector<vector<double> > limits,
     vector<int> dimensional_cells):
    wfi_(wfi),
    dimension_(dimension),
    limits_(limits),
    dimensional_cells_(dimensional_cells)
{
    initialize_mesh();
    initialize_connectivity();
}

void Weight_Function_Integration::Mesh::
initialize_mesh()
{
    vector<vector<double> > positions;
    
    // Check sizes
    Assert(dimensional_cells_.size() == dimension_);
    Assert(limits_.size() == dimension_);
    
    // Get total number of nodes
    dimensional_nodes_.resize(dimension_);
    number_of_background_nodes_ = 1;
    number_of_background_cells_ = 1;
    for (int d = 0; d < dimension_; ++d)
    {
        Assert(dimensional_cells_[d] >= 1);
        dimensional_nodes_[d] = dimensional_cells_[d] + 1;
        number_of_background_nodes_ *= dimensional_nodes_[d];
        number_of_background_cells_ *= dimensional_cells_[d];
    }
    
    // Get intervals between cells
    intervals_.resize(dimension_);
    for (int d = 0; d < dimension_; ++d)
    {
        intervals_[d] = (limits_[d][1] - limits_[d][0]) / static_cast<double>(dimensional_cells_[d]);
    }
    

    // Initialize nodes
    nodes_.resize(number_of_background_nodes_);
    switch (dimension_)
    {
    case 1:
    {
        int di = 0;
        for (int i = 0; i < dimensional_nodes_[0]; ++i)
        {
            double index = i;
            Node &node = nodes_[i];
            node.position = {limits_[di][0] + intervals_[di] * i};
        }
        break;
    }
    case 2:
    {
        int di = 0;
        int dj = 1;
        for (int i = 0; i < dimensional_nodes_[0]; ++i)
        {
            for (int j = 0; j < dimensional_nodes_[1]; ++i)
            {
                int index = j + dimensional_nodes_[1] * i;
                Node &node = nodes_[index];
                node.position = {limits_[di][0] + intervals_[di] * i,
                                 limits_[dj][0] + intervals_[dj] * j};
            }
        }
        break;
    }
    case 3:
    {
        int di = 0;
        int dj = 1;
        int dk = 2;
        for (int i = 0; i < dimensional_nodes_[0]; ++i)
        {
            for (int j = 0; j < dimensional_nodes_[1]; ++i)
            {
                for (int k = 0; k < dimensional_nodes_[2]; ++k)
                {
                    int index = k + dimensional_nodes_[2] * (j + dimensional_nodes_[1] * i);
                    Node &node = nodes_[index];
                    node.position = {limits_[di][0] + intervals_[di] * i,
                                     limits_[dj][0] + intervals_[dj] * j,
                                     limits_[dk][0] + intervals_[dk] * k};
                }
            }
            break;
        }
    }
    default:
        AssertMsg(false, "dimension (" + to_string(dimension_) +  ") not found");
    }

    // Initialize cells
    switch (dimension_)
    {
    case 1:
    {
        int di = 0;
        for (int i = 0; i < dimensional_cells_[di]; ++i)
        {
            int index = i;
            Cell &cell = cells_[index];
            
            // Set upper and lower limits
            cell.limits.resize(dimension_);
            for (int d = 0; d < dimension_; ++d)
            {
                int l0 = index;
                int l1 = l0 + 1;
                cell.limits[d] = {limits_[d][0] + intervals_[d] * l0,
                                  limits_[d][0] + intervals_[d] * l1};
            }
                
            // Set neighboring nodes and cells
            for (int ni = i; ni <= i + 1; ++ni)
            {
                int n_index = ni;
                Node &node = nodes_[n_index];
                node.neighboring_cells.push_back(index);
                cell.neighboring_nodes.push_back(n_index);
            }
        }
        break;
    }
    case 2:
    {
        int di = 0;
        int dj = 1;
        for (int i = 0; i < dimensional_cells_[di]; ++i)
        {
            for (int j = 0; j < dimensional_cells_[dj]; ++j)
            {
                int index = j + dimensional_cells_[dj] * i;
                vector<int> indices = {i, j};
                Cell &cell = cells_[index];
                
                // Set upper and lower limits
                cell.limits.resize(dimension_);
                for (int d = 0; d < dimension_; ++d)
                {
                    int l0 = indices[d];
                    int l1 = l0 + 1;
                    cell.limits[d] = {limits_[d][0] + intervals_[d] * l0,
                                      limits_[d][0] + intervals_[d] * l1};
                }
                
                // Set neighboring nodes and cells
                for (int ni = i; ni <= i + 1; ++ni)
                {
                    for (int nj = j; nj <= j + 1; ++nj)
                    {
                        int n_index = nj + dimensional_nodes_[dj] * ni;
                        Node &node = nodes_[n_index];
                        node.neighboring_cells.push_back(index);
                        cell.neighboring_nodes.push_back(n_index);
                    }
                }
            }
        }
        break;
    }
    case 3:
    {
        int di = 0;
        int dj = 1;
        int dk = 2;
        for (int i = 0; i < dimensional_cells_[di]; ++i)
        {
            for (int j = 0; j < dimensional_cells_[dj]; ++j)
            {
                for (int k = 0; k < dimensional_cells_[dk]; ++k)
                {
                    int index = k + dimensional_cells_[dk] * (j + dimensional_cells_[dj] * i);
                    vector<int> indices = {i, j, k};
                    Cell &cell = cells_[index];
                    
                    // Set neighboring nodes and cells
                    for (int ni = i; ni <= i + 1; ++ni)
                    {
                        for (int nj = j; nj <= j + 1; ++nj)
                        {
                            for (int nk = k; nk <= k + 1; ++nk)
                            {
                                int n_index = nk + dimensional_nodes_[dk] * (nj + dimensional_nodes_[dj] * ni);
                                Node &node = nodes_[n_index];
                                node.neighboring_cells.push_back(index);
                                cell.neighboring_nodes.push_back(n_index);
                            }
                        }
                    }
                
                    // Set upper and lower limits
                    cell.limits.resize(dimension_);
                    for (int d = 0; d < dimension_; ++d)
                    {
                        int l0 = indices[d];
                        int l1 = l0 + 1;
                        cell.limits[d] = {limits_[d][0] + intervals_[d] * l0,
                                          limits_[d][0] + intervals_[d] * l1};
                    }
                }
            }
        }
        break;
    }
    default:
        AssertMsg(false, "dimension (" + to_string(dimension_) + ") not found");
    }
    
    // Get KD tree
    vector<vector<double> > kd_positions(number_of_background_nodes_, vector<double>(dimension_));
    for (int i = 0; i < number_of_background_nodes_; ++i)
    {
        kd_positions[i] = nodes_[i].position;
    }
    kd_tree_ = make_shared<KD_Tree>(dimension_,
                                    number_of_background_nodes_,
                                    kd_positions);

    // Get maximum interval
    max_interval_ = *max_element(intervals_.begin(), intervals_.end());
}

void Weight_Function_Integration::Mesh::
initialize_connectivity()
{
    int number_of_points = wfi_.number_of_points_;

    // Get weight function connectivity
    for (int i = 0; i < number_of_points; ++i)
    {
        // Get weight function data
        shared_ptr<Weight_Function> weight = wfi_.weights_[i];
        double radius = get_inclusive_radius(weight->radius());
        vector<double> const position = weight->position();
        
        // Find nodes that intersect with the weight function
        vector<int> intersecting_nodes;
        vector<double> distances;
        int number_of_intersecting_nodes
            = kd_tree_->radius_search(radius,
                                      position,
                                      intersecting_nodes,
                                      distances);
        
        // Add cells for these nodes to the weight indices
        for (int j = 0; j < number_of_intersecting_nodes; ++j)
        {
            Node &node = nodes_[intersecting_nodes[j]];
            
            for (int c_index : node.neighboring_cells)
            {
                cells_[c_index].weight_indices.push_back(i);
            }
        }
    }
    
    // Get basis function connectivity
    for (int i = 0; i < number_of_points; ++i)
    {
        // Get basis function data
        shared_ptr<Basis_Function> basis = wfi_.bases_[i];
        double radius = get_inclusive_radius(basis->radius());
        vector<double> const position = basis->position();
        
        // Find nodes that intersect with the weight function
        vector<int> intersecting_nodes;
        vector<double> distances;
        int number_of_intersecting_nodes
            = kd_tree_->radius_search(radius,
                                      position,
                                      intersecting_nodes,
                                      distances);
        
        // Add cells for these nodes to the weight indices
        for (int j = 0; j < number_of_intersecting_nodes; ++j)
        {
            Node &node = nodes_[intersecting_nodes[j]];
            
            for (int c_index : node.neighboring_cells)
            {
                cells_[c_index].basis_indices.push_back(i);
            }
        }
    }

    // Remove duplicate basis/weight indices
    for (int i = 0; i < number_of_background_cells_; ++i)
    {
        Cell &cell = cells_[i];
        
        sort(cell.basis_indices.begin(), cell.basis_indices.end());
        cell.basis_indices.erase(unique(cell.basis_indices.begin(), cell.basis_indices.end()), cell.basis_indices.end());
        
        sort(cell.weight_indices.begin(), cell.weight_indices.end());
        cell.weight_indices.erase(unique(cell.weight_indices.begin(), cell.weight_indices.end()), cell.weight_indices.end());
    }
}

double Weight_Function_Integration::Mesh::
get_inclusive_radius(double radius) const
{
    // Move radius outward to account for edges
    switch (dimension_)
    {
    case 1:
        return radius;
    case 2:
        return sqrt(radius * radius + 0.25 * max_interval_ * max_interval_);
    case 3:
        return sqrt(radius * radius + 0.5 * max_interval_ * max_interval_);
    }
}

