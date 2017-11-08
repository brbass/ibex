#include "Integration_Mesh.hh"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <limits>

#include "Basis_Function.hh"
#include "Check.hh"
#include "KD_Tree.hh"
#include "Meshless_Function.hh"
#include "Meshless_Normalization.hh"
#include "Quadrature_Rule.hh"
#include "Weak_Spatial_Discretization.hh"
#include "Weight_Function.hh"

using namespace std;

Integration_Mesh::
Integration_Mesh(int dimension,
                 int number_of_points,
                 shared_ptr<Options> options,
                 vector<shared_ptr<Basis_Function> > const &bases,
                 vector<shared_ptr<Weight_Function> > const &weights):
    dimension_(dimension),
    number_of_points_(number_of_points),
    options_(options),
    bases_(bases),
    weights_(weights)
{
    // Initialize mesh
    initialize_mesh();
    initialize_connectivity();

    // Check to ensure sufficient number of cells
    if (number_of_nodes_ < number_of_points_
        && options_->adaptive_quadrature == false)
    {
        cerr << "more background integration cells or adaptive quadrature recommended" << endl;
        cerr << "number of points:\t" << number_of_points_ << endl;
        cerr << "number of cells:\t" << number_of_cells_ << endl;
    }
    
    // Get normalization information
    apply_basis_normalization_ = bases[0]->function()->depends_on_neighbors();
    apply_weight_normalization_ = weights[0]->function()->depends_on_neighbors();

    if (apply_basis_normalization_)
    {
        basis_normalization_ = bases[0]->function()->normalization();
    }
    if (apply_weight_normalization_)
    {
        weight_normalization_ = weights[0]->function()->normalization();
    }
    
    
}

void Integration_Mesh::
initialize_mesh()
{
    vector<vector<double> > positions;
    
    // Check sizes
    Assert(options_->dimensional_cells.size() == dimension_);
    Assert(options_->limits.size() == dimension_);
    
    // Get total number of nodes
    dimensional_nodes_.resize(dimension_);
    number_of_nodes_ = 1;
    number_of_cells_ = 1;
    for (int d = 0; d < dimension_; ++d)
    {
        Assert(options_->dimensional_cells[d] >= 1);
        dimensional_nodes_[d] = options_->dimensional_cells[d] + 1;
        number_of_nodes_ *= dimensional_nodes_[d];
        number_of_cells_ *= options_->dimensional_cells[d];
    }
    
    // Get intervals between cells
    intervals_.resize(dimension_);
    for (int d = 0; d < dimension_; ++d)
    {
        intervals_[d] = (options_->limits[d][1] - options_->limits[d][0]) / static_cast<double>(options_->dimensional_cells[d]);
    }
    
    // Initialize nodes
    nodes_.resize(number_of_nodes_);
    for (shared_ptr<Node> &node : nodes_)
    {
        node = make_shared<Node>();
    }
    switch (dimension_)
    {
    case 1:
    {
        int const di = 0;
        for (int i = 0; i < dimensional_nodes_[0]; ++i)
        {
            double index = i;
            shared_ptr<Node> node = nodes_[i];
            node->position = {options_->limits[di][0] + intervals_[di] * i};
        }
        break;
    }
    case 2:
    {
        int const di = 0;
        int const dj = 1;
        for (int i = 0; i < dimensional_nodes_[0]; ++i)
        {
            for (int j = 0; j < dimensional_nodes_[1]; ++j)
            {
                int index = j + dimensional_nodes_[1] * i;
                shared_ptr<Node> node = nodes_[index];
                node->position = {options_->limits[di][0] + intervals_[di] * i,
                                  options_->limits[dj][0] + intervals_[dj] * j};
            }
        }
        break;
    }
    case 3:
    {
        int const di = 0;
        int const dj = 1;
        int const dk = 2;
        for (int i = 0; i < dimensional_nodes_[0]; ++i)
        {
            for (int j = 0; j < dimensional_nodes_[1]; ++j)
            {
                for (int k = 0; k < dimensional_nodes_[2]; ++k)
                {
                    int index = k + dimensional_nodes_[2] * (j + dimensional_nodes_[1] * i);
                    shared_ptr<Node> node = nodes_[index];
                    node->position = {options_->limits[di][0] + intervals_[di] * i,
                                      options_->limits[dj][0] + intervals_[dj] * j,
                                      options_->limits[dk][0] + intervals_[dk] * k};
                }
            }
        }
        break;
    }
    default:
        AssertMsg(false, "dimension (" + to_string(dimension_) +  ") not found");
    }

    // Initialize cells
    cells_.resize(number_of_cells_);
    for (shared_ptr<Cell> &cell : cells_)
    {
        cell = make_shared<Cell>();
    }
    switch (dimension_)
    {
    case 1:
    {
        int di = 0;
        for (int i = 0; i < options_->dimensional_cells[di]; ++i)
        {
            int index = i;
            shared_ptr<Cell> cell = cells_[index];
            
            // Set upper and lower limits
            cell->limits.resize(dimension_);
            for (int d = 0; d < dimension_; ++d)
            {
                int l0 = index;
                int l1 = l0 + 1;
                cell->limits[d] = {options_->limits[d][0] + intervals_[d] * l0,
                                   options_->limits[d][0] + intervals_[d] * l1};
            }
                
            // Set neighboring nodes and cells
            for (int ni = i; ni <= i + 1; ++ni)
            {
                int n_index = ni;
                shared_ptr<Node> node = nodes_[n_index];
                node->neighboring_cells.push_back(index);
                cell->neighboring_nodes.push_back(n_index);
            }
        }
        break;
    }
    case 2:
    {
        int di = 0;
        int dj = 1;
        for (int i = 0; i < options_->dimensional_cells[di]; ++i)
        {
            for (int j = 0; j < options_->dimensional_cells[dj]; ++j)
            {
                int index = j + options_->dimensional_cells[dj] * i;
                vector<int> indices = {i, j};
                shared_ptr<Cell> cell = cells_[index];
                
                // Set upper and lower limits
                cell->limits.resize(dimension_);
                for (int d = 0; d < dimension_; ++d)
                {
                    int l0 = indices[d];
                    int l1 = l0 + 1;
                    cell->limits[d] = {options_->limits[d][0] + intervals_[d] * l0,
                                       options_->limits[d][0] + intervals_[d] * l1};
                }
                
                // Set neighboring nodes and cells
                for (int ni = i; ni <= i + 1; ++ni)
                {
                    for (int nj = j; nj <= j + 1; ++nj)
                    {
                        int n_index = nj + dimensional_nodes_[dj] * ni;
                        shared_ptr<Node> node = nodes_[n_index];
                        node->neighboring_cells.push_back(index);
                        cell->neighboring_nodes.push_back(n_index);
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
        for (int i = 0; i < options_->dimensional_cells[di]; ++i)
        {
            for (int j = 0; j < options_->dimensional_cells[dj]; ++j)
            {
                for (int k = 0; k < options_->dimensional_cells[dk]; ++k)
                {
                    int index = k + options_->dimensional_cells[dk] * (j + options_->dimensional_cells[dj] * i);
                    vector<int> indices = {i, j, k};
                    shared_ptr<Cell> cell = cells_[index];
                    
                    // Set neighboring nodes and cells
                    for (int ni = i; ni <= i + 1; ++ni)
                    {
                        for (int nj = j; nj <= j + 1; ++nj)
                        {
                            for (int nk = k; nk <= k + 1; ++nk)
                            {
                                int n_index = nk + dimensional_nodes_[dk] * (nj + dimensional_nodes_[dj] * ni);
                                shared_ptr<Node> node = nodes_[n_index];
                                node->neighboring_cells.push_back(index);
                                cell->neighboring_nodes.push_back(n_index);
                            }
                        }
                    }
                
                    // Set upper and lower limits
                    cell->limits.resize(dimension_);
                    for (int d = 0; d < dimension_; ++d)
                    {
                        int l0 = indices[d];
                        int l1 = l0 + 1;
                        cell->limits[d] = {options_->limits[d][0] + intervals_[d] * l0,
                                           options_->limits[d][0] + intervals_[d] * l1};
                    }
                }
            }
        }
        break;
    }
    default:
        AssertMsg(false, "dimension (" + to_string(dimension_) + ") not found");
    }

    // Initialize boundary surfaces
    switch (dimension_)
    {
    case 1:
    {
        number_of_surfaces_ = 2;
        surfaces_.resize(number_of_surfaces_);
        for (shared_ptr<Surface> &surface : surfaces_)
        {
            surface = make_shared<Surface>();
        }
        for (int i = 0; i < number_of_surfaces_; ++i)
        {
            surfaces_[i]->dimension = 0;
        }
        surfaces_[0]->normal = -1;
        surfaces_[1]->normal = 1;
        surfaces_[0]->neighboring_cell = 0;
        surfaces_[1]->neighboring_cell = number_of_cells_ - 1;
        nodes_[0]->neighboring_surfaces.push_back(0);
        nodes_[dimensional_nodes_[0] - 1]->neighboring_surfaces.push_back(1);
        break;
    }
    case 2:
    {
        int const di = 0;
        int const dj = 1;
        number_of_surfaces_ = 2 * (options_->dimensional_cells[di]
                                   + options_->dimensional_cells[dj]);
        surfaces_.clear();
        surfaces_.reserve(number_of_surfaces_);
        int index = 0;
        
        // Positive and negative (p = 0, 1) x boundaries
        for (int p = 0; p < 2; ++p)
        {
            int i = p == 0 ? 0 : options_->dimensional_cells[di] - 1;
            int ni = p == 0 ? 0 : dimensional_nodes_[di] - 1;
            double normal = p == 0 ? -1 : 1;
            for (int j = 0; j < options_->dimensional_cells[dj]; ++j)
            {
                shared_ptr<Surface> surface = make_shared<Surface>();
                surface->dimension = di;
                surface->normal = normal;
                surface->neighboring_cell = j + options_->dimensional_cells[dj] * i;
                surfaces_.push_back(surface);
                
                for (int nj = j; nj <= j + 1; ++nj)
                {
                    int n_index = nj + dimensional_nodes_[dj] * ni;
                    shared_ptr<Node> node = nodes_[n_index];
                    node->neighboring_surfaces.push_back(index);
                }

                index += 1;
            }
        }
        
        // Positive and negative (p = 0, 1) y boundaries
        for (int p = 0; p < 2; ++p)
        {
            int j = p == 0 ? 0 : options_->dimensional_cells[dj] - 1;
            int nj = p == 0 ? 0 : dimensional_nodes_[dj] - 1;
            double normal = p == 0 ? -1 : 1;
            for (int i = 0; i < options_->dimensional_cells[di]; ++i)
            {
                shared_ptr<Surface> surface = make_shared<Surface>();
                surface->dimension = dj;
                surface->normal = normal;
                surface->neighboring_cell = j + options_->dimensional_cells[dj] * i;
                surfaces_.push_back(surface);
                
                for (int ni = i; ni <= i + 1; ++ni)
                {
                    int n_index = nj + dimensional_nodes_[dj] * ni;
                    shared_ptr<Node> node = nodes_[n_index];
                    node->neighboring_surfaces.push_back(index);
                }
                
                index += 1;
            }
        }
        break;
    }
    case 3:
    {
        int const di = 0;
        int const dj = 1;
        int const dk = 2;
        number_of_surfaces_ = 2 * (options_->dimensional_cells[0] * options_->dimensional_cells[1]
                                   + options_->dimensional_cells[0] * options_->dimensional_cells[2]
                                   + options_->dimensional_cells[1] * options_->dimensional_cells[2]);
        surfaces_.clear();
        surfaces_.reserve(number_of_surfaces_);
        int index = 0;
        
        // Positive and negative (p = 0, 1) x boundaries
        for (int p = 0; p < 2; ++p)
        {
            int i = p == 0 ? 0 : options_->dimensional_cells[di] - 1;
            int ni = p == 0 ? 0 : dimensional_nodes_[di] - 1;
            double normal = p == 0 ? -1 : 1;
            for (int j = 0; j < options_->dimensional_cells[dj]; ++j)
            {
                for (int k = 0; k < options_->dimensional_cells[dk]; ++k)
                {
                    shared_ptr<Surface> surface = make_shared<Surface>();
                    surface->dimension = di;
                    surface->normal = normal;
                    surface->neighboring_cell = k + options_->dimensional_cells[dk] * (j + options_->dimensional_cells[dj]* i);
                    surfaces_.push_back(surface);
                    
                    for (int nj = j; nj <= j + 1; ++nj)
                    {
                        for (int nk = k; nk <= k + 1; ++nk)
                        {
                            int n_index = nk + dimensional_nodes_[dk] * (nj + dimensional_nodes_[dj] * ni);
                            shared_ptr<Node> node = nodes_[n_index];
                            node->neighboring_surfaces.push_back(index);
                        }
                    }

                    index += 1;
                }
            }
        }

        // Positive and negative (p = 0, 1) y boundaries
        for (int p = 0; p < 2; ++p)
        {
            int j = p == 0 ? 0 : options_->dimensional_cells[dj] - 1;
            int nj = p == 0 ? 0 : dimensional_nodes_[dj] - 1;
            double normal = p == 0 ? -1 : 1;
            for (int i = 0; i < options_->dimensional_cells[di]; ++i)
            {
                for (int k = 0; k < options_->dimensional_cells[dk]; ++k)
                {
                    shared_ptr<Surface> surface = make_shared<Surface>();
                    surface->dimension = dj;
                    surface->normal = normal;
                    surface->neighboring_cell = k + options_->dimensional_cells[dk] * (j + options_->dimensional_cells[dj]* i);
                    surfaces_.push_back(surface);
                    
                    for (int ni = i; ni <= i + 1; ++ni)
                    {
                        for (int nk = k; nk <= k + 1; ++nk)
                        {
                            int n_index = nk + dimensional_nodes_[dk] * (nj + dimensional_nodes_[dj] * ni);
                            shared_ptr<Node> node = nodes_[n_index];
                            node->neighboring_surfaces.push_back(index);
                        }
                    }

                    index += 1;
                }
            }
        }

        // Positive and negative (p = 0, 1) z boundaries
        for (int p = 0; p < 2; ++p)
        {
            int k = p == 0 ? 0 : options_->dimensional_cells[dk] - 1;
            int nk = p == 0 ? 0 : dimensional_nodes_[dk] - 1;
            double normal = p == 0 ? -1 : 1;
            for (int i = 0; i < options_->dimensional_cells[di]; ++i)
            {
                for (int j = 0; j < options_->dimensional_cells[dj]; ++j)
                {
                    shared_ptr<Surface> surface = make_shared<Surface>();
                    surface->dimension = dk;
                    surface->normal = normal;
                    surface->neighboring_cell = k + options_->dimensional_cells[dk] * (j + options_->dimensional_cells[dj]* i);
                    surfaces_.push_back(surface);
                    
                    for (int ni = i; ni <= i + 1; ++ni)
                    {
                        for (int nj = j; nj <= j + 1; ++nj)
                        {
                            int n_index = nk + dimensional_nodes_[dk] * (nj + dimensional_nodes_[dj] * ni);
                            shared_ptr<Node> node = nodes_[n_index];
                            node->neighboring_surfaces.push_back(index);
                        }
                    }
                    
                    index += 1;
                }
            }
        }
        break;
    }
    }
    Assert(surfaces_.size() == number_of_surfaces_);
    
    // Get KD tree
    vector<vector<double> > kd_positions(number_of_nodes_, vector<double>(dimension_));
    for (int i = 0; i < number_of_nodes_; ++i)
    {
        kd_positions[i] = nodes_[i]->position;
    }
    node_tree_ = make_shared<KD_Tree>(dimension_,
                                      number_of_nodes_,
                                      kd_positions);

    // Get maximum interval
    max_interval_ = *max_element(intervals_.begin(), intervals_.end());
}

void Integration_Mesh::
initialize_connectivity()
{
    // Get weight function connectivity
    for (int i = 0; i < number_of_points_; ++i)
    {
        // Get weight function data
        shared_ptr<Weight_Function> const weight = weights_[i];
        double const radius = get_inclusive_radius(weight->radius());
        vector<double> const position = weight->position();
        
        // Find nodes that intersect with the weight function
        vector<int> intersecting_nodes;
        vector<double> distances;
        int number_of_intersecting_nodes
            = node_tree_->radius_search(radius,
                                        position,
                                        intersecting_nodes,
                                        distances);
        
        // Add weight indices to cells and surfaces
        for (int j = 0; j < number_of_intersecting_nodes; ++j)
        {
            shared_ptr<Node> node = nodes_[intersecting_nodes[j]];
            
            for (int c_index : node->neighboring_cells)
            {
                cells_[c_index]->weight_indices.push_back(i);
            }
            
            for (int s_index : node->neighboring_surfaces)
            {
                surfaces_[s_index]->weight_indices.push_back(i);
            }
        }
    }
    
    // Get basis function connectivity
    if (!options_->identical_basis_functions)
    {
        for (int i = 0; i < number_of_points_; ++i)
        {
            // Get basis function data
            shared_ptr<Basis_Function> const basis = bases_[i];
            double radius = get_inclusive_radius(basis->radius());
            vector<double> const position = basis->position();

            // Find nodes that intersect with the basis function
            vector<int> intersecting_nodes;
            vector<double> distances;
            int number_of_intersecting_nodes
                = node_tree_->radius_search(radius,
                                            position,
                                            intersecting_nodes,
                                            distances);
        
            // Add basis indices to cells and surfaces
            for (int j = 0; j < number_of_intersecting_nodes; ++j)
            {
                shared_ptr<Node> node = nodes_[intersecting_nodes[j]];
            
                for (int c_index : node->neighboring_cells)
                {
                    cells_[c_index]->basis_indices.push_back(i);
                }

                for (int s_index : node->neighboring_surfaces)
                {
                    surfaces_[s_index]->basis_indices.push_back(i);
                }
            }
        }
    }
    
    // Remove duplicate volume indices
    for (int i = 0; i < number_of_cells_; ++i)
    {
        shared_ptr<Cell> cell = cells_[i];
        
        // Remove duplicate weight indices
        sort(cell->weight_indices.begin(), cell->weight_indices.end());
        cell->weight_indices.erase(unique(cell->weight_indices.begin(), cell->weight_indices.end()), cell->weight_indices.end());
        cell->number_of_weight_functions = cell->weight_indices.size();

        // Remove duplicate basis indices
        if (options_->identical_basis_functions)
        {
            cell->basis_indices = cell->weight_indices;
            cell->number_of_basis_functions = cell->number_of_weight_functions;
        }
        else
        {
            sort(cell->basis_indices.begin(), cell->basis_indices.end());
            cell->basis_indices.erase(unique(cell->basis_indices.begin(), cell->basis_indices.end()), cell->basis_indices.end());
            cell->number_of_basis_functions = cell->basis_indices.size();
        }
    }
    
    // Remove duplicate surface indices
    for (int i = 0; i < number_of_surfaces_; ++i)
    {
        shared_ptr<Surface> surface = surfaces_[i];

        // Remove duplicate weight indices
        sort(surface->weight_indices.begin(), surface->weight_indices.end());
        surface->weight_indices.erase(unique(surface->weight_indices.begin(), surface->weight_indices.end()), surface->weight_indices.end());
        surface->number_of_weight_functions = surface->weight_indices.size();

        // Remove duplicate basis indices
        if (options_->identical_basis_functions)
        {
            surface->basis_indices = surface->weight_indices;
            surface->number_of_basis_functions = surface->number_of_weight_functions;
        }
        else
        {
            sort(surface->basis_indices.begin(), surface->basis_indices.end());
            surface->basis_indices.erase(unique(surface->basis_indices.begin(), surface->basis_indices.end()), surface->basis_indices.end());
            surface->number_of_basis_functions = surface->basis_indices.size();
        }
    }

    // If applicable, change the number of integration ordinates
    if (options_->adaptive_quadrature)
    {
        // Get background cell integration ordinates
        for (int i = 0; i < number_of_cells_; ++i)
        {
            shared_ptr<Cell> cell = cells_[i];
            
            // Find minimum radius
            double min_radius = numeric_limits<double>::max();
            for (int j = 0; j < cell->number_of_weight_functions; ++j)
            {
                int const k = cell->weight_indices[j];
                shared_ptr<Weight_Function> const weight = weights_[k];
                double const radius = weight->radius();
                
                if (radius < min_radius)
                {
                    min_radius = radius;
                }
            }
            if (!options_->identical_basis_functions)
            {
                for (int j = 0; j < cell->number_of_basis_functions; ++j)
                {
                    int const k = cell->basis_indices[j];
                    shared_ptr<Basis_Function> const basis = bases_[k];
                    double const radius = basis->radius();
                    
                    if (radius < min_radius)
                    {
                        min_radius = radius;
                    }
                }
            }
            
            // Find min cell length
            double min_length = numeric_limits<double>::max();
            for (int d = 0; d < dimension_; ++d)
            {
                double length = cell->limits[d][1] - cell->limits[d][0];
                if (length < min_length)
                {
                    min_length = length;
                }
            }
            
            // Get expected number of integration ordinates
            int expected_number = ceil(min_length / min_radius * options_->minimum_radius_ordinates);
            
            // Compare to actual number of integration ordinates
            int const global_integration_ordinates = options_->integration_ordinates;
            if (expected_number > global_integration_ordinates)
            {
                cell->number_of_integration_ordinates = expected_number;
            }
            else
            {
                cell->number_of_integration_ordinates = global_integration_ordinates;
            }
        }
        
        // Get background surface integration ordinates
        for (int i = 0; i < number_of_surfaces_; ++i)
        {
            shared_ptr<Surface> surface = surfaces_[i];
            
            // Find minimum radius
            double min_radius = numeric_limits<double>::max();
            for (int j = 0; j < surface->number_of_weight_functions; ++j)
            {
                int const k = surface->weight_indices[j];
                shared_ptr<Weight_Function> const weight = weights_[k];
                double const radius = weight->radius();
                
                if (radius < min_radius)
                {
                    min_radius = radius;
                }
            }
            if (!options_->identical_basis_functions)
            {
                for (int j = 0; j < surface->number_of_basis_functions; ++j)
                {
                    int k = surface->basis_indices[j];
                    shared_ptr<Basis_Function> const basis = bases_[k];
                    double const radius = basis->radius();
                
                    if (radius < min_radius)
                    {
                        min_radius = radius;
                    }
                }
            }
            
            // Find min surface length
            shared_ptr<Cell> cell = cells_[surface->neighboring_cell];
            double min_length = numeric_limits<double>::max();
            for (int d = 0; d < dimension_; ++d)
            {
                if (d != surface->dimension)
                {
                    double length = cell->limits[d][1] - cell->limits[d][0];
                    if (length < min_length)
                    {
                        min_length = length;
                    }
                }
            }
            
            // Get expected number of integration ordinates
            int expected_number = ceil(min_length / min_radius * options_->minimum_radius_ordinates);
            
            // Compare to actual number of integration ordinates
            int const global_integration_ordinates = options_->integration_ordinates;
            if (expected_number > global_integration_ordinates)
            {
                surface->number_of_integration_ordinates = expected_number;
            }
            else
            {
                surface->number_of_integration_ordinates = global_integration_ordinates;
            }
        }
    }
    else
    {
        // Set number of integration ordinates for each cell and surface to global value
        for (int i = 0; i < number_of_cells_; ++i)
        {
            shared_ptr<Cell> cell = cells_[i];

            cell->number_of_integration_ordinates = options_->integration_ordinates;
        }
        for (int i = 0; i < number_of_surfaces_; ++i)
        {
            shared_ptr<Surface> surface = surfaces_[i];
        
            surface->number_of_integration_ordinates = options_->integration_ordinates;
        }
        
    }

}

double Integration_Mesh::
get_inclusive_radius(double radius) const
{
    // Move radius outward to account for surfaces
    switch (dimension_)
    {
    case 1:
        return radius;
    case 2:
        return sqrt(radius * radius + 0.25 * max_interval_ * max_interval_);
    case 3:
        return sqrt(radius * radius + 0.5 * max_interval_ * max_interval_);
    default:
        AssertMsg(false, "dimension cannot exceed 3");
        return -1.;
    }
}

void Integration_Mesh::
get_volume_quadrature(int i,
                      int &number_of_ordinates,
                      vector<vector<double> > &ordinates,
                      vector<double> &weights) const
{
    // Get limits of integration
    shared_ptr<Cell> const cell = cells_[i];
    vector<vector<double> > const &limits = cell->limits;
    int const number_of_integration_ordinates = cell->number_of_integration_ordinates;
    int const dx = 0;
    int const dy = 1;
    int const dz = 2;
    int const min = 0;
    int const max = 1;

    // Initialize temporary integration data
    Quadrature_Rule::Quadrature_Type quad_type = Quadrature_Rule::Quadrature_Type::GAUSS_LEGENDRE;
    vector<double> ordinates_x;
    vector<double> ordinates_y;
    vector<double> ordinates_z;

    // Get quadrature
    switch (dimension_)
    {
    case 1:
        Quadrature_Rule::cartesian_1d(quad_type,
                                      number_of_integration_ordinates,
                                      limits[dx][min],
                                      limits[dx][max],
                                      ordinates_x,
                                      weights);
        Quadrature_Rule::convert_to_position_1d(ordinates_x,
                                                ordinates);
        break;
    case 2:
        Quadrature_Rule::cartesian_2d(quad_type,
                                      quad_type,
                                      number_of_integration_ordinates,
                                      number_of_integration_ordinates,
                                      limits[dx][min],
                                      limits[dx][max],
                                      limits[dy][min],
                                      limits[dy][max],
                                      ordinates_x,
                                      ordinates_y,
                                      weights);
        Quadrature_Rule::convert_to_position_2d(ordinates_x,
                                                ordinates_y,
                                                ordinates);
        break;
    case 3:
        Quadrature_Rule::cartesian_3d(quad_type,
                                      quad_type,
                                      quad_type,
                                      number_of_integration_ordinates,
                                      number_of_integration_ordinates,
                                      number_of_integration_ordinates,
                                      limits[dx][min],
                                      limits[dx][max],
                                      limits[dy][min],
                                      limits[dy][max],
                                      limits[dz][min],
                                      limits[dz][max],
                                      ordinates_x,
                                      ordinates_y,
                                      ordinates_z,
                                      weights);
        Quadrature_Rule::convert_to_position_3d(ordinates_x,
                                                ordinates_y,
                                                ordinates_z,
                                                ordinates);
        break;
    }

    // Set number of ordinates
    number_of_ordinates = weights.size();
}

void Integration_Mesh::
get_surface_quadrature(int i,
                       int &number_of_ordinates,
                       vector<vector<double> > &ordinates,
                       vector<double> &weights) const
{
    // Get limits of integration
    shared_ptr<Surface> const surface = surfaces_[i];
    shared_ptr<Cell> const cell = cells_[surface->neighboring_cell];
    vector<vector<double> > const &limits = cell->limits;
    int const number_of_integration_ordinates = surface->number_of_integration_ordinates;
    int const dx = 0;
    int const dy = 1;
    int const dz = 2;
    int const min = 0;
    int const max = 1;

    // Get temporary integration data
    Quadrature_Rule::Quadrature_Type quad_type = Quadrature_Rule::Quadrature_Type::GAUSS_LEGENDRE;
    vector<double> ordinates_x;
    vector<double> ordinates_y;
    vector<double> ordinates_z;

    // Get quadrature
    switch (dimension_)
    {
    case 1:
        number_of_ordinates = 1;
        if (surface->normal < 0)
        {
            // At negative x boundary
            ordinates.assign(1, vector<double>(1, limits[dx][min]));
            weights.assign(1, 1.);
        }
        else
        {
            // At positive x boundary
            ordinates.assign(1, vector<double>(1, limits[dx][max]));
            weights.assign(1, 1.);
        }
        break;
    case 2:
        switch (surface->dimension)
        {
        case dx:
            Quadrature_Rule::cartesian_1d(quad_type,
                                          number_of_integration_ordinates,
                                          limits[dy][min],
                                          limits[dy][max],
                                          ordinates_y,
                                          weights);
            number_of_ordinates = weights.size();
            if (surface->normal < 0)
            {
                // At negative x boundary
                ordinates_x.assign(number_of_ordinates,
                                   limits[dx][min]);
            }
            else
            {
                // At positive x boundary
                ordinates_x.assign(number_of_ordinates,
                                   limits[dx][max]);
            }
            break;
        case dy:
            Quadrature_Rule::cartesian_1d(quad_type,
                                          number_of_integration_ordinates,
                                          limits[dx][min],
                                          limits[dx][max],
                                          ordinates_x,
                                          weights);
            number_of_ordinates = weights.size();
            if (surface->normal < 0)
            {
                // At negative x boundary
                ordinates_y.assign(number_of_ordinates,
                                   limits[dy][min]);
            }
            else
            {
                // At positive x boundary
                ordinates_y.assign(number_of_ordinates,
                                   limits[dy][max]);
            }
            break;
        }
        Quadrature_Rule::convert_to_position_2d(ordinates_x,
                                                ordinates_y,
                                                ordinates);
        break;
    case 3:
        switch (surface->dimension)
        {
        case dx:
            Quadrature_Rule::cartesian_2d(quad_type,
                                          quad_type,
                                          number_of_integration_ordinates,
                                          number_of_integration_ordinates,
                                          limits[dy][min],
                                          limits[dy][max],
                                          limits[dz][min],
                                          limits[dz][max],
                                          ordinates_y,
                                          ordinates_z,
                                          weights);
            number_of_ordinates = weights.size();
            if (surface->normal < 0)
            {
                // At negative x boundary
                ordinates_x.assign(number_of_ordinates,
                                   limits[dx][min]);
            }
            else
            {
                // At positive x boundary
                ordinates_x.assign(number_of_ordinates,
                                   limits[dx][max]);
            }
            break;
        case dy:
            Quadrature_Rule::cartesian_2d(quad_type,
                                          quad_type,
                                          number_of_integration_ordinates,
                                          number_of_integration_ordinates,
                                          limits[dx][min],
                                          limits[dx][max],
                                          limits[dz][min],
                                          limits[dz][max],
                                          ordinates_x,
                                          ordinates_z,
                                          weights);
            number_of_ordinates = weights.size();
            if (surface->normal < 0)
            {
                // At negative y boundary
                ordinates_y.assign(number_of_ordinates,
                                   limits[dy][min]);
            }
            else
            {
                // At positive y boundary
                ordinates_y.assign(number_of_ordinates,
                                   limits[dy][max]);
            }
            break;
        case dz:
            Quadrature_Rule::cartesian_2d(quad_type,
                                          quad_type,
                                          number_of_integration_ordinates,
                                          number_of_integration_ordinates,
                                          limits[dx][min],
                                          limits[dx][max],
                                          limits[dy][min],
                                          limits[dy][max],
                                          ordinates_x,
                                          ordinates_y,
                                          weights);
            number_of_ordinates = weights.size();
            if (surface->normal < 0)
            {
                // At negative z boundary
                ordinates_z.assign(number_of_ordinates,
                                   limits[dz][min]);
            }
            else
            {
                // At positive z boundary
                ordinates_z.assign(number_of_ordinates,
                                   limits[dz][max]);
            }
            break;
        }
        Quadrature_Rule::convert_to_position_3d(ordinates_x,
                                                ordinates_y,
                                                ordinates_z,
                                                ordinates);
        break;
    }
}

void Integration_Mesh::
get_basis_values(shared_ptr<Cell> const cell,
                 vector<double> const &position,
                 vector<vector<double> > const &basis_centers,
                 vector<double> &b_val) const
{
    // Initialize values 
    b_val.resize(cell->number_of_basis_functions);
    
    // Get values for basis functions at quadrature point
    for (int j = 0; j < cell->number_of_basis_functions; ++j)
    {
        shared_ptr<Meshless_Function> const func = bases_[cell->basis_indices[j]]->function()->base_function();
                
        b_val[j] = func->value(position);
    }

    // Normalize basis functions
    if (apply_basis_normalization_)
    {
        basis_normalization_->get_values(position,
                                         basis_centers,
                                         b_val,
                                         b_val);
    }
}

void Integration_Mesh::
get_volume_values(shared_ptr<Cell> const cell,
                  vector<double> const &position,
                  vector<vector<double> > const &basis_centers,
                  vector<vector<double> > const &weight_centers,
                  vector<double> &b_val,
                  vector<vector<double> > &b_grad,
                  vector<double> &w_val,
                  vector<vector<double> > &w_grad) const
{
    // Initialize values 
    b_val.resize(cell->number_of_basis_functions);
    b_grad.resize(cell->number_of_basis_functions);
    w_val.resize(cell->number_of_weight_functions);
    w_grad.resize(cell->number_of_weight_functions);

    // Get values for weight functions at quadrature point
    for (int j = 0; j < cell->number_of_weight_functions; ++j)
    {
        shared_ptr<Meshless_Function> const func = weights_[cell->weight_indices[j]]->function()->base_function();
        
        w_val[j] = func->value(position);
        w_grad[j] = func->gradient_value(position);
    }
    
    // Normalize weight functions
    if (apply_weight_normalization_)
    {
        weight_normalization_->get_gradient_values(position,
                                                   weight_centers,
                                                   w_val,
                                                   w_grad,
                                                   w_val,
                                                   w_grad);
    }
    
    if (options_->identical_basis_functions)
    {
        b_val = w_val;
        b_grad = w_grad;
    }
    else
    {
        // Get values for basis functions at quadrature point
        for (int j = 0; j < cell->number_of_basis_functions; ++j)
        {
            shared_ptr<Meshless_Function> const func = bases_[cell->basis_indices[j]]->function()->base_function();
                
            b_val[j] = func->value(position);
            b_grad[j] = func->gradient_value(position);
        }

        // Normalize basis functions
        if (apply_basis_normalization_)
        {
            basis_normalization_->get_gradient_values(position,
                                                      basis_centers,
                                                      b_val,
                                                      b_grad,
                                                      b_val,
                                                      b_grad);
        }
    }
}

void Integration_Mesh::
get_surface_values(shared_ptr<Surface> const surface,
                   vector<double> const &position,
                   vector<vector<double> > const &basis_centers,
                   vector<vector<double> > const &weight_centers,
                   vector<double> &b_val,
                   vector<double> &w_val) const
{
    // Initialize values
    b_val.resize(surface->number_of_basis_functions);
    w_val.resize(surface->number_of_weight_functions);
    
    // Get values for weight functions at quadrature point
    for (int j = 0; j < surface->number_of_weight_functions; ++j)
    {
        shared_ptr<Meshless_Function> const func = weights_[surface->weight_indices[j]]->function()->base_function();
                
        w_val[j] = func->value(position);
    }
    
    // Normalize weight functions
    if (apply_weight_normalization_)
    {
        weight_normalization_->get_values(position,
                                          weight_centers,
                                          w_val,
                                          w_val);
    }

    if (options_->identical_basis_functions)
    {
        b_val = w_val;
    }
    else
    {
        // Get values for basis functions at quadrature point
        for (int j = 0; j < surface->number_of_basis_functions; ++j)
        {
            shared_ptr<Meshless_Function> const func = bases_[surface->basis_indices[j]]->function()->base_function();
                
            b_val[j] = func->value(position);
        }
        if (apply_basis_normalization_)
        {
            basis_normalization_->get_values(position,
                                             basis_centers,
                                             b_val,
                                             b_val);
        }
    }
}

void Integration_Mesh::
get_cell_basis_indices(shared_ptr<Cell> const cell,
                       vector<vector<int> > &indices) const
{
    indices.assign(cell->number_of_weight_functions,
                   vector<int>(cell->number_of_basis_functions, -1));
    for (int i = 0; i < cell->number_of_weight_functions; ++i)
    {
        shared_ptr<Weight_Function> const weight = weights_[cell->weight_indices[i]];
        for (int j = 0; j < cell->number_of_basis_functions; ++j)
        {
            indices[i][j] = weight->local_basis_index(cell->basis_indices[j]);
        }
    }
}

void Integration_Mesh::
get_surface_basis_indices(shared_ptr<Surface> const surface,
                          vector<vector<int> > &indices) const
{
    indices.assign(surface->number_of_weight_functions,
                   vector<int>(surface->number_of_basis_functions, -1));
    for (int i = 0; i < surface->number_of_weight_functions; ++i)
    {
        shared_ptr<Weight_Function> const weight = weights_[surface->weight_indices[i]];
        for (int j = 0; j < surface->number_of_basis_functions; ++j)
        {
            indices[i][j] = weight->local_basis_index(surface->basis_indices[j]);
        }
    }
}

void Integration_Mesh::
get_weight_surface_indices(shared_ptr<Surface> const surface,
                           vector<int> &indices) const
{
    indices.assign(surface->number_of_weight_functions, Weight_Function::Errors::DOES_NOT_EXIST);
    for (int i = 0; i < surface->number_of_weight_functions; ++i)
    {
        shared_ptr<Weight_Function> const weight = weights_[surface->weight_indices[i]];
        indices[i] = weight->local_surface_index(surface->dimension,
                                                 surface->normal);
    }
}

void Integration_Mesh::
get_basis_centers(shared_ptr<Cell> const cell,
                  vector<vector<double> > &basis_positions) const
{
    int const number_of_basis_functions = cell->number_of_basis_functions;
    basis_positions.resize(number_of_basis_functions);
    for (int i = 0; i < number_of_basis_functions; ++i)
    {
        basis_positions[i] = bases_[cell->basis_indices[i]]->position();
    }
}

void Integration_Mesh::
get_basis_weight_centers(shared_ptr<Cell> const cell,
                         vector<vector<double> > &basis_positions,
                         vector<vector<double> > &weight_positions) const
{
    int const number_of_weight_functions = cell->number_of_weight_functions;
    weight_positions.resize(number_of_weight_functions);
    for (int i = 0; i < number_of_weight_functions; ++i)
    {
        weight_positions[i] = weights_[cell->weight_indices[i]]->position();
    }

    if (options_->identical_basis_functions)
    {
        basis_positions = weight_positions;
    }
    else
    {
        int const number_of_basis_functions = cell->number_of_basis_functions;
        basis_positions.resize(number_of_basis_functions);
        for (int i = 0; i < number_of_basis_functions; ++i)
        {
            basis_positions[i] = bases_[cell->basis_indices[i]]->position();
        }
    }
}

void Integration_Mesh::
get_basis_weight_centers(shared_ptr<Surface> const surface,
                         vector<vector<double> > &basis_positions,
                         vector<vector<double> > &weight_positions) const
{
    int const number_of_weight_functions = surface->number_of_weight_functions;
    weight_positions.resize(number_of_weight_functions);
    for (int i = 0; i < number_of_weight_functions; ++i)
    {
        weight_positions[i] = weights_[surface->weight_indices[i]]->position();
    }

    if (options_->identical_basis_functions)
    {
        basis_positions = weight_positions;
    }
    else
    {
        int const number_of_basis_functions = surface->number_of_basis_functions;
        basis_positions.resize(number_of_basis_functions);
        for (int i = 0; i < number_of_basis_functions; ++i)
        {
            basis_positions[i] = bases_[surface->basis_indices[i]]->position();
        }
    }
}

void Integration_Mesh::Options::
initialize_from_weak_options(std::shared_ptr<Weak_Spatial_Discretization_Options> weak_options)
{
    identical_basis_functions = (weak_options->identical_basis_functions == Weak_Spatial_Discretization_Options::Identical_Basis_Functions::TRUE);
    adaptive_quadrature = weak_options->adaptive_quadrature;
    minimum_radius_ordinates = weak_options->minimum_radius_ordinates;
    integration_ordinates = weak_options->integration_ordinates;
    limits = weak_options->limits;
    dimensional_cells = weak_options->dimensional_cells;
}
