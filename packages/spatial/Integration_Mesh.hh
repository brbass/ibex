#ifndef Integration_Mesh_hh
#define Integration_Mesh_hh

#include <memory>
#include <vector>

class Basis_Function;
class KD_Tree;
class Meshless_Normalization;
class Weak_Spatial_Discretization_Options;
class Weight_Function;

struct Integration_Mesh_Options
{
    bool identical_basis_functions = false;
    bool adaptive_quadrature = false;
    int minimum_radius_ordinates = 12;
    int integration_ordinates = 8;
    int maximum_integration_ordinates = 128;
    double boundary_tolerance = 1e-10;
    std::vector<std::vector<double> > limits;
    std::vector<int> dimensional_cells;
        
    void initialize_from_weak_options(std::shared_ptr<Weak_Spatial_Discretization_Options> weak_options);
};
    
class Integration_Mesh
{
public:

    // Cartesian cell
    struct Cell
    {
        // Spatial information
        std::vector<std::vector<double> > limits;

        // Connectivity
        std::vector<int> neighboring_nodes;
            
        // Intersections
        int number_of_basis_functions;
        int number_of_weight_functions;
        std::vector<int> basis_indices;
        std::vector<int> weight_indices;

        // Integration
        int number_of_integration_ordinates;
    };

    // Single point in the mesh
    struct Node
    {
        // Spatial information
        std::vector<double> position;

        // Connectivity
        std::vector<int> neighboring_cells;
        std::vector<int> neighboring_surfaces;
    };

    // Boundary surface in the cell
    struct Surface
    {
        // Spatial information
        int dimension;
        double normal;
            
        // Connectivity
        int neighboring_cell;

        // Intersections
        int number_of_basis_functions;
        int number_of_weight_functions;
        std::vector<int> basis_indices;
        std::vector<int> weight_indices;

        // Integration
        int number_of_integration_ordinates;
    };
        
    Integration_Mesh(int dimension,
                     int number_of_points,
                     std::shared_ptr<Integration_Mesh_Options> options,
                     std::vector<std::shared_ptr<Basis_Function> > const &bases,
                     std::vector<std::shared_ptr<Weight_Function> > const &weights);

    // Access data
    int dimension() const
    {
        return dimension_;
    }
    int number_of_points() const
    {
        return number_of_points_;
    }
    int number_of_cells() const
    {
        return number_of_cells_;
    }
    int number_of_surfaces() const
    {
        return number_of_surfaces_;
    }
    int number_of_nodes() const
    {
        return number_of_nodes_;
    }
    std::shared_ptr<Cell> const cell(int index) const
    {
        return cells_[index];
    }
    std::shared_ptr<Surface> const surface(int index) const
    {
        return surfaces_[index];
    }
    std::shared_ptr<Node> const node(int index) const
    {
        return nodes_[index];
    }

    // Quadrature methods
    void get_volume_quadrature(int i, // cell index
                               int &number_of_ordinates,
                               std::vector<std::vector<double> > &ordinates,
                               std::vector<double> &weights) const;
    void get_surface_quadrature(int i, // surface index
                                int &number_of_ordinates,
                                std::vector<std::vector<double> > &ordinates,
                                std::vector<double> &weights) const;

    // Get values at a single specified point
    void get_basis_values(std::shared_ptr<Cell> const cell,
                          std::vector<double> const &position,
                          std::vector<std::vector<double> > const &basis_centers,
                          std::vector<double> &b_val) const;
    void get_volume_values(std::shared_ptr<Cell> const cell,
                           std::vector<double> const &position,
                           std::vector<std::vector<double> > const &basis_centers,
                           std::vector<std::vector<double> > const &weight_centers,
                           std::vector<double> &b_val,
                           std::vector<std::vector<double> > &b_grad,
                           std::vector<double> &w_val,
                           std::vector<std::vector<double> > &w_grad) const;
    void get_surface_values(std::shared_ptr<Surface> const surface,
                            std::vector<double> const &position,
                            std::vector<std::vector<double> > const &basis_centers,
                            std::vector<std::vector<double> > const &weight_centers,
                            std::vector<double> &b_val,
                            std::vector<double> &w_val) const;
    
    // Indexing methods
    void get_cell_basis_indices(std::shared_ptr<Cell> const cell,
                                 std::vector<std::vector<int> > &basis_indices) const;
    void get_surface_basis_indices(std::shared_ptr<Surface> const surface,
                                   std::vector<std::vector<int> > &basis_indices) const;
    void get_weight_surface_indices(std::shared_ptr<Surface> const surface,
                                    std::vector<int> &surface_indices) const;

    // Get the positions of the centers basis/weight functions
    void get_basis_centers(std::shared_ptr<Cell> const cell,
                           std::vector<std::vector<double> > &basis_positions) const;
    void get_basis_weight_centers(std::shared_ptr<Cell> const cell,
                                  std::vector<std::vector<double> > &basis_positions,
                                  std::vector<std::vector<double> > &weight_positions) const;
    void get_basis_weight_centers(std::shared_ptr<Surface> const surface,
                                  std::vector<std::vector<double> > &basis_positions,
                                  std::vector<std::vector<double> > &weight_positions) const;

    // Get the cell index for a certain position
    int get_cell_at_position(std::vector<double> const &position) const;
    
private:

    // Initialization methods
    void initialize_mesh();
    void initialize_connectivity();
    double get_inclusive_radius(double radius) const;

    // Input data
    int dimension_;
    int number_of_points_;
    std::shared_ptr<Integration_Mesh_Options> options_;
    std::vector<std::shared_ptr<Basis_Function> > bases_;
    std::vector<std::shared_ptr<Weight_Function> > weights_;
    
    // Derived data
    int number_of_nodes_;
    int number_of_cells_;
    int number_of_surfaces_;
    bool apply_basis_normalization_;
    bool apply_weight_normalization_;
    double max_interval_;
    std::vector<int> dimensional_nodes_;
    std::vector<double> intervals_;
    std::shared_ptr<KD_Tree> node_tree_;
    std::shared_ptr<Meshless_Normalization> basis_normalization_;
    std::shared_ptr<Meshless_Normalization> weight_normalization_;
    std::vector<std::shared_ptr<Cell> > cells_;
    std::vector<std::shared_ptr<Node> > nodes_;
    std::vector<std::shared_ptr<Surface> > surfaces_;
};

#endif
