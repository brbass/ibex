#ifndef Weight_Function_Integration_hh
#define Weight_Function_Integration_hh

#include <memory>
#include <vector>

#include "Weight_Function.hh"

class Basis_Function;
class KD_Tree;
class Solid_Geometry;

class Weight_Function_Integration
{
public:
    
    class Mesh
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
        };
        
        Mesh(Weight_Function_Integration const &wfi,
             int dimension,
             std::vector<std::vector<double> > limits,
             std::vector<int> dimensional_cells);
        
    private:

        friend class Weight_Function_Integration;
        
        void initialize_mesh();
        void initialize_connectivity();
        double get_inclusive_radius(double radius) const;

        Weight_Function_Integration const &wfi_;
        int number_of_background_nodes_;
        int number_of_background_cells_;
        int number_of_background_surfaces_;
        int dimension_;
        double max_interval_;
        std::vector<std::vector<double> > limits_;
        std::vector<int> dimensional_cells_;
        std::vector<int> dimensional_nodes_;
        std::vector<double> intervals_;
        std::shared_ptr<KD_Tree> node_tree_;
        std::vector<Cell> cells_;
        std::vector<Node> nodes_;
        std::vector<Surface> surfaces_;
    };

    struct Material_Data
    {
        std::vector<double> sigma_t;
        std::vector<double> sigma_s;
        std::vector<double> nu;
        std::vector<double> sigma_f;
        std::vector<double> chi;
        std::vector<double> internal_source;
        std::vector<double> norm;
    };

    // Constructor
    Weight_Function_Integration(int number_of_points,
                                std::vector<std::shared_ptr<Basis_Function> > const &bases,
                                std::vector<std::shared_ptr<Weight_Function> > const &weights,
                                std::shared_ptr<Solid_Geometry> solid_geometry,
                                std::vector<std::vector<double> > limits,
                                std::vector<int> num_intervals);

    // Perform integration and put result into weight functions
    void perform_integration();
    
private:

    // Put volume, surface and material integrals into weight functions
    void put_integrals_into_weight(std::vector<Weight_Function::Integrals> const &integrals,
                                   std::vector<Material_Data> const &materials);
    
    // Perform all volume integrals
    void perform_volume_integration(std::vector<Weight_Function::Integrals> &integrals,
                                    std::vector<Material_Data> &materials) const;


    // Add weight function cell values to global integrals
    void add_volume_weight(Mesh::Cell const &cell,
                           double quad_weight,
                           std::vector<double> const &w_val,
                           std::vector<std::vector<double> > const &w_grad,
                           std::vector<Weight_Function::Integrals> &integrals) const;

    // Add basis/weight function cell values to global integrals
    void add_volume_basis_weight(Mesh::Cell const &cell,
                                 double quad_weight,
                                 std::vector<double> const &b_val,
                                 std::vector<std::vector<double> > const &b_grad,
                                 std::vector<double> const &w_val,
                                 std::vector<std::vector<double> > const &w_grad,
                                 std::vector<std::vector<int> > const &weight_basis_indices,
                                 std::vector<Weight_Function::Integrals> &integrals) const;

    // Add material cell values to global integrals
    void add_volume_material(Mesh::Cell const &cell,
                             double quad_weight,
                             std::vector<double> const &w_val,
                             std::vector<std::vector<double> > const &w_grad,
                             shared_ptr<Material> point_material,
                             std::vector<Material_Data> &materials) const;
    
    // Get values for a single quadrature point in a cell
    void get_volume_values(Mesh::Cell const &cell,
                           std::vector<double> const &position,
                           std::vector<double> &b_val,
                           std::vector<std::vector<double> > &b_grad,
                           std::vector<double> &w_val,
                           std::vector<std::vector<double> > &w_grad,
                           std::shared_ptr<Material> &point_material) const;

    // Get a Cartesian volume quadrature for the cell
    void get_volume_quadrature(int i, // cell index
                               int &number_of_ordinates,
                               std::vector<std::vector<double> > &ordinates,
                               std::vector<double> &weights) const;

    // Perform all surface integrals
    void perform_surface_integration(std::vector<Weight_Function::Integrals> &integrals) const;

    // Add basis/weight function surface values to global integrals
    void add_surface_basis_weight(Mesh::Surface const &surface,
                                  double quad_weight,
                                  std::vector<double> const &b_val,
                                  std::vector<double> const &w_val,
                                  std::vector<int> const &weight_surface_indices,
                                  std::vector<std::vector<int> > const &weight_basis_indices,
                                  std::vector<Weight_Function::Integrals> &integrals) const;
    
    // Get values for a single quadrature point on a surface
    void get_surface_values(Mesh::Surface const &surface,
                            std::vector<double> const &position,
                            std::vector<double> &b_val,
                            std::vector<double> &w_val) const;

    // Get a Cartesian surface quadrature for the surface
    void get_surface_quadrature(int i, // surface index
                                int &number_of_ordinates,
                                std::vector<std::vector<double> > &ordinates,
                                vector<double> &weights) const;
    
    // Get the local basis indices for the weight functions in a single cell
    void get_cell_basis_indices(Mesh::Cell const &cell,
                                 std::vector<std::vector<int> > &indices) const;
    
    // Get the local basis indices for the weight functions in a single surface
    void get_surface_basis_indices(Mesh::Surface const &surface,
                                   std::vector<std::vector<int> > &indices) const;

    // Get the local surface indices for the weight functions intersecting with a surface
    void get_weight_surface_indices(Mesh::Surface const &surface,
                                    std::vector<int> &indices) const;

    // Get a material from the material data
    void get_material(Material_Data const &material_data,
                      shared_ptr<Material> &material) const;
    
    // Initialize material data to zero
    void initialize_materials(std::vector<Material_Data> &materials) const;

    // Initialize integral data to zero
    void initialize_integrals(std::vector<Weight_Function::Integrals> &integrals) const;

    // Data
    Weight_Function::Options options_;
    int number_of_points_;
    std::vector<std::shared_ptr<Basis_Function> > bases_;
    std::vector<std::shared_ptr<Weight_Function> > weights_;
    std::shared_ptr<Solid_Geometry> solid_;
    std::shared_ptr<Mesh> mesh_;
    std::vector<std::vector<double> > limits_;
};

#endif
