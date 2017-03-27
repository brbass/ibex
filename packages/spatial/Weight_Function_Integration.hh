#ifndef Weight_Function_Integration_hh
#define Weight_Function_Integration_hh

#include <vector>

class Constructive_Solid_Geometry

class Weight_Function_Integration
{
public:

    struct Cell
    {
        // Data
        int index;
        std::vector<double> limits;
        
        // Connectivity
        std::vector<int> nodes;

        // Information
        std::vector<int> bases;
        std::vector<int> weights;
    };

    struct Node
    {
        // Data
        int index;
        vector<double> position;

        // Connectivity
        std::vector<int> cells;
        std::vector<int> surfaces;
    };

    struct Surface
    {
        // Data
        int index;
        
        // Connectivity
        std::vector<int> nodes;
    };
    
    
    Weight_Function_Integration(std::shared_ptr<Constructive_Solid_Geometry> solid_geometry);

    void get_weight_functions(int number_of_points,
                              std::vector<Weight_Function::Options> const &options,
                              std::vector<std::vector<int> > const &neighbors,
                              std::vector<std::shared_ptr<Basis_Function> > const &bases,
                              std::vector<std::shared_ptr<Weight_Function> > &weights);
    
private:

    void create_mesh();
    void create_mesh_1d();
    void create_mesh_2d();
    
    // Data
    std::vector<int> number_of_intervals_;
    std::vector<std::vector<double> > limits_;
    std::shared_ptr<Constructive_Solid_Geometry> solid_geometry_;
    
    // Geometry
    int dimension_;
    int number_of_cells_;
    int number_of_nodes_;
    int number_of_surfaces_;
    std::vector<Cell> cells_;
    std::vector<Node> nodes_;
    std::vector<Surface> surfaces_;
};

#endif
