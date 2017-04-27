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

        struct Cell
        {
            // Spatial information
            std::vector<std::vector<double> > limits;

            // Connectivity
            std::vector<int> neighboring_nodes;
            
            // Intersections
            std::vector<int> basis_indices;
            std::vector<int> weight_indices;
        };
        
        struct Node
        {
            // Spatial information
            std::vector<double> position;

            // Connectivity
            std::vector<int> neighboring_cells;
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
        int dimension_;
        double max_interval_;
        std::vector<std::vector<double> > limits_;
        std::vector<int> dimensional_cells_;
        std::vector<int> dimensional_nodes_;
        std::vector<double> intervals_;
        std::shared_ptr<KD_Tree> kd_tree_;
        std::vector<Cell> cells_;
        std::vector<Node> nodes_;
    };

    struct Material_Data
    {
        std::vector<double> sigma_t;
        std::vector<double> sigma_s;
        std::vector<double> nu;
        std::vector<double> sigma_f;
        std::vector<double> chi;
        std::vector<double> internal_source;
    };
    
    Weight_Function_Integration(int number_of_points,
                                std::vector<std::shared_ptr<Basis_Function> > const &bases,
                                std::vector<std::shared_ptr<Weight_Function> > const &weights,
                                std::shared_ptr<Solid_Geometry> solid_geometry,
                                std::vector<std::vector<double> > limits,
                                std::vector<int> num_intervals);

    void perform_integration();
    
private:
    
    Weight_Function::Options options_;
    int number_of_points_;
    std::vector<std::shared_ptr<Basis_Function> > bases_;
    std::vector<std::shared_ptr<Weight_Function> > weights_;
    std::shared_ptr<Solid_Geometry> solid_;
    std::shared_ptr<Mesh> mesh_;
    std::vector<std::vector<double> > limits_;
};

#endif
