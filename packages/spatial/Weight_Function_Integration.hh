#ifndef Weight_Function_Integration_hh
#define Weight_Function_Integration_hh

#include "Weight_Function.hh"

class KD_Tree;

class Weight_Function_Integration
{
public:
    
    class Mesh
    {
    public:

        class Cell
        {
            std::vector<int> basis_indices;
            std::vector<int> weight_indices;
            std::vector<std::vector<double> > limits;
            std::vector<int> neighboring_nodes;
        };
        
        class Node
        {
            std::vector<double> position;
            std::vector<int> neighboring_cells;
        };
        
        Mesh(Weight_Function_Integration const &wfi,
             int dimension,
             std::vector<std::vector<double> > limits,
             std::vector<int> num_intervals);
        
    private:

        std::shared_ptr<KD_Tree> get_kd_tree() const;

        Weight_Function_Integration const &wfi_;
        int number_of_local_nodes_;
        int number_of_global_nodes_;
        int number_of_global_cells_;
        int dimension_;
        std::vector<std::vector<double> > limits_;
        std::vector<int> dimensional_cells_;
        std::vector<int> dimensional_nodes_;
        std::vector<double> intervals_;
        std::shared_ptr<KD_Tree> kd_tree_;
        vector<Cell> cells_;
        vector<Node> nodes_;
    };
    
    Weight_Function_Integration(std::vector<std::shared_ptr<Basis_Function> > const &bases,
                                std::vector<std::shared_ptr<Weight_Function> > const &weights,
                                std::shared_ptr<Solid_Geometry> solid_geometry,
                                std::vector<std::vector<double> > limits,
                                std::vector<int> num_intervals);
    
private:
    
    Weight_Function::Options options_;
    std::vector<std::shared_ptr<Basis_Function> > bases_
    std::vector<std::shared_ptr<Weight_Function> > weights_;
    std::shared_ptr<Solid_Geometry> solid_;
    std::shared_ptr<Mesh> mesh_;
    std::vector<std::vector<double> > limits_;
};

#endif
