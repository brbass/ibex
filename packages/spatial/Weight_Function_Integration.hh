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
            std::vector<int> basis_indices_;
            std::vector<int> weight_indices_;
            
        };
        
        class Node
        {
            
        };
        
        Mesh(int dimension,
             std::vector<std::vector<double> > limits,
             std::vector<int> num_intervals);
        
    private:

        std::shared_ptr<KD_Tree> get_kd_tree() const;

        int number_of_points_;
        int number_of_cells_;
        int dimension_;
        std::vector<std::vector<double> > limits_;
        std::vector<int> num_intervals_;
        std::vector<double> intervals_;
        std::shared_ptr<KD_Tree> kd_tree_;
        std::vector<std::vector<double> > points_;
    };
    
    Weight_Function_Integration(std::vector<std::shared_ptr<Basis_Function> > const &bases,
                                std::vector<std::shared_ptr<Weight_Function> > const &weights,
                                std::shared_ptr<Solid_Geometry> solid_geometry,
                                std::vector<std::vector<double> > boundary_limits,
                                std::vector<int> num_intervals);
    
private:
    
    Weight_Function::Options options_;
    std::vector<std::shared_ptr<Basis_Function> > bases_
    std::vector<std::shared_ptr<Weight_Function> > weights_;
    std::shared_ptr<Solid_Geometry> solid_;
    Mesh mesh_;
    std::vector<std::vector<double> > limits_;
};

#endif
