#ifndef Boundary_Source_hh
#define Boundary_Source_hh

#include <memory>
#include <vector>

class Angular_Discretization;
class Energy_Discretization;
class XML_Node;

class Boundary_Source
{
public:

    Boundary_Source(int index,
                    std::shared_ptr<Angular_Discretization> angular_discretization,
                    std::shared_ptr<Energy_Discretization> energy_discretization,
                    std::vector<double> const &boundary_source,
                    std::vector<double> const &alpha);
    
    int index() const
    {
        return index_;
    }
    
    std::vector<double> const &boundary_source() const
    {
        return boundary_source_;
    }

    std::vector<double> const &alpha() const
    {
        return alpha_;
    }
    
    bool has_reflection() const;
    
    void check_class_invariants() const;
    
    void output(XML_Node output_node) const;

private:

    int index_;
    
    std::shared_ptr<Angular_Discretization> angular_discretization_;
    std::shared_ptr<Energy_Discretization> energy_discretization_;

    std::vector<double> boundary_source_;
    std::vector<double> alpha_;
};

#endif
