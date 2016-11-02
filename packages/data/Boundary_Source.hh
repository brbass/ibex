#ifndef Boundary_Source_hh
#define Boundary_Source_hh

#include <memory>
#include <vector>

#include "Angular_Discretization.hh"
#include "Energy_Discretization.hh"

using std::shared_ptr;
using std::vector;

class Boundary_Source
{
public:
    
    Boundary_Source(int index,
                    shared_ptr<Angular_Discretization> angular_discretization,
                    shared_ptr<Energy_Discretization> energy_discretization,
                    vector<double> const &boundary_source,
                    vector<double> const &alpha);
    
    int index() const
    {
        return index_;
    }
    
    vector<double> const &boundary_source() const
    {
        return boundary_source_;
    }

    vector<double> const &alpha() const
    {
        return alpha_;
    }
    
    bool has_reflection() const;
    
    void check_class_invariants() const;
    
    void output(pugi::xml_node &output_node) const;

private:

    int index_;
    
    shared_ptr<Angular_Discretization> angular_discretization_;
    shared_ptr<Energy_Discretization> energy_discretization_;

    vector<double> boundary_source_;
    vector<double> alpha_;
};

#endif
