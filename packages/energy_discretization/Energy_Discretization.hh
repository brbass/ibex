#ifndef Energy_Discretization_hh
#define Energy_Discretization_hh

#include <vector>

class XML_Node;

/* 
   Simple multigroup representation of energy
*/
class Energy_Discretization
{
public:

    // Constructors
    Energy_Discretization(int number_of_groups,
                          std::vector<double> const &energy_bounds);
    Energy_Discretization(int number_of_groups);
    
    // Number of energy groups
    int number_of_groups()
    {
        return number_of_groups_;
    }

    // Bounds of energy groups
    std::vector<double> const &energy_bounds() const
    {
        return energy_bounds_;
    }

    // Check class invariants
    void check_class_invariants() const;

    // Output data to XML file
    void output(XML_Node output_node) const;
    
private:
    
    int number_of_groups_;
    std::vector<double> energy_bounds_;

};

#endif
