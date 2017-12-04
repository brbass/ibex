#ifndef LDFE_Quadrature_hh
#define LDFE_Quadrature_hh

#include "Angular_Discretization.hh"

/*
  Quadrature based on linear discontinuous finite elements
  
  Ref: Discrete-ordinates quadrature sets based on linear discontinuous finite elements
  Authors: Jarrell and Adams
*/
class LDFE_Quadrature : public Angular_Discretization
{
public:
    
    // Constructor
    LDFE_Quadrature(int dimension,
                    int number_of_scattering_moments,
                    int rule);

    // Return direction
    virtual std::vector<double> const &direction(int ord) const override
    {
        return directions_[ord];
    }
    
    // Return all ordinates
    virtual std::vector<double> const &ordinates() const override
    {
        return ordinates_;
    }

    // Return weights
    virtual std::vector<double> const &weights() const override
    {
        return weights_;
    }
    
    // X component of ordinates
    virtual std::vector<double> const &mu() const
    {
        return mu_;
    }

    // Y component of ordinates
    virtual std::vector<double> const &eta() const
    {
        return mu_;
    }
    
    // Z component of ordinates
    virtual std::vector<double> const &xi() const
    {
        return mu_;
    }
    
    // Check class invariants
    virtual void check_class_invariants() const override;

    // Output data to XML file
    virtual void output(XML_Node output_node) const override;
    
    // Get reflected directions
    virtual int reflect_ordinate(int o,
                                 std::vector<double> const &normal) const override;

private:

    void get_octant_values(int rule,
                           std::vector<double> &values) const;
    void initialize_quadrature();
    
    int rule_;
    double reflection_tolerance_;

    std::vector<double> mu_;
    std::vector<double> eta_;
    std::vector<double> xi_;
    std::vector<double> ordinates_;
    std::vector<double> weights_;
    std::vector<std::vector<double> > directions_;
};

#endif
