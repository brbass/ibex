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
    virtual vector<double> const &direction(int ord) const override
    {
        return directions_[ord];
    }
    
    // Return all ordinates
    virtual vector<double> const &ordinates() const override
    {
        return ordinates_;
    }

    // Return weights
    virtual vector<double> const &weights() const override
    {
        return weights_;
    }
    
    // X component of ordinates
    virtual vector<double> const &mu() const
    {
        return mu_;
    }

    // Y component of ordinates
    virtual vector<double> const &eta() const
    {
        return mu_;
    }
    
    // Z component of ordinates
    virtual vector<double> const &xi() const
    {
        return mu_;
    }
    
    // Check class invariants
    virtual void check_class_invariants() const override;

    // Output data to XML file
    virtual void output(pugi::xml_node &output_node) const override;

    // Get reflected directions
    virtual int reflect_ordinate(int o,
                                 vector<double> const &normal) const override;

private:

    void initialize_quadrature();
    void initialize_1();
    void initialize_2();
    void initialize_3();
    
    int rule_;

    double reflection_tolerance_;

    vector<double> mu_;
    vector<double> eta_;
    vector<double> xi_;
    vector<double> ordinates_;
    vector<double> weights_;
    vector<vector<double> > directions_;
};

#endif
