#ifndef Angular_Discretization_hh
#define Angular_Discretization_hh

#include <vector>

class XML_Node;

/*
  Pure virtual class for an S_N angular discretization
*/
class Angular_Discretization
{
public:

    // Constructor
    Angular_Discretization(int dimension,
                           int number_of_scattering_moments,
                           int number_of_ordinates);

    // Number of spatial dimensions
    virtual int dimension()
    {
        return dimension_;
    }
    
    // Number of angular moments
    virtual int number_of_moments()
    {
        return number_of_moments_;
    }

    // Number of scattering moments
    virtual int number_of_scattering_moments()
    {
        return number_of_scattering_moments_;
    }
    
    // Number of discrete ordinates
    virtual int number_of_ordinates()
    {
        return number_of_ordinates_;
    }

    // Angular normalization factor: 2 for 1D and otherwise 4\pi
    virtual double angular_normalization()
    {
        return angular_normalization_;
    }

    // Return a direction for a single ordinate
    virtual std::vector<double> const &direction(int ord) const = 0;
    
    // Return std::vector of ordinates
    virtual std::vector<double> const &ordinates() const = 0;

    // Return std::vector of weights
    virtual std::vector<double> const &weights() const = 0;

    // Find a value of the moment in a particular ordinate direction
    virtual double moment(int mom,
                          int ord);
    
    // Return base index of spherical harmonics function
    virtual std::vector<int> const &scattering_indices() const
    {
        return l_indices_;
    }

    // Output data to XML file
    virtual void output(XML_Node output_node) const = 0;

    virtual int reflect_ordinate(int o,
                                 std::vector<double> const &normal) const = 0;

    virtual void check_class_invariants() const = 0;
    
protected:
    
    int dimension_;
    int number_of_scattering_moments_;
    int number_of_moments_;
    int number_of_ordinates_;
    double angular_normalization_;

private:

    virtual void initialize_moment_data();
    
    std::vector<int> l_indices_;
    std::vector<int> m_indices_;
};

#endif
