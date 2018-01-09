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
                          int ord) const;
    
    // Return base index of spherical harmonics function
    virtual std::vector<int> const &scattering_indices() const
    {
        return l_indices_;
    }

    // Return the spherical harmonic indices
    virtual std::vector<int> const &harmonic_degrees() const
    {
        return l_indices_;
    }
    virtual std::vector<int> const &harmonic_orders() const
    {
        return m_indices_;
    }
    
    // Output data to XML file
    virtual void output(XML_Node output_node) const = 0;

    // Get reflected direction for ordinate direction
    virtual int reflect_ordinate(int o,
                                 std::vector<double> const &normal) const = 0;

    // Check class member sizes
    virtual void check_class_invariants() const = 0;
    
    // Perform moment-to-discrete operation
    virtual void moment_to_discrete(std::vector<double> &data) const;

    // Perform discrete-to-moment operation
    virtual void discrete_to_moment(std::vector<double> &data) const;

    // Get list of integration coefficients for a manufactured solution
    // (2l'+1)/(4\pi) \int \Omega Y_m1 Y_m2 d\Omega
    // Indexing:
    //     size[source moment]
    //     indices[source moment][moment num]
    //     coefficients[source moment][moment num][dimension]
    virtual void manufactured_coefficients(std::vector<int> &size,
                                           std::vector<std::vector<int> > &indices,
                                           std::vector<std::vector<double> > &coefficients) const;
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
