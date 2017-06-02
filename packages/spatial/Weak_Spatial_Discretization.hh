#ifndef Weak_Spatial_Discretization_hh
#define Weak_Spatial_Discretization_hh

#include "Spatial_Discretization.hh"
#include "Weight_Function.hh"

class Basis_Function;
class KD_Tree;

class Weak_Spatial_Discretization : public Spatial_Discretization
{
public:

    // Constructor
    Weak_Spatial_Discretization(std::vector<std::shared_ptr<Basis_Function> > &bases,
                                std::vector<std::shared_ptr<Weight_Function> > &weights,
                                std::shared_ptr<Dimensional_Moments> dimensional_moments,
                                std::shared_ptr<Weak_Spatial_Discretization_Options> options,
                                std::shared_ptr<KD_Tree> kd_tree = std::shared_ptr<KD_Tree>());
    
    // Point functions
    virtual bool has_reflection() const override
    {
        return has_reflection_;
    }
    virtual int number_of_points() const override
    {
        return number_of_points_;
    }
    virtual int number_of_boundary_points() const override
    {
        return number_of_boundary_bases_;
    }
    virtual int dimension() const override
    {
        return dimension_;
    }
    virtual int number_of_nodes() const override
    {
        return number_of_nodes_;
    }
    virtual int number_of_dimensional_moments() const override
    {
        return number_of_dimensional_moments_;
    }
    virtual std::shared_ptr<Dimensional_Moments> dimensional_moments() const override
    {
        return dimensional_moments_;
    }
    virtual std::shared_ptr<Point> point(int point_index) const override
    {
        return weights_[point_index];
    }
    virtual void output(XML_Node output_node) const override;
    virtual void check_class_invariants() const override;

    // Weight_Function functions
    virtual shared_ptr<Weak_Spatial_Discretization_Options> options() const
    {
        return options_;
    }
    virtual int number_of_boundary_weights() const
    {
        return number_of_boundary_weights_;
    }
    virtual std::vector<int> const &number_of_basis_functions() const
    {
        return number_of_basis_functions_;
    }
    virtual std::shared_ptr<Weight_Function> weight(int point_index) const
    {
        return weights_[point_index];
    }
    virtual std::shared_ptr<Basis_Function> basis(int point_index) const
    {
        return bases_[point_index];
    }
    virtual std::shared_ptr<Weight_Function> boundary_weight(int boundary_index) const
    {
        return boundary_weights_[boundary_index];
    }
    virtual std::shared_ptr<Basis_Function> boundary_basis(int boundary_index) const
    {
        return boundary_bases_[boundary_index];
    }

    // Functions to get values given expansion coefficients
    
    // Get the nearest weight function to a point: this can be used to find the basis functions applicable to a point
    virtual int nearest_point(std::vector<double> const &position) const;
    
    // Get the basis expansion values at the centers of the weight functions
    virtual double collocation_value(int i,
                                     std::vector<double> const &coefficients) const;
    virtual void collocation_values(std::vector<double> const &coefficients,
                                    std::vector<double> &values) const;
    
    // Get the weighted basis expansion values for each weight function
    virtual double weighted_collocation_value(int i,
                                              std::vector<double> const &coefficients) const;
    virtual void weighted_collocation_values(std::vector<double> const &coefficients,
                                             std::vector<double> &values) const;
    
    // Get expansion values at arbitrary points given the coefficients
    virtual double expansion_value(int i,
                                   std::vector<double> const &position,
                                   std::vector<double> const &coefficients) const;
    virtual double expansion_value(std::vector<double> const &position,
                                   std::vector<double> const &coefficients) const;
    
    // Get expansion values at arbitrary points given the group-dependent coefficients
    // Note that since the point is the last index, the groups could refer to energy groups plus moments
    virtual std::vector<double> expansion_values(int i,
                                                 int number_of_groups,
                                                 std::vector<double> const &position,
                                                 std::vector<double> const &coefficients) const;
    virtual std::vector<double> expansion_values(int number_of_groups,
                                                 std::vector<double> const &position,
                                                 std::vector<double> const &coefficients) const;

private:

    // Data
    bool has_reflection_;
    bool include_supg_;
    int number_of_points_;
    int number_of_boundary_weights_;
    int number_of_boundary_bases_;
    int dimension_;
    int number_of_nodes_;
    std::shared_ptr<Weak_Spatial_Discretization_Options> options_;
    std::vector<int> number_of_basis_functions_;
    std::vector<std::shared_ptr<Weight_Function> > weights_;
    std::vector<std::shared_ptr<Weight_Function> > boundary_weights_;
    std::vector<std::shared_ptr<Basis_Function> > bases_;
    std::vector<std::shared_ptr<Basis_Function> > boundary_bases_;
    std::shared_ptr<Dimensional_Moments> dimensional_moments_;
    std::shared_ptr<KD_Tree> kd_tree_;
};

struct Weak_Spatial_Discretization_Options
{
    // Main value to set: method of weighting
    enum class Weighting
    {
        POINT,
        WEIGHT,
        FLUX
    };
    std::shared_ptr<Conversion<Weighting, std::string> > weighting_conversion() const;
        
    // Total cross section method: changed to ISOTROPIC
    // unless Weighting::FLUX is used
    enum class Total
    {
        ISOTROPIC,
        MOMENT
    };
    std::shared_ptr<Conversion<Total, std::string> > total_conversion() const;

    // Tau scaling according to distance from boundary
    enum class Tau_Scaling
    {
        NONE,
        FUNCTIONAL, // 1 - b(boundary) - b(center)
        LINEAR, // pos_boundary / radius
        ABSOLUTE // 0 if on boundary
    };
    std::shared_ptr<Conversion<Tau_Scaling, std::string> > tau_scaling_conversion() const;

    // Galerkin: should be set to true or false by time options are passed to Weight_Function
    enum class Identical_Basis_Functions
    {
        AUTO, // Set in creation methods
        TRUE,
        FALSE
    };
    std::shared_ptr<Conversion<Identical_Basis_Functions, std::string> > identical_basis_functions_conversion() const;
        
    // External integration parameters: must be set for external calculation
    bool external_integral_calculation = true;
    int integration_ordinates = 8; // Dimensional integration quadrature
    std::vector<std::vector<double> > limits;
    std::shared_ptr<Solid_Geometry> solid;
    std::vector<int> dimensional_cells;
    
    // Parameters for the user to set
    bool include_supg = false;
    Identical_Basis_Functions identical_basis_functions = Identical_Basis_Functions::FALSE;
    Weighting weighting = Weighting::WEIGHT; 
    Total total = Total::ISOTROPIC;
    Tau_Scaling tau_scaling = Tau_Scaling::LINEAR;
    std::function<double(int /*moment*/,
                         int /*group*/,
                         std::vector<double> const & /*position*/)> flux;

    // Automatically set parameters
    bool input_finalized = false;
    bool normalized = true;

    // Check input and set automatic parameters
    void finalize_input();
};

#endif
