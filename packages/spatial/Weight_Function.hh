#ifndef Weight_Function_hh
#define Weight_Function_hh

#include "Point.hh"

#include <functional>

class Basis_Function;
class Cartesian_Plane;
class Meshless_Function;
class Solid_Geometry;

class Weight_Function : public Point
{
public:
    
    struct Options
    {
    public:

        // Main value to set: method of weighting
        enum class Weighting
        {
            POINT,
            WEIGHT,
            FLUX
        };
        
        // Determines whether to add the dimensional material moments
        enum class Output
        {
            STANDARD,
            SUPG
        };
        
        // Total cross section method: changed to ISOTROPIC
        // unless Weighting::FLUX is used
        enum class Total
        {
            ISOTROPIC,
            MOMENT
        };

        // Parameters that are automatically set
        bool include_supg = false;
        bool normalized = true;
        double tau; // SUPG parameter (tau_const / shape)
        bool outside_integral_calculation = false;

        // Parameters for the user to set
        int integration_ordinates = 32;
        double tau_const = 1; // Constant in front of 1/shape
        Weighting weighting = Weighting::WEIGHT;
        Total total = Total::ISOTROPIC;
        Output output = Output::STANDARD;
        std::function<double(int /*moment*/,
                             int /*group*/,
                             std::vector<double> const & /*position*/)> flux;
    };
    
    // Constructor
    Weight_Function(int index,
                    int dimension,
                    Options options,
                    std::shared_ptr<Meshless_Function> meshless_function,
                    std::vector<std::shared_ptr<Basis_Function> > basis_functions,
                    std::shared_ptr<Solid_Geometry> solid_geometry,
                    std::vector<std::shared_ptr<Cartesian_Plane> > boundary_surfaces);
    
    // Point functions
    virtual int index() const override
    {
        return index_;
    }
    virtual int dimension() const override
    {
        return dimension_;
    }
    virtual int number_of_nodes() const override
    {
        return 1;
    }
    virtual Point_Type point_type() const override
    {
        return point_type_;
    }
    virtual std::vector<double> const &position() const override
    {
        return position_;
    }
    virtual std::shared_ptr<Material> material() const override
    {
        return material_;
    }
    virtual void output(XML_Node output_node) const override;
    virtual void check_class_invariants() const override;

    // Weight_Function functions
    virtual int number_of_basis_functions() const
    {
        return number_of_basis_functions_;
    }
    virtual int number_of_boundary_surfaces() const
    {
        return number_of_boundary_surfaces_;
    }
    virtual int number_of_dimensional_moments() const
    {
        return number_of_dimensional_moments_;
    }
    virtual double radius() const
    {
        return radius_;
    }
    virtual Options options() const
    {
        return options_;
    }
    virtual std::vector<int> const &basis_function_indices() const
    {
        return basis_function_indices_;
    }
    virtual std::shared_ptr<Meshless_Function> function() const
    {
        return meshless_function_;
    }
    virtual std::shared_ptr<Basis_Function> basis_function(int i) const
    {
        return basis_functions_[i];
    }
    virtual std::shared_ptr<Solid_Geometry> solid_geometry() const
    {
        return solid_geometry_;
    }
    virtual std::shared_ptr<Cartesian_Plane> boundary_surface(int i) const
    {
        return weighted_boundary_surfaces_[i];
    }
    
    // Collocation values
    virtual std::vector<double> const &v_b()
    {
        return v_b_;
    }
    virtual std::vector<double> const &v_db()
    {
        return v_db_;
    }
    
    // Integral values
    virtual std::vector<double> const &is_w()
    {
        return is_w_;
    }
    virtual std::vector<double> const &is_b_w()
    {
        return is_b_w_;
    }
    virtual std::vector<double> const &iv_w()
    {
        return iv_w_;
    }
    virtual std::vector<double> const &iv_dw()
    {
        return iv_dw_;
    }
    virtual std::vector<double> const &iv_b_w()
    {
        return iv_b_w_;
    }
    virtual std::vector<double> const &iv_b_dw()
    {
        return iv_b_dw_;
    }
    virtual std::vector<double> const &iv_db_w()
    {
        return iv_db_w_;
    }
    virtual std::vector<double> const &iv_db_dw()
    {
        return iv_db_dw_;
    }
    
    // Quadrature methods
    virtual bool get_full_quadrature(std::vector<std::vector<double> > &ordinates,
                                     std::vector<double> &weights) const;
    virtual bool get_basis_quadrature(int i,
                                      std::vector<std::vector<double> > &ordinates,
                                      std::vector<double> &weights) const;
    virtual bool get_full_surface_quadrature(int s,
                                             std::vector<std::vector<double> > &ordinates,
                                             std::vector<double> &weights) const;
    virtual bool get_basis_surface_quadrature(int i,
                                              int s,
                                              std::vector<std::vector<double> > &ordinates,
                                              std::vector<double> &weights) const;
    
private:

    // Specific quadrature methods
    virtual bool get_full_quadrature_1d(std::vector<std::vector<double> > &ordinates,
                                        std::vector<double> &weights) const;
    virtual bool get_full_quadrature_2d(std::vector<std::vector<double> > &ordinates,
                                        std::vector<double> &weights) const;
    virtual bool get_basis_quadrature_1d(int i,
                                         std::vector<std::vector<double> > &ordinates,
                                         std::vector<double> &weights) const;
    virtual bool get_basis_quadrature_2d(int i,
                                         std::vector<std::vector<double> > &ordinates,
                                         std::vector<double> &weights) const;
    virtual bool get_full_surface_quadrature_2d(int s,
                                                std::vector<std::vector<double> > &ordinates,
                                                std::vector<double> &weights) const;
    virtual bool get_basis_surface_quadrature_2d(int i,
                                                 int s,
                                                 std::vector<std::vector<double> > &ordinates,
                                                 std::vector<double> &weights) const;
    
    // Integration methods
    virtual void set_options_and_limits();
    virtual void calculate_values();
    virtual void calculate_integrals();
    virtual void calculate_material();
    virtual void calculate_boundary_source();
    
    // Specific integration methods
    virtual void calculate_standard_point_material();
    virtual void calculate_standard_weight_material();
    virtual void calculate_supg_point_material();
    virtual void calculate_supg_weight_material();
    virtual void calculate_weight_boundary_source();
    
    // Point data
    int index_;
    int dimension_;
    Point_Type point_type_;
    std::vector<double> position_;
    std::shared_ptr<Material> material_;

    // Weight_Function data
    int number_of_basis_functions_;
    int number_of_boundary_surfaces_;
    int number_of_dimensional_moments_;
    double radius_;
    Options options_;
    std::vector<int> basis_function_indices_;
    std::shared_ptr<Meshless_Function> meshless_function_;
    std::vector<std::shared_ptr<Basis_Function> > basis_functions_;
    std::shared_ptr<Solid_Geometry> solid_geometry_;
    std::vector<std::shared_ptr<Cartesian_Plane> > boundary_surfaces_;
    std::vector<std::shared_ptr<Cartesian_Plane> > weighted_boundary_surfaces_;
    
    // Calculated data
    std::vector<double> min_boundary_limits_;
    std::vector<double> max_boundary_limits_;

    // Values
    std::vector<double> v_b_; // basis function at weight center
    std::vector<double> v_db_; // derivative of basis function at weight pcenter 
    
    // Surface integrals
    std::vector<double> is_w_; // weight function: s
    std::vector<double> is_b_w_; // weight/basis functions: s->i

    // Volume integrals
    std::vector<double> iv_w_; // weight function: none
    std::vector<double> iv_dw_; // derivative of weight function: d
    std::vector<double> iv_b_w_; // basis function and weight function: i
    std::vector<double> iv_b_dw_; // basis function and derivative of weight function: dw->i
    std::vector<double> iv_db_w_; // weight function and derivative of basis function: db->i
    std::vector<double> iv_db_dw_; // derivative of basis and weight functions: db->dw->i
};

#endif
