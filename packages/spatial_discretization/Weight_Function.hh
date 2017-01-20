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
    
    enum class Material_Type
    {
        STANDARD,
        SUPG,
        NONE
    };

    struct Material_Options
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

        // Scattering cross section method: changed to SCATTERING_MOMENTS
        // unless Weighting::FLUX is used
        enum class Scattering
        {
            SCATTERING_MOMENTS,
            MOMENTS
        };

        Weighting weighting = Weighting::WEIGHT;
        Total total = Total::ISOTROPIC;
        Scattering scattering = Scattering::SCATTERING_MOMENTS;
        Output output = Output::STANDARD;
        std::function<double(int /*moment*/,
                             int /*group*/,
                             std::vector<double> const & /*position*/)> flux;
    };
    
    // Constructor
    Weight_Function(int index,
                    int dimension,
                    int integration_ordinates,
                    Material_Options material_options,
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
    virtual Point_Type point_type() const override
    {
        return point_type_;
    }
    virtual std::shared_ptr<Material> material() const override
    {
        return material_[0];
    }
    virtual std::vector<std::shared_ptr<Material> > supg_material() const
    {
        return material_;
    }
    virtual std::vector<double> const &position() const override
    {
        return position_;
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
        return boundary_surfaces_[i];
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
    virtual void calculate_integrals();
    virtual void calculate_material();
    
    // Point data
    int index_;
    int dimension_;
    Point_Type point_type_;
    std::vector<double> position_;

    // Weight_Function data
    int integration_ordinates_;
    int number_of_basis_functions_;
    int number_of_boundary_surfaces_;
    double radius_;
    Material_Options material_options_;
    std::shared_ptr<Meshless_Function> meshless_function_;
    std::vector<std::shared_ptr<Basis_Function> > basis_functions_;
    std::shared_ptr<Solid_Geometry> solid_geometry_;
    std::vector<std::shared_ptr<Cartesian_Plane> > boundary_surfaces_;

    // Calculated data
    std::vector<double> min_boundary_limits_;
    std::vector<double> max_boundary_limits_;
    
    // Integrals
    std::vector<double> is_b_w_;
    std::vector<double> iv_b_w_;
    std::vector<double> iv_b_dw_;
    std::vector<double> iv_db_w_;
    std::vector<double> iv_db_dw_;

    // Material
    std::vector<std::shared_ptr<Material> > material_;
};

#endif
