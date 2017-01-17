#ifndef Weight_Function_hh
#define Weight_Function_hh

#include "Point.hh"

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
    
    // Constructor
    Weight_Function(int index,
                    int dimension,
                    int integration_ordinates,
                    std::shared_ptr<Meshless_Function> meshless_function,
                    std::vector<std::shared_ptr<Basis_Function> > basis_functions,
                    std::shared_ptr<Solid_Geometry> solid_geometry,
                    std::vector<std::shared_ptr<Cartesian_Plane> > boundary_surfaces);

    // Data access
    virtual int index() const
    {
        return index_;
    }
    virtual int dimension() const
    {
        return dimension_;
    }
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
    virtual void get_full_quadrature_1d(std::vector<double> &ordinates,
                                        std::vector<double> &weights) const;
    virtual void get_basis_quadrature_1d(int i,
                                         std::vector<double> &ordinates,
                                         std::vector<double> &weights) const;
    
private:
    
    // Integration methods
    virtual void calculate_integrals_1d() const;
    
    // Data
    int index_;
    int dimension_;
    int integration_ordinates_;
    int number_of_basis_functions_;
    int number_of_boundary_surfaces_;
    std::shared_ptr<Meshless_Function> meshless_function_;
    std::vector<std::shared_ptr<Basis_Function> > basis_functions_;
    std::shared_ptr<Solid_Geometry> solid_geometry_;
    std::vector<std::shared_ptr<Cartesian_Plane> > boundary_surfaces_;

    // Calculated data
    
};

#endif
