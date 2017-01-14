#ifndef Weight_Function_hh
#define Weight_Function_hh

#include "Point.hh"

class Basis_Function;
class Cartesian_Plane;
class RBF_Function;

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
                    shared_ptr<Meshless_Function> meshless_function,
                    vector<shared_ptr<Basis_Function> > basis_functions,
                    shared_ptr<Solid_Geometry> solid_geometry,
                    vector<shared_ptr<Cartesian_Plane> > boundary_surfaces);

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
    shared_ptr<Meshless_Function> function() const
    {
        return meshless_function_;
    }
    shared_ptr<Basis_Function> basis_function(int i) const
    {
        return basis_functions_[i];
    }
    shared_ptr<Solid_Geometry> solid_geometry() const
    {
        return solid_geometry_;
    }
    shared_ptr<Cartesian_Plane> boundary_surface(int i) const
    {
        return boundary_surfaces_[i];
    }
    
    // Quadrature methods
    virtual void get_full_quadrature(vector<double> &quadrature,
                                     vector<double> &weight) const;
    virtual void get_basis_quadrature(int i,
                                      vector<double> &quadrature,
                                      vector<double> &weight) const;
    
private:
    
    // Integration methods
    virtual void calculate_integrals();
    
    // Data
    int index_;
    int dimension_;
    int integration_ordinates_;
    int number_of_basis_functions_;
    int number_of_boundary_surfaces_;
    shared_ptr<Meshless_Function> meshless_function_;
    vector<shared_ptr<Basis_Function> > basis_functions_;
    shared_ptr<Solid_Geometry> solid_geometry_;
    vector<shared_ptr<Cartesian_Plane> > boundary_surfaces_;

    // Calculated data
    
}

#endif
